#include "../copyright.h"
/*==========================================================================================*/
/*! \file bvals_particle_self_gravity.c
 *  \brief Sets boundary conditions (quantities in ghost zones) for the
 *   particle gravitational potential on each edge of a Grid. 
 *
 * PURPOSE: Sets boundary conditions (quantities in ghost zones) for the
 *   particle gravitational potential on each edge of a Grid.  See comments at
 *   start of bvals_mhd.c for more details.
 * The only BC functions implemented here are for:
 *- 1 = reflecting, 4 = periodic, and MPI boundaries
 *
 * CONTAINS PUBLIC FUNCTIONS: 
 * - bvals_particle_self_gravity()      - calls appropriate functions to set ghost cells
 * - bvals_particle_self_gravity_init() - sets function pointers used by bvals_grav()
 * - bvals_grav_fun()                   - enrolls a pointer to a user-defined BC function */
/*==========================================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include "../defs.h"
#include "../athena.h"
#include "prototypes.h"
#include "../prototypes.h"
#include "../globals.h"

/* The functions in this file will only work with PARTICLE_SELF_GRAVITY */
#ifdef PARTICLE_SELF_GRAVITY

/*! \struct Remap
 *  \brief Define structure which holds variables remapped 
 *  by shearing sheet BCs */
typedef struct Remap_s{
  Real par_Phi;
}Remap;

/* The memory for all the arrays below is allocated in bvals_particle_self_gravity_shear_init */
/* Arrays of ghost zones containing remapped conserved quantities */
static Remap ***GhstZns=NULL, ***GhstZnsBuf=NULL;
/* 1D vectors for reconstruction in conservative remap step */
static Real *U=NULL, *Flx=NULL;
/* temporary vector needed for 3rd order reconstruction in ghost zones */
#if defined(THIRD_ORDER_CHAR) || defined(THIRD_ORDER_PRIM)
static Real *Uhalf=NULL;
#endif

#ifdef MPI_PARALLEL
/* MPI send and receive buffers */
static double **send_buf = NULL, **recv_buf = NULL;
static double *send_bufss = NULL, *recv_bufss = NULL;
static MPI_Request *recv_rq, *send_rq;
#endif /* MPI_PARALLEL */

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 *   reflect_Phi_???()  - apply reflecting BCs at boundary ???
 *   periodic_Phi_???() - apply periodic BCs at boundary ???
 *   pack_Phi_???()   - pack data for MPI non-blocking send at ??? boundary
 *   unpack_Phi_???() - unpack data for MPI non-blocking receive at ??? boundary
 *   RemapFluxPSG()   - Does conservative remap of fluxes
 *   ShearingSheet_part_grav_??() - Shearing sheet boundaries for potential 
 *============================================================================*/

void RemapFluxPSG(const Real *U,const Real eps,const int ji,const int jo, Real *F);
void ShearingSheet_par_grav_ix1(DomainS *pD);
void ShearingSheet_par_grav_ox1(DomainS *pD);

static void reflect_Phi_ix1(GridS *pG);
static void reflect_Phi_ox1(GridS *pG);
static void reflect_Phi_ix2(GridS *pG);
static void reflect_Phi_ox2(GridS *pG);
static void reflect_Phi_ix3(GridS *pG);
static void reflect_Phi_ox3(GridS *pG);

static void periodic_Phi_ix1(GridS *pG);
static void periodic_Phi_ox1(GridS *pG);
static void periodic_Phi_ix2(GridS *pG);
static void periodic_Phi_ox2(GridS *pG);
static void periodic_Phi_ix3(GridS *pG);
static void periodic_Phi_ox3(GridS *pG);

static void obc_fft_Phi_ix1(GridS *pG);
static void obc_fft_Phi_ox1(GridS *pG);
static void obc_fft_Phi_ix2(GridS *pG);
static void obc_fft_Phi_ox2(GridS *pG);
static void obc_fft_Phi_ix3(GridS *pG);
static void obc_fft_Phi_ox3(GridS *pG);

static void ProlongateLater(GridS *pG);

#ifdef MPI_PARALLEL
static void pack_Phi_ix1(GridS *pG);
static void pack_Phi_ox1(GridS *pG);
static void pack_Phi_ix2(GridS *pG);
static void pack_Phi_ox2(GridS *pG);
static void pack_Phi_ix3(GridS *pG);
static void pack_Phi_ox3(GridS *pG);

static void unpack_Phi_ix1(GridS *pG);
static void unpack_Phi_ox1(GridS *pG);
static void unpack_Phi_ix2(GridS *pG);
static void unpack_Phi_ox2(GridS *pG);
static void unpack_Phi_ix3(GridS *pG);
static void unpack_Phi_ox3(GridS *pG);
#endif /* MPI_PARALLEL */

/*=========================== PUBLIC FUNCTIONS ==============================*/
/*----------------------------------------------------------------------------*/
/*! \fn void bvals_particle_self_gravity(DomainS *pD)
 *  \brief Calls appropriate functions to set ghost zones.  
 *
 *   The function
 *   pointers (*_GBCFun) are set during initialization by 
 *   bvals_particle_self_gravity_init() to be one of the 
 *   functions corresponding to reflecting or periodic.  If the
 *   left- or right-Grid ID numbers are >= 1 (neighboring grids exist), then MPI
 *   calls are used.
 *
 * Order for updating boundary conditions must always be x1-x2-x3 in order to
 * fill the corner cells properly
 */

void bvals_particle_self_gravity(DomainS *pD)
{
  GridS *pGrid = (pD->Grid);
#ifdef SHEARING_BOX
  int myL,myM,myN;
#endif
#ifdef MPI_PARALLEL
  int cnt1, cnt2, cnt3, cnt, ierr, mIndex;
#endif /* MPI_PARALLEL */

/*--- Step 1. ------------------------------------------------------------------
 * Boundary Conditions in x1-direction */

  if (pGrid->Nx[0] > 1){

#ifdef MPI_PARALLEL

    cnt = nghost*(pGrid->Nx[1])*(pGrid->Nx[2]);

/* MPI blocks to both left and right */
    if (pGrid->rx1_Gid >= 0 && pGrid->lx1_Gid >= 0) {
      /* Post non-blocking receives for data from L and R Grids */
      ierr = MPI_Irecv(&(recv_buf[0][0]),cnt,MPI_DOUBLE,pGrid->lx1_Gid,LtoR_tag,
        pD->Comm_Domain, &(recv_rq[0]));
      ierr = MPI_Irecv(&(recv_buf[1][0]),cnt,MPI_DOUBLE,pGrid->rx1_Gid,RtoL_tag,
        pD->Comm_Domain, &(recv_rq[1]));

      /* pack and send data L and R */
      pack_Phi_ix1(pGrid);
      ierr = MPI_Isend(&(send_buf[0][0]),cnt,MPI_DOUBLE,pGrid->lx1_Gid,RtoL_tag,
        pD->Comm_Domain, &(send_rq[0]));

      pack_Phi_ox1(pGrid);
      ierr = MPI_Isend(&(send_buf[1][0]),cnt,MPI_DOUBLE,pGrid->rx1_Gid,LtoR_tag,
        pD->Comm_Domain, &(send_rq[1]));


      /* check non-blocking sends have completed. */
      ierr = MPI_Waitall(2, send_rq, MPI_STATUS_IGNORE);

      /* check non-blocking receives and unpack data in any order. */
      ierr = MPI_Waitany(2,recv_rq,&mIndex,MPI_STATUS_IGNORE);
      if (mIndex == 0) unpack_Phi_ix1(pGrid);
      if (mIndex == 1) unpack_Phi_ox1(pGrid);
      ierr = MPI_Waitany(2,recv_rq,&mIndex,MPI_STATUS_IGNORE);
      if (mIndex == 0) unpack_Phi_ix1(pGrid);
      if (mIndex == 1) unpack_Phi_ox1(pGrid);

    }

/* Physical boundary on left, MPI block on right */
    if (pGrid->rx1_Gid >= 0 && pGrid->lx1_Gid < 0) {

      /* Post non-blocking receive for data from R Grid */
      ierr = MPI_Irecv(&(recv_buf[1][0]),cnt,MPI_DOUBLE,pGrid->rx1_Gid,RtoL_tag,
        pD->Comm_Domain, &(recv_rq[1]));

      /* pack and send data R */
      pack_Phi_ox1(pGrid);
      ierr = MPI_Isend(&(send_buf[1][0]),cnt,MPI_DOUBLE,pGrid->rx1_Gid,LtoR_tag,
        pD->Comm_Domain, &(send_rq[1]));

      /* set physical boundary */
      (*(pD->ix1_GBCFun))(pGrid);

      /* check non-blocking send has completed. */
      ierr = MPI_Wait(&(send_rq[1]), MPI_STATUS_IGNORE);

      /* wait on non-blocking receive from R and unpack data */
      ierr = MPI_Wait(&(recv_rq[1]), MPI_STATUS_IGNORE);
      unpack_Phi_ox1(pGrid);

    }

/* MPI block on left, Physical boundary on right */
    if (pGrid->rx1_Gid < 0 && pGrid->lx1_Gid >= 0) {

      /* Post non-blocking receive for data from L grid */
      ierr = MPI_Irecv(&(recv_buf[0][0]),cnt,MPI_DOUBLE,pGrid->lx1_Gid,LtoR_tag,
        pD->Comm_Domain, &(recv_rq[0]));

      /* pack and send data L */
      pack_Phi_ix1(pGrid);
      ierr = MPI_Isend(&(send_buf[0][0]),cnt,MPI_DOUBLE,pGrid->lx1_Gid,RtoL_tag,
        pD->Comm_Domain, &(send_rq[0]));

      /* set physical boundary */
      (*(pD->ox1_GBCFun))(pGrid);

      /* check non-blocking send has completed. */
      ierr = MPI_Wait(&(send_rq[0]), MPI_STATUS_IGNORE);

      /* wait on non-blocking receive from L and unpack data */
      ierr = MPI_Wait(&(recv_rq[0]), MPI_STATUS_IGNORE);
      unpack_Phi_ix1(pGrid);
    
    }
#endif /* MPI_PARALLEL */

/* Physical boundaries on both left and right */
    if (pGrid->rx1_Gid < 0 && pGrid->lx1_Gid < 0) {
      (*(pD->ix1_GBCFun))(pGrid);
      (*(pD->ox1_GBCFun))(pGrid);
    }

  }

/*--- Step 2. ------------------------------------------------------------------
 * Boundary Conditions in x2-direction */

  if (pGrid->Nx[1] > 1){

#ifdef MPI_PARALLEL

    cnt = (pGrid->Nx[0] + 2*nghost)*nghost*(pGrid->Nx[2]);

/* MPI blocks to both left and right */
    if (pGrid->rx2_Gid >= 0 && pGrid->lx2_Gid >= 0) {

      /* Post non-blocking receives for data from L and R Grids */
      ierr = MPI_Irecv(&(recv_buf[0][0]),cnt,MPI_DOUBLE,pGrid->lx2_Gid,LtoR_tag,
        pD->Comm_Domain, &(recv_rq[0]));
      ierr = MPI_Irecv(&(recv_buf[1][0]),cnt,MPI_DOUBLE,pGrid->rx2_Gid,RtoL_tag,
        pD->Comm_Domain, &(recv_rq[1]));

      /* pack and send data L and R */
      pack_Phi_ix2(pGrid);
      ierr = MPI_Isend(&(send_buf[0][0]),cnt,MPI_DOUBLE,pGrid->lx2_Gid,RtoL_tag,
        pD->Comm_Domain, &(send_rq[0]));

      pack_Phi_ox2(pGrid);
      ierr = MPI_Isend(&(send_buf[1][0]),cnt,MPI_DOUBLE,pGrid->rx2_Gid,LtoR_tag,
        pD->Comm_Domain, &(send_rq[1]));

      /* check non-blocking sends have completed. */
      ierr = MPI_Waitall(2, send_rq, MPI_STATUS_IGNORE);

      /* check non-blocking receives and unpack data in any order. */
      ierr = MPI_Waitany(2,recv_rq,&mIndex,MPI_STATUS_IGNORE);
      if (mIndex == 0) unpack_Phi_ix2(pGrid);
      if (mIndex == 1) unpack_Phi_ox2(pGrid);
      ierr = MPI_Waitany(2,recv_rq,&mIndex,MPI_STATUS_IGNORE);
      if (mIndex == 0) unpack_Phi_ix2(pGrid);
      if (mIndex == 1) unpack_Phi_ox2(pGrid);

    }

/* Physical boundary on left, MPI block on right */
    if (pGrid->rx2_Gid >= 0 && pGrid->lx2_Gid < 0) {

      /* Post non-blocking receive for data from R Grid */
      ierr = MPI_Irecv(&(recv_buf[1][0]),cnt,MPI_DOUBLE,pGrid->rx2_Gid,RtoL_tag,
        pD->Comm_Domain, &(recv_rq[1]));

      /* pack and send data R */
      pack_Phi_ox2(pGrid);
      ierr = MPI_Isend(&(send_buf[1][0]),cnt,MPI_DOUBLE,pGrid->rx2_Gid,LtoR_tag,
        pD->Comm_Domain, &(send_rq[1]));

      /* set physical boundary */
      (*(pD->ix2_GBCFun))(pGrid);

      /* check non-blocking send has completed. */
      ierr = MPI_Wait(&(send_rq[1]), MPI_STATUS_IGNORE);

      /* wait on non-blocking receive from R and unpack data */
      ierr = MPI_Wait(&(recv_rq[1]), MPI_STATUS_IGNORE);
      unpack_Phi_ox2(pGrid);

    }

/* MPI block on left, Physical boundary on right */
    if (pGrid->rx2_Gid < 0 && pGrid->lx2_Gid >= 0) {

      /* Post non-blocking receive for data from L grid */
      ierr = MPI_Irecv(&(recv_buf[0][0]),cnt,MPI_DOUBLE,pGrid->lx2_Gid,LtoR_tag,
        pD->Comm_Domain, &(recv_rq[0]));

      /* pack and send data L */
      pack_Phi_ix2(pGrid);
      ierr = MPI_Isend(&(send_buf[0][0]),cnt,MPI_DOUBLE,pGrid->lx2_Gid,RtoL_tag,
        pD->Comm_Domain, &(send_rq[0]));

      /* set physical boundary */
      (*(pD->ox2_GBCFun))(pGrid);

      /* check non-blocking send has completed. */
      ierr = MPI_Wait(&(send_rq[0]), MPI_STATUS_IGNORE);

      /* wait on non-blocking receive from L and unpack data */
      ierr = MPI_Wait(&(recv_rq[0]), MPI_STATUS_IGNORE);
      unpack_Phi_ix2(pGrid);

    }
#endif /* MPI_PARALLEL */

/* Physical boundaries on both left and right */
    if (pGrid->rx2_Gid < 0 && pGrid->lx2_Gid < 0) {
      (*(pD->ix2_GBCFun))(pGrid);
      (*(pD->ox2_GBCFun))(pGrid);
    }

/* shearing sheet BCs */
#ifdef SHEARING_BOX
    get_myGridIndex(pD, myID_Comm_world, &myL, &myM, &myN);
    if (myL == 0) {
      ShearingSheet_par_grav_ix1(pD);
    }
    if (myL == (pD->NGrid[0]-1)) {
      ShearingSheet_par_grav_ox1(pD);
    }
#endif

  }

/*--- Step 3. ------------------------------------------------------------------
 * Boundary Conditions in x3-direction */

  if (pGrid->Nx[2] > 1){

#ifdef MPI_PARALLEL

    cnt = (pGrid->Nx[0] + 2*nghost)*(pGrid->Nx[1] + 2*nghost)*nghost;

/* MPI blocks to both left and right */
    if (pGrid->rx3_Gid >= 0 && pGrid->lx3_Gid >= 0) {

      /* Post non-blocking receives for data from L and R Grids */
      ierr = MPI_Irecv(&(recv_buf[0][0]),cnt,MPI_DOUBLE,pGrid->lx3_Gid,LtoR_tag,
        pD->Comm_Domain, &(recv_rq[0]));
      ierr = MPI_Irecv(&(recv_buf[1][0]),cnt,MPI_DOUBLE,pGrid->rx3_Gid,RtoL_tag,
        pD->Comm_Domain, &(recv_rq[1]));

      /* pack and send data L and R */
      pack_Phi_ix3(pGrid);
      ierr = MPI_Isend(&(send_buf[0][0]),cnt,MPI_DOUBLE,pGrid->lx3_Gid,RtoL_tag,
        pD->Comm_Domain, &(send_rq[0]));

      pack_Phi_ox3(pGrid);
      ierr = MPI_Isend(&(send_buf[1][0]),cnt,MPI_DOUBLE,pGrid->rx3_Gid,LtoR_tag,
        pD->Comm_Domain, &(send_rq[1]));

      /* check non-blocking sends have completed. */
      ierr = MPI_Waitall(2, send_rq, MPI_STATUS_IGNORE);

      /* check non-blocking receives and unpack data in any order. */
      ierr = MPI_Waitany(2,recv_rq,&mIndex,MPI_STATUS_IGNORE);
      if (mIndex == 0) unpack_Phi_ix3(pGrid);
      if (mIndex == 1) unpack_Phi_ox3(pGrid);
      ierr = MPI_Waitany(2,recv_rq,&mIndex,MPI_STATUS_IGNORE);
      if (mIndex == 0) unpack_Phi_ix3(pGrid);
      if (mIndex == 1) unpack_Phi_ox3(pGrid);

    }

/* Physical boundary on left, MPI block on right */
    if (pGrid->rx3_Gid >= 0 && pGrid->lx3_Gid < 0) {

      /* Post non-blocking receive for data from R Grid */
      ierr = MPI_Irecv(&(recv_buf[1][0]),cnt,MPI_DOUBLE,pGrid->rx3_Gid,RtoL_tag,
        pD->Comm_Domain, &(recv_rq[1]));

      /* pack and send data R */
      pack_Phi_ox3(pGrid);
      ierr = MPI_Isend(&(send_buf[1][0]),cnt,MPI_DOUBLE,pGrid->rx3_Gid,LtoR_tag,
        pD->Comm_Domain, &(send_rq[1]));

      /* set physical boundary */
      (*(pD->ix3_GBCFun))(pGrid);

      /* check non-blocking send has completed. */
      ierr = MPI_Wait(&(send_rq[1]), MPI_STATUS_IGNORE);

      /* wait on non-blocking receive from R and unpack data */
      ierr = MPI_Wait(&(recv_rq[1]), MPI_STATUS_IGNORE);
      unpack_Phi_ox3(pGrid);

    }

/* MPI block on left, Physical boundary on right */
    if (pGrid->rx3_Gid < 0 && pGrid->lx3_Gid >= 0) {

      /* Post non-blocking receive for data from L grid */
      ierr = MPI_Irecv(&(recv_buf[0][0]),cnt,MPI_DOUBLE,pGrid->lx3_Gid,LtoR_tag,
        pD->Comm_Domain, &(recv_rq[0]));

      /* pack and send data L */
      pack_Phi_ix3(pGrid);
      ierr = MPI_Isend(&(send_buf[0][0]),cnt,MPI_DOUBLE,pGrid->lx3_Gid,RtoL_tag,
        pD->Comm_Domain, &(send_rq[0]));

      /* set physical boundary */
      (*(pD->ox3_GBCFun))(pGrid);

      /* check non-blocking send has completed. */
      ierr = MPI_Wait(&(send_rq[0]), MPI_STATUS_IGNORE);

      /* wait on non-blocking receive from L and unpack data */
      ierr = MPI_Wait(&(recv_rq[0]), MPI_STATUS_IGNORE);
      unpack_Phi_ix3(pGrid);
    }
#endif /* MPI_PARALLEL */

/* Physical boundaries on both left and right */
    if (pGrid->rx3_Gid < 0 && pGrid->lx3_Gid < 0) {
      (*(pD->ix3_GBCFun))(pGrid);
      (*(pD->ox3_GBCFun))(pGrid);
    }

  }

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn void bvals_particle_self_gravity_init(MeshS *pM) 
 *  \brief Sets function pointers for physical boundaries during
 *   initialization, allocates memory for send/receive buffers with MPI
 */

void bvals_particle_self_gravity_init(MeshS *pM)
{
  GridS *pG;
  DomainS *pD;
  int i,nl,nd,irefine;
#ifdef MPI_PARALLEL
  int myL,myM,myN,l,m,n,nx1t,nx2t,nx3t,size;
  int x1cnt=0, x2cnt=0, x3cnt=0; /* Number of words passed in x1/x2/x3-dir. */
#endif /* MPI_PARALLEL */

/* Cycle through all the Domains that have active Grids on this proc */
  for (nl=0; nl<(pM->NLevels); nl++){
  for (nd=0; nd<(pM->DomainsPerLevel[nl]); nd++){
  if (pM->Domain[nl][nd].Grid != NULL) {
    pD = (DomainS*)&(pM->Domain[nl][nd]);  /* ptr to Domain */
    pG = pM->Domain[nl][nd].Grid;          /* ptr to Grid */
    irefine = 1;
    for (i=1;i<=nl;i++) irefine *= 2;   /* C pow fn only takes doubles !! */
#ifdef MPI_PARALLEL
/* get (l,m,n) coordinates of Grid being updated on this processor */
    get_myGridIndex(pD, myID_Comm_world, &myL, &myM, &myN);
#endif /* MPI_PARALLEL */

/* Set function pointers for physical boundaries in x1-direction */

    if(pG->Nx[0] > 1) {

/*---- ix1 boundary ----------------------------------------------------------*/

      if(pD->ix1_GBCFun == NULL){

/* Domain boundary is in interior of root */
        if(pD->Disp[0] != 0) {
          pD->ix1_GBCFun = ProlongateLater;

/* Domain is at L-edge of root Domain */
        } else {
#ifndef PARTICLE_SELF_GRAVITY_USING_FFT_OBC
          switch(pM->BCFlag_ix1){

          case 1: /* Reflecting */
            pD->ix1_GBCFun = reflect_Phi_ix1;
          break;

          case 2: /* Outflow */
            ath_error("[bvals_grav_init]: BCFlag_ix1 = 2 not implemented\n");
          break;

          case 4: /* Periodic */
            pD->ix1_GBCFun = periodic_Phi_ix1;
#ifdef MPI_PARALLEL
            if(pG->lx1_Gid < 0 && pD->NGrid[0] > 1){
              pG->lx1_Gid = pD->GData[myN][myM][pD->NGrid[0]-1].ID_Comm_Domain;
	    }
#endif /* MPI_PARALLEL */
          break;

          case 5: /* Reflecting, B_normal!=0 */
            pD->ix1_GBCFun = reflect_Phi_ix1;
          break;

          default:
            ath_error("[bvals_grav_init]: BCFlag_ix1 = %d unknown\n",
            pM->BCFlag_ix1);
          }
#else /* PARTICLE_SELF_GRAVITY_USING_FFT_OBC */
        pD->ix1_GBCFun = obc_fft_Phi_ix1;
#endif
        }
      }

/*---- ox1 boundary ----------------------------------------------------------*/

      if(pD->ox1_GBCFun == NULL){

/* Domain boundary is in interior of root */
        if((pD->Disp[0] + pD->Nx[0])/irefine != pM->Nx[0]) {
          pD->ox1_GBCFun = ProlongateLater;

/* Domain is at R-edge of root Domain */
        } else {
#ifndef PARTICLE_SELF_GRAVITY_USING_FFT_OBC
          switch(pM->BCFlag_ox1){

          case 1: /* Reflecting */
            pD->ox1_GBCFun = reflect_Phi_ox1;
          break;

          case 2: /* Outflow */
            ath_error("[bvals_grav_init]: BCFlag_ox1 = 2 not implemented\n");
          break;

          case 4: /* Periodic */
            pD->ox1_GBCFun = periodic_Phi_ox1;
#ifdef MPI_PARALLEL
            if(pG->rx1_Gid < 0 && pD->NGrid[0] > 1){
              pG->rx1_Gid = pD->GData[myN][myM][0].ID_Comm_Domain;
            }
#endif /* MPI_PARALLEL */
          break;

          case 5: /* Reflecting, B_normal!=0 */
            pD->ox1_GBCFun = reflect_Phi_ox1;
          break;

          default:
            ath_error("[bvals_grav_init]: BCFlag_ox1 = %d unknown\n",
            pM->BCFlag_ox1);
          }
#else /* PARTICLE_SELF_GRAVITY_USING_FFT_OBC */
          pD->ox1_GBCFun = obc_fft_Phi_ox1;
#endif
        }
      }
    }

/* Set function pointers for physical boundaries in x2-direction */

    if(pG->Nx[1] > 1) {

/*---- ix2 boundary ----------------------------------------------------------*/

      if(pD->ix2_GBCFun == NULL){

/* Domain boundary is in interior of root */
        if(pD->Disp[1] != 0) {
          pD->ix2_GBCFun = ProlongateLater;

/* Domain is at L-edge of root Domain */
        } else {
#ifndef PARTICLE_SELF_GRAVITY_USING_FFT_OBC
          switch(pM->BCFlag_ix2){

          case 1: /* Reflecting */
            pD->ix2_GBCFun = reflect_Phi_ix2;
          break;

          case 2: /* Outflow */
            ath_error("[bvals_grav_init]: BCFlag_ix2 = 2 not implemented\n");
          break;

          case 4: /* Periodic */
            pD->ix2_GBCFun = periodic_Phi_ix2;
#ifdef MPI_PARALLEL
            if(pG->lx2_Gid < 0 && pD->NGrid[1] > 1){
              pG->lx2_Gid = pD->GData[myN][pD->NGrid[1]-1][myL].ID_Comm_Domain;
            }
#endif /* MPI_PARALLEL */
          break;

          case 5: /* Reflecting, B_normal!=0 */
            pD->ix2_GBCFun = reflect_Phi_ix2;
          break;

          default:
            ath_error("[bvals_grav_init]: BCFlag_ix2 = %d unknown\n",
            pM->BCFlag_ix2);
          }
#else /* PARTICLE_SELF_GRAVITY_USING_FFT_OBC */
          pD->ix2_GBCFun = obc_fft_Phi_ix2;
#endif
        }
      }

/*---- ox2 boundary ----------------------------------------------------------*/

    if(pD->ox2_GBCFun == NULL){

/* Domain boundary is in interior of root */
        if((pD->Disp[1] + pD->Nx[1])/irefine != pM->Nx[1]) {
          pD->ox2_GBCFun = ProlongateLater;

/* Domain is at R-edge of root Domain */
        } else {
#ifndef PARTICLE_SELF_GRAVITY_USING_FFT_OBC
          switch(pM->BCFlag_ox2){

          case 1: /* Reflecting */
            pD->ox2_GBCFun = reflect_Phi_ox2;
          break;

          case 2: /* Outflow */
            ath_error("[bvals_grav_init]: BCFlag_ox2 = 2 not implemented\n");
          break;

          case 4: /* Periodic */
            pD->ox2_GBCFun = periodic_Phi_ox2;
#ifdef MPI_PARALLEL
            if(pG->rx2_Gid < 0 && pD->NGrid[1] > 1){
              pG->rx2_Gid = pD->GData[myN][0][myL].ID_Comm_Domain;
            }
#endif /* MPI_PARALLEL */
          break;

          case 5: /* Reflecting, B_normal!=0 */
            pD->ox2_GBCFun = reflect_Phi_ox2;
          break;

          default:
            ath_error("[bvals_grav_init]: BCFlag_ox2 = %d unknown\n",
            pM->BCFlag_ox2);
          }
#else /* PARTICLE_SELF_GRAVITY_USING_FFT_OBC */
          pD->ox2_GBCFun = obc_fft_Phi_ox2;
#endif
        }
      }
    }

/* Set function pointers for physical boundaries in x3-direction */

    if(pG->Nx[2] > 1) {

/*---- ix3 boundary ----------------------------------------------------------*/

      if(pD->ix3_GBCFun == NULL){

/* Domain boundary is in interior of root */
        if(pD->Disp[2] != 0) {
          pD->ix3_GBCFun = ProlongateLater;

/* Domain is at L-edge of root Domain */
        } else {
#ifndef PARTICLE_SELF_GRAVITY_USING_FFT_OBC
          switch(pM->BCFlag_ix3){

          case 1: /* Reflecting */
            pD->ix3_GBCFun = reflect_Phi_ix3;
          break;

          case 2: /* Outflow */
            pD->ix3_GBCFun = obc_fft_Phi_ix3;
          break;

          case 4: /* Periodic */
            pD->ix3_GBCFun = periodic_Phi_ix3;
#ifdef MPI_PARALLEL
            if(pG->lx3_Gid < 0 && pD->NGrid[2] > 1){
              pG->lx3_Gid = pD->GData[pD->NGrid[2]-1][myM][myL].ID_Comm_Domain;
            }
#endif /* MPI_PARALLEL */
          break;

          case 5: /* Reflecting, B_normal!=0 */
            pD->ix3_GBCFun = reflect_Phi_ix3;
          break;

          default:
            ath_error("[bvals_particle_self_gravity_init]: BCFlag_ix3 = %d unknown\n",
            pM->BCFlag_ix3);
          }
#else /* PARTICLE_SELF_GRAVITY_USING_FFT_OBC */
          pD->ix3_GBCFun = obc_fft_Phi_ix3;
#endif
        }
      }

/*---- ox3 boundary ----------------------------------------------------------*/

    if(pD->ox3_GBCFun == NULL){

/* Domain boundary is in interior of root */
        if((pD->Disp[2] + pD->Nx[2])/irefine != pM->Nx[2]) {
          pD->ox3_GBCFun = ProlongateLater;

/* Domain is at R-edge of root Domain */
        } else {
#ifndef PARTICLE_SELF_GRAVITY_USING_FFT_OBC
          switch(pM->BCFlag_ox3){

          case 1: /* Reflecting */
            pD->ox3_GBCFun = reflect_Phi_ox3;
          break;

          case 2: /* Outflow */
            pD->ox3_GBCFun = obc_fft_Phi_ox3;
          break;

          case 4: /* Periodic */
            pD->ox3_GBCFun = periodic_Phi_ox3;
#ifdef MPI_PARALLEL
            if(pG->rx3_Gid < 0 && pD->NGrid[2] > 1){
              pG->rx3_Gid = pD->GData[0][myM][myL].ID_Comm_Domain;
            }
#endif /* MPI_PARALLEL */
          break;

          case 5: /* Reflecting, B_normal!=0 */
            pD->ox3_GBCFun = reflect_Phi_ox3;
          break;

          default:
            ath_error("[bvals_particle_self_gravity_init]: BCFlag_ox3 = %d unknown\n",
            pM->BCFlag_ox3);
          }
#else /* PARTICLE_SELF_GRAVITY_USING_FFT_OBC */
          pD->ox3_GBCFun = obc_fft_Phi_ox3;
#endif
        }
      }
    }

/* Figure out largest size needed for send/receive buffers with MPI ----------*/

#ifdef MPI_PARALLEL

    for (n=0; n<(pD->NGrid[2]); n++){
    for (m=0; m<(pD->NGrid[1]); m++){
      for (l=0; l<(pD->NGrid[0]); l++){

/* x1cnt is surface area of x1 faces */
        if(pD->NGrid[0] > 1){
          nx2t = pD->GData[n][m][l].Nx[1];
          if(nx2t > 1) nx2t += 1;

          nx3t = pD->GData[n][m][l].Nx[2];
          if(nx3t > 1) nx3t += 1;

          if(nx2t*nx3t > x1cnt) x1cnt = nx2t*nx3t;
        }

/* x2cnt is surface area of x2 faces */
        if(pD->NGrid[1] > 1){
          nx1t = pD->GData[n][m][l].Nx[0];
          if(nx1t > 1) nx1t += 2*nghost;

          nx3t = pD->GData[n][m][l].Nx[2];
          if(nx3t > 1) nx3t += 1;

          if(nx1t*nx3t > x2cnt) x2cnt = nx1t*nx3t;
        }

/* x3cnt is surface area of x3 faces */
        if(pD->NGrid[2] > 1){
          nx1t = pD->GData[n][m][l].Nx[0];
          if(nx1t > 1) nx1t += 2*nghost;

          nx2t = pD->GData[n][m][l].Nx[1];
          if(nx2t > 1) nx2t += 2*nghost;

          if(nx1t*nx2t > x3cnt) x3cnt = nx1t*nx2t;
        }
      }
    }}
#endif /* MPI_PARALLEL */

  }}}  /* End loop over all Domains with active Grids -----------------------*/

#ifdef MPI_PARALLEL
/* Allocate memory for send/receive buffers and MPI_Requests */

  size = x1cnt > x2cnt ? x1cnt : x2cnt;
  size = x3cnt >  size ? x3cnt : size;

  size *= nghost; /* Multiply by the third dimension */

  if (size > 0) {
    if((send_buf = (double**)calloc_2d_array(2,size,sizeof(double))) == NULL)
      ath_error("[bvals_init]: Failed to allocate send buffer\n");

    if((recv_buf = (double**)calloc_2d_array(2,size,sizeof(double))) == NULL)
      ath_error("[bvals_init]: Failed to allocate recv buffer\n");
  }

  if((recv_rq = (MPI_Request*) calloc_1d_array(2,sizeof(MPI_Request))) == NULL)
    ath_error("[bvals_init]: Failed to allocate recv MPI_Request array\n");
  if((send_rq = (MPI_Request*) calloc_1d_array(2,sizeof(MPI_Request))) == NULL)
    ath_error("[bvals_init]: Failed to allocate send MPI_Request array\n");

#endif /* MPI_PARALLEL */

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn void bvals_grav_fun(DomainS *pD, enum BCDirection dir, VGFun_t prob_bc)
 *  \brief Sets function pointers for user-defined BCs in problem 
 */

void bvals_grav_fun(DomainS *pD, enum BCDirection dir, VGFun_t prob_bc)
{
  switch(dir){
  case left_x1:
    pD->ix1_GBCFun = prob_bc;
    break;
  case right_x1:
    pD->ox1_GBCFun = prob_bc;
    break;
  case left_x2:
    pD->ix2_GBCFun = prob_bc;
    break;
  case right_x2:
    pD->ox2_GBCFun = prob_bc;
    break;
  case left_x3:
    pD->ix3_GBCFun = prob_bc;
    break;
  case right_x3:
    pD->ox3_GBCFun = prob_bc;
    break;
  default:
    ath_perr(-1,"[bvals_grav_fun]: Unknown direction = %d\n",dir);
    exit(EXIT_FAILURE);
  }
  return;
}

/*=========================== PRIVATE FUNCTIONS ==============================*/
/* Following are the functions:
 *   reflecting_???
 *   periodic_???
 *   send_???
 *   receive_???
 * where ???=[ix1,ox1,ix2,ox2,ix3,ox3]
 */

/*----------------------------------------------------------------------------*/
/*! \fn static void reflect_Phi_ix1(GridS *pGrid)
 *  \brief REFLECTING boundary conditions, Inner x1 boundary (ibc_x1=1)
 */

static void reflect_Phi_ix1(GridS *pGrid)
{
  int is = pGrid->is;
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k;

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {
        pGrid->par_Phi[k][j][is-i] = pGrid->par_Phi[k][j][is+(i-1)];
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void reflect_Phi_ox1(GridS *pGrid)
 *  \brief REFLECTING boundary conditions, Outer x1 boundary (obc_x1=1)
 */

static void reflect_Phi_ox1(GridS *pGrid)
{
  int ie = pGrid->ie;
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k;

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {
        pGrid->par_Phi[k][j][ie+i] = pGrid->par_Phi[k][j][ie-(i-1)];
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void reflect_Phi_ix2(GridS *pGrid)
 *  \brief REFLECTING boundary conditions, Inner x2 boundary (ibc_x2=1)
 */

static void reflect_Phi_ix2(GridS *pGrid)
{
  int is = pGrid->is, ie = pGrid->ie;
  int js = pGrid->js;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k;

  for (k=ks; k<=ke; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pGrid->par_Phi[k][js-j][i]    =  pGrid->par_Phi[k][js+(j-1)][i];
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void reflect_Phi_ox2(GridS *pGrid)
 *  \brief REFLECTING boundary conditions, Outer x2 boundary (obc_x2=1)
 */

static void reflect_Phi_ox2(GridS *pGrid)
{
  int is = pGrid->is, ie = pGrid->ie;
  int je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k;

  for (k=ks; k<=ke; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pGrid->par_Phi[k][je+j][i] = pGrid->par_Phi[k][je-(j-1)][i];
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void reflect_Phi_ix3(GridS *pGrid)
 *  \brief REFLECTING boundary conditions, Inner x3 boundary (ibc_x3=1)
 */

static void reflect_Phi_ix3(GridS *pGrid)
{
  int is = pGrid->is, ie = pGrid->ie;
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks;
  int i,j,k;

  for (k=1; k<=nghost; k++) {
    for (j=js-nghost; j<=je+nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pGrid->par_Phi[ks-k][j][i] = pGrid->par_Phi[ks+(k-1)][j][i];
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void reflect_Phi_ox3(GridS *pGrid)
 *  \brief REFLECTING boundary conditions, Outer x3 boundary (obc_x3=1)
 */

static void reflect_Phi_ox3(GridS *pGrid)
{
  int is = pGrid->is, ie = pGrid->ie;
  int js = pGrid->js, je = pGrid->je;
  int ke = pGrid->ke;
  int i,j,k;

  for (k=1; k<=nghost; k++) {
    for (j=js-nghost; j<=je+nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pGrid->par_Phi[ke+k][j][i] = pGrid->par_Phi[ke-(k-1)][j][i];
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void periodic_Phi_ix1(GridS *pGrid)
 *  \brief PERIODIC boundary conditions, Inner x1 boundary (ibc_x1=4)
 */

static void periodic_Phi_ix1(GridS *pGrid)
{
  int is = pGrid->is, ie = pGrid->ie;
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k;

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {
        pGrid->par_Phi[k][j][is-i] = pGrid->par_Phi[k][j][ie-(i-1)];
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void periodic_Phi_ox1(GridS *pGrid)
 *  \brief PERIODIC boundary conditions (cont), Outer x1 boundary (obc_x1=4)
 */

static void periodic_Phi_ox1(GridS *pGrid)
{
  int is = pGrid->is, ie = pGrid->ie;
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k;

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {
        pGrid->par_Phi[k][j][ie+i] = pGrid->par_Phi[k][j][is+(i-1)];
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void periodic_Phi_ix2(GridS *pGrid)
 *  \brief PERIODIC boundary conditions (cont), Inner x2 boundary (ibc_x2=4)
 */

static void periodic_Phi_ix2(GridS *pGrid)
{
  int is = pGrid->is, ie = pGrid->ie;
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k;

  for (k=ks; k<=ke; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pGrid->par_Phi[k][js-j][i] = pGrid->par_Phi[k][je-(j-1)][i];
      }
    }
  }
  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void periodic_Phi_ox2(GridS *pGrid)
 *  \brief PERIODIC boundary conditions (cont), Outer x2 boundary (obc_x2=4)
 */

static void periodic_Phi_ox2(GridS *pGrid)
{
  int is = pGrid->is, ie = pGrid->ie;
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k;

  for (k=ks; k<=ke; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pGrid->par_Phi[k][je+j][i] = pGrid->par_Phi[k][js+(j-1)][i];
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void periodic_Phi_ix3(GridS *pGrid)
 *  \brief PERIODIC boundary conditions (cont), Inner x3 boundary (ibc_x3=4)
 */

static void periodic_Phi_ix3(GridS *pGrid)
{
  int is = pGrid->is, ie = pGrid->ie;
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k;

  for (k=1; k<=nghost; k++) {
    for (j=js-nghost; j<=je+nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pGrid->par_Phi[ks-k][j][i] = pGrid->par_Phi[ke-(k-1)][j][i];
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void periodic_Phi_ox3(GridS *pGrid)
 *  \brief PERIODIC boundary conditions (cont), Outer x3 boundary (obc_x3=4)
 */

static void periodic_Phi_ox3(GridS *pGrid)
{
  int is = pGrid->is, ie = pGrid->ie;
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k;

  for (k=1; k<=nghost; k++) {
    for (j=js-nghost; j<=je+nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pGrid->par_Phi[ke+k][j][i] = pGrid->par_Phi[ks+(k-1)][j][i];
      }
    }
  }

  return;
}


/*----------------------------------------------------------------------------
 * OPEN BOUNDARY CONDITION FUNCTIONS
 *----------------------------------------------------------------------------*/

/* For linear extrapolation, just comment out the following line */
#define QUADRATIC_EXTRAPOLATION

/*----------------------------------------------------------------------------*/
/* OPEN (VACUUM) boundary conditions, Inner x1 boundary
 *
 * This BC uses a linear (quadratic) extrapolation of Phi normal to each face.
 * NOTE:  This version requires AT LEAST 2(3) active zones in each direction! */
static void obc_fft_Phi_ix1(GridS *pG)
{
  int i,is=pG->is;
  int j,js=pG->js,je=pG->je;
  int k,ks=pG->ks,ke=pG->ke;
  
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {
#ifdef QUADRATIC_EXTRAPOLATION
        pG->par_Phi[k][j][is-i] = 3.0*pG->par_Phi[k][j][is-i+1] - 3.0*pG->par_Phi[k][j][is-i+2] + pG->par_Phi[k][j][is-i+3];
#else  /* LINEAR EXTRAPOLATION */
        pG->par_Phi[k][j][is-i] = 2.0*pG->par_Phi[k][j][is-i+1] - pG->par_Phi[k][j][is-i+2];
#endif
      }
    }
  }
}

/*----------------------------------------------------------------------------*/
/* OPEN (VACUUM) boundary conditions, Outer x1 boundary
 */
static void obc_fft_Phi_ox1(GridS *pG)
{
  int i,ie=pG->ie;
  int j,js=pG->js,je=pG->je;
  int k,ks=pG->ks,ke=pG->ke;
  
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {
#ifdef QUADRATIC_EXTRAPOLATION
        pG->par_Phi[k][j][ie+i] = 3.0*pG->par_Phi[k][j][ie+i-1] - 3.0*pG->par_Phi[k][j][ie+i-2] + pG->par_Phi[k][j][ie+i-3];
#else  /* LINEAR EXTRAPOLATION */
        pG->par_Phi[k][j][ie+i] = 2.0*pG->par_Phi[k][j][ie+i-1] - pG->par_Phi[k][j][ie+i-2];
#endif
      }
    }
  }
}

/*----------------------------------------------------------------------------*/
/* OPEN (VACUUM) boundary conditions, Innter x2 boundary
 */
static void obc_fft_Phi_ix2(GridS *pG)
{
  int i,is=pG->is,ie=pG->ie;
  int j,js=pG->js;
  int k,ks=pG->ks,ke=pG->ke;
  
  for (k=ks; k<=ke; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
#ifdef QUADRATIC_EXTRAPOLATION
        pG->par_Phi[k][js-j][i] = 3.0*pG->par_Phi[k][js-j+1][i] - 3.0*pG->par_Phi[k][js-j+2][i] + pG->par_Phi[k][js-j+3][i];
#else  /* LINEAR EXTRAPOLATION */
        pG->par_Phi[k][js-j][i] = 2.0*pG->par_Phi[k][js-j+1][i] - pG->par_Phi[k][js-j+2][i];
#endif
      }
    }
  }
}

/*----------------------------------------------------------------------------*/
/* OPEN (VACUUM) boundary conditions, Outer x2 boundary
 */
static void obc_fft_Phi_ox2(GridS *pG)
{
  int i,is=pG->is,ie=pG->ie;
  int j,je=pG->je;
  int k,ks=pG->ks,ke=pG->ke;
  
  for (k=ks; k<=ke; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
#ifdef QUADRATIC_EXTRAPOLATION
        pG->par_Phi[k][je+j][i] = 3.0*pG->par_Phi[k][je+j-1][i] - 3.0*pG->par_Phi[k][je+j-2][i] + pG->par_Phi[k][je+j-3][i];
#else  /* LINEAR EXTRAPOLATION */
        pG->par_Phi[k][je+j][i] = 2.0*pG->par_Phi[k][je+j-1][i] - pG->par_Phi[k][je+j-2][i];
#endif
      }
    }
  }
}

/*----------------------------------------------------------------------------*/
/* OPEN (VACUUM) boundary conditions, Inner x3 boundary
 */
static void obc_fft_Phi_ix3(GridS *pG)
{
  int i,is=pG->is,ie=pG->ie;
  int j,js=pG->js,je=pG->je;
  int k,ks=pG->ks;
  
  for (k=1; k<=nghost; k++) {
    for (j=js-nghost; j<=je+nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
#ifdef QUADRATIC_EXTRAPOLATION
        pG->par_Phi[ks-k][j][i] = 3.0*pG->par_Phi[ks-k+1][j][i] - 3.0*pG->par_Phi[ks-k+2][j][i] + pG->par_Phi[ks-k+3][j][i];
#else  /* LINEAR EXTRAPOLATION */
        pG->par_Phi[ks-k][j][i] = 2.0*pG->par_Phi[ks-k+1][j][i] - pG->par_Phi[ks-k+2][j][i];
#endif
      }
    }
  }
}

/*----------------------------------------------------------------------------*/
/* OPEN (VACUUM) boundary conditions, Outer x3 boundary
 */
static void obc_fft_Phi_ox3(GridS *pG)
{
  int i,is=pG->is,ie=pG->ie;
  int j,js=pG->js,je=pG->je;
  int k,ke=pG->ke;
  
  for (k=1; k<=nghost; k++) {
    for (j=js-nghost; j<=je+nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
#ifdef QUADRATIC_EXTRAPOLATION
        pG->par_Phi[ke+k][j][i] = 3.0*pG->par_Phi[ke+k-1][j][i] - 3.0*pG->par_Phi[ke+k-2][j][i] + pG->par_Phi[ke+k-3][j][i];
#else  /* LINEAR EXTRAPOLATION */
        pG->par_Phi[ke+k][j][i] = 2.0*pG->par_Phi[ke+k-1][j][i] - pG->par_Phi[ke+k-2][j][i];
#endif
      }
    }
  }
}

/*----------------------------------------------------------------------------*/
/*! \fn static void ProlongateLater(GridS *pGrid)
 *  \brief PROLONGATION boundary conditions.  
 *
 * Nothing is actually done here, the
 * prolongation is actually handled in ProlongateGhostZones in main loop, so
 * this is just a NoOp Grid function.  */

static void ProlongateLater(GridS *pGrid)
{
  return;
}

#ifdef MPI_PARALLEL  /* This ifdef wraps the next 12 funs; ~550 lines */
/*----------------------------------------------------------------------------*/
/*! \fn static void pack_Phi_ix1(GridS *pG)
 *  \brief PACK boundary conditions for MPI_Isend, Inner x1 boundary */

static void pack_Phi_ix1(GridS *pG)
{
  int is = pG->is, ie = pG->ie;
  int js = pG->js, je = pG->je;
  int ks = pG->ks, ke = pG->ke;
  int i,j,k;
  double *pSnd;
  pSnd = (double*)&(send_buf[0][0]);

/* Pack only Phi into send buffer */
  for (k=ks; k<=ke; k++){
    for (j=js; j<=je; j++){
      for (i=is; i<=is+(nghost-1); i++){
        *(pSnd++) = pG->par_Phi[k][j][i];
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void pack_Phi_ox1(GridS *pG)
 *  \brief PACK boundary conditions for MPI_Isend, Outer x1 boundary */

static void pack_Phi_ox1(GridS *pG)
{
  int is = pG->is, ie = pG->ie;
  int js = pG->js, je = pG->je;
  int ks = pG->ks, ke = pG->ke;
  int i,j,k;
  double *pSnd;
  pSnd = (double*)&(send_buf[1][0]);

/* Pack only Phi into send buffer */
  for (k=ks; k<=ke; k++){
    for (j=js; j<=je; j++){
      for (i=ie-(nghost-1); i<=ie; i++){
        *(pSnd++) = pG->par_Phi[k][j][i];
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void pack_Phi_ix2(GridS *pG)
 *  \brief PACK boundary conditions for MPI_Isend, Inner x2 boundary */

static void pack_Phi_ix2(GridS *pG)
{
  int is = pG->is, ie = pG->ie;
  int js = pG->js, je = pG->je;
  int ks = pG->ks, ke = pG->ke;
  int i,j,k;
  double *pSnd;
  pSnd = (double*)&(send_buf[0][0]);

/* Pack only Phi into send buffer */
  for (k=ks; k<=ke; k++){
    for (j=js; j<=js+(nghost-1); j++){
      for (i=is-nghost; i<=ie+nghost; i++){
        *(pSnd++) = pG->par_Phi[k][j][i];
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void pack_Phi_ox2(GridS *pG)
 *  \brief PACK boundary conditions for MPI_Isend, Outer x2 boundary */

static void pack_Phi_ox2(GridS *pG)
{
  int is = pG->is, ie = pG->ie;
  int js = pG->js, je = pG->je;
  int ks = pG->ks, ke = pG->ke;
  int i,j,k;
  double *pSnd;
  pSnd = (double*)&(send_buf[1][0]);

/* Pack only Phi into send buffer */

  for (k=ks; k<=ke; k++){
    for (j=je-(nghost-1); j<=je; j++){
      for (i=is-nghost; i<=ie+nghost; i++){
        *(pSnd++) = pG->par_Phi[k][j][i];
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void pack_Phi_ix3(GridS *pG)
 *  \brief PACK boundary conditions for MPI_Isend, Inner x3 boundary */

static void pack_Phi_ix3(GridS *pG)
{
  int is = pG->is, ie = pG->ie;
  int js = pG->js, je = pG->je;
  int ks = pG->ks, ke = pG->ke;
  int i,j,k;
  double *pSnd;
  pSnd = (double*)&(send_buf[0][0]);

/* Pack only Phi into send buffer */

  for (k=ks; k<=ks+(nghost-1); k++){
    for (j=js-nghost; j<=je+nghost; j++){
      for (i=is-nghost; i<=ie+nghost; i++){
        *(pSnd++) = pG->par_Phi[k][j][i];
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void pack_Phi_ox3(GridS *pG)
 *  \brief PACK boundary conditions for MPI_Isend, Outer x3 boundary */

static void pack_Phi_ox3(GridS *pG)
{
  int is = pG->is, ie = pG->ie;
  int js = pG->js, je = pG->je;
  int ks = pG->ks, ke = pG->ke;
  int i,j,k;
  double *pSnd;
  pSnd = (double*)&(send_buf[1][0]);

/* Pack only Phi into send buffer */

  for (k=ke-(nghost-1); k<=ke; k++){
    for (j=js-nghost; j<=je+nghost; j++){
      for (i=is-nghost; i<=ie+nghost; i++){
        *(pSnd++) = pG->par_Phi[k][j][i];
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void unpack_Phi_ix1(GridS *pG)
 *  \brief UNPACK boundary conditions after MPI_Irecv, Inner x1 boundary */

static void unpack_Phi_ix1(GridS *pG)
{
  int is = pG->is;
  int js = pG->js, je = pG->je;
  int ks = pG->ks, ke = pG->ke;
  int i,j,k;
  double *pRcv;
  pRcv = (double*)&(recv_buf[0][0]);

/* Manually unpack the data from the receive buffer */

  for (k=ks; k<=ke; k++){
//    for (j=js; j<=js; j++){
//MODIFIED BY HAO GONG
    for (j=js; j<=je; j++){
      for (i=is-nghost; i<=is-1; i++){
        pG->par_Phi[k][j][i] = *(pRcv++);
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void unpack_Phi_ox1(GridS *pG)
 *  \brief UNPACK boundary conditions after MPI_Irecv, Outer x1 boundary */

static void unpack_Phi_ox1(GridS *pG)
{
  int ie = pG->ie;
  int js = pG->js, je = pG->je;
  int ks = pG->ks, ke = pG->ke;
  int i,j,k;
  double *pRcv;
  pRcv = (double*)&(recv_buf[1][0]);

/* Manually unpack the data from the receive buffer */

  for (k=ks; k<=ke; k++){
    for (j=js; j<=je; j++){
      for (i=ie+1; i<=ie+nghost; i++){
        pG->par_Phi[k][j][i] = *(pRcv++);
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void unpack_Phi_ix2(GridS *pG)
 *  \brief UNPACK boundary conditions after MPI_Irecv, Inner x2 boundary */

static void unpack_Phi_ix2(GridS *pG)
{
  int is = pG->is, ie = pG->ie;
  int js = pG->js;
  int ks = pG->ks, ke = pG->ke;
  int i,j,k;
  double *pRcv;
  pRcv = (double*)&(recv_buf[0][0]);

/* Manually unpack the data from the receive buffer */

  for (k=ks; k<=ke; k++){
    for (j=js-nghost; j<=js-1; j++){
      for (i=is-nghost; i<=ie+nghost; i++){
        pG->par_Phi[k][j][i] = *(pRcv++);
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void unpack_Phi_ox2(GridS *pG)
 *  \brief UNPACK boundary conditions after MPI_Irecv, Outer x2 boundary */

static void unpack_Phi_ox2(GridS *pG)
{
  int is = pG->is, ie = pG->ie;
  int je = pG->je;
  int ks = pG->ks, ke = pG->ke;
  int i,j,k;
  double *pRcv;
  pRcv = (double*)&(recv_buf[1][0]);

/* Manually unpack the data from the receive buffer */

  for (k=ks; k<=ke; k++){
    for (j=je+1; j<=je+nghost; j++){
      for (i=is-nghost; i<=ie+nghost; i++){
        pG->par_Phi[k][j][i] = *(pRcv++);
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void unpack_Phi_ix3(GridS *pG)
 *  \brief UNPACK boundary conditions after MPI_Irecv, Inner x3 boundary */

static void unpack_Phi_ix3(GridS *pG)
{
  int is = pG->is, ie = pG->ie;
  int js = pG->js, je = pG->je;
  int ks = pG->ks;
  int i,j,k;
  double *pRcv;
  pRcv = (double*)&(recv_buf[0][0]);

/* Manually unpack the data from the receive buffer */

  for (k=ks-nghost; k<=ks-1; k++){
    for (j=js-nghost; j<=je+nghost; j++){
      for (i=is-nghost; i<=ie+nghost; i++){
        pG->par_Phi[k][j][i] = *(pRcv++);
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void unpack_Phi_ox3(GridS *pG)
 *  \brief UNPACK boundary conditions after MPI_Irecv, Outer x3 boundary */

static void unpack_Phi_ox3(GridS *pG)
{
  int is = pG->is, ie = pG->ie;
  int js = pG->js, je = pG->je;
  int ke = pG->ke;
  int i,j,k;
  double *pRcv;
  pRcv = (double*)&(recv_buf[1][0]);

/* Manually unpack the data from the receive buffer */

  for (k=ke+1; k<=ke+nghost; k++){
    for (j=js-nghost; j<=je+nghost; j++){
      for (i=is-nghost; i<=ie+nghost; i++){
        pG->par_Phi[k][j][i] = *(pRcv++);
      }
    }
  }

  return;
}

#endif /* MPI_PARALLEL */

#ifdef SHEARING_BOX
/*----------------------------------------------------------------------------*/
/*! \fn void ShearingSheet_par_grav_ix1(DomainS *pD)
 *  \brief 3D shearing-sheet BCs for gravitational potential in x1.  
 *
 * It applies a remap
 * in Y after the ghost cells have been set by the usual periodic BCs in X and
 * Y */
/*----------------------------------------------------------------------------*/

void ShearingSheet_par_grav_ix1(DomainS *pD)
{
  GridS *pG = pD->Grid;
  int is = pG->is, ie = pG->ie;
  int js = pG->js, je = pG->je;
  int ks = pG->ks, ke = pG->ke;
  int i,ii,j,k,ku,n,joffset,jremap;
  Real xmin,xmax,Lx,Ly,qomL,yshear,deltay,epsi;
#ifdef MPI_PARALLEL
  int my_iproc,my_jproc,my_kproc,cnt,jproc,joverlap,Ngrids;
  int ierr,sendto_id,getfrom_id;
  double *pSnd,*pRcv;
  Remap *pRemap;
  MPI_Request rq;
#endif

/*--- Step 1. ------------------------------------------------------------------
 * Compute the distance the computational domain has sheared in y */

  xmin = pD->RootMinX[0];
  xmax = pD->RootMaxX[0];
  Lx = xmax - xmin;

  xmin = pD->RootMinX[1];
  xmax = pD->RootMaxX[1];
  Ly = xmax - xmin;

  qomL = qshear*Omega_0*Lx;
  yshear = qomL*(pG->time);

/* Split this into integer and fractional peices of the Domain in y.  Ignore
 * the integer piece because the Grid is periodic in y */

  deltay = fmod(yshear, Ly);

/* further decompose the fractional peice into integer and fractional pieces of
 * a grid cell.  Note 0.0 <= epsi < 1.0.  If Domain has MPI decomposition in Y,
 * then it is possible that:  pD->Nx2 > joffset > pG->Nx2   */

  joffset = (int)(deltay/pG->dx2);
  epsi = (fmod(deltay,pG->dx2))/pG->dx2;

/*--- Step 2. ------------------------------------------------------------------
 * Copy data into GhstZns array.  Note i and j indices are switched.
 * Steps 2-10 are for 3D or 2D xy runs.  Step 11 handles 2D xz separately */

  if (pG->Nx[2] > 1 || ShBoxCoord==xy) {  /* this if ends at end of step 10 */

  if (pG->Nx[2] > 1) ku=ke+1; else ku=ke;
  for(k=ks; k<=ku; k++) {
    for(j=js-nghost; j<=je+nghost; j++){
      for(i=0; i<nghost; i++){
        ii = is-nghost+i;
        GhstZns[k][i][j].par_Phi = pG->par_Phi[k][j][ii];
      }
    }
  }

/*--- Step 3. ------------------------------------------------------------------
 * Copy GhstZns into buffer, at the same time apply a conservative remap of
 * solution over the fractional part of grid cell */

  for(k=ks; k<=ku; k++) {
    for(i=0; i<nghost; i++){

      for (j=js-nghost; j<=je+nghost; j++) U[j] = GhstZns[k][i][j].par_Phi;
      RemapFluxPSG(U,epsi,js,je+1,Flx);
      for(j=js; j<=je; j++){
        GhstZnsBuf[k][i][j].par_Phi = GhstZns[k][i][j].par_Phi - (Flx[j+1]-Flx[j]);
      }
    }
  }

/*--- Step 4. ------------------------------------------------------------------
 * If no MPI decomposition in Y, apply shift over integer number of
 * grid cells during copy from buffer back into GhstZns.  */

  if (pD->NGrid[1] == 1) {

    for(k=ks; k<=ku; k++) {
      for(j=js; j<=je; j++){
        jremap = j - joffset;
        if (jremap < (int)js) jremap += pG->Nx[1];

        for(i=0; i<nghost; i++){
          GhstZns[k][i][j].par_Phi = GhstZnsBuf[k][i][jremap].par_Phi;
        }
      }
    }

#ifdef MPI_PARALLEL
  } else {

/*--- Step 5. ------------------------------------------------------------------
 * If Domain contains MPI decomposition in Y, then MPI calls are required for
 * the cyclic shift needed to apply shift over integer number of grid cells
 * during copy from buffer back into GhstZns.  */

    get_myGridIndex(pD, myID_Comm_world, &my_iproc, &my_jproc, &my_kproc);

/* Find integer and fractional number of grids over which offset extends.
 * This assumes every grid has same number of cells in x2-direction! */
    Ngrids = (int)(joffset/pG->Nx[1]);
    joverlap = joffset - Ngrids*pG->Nx[1];

/*--- Step 5a. -----------------------------------------------------------------
 * Find ids of processors that data in [je-(joverlap-1):je] is sent to, and
 * data in [js:js+(joverlap-1)] is received from.  Only execute if joverlap>0  */

    if (joverlap != 0) {

      jproc = my_jproc + (Ngrids + 1);
      if (jproc > (pD->NGrid[1]-1)) jproc -= pD->NGrid[1]; 
      sendto_id = pD->GData[my_kproc][jproc][my_iproc].ID_Comm_Domain;

      jproc = my_jproc - (Ngrids + 1);
      if (jproc < 0) jproc += pD->NGrid[1]; 
      getfrom_id = pD->GData[my_kproc][jproc][my_iproc].ID_Comm_Domain;

/*--- Step 5b. -----------------------------------------------------------------
 * Pack send buffer and send data in [je-(joverlap-1):je] from GhstZnsBuf */

      cnt = nghost*joverlap*(ku-ks+1);
/* Post a non-blocking receive for the input data */
      ierr = MPI_Irecv(recv_bufss, cnt, MPI_DOUBLE, getfrom_id,
                      shearing_sheet_grav_ix1_tag, pD->Comm_Domain, &rq);

      pSnd = send_bufss;
      for (k=ks; k<=ku; k++) {
        for (j=je-(joverlap-1); j<=je; j++) {
          for(i=0; i<nghost; i++){
            /* Get a pointer to the Remap structure */
            pRemap = &(GhstZnsBuf[k][i][j]);
            *(pSnd++) = pRemap->par_Phi;
          }
        }
      }
      ierr = MPI_Send(send_bufss, cnt, MPI_DOUBLE, sendto_id,
                   shearing_sheet_grav_ix1_tag, pD->Comm_Domain);

/*--- Step 5c. -----------------------------------------------------------------
 * unpack data sent from [je-(joverlap-1):je], and remap into cells in
 * [js:js+(joverlap-1)] in GhstZns */

      ierr = MPI_Wait(&rq, MPI_STATUS_IGNORE);

      pRcv = recv_bufss;
      for (k=ks; k<=ku; k++) {
        for (j=js; j<=js+(joverlap-1); j++) {
          for(i=0; i<nghost; i++){
            pRemap = &(GhstZns[k][i][j]);
            pRemap->par_Phi = *(pRcv++);
          }
        }
      }

    }

/*--- Step 5d. -----------------------------------------------------------------
 * If shear is less one full Grid, remap cells which remain on same processor
 * from GhstZnsBuf into GhstZns.  Cells in [js:je-joverlap] are shifted by
 * joverlap into [js+joverlap:je] */

    if (Ngrids == 0) {

      for(k=ks; k<=ku; k++) {
        for(j=js+joverlap; j<=je; j++){
          jremap = j-joverlap;
          for(i=0; i<nghost; i++){
            GhstZns[k][i][j].par_Phi = GhstZnsBuf[k][i][jremap].par_Phi;
          }
        }
      }

/*--- Step 5e. -----------------------------------------------------------------
 * If shear is more than one Grid, pack and send data from [js:je-joverlap]
 * from GhstZnsBuf (this step replaces 5d) */

    } else {

/* index of sendto and getfrom processors in GridArray are -/+1 from Step 5a */

      jproc = my_jproc + Ngrids;
      if (jproc > (pD->NGrid[1]-1)) jproc -= pD->NGrid[1];
      sendto_id = pD->GData[my_kproc][jproc][my_iproc].ID_Comm_Domain;

      jproc = my_jproc - Ngrids;
      if (jproc < 0) jproc += pD->NGrid[1];
      getfrom_id = pD->GData[my_kproc][jproc][my_iproc].ID_Comm_Domain;

      cnt = nghost*(pG->Nx[1]-joverlap)*(ku-ks+1);
/* Post a non-blocking receive for the input data from the left grid */
      ierr = MPI_Irecv(recv_bufss, cnt, MPI_DOUBLE, getfrom_id,
                      shearing_sheet_grav_ix1_tag, pD->Comm_Domain, &rq);

      pSnd = send_bufss;
      for (k=ks; k<=ku; k++) {
        for (j=js; j<=je-joverlap; j++) {
          for(i=0; i<nghost; i++){
            /* Get a pointer to the Remap structure */
            pRemap = &(GhstZnsBuf[k][i][j]);
            *(pSnd++) = pRemap->par_Phi;
          }
        }
      }
      ierr = MPI_Send(send_bufss, cnt, MPI_DOUBLE, sendto_id,
                     shearing_sheet_grav_ix1_tag, pD->Comm_Domain);

/* unpack data sent from [js:je-overlap], and remap into cells in
 * [js+joverlap:je] in GhstZns */
      ierr = MPI_Wait(&rq, MPI_STATUS_IGNORE);

      pRcv = recv_bufss;
      for (k=ks; k<=ku; k++) {
        for (j=js+joverlap; j<=je; j++) {
          for(i=0; i<nghost; i++){
            /* Get a pointer to the Remap structure */
            pRemap = &(GhstZns[k][i][j]);
            pRemap->par_Phi = *(pRcv++);
          }
        }
      }
    } /* end of step 5e - shear is more than one Grid */

#endif /* MPI_PARALLEL */
  } /* end of step 5 - MPI decomposition in Y */

/*--- Step 6. ------------------------------------------------------------------
 * Now copy remapped variables back into ghost cells */

  for(k=ks; k<=ke; k++) {
    for(j=js; j<=je; j++){
      for(i=0; i<nghost; i++){
        pG->par_Phi[k][j][is-nghost+i] = GhstZns[k][i][j].par_Phi;
      }
    }
  }

/*--- Step 8. ------------------------------------------------------------------
 * With no MPI decomposition in Y, apply periodic BCs in Y (similar to
 * periodic_ix2() and periodic_ox2() in set_bavls_mhd.c) */

  if (pD->NGrid[1] == 1) {

    for(k=ks; k<=ke; k++) {
      for(j=1; j<=nghost; j++){
        for(i=is-nghost; i<is; i++){
          pG->par_Phi[k][js-j][i] = pG->par_Phi[k][je-(j-1)][i];
          pG->par_Phi[k][je+j][i] = pG->par_Phi[k][js+(j-1)][i];
        }
      }
    }

#ifdef MPI_PARALLEL
  } else {

/*--- Step 9. ------------------------------------------------------------------
 * With MPI decomposition in Y, use MPI calls to handle periodic BCs in Y (like
 * send_ox2/receive_ix1 and send_ix1/receive_ox2 pairs in set_bvals_mhd.c */


/* Post a non-blocking receive for the input data from the left grid */
    cnt = nghost*nghost*(ku-ks+1);
    ierr = MPI_Irecv(recv_bufss, cnt, MPI_DOUBLE, pG->lx2_id,
                    shearing_sheet_grav_ix1_tag, pD->Comm_Domain, &rq);

    pSnd = send_bufss;
    for (k=ks; k<=ku; k++){
      for (j=je-nghost+1; j<=je; j++){
        for (i=is-nghost; i<is; i++){
          *(pSnd++) = pG->par_Phi[k][j][i];
        }
      }
    }

/* send contents of buffer to the neighboring grid on R-x2 */
    ierr = MPI_Send(send_bufss, cnt, MPI_DOUBLE, pG->rx2_id,
                   shearing_sheet_grav_ix1_tag, pD->Comm_Domain);

/* Wait to receive the input data from the left grid */
    ierr = MPI_Wait(&rq, MPI_STATUS_IGNORE);

    pRcv = recv_bufss;
    for (k=ks; k<=ku; k++){
      for (j=js-nghost; j<=js-1; j++){
        for (i=is-nghost; i<is; i++){
          pG->par_Phi[k][j][i] = *(pRcv++);
        }
      }
    }

/* Post a non-blocking receive for the input data from the right grid */
    ierr = MPI_Irecv(recv_bufss, cnt, MPI_DOUBLE, pG->rx2_id,
                    shearing_sheet_grav_ix1_tag, pD->Comm_Domain, &rq);

    pSnd = send_bufss;
    for (k=ks; k<=ku; k++){
      for (j=js; j<=js+nghost-1; j++){
        for (i=is-nghost; i<is; i++){
          *(pSnd++) = pG->par_Phi[k][j][i];
        }
      }
    }

/* send contents of buffer to the neighboring grid on L-x2 */
    ierr = MPI_Send(send_bufss, cnt, MPI_DOUBLE, pG->lx2_id,
                   shearing_sheet_grav_ix1_tag, pD->Comm_Domain);

/* Wait to receive the input data from the left grid */
    ierr = MPI_Wait(&rq, MPI_STATUS_IGNORE);

    pRcv = recv_bufss;
    for (k=ks; k<=ku; k++){
      for (j=je+1; j<=je+nghost; j++){
        for (i=is-nghost; i<is; i++){
          pG->par_Phi[k][j][i] = *(pRcv++);
        }
      }
    }
#endif /* MPI_PARALLEL */

  } /* end of step 9 - periodic BC in Y with MPI */

  } /* end of if */

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn void ShearingSheet_par_grav_ox1(DomainS *pD)
 *  \brief 3D shearing-sheet BCs for gravitational potential in x1.  
 *
 * It applies a remap
 * in Y after the ghost cells have been set by the usual periodic BCs in X and
 * Y */
/*----------------------------------------------------------------------------*/

void ShearingSheet_par_grav_ox1(DomainS *pD)
{
  GridS *pG = pD->Grid;
  int is = pG->is, ie = pG->ie;
  int js = pG->js, je = pG->je;
  int ks = pG->ks, ke = pG->ke;
  int i,ii,j,k,ku,n,joffset,jremap;
  Real xmin,xmax,Lx,Ly,qomL,yshear,deltay,epso;
#ifdef MPI_PARALLEL
  int my_iproc,my_jproc,my_kproc,cnt,jproc,joverlap,Ngrids;
  int ierr,sendto_id,getfrom_id;
  double *pSnd,*pRcv;
  Remap *pRemap;
  MPI_Request rq;
#endif

/*--- Step 1. ------------------------------------------------------------------
 * Compute the distance the computational domain has sheared in y */

  xmin = pD->RootMinX[0];
  xmax = pD->RootMaxX[0];
  Lx = xmax - xmin;

  xmin = pD->RootMinX[1];
  xmax = pD->RootMaxX[1];
  Ly = xmax - xmin;

  qomL = qshear*Omega_0*Lx;
  yshear = qomL*(pG->time);

/* Split this into integer and fractional peices of the Domain in y.  Ignore
 * the integer piece because the Grid is periodic in y */

  deltay = fmod(yshear, Ly);

/* further decompose the fractional peice into integer and fractional pieces of
 * a grid cell.  Note 0.0 <= epsi < 1.0.  If Domain has MPI decomposition in Y,
 * then it is possible that:  pD->Nx2 > joffset > pG->Nx2   */

  joffset = (int)(deltay/pG->dx2);
  epso = -(fmod(deltay,pG->dx2))/pG->dx2;

/*--- Step 2. ------------------------------------------------------------------
 * Copy data into GhstZns array.  Note i and j indices are switched.
 * Steps 2-10 are for 3D runs.  Step 11 handles 2D separately */

  if (pG->Nx[2] > 1 || ShBoxCoord==xy) {  /* this if ends at end of step 10 */

  if (pG->Nx[2] > 1) ku=ke+1; else ku=ke;

  for(k=ks; k<=ku; k++) {
    for(j=js-nghost; j<=je+nghost; j++){
      for(i=0; i<nghost; i++){
        ii = ie+1+i;
        GhstZns[k][i][j].par_Phi = pG->par_Phi[k][j][ii];
      }
    }
  }

/*--- Step 3. ------------------------------------------------------------------
 * Copy GhstZns into buffer, at the same time apply a conservative remap of
 * solution over the fractional part of grid cell */

  for(k=ks; k<=ku; k++) {
    for(i=0; i<nghost; i++){
      for (j=js-nghost; j<=je+nghost; j++) U[j] = GhstZns[k][i][j].par_Phi;
        RemapFluxPSG(U,epso,js,je+1,Flx);
      for(j=js; j<=je; j++){
        GhstZnsBuf[k][i][j].par_Phi = GhstZns[k][i][j].par_Phi - (Flx[j+1]-Flx[j]);
      }
    }
  }

/*--- Step 4. ------------------------------------------------------------------
 * If no MPI decomposition in Y, apply shift over integer number of
 * grid cells during copy from buffer back into GhstZns.  */

  if (pD->NGrid[1] == 1) {

    for(k=ks; k<=ku; k++) {
      for(j=js; j<=je; j++){
        jremap = j + joffset;
        if (jremap > (int)je) jremap -= pG->Nx[1];

        for(i=0; i<nghost; i++){
          GhstZns[k][i][j].par_Phi = GhstZnsBuf[k][i][jremap].par_Phi;
        }

      }
    }

#ifdef MPI_PARALLEL
  } else {

/*--- Step 5. ------------------------------------------------------------------
 * If Domain contains MPI decomposition in Y, then MPI calls are required for
 * the cyclic shift needed to apply shift over integer number of grid cells
 * during copy from buffer back into GhstZns.  */

    get_myGridIndex(pD, myID_Comm_world, &my_iproc, &my_jproc, &my_kproc);

/* Find integer and fractional number of grids over which offset extends.
 * This assumes every grid has same number of cells in x2-direction! */
    Ngrids = (int)(joffset/pG->Nx[1]);
    joverlap = joffset - Ngrids*pG->Nx[1];

/*--- Step 5a. -----------------------------------------------------------------
 * Find ids of processors that data in [js:js+(joverlap-1)] is sent to, and
 * data in [je-(overlap-1):je] is received from.  Only execute if joverlap>0  */

    if (joverlap != 0) {

      jproc = my_jproc - (Ngrids + 1);
      if (jproc < 0) jproc += pD->NGrid[1]; 
      sendto_id = pD->GData[my_kproc][jproc][my_iproc].ID_Comm_Domain;

      jproc = my_jproc + (Ngrids + 1);
      if (jproc > (pD->NGrid[1]-1)) jproc -= pD->NGrid[1]; 
      getfrom_id = pD->GData[my_kproc][jproc][my_iproc].ID_Comm_Domain;

/*--- Step 5b. -----------------------------------------------------------------
 * Pack send buffer and send data in [js:js+(joverlap-1)] from GhstZnsBuf */

      cnt = nghost*joverlap*(ku-ks+1);
/* Post a non-blocking receive for the input data */
      ierr = MPI_Irecv(recv_bufss, cnt, MPI_DOUBLE, getfrom_id,
                      shearing_sheet_grav_ox1_tag, pD->Comm_Domain, &rq);

      pSnd = send_bufss;
      for (k=ks; k<=ku; k++) {
        for (j=js; j<=js+(joverlap-1); j++) {
          for(i=0; i<nghost; i++){
            /* Get a pointer to the Remap structure */
            pRemap = &(GhstZnsBuf[k][i][j]);
            *(pSnd++) = pRemap->par_Phi;
          }
        }
      }
      ierr = MPI_Send(send_bufss, cnt, MPI_DOUBLE, sendto_id,
                   shearing_sheet_grav_ox1_tag, pD->Comm_Domain);


/*--- Step 5c. -----------------------------------------------------------------
 * unpack data sent from [js:js+(joverlap-1)], and remap into cells in
 * [je-(joverlap-1):je] in GhstZns
 */

      ierr = MPI_Wait(&rq, MPI_STATUS_IGNORE);

      pRcv = recv_bufss;
      for (k=ks; k<=ku; k++) {
        for (j=je-(joverlap-1); j<=je; j++) {
          for(i=0; i<nghost; i++){
            /* Get a pointer to the Remap structure */
            pRemap = &(GhstZns[k][i][j]);
            pRemap->par_Phi = *(pRcv++);
          }
        }
      }

    }

/*--- Step 5d. -----------------------------------------------------------------
 * If shear is less one full Grid, remap cells which remain on same processor
 * from GhstZnsBuf into GhstZns.  Cells in [js+joverlap:je] are shifted by
 * joverlap into [js:je-joverlap] */

    if (Ngrids == 0) {

      for(k=ks; k<=ku; k++) {
        for(j=js; j<=je-joverlap; j++){
          jremap = j+joverlap;
          for(i=0; i<nghost; i++){
            GhstZns[k][i][j].par_Phi = GhstZnsBuf[k][i][jremap].par_Phi;
          }
        }
      }

/*--- Step 5e. -----------------------------------------------------------------
 * If shear is more than one Grid, pack and send data from [js+joverlap:je]
 * from GhstZnsBuf (this step replaces 5d) */

    } else {

/* index of sendto and getfrom processors in GridArray are -/+1 from Step 5a */

      jproc = my_jproc - Ngrids;
      if (jproc < 0) jproc += pD->NGrid[1];
      sendto_id = pD->GData[my_kproc][jproc][my_iproc].ID_Comm_Domain;

      jproc = my_jproc + Ngrids;
      if (jproc > (pD->NGrid[1]-1)) jproc -= pD->NGrid[1];
      getfrom_id = pD->GData[my_kproc][jproc][my_iproc].ID_Comm_Domain;

      cnt = nghost*(pG->Nx[1]-joverlap)*(ku-ks+1);
/* Post a non-blocking receive for the input data from the left grid */
      ierr = MPI_Irecv(recv_bufss, cnt, MPI_DOUBLE, getfrom_id,
                      shearing_sheet_grav_ox1_tag, pD->Comm_Domain, &rq);

      pSnd = send_bufss;
      for (k=ks; k<=ku; k++) {
        for (j=js+joverlap; j<=je; j++) {
          for(i=0; i<nghost; i++){
            /* Get a pointer to the Remap structure */
            pRemap = &(GhstZnsBuf[k][i][j]);
            *(pSnd++) = pRemap->par_Phi;
          }
        }
      }
      ierr = MPI_Send(send_bufss, cnt, MPI_DOUBLE, sendto_id,
                   shearing_sheet_grav_ox1_tag, pD->Comm_Domain);

/* unpack data sent from [js+joverlap:je], and remap into cells in
 * [js:je-joverlap] in GhstZns */

      ierr = MPI_Wait(&rq, MPI_STATUS_IGNORE);

      pRcv = recv_bufss;
      for (k=ks; k<=ku; k++) {
        for (j=js; j<=je-joverlap; j++) {
          for(i=0; i<nghost; i++){
            /* Get a pointer to the Remap structure */
            pRemap = &(GhstZns[k][i][j]);
            pRemap->par_Phi = *(pRcv++);
          }
        }
      }
    } /* end of step 5e - shear is more than one Grid */

#endif /* MPI_PARALLEL */
  } /* end of step 5 - MPI decomposition in Y */

/*--- Step 6. ------------------------------------------------------------------
 * Now copy remapped variables back into ghost cells */

  for(k=ks; k<=ke; k++) {
    for(j=js; j<=je; j++){
      for(i=0; i<nghost; i++){
        pG->par_Phi[k][j][ie+1+i] = GhstZns[k][i][j].par_Phi;
      }
    }
  }

/*--- Step 8. ------------------------------------------------------------------
 * With no MPI decomposition in Y, apply periodic BCs in Y (similar to
 * periodic_ix2() and periodic_ox2() in set_bavls_mhd.c) */

  if (pD->NGrid[1] == 1) {

    for(k=ks; k<=ke; k++) {
      for(j=1; j<=nghost; j++){
        for(i=ie+1; i<=ie+nghost; i++){
          pG->par_Phi[k][js-j][i] = pG->par_Phi[k][je-(j-1)][i];
          pG->par_Phi[k][je+j][i] = pG->par_Phi[k][js+(j-1)][i];
        }
      }
    }
#ifdef MPI_PARALLEL
  } else {

/*--- Step 9. ------------------------------------------------------------------
 * With MPI decomposition in Y, use MPI calls to handle periodic BCs in Y (like
 * send_ox2/receive_ix1 and send_ix1/receive_ox2 pairs in set_bvals_mhd.c */


/* Post a non-blocking receive for the input data from the left grid */
    cnt = nghost*nghost*(ku-ks + 1);
    ierr = MPI_Irecv(recv_bufss, cnt, MPI_DOUBLE, pG->lx2_id,
                    shearing_sheet_grav_ox1_tag, pD->Comm_Domain, &rq);

    pSnd = send_bufss;
    for (k=ks; k<=ku; k++){
      for (j=je-nghost+1; j<=je; j++){
        for (i=ie+1; i<=ie+nghost; i++){
          *(pSnd++) = pG->par_Phi[k][j][i];
        }
      }
    }

/* send contents of buffer to the neighboring grid on R-x2 */
    ierr = MPI_Send(send_bufss, cnt, MPI_DOUBLE, pG->rx2_id,
                   shearing_sheet_grav_ox1_tag, pD->Comm_Domain);

/* Wait to receive the input data from the left grid */
    ierr = MPI_Wait(&rq, MPI_STATUS_IGNORE);

    pRcv = recv_bufss;
    for (k=ks; k<=ku; k++){
      for (j=js-nghost; j<=js-1; j++){
        for (i=ie+1; i<=ie+nghost; i++){
          pG->par_Phi[k][j][i] = *(pRcv++);
        }
      }
    }

/* Post a non-blocking receive for the input data from the right grid */
    ierr = MPI_Irecv(recv_bufss, cnt, MPI_DOUBLE, pG->rx2_id,
                    shearing_sheet_grav_ox1_tag, pD->Comm_Domain, &rq);

    pSnd = send_bufss;
    for (k=ks; k<=ku; k++){
      for (j=js; j<=js+nghost-1; j++){
        for (i=ie+1; i<=ie+nghost; i++){
          *(pSnd++) = pG->par_Phi[k][j][i];
        }
      }
    }

/* send contents of buffer to the neighboring grid on L-x2 */
    ierr = MPI_Send(send_bufss, cnt, MPI_DOUBLE, pG->lx2_id,
                   shearing_sheet_grav_ox1_tag, pD->Comm_Domain);

/* Wait to receive the input data from the left grid */
    ierr = MPI_Wait(&rq, MPI_STATUS_IGNORE);

    pRcv = recv_bufss;
    for (k=ks; k<=ku; k++){
      for (j=je+1; j<=je+nghost; j++){
        for (i=ie+1; i<=ie+nghost; i++){
          pG->par_Phi[k][j][i] = *(pRcv++);
        }
      }
    }
#endif /* MPI_PARALLEL */

  } /* end of step 9 - periodic BC in Y with MPI */

  } /* end of if */

  return;
}

#endif /* Shearing Box */

#if defined(SHEARING_BOX) 
/*----------------------------------------------------------------------------*/
/*! \fn void bvals_shear_init(MeshS *pM)
 *  \brief Allocates memory for temporary arrays/buffers
 */

void bvals_particle_self_gravity_shear_init(MeshS *pM)
{
  GridS *pG;
  int nl,nd,nx1,nx2,nx3,max1=0,max2=0,max3=0;
#ifdef MPI_PARALLEL
  int size1=0,size2=0,size;
#endif

/* Loop over all Grids on this processor to find maximum size of arrays */

  for (nl=0; nl<(pM->NLevels); nl++){
    for (nd=0; nd<(pM->DomainsPerLevel[nl]); nd++){
      if (pM->Domain[nl][nd].Grid != NULL) { /* there is a Grid on this proc */
        pG=pM->Domain[nl][nd].Grid;          /* set pointer to Grid */

        nx1 = pG->Nx[0] + 2*nghost;
        nx2 = pG->Nx[1] + 2*nghost;
        nx3 = pG->Nx[2] + 2*nghost;
        max1 = MAX(max1,nx1);
        max2 = MAX(max2,nx2);
        max3 = MAX(max3,nx3);
      }
    }
  }

/* Allocate memory for temporary arrays and vectors */

  if((GhstZns=(Remap***)calloc_3d_array(max3,nghost,max2,sizeof(Remap)))==NULL)
    ath_error("[bvals_particle_self_gravity_shear_init]: malloc returned a NULL pointer\n");

  if((GhstZnsBuf=(Remap***)calloc_3d_array(max3,nghost,max2,sizeof(Remap))) ==
    NULL) ath_error("[bvals_particle_self_gravity_shear_init]: malloc returned a NULL pointer\n");

  if((U = (Real*)malloc(max2*sizeof(Real))) == NULL)
    ath_error("[bvals_particle_self_gravity_shear_init]: malloc returned a NULL pointer\n");

  if((Flx = (Real*)malloc(max2*sizeof(Real))) == NULL)
    ath_error("[bvals_particle_self_gravity_shear_init]: malloc returned a NULL pointer\n");

#if defined(THIRD_ORDER_CHAR) || defined(THIRD_ORDER_PRIM)
  if ((Uhalf = (Real*)malloc(max2*sizeof(Real))) == NULL)
    ath_error("[bvals_particle_self_gravity_shear_init]: malloc returned a NULL pointer\n");
#endif

/* allocate memory for send/receive buffers in MPI parallel calculations */

#ifdef MPI_PARALLEL
  size1 = nghost*pG->Nx[1]*(pG->Nx[2]+1);
  size = MAX(size1,size2);

  if((send_bufss = (double*)malloc(size*sizeof(double))) == NULL)
    ath_error("[bvals_particle_self_gravity_shear_init]: Failed to allocate send buffer\n");

  if((recv_bufss = (double*)malloc(size*sizeof(double))) == NULL)
    ath_error("[bvals_particle_self_gravity_shear_init]: Failed to allocate receive buffer\n");
#endif /* MPI_PARALLEL */

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn void bvals_shear_destruct(void)
 *  \brief Free temporary arrays
 */

void bvals_particle_self_gravity_shear_destruct(void)
{
  if (GhstZns    != NULL) free_3d_array(GhstZns);
  if (GhstZnsBuf != NULL) free_3d_array(GhstZnsBuf);
#if defined(THIRD_ORDER_CHAR) || defined(THIRD_ORDER_PRIM)
  if (Uhalf != NULL) free(Uhalf);
#endif
  if (U   != NULL) free(U);
  if (Flx != NULL) free(Flx);
#ifdef MPI_PARALLEL
  if (send_buf != NULL) free(send_buf);
  if (recv_buf != NULL) free(recv_buf);
  if (send_bufss != NULL) free(send_bufss);
  if (recv_bufss != NULL) free(recv_bufss);
#endif /* MPI_PARALLEL */

  return;
}

/*------------------------------------------------------------------------------
 * RemapFluxPSG: computes "fluxes" of particle potential for conservative remap
 * Input Arguments:
 *   U = 1D vector of conserved variable at cell centers along 1-D slice
 *   eps = fraction of a cell to be remapped
 * Output Arguments:
 *   Flux = fluxes of conserved variable at interfaces over [jinner:jouter]
 */

#if defined(SECOND_ORDER_CHAR) || defined (SECOND_ORDER_PRIM)
/*----------------------------------------------------------------------------*/
/*! \fn void RemapFluxPSG(const Real *U, const Real eps,
 *             const int jinner, const int jouter, Real *Flux)
 *  \brief Second order reconstruction for conservative remap.
 *   using piecewise linear reconstruction and min/mod limiters
 */

void RemapFluxPSG(const Real *U, const Real eps,
               const int jinner, const int jouter, Real *Flux)
{
  int j,jl,ju;
  Real dUc,dUl,dUr,dUm,lim_slope;

/* jinner,jouter are index range over which flux must be returned.  Set loop
 * limits depending on direction of upwind differences  */

  if (eps > 0.0) { /* eps always > 0 for inner i boundary */
    jl = jinner-1;
    ju = jouter-1;
  } else {         /* eps always < 0 for outer i boundary */
    jl = jinner;
    ju = jouter;
  }

  for (j=jl; j<=ju; j++) {
      dUc = U[j+1] - U[j-1];
      dUl = U[j  ] - U[j-1];
      dUr = U[j+1] - U[j  ];

      dUm = 0.0;
      if (dUl*dUr > 0.0) {
        lim_slope = MIN(fabs(dUl),fabs(dUr));
        dUm = SIGN(dUc)*MIN(0.5*fabs(dUc),2.0*lim_slope);
      }
 
    if (eps > 0.0) { /* eps always > 0 for inner i boundary */
      Flux[j+1] = eps*(U[j] + 0.5*(1.0 - eps)*dUm);
    } else {         /* eps always < 0 for outer i boundary */
      Flux[j  ] = eps*(U[j] - 0.5*(1.0 + eps)*dUm);
    }
  }

  return;
}

#endif /* SECOND_ORDER */

#if defined(THIRD_ORDER_CHAR) || defined(THIRD_ORDER_PRIM)
/*----------------------------------------------------------------------------*/
/*! \fn void RemapFluxPSG(const Real *U, const Real eps,
 *             const int jinner, const int jouter, Real *Flux)
 *  \brief third order reconstruction for conservative remap 
 *   using Colella & Sekora extremum preserving algorithm (PPME)
 */

void RemapFluxPSG(const Real *U, const Real eps,
               const int jinner, const int jouter, Real *Flux)
{
  int j,jl,ju;
  Real d2Uc,d2Ul,d2Ur,d2U,d2Ulim,lim_slope,Ulv,Urv,dU,U6,qa,qb,qc,qx;

/* jinner,jouter are index range over which flux must be returned.  Set loop
 * limits depending on direction of upwind differences  */

  if (eps > 0.0) { /* eps always > 0 for inner i boundary */
    jl = jinner-1;
    ju = jouter-1;
  } else {         /* eps always < 0 for outer i boundary */
    jl = jinner;
    ju = jouter;
  }

  for (j=jl; j<=ju+1; j++) {
    Uhalf[j]=(7.0*(U[j-1]+U[j]) - (U[j-2]+U[j+1]))/12.0;
    d2Uc = 3.0*(U[j-1] - 2.0*Uhalf[j] + U[j]);
    d2Ul = (U[j-2] - 2.0*U[j-1] + U[j  ]);
    d2Ur = (U[j-1] - 2.0*U[j  ] + U[j+1]);
    d2Ulim = 0.0;
    lim_slope = MIN(fabs(d2Ul),fabs(d2Ur));
    if (d2Uc > 0.0 && d2Ul > 0.0 && d2Ur > 0.0) {
      d2Ulim = SIGN(d2Uc)*MIN(1.25*lim_slope,fabs(d2Uc));
    }
    if (d2Uc < 0.0 && d2Ul < 0.0 && d2Ur < 0.0) {
      d2Ulim = SIGN(d2Uc)*MIN(1.25*lim_slope,fabs(d2Uc));
    }
    Uhalf[j] = 0.5*((U[j-1]+U[j]) - d2Ulim/3.0);
  }

  for (j=jl; j<=ju; j++) {
    Ulv = Uhalf[j  ];
    Urv = Uhalf[j+1];

    qa = (Urv-U[j])*(U[j]-Ulv);
    qb = (U[j-1]-U[j])*(U[j]-U[j+1]);
    if (qa <= 0.0 && qb <= 0.0) {
      qc = 6.0*(U[j] - 0.5*(Ulv+Urv));
      d2U  = -2.0*qc;
      d2Uc = (U[j-1] - 2.0*U[j  ] + U[j+1]);
      d2Ul = (U[j-2] - 2.0*U[j-1] + U[j  ]);
      d2Ur = (U[j  ] - 2.0*U[j+1] + U[j+2]);
      d2Ulim = 0.0;
      lim_slope = MIN(fabs(d2Ul),fabs(d2Ur));
      lim_slope = MIN(fabs(d2Uc),lim_slope);
      if (d2Uc > 0.0 && d2Ul > 0.0 && d2Ur > 0.0 && d2U > 0.0) {
        d2Ulim = SIGN(d2U)*MIN(1.25*lim_slope,fabs(d2U));
      }
      if (d2Uc < 0.0 && d2Ul < 0.0 && d2Ur < 0.0 && d2U < 0.0) {
        d2Ulim = SIGN(d2U)*MIN(1.25*lim_slope,fabs(d2U));
      }
      if (d2U == 0.0) {
        Ulv = U[j];
        Urv = U[j];
      } else {
        Ulv = U[j] + (Ulv - U[j])*d2Ulim/d2U;
        Urv = U[j] + (Urv - U[j])*d2Ulim/d2U;
      }
    }

    qa = (Urv-U[j])*(U[j]-Ulv);
    qb = Urv-Ulv;
    qc = 6.0*(U[j] - 0.5*(Ulv+Urv));
    if (qa <= 0.0) {
      Ulv = U[j];
      Urv = U[j];
    } else if ((qb*qc) > (qb*qb)) {
      Ulv = 3.0*U[j] - 2.0*Urv;
    } else if ((qb*qc) < -(qb*qb)) {
      Urv = 3.0*U[j] - 2.0*Ulv;
    }

    dU = Urv - Ulv;
    U6 = 6.0*(U[j] - 0.5*(Ulv + Urv));

    if (eps > 0.0) { /* eps always > 0 for inner i boundary */
      qx = TWO_3RDS*eps;
      Flux[j+1] = eps*(Urv - 0.75*qx*(dU - (1.0 - qx)*U6));

    } else {         /* eps always < 0 for outer i boundary */
      qx = -TWO_3RDS*eps;
      Flux[j  ] = eps*(Ulv + 0.75*qx*(dU + (1.0 - qx)*U6));
    }
  }

  return;
}

#endif /* THIRD_ORDER */

#endif /* SHEARING_BOX */

#endif /* PARTICLE_SELF_GRAVITY */
