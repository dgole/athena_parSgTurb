#include "../copyright.h"
/*===========================================================================*/
/*! \file self_gravity_particle.c
 *  \brief Provides routines for calculating the particle self-gravity
 *         potential and force.
 *
 * PURPOSE: provide routines for calculating particle self=gravity
 *          potential and force
 * 
 * CONTAINS PUBLIC FUNCTIONS:
 * - Get_Force_SG();  - contains methods to calculate self-gravity force
 *                      and map force to particle location
 *
 * PRIVATE FUNCTION PROTOTYPES:
 *   Calculate_Force_SG()             - calculate the force by differencing potential
 *   interp_force()                   - interpolates force from grid to particle location
 *   distr_mass()                     - maps particle mass to grid density 
 *   selfg_particle_fft_3d_init       - initialzation for periodic self-gravity solver
 *   selfg_particle_fft_obc_3d_init   - initialzation for outflow self-gravity solver
 *   selfg_particle_fft_disk_3d_init  - initialzation for disk (x3 outflow) self-gravity solver
 *
 *   Written by Jacob B. Simon						      */
/*============================================================================*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../defs.h"
#include "../athena.h"
#include "../prototypes.h"
#include "prototypes.h"
#include "particle.h"
#include "../globals.h"

#ifdef PARTICLE_SELF_GRAVITY         /* endif at the end of the file */

#ifndef FFT_ENABLED
#error self gravity with FFT requires configure --enable-fft
#endif /* FFT_ENABLED */

/* plans for forward and backward FFTs; work space for FFTW */
static struct ath_3d_fft_plan *fplan3d, *bplan3d;
static ath_fft_data *work=NULL, *work2=NULL;

#ifdef STATIC_MESH_REFINEMENT       /* SMR doesn't work with particle self-gravity */
#error particle self gravity with FFT not yet implemented to work with SMR
#endif


/*====================================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 *   Calculate_Force_SG()             - calculate the force by differencing potential
 *   interp_force()                   - interpolates force from grid to particle location
 *   distr_mass()                     - maps particle mass to grid density 
 *   selfg_particle_fft_3d_init       - initialzation for periodic self-gravity solver
 *   selfg_particle_fft_obc_3d_init   - initialzation for outflow self-gravity solver
 *====================================================================================*/
//void Calculate_Force_SG(GridS *pG);
int interp_force(GridS *pG, Real weight[3][3][3], int is, int js, int ks, Real *gxp, Real *gyp, Real *gzp);
void distr_mass(DomainS *pD, Real mass, Real weight[3][3][3], int is, int js, int ks);
void selfg_particle_fft_3d_init(MeshS *pM);
void selfg_particle_fft_obc_3d_init(MeshS *pM);

/*----------------------------------------------------------------------------*/
/*! \fn void Get_Force_SG(GridS *pG, Real x1, Real x2, Real x3, Real3Vect cell1)
 *  \brief Gets force from self-gravity and maps it to particle location.
 *
 * Input: pG, location of particles (x1, x2, x3), cell1
 * Output: Real3Vect containing force due to self-gravity
 */

Real3Vect Get_Force_SG(GridS *pG, Real x1, Real x2, Real x3,
                       Real3Vect cell1)
{
  int is, js, ks;
  Real3Vect fg;
  Real gxp, gyp, gzp;
  Real weight[3][3][3];

  fg.x1 = 0.; fg.x2 = 0.0; fg.x3 = 0.0;

/*
 * Interpolate the grid values of the forces to the particle location
 */
  getweight(pG, x1, x2, x3, cell1, weight, &is, &js, &ks);  
  if (interp_force(pG, weight, is, js, ks, &gxp, &gyp, &gzp) == 0) 
  {
/* Here is the force to be applied to the particles */
    fg.x1 = gxp;
    fg.x2 = gyp;
    fg.x3 = gzp;
  } else {
    fg.x1 = 0.0;
    fg.x2 = 0.0;
    fg.x3 = 0.0;
  }

  return fg;


}

Real3Vect Get_Force_SG_Direct(GridS *pG, Real x1, Real x2, Real x3, int nparticle_root)
{
  int p;
  GrainS *gr;
  Real3Vect fg;
  Real GRAV, mass, r3, soft_length;
  GRAV = four_pi_G_par/(4.*PI);

//  soft_length = 0.5*MIN(pG->dx1,MIN(pG->dx2,pG->dx3));
  soft_length = 0.;

  fg.x1 = 0.; fg.x2 = 0.; fg.x3 = 0.;
  p = 0;
  while (p<pG->nparticle) {
    if (p != nparticle_root) {
       gr = &(pG->particle[p]);
       mass = grproperty[gr->property].m;
       r3 = pow(sqrt(SQR(gr->x1-x1)+SQR(gr->x2-x2)+SQR(gr->x3-x3)+SQR(soft_length)),3.);
       fg.x1 += GRAV*mass/r3*(gr->x1-x1);
       fg.x2 += GRAV*mass/r3*(gr->x2-x2);
       fg.x3 += GRAV*mass/r3*(gr->x3-x3);
    }
    p++;
  }
/* Here is the force to be applied to the particles */

  return fg;

}

/*----------------------------------------------------------------------------*/
/*! \fn void Calculate_Force_SG(GridS *pG)
 *  \brief Calculates force by differencing potential
 *
 * Input: pG
 */

void Calculate_Force_SG(DomainS *pD)
{

  GridS *pG = (pD->Grid);
/* Here, we take par_Phi and calculate the force at the grid center */
  int i, is = pG->is, ie = pG->ie;
  int j, js = pG->js, je = pG->je;
  int k, ks = pG->ks, ke = pG->ke;


/* Calculate x,y,z forces from face centered potential
 * Divide by 4 PI G as is done for the normal fluid self-grav force calculation
 */ 

  for (k=klp+2; k<=kup-1; k++) {
  for (j=jlp+2; j<=jup-1; j++) {
  for (i=ilp+2; i<=iup-1; i++) { 
    pG->par_gx[k][j][i] = -(pG->par_Phi[k][j][i+1]-pG->par_Phi[k][j][i-1])/(2.*pG->dx1);
    pG->par_gy[k][j][i] = -(pG->par_Phi[k][j+1][i]-pG->par_Phi[k][j-1][i])/(2.*pG->dx2);
    pG->par_gz[k][j][i] = -(pG->par_Phi[k+1][j][i]-pG->par_Phi[k-1][j][i])/(2.*pG->dx3);
  }}}

  return;

}

/*----------------------------------------------------------------------------*/
/*! \fn int interp_force(GridS *pG Real weight[3][3][3], int is, int js, int ks,
 *                         Real *gxp, Real *gyp, Real *gzp)
 *  \brief Interpolates force from grid cell center to location of particles
 *
 * Input: pG, weighting array (weight), grid cells location (is, js, ks),
 * Output: Returns integer determining success of routine.  Also, 
 *         gxp, gyp, gzp are forces at particle location.       
 */

int interp_force(GridS *pG, Real weight[3][3][3], int is, int js, int ks,
                  Real *gxp, Real *gyp, Real *gzp)
{

  int n0,i,j,k,i0,j0,k0,i1,j1,k1,i2,j2,k2;
  Real gx_tmp,gy_tmp,gz_tmp;           /* density and velocity of the fluid */
  Real totwei, totwei1;         /* total weight (in case of edge cells) */

  /* linear interpolation */
  gx_tmp = 0.0;  gy_tmp = 0.0;  gz_tmp = 0.0; 
  totwei = 0.0;         totwei1 = 1.0;

  /* Interpolate graviational force */
  /* Note: in lower dimensions only wei[0] is non-zero */
  n0 = ncell-1;
  k1 = MAX(ks, klp);    k2 = MIN(ks+n0, kup);
  j1 = MAX(js, jlp);    j2 = MIN(js+n0, jup);
  i1 = MAX(is, ilp);    i2 = MIN(is+n0, iup);

  for (k=k1; k<=k2; k++) {
    k0=k-k1;
    for (j=j1; j<=j2; j++) {
      j0=j-j1;
      for (i=i1; i<=i2; i++) {
        i0=i-i1;

        gx_tmp += weight[k0][j0][i0] * pG->par_gx[k][j][i];
        gy_tmp += weight[k0][j0][i0] * pG->par_gy[k][j][i];
        gz_tmp += weight[k0][j0][i0] * pG->par_gz[k][j][i];

        totwei += weight[k0][j0][i0];
      }
    }
  }
  if (totwei < TINY_NUMBER) /* particle lies out of the grid, warning! */
    return -1;

  totwei1 = 1.0/totwei;
  *gxp = gx_tmp*totwei1;     *gyp = gy_tmp*totwei1;       *gzp = gz_tmp*totwei1;

  return 0;

}

#ifdef PARTICLE_SELF_GRAVITY_USING_FFT

/*----------------------------------------------------------------------------*/
/*! \fn void selfg_particle_fft_3d(DomainS *pD)
 *  \brief Calculate potential for 3D, periodic boundaries using FFT
 * 
 * Input: DomainS pD
 *
 * Notes: Only works for uniform grid, periodic boundary conditions
 */

void selfg_particle_fft_3d(DomainS *pD)
{
  GridS *pG = (pD->Grid);
  GrainS *gr;
  int i, is = pG->is, ie = pG->ie;
  int j, js = pG->js, je = pG->je;
  int k, ks = pG->ks, ke = pG->ke;
  int isw,jsw,ksw;
  long p;
  Real dx1sq=(pG->dx1*pG->dx1),dx2sq=(pG->dx2*pG->dx2),dx3sq=(pG->dx3*pG->dx3);
  Real dkx,dky,dkz,pcoeff;
  Real weight[3][3][3];
  Real3Vect cell1;

#ifdef SHEARING_BOX
  Real qomt,Lx,Ly,dt;
  Real kxtdx;
  Real xmin,xmax;
  int ip,jp;
  int nx3=pG->Nx[2]+2*nghost;
  int nx2=pG->Nx[1]+2*nghost;
  int nx1=pG->Nx[0]+2*nghost;
  Real ***RollDen, ***UnRollPhi;

  if((RollDen=(Real***)calloc_3d_array(nx3,nx1,nx2,sizeof(Real)))==NULL)
    ath_error("[selfg_particle_fft_3d]: malloc returned a NULL pointer\n");
  if((UnRollPhi=(Real***)calloc_3d_array(nx3,nx1,nx2,sizeof(Real)))==NULL)
    ath_error("[selfg_particle_fft_3d]: malloc returned a NULL pointer\n");

  xmin = pD->RootMinX[0];
  xmax = pD->RootMaxX[0];
  Lx = xmax - xmin;

  xmin = pD->RootMinX[1];
  xmax = pD->RootMaxX[1];
  Ly = xmax - xmin;

  dt = pG->time-((int)(qshear*Omega_0*pG->time*Lx/Ly))*Ly/(qshear*Omega_0*Lx);
  qomt = qshear*Omega_0*dt;
#endif

/* Here, we map the mass of the particles to a grid density */

  /* convenient expressions */
  if (pG->Nx[0] > 1)  cell1.x1 = 1.0/pG->dx1;
  else                cell1.x1 = 0.0;

  if (pG->Nx[1] > 1)  cell1.x2 = 1.0/pG->dx2;
  else                cell1.x2 = 0.0;

  if (pG->Nx[2] > 1)  cell1.x3 = 1.0/pG->dx3;
  else                cell1.x3 = 0.0;

  /* initialization */
  for (k=klp; k<=kup; k++)
    for (j=jlp; j<=jup; j++)
      for (i=ilp; i<=iup; i++) {
        pG->Coup[k][j][i].grid_d = 0.0;
      }

  for (p=0; p<pG->nparticle; p++)
  {/* loop over all particles */
    gr = &(pG->particle[p]);

/* Get weight */
    getweight(pG, gr->x1, gr->x2, gr->x3, cell1, weight, &isw, &jsw, &ksw);
   
/* Calculate particle density on the grid */ 
    distr_mass(pD,grproperty[gr->property].m,weight,isw,jsw,ksw);
  } 

/* deposit ghost zone values into the boundary zones */
  exchange_gpcouple(pD, 0);

#ifdef SHEARING_BOX
  for (k=ks-nghost; k<=ke+nghost; k++){
  for (j=js-nghost; j<=je+nghost; j++){
    for (i=is-nghost; i<=ie+nghost; i++){
      RollDen[k][i][j] = pG->Coup[k][j][i].grid_d;
    }
  }}
#endif

/* Forward FFT of 4\piG*d */

/* For shearing-box, need to roll density to the nearest periodic point */
#ifdef SHEARING_BOX
  RemapVar(pD,RollDen,-dt);
#endif
  for (k=ks; k<=ke; k++){
  for (j=js; j<=je; j++){
    for (i=is; i<=ie; i++){
      work[F3DI(i-is,j-js,k-ks,pG->Nx[0],pG->Nx[1],pG->Nx[2])][0] = 
#ifdef SHEARING_BOX
        RollDen[k][i][j];
#else
        pG->Coup[k][j][i].grid_d;
#endif
      work[F3DI(i-is,j-js,k-ks,pG->Nx[0],pG->Nx[1],pG->Nx[2])][1] = 0.0;
    }
  }}

   ath_3d_fft(fplan3d, work);
     
/* Compute potential in Fourier space.  Multiple loops are used to avoid divide
 * by zero at i=is,j=js,k=ks, and to avoid if statement in loop   */
/* To compute kx,ky,kz, note that indices relative to whole Domain are needed */

  dkx = 2.0*PI/(double)(pD->Nx[0]);
  dky = 2.0*PI/(double)(pD->Nx[1]);
  dkz = 2.0*PI/(double)(pD->Nx[2]);

#ifdef SHEARING_BOX
  ip=KCOMP(0,pG->Disp[0],pD->Nx[0]);
  jp=KCOMP(0,pG->Disp[1],pD->Nx[1]);
  kxtdx  = (ip+qomt*Lx/Ly*jp)*dkx;
#endif

  if ((pG->Disp[2])==0 && (pG->Disp[1])==0 && (pG->Disp[0])==0) {
    work[F3DI(0,0,0,pG->Nx[0],pG->Nx[1],pG->Nx[2])][0] = 0.0;
    work[F3DI(0,0,0,pG->Nx[0],pG->Nx[1],pG->Nx[2])][1] = 0.0;
  } else {
#ifdef SHEARING_BOX
    pcoeff = 1.0/(((2.0*cos( kxtdx           )-2.0)/dx1sq) +
                  ((2.0*cos((pG->Disp[1])*dky)-2.0)/dx2sq) +
                  ((2.0*cos((pG->Disp[2])*dkz)-2.0)/dx3sq));
#else
    pcoeff = 1.0/(((2.0*cos((pG->Disp[0])*dkx)-2.0)/dx1sq) +
                  ((2.0*cos((pG->Disp[1])*dky)-2.0)/dx2sq) +
                  ((2.0*cos((pG->Disp[2])*dkz)-2.0)/dx3sq));
#endif
    work[F3DI(0,0,0,pG->Nx[0],pG->Nx[1],pG->Nx[2])][0] *= pcoeff;
    work[F3DI(0,0,0,pG->Nx[0],pG->Nx[1],pG->Nx[2])][1] *= pcoeff;
  }


  for (k=ks+1; k<=ke; k++){
#ifdef SHEARING_BOX
    pcoeff = 1.0/(((2.0*cos( kxtdx                    )-2.0)/dx1sq) +
                  ((2.0*cos((        pG->Disp[1] )*dky)-2.0)/dx2sq) +
                  ((2.0*cos(( (k-ks)+pG->Disp[2] )*dkz)-2.0)/dx3sq));
#else
    pcoeff = 1.0/(((2.0*cos((        pG->Disp[0] )*dkx)-2.0)/dx1sq) +
                  ((2.0*cos((        pG->Disp[1] )*dky)-2.0)/dx2sq) +
                  ((2.0*cos(( (k-ks)+pG->Disp[2] )*dkz)-2.0)/dx3sq));
#endif
    work[F3DI(0,0,k-ks,pG->Nx[0],pG->Nx[1],pG->Nx[2])][0] *= pcoeff;
    work[F3DI(0,0,k-ks,pG->Nx[0],pG->Nx[1],pG->Nx[2])][1] *= pcoeff;
  }

  for (j=js+1; j<=je; j++){
    for (k=ks; k<=ke; k++){
#ifdef SHEARING_BOX
      jp=KCOMP(j-js ,pG->Disp[1],pD->Nx[1]);
      kxtdx  = (ip+qomt*Lx/Ly*jp)*dkx;
      pcoeff = 1.0/(((2.0*cos( kxtdx                    )-2.0)/dx1sq) +
                    ((2.0*cos(( (j-js)+pG->Disp[1] )*dky)-2.0)/dx2sq) +
                    ((2.0*cos(( (k-ks)+pG->Disp[2] )*dkz)-2.0)/dx3sq));
#else
      pcoeff = 1.0/(((2.0*cos((        pG->Disp[0] )*dkx)-2.0)/dx1sq) +
                    ((2.0*cos(( (j-js)+pG->Disp[1] )*dky)-2.0)/dx2sq) +
                    ((2.0*cos(( (k-ks)+pG->Disp[2] )*dkz)-2.0)/dx3sq));
#endif
      work[F3DI(0,j-js,k-ks,pG->Nx[0],pG->Nx[1],pG->Nx[2])][0] *= pcoeff;
      work[F3DI(0,j-js,k-ks,pG->Nx[0],pG->Nx[1],pG->Nx[2])][1] *= pcoeff;
    }
  }

  for (i=is+1; i<=ie; i++){
  for (j=js; j<=je; j++){
    for (k=ks; k<=ke; k++){
#ifdef SHEARING_BOX
      ip=KCOMP(i-is ,pG->Disp[0],pD->Nx[0]);
      jp=KCOMP(j-js ,pG->Disp[1],pD->Nx[1]);
      kxtdx  = (ip+qomt*Lx/Ly*jp)*dkx;
      pcoeff = 1.0/(((2.0*cos( kxtdx                    )-2.0)/dx1sq) +
                    ((2.0*cos(( (j-js)+pG->Disp[1] )*dky)-2.0)/dx2sq) +
                    ((2.0*cos(( (k-ks)+pG->Disp[2] )*dkz)-2.0)/dx3sq));
#else
      pcoeff = 1.0/(((2.0*cos(( (i-is)+pG->Disp[0] )*dkx)-2.0)/dx1sq) +
                    ((2.0*cos(( (j-js)+pG->Disp[1] )*dky)-2.0)/dx2sq) +
                    ((2.0*cos(( (k-ks)+pG->Disp[2] )*dkz)-2.0)/dx3sq));
#endif
      work[F3DI(i-is,j-js,k-ks,pG->Nx[0],pG->Nx[1],pG->Nx[2])][0] *= pcoeff;
      work[F3DI(i-is,j-js,k-ks,pG->Nx[0],pG->Nx[1],pG->Nx[2])][1] *= pcoeff;
    }
  }}

/* Backward FFT and set potential in real space.  Normalization of Phi is over
 * total number of cells in Domain */

  ath_3d_fft(bplan3d, work);

  for (k=ks; k<=ke; k++){
  for (j=js; j<=je; j++){
    for (i=is; i<=ie; i++){
#ifdef SHEARING_BOX
      UnRollPhi[k][i][j] = 
       four_pi_G_par*work[F3DI(i-is,j-js,k-ks,pG->Nx[0],pG->Nx[1],pG->Nx[2])][0]
        / bplan3d->gcnt;
#else
      pG->par_Phi[k][j][i] =
       four_pi_G_par*work[F3DI(i-is,j-js,k-ks,pG->Nx[0],pG->Nx[1],pG->Nx[2])][0]
        / bplan3d->gcnt;
#endif
    }
  }}

#ifdef SHEARING_BOX
  RemapVar(pD,UnRollPhi,dt);

  for (k=ks; k<=ke; k++){
    for (j=js; j<=je; j++){
      for (i=is; i<=ie; i++){
         pG->par_Phi[k][j][i] = UnRollPhi[k][i][j];
      }
    }
  }

  free_3d_array(RollDen);
  free_3d_array(UnRollPhi);
#endif

  return;
}

#endif /* PARTICLE_SELF_GRAVITY_USING_FFT */

#ifdef PARTICLE_SELF_GRAVITY_USING_FFT_OBC

/*----------------------------------------------------------------------------*/
/*! \fn void selfg_particle_fft_obc_3d(DomainS *pD)
 *  \brief Calculate potential for 3D, outflow boundaries using FFT
 * 
 * Input: DomainS pD
 *
 * Notes: Only works for uniform grid, outflow boundary conditions
 */

void selfg_particle_fft_obc_3d(DomainS *pD)
{
  GridS *pG = (pD->Grid);
  GrainS *gr;
  int i,ioff,ip, is = pG->is, ie = pG->ie;
  int j,joff,jp, js = pG->js, je = pG->je;
  int k,koff,kp, ks = pG->ks, ke = pG->ke;
  int isw,jsw,ksw;
  long p;
  Real pcoeff,offset;
  int Nx1=pD->Nx[0], Nx2=pD->Nx[1], Nx3=pD->Nx[2];
  int hNx1=Nx1/2, hNx2=Nx2/2, hNx3=Nx3/2;
  Real idx1sq=1.0/SQR(pG->dx1),idx2sq=1.0/SQR(pG->dx2),idx3sq=1.0/SQR(pG->dx3);
  Real dkx=2.0*PI/(double)(Nx1),dky=2.0*PI/(double)(Nx2),dkz=2.0*PI/(double)(Nx3);
  Real weight[3][3][3];
  Real3Vect cell1;

/* Here, we map the mass of the particles to a grid density */

  /* convenient expressions */
  if (pG->Nx[0] > 1)  cell1.x1 = 1.0/pG->dx1;
  else                cell1.x1 = 0.0;

  if (pG->Nx[1] > 1)  cell1.x2 = 1.0/pG->dx2;
  else                cell1.x2 = 0.0;

  if (pG->Nx[2] > 1)  cell1.x3 = 1.0/pG->dx3;
  else                cell1.x3 = 0.0;

  /* initialization */
  for (k=klp; k<=kup; k++)
    for (j=jlp; j<=jup; j++)
      for (i=ilp; i<=iup; i++) {
        pG->Coup[k][j][i].grid_d = 0.0;
      }
  for (p=0; p<pG->nparticle; p++)
  {/* loop over all particles */
    gr = &(pG->particle[p]);

/* Get weight */
    getweight(pG, gr->x1, gr->x2, gr->x3, cell1, weight, &isw, &jsw, &ksw);

/* Calculate particle density on the grid */
    distr_mass(pD,grproperty[gr->property].m,weight,isw,jsw,ksw);
  }

/* deposit ghost zone values into the boundary zones */
  exchange_gpcouple(pD, 0);

  /* Loop over the index offsets (0 for even, 1 for odd) */
  for (koff=0; koff<=1; koff++) {
  for (joff=0; joff<=1; joff++) {
    for (ioff=0; ioff<=1; ioff++) {

      /* STEP 1: Forward FFT of 4\piG*d */
      /* Add gas density into work array 0. */
      for (k=ks; k<=ke; k++) {
        for (j=js; j<=je; j++) {
          for (i=is; i<=ie; i++) {
            work[F3DI(i-is,j-js,k-ks,pG->Nx[0],pG->Nx[1],pG->Nx[2])][0] = pG->Coup[k][j][i].grid_d;
          }
        }
      }

      /* Copy work array 0 into work array 1, then multiply by complex offsets.
       * To compute offsets, note that indices relative to whole Domain are
       * needed. */
      for (k=ks; k<=ke; k++) {
        for (j=js; j<=je; j++) {
          for (i=is; i<=ie; i++) {
            work[F3DI(i-is,j-js,k-ks,pG->Nx[0],pG->Nx[1],pG->Nx[2])][1] =
              work[F3DI(i-is,j-js,k-ks,pG->Nx[0],pG->Nx[1],pG->Nx[2])][0];

            offset = 0.5*((i-is+pG->Disp[0])*ioff*dkx
                        + (j-js+pG->Disp[1])*joff*dky
                        + (k-ks+pG->Disp[2])*koff*dkz);
            work[F3DI(i-is,j-js,k-ks,pG->Nx[0],pG->Nx[1],pG->Nx[2])][0] *=  cos(offset);
            work[F3DI(i-is,j-js,k-ks,pG->Nx[0],pG->Nx[1],pG->Nx[2])][1] *= -sin(offset);
          }
        }
      }

      /* Forward FFT */
      ath_3d_fft(fplan3d, work);

      /* STEP 2:  Compute potential in Fourier space.  Multiple loops are used
       * to avoid divide by zero at i=is,j=js,k=ks, and to avoid if statement in
       * loop.  To compute kx,ky,kz, note that indices relative to whole Domain
       * are needed. */
      if ((pG->Disp[0])==0 && (pG->Disp[1])==0 && (pG->Disp[2])==0
           && ioff==0 && joff==0 && koff==0) {
        work[F3DI(0,0,0,pG->Nx[0],pG->Nx[1],pG->Nx[2])][0] = 0.0;
        work[F3DI(0,0,0,pG->Nx[0],pG->Nx[1],pG->Nx[2])][1] = 0.0;
      } else {
        pcoeff = -0.5/((1.0-cos((pG->Disp[0] + 0.5*ioff)*dkx))*idx1sq +
                       (1.0-cos((pG->Disp[1] + 0.5*joff)*dky))*idx2sq +
                       (1.0-cos((pG->Disp[2] + 0.5*koff)*dkz))*idx3sq);
        work[F3DI(0,0,0,pG->Nx[0],pG->Nx[1],pG->Nx[2])][0] *= pcoeff;
        work[F3DI(0,0,0,pG->Nx[0],pG->Nx[1],pG->Nx[2])][1] *= pcoeff;
      }

      for (k=ks+1; k<=ke; k++) {
        pcoeff = -0.5/((1.0-cos((         pG->Disp[0] + 0.5*ioff)*dkx))*idx1sq +
                       (1.0-cos((         pG->Disp[1] + 0.5*joff)*dky))*idx2sq +
                       (1.0-cos(((k-ks) + pG->Disp[2] + 0.5*koff)*dkz))*idx3sq);
        work[F3DI(0,0,k-ks,pG->Nx[0],pG->Nx[1],pG->Nx[2])][0] *= pcoeff;
        work[F3DI(0,0,k-ks,pG->Nx[0],pG->Nx[1],pG->Nx[2])][1] *= pcoeff;
      }

      for (j=js+1; j<=je; j++) {
        for (k=ks; k<=ke; k++) {
          pcoeff = -0.5/((1.0-cos((         pG->Disp[0] + 0.5*ioff)*dkx))*idx1sq +
                         (1.0-cos(((j-js) + pG->Disp[1] + 0.5*joff)*dky))*idx2sq +
                         (1.0-cos(((k-ks) + pG->Disp[2] + 0.5*koff)*dkz))*idx3sq);
          work[F3DI(0,j-js,k-ks,pG->Nx[0],pG->Nx[1],pG->Nx[2])][0] *= pcoeff;
          work[F3DI(0,j-js,k-ks,pG->Nx[0],pG->Nx[1],pG->Nx[2])][1] *= pcoeff;
        }
      }

      for (i=is+1; i<=ie; i++) {
        for (j=js; j<=je; j++) {
          for (k=ks; k<=ke; k++) {
            pcoeff = -0.5/((1.0-cos(((i-is) + pG->Disp[0] + 0.5*ioff)*dkx))*idx1sq +
                           (1.0-cos(((j-js) + pG->Disp[1] + 0.5*joff)*dky))*idx2sq +
                           (1.0-cos(((k-ks) + pG->Disp[2] + 0.5*koff)*dkz))*idx3sq);
            work[F3DI(i-is,j-js,k-ks,pG->Nx[0],pG->Nx[1],pG->Nx[2])][0] *= pcoeff;
            work[F3DI(i-is,j-js,k-ks,pG->Nx[0],pG->Nx[1],pG->Nx[2])][1] *= pcoeff;
          }
        }
      }


      /* STEP 3:  Backward FFT and set potential in real space */
      ath_3d_fft(bplan3d, work);

      /* Multiply by complex offsets and add real part to Phi.  To compute
       * offsets, note that indices relative to whole Domain are needed. */
      for (k=ks; k<=ke; k++) {
        for (j=js; j<=je; j++) {
          for (i=is; i<=ie; i++) {
            offset = 0.5*((i-is+pG->Disp[0])*ioff*dkx
                        + (j-js+pG->Disp[1])*joff*dky
                        + (k-ks+pG->Disp[2])*koff*dkz);
            pG->par_Phi[k][j][i] +=
              cos(offset)*work[F3DI(i-is,j-js,k-ks,pG->Nx[0],pG->Nx[1],pG->Nx[2])][0]
            - sin(offset)*work[F3DI(i-is,j-js,k-ks,pG->Nx[0],pG->Nx[1],pG->Nx[2])][1];
          }
        }
      }
    }
  }}

  /* Finally, normalize the transforms. NOTE: There is an 8.0 here because
   * formally the transforms are performed over the extended domain of size
   * (2Nx)(2Ny)(2Nz). */
  pcoeff = four_pi_G_par/(8.0*bplan3d->gcnt);
  for (k=ks; k<=ke; k++){
    for (j=js; j<=je; j++){
      for (i=is; i<=ie; i++){
        pG->par_Phi[k][j][i] *= pcoeff;
      }
    }
  }

  return;
}

#endif /* PARTICLE_SELF_GRAVITY_USING_FFT_OBC */

#ifdef PARTICLE_SELF_GRAVITY_USING_FFT_DISK

/*----------------------------------------------------------------------------*/
/*! \fn void selfg_particle_fft_disk_3d(DomainS *pD)
 *  \brief Calculate potential for 3D disk boundaries using FFT
 *         Periodic boundary conditions in x1 and x2; open bc in x3
 * 
 * Input: DomainS pD
 *
 * Notes: Only works for uniform grid, outflow x3, periodic x1, x2
 */

void selfg_particle_fft_disk_3d(DomainS *pD)
{
  GridS *pG = (pD->Grid);
  GrainS *gr;
  int i, is = pG->is, ie = pG->ie;
  int j, js = pG->js, je = pG->je;
  int k, ks = pG->ks, ke = pG->ke;
  int ip, jp;
  int isw,jsw,ksw;
  long p;
  Real kxtdx,kydy;
  Real dkx,dky,dkz;
  Real dx1sq=(pG->dx1*pG->dx1),dx2sq=(pG->dx2*pG->dx2),dx3sq=(pG->dx3*pG->dx3);
  Real xmin,xmax;
  Real Lperp,den;
  Real ***Acoeff=NULL,***Bcoeff=NULL;
  Real weight[3][3][3];
  Real3Vect cell1;

#ifdef SHEARING_BOX
  int nx3=pG->Nx[2]+2*nghost;
  int nx2=pG->Nx[1]+2*nghost;
  int nx1=pG->Nx[0]+2*nghost;
  Real ***RollDen=NULL, ***UnRollPhi=NULL;
  Real Lx,Ly,qomt,dt;

  if((RollDen=(Real***)calloc_3d_array(nx3,nx1,nx2,sizeof(Real)))==NULL)
    ath_error("[selfg_fft_disk]: malloc returned a NULL pointer\n");
  if((UnRollPhi=(Real***)calloc_3d_array(nx3,nx1,nx2,sizeof(Real)))==NULL)
    ath_error("[selfg_fft_disk]: malloc returned a NULL pointer\n");

  xmin = pD->RootMinX[0];
  xmax = pD->RootMaxX[0];
  Lx = xmax - xmin;

  xmin = pD->RootMinX[1];
  xmax = pD->RootMaxX[1];
  Ly = xmax - xmin;

  dt = pG->time-((int)(qshear*Omega_0*pG->time*Lx/Ly))*Ly/(qshear*Omega_0*Lx);
  qomt = qshear*Omega_0*dt;
#endif

/* Here, we map the mass of the particles to a grid density */

  /* convenient expressions */
  if (pG->Nx[0] > 1)  cell1.x1 = 1.0/pG->dx1;
  else                cell1.x1 = 0.0;

  if (pG->Nx[1] > 1)  cell1.x2 = 1.0/pG->dx2;
  else                cell1.x2 = 0.0;

  if (pG->Nx[2] > 1)  cell1.x3 = 1.0/pG->dx3;
  else                cell1.x3 = 0.0;

  /* initialization */
  for (k=klp; k<=kup; k++)
    for (j=jlp; j<=jup; j++)
      for (i=ilp; i<=iup; i++) {
        pG->Coup[k][j][i].grid_d = 0.0;
      }

  for (p=0; p<pG->nparticle; p++)
  {/* loop over all particles */
    gr = &(pG->particle[p]);

/* Get weight */
    getweight(pG, gr->x1, gr->x2, gr->x3, cell1, weight, &isw, &jsw, &ksw);

/* Calculate particle density on the grid */
    distr_mass(pD,grproperty[gr->property].m,weight,isw,jsw,ksw);
  }

/* deposit ghost zone values into the boundary zones */
  exchange_gpcouple(pD, 0);

/* allocates memory for Acoeff and Bcoeff arrays */
  if ((Acoeff = (Real***)calloc_3d_array(pG->Nx[0],pG->Nx[1],pG->Nx[2],sizeof(Real))) == NULL)
    ath_error("[selfg_fft_disk]: malloc returned a NULL pointer\n");

  if ((Bcoeff = (Real***)calloc_3d_array(pG->Nx[0],pG->Nx[1],pG->Nx[2],sizeof(Real))) == NULL)
    ath_error("[selfg_fft_disk]: malloc returned a NULL pointer\n");

/* To compute kx,ky,kz, note indices relative to whole Domain are needed */
  dkx = 2.0*PI/(double)(pD->Nx[0]);
  dky = 2.0*PI/(double)(pD->Nx[1]);
  dkz = 2.0*PI/(double)(pD->Nx[2]);

/* This is size of whole Domain perpendicular to the plane (=disk thickness)*/
  xmin = pD->RootMinX[2];
  xmax = pD->RootMaxX[2];
  Lperp = xmax-xmin;

/* Compute potential coeffs in k space. Zero wavenumber is special
   case; need to avoid divide by zero */
  for (i=is; i<=ie; i++){
    for (j=js; j<=je; j++){
      for (k=ks; k<=ke; k++){
        ip=KCOMP(i-is,pG->Disp[0],pD->Nx[0]);
        jp=KCOMP(j-js,pG->Disp[1],pD->Nx[1]);
#ifdef SHEARING_BOX
        kxtdx = (ip+qomt*Lx/Ly*jp)*dkx;
#else
        kxtdx = ip*dkx;
#endif
        kydy = jp*dky;
        if (((k-ks)+pG->Disp[2])==0 && ((j-js)+pG->Disp[1])==0 && ((i-is)+pG->Disp[0])==0)
          Acoeff[0][0][0] = 0.0;
        else{
          Acoeff[i-is][j-js][k-ks] = 0.5*
            (1.0-exp(-sqrt(SQR(kxtdx)/dx1sq+SQR(kydy)/dx2sq)*Lperp))/
            (((2.0*cos(  kxtdx                 )-2.0)/dx1sq) +
             ((2.0*cos(  kydy                  )-2.0)/dx2sq) +
             ((2.0*cos(((k-ks)+pG->Disp[2])*dkz)-2.0)/dx3sq));
        }
        Bcoeff[i-is][j-js][k-ks] = 0.5*
          (1.0+exp(-sqrt(SQR(kxtdx)/dx1sq+SQR(kydy)/dx2sq)*Lperp))/
          (((2.0*cos(      kxtdx                 )-2.0)/dx1sq) +
           ((2.0*cos(      kydy                  )-2.0)/dx2sq) +
           ((2.0*cos((0.5+(k-ks)+pG->Disp[2])*dkz)-2.0)/dx3sq));
      }
    }
  }

/* Copy current potential into old */

#ifdef SHEARING_BOX
  for (k=ks-nghost; k<=ke+nghost; k++){
    for (j=js-nghost; j<=je+nghost; j++){
      for (i=is-nghost; i<=ie+nghost; i++){
        RollDen[k][i][j] = pG->Coup[k][j][i].grid_d;
      }
    }
  }
#endif

#ifdef SHEARING_BOX
  RemapVar(pD,RollDen,-dt);
#endif

/* Fill arrays of 4\piG*d and 4\piG*d *exp(-i pi x2/Lperp) */


  for (k=ks; k<=ke; k++){
    for (j=js; j<=je; j++){
      for (i=is; i<=ie; i++){
#ifdef SHEARING_BOX
        den=RollDen[k][i][j];
#else
        den=pG->Coup[k][j][i].grid_d;
#endif
        work[F3DI(i-is,j-js,k-ks,pG->Nx[0],pG->Nx[1],pG->Nx[2])][0] = den;
      }
    }
  }

  for (k=ks; k<=ke; k++){
    for (j=js; j<=je; j++){
      for (i=is; i<=ie; i++){
        work[F3DI(i-is,j-js,k-ks,pG->Nx[0],pG->Nx[1],pG->Nx[2])][0] *=four_pi_G_par;
        work[F3DI(i-is,j-js,k-ks,pG->Nx[0],pG->Nx[1],pG->Nx[2])][1] = 0.0;

        work2[F3DI(i-is,j-js,k-ks,pG->Nx[0],pG->Nx[1],pG->Nx[2])][0] =
          cos(0.5*((k-ks)+pG->Disp[2])*dkz)*
              work[F3DI(i-is,j-js,k-ks,pG->Nx[0],pG->Nx[1],pG->Nx[2])][0];
        work2[F3DI(i-is,j-js,k-ks,pG->Nx[0],pG->Nx[1],pG->Nx[2])][1] =
         -sin(0.5*((k-ks)+pG->Disp[2])*dkz)*
              work[F3DI(i-is,j-js,k-ks,pG->Nx[0],pG->Nx[1],pG->Nx[2])][0];
      }
    }
  }

/* Forward FFT of 4\piG*d and 4\piG*d *exp(-i pi x2/Lperp) */

  ath_3d_fft(fplan3d, work);
  ath_3d_fft(fplan3d, work2);

/* Compute potential in Fourier space, using pre-computed coefficients */

  for (i=is; i<=ie; i++){
    for (j=js; j<=je; j++){
      for (k=ks; k<=ke; k++){
        work[F3DI(i-is,j-js,k-ks,pG->Nx[0],pG->Nx[1],pG->Nx[2])][0] *=
          Acoeff[i-is][j-js][k-ks];
        work[F3DI(i-is,j-js,k-ks,pG->Nx[0],pG->Nx[1],pG->Nx[2])][1] *=
          Acoeff[i-is][j-js][k-ks];
        work2[F3DI(i-is,j-js,k-ks,pG->Nx[0],pG->Nx[1],pG->Nx[2])][0] *=
          Bcoeff[i-is][j-js][k-ks];
        work2[F3DI(i-is,j-js,k-ks,pG->Nx[0],pG->Nx[1],pG->Nx[2])][1] *=
          Bcoeff[i-is][j-js][k-ks];
      }
    }
  }


/* Backward FFT and set potential in real space.  Normalization of Phi is over
 * total number of cells in Domain */

  ath_3d_fft(bplan3d, work);
  ath_3d_fft(bplan3d, work2);

  for (k=ks; k<=ke; k++){
    for (j=js; j<=je; j++){
      for (i=is; i<=ie; i++){
#ifdef SHEARING_BOX
        UnRollPhi[k][i][j] =
#else
        pG->par_Phi[k][j][i] =
#endif
                           (work[F3DI(i-is,j-js,k-ks,pG->Nx[0],pG->Nx[1],pG->Nx[2])][0]
                         + cos(0.5*((k-ks)+pG->Disp[2])*dkz)*
                           work2[F3DI(i-is,j-js,k-ks,pG->Nx[0],pG->Nx[1],pG->Nx[2])][0]
                         - sin(0.5*((k-ks)+pG->Disp[2])*dkz)*
                           work2[F3DI(i-is,j-js,k-ks,pG->Nx[0],pG->Nx[1],pG->Nx[2])][1])/
                           bplan3d->gcnt;
      }
    }
  }

#ifdef SHEARING_BOX
  RemapVar(pD,UnRollPhi,dt);

  for (k=ks; k<=ke; k++){
    for (j=js; j<=je; j++){
      for (i=is; i<=ie; i++){
         pG->par_Phi[k][j][i] = UnRollPhi[k][i][j];
      }
    }
  }

  free_3d_array(RollDen);
  free_3d_array(UnRollPhi);
#endif
  free_3d_array(Acoeff);
  free_3d_array(Bcoeff);

  return;
}

#endif /* PARTICLE_SELF_GRAVITY_USING_FFT_DISK */

/*----------------------------------------------------------------------------*/
/*! \fn void distr_mass(DomainS *pD, Real mass, Real weight[3][3][3], 
 *                              int is, int js, int ks)
 *  \brief Maps particle mass onto a grid density
 *
 * Input: pG, mass of particle (mass), weighting array (weight), grid cell
 *        location (is,js,ks)
 *
 */
void distr_mass(DomainS *pD, Real mass, Real weight[3][3][3], int is, int js, int ks)
{
  GridS *pG = pD->Grid;
  int n0,i,j,k,i0,j0,k0,i1,j1,k1,i2,j2,k2;

  n0 = ncell-1;
  k1 = MAX(ks, klp);    k2 = MIN(ks+n0, kup);
  j1 = MAX(js, jlp);    j2 = MIN(js+n0, jup);
  i1 = MAX(is, ilp);    i2 = MIN(is+n0, iup);

  for (k=k1; k<=k2; k++) {
    k0 = k-k1;
    for (j=j1; j<=j2; j++) {
      j0 = j-j1;
      for (i=i1; i<=i2; i++) {
        i0 = i-i1;
        pG->Coup[k][j][i].grid_d += weight[k0][j0][i0] * mass;
        pG->Coup[k][j][i].grid_v1 += 0.;      
        pG->Coup[k][j][i].grid_v2 += 0.;      
        pG->Coup[k][j][i].grid_v3 += 0.;      
      }
    }
  }

  return;
}

#ifdef PARTICLE_SELF_GRAVITY_USING_FFT

/*----------------------------------------------------------------------------*/
/*! \fn void selfg_particle_fft_3d_init(MeshS *pM)
 *  \brief Initializes plans for forward/backward FFTs, and allocates memory 
 *  needed by FFTW.
 */

void selfg_particle_fft_3d_init(MeshS *pM)
{
  DomainS *pD;
  int nl,nd;
  for (nl=0; nl<(pM->NLevels); nl++){
    for (nd=0; nd<(pM->DomainsPerLevel[nl]); nd++){
      if (pM->Domain[nl][nd].Grid != NULL){
        pD = (DomainS*)&(pM->Domain[nl][nd]);
        fplan3d = ath_3d_fft_quick_plan(pD, NULL, ATH_FFT_FORWARD);
        bplan3d = ath_3d_fft_quick_plan(pD, NULL, ATH_FFT_BACKWARD);
        work = ath_3d_fft_malloc(fplan3d);
      }
    }
  }
}

#endif /* PARTICLE_SELF_GRAVITY_USING_FFT */

#ifdef PARTICLE_SELF_GRAVITY_USING_FFT_OBC

/*----------------------------------------------------------------------------*/
/*! \fn void selfg_particle_fft_obc_3d_init(MeshS *pM)
 *  \brief Initializes plans for forward/backward FFTs, and allocates memory 
 *   needed by FFTW.
 */
void selfg_particle_fft_obc_3d_init(MeshS *pM)
{
  DomainS *pD;
  int nl,nd;
  for (nl=0; nl<(pM->NLevels); nl++) {
    for (nd=0; nd<(pM->DomainsPerLevel[nl]); nd++) {
      if (pM->Domain[nl][nd].Grid != NULL) {
        pD = (DomainS*)&(pM->Domain[nl][nd]);
        fplan3d = ath_3d_fft_quick_plan(pD, NULL, ATH_FFT_FORWARD);
        bplan3d = ath_3d_fft_quick_plan(pD, NULL, ATH_FFT_BACKWARD);
        work = ath_3d_fft_malloc(fplan3d);
      }
    }
  }
}

#endif /* PARTICLE_SELF_GRAVITY_USING_FFT_OBC */

#ifdef PARTICLE_SELF_GRAVITY_USING_FFT_DISK

/*----------------------------------------------------------------------------*/
/*! \fn void selfg_particle_fft_disk_3d_init(MeshS *pM)
 *  \brief Initializes plans for forward/backward FFTs, and allocates memory 
 *  needed by FFTW.
 */

void selfg_particle_fft_disk_3d_init(MeshS *pM)
{
  DomainS *pD;
  int nl,nd;
  for (nl=0; nl<(pM->NLevels); nl++){
    for (nd=0; nd<(pM->DomainsPerLevel[nl]); nd++){
      if (pM->Domain[nl][nd].Grid != NULL){
        pD = (DomainS*)&(pM->Domain[nl][nd]);
        fplan3d = ath_3d_fft_quick_plan(pD, NULL, ATH_FFT_FORWARD);
        bplan3d = ath_3d_fft_quick_plan(pD, NULL, ATH_FFT_BACKWARD);
        work = ath_3d_fft_malloc(fplan3d);
        work2 = ath_3d_fft_malloc(fplan3d);
      }
    }
  }
}

#endif /* PARTICLE_SELF_GRAVITY_USING_FFT_DISK */

/*----------------------------------------------------------------------------*/
/*! \fn void particle_self_gravity_init(MeshS *pM)
 *  \brief Runs some checks and initialization routines
 */
void particle_self_gravity_init(MeshS *pM)
{
  int dim = 0;

/* Calculate the dimensions  */
  if(pM->Nx[0] > 1) dim++;
  if(pM->Nx[1] > 1) dim++;
  if(pM->Nx[2] > 1) dim++;

/* test that user set values for constants */
  if (four_pi_G_par < 0.0)
    ath_error("[particle_self_gravity_init] four_pi_G_par must be set >0 in prob generator\n");

/* Return function pointer based on dimensions and algorithm */

#ifdef PARTICLE_SELF_GRAVITY_USING_FFT
  selfg_particle_fft_3d_init(pM);
#endif
#ifdef PARTICLE_SELF_GRAVITY_USING_FFT_OBC
  selfg_particle_fft_obc_3d_init(pM);
#endif
#ifdef PARTICLE_SELF_GRAVITY_USING_FFT_DISK
  selfg_particle_fft_disk_3d_init(pM);
#endif

  return;

}

/*----------------------------------------------------------------------------*/
/*! \fn void particle_self_gravity_destruct(void)
 *  \brief Place holder for freeing up stuff, but nothing is currently done
 *         here.
 */
void particle_self_gravity_destruct(void)
{
  return;
}

#endif /* PARTICLE_SELF_GRAVITY */
