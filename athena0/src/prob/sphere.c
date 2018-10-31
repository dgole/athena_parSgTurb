#include "copyright.h"
/*============================================================================*/
/*! \file sphere.c
 *  \brief Problem generator for spherical collapse of 
 *   uniform density sphere. Gravity test. */
/*============================================================================*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"

/*============================================================================* 
 * PRIVATE FUNCTION PROTOTYPES:
 *============================================================================*/

static Real dsphere, damb, p0, bx0, by0, bz0, radius;

void problem(DomainS *pDomain)
{
  GridS *pGrid = (pDomain->Grid);
  int i,is,ie,j,js,je,k,ks,ke,nx1,nx2,nx3;
  int ki, kj, kk, count, err, ndim;
  Real x1,x2,x3,x1min,x1max,x2min,x2max,x3min,x3max;
  Real x1len,x2len,x3len,dx1,dx2,dx3,xc,yc,zc,r2;
  Real huge=1.0e60,radius=0.0;

#ifdef MPI_PARALLEL
  int myid = myID_Comm_world;
#else
  int myid = 0;
#endif

  damb          = par_getd("problem","damb"     );
  dsphere       = par_getd("problem","dsphere"  );
#ifndef ISOTHERMAL
  p0       = par_getd("problem","p0"  );
#endif
#ifdef MHD
  bx0      = par_getd("problem","bx0");
  by0      = par_getd("problem","by0");
  bz0      = par_getd("problem","bz0");
#endif
  xc      = par_getd("problem","xc" );
  yc      = par_getd("problem","yc" );
  zc      = par_getd("problem","zc" );
  radius   = par_getd("problem","radius");
#ifdef SELF_GRAVITY
  four_pi_G= par_getd("problem","four_pi_G");
#endif
  if (myid == 0) {
    fprintf(stdout,"[sphere]: damb           = %13.5e\n",damb   );
    fprintf(stdout,"[sphere]: dsphere        = %13.5e\n",dsphere);
#ifndef ISOTHERMAL
    fprintf(stdout,"[sphere]: p0        = %13.5e\n",p0);
#endif
#ifdef MHD
    fprintf(stdout,"[sphere]: bx0       = %13.5e\n",bx0);
    fprintf(stdout,"[sphere]: by0       = %13.5e\n",by0);
    fprintf(stdout,"[sphere]: bz0       = %13.5e\n",bz0);
#endif  
    fprintf(stdout,"[sphere]: xc        = %13.5e\n",xc);
    fprintf(stdout,"[sphere]: yc        = %13.5e\n",yc);
    fprintf(stdout,"[sphere]: zc        = %13.5e\n",zc);
    fprintf(stdout,"[sphere]: radius    = %13.5e\n",radius);
#ifdef SELF_GRAVITY
    fprintf(stdout,"[sphere]: four_pi_G = %13.5e\n",four_pi_G);
#endif
  }

  is = pGrid->is; ie = pGrid->ie;
  js = pGrid->js; je = pGrid->je;
  ks = pGrid->ks; ke = pGrid->ke;
  nx1 = (ie-is)+1; 
  nx2 = (je-js)+1;
  nx3 = (ke-ks)+1;
  ndim = 1;
  if (is == ie) {
    ath_error("[sphere]: This problem can only be run with Nx1>1\n");
  }
  nx1 += 2*nghost;
  if (nx2 > 1) {
    ndim++;
    nx2 += 2*nghost;
  }
  if (nx3 > 1) {
    ndim++;
    nx3 += 2*nghost;
  }
  if (myid == 0) {
    fprintf(stdout,"[sphere]: ndim      = %13d\n",ndim);
    fprintf(stdout,"[sphere]: nx1,2,3   = %5d %5d %5d\n",nx1,nx2,nx3);
    fprintf(stdout,"[sphere]: is,js,ks  = %5d %5d %5d\n",is,js,ks);
    fprintf(stdout,"[sphere]: ie,je,ke  = %5d %5d %5d\n",ie,je,ke);
  }

/* ONLY FOR CONSTANT GRIDS */
  cc_pos(pGrid,is,js,ks,&x1min,&x2min,&x3min);
  cc_pos(pGrid,ie,je,ke,&x1max,&x2max,&x3max);
  dx1   = pGrid->dx1;
  dx2   = pGrid->dx2;
  dx3   = pGrid->dx3;
  x1min = par_getd("domain1","x1min");
  x1max = par_getd("domain1","x1max");
  x2min = par_getd("domain1","x2min");
  x2max = par_getd("domain1","x2max");
  x3min = par_getd("domain1","x3min");
  x3max = par_getd("domain1","x3max");
  x1len =  x1max-x1min;
  x2len =  x2max-x2min;
  x3len =  x3max-x3min;
  if (myid == 0) {
    fprintf(stdout,"[sphere]: x1min = %13.5e, x2min = %13.5e, x3min = %13.5e\n",x1min,x2min,x3min);
    fprintf(stdout,"[sphere]: x1max = %13.5e, x2max = %13.5e, x3max = %13.5e\n",x1max,x2max,x3max);
    fprintf(stdout,"[sphere]: dx1   = %13.5e, dx2   = %13.5e, dx3   = %13.5e\n",dx1,dx2,dx3);
    fprintf(stdout,"[sphere]: x1len = %13.5e, x2len = %13.5e, x3len = %13.5e\n",x1len,x2len,x3len);
  }

  grav_mean_rho = damb;

  if (ndim == 3) { /* this is the 3D case */
    for (k=ks; k<=ke; k++) {
      for (j=js; j<=je; j++) {
        for (i=is; i<=ie; i++) {
          /* these are cell centers */
          cc_pos(pGrid,i,j,k,&x1,&x2,&x3);
          r2  = sqrt(SQR(x1-xc)+SQR(x2-yc)+SQR(x3-zc));          
          if (r2 <= radius) {
            pGrid->U[k][j][i].d = dsphere;
          } else {
            pGrid->U[k][j][i].d = damb;
          }     
          pGrid->U[k][j][i].M1 = 0.0;
          pGrid->U[k][j][i].M2 = 0.0;
          pGrid->U[k][j][i].M3 = 0.0;
#ifdef MHD
          pGrid->B1i[k][j][i]  = bx0;
          pGrid->B2i[k][j][i]  = by0;
          pGrid->B3i[k][j][i]  = bz0;
#endif
        } /* i */
      } /* j */
    } /* k */
    if (myid == 0) 
      fprintf(stdout,"[sphere]: 3D setup finished\n");
  }
 
/* boundary conditions on interface B */
#ifdef MHD
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      pGrid->B1i[k][j][ie+1] = bx0;
    }
  }
  for (k=ks; k<=ke; k++) {
    for (i=is; i<=ie; i++) {
      pGrid->B2i[k][je+1][i] = by0;
    }
  }
  if (ndim == 3) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        pGrid->B3i[ke+1][j][i] = bz0;
      }
    }
  }
#endif

/* initialize total energy and cell-centered B */

  if (ndim == 3) {
    /* 3D case */
    for (k=ks; k<=ke; k++) {
      for (j=js; j<=je; j++) {
        for (i=is; i<=ie; i++) {
#ifdef MHD
          pGrid->U[k][j][i].B1c= 0.5*(pGrid->B1i[k][j][i]+pGrid->B1i[k][j][i+1]);
          pGrid->U[k][j][i].B2c= 0.5*(pGrid->B2i[k][j][i]+pGrid->B2i[k][j+1][i]);
          pGrid->U[k][j][i].B3c= 0.5*(pGrid->B3i[k][j][i]+pGrid->B3i[k+1][j][i]);
#endif
#ifndef ISOTHERMAL
          pGrid->U[k][j][i].E = p0/Gamma_1
#ifdef MHD
              + 0.5*(SQR(pGrid->U[k][j][i].B1c) + SQR(pGrid->U[k][j][i].B2c)
                   + SQR(pGrid->U[k][j][i].B3c))
#endif
              + 0.5*(SQR(pGrid->U[k][j][i].M1) + SQR(pGrid->U[k][j][i].M2)
                   + SQR(pGrid->U[k][j][i].M3))/pGrid->U[k][j][i].d;
#endif
        }
      }
    }
  }

  return;
}

/*==============================================================================
 * PROBLEM USER FUNCTIONS:
 * problem_write_restart() - writes problem-specific user data to restart files
 * problem_read_restart()  - reads problem-specific user data from restart files
 * get_usr_expr()          - sets pointer to expression for special output data
 * Userwork_in_loop        - problem specific work IN     main loop
 * Userwork_after_loop     - problem specific work AFTER  main loop
 *----------------------------------------------------------------------------*/

void problem_write_restart(MeshS *pM, FILE *fp)
{
  return;
}

/* problem restart needs to enroll cooling function and boundary conditions. 
 */

void problem_read_restart(MeshS *pM, FILE *fp)
{
  return;
}

ConsFun_t get_usr_expr(const char *expr)
{
  return NULL;
}

VOutFun_t get_usr_out_fun(const char *name){
  return NULL;
}

void Userwork_in_loop(MeshS *pM)
{
  return;
}

void Userwork_after_loop(MeshS *pM)
{
  return;
}
