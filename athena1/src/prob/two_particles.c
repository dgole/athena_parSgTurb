#include "copyright.h"
/*============================================================================*/
/*! \file particle_sphere.c
 *  \brief Problem generator for spherical collapse of 
 *   uniform density sphere of particles. Particle Gravity test. */
/*============================================================================*/
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"
#include "particles/particle.h"

#ifndef PARTICLES
#error : The particle_sphere problem requires particles to be enabled.
#endif /* PARTICLES */

#ifndef ISOTHERMAL
#error : The particle_sphere problem requires isothermal equation of state.
#endif /* ISOTHERMAL */

#ifdef MHD
#error: The particle_sphere problem does not allow MHD
#endif /* MHD */ 

/*------------------------ filewide global variables -------------------------*/
/* Sphere variables */
static Real radius, mratio;
/* particle number variables */
int Npar,Npar3,downsamp,Nx;
/* output filename */
char name[50];

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 * property_mybin()  - particle property selection function
 *============================================================================*/
static double ran2(long int *idum);
static int property_mybin(const GrainS *gr, const GrainAux *grsub);

/*=========================== PUBLIC FUNCTIONS =================================
 *============================================================================*/
/*----------------------------------------------------------------------------*/
/* problem:   */

void problem(DomainS *pDomain)
{
  GridS *pGrid = pDomain->Grid;
  int i,j,k,ip,jp,kp,ixs,jxs,kxs;
  int ipert;
  long p;
  long iseed = -1; /* Initialize on the first call to ran2 */
  Real x1,x2,x3,x1l,x1u,x2l,x2u,x3l,x3u,x1p,x2p,x3p;
  Real x1min,x1max,x2min,x2max,x3min,x3max,Lx,Ly,Lz;
  Real r2,xc,yc,zc,dgas,offset;
  double rval;

  if (pDomain->Nx[2] == 1) {
    ath_error("[particle_sphere]: particle_sphere only works for 3D problem.\n");
  }

/* Initialize boxsize */
  x1min = pDomain->RootMinX[0];
  x1max = pDomain->RootMaxX[0];
  Lx = x1max - x1min;

  x2min = pDomain->RootMinX[1];
  x2max = pDomain->RootMaxX[1];
  Ly = x2max - x2min;

  x3min = pDomain->RootMinX[2];
  x3max = pDomain->RootMaxX[2];
  Lz = x3max - x3min;

  Nx = pDomain->Nx[0]; /* used for particle selection function */

/* Ensure a different initial random seed for each process in an MPI calc. */
  ixs = pGrid->Disp[0];
  jxs = pGrid->Disp[1];
  kxs = pGrid->Disp[2];
  iseed = -1 - (ixs + pDomain->Nx[0]*(jxs + pDomain->Nx[1]*kxs));

/* Read in perturbation mode.  If ipert != 1, no randomization is done.  
 * This is the default. */
  ipert = par_geti_def("problem","ipert", 0);
/* Read in offset from grid cell center */
  offset = par_getd_def("problem","offset",offset);

/* Read initial conditions */
  dgas = par_getd_def("problem","dgas",1.0);
  xc      = par_getd_def("problem","xc",0.);
  yc      = par_getd_def("problem","yc",0.);
  zc      = par_getd_def("problem","zc",0.);
  radius   = par_getd_def("problem","radius",0.5);
  mratio = par_getd_def("problem","mratio",2.0); /* total mass fraction */
#ifdef PARTICLE_SELF_GRAVITY
  four_pi_G_par = par_getd_def("problem","four_pi_G_par",1.0);
#endif

  /* particle number */
//  if (npartypes != 1)
//    ath_error("[particle_sphere]: This test only allows ONE particle species!\n");

  Npar  = (int)(pow(par_geti("particle","parnumcell"),1.0/3.0));
  Npar3 = Npar*SQR(Npar);

  pGrid->nparticle         = 2;
  grproperty[0].num = 1;
  grproperty[0].m = 1.;


  /* set the down sampling of the particle list output
   * (by default, output 1 particle per cell)
   */
  downsamp = 0;

  /* particle stopping time */
  tstop0[0] = par_getd("problem","tstop"); /* in code unit */
  if (par_geti("particle","tsmode") != 3)
    ath_error("[particle_sphere]: This test works only for fixed stopping time!\n");


/* Now set initial conditions for the gas */
  for (k=pGrid->ks; k<=pGrid->ke; k++) {
  for (j=pGrid->js; j<=pGrid->je; j++) {
  for (i=pGrid->is; i<=pGrid->ie; i++) {
    pGrid->U[k][j][i].d = dgas;
    pGrid->U[k][j][i].M1 = 0.;
    pGrid->U[k][j][i].M2 = 0.;
    pGrid->U[k][j][i].M3 = 0.;
  }}}

  pGrid->particle[0].property = 0;
  pGrid->particle[1].property = 0;
  pGrid->particle[0].pos = 1;
  pGrid->particle[1].pos = 1;
  pGrid->particle[0].my_id = 0;
  pGrid->particle[1].my_id = 1;
  pGrid->particle[0].x1 = -0.5;
  pGrid->particle[0].x2 =  0.0;
  pGrid->particle[0].x3 =  0.0;
  pGrid->particle[1].x1 =  0.5;
  pGrid->particle[1].x2 =  0.0;
  pGrid->particle[1].x3 =  0.0;
  pGrid->particle[0].v1 =  0.0;
  pGrid->particle[0].v2 =  0.0;
  pGrid->particle[0].v3 =  0.0;
  pGrid->particle[1].v1 =  0.0;
  pGrid->particle[1].v2 =  0.0;
  pGrid->particle[1].v3 =  0.0;
 
  return;
}

/*==============================================================================
 * PROBLEM USER FUNCTIONS:
 * problem_write_restart() - writes problem-specific user data to restart files
 * problem_read_restart()  - reads problem-specific user data from restart files
 * get_usr_expr()          - sets pointer to expression for special output data
 * get_usr_out_fun()       - returns a user defined output function pointer
 * Userwork_in_loop        - problem specific work IN     main loop
 * Userwork_after_loop     - problem specific work AFTER  main loop
 *----------------------------------------------------------------------------*/

void problem_write_restart(MeshS *pM, FILE *fp)
{
  return;
}

void problem_read_restart(MeshS *pM, FILE *fp)
{
  Nx = pM->Nx[0]; /* used for particle selection function */
  Npar  = (int)(sqrt(par_geti("particle","parnumcell")));
  Npar3 = SQR(Npar)*Npar;

  return;
}

ConsFun_t get_usr_expr(const char *expr)
{
  return NULL;
}

VOutFun_t get_usr_out_fun(const char *name)
{
  return NULL;
}

/*! \fn void gasvshift(const Real x1, const Real x2, const Real x3,
 *                                  Real *u1, Real *u2, Real *u3)
 *  \brief Gas velocity shift*/
void gasvshift(const Real x1, const Real x2, const Real x3,
                                    Real *u1, Real *u2, Real *u3)
{
  return;
}

#ifdef PARTICLES
PropFun_t get_usr_par_prop(const char *name)
{
  if (strcmp(name,"downsamp")==0) return property_mybin;
  return NULL;
}

void Userforce_particle(Real3Vect *ft, const Real x1, const Real x2, const Real x3,
                                    const Real v1, const Real v2, const Real v3)
{
  return;
}
#endif

void Userwork_in_loop(MeshS *pM)
{
  return;
}

/*---------------------------------------------------------------------------
 * Userwork_after_loop: computes L1-error in linear waves,
 * ASSUMING WAVE HAS PROPAGATED AN INTEGER NUMBER OF PERIODS
 * Must set parameters in input file appropriately so that this is true
 */

void Userwork_after_loop(MeshS *pM)
{
  return;
}
 
static int property_mybin(const GrainS *gr, const GrainAux *grsub)
{
  long a,b,c,d,e,ds,sp;

  sp = MAX(downsamp/Npar3,1);  /* spacing in cells */
  ds = Npar3*sp;               /* actual dowmsampling */

  a = gr->my_id/ds;
  b = gr->my_id - a*ds;

  c = gr->my_id/(Npar3*Nx);    /* column number */
  d = c/sp;
  e = c-sp*d;

  if ((e == 0) && (b == 0) && (gr->pos == 1))
    return 1;
  else
    return 0;
}

/*------------------------------------------------------------------------------
 */

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define RNMX (1.0-DBL_EPSILON)

/*! \fn double ran2(long int *idum)
 *  \brief Extracted from the Numerical Recipes in C (version 2) code.  Modified
 *   to use doubles instead of floats. -- T. A. Gardiner -- Aug. 12, 2003 
 *
 * Long period (> 2 x 10^{18}) random number generator of L'Ecuyer
 * with Bays-Durham shuffle and added safeguards.  Returns a uniform
 * random deviate between 0.0 and 1.0 (exclusive of the endpoint
 * values).  Call with idum = a negative integer to initialize;
 * thereafter, do not alter idum between successive deviates in a
 * sequence.  RNMX should appriximate the largest floating point value
 * that is less than 1.
 */

double ran2(long int *idum)
{
  int j;
  long int k;
  static long int idum2=123456789;
  static long int iy=0;
  static long int iv[NTAB];
  double temp;

  if (*idum <= 0) { /* Initialize */
    if (-(*idum) < 1) *idum=1; /* Be sure to prevent idum = 0 */
    else *idum = -(*idum);
    idum2=(*idum);
    for (j=NTAB+7;j>=0;j--) { /* Load the shuffle table (after 8 warm-ups) */
      k=(*idum)/IQ1;
      *idum=IA1*(*idum-k*IQ1)-k*IR1;
      if (*idum < 0) *idum += IM1;
      if (j < NTAB) iv[j] = *idum;
    }
    iy=iv[0];
  }
  k=(*idum)/IQ1;                 /* Start here when not initializing */
  *idum=IA1*(*idum-k*IQ1)-k*IR1; /* Compute idum=(IA1*idum) % IM1 without */
  if (*idum < 0) *idum += IM1;   /* overflows by Schrage's method */
  k=idum2/IQ2;
  idum2=IA2*(idum2-k*IQ2)-k*IR2; /* Compute idum2=(IA2*idum) % IM2 likewise */
  if (idum2 < 0) idum2 += IM2;
  j=(int)(iy/NDIV);              /* Will be in the range 0...NTAB-1 */
  iy=iv[j]-idum2;                /* Here idum is shuffled, idum and idum2 */
  iv[j] = *idum;                 /* are combined to generate output */
  if (iy < 1) iy += IMM1;
  if ((temp=AM*iy) > RNMX) return RNMX; /* No endpoint values */
  else return temp;
}

#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV
#undef RNMX
