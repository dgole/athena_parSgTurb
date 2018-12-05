#include "copyright.h"
/*============================================================================*/
/*! \file par_strat3d.c
 *  \brief Problem generator for non-linear streaming instability in stratified
 *   disks.
 *
 * PURPOSE: Problem generator for non-linear streaming instability in stratified
 *   disks. This code works in 3D ONLY. Isothermal eos is assumed, and the value
 *   etavk/iso_sound is fixed. MPI domain decomposition in x is allowed, but
 *   not in z.
 *
 * Perturbation modes:
 * -  ipert = 0: multi-nsh equilibtium
 * -  ipert = 1: white noise within the entire grid
 * -  ipert = 2: non-nsh velocity
 *
 *  Should be configured using --enable-shearing-box and --with-eos=isothermal.
 *  FARGO is recommended.
 *
 * Reference:
 * - Johansen & Youdin, 2007, ApJ, 662, 627
 * - Bai & Stone, 2009, in preparation */
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
#include "particles/prototypes.h"

#ifndef SHEARING_BOX
#error : The streaming3d problem requires shearing-box to be enabled.
#endif /* SHEARING_BOX */

#ifndef PARTICLES
#error : The streaming3d problem requires particles to be enabled.
#endif /* PARTICLES */

#ifndef ISOTHERMAL
#error : The streaming3d problem requires isothermal equation of state.
#endif /* ISOTHERMAL */

/* Uncomment the following define to drive the flow in an impulsive manner
   as was done originally.  Restarts for this mode not yet implemented! */
#define IMPULSIVE_DRIVING

/* KEEP SEMI-COLONS OUT OF THESE PRE-PROCESSOR DIRECTIVES! */
/* FFT indexing Nfast=k, Nmid=j, Nslow=i (opposite to Athena)
 * For OFST, i,j,k,nx2,nx3 reference the local grid */
#define OFST(i, j, k) ((k) + nx3*((j) + nx2*(i)))
/* KWVM: magnitude of wavenumber k in units of dkx */
#define KWVM(i, j, k) (sqrt(SQR(KCOMP(i,gis,gnx1))+ \
                            SQR(KCOMP(j,gjs,gnx2))+SQR(KCOMP(k,gks,gnx3))))

/* FFTW - Variables, Plan, etc. */
/* These are made static global variables so that they need not be
   allocated AND destroyed with each call to pspect! */
static struct ath_3d_fft_plan *plan;
/* Between calls to generate(), these have unshifted, unnormalized
 * velocity perturbations. */
static ath_fft_data *fv1=NULL, *fv2=NULL, *fv3=NULL;

/* Amplitude to stir velocities and input alpha */
static Real amp_force, alpha_in;
/* Entire size of mesh in all three dimensions */
Real Lx_all,Ly_all,Lz_all;
/* Forcing flag */
static int kforce;
/* Normalized, shifted velocity perturbations */
static Real ***dv1=NULL, ***dv2=NULL, ***dv3=NULL;
/* Cutoff wavenumbers, G&O spect peak, power law spect exponent, 2 pi/L */
static Real klow,khigh,kpeak,expo,dkx;
/* Energy injection rate, last planned driving time, driving interval.
 * If not using impulsive driving, then the time quantities above are for
 * computing a new spectrum, not driving */
static Real dedt,tdrive,dtdrive;
/* Driving properties */
static int ispect,idrive;
/* Number of cells in local grid, number of cells in global grid */
static int nx1,nx2,nx3,gnx1,gnx2,gnx3;
/* Starting and ending indices for global grid */
static int gis,gie,gjs,gje,gks,gke;
/* Seed for random number generator */
long int rseed;
/* Initial density (will be average density throughout simulation) */
static const Real rhobar = 1.0;

/* Functions appear in this file in the same order that they appear in the
 * prototypes below */

/* Function prototypes for generating velocity perturbations */
static void pspect(ath_fft_data *ampl);
static void project();
static inline void transform();
static inline void generate();
static void perturb(GridS *pGrid, Real dt);

Real vsc1,vsc2;
int ipert;
/* domain size variables */
Real x1min,x1max,x2min,x2max,x3min,x3max,Lx,Ly,Lz,Lg;
long Npar;
/* output variables */
long ntrack;   /* number of tracer particles */
long nlis;     /* number of output particles for list output */
int mytype;    /* specific particle type for output */
Real dpar_thresh; /* threshold particle density */

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 * ran2()            - random number generator
 * Normal()          - normal distribution generator
 * Erf()             - error function
 * UnstratifiedDisk() - tidal potential in 3D shearing box
 * VertGrav()         - potential for vertical component of gravity
 * MultiNSH()        - multiple component NSH equilibrium solver
 * ShearingBoxPot()  - shearing box tidal gravitational potential
 * strat_ix3         - vertical outflow boundary for bottom of grid
 * strat_ox3         - vertical outflow boundary for top of grid
 * hst_rho_Vx_dVy()  - total Reynolds stress for history dump
 * expr_*()          - computes new output variables
 * hst_*             - adds new history variables
 * property_???()    - particle property selection function
  * output_1d()      - dumps horizontally averaged quantities to a text file
 *============================================================================*/

double ran2(long int *idum);
double Normal(long int *idum);
Real Erf(Real z);

void MultiNSH(int n, Real *tstop, Real *mratio, Real etavk,
                     Real *uxNSH, Real *uyNSH, Real *wxNSH, Real *wyNSH);
static Real hst_rho_Vx_dVy(const GridS *pG,const int i,const int j,const int k);
static Real UnstratifiedDisk(const Real x1, const Real x2, const Real x3);
static Real VertGrav(const Real x1, const Real x2, const Real x3);
static void strat_ix3(GridS *pG);
static void strat_ox3(GridS *pG);
static int property_limit(const GrainS *gr, const GrainAux *grsub);
static int property_trace(const GrainS *gr, const GrainAux *grsub);
static int property_type(const GrainS *gr, const GrainAux *grsub);
static int property_dense(const GrainS *gr, const GrainAux *grsub);
static Real expr_KE(const GridS *pG, const int i, const int j, const int k);
static Real mass_cons(MeshS *pM);
#ifdef ADIABATIC
static Real hst_E_total(const GridS *pG, const int i, const int j, const int k);
#endif
static void output_1d(MeshS *pM, OutputS *pOut);

/* Function prototypes for initializing and interfacing with Athena */
static void initialize(GridS *pGrid, DomainS *pD);

extern Real expr_dpar(const GridS *pG, const int i, const int j, const int k);
extern Real expr_V1par(const GridS *pG, const int i, const int j, const int k);
extern Real expr_V2par(const GridS *pG, const int i, const int j, const int k);
extern Real expr_V3par(const GridS *pG, const int i, const int j, const int k);
extern Real expr_V2(const GridS *pG, const int i, const int j, const int k);
/* top and bottom of root Domain, shared with outputs, etc. */
static Real ztop, zbtm;
/* Total mass in grid */
static Real mass;
/* Apply a density floor - useful for large |z| regions */
static Real D_FLOOR = 1.e-4;
/* Flag to determine whether or not to employ outflow boundaries
   in z
   zbc_out = 1 is outflow
   zbc_out = 0 is periodic or reflecting
*/
static int zbc_out = 0;

/*=========================== PUBLIC FUNCTIONS =================================
 *============================================================================*/
/*----------------------------------------------------------------------------*/
/* problem:   */

void problem(DomainS *pDomain)
{
  GridS *pGrid = pDomain->Grid;
  int i,j,k,ks,pt,tsmode;
  int ixs,jxs,kxs;
  long p,q;
  Real ScaleHg,tsmin,tsmax,tscrit,amin,amax,Hparmin,Hparmax;
  Real *ep,*ScaleHpar,epsum,mratio,pwind,rhoaconv,etavk;
  Real *epsilon,*uxNSH,*uyNSH,**wxNSH,**wyNSH;
  Real rhog,h,x1,x2,x3,t,x1p,x2p,x3p,zmin,zmax,dx3_1,b;
  Real dVol,mtot=0.0,my_mtot;
#ifdef MPI_PARALLEL
  long int iseed = myID_Comm_world; /* Initialize on the first call to ran2 */
  int ierr;
#endif
#ifdef PARTICLE_SELF_GRAVITY
  four_pi_G_par = par_getd_def("problem","four_pi_G_par",1.0);
#endif

  if (pDomain->Nx[2] == 1) {
    ath_error("[par_strat3d]: par_strat3d only works for 3D problem.\n");
  }

#ifdef MPI_PARALLEL
  if (pDomain->NGrid[2] > 2) {
    ath_error(
  "[par_strat3d]: The z-domain can not be decomposed into more than 2 grids\n");
  }
#endif

/* Ensure a different initial random seed for each process in an MPI calc. */
  ixs = pGrid->Disp[0];
  jxs = pGrid->Disp[1];
  kxs = pGrid->Disp[2];
  rseed = -1 - (ixs + pDomain->Nx[0]*(jxs + pDomain->Nx[1]*kxs));
  initialize(pGrid, pDomain);
  tdrive = 0.0;

/* Initialize boxsize */
  x1min = pGrid->MinX[0];
  x1max = pGrid->MaxX[0];
  Lx = x1max - x1min;

  x2min = pGrid->MinX[1];
  x2max = pGrid->MaxX[1];
  Ly = x2max - x2min;

  x3min = par_getd("domain1","x3min");
  x3max = par_getd("domain1","x3max");
  Lz = x3max - x3min;

  Lx_all = par_getd("domain1","x1max") - par_getd("domain1","x1min");
  Ly_all = par_getd("domain1","x2max") - par_getd("domain1","x2min");
  Lz_all = par_getd("domain1","x3max") - par_getd("domain1","x3min");

  Lg = nghost*pGrid->dx3; /* size of the ghost zone */

  ks = pGrid->ks;

/* Initialize boxsize */
  ztop = x3max;
  zbtm = x3min;

/* Read initial conditions and parameters */
  Omega_0 = par_getd("problem","omega");
  qshear = par_getd_def("problem","qshear",1.5);
  ipert = par_geti_def("problem","ipert",1);
  vsc1 = par_getd_def("problem","vsc1",0.05); /* in unit of iso_sound (N.B.!) */
  vsc2 = par_getd_def("problem","vsc2",0.0);

  alpha_in = par_getd_def("problem","alpha_in",1.e-4);
  kforce   = par_geti_def("problem","kforce",0);

  vsc1 = vsc1 * Iso_csound;
  vsc2 = vsc2 * Iso_csound;

  ScaleHg = Iso_csound/Omega_0;

  /* particle number */
  Npar  = (long)(par_geti("particle","parnumgrid"));

  pGrid->nparticle = Npar*npartypes;
  for (i=0; i<npartypes; i++)
    grproperty[i].num = Npar;

  if (pGrid->nparticle+2 > pGrid->arrsize)
    particle_realloc(pGrid, pGrid->nparticle+2);

  ep = (Real*)calloc_1d_array(npartypes, sizeof(Real));
  ScaleHpar = (Real*)calloc_1d_array(npartypes, sizeof(Real));

  epsilon = (Real*)calloc_1d_array(npartypes, sizeof(Real));
  wxNSH   = (Real**)calloc_2d_array(pGrid->Nx[2]+1, npartypes,sizeof(Real));
  wyNSH   = (Real**)calloc_2d_array(pGrid->Nx[2]+1, npartypes,sizeof(Real));
  uxNSH   = (Real*)calloc_1d_array(pGrid->Nx[2]+1, sizeof(Real));
  uyNSH   = (Real*)calloc_1d_array(pGrid->Nx[2]+1, sizeof(Real));

  /* particle stopping time */
  tsmode = par_geti("particle","tsmode");
  if (tsmode == 3) {/* fixed stopping time */
    tsmin = par_getd("problem","tsmin"); /* in code unit */
    tsmax = par_getd("problem","tsmax");
    tscrit= par_getd("problem","tscrit");

    for (i=0; i<npartypes; i++) {
      tstop0[i] = tsmin*exp(i*log(tsmax/tsmin)/MAX(npartypes-1,1.0));
      grproperty[i].rad = tstop0[i];
      /* use fully implicit integrator for well coupled particles */
      if (tstop0[i] < tscrit) grproperty[i].integrator = 3;
    }
  }
  else {
    amin = par_getd("problem","amin");
    amax = par_getd("problem","amax");

    for (i=0; i<npartypes; i++)
      grproperty[i].rad = amin*exp(i*log(amax/amin)/MAX(npartypes-1,1.0));

    if (tsmode <= 2) {/* Epstein/General regime */
      /* conversion factor for rhoa */
      rhoaconv = par_getd_def("problem","rhoaconv",1.0);

      for (i=0; i<npartypes; i++)
        grrhoa[i]=grproperty[i].rad*rhoaconv;
    }

    if (tsmode == 1)  /* General drag formula */
      alamcoeff = par_getd("problem","alamcoeff");
  }

  /* particle scale height */
  Hparmax = par_getd("problem","hparmax"); /* in unit of gas scale height */
  Hparmin = par_getd("problem","hparmin");
  for (i=0; i<npartypes; i++)
    ScaleHpar[i] = Hparmax*
                   exp(-i*log(Hparmax/Hparmin)/MAX(npartypes-1,1.0));

#ifdef FEEDBACK
  mratio = par_getd_def("problem","mratio",0.0); /* total mass fraction */
  pwind = par_getd_def("problem","pwind",0.0);   /* power law index */
  if (mratio < 0.0)
    ath_error("[par_strat2d]: mratio must be positive!\n");

  epsum = 0.0;
  for (i=0; i<npartypes; i++)
  {
    ep[i] = pow(grproperty[i].rad,pwind);	epsum += ep[i];
  }

  for (i=0; i<npartypes; i++)
  {
    ep[i] = mratio*ep[i]/epsum;
    grproperty[i].m = sqrt(2.0*PI)*ScaleHg/Lz*ep[i]*
                                   pGrid->Nx[0]*pGrid->Nx[1]*pGrid->Nx[2]/Npar;
  }
#else
  mratio = 0.0;
  for (i=0; i<npartypes; i++)
    ep[i] = 0.0;
#endif

  /* NSH equilibrium */
  for (k=pGrid->ks; k<=pGrid->ke+1; k++) {

    h = pGrid->MinX[2] + (k-pGrid->ks)*pGrid->dx3;
    q = k - ks;
    etavk = fabs(vsc1+vsc2*SQR(h));

    for (i=0; i<npartypes; i++) {
      epsilon[i] = ep[i]/ScaleHpar[i]*exp(-0.5*SQR(h/ScaleHg)
         *(SQR(1.0/ScaleHpar[i])-1.0))/erf(Lz/(sqrt(8.0)*ScaleHpar[i]*ScaleHg));

      if (tsmode != 3)
        tstop0[i] = get_ts(pGrid,i,exp(-0.5*SQR(h/ScaleHg)),Iso_csound,etavk);
    }

    MultiNSH(npartypes, tstop0, epsilon, etavk,
                              &uxNSH[q], &uyNSH[q], wxNSH[q], wyNSH[q]);
  }

/* Now set initial conditions for the gas */
  for (k=pGrid->ks; k<=pGrid->ke; k++) {
  for (j=pGrid->js; j<=pGrid->je; j++) {
  for (i=pGrid->is; i<=pGrid->ie; i++) {
    cc_pos(pGrid,i,j,k,&x1,&x2,&x3);

    rhog = exp(-0.5*SQR(x3/ScaleHg));
    pGrid->U[k][j][i].d = rhog;

    if (ipert != 1) {/* NSH velocity */
      pGrid->U[k][j][i].M1 = 0.5*rhog*(uxNSH[k-ks]+uxNSH[k-ks+1]);
      pGrid->U[k][j][i].M2 = 0.5*rhog*(uyNSH[k-ks]+uyNSH[k-ks+1]);
    } else {
      pGrid->U[k][j][i].M1 = 0.0;
      pGrid->U[k][j][i].M2 = 0.0;
    }

    pGrid->U[k][j][i].M3 = 0.0;
#ifndef FARGO
    pGrid->U[k][j][i].M2 -= qshear*rhog*Omega_0*x1;
#endif

  }}}

/* Now set initial conditions for the particles */
  p = 0;
  dx3_1 = 1.0/pGrid->dx3;
  zmin = pGrid->MinX[2];
  zmax = pGrid->MaxX[2];

  for (q=0; q<Npar; q++) {

    for (pt=0; pt<npartypes; pt++) {

      x1p = x1min + Lx*ran2(&iseed);
      x2p = x2min + Ly*ran2(&iseed);
      x3p = ScaleHpar[pt]*ScaleHg*Normal(&iseed);
      while ((x3p >= zmax) || (x3p < zmin))
        x3p = ScaleHpar[pt]*ScaleHg*Normal(&iseed);

      pGrid->particle[p].property = pt;
      pGrid->particle[p].x1 = x1p;
      pGrid->particle[p].x2 = x2p;
      pGrid->particle[p].x3 = x3p;

      if (ipert != 1) {/* NSH velocity */

        cellk(pGrid, x3p, dx3_1, &k, &b);
        k = k-pGrid->ks;  b = b - pGrid->ks;

        pGrid->particle[p].v1 = (k+1-b)*wxNSH[k][pt]+(b-k)*wxNSH[k+1][pt];
        pGrid->particle[p].v2 = (k+1-b)*wyNSH[k][pt]+(b-k)*wyNSH[k+1][pt];

      } else {

        pGrid->particle[p].v1 = 0.0;
        pGrid->particle[p].v2 = vsc1+vsc2*SQR(x2p);

      }

      pGrid->particle[p].v3 = 0.0;
#ifndef FARGO
      pGrid->particle[p].v2 -= qshear*Omega_0*x1p;
#endif

      pGrid->particle[p].pos = 1; /* grid particle */
      pGrid->particle[p].my_id = p;
#ifdef MPI_PARALLEL
      pGrid->particle[p].init_id = myID_Comm_world;
#endif
      p++;
  }}

/* enroll gravitational potential function */

  StaticGravPot = VertGrav;
  ShearingBoxPot = UnstratifiedDisk;

/* Enroll vertically stratified outflow boundaries */
  if (zbc_out == 1) {
    set_bvals_particle_fun(left_x3, strat_ix3);
    set_bvals_particle_fun(right_x3, strat_ox3);
  }

  dump_history_enroll(hst_rho_Vx_dVy, "<rho Vx dVy>");

  /* set the # of the particles in list output
   * (by default, output 1 particle per cell)
   */
  nlis = par_geti_def("problem","nlis",pGrid->Nx[0]*pGrid->Nx[1]*pGrid->Nx[2]);

  /* set the number of particles to keep track of */
  ntrack = par_geti_def("problem","ntrack",2000);

  /* set the threshold particle density */
  dpar_thresh = par_geti_def("problem","dpar_thresh",10.0);

/* Calculate total mass on mesh */
  dVol = 1.0;
  if (pGrid->dx1 > 0.0) dVol *= pGrid->dx1;
  if (pGrid->dx2 > 0.0) dVol *= pGrid->dx2;
  if (pGrid->dx3 > 0.0) dVol *= pGrid->dx3;

  for (k=pGrid->ks; k<=pGrid->ke; k++) {
  for (j=pGrid->js; j<=pGrid->je; j++) {
  for (i=pGrid->is; i<=pGrid->ie; i++) {
    mtot += pGrid->U[k][j][i].d;
  }}}
    mtot *= dVol;

#ifdef MPI_PARALLEL
    ierr = MPI_Reduce(&(mtot), &(my_mtot), 1,
           MPI_DOUBLE, MPI_SUM, 0, pDomain->Comm_Domain);
#else
    my_mtot = mtot;
#endif

    if(myID_Comm_world==0){  /* I'm the parent */

      dVol = (pDomain->MaxX[0]-pDomain->MinX[0])
            *(pDomain->MaxX[1]-pDomain->MinX[1])*(pDomain->MaxX[2]-pDomain->MinX[2]);
      my_mtot /= dVol;

      /* calculate the mass conversion ratio */
      mass = my_mtot;
    }

#ifdef MPI_PARALLEL
  ierr = MPI_Bcast(&(mass), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif

  /* Set the initial perturbations.  Note that we're putting in too much
   * energy this time.  This is okay since we're only interested in the
   * saturated state. */
//  generate();
//  perturb(pGrid, dtdrive);

  /* finalize */
  free(ep);  free(ScaleHpar);
  free(epsilon);
  free_2d_array(wxNSH);  free_2d_array(wyNSH);
  free(uxNSH);           free(uyNSH);

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
  DomainS *pD = (DomainS*)&(pM->Domain[0][0]);
  GridS *pG = pD->Grid;
  Real x1,x2,x3;
  Real dVol,mtot=0.0,my_mtot;
  int i,j,k;
#ifdef MPI_PARALLEL
  int ierr;
#endif
#ifdef PARTICLE_SELF_GRAVITY
  four_pi_G_par = par_getd_def("problem","four_pi_G_par",1.0);
#endif

  Omega_0 = par_getd("problem","omega");
  qshear = par_getd_def("problem","qshear",1.5);
  ipert = par_geti_def("problem","ipert",1);

  x1min = pG->MinX[0];
  x1max = pG->MaxX[0];
  Lx = x1max - x1min;

  x2min = pG->MinX[1];
  x2max = pG->MaxX[1];
  Ly = x2max - x2min;

  x3min = pM->RootMinX[2];
  x3max = pM->RootMaxX[2];
  Lz = x3max - x3min;

  Lx_all = par_getd("domain1","x1max") - par_getd("domain1","x1min");
  Ly_all = par_getd("domain1","x2max") - par_getd("domain1","x2min");
  Lz_all = par_getd("domain1","x3max") - par_getd("domain1","x3min");

/* Initialize boxsize */
  ztop = x3max;
  zbtm = x3min;

  Lg = nghost*pG->dx3; /* size of the ghost zone */

  vsc1 = par_getd_def("problem","vsc1",0.05); /* in unit of iso_sound (N.B.!) */
  vsc2 = par_getd_def("problem","vsc2",0.0);

  vsc1 = vsc1 * Iso_csound;
  vsc2 = vsc2 * Iso_csound;

  alpha_in = par_getd_def("problem","alpha_in",1.e-4);
  kforce   = par_geti_def("problem","kforce",0);

  Npar  = (int)(sqrt(par_geti("particle","parnumgrid")));
  nlis = par_geti_def("problem","nlis",pG->Nx[0]*pG->Nx[1]*pG->Nx[2]);
  ntrack = par_geti_def("problem","ntrack",2000);

/* enroll gravitational potential function */

  StaticGravPot = VertGrav;
  ShearingBoxPot = UnstratifiedDisk;

  dump_history_enroll(hst_rho_Vx_dVy, "<rho Vx dVy>");

/* If using outflow boundaries, have to enroll them here too */

  if (zbc_out == 1) {
    set_bvals_particle_fun(left_x3, strat_ix3);
    set_bvals_particle_fun(right_x3, strat_ox3);
  }

/* Calculate total mass on mesh */
  dVol = 1.0;
  if (pG->dx1 > 0.0) dVol *= pG->dx1;
  if (pG->dx2 > 0.0) dVol *= pG->dx2;
  if (pG->dx3 > 0.0) dVol *= pG->dx3;

  for (k=pG->ks; k<=pG->ke; k++) {
  for (j=pG->js; j<=pG->je; j++) {
  for (i=pG->is; i<=pG->ie; i++) {
    mtot += pG->U[k][j][i].d;
  }}}
    mtot *= dVol;

#ifdef MPI_PARALLEL
    ierr = MPI_Reduce(&(mtot), &(my_mtot), 1,
           MPI_DOUBLE, MPI_SUM, 0, pD->Comm_Domain);
#else
    my_mtot = mtot;
#endif

    if(myID_Comm_world==0){  /* I'm the parent */

      dVol = (pD->MaxX[0]-pD->MinX[0])
            *(pD->MaxX[1]-pD->MinX[1])*(pD->MaxX[2]-pD->MinX[2]);
      my_mtot /= dVol;

      /* calculate the mass conversion ratio */
      mass = my_mtot;
    }

#ifdef MPI_PARALLEL
  ierr = MPI_Bcast(&(mass), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif

   return;
}

/* Enforce mass conservation */
static Real mass_cons(MeshS *pM)
{
  GridS *pG;
  DomainS *pD;
  int i,j,k,ierr;
  Real dVol,mtot=0.0,my_mtot,mratio;

  /* Loop over the root domain */
  if (pM->Domain[0][0].Grid != NULL){
    pG = pM->Domain[0][0].Grid;
    pD = (DomainS*)&(pM->Domain[0][0]);

    dVol = 1.0;
    if (pG->dx1 > 0.0) dVol *= pG->dx1;
    if (pG->dx2 > 0.0) dVol *= pG->dx2;
    if (pG->dx3 > 0.0) dVol *= pG->dx3;

    for (k=pG->ks; k<=pG->ke; k++) {
    for (j=pG->js; j<=pG->je; j++) {
    for (i=pG->is; i<=pG->ie; i++) {
      mtot += pG->U[k][j][i].d;
    }}}
    mtot *= dVol;

#ifdef MPI_PARALLEL
    ierr = MPI_Reduce(&(mtot), &(my_mtot), 1,
           MPI_DOUBLE, MPI_SUM, 0, pD->Comm_Domain);
#else
    my_mtot = mtot;
#endif

    if(myID_Comm_world==0){  /* I'm the parent */

      dVol = (pD->MaxX[0]-pD->MinX[0])
            *(pD->MaxX[1]-pD->MinX[1])*(pD->MaxX[2]-pD->MinX[2]);
      my_mtot /= dVol;

      /* calculate the mass conversion ratio */
      mratio = mass / my_mtot;
    }
  }

#ifdef MPI_PARALLEL
  ierr = MPI_Bcast(&(mratio), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif

  return mratio;
}

static Real hst_rho_Vx_dVy(const GridS *pG,
                           const int i, const int j, const int k)
{
  Real x1,x2,x3;
  cc_pos(pG,i,j,k,&x1,&x2,&x3);
#ifdef FARGO
  return pG->U[k][j][i].M1*pG->U[k][j][i].M2/pG->U[k][j][i].d;
#else
  return pG->U[k][j][i].M1*(pG->U[k][j][i].M2/pG->U[k][j][i].d
                            + qshear*Omega_0*x1);
#endif
}

/*----------------------------------------------------------------------------*/
/*! \fn static Real expr_KE(const GridS *pG, const int i, const int j,
 *                          const int k)
 *  \brief Computes dens*(Vx^2+Vy^2+Vz^2)/2
 */
static Real expr_KE(const GridS *pG, const int i, const int j, const int k)
{
  Real x1,x2,x3,Vy,Vx,Vz;
  cc_pos(pG,i,j,k,&x1,&x2,&x3);
#ifdef FARGO
  Vy = (pG->U[k][j][i].M2/pG->U[k][j][i].d);
#else
  Vy = (pG->U[k][j][i].M2/pG->U[k][j][i].d + qshear*Omega_0*x1);
#endif
  Vx = pG->U[k][j][i].M1/pG->U[k][j][i].d;
  Vz = pG->U[k][j][i].M3/pG->U[k][j][i].d;

  return pG->U[k][j][i].d*(Vx*Vx + Vy*Vy + Vz*Vz)/2.0;

}

#ifdef ADIABATIC
/*! \fn static Real hst_E_total(const GridS *pG, const int i, const int j,
 *                              const int k)
 *  \brief total energy (including tidal potential). */
static Real hst_E_total(const GridS *pG, const int i, const int j, const int k)
{
  Real x1,x2,x3,phi;
  cc_pos(pG,i,j,k,&x1,&x2,&x3);
  phi = UnstratifiedDisk(x1, x2, x3);

  return pG->U[k][j][i].E + pG->U[k][j][i].d*phi;
}
#endif /* ADIABATIC */

/* Get_user_expression computes dVy */
ConsFun_t get_usr_expr(const char *expr)
{
  if(strcmp(expr,"KE")==0) return expr_KE;
  return NULL;
}

VOutFun_t get_usr_out_fun(const char *name)
{
  if(strcmp(name,"1d")==0) return output_1d;
  return NULL;
}


#ifdef PARTICLES
PropFun_t get_usr_par_prop(const char *name)
{
  if (strcmp(name,"limit")==0)    return property_limit;
  if (strcmp(name,"trace")==0)    return property_trace;
  if (strcmp(name,"type")==0)     return property_type;
  if (strcmp(name,"dense")==0)    return property_dense;

  return NULL;
}

void gasvshift(const Real x1, const Real x2, const Real x3,
                                    Real *u1, Real *u2, Real *u3)
{
  return;
}

void Userforce_particle(Real3Vect *ft, const Real x1, const Real x2, const Real x3,
                                    const Real v1, const Real v2, const Real v3)
{
  Real z,fac;

  ft->x1 -= 2.0*(vsc1 + vsc2*SQR(x3))*Omega_0;

  if (zbc_out == 1) {
    ft->x3 -= SQR(Omega_0)*x3;
  } else {
    if(x3 > x3max)
      z = x3-Lz;
    else if (x3 < x3min)
      z = x3+Lz;
    else
      z = x3;
    fac = Lg/(0.5*Lz+Lg-fabs(z));
    ft->x3 -= SQR(Omega_0)*z*(1.0-SQR(fac)*fac); /* 3rd order sharp */
//  ft->x3 -= SQR(Omega_0)*z*(1.0-SQR(fac));  /* 2nd order sharp */
  }

  return;
}
#endif

void Userwork_in_loop(MeshS *pM)
{
  GridS *pG;
  DomainS *pD;
  int nl,nd,i,j,k;
  Real mratio,newtime;
  Real x1,x2,x3,kx,ky,kz,vx,vy,vz;
  Real phixx,phixy,phixz,phiyx,phiyy,phiyz,phizx,phizy,phizz;
  long iseedxx,iseedxy,iseedxz,iseedyx,iseedyy,iseedyz,iseedzx,iseedzy,iseedzz;

/* Enforce mass conservation */
  mratio = mass_cons(pM);

  for (nl=0; nl<(pM->NLevels); nl++){
    for (nd=0; nd<(pM->DomainsPerLevel[nl]); nd++){
      if (pM->Domain[nl][nd].Grid != NULL){

        pG = pM->Domain[nl][nd].Grid;
        for (k=pG->ks; k<=pG->ke; k++) {
          for (j=pG->js; j<=pG->je; j++) {
            for (i=pG->is; i<=pG->ie; i++) {
              pG->U[k][j][i].d = MAX(pG->U[k][j][i].d * mratio, D_FLOOR);
            }}}

        if (isnan(pG->dt)) ath_error("Time step is NaN!");

        if (kforce == 1) {
          if (idrive == 0) {  /* driven turbulence */
            /* Integration has already been done, but time not yet updated */
            newtime = pG->time + pG->dt;

#ifndef IMPULSIVE_DRIVING
            /* Drive on every time step */
            perturb(pG, pG->dt);
#endif /* IMPULSIVE_DRIVING */

            if (newtime >= (tdrive+dtdrive)) {
              /* If we start with large time steps so that tdrive would get way
               * behind newtime, this makes sure we don't keep generating after
               * dropping down to smaller time steps */
              while ((tdrive+dtdrive) <= newtime) tdrive += dtdrive;

#ifdef IMPULSIVE_DRIVING
              /* Only drive at intervals of dtdrive */
              perturb(pG, dtdrive);
#endif /* IMPULSIVE_DRIVING */

              /* Compute new spectrum after dtdrive.  Putting this after perturb()
               * means we won't be applying perturbations from a new power spectrum
               * just before writing outputs.  At the very beginning, we'll go a
               * little longer before regenerating, but the energy injection rate
               * was off on the very first timestep anyway.  When studying driven
               * turbulence, all we care about is the saturated state. */
               generate();
            }
          }
        }

        if (kforce == 0) {

          // try driving in one fourth of the box
          //float factor = 4.0;
          //amp_force = factor*4.*sqrt(0.22*alpha_in*(1./Omega_0)/(Lx_all*Ly_all*Lz_all));

          amp_force = 4.*sqrt(0.22*alpha_in*(1./Omega_0)/(Lx_all*Ly_all*Lz_all));

          kx = 2.*PI/(Lx_all/factor);
          ky = 2.*PI/(Ly_all/factor);
          kz = 2.*PI/(Lz_all/factor);

          iseedxx = pM->nstep*1;
          iseedxy = pM->nstep*3;
          iseedxz = pM->nstep*5;

          iseedyx = pM->nstep*2;
          iseedyy = pM->nstep*4;
          iseedyz = pM->nstep*6;

          iseedzx = pM->nstep*9;
          iseedzy = pM->nstep*8;
          iseedzz = pM->nstep*7;

          phixx = ran2(&iseedxx)*2.*PI;
          phixy = ran2(&iseedxy)*2.*PI;
          phixz = ran2(&iseedxz)*2.*PI;

          phiyx = ran2(&iseedyx)*2.*PI;
          phiyy = ran2(&iseedyy)*2.*PI;
          phiyz = ran2(&iseedyz)*2.*PI;

          phizx = ran2(&iseedzx)*2.*PI;
          phizy = ran2(&iseedzy)*2.*PI;
          phizz = ran2(&iseedzz)*2.*PI;

          for (k=pG->ks; k<=pG->ke; k++) {
          for (j=pG->js; j<=pG->je; j++) {
          for (i=pG->is; i<=pG->ie; i++) {

            cc_pos(pG,i,j,k,&x1,&x2,&x3);

            /*if (fabs(x3)<(Lz_all/(factor*2.))) {
              vx = amp_force*(sin(kx*x1+phixx)*sin(ky*x2+phiyx)*sin(kz*x3+phizx));
              vy = amp_force*(sin(kx*x1+phixy)*sin(ky*x2+phiyy)*sin(kz*x3+phizy));
              vz = amp_force*(sin(kx*x1+phixz)*sin(ky*x2+phiyz)*sin(kz*x3+phizz));
              pG->U[k][j][i].M1 += pG->U[k][j][i].d*vx*pG->dt;
              pG->U[k][j][i].M2 += pG->U[k][j][i].d*vy*pG->dt;
              pG->U[k][j][i].M3 += pG->U[k][j][i].d*vz*pG->dt;
            }*/
            vx = amp_force*(sin(kx*x1+phixx)*sin(ky*x2+phiyx)*sin(kz*x3+phizx));
            vy = amp_force*(sin(kx*x1+phixy)*sin(ky*x2+phiyy)*sin(kz*x3+phizy));
            vz = amp_force*(sin(kx*x1+phixz)*sin(ky*x2+phiyz)*sin(kz*x3+phizz));
            pG->U[k][j][i].M1 += pG->U[k][j][i].d*vx*pG->dt;
            pG->U[k][j][i].M2 += pG->U[k][j][i].d*vy*pG->dt;
            pG->U[k][j][i].M3 += pG->U[k][j][i].d*vz*pG->dt;

          }}}
        }

  }}}

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


/*=========================== PRIVATE FUNCTIONS ==============================*/
/*--------------------------------------------------------------------------- */
/*----------------------------------------------------------------------------*/
/*! \fn static Real UnstratifiedDisk(const Real x1, const Real x2,const Real x3)
 *  \brief tidal potential in 3D shearing box */
static Real UnstratifiedDisk(const Real x1, const Real x2, const Real x3)
{
  Real phi=0.0;
#ifndef FARGO
  phi -= qshear*Omega_0*Omega_0*x1*x1;
#endif
  return phi;
}

/*----------------------------------------------------------------------------*/
/*! \fn static Real VertGrav(const Real x1, const Real x2, const Real x3)
 *  \brief potential for vertical component of gravity */
static Real VertGrav(const Real x1, const Real x2, const Real x3)
{
  Real phi=0.0,z,z0;

/* If outflow boundaries are used in z, we just use the normal
   z potential.  Otherwise, we ensure periodicity and also
   smooth the potential near the vertical boundaries */

  if (zbc_out == 1) {
    z = x3;
    phi += 0.5*Omega_0*Omega_0*z*z;
  } else {
    if(x3 > ztop)
      z=x3-ztop+zbtm;
    else if (x3 < zbtm)
      z=x3-zbtm+ztop;
    else
      z=fabs(x3);

      phi += 0.5*SQR(Omega_0*z);

      /* smooth the potential at domain edges */
      z0 = 0.5*Lz+Lg;
      phi -= SQR(Omega_0*Lg)*Lg*(2.0*z-z0)/(2.0*SQR(z0-z)); /* 3rd order sharp */
//  phi -= SQR(Omega_0*Lg)*(z/(z0-z)+log(z0/(z0-z))); /* 2nd order sharp */
  }
  return phi;
}

/*! \fn static void strat_ix3(GridS *pG)
 *  \brief  Here is the lower z outflow boundary.
            The basic idea is that the pressure and density
            are exponentially extrapolated in the ghost zones
            assuming a constant temperature there (i.e., an
            isothermal atmosphere). The z velocity (NOT the
            momentum) are set to zero in the ghost zones in the
            case of the last lower physical zone having an inward
            flow.  All other variables are extrapolated into the
            ghost zones with zero slope.
*/

static void strat_ix3(GridS *pG)
{
  int ks = pG->ks;
  int ie = pG->ie;
  int je = pG->je;
  int i,j,k,il,iu,jl,ju; /* i-lower/upper;  j-lower/upper */
  Real x1,x2,x3;
  Real press,pressks,csq;
  static Real x3b;

  x3b = zbtm+0.5*pG->dx3;

  if (pG->Nx[0] > 1){
    iu = pG->ie + nghost;
    il = pG->is - nghost;
  } else {
    iu = pG->ie;
    il = pG->is;
  }
  if (pG->Nx[1] > 1){
    ju = pG->je + nghost;
    jl = pG->js - nghost;
  } else {
    ju = pG->je;
    jl = pG->js;
  }

  for (k=1; k<=nghost; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        cc_pos(pG,i,j,ks-k,&x1,&x2,&x3);
/* First calculate the effective gas temperature in the last physical zone */
#ifdef ADIABATIC
        pressks = pG->U[ks][j][i].E - 0.5*(SQR(pG->U[ks][j][i].M1)
                      + SQR(pG->U[ks][j][i].M2)
                      + SQR(pG->U[ks][j][i].M3))/pG->U[ks][j][i].d;
        pressks *= Gamma_1;
        pressks = MAX(pressks,TINY_NUMBER);
        csq = pressks/pG->U[ks][j][i].d;
#else
        csq = Iso_csound2;
#endif /* ADIABATIC */
/* Now extrapolate the density to balance gravity assuming a constant temperature in the ghost zones */
        pG->U[ks-k][j][i].d = pG->U[ks][j][i].d*exp(-(x3*x3-x3b*x3b)/(2.0*csq/(Omega_0*Omega_0)));
/* Copy the velocities, but not the momenta --- important because of the density extrapolation above */
        pG->U[ks-k][j][i].M1 = pG->U[ks][j][i].M1/pG->U[ks][j][i].d*pG->U[ks-k][j][i].d;
        pG->U[ks-k][j][i].M2 = pG->U[ks][j][i].M2/pG->U[ks][j][i].d*pG->U[ks-k][j][i].d;
/* If there's inflow into the grid, set the normal velocity to zero */
        if (pG->U[ks][j][i].M3 >= 0.0) {
          pG->U[ks-k][j][i].M3 = 0.0;
        } else {
          pG->U[ks-k][j][i].M3 = pG->U[ks][j][i].M3/pG->U[ks][j][i].d*pG->U[ks-k][j][i].d;
        }
#ifdef ADIABATIC
        press = pG->U[ks-k][j][i].d*csq;
        pG->U[ks-k][j][i].E = press/Gamma_1
        + 0.5*(SQR(pG->U[ks-k][j][i].M1) + SQR(pG->U[ks-k][j][i].M2)
             + SQR(pG->U[ks-k][j][i].M3))/pG->U[ks-k][j][i].d;
#endif /* ADIABATIC */
      }
    }
  }

  return;

}

/*! \fn static void strat_ox3(GridS *pG)
 *  \brief  Here is the upper z outflow boundary.
            The basic idea is that the pressure and density
            are exponentially extrapolated in the ghost zones
            assuming a constant temperature there (i.e., an
            isothermal atmosphere). The z velocity (NOT the
            momentum) are set to zero in the ghost zones in the
            case of the last upper physical zone having an inward
            flow.  All other variables are extrapolated into the
            ghost zones with zero slope.
*/

static void strat_ox3(GridS *pG)
{
  int ke = pG->ke;
  int ie = pG->ie;
  int je = pG->je;
  int i,j,k,il,iu,jl,ju; /* i-lower/upper;  j-lower/upper */
  Real x1,x2,x3;
  Real press,presske,csq;
  static Real x3t;

  x3t = ztop-0.5*pG->dx3;

  if (pG->Nx[0] > 1){
    iu = pG->ie + nghost;
    il = pG->is - nghost;
  } else {
    iu = pG->ie;
    il = pG->is;
  }
  if (pG->Nx[1] > 1){
    ju = pG->je + nghost;
    jl = pG->js - nghost;
  } else {
    ju = pG->je;
    jl = pG->js;
  }

  for (k=1; k<=nghost; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        cc_pos(pG,i,j,ke+k,&x1,&x2,&x3);
#ifdef ADIABATIC
        presske = pG->U[ke][j][i].E - 0.5*(SQR(pG->U[ke][j][i].M1)
                      + SQR(pG->U[ke][j][i].M2)
                      + SQR(pG->U[ke][j][i].M3))/pG->U[ke][j][i].d;
                      + SQR(pG->U[ke][j][i].M2)
                      + SQR(pG->U[ke][j][i].M3))/pG->U[ke][j][i].d;
        presske *= Gamma_1;
        presske = MAX(presske,TINY_NUMBER);
        csq = presske/pG->U[ke][j][i].d;
#else
        csq = Iso_csound2;
#endif /* ADIABATIC */
/* Now extrapolate the density to balance gravity assuming a constant temperature in the ghost zones */
        pG->U[ke+k][j][i].d = pG->U[ke][j][i].d*exp(-(x3*x3-x3t*x3t)/(2.0*csq/(Omega_0*Omega_0)));
/* Copy the velocities, but not the momenta --- important because of the density extrapolation above */
        pG->U[ke+k][j][i].M1 = pG->U[ke][j][i].M1/pG->U[ke][j][i].d*pG->U[ke+k][j][i].d;
        pG->U[ke+k][j][i].M2 = pG->U[ke][j][i].M2/pG->U[ke][j][i].d*pG->U[ke+k][j][i].d;
/* If there's inflow into the grid, set the normal velocity to zero */
        if (pG->U[ke][j][i].M3 <= 0.0) {
          pG->U[ke+k][j][i].M3 = 0.0;
        } else {
          pG->U[ke+k][j][i].M3 = pG->U[ke][j][i].M3/pG->U[ke][j][i].d*pG->U[ke+k][j][i].d;
        }
#ifdef ADIABATIC
        press = pG->U[ke+k][j][i].d*csq;
        pG->U[ke+k][j][i].E = press/Gamma_1
        + 0.5*(SQR(pG->U[ke+k][j][i].M1) + SQR(pG->U[ke+k][j][i].M2)
             + SQR(pG->U[ke+k][j][i].M3))/pG->U[ke+k][j][i].d;
#endif /* ADIABATIC */
      }
    }
  }

  return;

}

/*! \fn static int property_limit(const GrainS *gr, const GrainAux *grsub)
 *  \brief user defined particle selection function (1: true; 0: false) */
static int property_limit(const GrainS *gr, const GrainAux *grsub)
{
  if ((gr->pos == 1) && (gr->my_id<nlis))
    return 1;
  else
    return 0;
}

/*! \fn static int property_trace(const GrainS *gr, const GrainAux *grsub)
 *  \brief user defined particle selection function (1: true; 0: false) */
static int property_trace(const GrainS *gr, const GrainAux *grsub)
{
  if ((gr->pos == 1) && (gr->my_id<ntrack))
    return 1;
  else
    return 0;
}

/*! \fn static int property_type(const GrainS *gr, const GrainAux *grsub)
 *  \brief user defined particle selection function (1: true; 0: false) */
static int property_type(const GrainS *gr, const GrainAux *grsub)
{
  if ((gr->pos == 1) && (gr->property == mytype))
    return 1;
  else
    return 0;
}

/*! \fn static int property_dense(const GrainS *gr, const GrainAux *grsub)
 *  \brief user defined particle selection function (1: true; 0: false) */
static int property_dense(const GrainS *gr, const GrainAux *grsub)
{
  if ((gr->pos == 1) && (grsub->dpar > dpar_thresh))
    return 1;
  else
    return 0;
}

/*! \fn static void output_1d(MeshS *pM, OutputS *pOut)
 *  \brief output routine to calculate 1D horizontally
    averaged quantities.  Currently, only outputs at lowest
    refinement level */

static void output_1d(MeshS *pM, OutputS *pOut)
{
  GridS *pGrid;
  DomainS *pD;
  int i,j,k;
  int tot1d,i1d,nzmx,my_nz,kg,kdisp;
  int dnum = pOut->num,nl,nd;
  static int FIRST = 0;
  double darea,**out1d;
  double x1,x2,x3,Lx,Ly,press;
  static double *out_x3;

  FILE *p_1dfile;
  char *fname;
  double area_rat; /* (Grid Volume)/(dx1*dx2*dx3) */

#ifdef MPI_PARALLEL
  double *my_out1d;
  double *g_out1d;
  int zproc;
  int ierr,myID_Comm_Domain;
#endif

  tot1d=7;
#ifdef ADIABATIC
  tot1d=tot1d+3;
#endif /* ADIABATIC */

  Lx = pM->RootMaxX[0] - pM->RootMinX[0];
  Ly = pM->RootMaxX[1] - pM->RootMinX[1];
  nzmx = pM->Nx[2];

/* At level=0, there is only one domain */

  pGrid = pM->Domain[0][0].Grid;
  int is = pGrid->is, ie = pGrid->ie;
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  pD = (DomainS*)&(pM->Domain[0][0]);

#ifdef MPI_PARALLEL
  int nproc = pD->NGrid[0]*pD->NGrid[1]*pD->NGrid[2];
#endif

#ifdef MPI_PARALLEL
  ierr = MPI_Comm_rank(pD->Comm_Domain, &myID_Comm_Domain);
  if(ierr != MPI_SUCCESS)
    ath_error("[change_rundir]: MPI_Comm_rank error = %d\n",ierr);
#endif
  if (FIRST == 0){
#ifdef MPI_PARALLEL
    if (myID_Comm_Domain == 0) {
#endif
      out_x3 = (double *) calloc_1d_array(nzmx,sizeof(double));
#ifdef MPI_PARALLEL
    }
#endif
  }

  out1d = (double **) calloc_2d_array(nzmx,tot1d,sizeof(double));
#ifdef MPI_PARALLEL
  my_out1d = (double *) calloc_1d_array(nzmx,sizeof(double));
  g_out1d = (double *) calloc_1d_array(nzmx,sizeof(double));
#endif
  for (k=0; k<nzmx; k++) {
    for (i1d=0; i1d<tot1d; i1d++) {
      out1d[k][i1d] = 0.0;
    }
  }
  kdisp=pGrid->Disp[2];

/* First calculate the x3 coordinate and save it to be dumped
   by root in every 1d file */
  if (FIRST == 0) {
#ifdef MPI_PARALLEL
  if (myID_Comm_Domain == 0) {
#endif
    for (k=0; k<nzmx; k++) {
      x3 = pM->RootMinX[2] + (k + 0.5)*pGrid->dx3;
      out_x3[k] = x3;
    }
#ifdef MPI_PARALLEL
  }
#endif
  }

/* Compute 1d averaged variables */
  for (k=ks; k<=ke; k++) {
    kg=k+kdisp-nghost;
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        i1d=0;
        out1d[kg][i1d] += pGrid->U[k][j][i].d;
        i1d++;
#ifdef ISOTHERMAL
        out1d[kg][i1d] += pGrid->U[k][j][i].d*Iso_csound2;
#else
        press           = MAX(Gamma_1*(pGrid->U[k][j][i].E - expr_KE(pGrid,i,j,k)
                                ),TINY_NUMBER);
        out1d[kg][i1d] += press;
#endif
#ifdef ADIABATIC
        i1d++;
        out1d[kg][i1d] += press/(Gamma_1*pGrid->U[k][j][i].d);
        i1d++;
        out1d[kg][i1d] += pGrid->U[k][j][i].E;
        i1d++;
        out1d[kg][i1d] += hst_E_total(pGrid,i,j,k);
#endif
        i1d++;
        out1d[kg][i1d] += 0.5*SQR(pGrid->U[k][j][i].M1)/pGrid->U[k][j][i].d;
        i1d++;
        cc_pos(pGrid,i,j,k,&x1,&x2,&x3);
#ifdef FARGO
        out1d[kg][i1d] += 0.5*SQR(pGrid->U[k][j][i].M2)/pGrid->U[k][j][i].d;
#else
        out1d[kg][i1d] += 0.5*pGrid->U[k][j][i].d*SQR(pGrid->U[k][j][i].M2/pGrid->U[k][j][i].d + qshear*Omega_0*x1);
#endif
        i1d++;
        out1d[kg][i1d] += 0.5*SQR(pGrid->U[k][j][i].M3)/pGrid->U[k][j][i].d;
        i1d++;
        out1d[kg][i1d] += expr_KE(pGrid,i,j,k);
        i1d++;
        out1d[kg][i1d] += hst_rho_Vx_dVy(pGrid,i,j,k);
      }
    }
  }

  /* Calculate the (Grid Volume) / (Grid Cell Volume) Ratio */
  area_rat = Lx*Ly/(pGrid->dx1*pGrid->dx2);

/* The parent sums the scal[] array.
 * Note that this assumes (dx1,dx2,dx3) = const. */

#ifdef MPI_PARALLEL
  for(i1d=0; i1d<tot1d; i1d++){
    for (k=0; k<nzmx; k++) {
      my_out1d[k] = out1d[k][i1d];
    }
    ierr = MPI_Reduce(my_out1d, g_out1d, nzmx,
                      MPI_DOUBLE, MPI_SUM, 0, pD->Comm_Domain);
    if(ierr)
      ath_error("[output_1d]: MPI_Reduce call returned error = %d\n",ierr);
    for (k=0; k<nzmx; k++) {
      out1d[k][i1d] = g_out1d[k];
    }
  }
#endif

/* For parallel calculations, only the parent computes the average
 * and writes the output. */
#ifdef MPI_PARALLEL
  if(myID_Comm_Domain == 0){ /* I'm the parent */
#endif

  darea = 1.0/(double)area_rat;
  for (k=0; k<nzmx; k++) {
    for (i1d=0; i1d<tot1d; i1d++) {
      out1d[k][i1d] *= darea;
    }
  }

/* Generate filename */
#ifdef MPI_PARALLEL
  fname = ath_fname("../",pM->outfilename,NULL,NULL,num_digit,dnum,NULL,"1d");
#else
  fname = ath_fname(NULL,pM->outfilename,NULL,NULL,num_digit,dnum,NULL,"1d");
#endif
  if (fname == NULL) {
    ath_error("[output_1d]: Error constructing output filename\n");
    return;
  }

/* open filename */
  p_1dfile = fopen(fname,"w");
  if (p_1dfile == NULL) {
    ath_error("[output_1d]: Unable to open 1d average file %s\n",fname);
    return;
  }

/* Write out data */

  for (k=0; k<nzmx; k++) {
#ifdef ISOTHERMAL
    if (k == 0) {
      fprintf(p_1dfile,"# x3     dens  pressure    KEx         KEy         KEz         KE          Reynolds\n");
    }
    fprintf(p_1dfile,"%G %G %G %G %G %G %G %G\n",out_x3[k],out1d[k][0],out1d[k][1],out1d[k][2],out1d[k][3],out1d[k][4],
            out1d[k][5],out1d[k][6]);
#else
    if (k == 0) {
      fprintf(p_1dfile,"# x3     dens    pressure    temperature  E     Etot     KEx         KEy         KEz         KE          Reynolds\n");
    }
    fprintf(p_1dfile,"%G %G %G %G %G %G %G %G %G %G %G\n",out_x3[k],out1d[k][0],out1d[k][1],out1d[k][2],out1d[k][3],out1d[k][4],
            out1d[k][5],out1d[k][6],out1d[k][7],out1d[k][8],out1d[k][9]);
#endif /* ISOTHERMAL */
  }

  fclose(p_1dfile);
  free(fname);
#ifdef MPI_PARALLEL
  }
#endif

  free_2d_array(out1d); /* Free the memory we malloc'd */
#ifdef MPI_PARALLEL
  free_1d_array(my_out1d); /* Free the memory we malloc'd */
  free_1d_array(g_out1d); /* Free the memory we malloc'd */
#endif
  if (FIRST == 0) {
    FIRST = 1;
  }

return;
}

/*----------------------------------------------------------------------------*/
/*! \fn void MultiNSH(int n, Real *tstop, Real *mratio, Real etavk,
 *                   Real *uxNSH, Real *uyNSH, Real *wxNSH, Real *wyNSH)
 *  \brief Multi-species NSH equilibrium
 *
 * Input: # of particle types (n), dust stopping time and mass ratio array, and
 *        drift speed etavk.
 * Output: gas NSH equlibrium velocity (u), and dust NSH equilibrium velocity
 *         array (w).
 */
void MultiNSH(int n, Real *tstop, Real *mratio, Real etavk,
                     Real *uxNSH, Real *uyNSH, Real *wxNSH, Real *wyNSH)
{
  int i,j;
  Real *Lambda1,**Lam1GamP1, **A, **B, **Tmp;

  Lambda1 = (Real*)calloc_1d_array(n, sizeof(Real));     /* Lambda^{-1} */
  Lam1GamP1=(Real**)calloc_2d_array(n, n, sizeof(Real)); /* Lambda1*(1+Gamma) */
  A       = (Real**)calloc_2d_array(n, n, sizeof(Real));
  B       = (Real**)calloc_2d_array(n, n, sizeof(Real));
  Tmp     = (Real**)calloc_2d_array(n, n, sizeof(Real));

  /* definitions */
  for (i=0; i<n; i++){
    for (j=0; j<n; j++)
      Lam1GamP1[i][j] = mratio[j];
    Lam1GamP1[i][i] += 1.0;
    Lambda1[i] = 1.0/(tstop[i]+1.0e-16);
    for (j=0; j<n; j++)
      Lam1GamP1[i][j] *= Lambda1[i];
  }

  /* Calculate A and B */
  MatrixMult(Lam1GamP1, Lam1GamP1, n,n,n, Tmp);
  for (i=0; i<n; i++) Tmp[i][i] += 1.0;
  InverseMatrix(Tmp, n, B);
  for (i=0; i<n; i++)
  for (j=0; j<n; j++)
    B[i][j] *= Lambda1[j];
  MatrixMult(Lam1GamP1, B, n,n,n, A);

  /* Obtain NSH velocities */
  *uxNSH = 0.0;  *uyNSH = 0.0;
  for (i=0; i<n; i++){
    wxNSH[i] = 0.0;
    wyNSH[i] = 0.0;
    for (j=0; j<n; j++){
      wxNSH[i] -= B[i][j];
      wyNSH[i] -= A[i][j];
    }
    wxNSH[i] *= 2.0*etavk;
    wyNSH[i] *= etavk;
    *uxNSH -= mratio[i]*wxNSH[i];
    *uyNSH -= mratio[i]*wyNSH[i];
    wyNSH[i] += etavk;
  }

  free(Lambda1);
  free_2d_array(A);         free_2d_array(B);
  free_2d_array(Lam1GamP1); free_2d_array(Tmp);

  return;
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

/*--------------------------------------------------------------------------- */
/*! \fn double Normal(long int *idum)
 *  \brief Normal distribution random number generator */
double Normal(long int *idum)
{
  double Y,X1,X2;

  X1 = ran2(idum);
  X2 = ran2(idum);

  Y = sqrt(-2.0*log(X1+TINY_NUMBER))*cos(2*PI*X2);

  return Y;
}



/* ========================================================================== */

/*! \fn static void pspect(ath_fft_data *ampl)
 *  \brief computes component of velocity with specific power
 *  spectrum in Fourier space determined by ispect
 *
 *  Velocity power spectrum returned in ampl
 *  - klow   = multiple of 2 pi/L for cut-off at low  wavenumbers
 *  - khigh  = multiple of 2 pi/L for cut-off at high wavenumbers
 *  - expo   = exponent of power law
 *  - ispect = integer flag which specifies spectrum
 *
 *  Note that the fourier amplitudes are stored in an array with no
 *  ghost zones
 */
static void pspect(ath_fft_data *ampl)
{
  int i,j,k;
  double q1,q2,q3;

  /* set random amplitudes with gaussian deviation */
  for (k=0; k<nx3; k++) {
    for (j=0; j<nx2; j++) {
      for (i=0; i<nx1; i++) {
        q1 = ran2(&rseed);
        q2 = ran2(&rseed);
        q3 = sqrt(-2.0*log(q1+1.0e-20))*cos(2.0*PI*q2);
        q1 = ran2(&rseed);
        ampl[OFST(i,j,k)][0] = q3*cos(2.0*PI*q1);
        ampl[OFST(i,j,k)][1] = q3*sin(2.0*PI*q1);
      }
    }
  }

  /* set power spectrum
   *   ispect=1: power law - original form
   *   ispect=2: form from Gammie&Ostriker
   */
  for (k=0; k<nx3; k++) {
    for (j=0; j<nx2; j++) {
      for (i=0; i<nx1; i++) {
        /* compute k/dkx */
        q3 = KWVM(i,j,k);
        if ((q3 > klow) && (q3 < khigh)) {
          q3 *= dkx; /* multiply by 2 pi/L */
          if (ispect == 1) {
            /* decreasing power law */
            ampl[OFST(i,j,k)][0] /= pow(q3,expo);
            ampl[OFST(i,j,k)][1] /= pow(q3,expo);
          } else if (ispect == 2) {
            /* G&O form */
            ampl[OFST(i,j,k)][0] *= pow(q3,3.0)*exp(-4.0*q3/kpeak);
            ampl[OFST(i,j,k)][1] *= pow(q3,3.0)*exp(-4.0*q3/kpeak);
          }
        } else {
          /* introduce cut-offs at klow and khigh */
          ampl[OFST(i,j,k)][0] = 0.0;
          ampl[OFST(i,j,k)][1] = 0.0;
        }
      }
    }
  }
  ampl[0][0] = 0.0;
  ampl[0][1] = 0.0;

  return;
}

/* ========================================================================== */

/*! \fn static void project()
 *  \brief Makes velocity perturbations divergence free
 */
static void project()
{
  int i,j,k,m,ind;
  double kap[3], kapn[3], mag;
  ath_fft_data dot;

  /* Project off non-solenoidal component of velocity */
  for (k=0; k<nx3; k++) {
    kap[2] = sin(2.0*PI*(gks+k)/gnx3);
    for (j=0; j<nx2; j++) {
      kap[1] = sin(2.0*PI*(gjs+j)/gnx2);
      for (i=0; i<nx1; i++) {
        if (((gis+i)+(gjs+j)+(gks+k)) != 0) {
          kap[0] = sin(2.0*PI*(gis+i)/gnx1);
          ind = OFST(i,j,k);

          /* make kapn a unit vector */
          mag = sqrt(SQR(kap[0]) + SQR(kap[1]) + SQR(kap[2]));
          for (m=0; m<3; m++) kapn[m] = kap[m] / mag;

          /* find fv_0 dot kapn */
          dot[0] = fv1[ind][0]*kapn[0]+fv2[ind][0]*kapn[1]+fv3[ind][0]*kapn[2];
          dot[1] = fv1[ind][1]*kapn[0]+fv2[ind][1]*kapn[1]+fv3[ind][1]*kapn[2];

          /* fv = fv_0 - (fv_0 dot kapn) * kapn */
          fv1[ind][0] -= dot[0]*kapn[0];
          fv2[ind][0] -= dot[0]*kapn[1];
          fv3[ind][0] -= dot[0]*kapn[2];

          fv1[ind][1] -= dot[1]*kapn[0];
          fv2[ind][1] -= dot[1]*kapn[1];
          fv3[ind][1] -= dot[1]*kapn[2];
        }
      }
    }
  }

  return;
}

/* ========================================================================== */

/*! \fn static inline void transform()
 *  \brief Generate velocities from fourier transform
 */
static inline void transform()
{
  /* Transform velocities from k space to physical space */
  ath_3d_fft(plan, fv1);
  ath_3d_fft(plan, fv2);
  ath_3d_fft(plan, fv3);

  /* Should technically renormalize (divide by gnx1*gnx2*gnx3) here, but
   * since we're going to renormalize to get the desired energy injection
   * rate anyway, there's no point */

  return;
}

/* ========================================================================== */

/*! \fn static inline void generate()
 *  \brief Generate the velocity perturbations
 */
static inline void generate()
{
  /* Generate new perturbations following appropriate power spectrum */
  pspect(fv1);
  pspect(fv2);
  pspect(fv3);

  /* Require div V = 0 */
  project();

  /* Transform perturbations to real space, but don't normalize until
   * just before we apply them in perturb() */
  transform();

  return;
}

/* ========================================================================== */

/* ========================================================================== */

/*! \fn static void perturb(Grid *pGrid, Real dt)
 *  \brief  Shifts velocities so no net momentum change, normalizes to keep
 *  dedt fixed, and then sets velocities
 */
static void perturb(GridS *pGrid, Real dt)
{
  int i, is=pGrid->is, ie = pGrid->ie;
  int j, js=pGrid->js, je = pGrid->je;
  int k, ks=pGrid->ks, ke = pGrid->ke;
  int ind, mpierr;
  Real dvol, aa, b, c, s, de, qa, v1, v2, v3;
  Real t0, t0ij, t0i, t1, t1ij, t1i;
  Real t2, t2ij, t2i, t3, t3ij, t3i;
  Real m[4], gm[4];

  /* Set the velocities in real space */
  dvol = 1.0/((Real)(gnx1*gnx2*gnx3));
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        ind = OFST(i-is,j-js,k-ks);
        dv1[k][j][i] = fv1[ind][0]*dvol;
        dv2[k][j][i] = fv2[ind][0]*dvol;
        dv3[k][j][i] = fv3[ind][0]*dvol;
      }
    }
  }

  /* Calculate net momentum pertubation components t1, t2, t3 */
  t0 = 0.0;  t1 = 0.0;  t2 = 0.0;  t3 = 0.0;
  for (k=ks; k<=ke; k++) {
    t0ij = 0.0;  t1ij = 0.0;  t2ij = 0.0;  t3ij = 0.0;
    for (j=js; j<=je; j++) {
      t0i = 0.0;  t1i = 0.0;  t2i = 0.0;  t3i = 0.0;
      for (i=is; i<=ie; i++) {
        t0i += pGrid->U[k][j][i].d;

        /* The net momentum perturbation */
        t1i += pGrid->U[k][j][i].d * dv1[k][j][i];
        t2i += pGrid->U[k][j][i].d * dv2[k][j][i];
        t3i += pGrid->U[k][j][i].d * dv3[k][j][i];
      }
      t0ij += t0i;  t1ij += t1i;  t2ij += t2i;  t3ij += t3i;
    }
    t0 += t0ij;  t1 += t1ij;  t2 += t2ij;  t3 += t3ij;
  }

#ifdef MPI_PARALLEL
  /* Sum the perturbations over all processors */
  m[0] = t0;  m[1] = t1;  m[2] = t2;  m[3] = t3;
  mpierr = MPI_Allreduce(m, gm, 4, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  if (mpierr) ath_error("[normalize]: MPI_Allreduce error = %d\n", mpierr);
  t0 = gm[0];  t1 = gm[1];  t2 = gm[2];  t3 = gm[3];
#endif /* MPI_PARALLEL */

  /* Subtract the mean velocity perturbation so that the net momentum
   * perturbation is zero. */
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        dv1[k][j][i] -= t1/t0;
        dv2[k][j][i] -= t2/t0;
        dv3[k][j][i] -= t3/t0;
      }
    }
  }

  /* Calculate unscaled energy of perturbations */
  t1 = 0.0;  t2 = 0.0;
  for (k=ks; k<=ke; k++) {
    t1ij = 0.0;  t2ij = 0.0;
    for (j=js; j<=je; j++) {
      t1i = 0.0;  t2i = 0.0;
      for (i=is; i<=ie; i++) {
        /* Calculate velocity pertubation at cell center from
         * perturbations at cell faces */
        v1 = dv1[k][j][i];
        v2 = dv2[k][j][i];
        v3 = dv3[k][j][i];

        t1i += (pGrid->U[k][j][i].d)*(SQR(v1) + SQR(v2) + SQR(v3));
        t2i +=  (pGrid->U[k][j][i].M1)*v1 + (pGrid->U[k][j][i].M2)*v2 +
                     (pGrid->U[k][j][i].M3)*v3;
      }
      t1ij += t1i;  t2ij += t2i;
    }
    t1 += t1ij;  t2 += t2ij;
  }

#ifdef MPI_PARALLEL
  /* Sum the perturbations over all processors */
  m[0] = t1;  m[1] = t2;
  mpierr = MPI_Allreduce(m, gm, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  if (mpierr) ath_error("[normalize]: MPI_Allreduce error = %d\n", mpierr);
  t1 = gm[0];  t2 = gm[1];
#endif /* MPI_PARALLEL */

  /* Rescale to give the correct energy injection rate */
  dvol = pGrid->dx1*pGrid->dx2*pGrid->dx3;
  if (idrive == 0) {
    /* driven turbulence */
    de = dedt*dt;
  } else {
    /* decaying turbulence (all in one shot) */
    de = dedt;
  }
  aa = 0.5*t1;
  aa = MAX(aa,1.0e-20);
  b = t2;
  c = -de/dvol;
  if(b >= 0.0)
    s = (-2.0*c)/(b + sqrt(b*b - 4.0*aa*c));
  else
    s = (-b + sqrt(b*b - 4.0*aa*c))/(2.0*aa);

  if (isnan(s)) ath_error("[perturb]: s is NaN!\n");

  /* Apply momentum pertubations */
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        qa = s*pGrid->U[k][j][i].d;
        pGrid->U[k][j][i].M1 += qa*dv1[k][j][i];
        pGrid->U[k][j][i].M2 += qa*dv2[k][j][i];
        pGrid->U[k][j][i].M3 += qa*dv3[k][j][i];
      }
    }
  }

  return;
}

/* ========================================================================== */
/*! \fn static void initialize(Grid *pGrid, Domain *pD)
 *  \brief  Allocate memory and initialize FFT plans */
static void initialize(GridS *pGrid, DomainS *pD)
{
  int i, is=pGrid->is, ie = pGrid->ie;
  int j, js=pGrid->js, je = pGrid->je;
  int k, ks=pGrid->ks, ke = pGrid->ke;
  int nbuf, mpierr, nx1gh, nx2gh, nx3gh;
  float kwv, kpara, kperp;
  char donedrive = 0;

/* -----------------------------------------------------------
 * Variables within this block are stored globally, and used
 * within preprocessor macros.  Don't create variables with
 * these names within your function if you are going to use
 * OFST(), KCOMP(), or KWVM() within the function! */

  /* Get local grid size */
  nx1 = (ie-is+1);
  nx2 = (je-js+1);
  nx3 = (ke-ks+1);

  /* Get global grid size */
  gnx1 = pD->Nx[0];
  gnx2 = pD->Nx[1];
  gnx3 = pD->Nx[2];

  /* Get extents of local FFT grid in global coordinates */
  gis=is+pGrid->Disp[0];  gie=ie+pGrid->Disp[0];
  gjs=js+pGrid->Disp[1];  gje=je+pGrid->Disp[1];
  gks=ks+pGrid->Disp[2];  gke=ke+pGrid->Disp[2];
/* ----------------------------------------------------------- */

  /* Get size of arrays with ghost cells */
  nx1gh = nx1 + 2*nghost;
  nx2gh = nx2 + 2*nghost;
  nx3gh = nx3 + 2*nghost;

  /* Get input parameters */

  /* interval for generating new driving spectrum; also interval for
   * driving when IMPULSIVE_DRIVING is used */
  dtdrive = par_getd("problem","dtdrive");
#ifdef MHD
  /* magnetic field strength */
  beta = par_getd("problem","beta");
  /* beta = isothermal pressure/magnetic pressure */
  B0 = sqrt(2.0*Iso_csound2*rhobar/beta);
#endif /* MHD */
  /* energy injection rate */
  dedt = par_getd("problem","dedt");

  /* parameters for spectrum */
  ispect = par_geti("problem","ispect");
  if (ispect == 1) {
    expo = par_getd("problem","expo");
  } else if (ispect == 2) {
    kpeak = par_getd("problem","kpeak")*2.0*PI;
  } else {
    ath_error("Invalid value for ispect\n");
  }
  /* Cutoff wavenumbers of spectrum */
  klow = par_getd("problem","klow"); /* in integer units */
  khigh = par_getd("problem","khigh"); /* in integer units */
  dkx = 2.0*PI/(pGrid->dx1*gnx1); /* convert k from integer */

  /* Driven or decaying */
  idrive = par_geti("problem","idrive");
  if ((idrive < 0) || (idrive > 1)) ath_error("Invalid value for idrive\n");
  /* If restarting with decaying turbulence, no driving necessary. */
  if ((idrive == 1) && (pGrid->time > 0.)) {
    donedrive = 1;
  }

  if (donedrive == 0) {
    /* Allocate memory for components of velocity perturbation */
    if ((dv1=(Real***)calloc_3d_array(nx3gh,nx2gh,nx1gh,sizeof(Real)))==NULL) {
      ath_error("[problem]: Error allocating memory for vel pert\n");
    }
    if ((dv2=(Real***)calloc_3d_array(nx3gh,nx2gh,nx1gh,sizeof(Real)))==NULL) {
      ath_error("[problem]: Error allocating memory for vel pert\n");
    }
    if ((dv3=(Real***)calloc_3d_array(nx3gh,nx2gh,nx1gh,sizeof(Real)))==NULL) {
      ath_error("[problem]: Error allocating memory for vel pert\n");
    }
  }

  /* Initialize the FFT plan */
  plan = ath_3d_fft_quick_plan(pD, NULL, ATH_FFT_BACKWARD);

  /* Allocate memory for FFTs */
  if (donedrive == 0) {
    fv1 = ath_3d_fft_malloc(plan);
    fv2 = ath_3d_fft_malloc(plan);
    fv3 = ath_3d_fft_malloc(plan);
  }

  /* Enroll outputs */
//  dump_history_enroll(hst_dEk,"<dE_K>");
//  dump_history_enroll(hst_dEb,"<dE_B>");

  return;
}


/*--------------------------------------------------------------------------- */
/*! \fn Real Erf(Real z)
 *  \brief Error function  */
Real Erf(Real z)
{
  /* arrays of the error function y=erf(x) */
  static double x[101]={
        0.000000e+000,  3.783387e-003,  7.709914e-003,  1.178500e-002,  1.601425e-002,  2.040352e-002,
        2.495885e-002,  2.968653e-002,  3.459307e-002,  3.968525e-002,  4.497008e-002,  5.045486e-002,
        5.614715e-002,  6.205480e-002,  6.818596e-002,  7.454909e-002,  8.115295e-002,  8.800667e-002,
        9.511969e-002,  1.025018e-001,  1.101632e-001,  1.181145e-001,  1.263667e-001,  1.349310e-001,
        1.438193e-001,  1.530440e-001,  1.626176e-001,  1.725534e-001,  1.828652e-001,  1.935671e-001,
        2.046738e-001,  2.162008e-001,  2.281639e-001,  2.405796e-001,  2.534651e-001,  2.668380e-001,
        2.807169e-001,  2.951209e-001,  3.100699e-001,  3.255844e-001,  3.416859e-001,  3.583966e-001,
        3.757395e-001,  3.937386e-001,  4.124186e-001,  4.318054e-001,  4.519256e-001,  4.728071e-001,
        4.944786e-001,  5.169701e-001,  5.403124e-001,  5.645379e-001,  5.896800e-001,  6.157732e-001,
        6.428537e-001,  6.709587e-001,  7.001271e-001,  7.303990e-001,  7.618162e-001,  7.944220e-001,
        8.282614e-001,  8.633812e-001,  8.998296e-001,  9.376570e-001,  9.769156e-001,  1.017659e+000,
        1.059945e+000,  1.103830e+000,  1.149376e+000,  1.196644e+000,  1.245701e+000,  1.296614e+000,
        1.349454e+000,  1.404292e+000,  1.461205e+000,  1.520272e+000,  1.581573e+000,  1.645193e+000,
        1.711221e+000,  1.779746e+000,  1.850864e+000,  1.924673e+000,  2.001274e+000,  2.080774e+000,
        2.163281e+000,  2.248909e+000,  2.337778e+000,  2.430008e+000,  2.525728e+000,  2.625070e+000,
        2.728170e+000,  2.835170e+000,  2.946219e+000,  3.061469e+000,  3.181080e+000,  3.305216e+000,
        3.434048e+000,  3.567755e+000,  3.706521e+000,  3.850536e+000,  4.000000e+000
  };
  static double y[101]={
        0.00000000e+000,  4.26907434e-003,  8.69953340e-003,  1.32973284e-002,  1.80686067e-002,  2.30197153e-002,
        2.81572033e-002,  3.34878242e-002,  3.90185379e-002,  4.47565113e-002,  5.07091186e-002,  5.68839404e-002,
        6.32887618e-002,  6.99315688e-002,  7.68205444e-002,  8.39640613e-002,  9.13706742e-002,  9.90491090e-002,
        1.07008250e-001,  1.15257124e-001,  1.23804883e-001,  1.32660778e-001,  1.41834139e-001,  1.51334337e-001,
        1.61170754e-001,  1.71352743e-001,  1.81889576e-001,  1.92790394e-001,  2.04064148e-001,  2.15719527e-001,
        2.27764884e-001,  2.40208149e-001,  2.53056730e-001,  2.66317410e-001,  2.79996226e-001,  2.94098338e-001,
        3.08627885e-001,  3.23587825e-001,  3.38979770e-001,  3.54803790e-001,  3.71058224e-001,  3.87739454e-001,
        4.04841688e-001,  4.22356710e-001,  4.40273635e-001,  4.58578645e-001,  4.77254725e-001,  4.96281391e-001,
        5.15634428e-001,  5.35285634e-001,  5.55202571e-001,  5.75348359e-001,  5.95681482e-001,  6.16155658e-001,
        6.36719759e-001,  6.57317799e-001,  6.77889021e-001,  6.98368078e-001,  7.18685336e-001,  7.38767318e-001,
        7.58537287e-001,  7.77916009e-001,  7.96822665e-001,  8.15175962e-001,  8.32895397e-001,  8.49902691e-001,
        8.66123358e-001,  8.81488386e-001,  8.95935967e-001,  9.09413237e-001,  9.21877939e-001,  9.33299942e-001,
        9.43662512e-001,  9.52963249e-001,  9.61214608e-001,  9.68443923e-001,  9.74692870e-001,  9.80016358e-001,
        9.84480847e-001,  9.88162149e-001,  9.91142807e-001,  9.93509200e-001,  9.95348535e-001,  9.96745927e-001,
        9.97781755e-001,  9.98529475e-001,  9.99054014e-001,  9.99410828e-001,  9.99645625e-001,  9.99794704e-001,
        9.99885782e-001,  9.99939162e-001,  9.99969080e-001,  9.99985060e-001,  9.99993164e-001,  9.99997050e-001,
        9.99998805e-001,  9.99999548e-001,  9.99999841e-001,  9.99999948e-001,  9.99999985e-001
  };
  double z1, g;
  int il, iu, im;
  z1 = fabs(z);
  /* look for the location of z1 in the x array */
  il=0; iu=100;
  while(iu-il>1) {
    im = (iu+il)/2;
    if (x[im] > z1) iu = im;
    else il = im;
  }
  /* linear interpolation */
  g = (y[iu]*(z1-x[il])+y[il]*(x[iu]-z1))/(x[iu]-x[il]);

  g = MIN(g,1.0);
  if (z<0.0) g = -g;
  return g;
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
