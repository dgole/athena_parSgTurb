#ifndef DEFINITIONS_H
#define DEFINITIONS_H 
/*  
 *  WARNING! This file has been automatically generated by configure.
 *  Any changes to it will be overwritten the next time configure is run.
 */
/*==============================================================================
 * FILE: defs.h.in
 *
 * PURPOSE: Template file for defs.h.  When 'configure' is run, a new defs.h
 *   file will be created (overwriting the last) from this template in which
 *   various cpp macros are defined selected from the options available here.  
 *
 * TO BY-PASS CONFIGURE: copy this file into defs.h, and edit the cpp macros
 *   by hand.
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*  macros which define physics and algorithm
 *  (user modified via configure) */

/* Version identifier when configure was last run */
#define CONFIGURE_DATE "Thu Nov 15 13:38:53 CST 2018"

/* Problem generator with which Athena is compiled */
#define A_PROBLEM "par_strat3d_turb"

/* HYDRO or MHD */
#define HYDRO

/* ADIABATIC or ISOTHERMAL or ... */
#define ISOTHERMAL
#define EOS_STR "ISOTHERMAL"

#if defined(ISOTHERMAL) /* || defined PIECEWISE_POLYTROPE ... */
#define BAROTROPIC        /* P = P(rho) */
#elif !defined(ADIABATIC)  /* P = P(eth) = (gamma - 1)*eth */
#define GENERAL_EOS       /* P = P(rho,eth) */
#endif

/* Coordinate system: CARTESIAN or CYLINDRICAL */
#define CARTESIAN

/* Number of passively advected scalars */
#define NSCALARS 0

/* Self-gravity */
#define NO_SELF_GRAVITY
#define SELF_GRAVITY_NONE

/* Particles */
#define PARTICLES
#define FEEDBACK

/* Particle self-gravity */
#define PARTICLE_SELF_GRAVITY
#define PARTICLE_SELF_GRAVITY_USING_FFT_DISK

/* implicit cooling */
#define NO_COOLING

/* resistivity, viscosity, and thermal conduction */
#define NO_RESISTIVITY
#define NO_VISCOSITY
#define NO_THERMAL_CONDUCTION
#define NO_STS

/* special relativity */
#define NEWTONIAN

/* order of spatial reconstruction: FIRST_ORDER,
 * SECOND_ORDER_CHAR, SECOND_ORDER_PRIM, THIRD_ORDER_CHAR, THIRD_ORDER_PRIM */
#define THIRD_ORDER_PRIM

/* flux type
 * ROE_FLUX, HLLE_FLUX, HLLC_FLUX, HLLD_FLUX, FORCE_FLUX, EXACT_FLUX,
 * TWO_SHOCK_FLUX */
#define HLLC_FLUX

/* unsplit integrator:
 * CTU_INTEGRATOR or VL_INTEGRATOR */
#define CTU_INTEGRATOR

/* Real: DOUBLE_PREC or SINGLE_PREC */
#define DOUBLE_PREC

/* debug mode: DEBUG or OPTIMIZE */
#define OPTIMIZE

/* Write ghost cells in outputs: WRITE_GHOST_CELLS or NO_WRITE_GHOST_CELLS */
#define NO_WRITE_GHOST_CELLS

/* MPI parallelism: MPI_PARALLEL or NO_MPI_PARALLEL */
#define MPI_PARALLEL

/* H-correction: H_CORRECTION or NO_H_CORRECTION */
#define NO_H_CORRECTION

/* FFT mode: FFT_ENABLED or NO_FFT */
#define FFT_ENABLED

/* shearing-box: SHEARING_BOX or NO_SHEARING_BOX */
#define SHEARING_BOX

/* fargo: FARGO or NO_FARGO */
#define FARGO

/* rotating-frame: ROTATING_FRAME or NO_ROTATING_FRAME */
#define NO_ROTATING_FRAME

/* l1_inflow: L1_INFLOW or NO_L1_INFLOW */
#define NO_L1_INFLOW

/* Mesh Refinement mode: STATIC_MESH_REFINEMENT or NO_MESH_REFINEMENT */
#define NO_MESH_REFINEMENT

/* First order flux correction in VL integrator:
 * FIRST_ORDER_FLUX_CORRECTION or NO_FIRST_ORDER_FLUX_CORRECTION */
#define NO_FIRST_ORDER_FLUX_CORRECTION

/*----------------------------------------------------------------------------*/
/* macros associated with numerical algorithm (rarely modified) */

/* nghost = Number of Ghost Cells 
 * num_digit = Number of digits in data dump file
 * MAXLEN = maximum line length in input parameter file
 */

/* Number of ghost cells must be 5 with particles and 3rd order */
enum {
#ifdef PARTICLES 
#if defined(THIRD_ORDER_CHAR) || defined(THIRD_ORDER_PRIM)
  nghost = 5,
#else
  nghost = 4,
#endif
#else
  nghost = 4,
#endif
  num_digit = 4
};
#define MAXLEN 256

/*----------------------------------------------------------------------------*/
/* general purpose macros (never modified) */
#ifndef MIN
#define MIN(a,b) ( ((a) < (b)) ? (a) : (b) )
#endif
#ifndef MAX
#define MAX(a,b) ( ((a) > (b)) ? (a) : (b) )
#endif
#define SIGN(a) ( ((a) < 0.) ? -1. : 1. )
#define SQR(x) ( (x)*(x) )
#define CUBE(x) ( (x)*(x)*(x) )
#define STR(x) #x
#define SQRT2 1.4142135623730951
#define ONE_OVER_SQRT2 0.7071067811865475
#define PI       3.14159265358979323846
#define ONE_3RD  0.3333333333333333
#define TWO_3RDS 0.6666666666666667
#define FOUR_3RDS 1.333333333333333
#define TINY_NUMBER 1.0e-20
#define HUGE_NUMBER 1.0e+20

/*----------------------------------------------------------------------------*/
/* computed macros based on above choices (never modified) */

#ifdef BAROTROPIC /* BAROTROPIC EOS */
#ifdef HYDRO
 enum {NWAVE = 4, NVAR = 4 + NSCALARS};
#endif
#ifdef MHD
 enum {NWAVE = 6, NVAR = 7 + NSCALARS};
#endif
#else /* ADIABATIC or other EOS */
#ifdef HYDRO
 enum {NWAVE = 5, NVAR = 5 + NSCALARS};
#endif
#ifdef MHD
 enum {NWAVE = 7, NVAR = 8 + NSCALARS};
#endif
#endif /* EOS */

/*----------------------------------------------------------------------------*/

#ifdef MPI_PARALLEL
/* Integer constants to identify MPI messages sent in various parts of code */
enum {LtoR_tag,
      RtoL_tag,
      boundary_particle_tag,
      boundary_feedback_tag,
      shearing_sheet_ix1_tag,
      shearing_sheet_ox1_tag,
#if defined(SELF_GRAVITY) || defined(PARTICLE_SELF_GRAVITY) 
      shearing_sheet_grav_ix1_tag,
      shearing_sheet_grav_ox1_tag,
      remapvar_tag,
#endif
      remapFlx_tag,
      fargo_tag,
      ch_rundir0_tag,
      ch_rundir1_tag
};
#endif /* MPI_PARALLEL */

#ifdef SHEARING_BOX
/* integer constants to denote direction of 2D slice in shearing box */
enum SS2DCoord {xy, xz};
#endif

#endif /* DEFINITIONS_H */
