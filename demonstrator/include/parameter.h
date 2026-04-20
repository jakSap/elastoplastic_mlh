//
// Created by Johannes Martin on 21.09.22.
//

#ifndef DEMONSTRATOR_PARAMETER_H
#define DEMONSTRATOR_PARAMETER_H

/// possible values 2 or 3 for 2D or 3D simulations
#define DIM 2

/// define if periodic boundaries should be employed
#define PERIODIC_BOUNDARIES 0

/// define if timestep is adaptive
#define ADAPTIVE_TIMESTEP 1

/// define Courant-Friedrichs-Levy number, should be smaller than 1
#define CFL .2           // TODO: move to config

/// maximum number of interactions for each particle
#define MAX_NUM_INTERACTIONS 200

/// use a per-particle, iteratively determined smoothing length (Hopkins 2015,
/// see Martin's master thesis section 3.4.2). When 0, the global
/// `config.kernelSize` is used everywhere as before. Even with VARIABLE_SML 0
/// the per-particle `sml[i]` array is allocated and initialized to
/// `config.kernelSize` so the new IO path is exercised.
#define VARIABLE_SML 1

/// target effective neighbor number N_NN for the smoothing length iteration.
/// Defined as N_NN = V_kern(h_i) * sum_j W(r_ij, h_i) where
/// V_kern is pi h^2 in 2D and (4/3) pi h^3 in 3D.
/// For the single_block 2D testcase (delta_p = 0.005, kernelSize = 0.02),
/// this evaluates to ~50, which is also a reasonable default for 2D MFM.
#define NNN_TARGET 20

/// numerical tolerance for the Newton-Raphson root finding on f_SML
#define SML_TOL 1e-3

/// maximum number of Newton-Raphson iterations
#define SML_MAX_ITER 20

/// maximum allowed relative change of h per NR iteration (clamps overshoot)
#define SML_MAX_FACTOR 1.2

/// hard lower bound on h_i, as a factor of config.kernelSize.
/// Clamps every Newton step so a sick particle cannot collapse h to 0.
#define SML_H_MIN_FACTOR 0.1

/// hard upper bound on h_i, as a factor of config.kernelSize.
/// Clamps every Newton step so a sick particle cannot run h away to infinity.
/// Must stay <= the NNS search radius budget or neighbors will be missed.
#define SML_H_MAX_FACTOR 4.0

/// fraction of particles allowed to end the Newton iteration in a bad state
/// (unconverged OR clamped at hMin/hMax) before a WARN is logged.
#define SML_WARN_FRACTION 0.01

/// fraction of particles in a bad state that triggers a hard abort.
/// Above this the simulation state is considered unrecoverable.
#define SML_PANIC_FRACTION 0.1

/* Define which equation of state to use.
For now, ideal gas (=0) and murnaghan EOS (=1) are supported: */
// #define EOS 0 // Ideal gas EOS
#define EOS 1 // Murnaghan EOS
// #define EOS 2 // Tillotson EOS

/* Simulate elastic dynamics*/
#define ELASTIC 1
#define SHEAR_MODULUS 0.22
/// Enable von Mises plasticity (1 = on, 0 = off).
/// When enabled, the deviatoric stress is clamped to the yield surface
/// after each integration half-step via a radial return.
#define PLASTICITY 0
/// Von Mises yield stress Y0 (only used when PLASTICITY is 1).
#define YIELD_STRESS 0.1
// #define SHEAR_MODULUS 0
/** maximum interactions with ghost particles
 *  ignored when `PERIODIC_BOUNDARIES` is not set
**/
#define MAX_NUM_GHOST_INTERACTIONS 150

/// flag for slope limiting, 0: no slope limiting
#define SLOPE_LIMITING 1

/// slope limiting parameter, ignored when `SLOPE_LIMITING` is false
#define BETA 4.             // TODO: move to config

/// use pairwise limiter
#define PAIRWISE_LIMITER 1
#define PSI_1 .5            // TODO: move to config
#define PSI_2 .25           // TODO: move to config

/// meshless finite mass method instead of meshless finite volume
#define MESHLESS_FINITE_MASS 1


#define USE_MATID 0
// Use HLLC or HLL solver for EOS != ideal gas
#define USE_HLLC 1

// Use HLL solver
#define HLLC_general_EOS 1
#define USE_HLL 0

// Use Roe Average for HLL solver. Otherwise direct estimate is used.
#define USE_ROE 1

// Convert HLL(C) flux to primitive vars for debugging
#define DEBUG_HLL 1

/// enforcing flux symmetry by only calculating on side
#define ENFORCE_FLUX_SYM 1

/// define if particles should move, otherwise a fixed grid is used
#define MOVE_PARTICLES 1

/** define debug level to enable additional output:
 * 0: no debug additions
 * 1: additional checks
 * 2: dump NNL and ghosts to files (this should not be used for large amounts of particles)
**/
#define DEBUG_LVL 1

/// use first order quadrature point for Riemann problems
#define FIRST_ORDER_QUAD_POINT 1

/// define if code should run as SPH, which ignores most of the directives above
#define RUNSPH 0

/// define if arificial viscosity should be employed
#define ARTVISC 1

/// artificial viscosity parameters
#define ALPHA_VISC 1
#define BETA_VISC 2
#define EPSMU .01

/// define verbosity for each VERBOSITY_PARTICLES particles
// TODO: use this flag when debug level 1
#define VERBOSITY_PARTICLES 10

/// deprecated, ENFORCE_FLUX_SYM should be set to 1
// define how much tolerance of flux antisymmetry is allowed in checkFluxSymmetry
#define FLUX_SYM_TOL 1e-20

// floor for internal energy u, set to negative value to disable
#define ENERGY_FLOOR -1

// Output per-particle condition number of the gradient estimation matrix E
#define OUTPUT_CONDITION_NUMBER 1

// Use per-particle local reference density rho0 for Murnaghan EOS.
// 1: rho0 is stored per particle and initialized from initial density (local rho0)
// 0: rho0 from the global EOS parameter is used for all particles
#define LOCAL_RHO0 0

// using pressure floor to avoid predicting negative pressures
//#define PRESSURE_FLOOR 1e-8

// Use density floor
#define DENSITY_FLOOR .01

#endif //DEMONSTRATOR_PARAMETER_H
