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

/// maximum number of interactions for each particle.
/// With hMaxFactor=4.0 (hMax=0.4), a tensile particle at the edge of
/// initial-density material (rho=1.0, n=400/unit^2) has up to
/// pi*0.16*400 ~ 200 neighbors; 400 gives safe headroom.
#define MAX_NUM_INTERACTIONS 400

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
Supported: ideal gas (=0), Murnaghan (=1), Tillotson (=2). */
// #define EOS 0 // Ideal gas EOS
// #define EOS 1 // Murnaghan EOS
#define EOS 2 // Tillotson EOS

/* Simulate elastic dynamics*/
#define ELASTIC 1
// Shear modulus mu now lives per-material in the EOS material table
// (see EquationOfState::MurnaghanMaterial.mu); look it up via
// MeshlessEOS->EOSShearModulus(matId[i]).
/// Enable von Mises plasticity (1 = on, 0 = off).
/// When enabled, the deviatoric stress is clamped to the yield surface
/// after each integration half-step via a radial return.
#define PLASTICITY 0
/// Von Mises yield stress Y0 (only used when PLASTICITY is 1).
#define YIELD_STRESS 0.1

// Enable tensile correction. This activates computation of f_ab
#define TENSILE_CORRECTION 1
// Easiest form of tensile correction: Reduce effective Pressure after Riemann problem.
// Using f_ab calculated as in Monaghan et al, 2000, implemented as in GIZMO.
#define TENSILE_CORRECTION_1 1
// Two Parameters for the tensile correction facor, as in Monaghan, 2000
#define TENSILE_CORRECTION_PREFAC 0.2
#define TENSILE_CORRECTION_POWER 4.

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


#define USE_MATID 1
#if EOS == 1 && !USE_MATID
#error "Murnaghan EOS (EOS == 1) requires USE_MATID == 1: per-particle matId drives the material table."
#endif
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

// Hopkins/GIZMO surface volume closure correction (Reinhardt & Stadel 2017,
// arXiv:1701.08296, Fig. 3 calibration; GIZMO flag HYDRO_KERNEL_SURFACE_VOLCORR).
// When 1, each particle's kernel-sum density is divided by a per-particle
// FaceClosureError(FCE), computed from the asymmetry of its kernel neighbour
// distribution:
//   xi_i  = |sum_j W_ij r_ij| / (h_i * Omega_i)         (dimensionless)
//   FCE_i = clamp(SURFACE_VOLCORR_A - SURFACE_VOLCORR_B*xi_i,
//                 SURFACE_VOLCORR_FLOOR, 1.0)
//   rho_i = m_i * Omega_i / FCE_i
// FCE = 1 in the bulk (symmetric stencil), drops to ~0.65 at a flat free
// surface, and the floor 0.344 caps the maximum bump (~ x3) at 2D corners.
// The corrected (FCE/Omega) is also used as the volume V_i in compEffectiveFace
// and sumVolume, mirroring GIZMO's downstream use of `Density` as 1/V.
// Set to 0 to disable the correction entirely (rho_i = m_i*Omega_i).
#define SURFACE_VOLCORR 1
#define SURFACE_VOLCORR_A 1.0259
#define SURFACE_VOLCORR_B 2.52444
#define SURFACE_VOLCORR_FLOOR 0.344301

#if SURFACE_VOLCORR && PERIODIC_BOUNDARIES
#error "SURFACE_VOLCORR is only implemented for the non-periodic path. The closure asymmetry sum currently ignores ghost neighbours, which would spuriously trigger the correction at periodic boundaries. Set SURFACE_VOLCORR to 0 for periodic runs."
#endif

// GIZMO HYDRO_EXPLICITLY_INTEGRATE_VOLUME port (Monaghan 2000 continuity).
// When 1, the kernel-sum density rho_kern produced by compDensity() is
// replaced everywhere downstream (pressure, gradients, Riemann, fluxes,
// stress integration, output) by an explicitly-integrated density rhoExplicit
// evolving by
//     d ln rhoExplicit / dt = - div v
// The integrated value is relaxed back toward rho_kern in log-space on the
// timescale t_relax ~ (1/EXPLICIT_VOL_RELAX_COEF) * L_grad / c_eff with
//     L_grad = max(rho/|grad rho|, sml)
//     c_eff  = min(c_sound, sqrt(mu/rho))  (deviatoric-wave speed if ELASTIC)
// The advection step's argument is clamped to +/-EXPLICIT_VOL_DIVV_CLAMP so a
// single rogue gradient cannot collapse or blow up the density.
// Two half-step kicks bracket the flux/updateState block in MeshlessScheme::run(),
// mirroring GIZMO's drift-kick-drift structure (each half-kick uses the
// gradients available at that point: pre-update for kick A, post-update for B).
// Bisection on the colliding-rings test (MA-Obsidian/2026-05-25-bisection-results.md)
// identified this as the single highest-impact missing GIZMO feature for elastic
// rebound; without it the rings stick instead of bouncing.
#define EXPLICIT_VOL_INTEGRATION 0
/// Maximum |div v * dt/2| inside the half-step exponential. Caps the per-step
/// density change at exp(EXPLICIT_VOL_DIVV_CLAMP) ~ 4.5 (matches GIZMO).
#define EXPLICIT_VOL_DIVV_CLAMP 1.5
/// Coefficient in delta = COEF * (dt/2) * c_eff / L_grad. 0.1 gives
/// t_relax ~ 10 * L_grad / c_eff, the GIZMO default.
#define EXPLICIT_VOL_RELAX_COEF 0.1
/// Whether the log-space relaxation toward the kernel-sum density actually
/// fires. GIZMO's kicks.c:317 computes the relaxed value `qn` but then sets
/// Density_ExplicitInt = exp(q0) -- i.e. it discards the relaxation and keeps
/// the PURE-ADVECTION (Monaghan continuity) density. The colliding-rings
/// baseline that bounces (ring_sep=18.4) runs with this advection-only form.
/// Enabling the relaxation (qn) drags rhoExplicit toward the FCE-bumped,
/// surface-noisy kernel density every step; at the free surface this injects
/// noise into the evolved density and triggers the rim instability that
/// crashed the first port at t~51 (Surface_issiues/experiment_setup_20260602.md).
/// 0 = advection only (exp(q0), matches GIZMO -- DEFAULT); 1 = relax (exp(qn),
/// legacy demonstrator behaviour). When 0 the relaxation block (Lgrad / cEff /
/// rhoGrad) is dead code, so the stale-gradient concerns there are moot.
#define EXPLICIT_VOL_RELAX 0
/// Whether the working `rho` stays frozen at the once-per-step override value
/// during the hydro pass. GIZMO freezes Density throughout density->pressure->
/// gradients->Riemann->flux (the swap to Density_ExplicitInt happens once in
/// density.c:982-987 and is finalised once at kicks.c:395, mode==1); the
/// advection only accumulates into Density_ExplicitInt across the two half-kicks.
/// 1 = freeze: half-step kick A leaves `rho` untouched so pressure, gradients,
///     faces and Riemann reconstruction all see a single self-consistent density;
///     `rho` is finalised from rhoExplicit only at end-of-step (kick B). DEFAULT.
/// 0 = legacy: re-sync rho = rhoExplicit inside every half-kick (desyncs P and
///     rho for the flux pass, since pressure was computed before kick A).
#define EXPLICIT_VOL_FREEZE_RHO 1

#if EXPLICIT_VOL_INTEGRATION && PERIODIC_BOUNDARIES
#error "EXPLICIT_VOL_INTEGRATION currently only wires the non-periodic path. The periodic gradient-recompute hooks would need additional ghost updates of rhoExplicit/rhoKernel. Set EXPLICIT_VOL_INTEGRATION to 0 for periodic runs."
#endif

#if EXPLICIT_VOL_INTEGRATION && !ELASTIC
#error "EXPLICIT_VOL_INTEGRATION currently only wires the ELASTIC branch of MeshlessScheme::run(). The non-ELASTIC path uses updateStateAndPosition() with no post-update gradient recompute, so the half-step splitting has no second-kick anchor. Wire a single full-step kick before updateStateAndPosition and lift this guard if you need !ELASTIC."
#endif

// Diagnostic: dump per-pair Riemann-flux contributions to the first particle
// whose conditionNumber crosses DIAG_COND_TRIGGER. The target particle is
// selected dynamically (highest cond above the threshold on the first step
// where any particle exceeds it), then tracked for DIAG_WINDOW_STEPS further
// calls to collectFluxes.
// DIAG_COND_ENABLE: set to 1 to enable, 0 to disable.
// DIAG_COND_TRIGGER: threshold as a double literal (not used in #if).
#define DIAG_COND_ENABLE     1
#define DIAG_COND_TRIGGER    1000.
#define DIAG_WINDOW_STEPS    100

// using pressure floor to avoid predicting negative pressures
//#define PRESSURE_FLOOR 1e-8

// Use density floor
#define DENSITY_FLOOR .01

// Hopkins 2015 MFM <-> SPH gradient blend, controlled per particle by
// conditionNumber[i]. The output gradient is
//   grad = (1 - w) * grad_MFM + w * grad_SPH
// with w determined by a linear ramp:
//   cond <= COND_MAX_FOR_GRADIENT  -> w = 0  (pure MFM, bulk)
//   cond >= COND_BLEND_HI          -> w = 1  (pure SPH, sick particle)
//   in between                     -> w = (cond - lo) / (hi - lo)
// When COND_BLEND_HI <= COND_MAX_FOR_GRADIENT, the ramp degenerates to a
// hard switch at COND_MAX_FOR_GRADIENT (Attempt-6 behaviour). The SPH side
// uses the same `cubicSpline`/`dWdr` and `h_i` as MFM. Per-particle, all
// gradient quantities (rho, P, v, S) blend with the same w.
// Set COND_MAX_FOR_GRADIENT to a negative value to disable the fallback
// entirely (gradient stays MFM regardless of conditioning).
#define COND_MAX_FOR_GRADIENT 1.0e7
// Negative -> blend disabled, hard switch at COND_MAX_FOR_GRADIENT.
// (A7 with LO=100/HI=1000 and A8 with LO=10/HI=100 both regressed
// vs Attempt 6's hard switch — see sph_fallback_attempts.md.)
#define COND_BLEND_HI         -1.

// Attempt-5 neighbor-side filters (skip a neighbor whose own E is sick or
// whose stencil is sparse). DISABLED — Attempt 5 in
// MA-Obsidian/DEBUG/colliding_rings_negative_density_reconstruction.md
// found these filters made per-pair extremes worse, because the propagation
// of sick lab-velocity values through Riemann fluxes is a separate physics
// problem (free-surface tensile instability), not a gradient-estimator one.
// Kept for future experimentation; set COND_MAX_NEIGHBOR_FOR_GRADIENT > 0
// or MIN_NOI_FOR_NEIGHBOR_USE > 0 to re-enable.
#define COND_MAX_NEIGHBOR_FOR_GRADIENT -1.
#define MIN_NOI_FOR_NEIGHBOR_USE 0

#endif //DEMONSTRATOR_PARAMETER_H
