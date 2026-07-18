//
// Created by Johannes Martin on 21.09.22.
//

#ifndef DEMONSTRATOR_PARAMETER_H
#define DEMONSTRATOR_PARAMETER_H

/// possible values 2 or 3 for 2D or 3D simulations
#define DIM 2

// Order of the least-squares polynomial gradient estimator:
//   1 = linear-exact  (classic Hopkins-2015 MFM first moment matrix)
//   2 = quadratic-exact (augmented with the symmetric Hessian terms)
//   3 = cubic-exact     (also the third-derivative terms)
#define GRADIENT_ORDER 1

// The augmented moment-matrix machinery (separate psijTildeGrad array, PxP
// solve, per-particle fallback) is shared by order 2 and 3; only the basis
// size GRAD_MAT_DIM and the polynomial basis grow with the order.
#define SECOND_ORDER_GRADIENTS (GRADIENT_ORDER >= 2)

// Size of the least-squares moment matrix = number of polynomial basis terms
// (no constant term). GRAD_MAT_DIM also sizes the LAPACK scratch buffers in
// Helper, so it is always defined.
//   order 1: DIM                                             (2D: 2,  3D: 3)
//   order 2: + DIM*(DIM+1)/2          symmetric Hessian      (2D: 5,  3D: 9)
//   order 3: + DIM*(DIM+1)*(DIM+2)/6  symmetric 3rd deriv    (2D: 9,  3D: 19)
#if GRADIENT_ORDER == 1
#define GRAD_MAT_DIM (DIM)
#elif GRADIENT_ORDER == 2
#define GRAD_MAT_DIM (DIM + (DIM*(DIM+1))/2)
#else
#define GRAD_MAT_DIM (DIM + (DIM*(DIM+1))/2 + (DIM*(DIM+1)*(DIM+2))/6)
#endif

// Conditioning guard for the augmented least-squares gradient fit
// (inverseMatrixChecked). A one-sided free-surface stencil makes the moment
// matrix near-singular without making it exactly singular, so LAPACK returns a
// garbage inverse and the recovered gradient explodes -- most severely at
// order 3, where the linear and cubic basis terms share odd parity and go
// collinear on such stencils. When the Frobenius condition estimate
// kappa_F = ||M||_F * ||M^-1||_F exceeds this threshold the fit is rejected and
// the estimator falls back to the first-order weights for that particle. This
// is the primary tuning knob: lower -> fall back sooner (safer, closer to
// first order), higher -> keep the high-order fit on more marginal stencils.
#define GRAD_MAT_COND_MAX 1.0e5

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
#define MAX_NUM_INTERACTIONS 300

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
#define NNN_TARGET 32

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

/// MeshlessScheme::run() searches for neighbours ONCE per step at a fixed
/// nnsRadius = max(kernelSize, hMax_from_previous_step), then hands that
/// candidate list (nnl/noi) to updateAllSmoothingLengths(), whose Newton
/// solve only reweights the particles already in the list -- it cannot
/// discover new ones as h grows. A particle whose h has to grow a lot this
/// step (e.g. a free-surface corner with few real neighbours) can converge
/// to an h that exceeds nnsRadius, leaving it stuck with an incomplete
/// neighbour list (wrong Omega/FCE/density/gradients), even though the
/// Newton iteration itself reports converged. If h_max after the solve
/// exceeds nnsRadius, MeshlessScheme::run() redoes the whole NNS+solve pass
/// at nnsRadius = h_max * NNS_SEARCH_MARGIN, repeating until the search
/// radius covers where h actually converged (or NNS_MAX_RETRIES is hit).
#define NNS_SEARCH_MARGIN 1.5

/// safety cap on how many times MeshlessScheme::run() will redo the NNS+SML
/// solve pass per step chasing NNS_SEARCH_MARGIN headroom. Hitting this is a
/// sign of a genuinely pathological configuration (h growing without bound),
/// not just a slow-to-settle corner -- treated as fatal, like SML_PANIC_FRACTION.
#define NNS_MAX_RETRIES 5

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
#define EOS 2

/* Simulate elastic dynamics*/
#define ELASTIC 1
// Shear modulus mu now lives per-material in the EOS material table
// (see EquationOfState::MurnaghanMaterial.mu); look it up via
// MeshlessEOS->EOSShearModulus(matId[i]).
/// Legacy von Mises plasticity (1 = on, 0 = off). When enabled (and no
/// per-material plasticity model below is selected) the deviatoric stress is
/// clamped to a constant yield surface YIELD_STRESS via a radial return after
/// each stress half-step. Kept for backward compatibility; new work should use
/// the per-material *_PLASTICITY models below.
#define PLASTICITY 0
/// Constant von Mises yield stress Y0 used by the legacy PLASTICITY path.
#define YIELD_STRESS 0.1

// --- Solid failure suite (ported from miluphcuda; published model equations) ---

/// Grady-Kipp fragmentation/damage model, following Benz & Asphaug (1995).
/// Weibull-distributed activation-strain flaws are read from the IC HDF5
/// (/numFlaws, /flaws). Damage always reduces (negative) pressure; it reduces
/// the deviatoric stress only if DAMAGE_ACTS_ON_S is set.
#define FRAGMENTATION 0
#define DAMAGE_ACTS_ON_S 0
/// Maximum number of activation thresholds stored per particle. Only relevant
/// for FRAGMENTATION; keep at 1 otherwise so the flaw array stays tiny.
#define MAX_NUM_FLAWS 1
/// Skip damage evolution while simulated time < DAMAGE_START_TIME. The kernel-sum
/// density is truncated (rho < rho0, hence spurious tension everywhere) for the
/// first step or two until the smoothing length converges; without this gate that
/// startup transient latches every Weibull flaw at once. Damage models in general
/// require a relaxed IC. Default 0.0 = no delay.
#define DAMAGE_START_TIME 0.0

/// Per-material plasticity yield models (mutually exclusive). The yield
/// strength Y(P, damage, u) is supplied by EquationOfState::EOSYieldStrength
/// and fed to the same radial return as the legacy path.
/// Constant von Mises yield strength (per-material Y0).
#define VON_MISES_PLASTICITY 0
/// Mohr-Coulomb: Y = min(Y0 + mu_fric * P, Y_M).
#define MOHR_COULOMB_PLASTICITY 0
/// Drucker-Prager cone fitted to Mohr-Coulomb (cohesion + friction angle).
#define DRUCKER_PRAGER_PLASTICITY 0
/// Pressure-dependent strength following Collins et al. (2004), Jutzi (2015):
/// damage blends intact and damaged yield curves.
#define COLLINS_PLASTICITY 0
/// Simplified Collins (Lundborg strength representation only).
#define COLLINS_PLASTICITY_SIMPLE 0
/// Degrade the Collins yield strength with melt energy (uses u vs u_melt).
#define COLLINS_PLASTICITY_INCLUDE_MELT_ENERGY 0

/// Number of per-material plasticity models selected; must be 0 or 1.
#define PLASTICITY_MODEL_COUNT (VON_MISES_PLASTICITY + MOHR_COULOMB_PLASTICITY \
                                + DRUCKER_PRAGER_PLASTICITY + COLLINS_PLASTICITY)
/// True when any radial-return plasticity (legacy or per-material) is active.
#define PLASTICITY_ANY (PLASTICITY || PLASTICITY_MODEL_COUNT)

#if PLASTICITY_MODEL_COUNT > 1
#error "Select at most one of VON_MISES/MOHR_COULOMB/DRUCKER_PRAGER/COLLINS_PLASTICITY."
#endif
#if (COLLINS_PLASTICITY_SIMPLE || COLLINS_PLASTICITY_INCLUDE_MELT_ENERGY) && !COLLINS_PLASTICITY
#error "COLLINS_PLASTICITY_SIMPLE / _INCLUDE_MELT_ENERGY require COLLINS_PLASTICITY."
#endif
#if PLASTICITY_MODEL_COUNT && !ELASTIC
#error "Per-material plasticity requires ELASTIC == 1."
#endif
#if PLASTICITY_MODEL_COUNT && EOS != 2
#error "Per-material plasticity (EOSYieldStrength) is only implemented for the Tillotson EOS (EOS == 2)."
#endif
#if FRAGMENTATION && (!ELASTIC || EOS != 2)
#error "FRAGMENTATION requires ELASTIC == 1 and the Tillotson EOS (EOS == 2)."
#endif
#if COLLINS_PLASTICITY && !FRAGMENTATION
#error "COLLINS_PLASTICITY blends intact/damaged curves and requires FRAGMENTATION."
#endif
#if MAX_NUM_FLAWS < 1
#error "MAX_NUM_FLAWS must be >= 1."
#endif

// Kernel for the MFM path, GIZMO numbering (GIZMO kernel.h KERNEL_FUNCTION):
// 3 = cubic spline (default), 6 = Wendland C2, 7 = Wendland C4, 9 = Wendland C6.
#define KERNEL_FUNCTION 6

// Enable tensile correction. This activates computation of f_ab
#define TENSILE_CORRECTION 1
// Easiest form of tensile correction: Reduce effective Pressure after Riemann problem.
// Using f_ab calculated as in Monaghan et al, 2000, implemented as in GIZMO.
#define TENSILE_CORRECTION_1 1
// Damp the deviatoric stress of a negative-pressure side by (1 - f), GIZMO's
// simple sign-test branch in solids/elastic_stress_tensor_force.h ("#if 0" branch).
#define TENSILE_CORRECTION_2 0
// Damp only tensile (positive) principal components of the deviatoric stress by
// (1 - f), GIZMO's default eigenvalue branch in solids/elastic_stress_tensor_force.h.
#define TENSILE_CORRECTION_3 0
// Two Parameters for the tensile correction facor, as in Monaghan, 2000
#define TENSILE_CORRECTION_PREFAC 0.2
#define TENSILE_CORRECTION_POWER 4.
// Output per-particle average f_ab to h5 dump (requires TENSILE_CORRECTION)
#if TENSILE_CORRECTION
#define OUTPUT_FAB 0
#endif

#if TENSILE_CORRECTION_2 && TENSILE_CORRECTION_3
#error "TENSILE_CORRECTION_2 and TENSILE_CORRECTION_3 are mutually exclusive stress-damping branches (GIZMO's #if 0/#else)."
#endif
#if TENSILE_CORRECTION_3 && DIM == 3
#error "TENSILE_CORRECTION_3 is only implemented for DIM == 2 (analytic 2x2 eigendecomposition)."
#endif

// ---------------------------------------------------------------------------
// Angular-momentum conservation (20260718 study, testcases/angular_momentum_study).
// MFM pair momentum fluxes are generally NOT parallel to the pair's line of
// centers (effective faces Aij are non-central, the deviatoric stress flux
// S*Fhat is non-radial, and the transverse shear-wave HLL dissipation acts
// along dv_t), so each antisymmetric pair exchange exerts a net torque
// -(r_i - r_j) x F_ij and total L_z drains at collision contact. Selector:
//   0 = baseline MFM (no fix)
//   1 = RADIAL_FLUX: project each pair momentum flux onto the line of centers.
//       Exact L conservation by construction; transverse (shear) momentum
//       exchange between pairs is dropped (contact becomes normal-force-like).
//   2 = RADIAL_FACE: project the effective face Aij onto the line of centers
//       (fixes only the geometric non-centrality of the HLLC flux).
//   3 = NO_TDISS: drop the transverse shear-wave HLL dissipation wt_t*dv_t in
//       the GIZMO elastic stress flux (the purely numerical non-central term);
//       longitudinal dissipation is kept unchanged.
//   4 = SPIN: bookkeeping only -- accumulate each pair's residual torque into
//       a per-particle spin ledger so that L_orbital + L_spin is conserved
//       exactly; the velocity field itself stays baseline.
//   5 = GLOBAL_CORR: after each velocity update, apply the minimal-kinetic-
//       energy rigid-rotation velocity correction (about the center of mass)
//       that restores the pre-step L_z exactly. Conserves P exactly; the KE
//       change is compensated from internal energy.
//   6 = LOCAL_CORR: attribute half of each pair's residual torque to each
//       partner as a torque debt, then cancel every particle's debt by a small
//       rotation kick applied to its own neighborhood (mass-centered, so P is
//       conserved exactly and the injected L_z equals the debt exactly). KE
//       change compensated from internal energy. Local version of 5.
//   7 = RADIAL_FACE + NO_TDISS combined (fix both numerical non-central
//       pieces, keep the physical deviatoric stress flux untouched).
#define AM_METHOD 0

#define AM_RADIAL_FLUX  (AM_METHOD == 1)
#define AM_RADIAL_FACE  (AM_METHOD == 2 || AM_METHOD == 7)
#define AM_NO_TDISS     (AM_METHOD == 3 || AM_METHOD == 7)
#define AM_SPIN         (AM_METHOD == 4)
#define AM_GLOBAL_CORR  (AM_METHOD == 5)
#define AM_LOCAL_CORR   (AM_METHOD == 6)
// methods that need the per-particle residual-torque accumulator tqF
#define AM_TORQUE_TRACK (AM_SPIN || AM_LOCAL_CORR || AM_GLOBAL_CORR)

#if AM_METHOD && DIM != 2
#error "AM_METHOD > 0 is only implemented for DIM == 2 (torque about z)."
#endif
#if (AM_SPIN || AM_GLOBAL_CORR || AM_LOCAL_CORR) && !ELASTIC
#error "AM_METHOD 4-6 hook into the ELASTIC updateState path only."
#endif

// GIZMO-faithful elastic coupling: keep the HLLC Riemann problem purely isotropic
// (pressure only + dummy-pressure shift) and add the deviatoric stress as a separate
// SPH stress flux with longitudinal+transverse (shear-wave) HLL dissipation and
// per-eigenvalue tensile correction (GIZMO solids/elastic_stress_tensor_force.h).
// (consistency guards live below, after USE_HLLC is defined)
#define GIZMO_ELASTIC_FLUX 1

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

// GIZMO_ELASTIC_FLUX consistency guards (USE_HLLC/ELASTIC are now defined)
#if GIZMO_ELASTIC_FLUX && !(USE_HLLC && ELASTIC)
#error "GIZMO_ELASTIC_FLUX requires USE_HLLC == 1 and ELASTIC == 1."
#endif
// 3D MFM is supported: the GIZMO elastic stress flux and tensile-principal damping use the
// dimension-agnostic Helper::eigenDecompositionSym (3x3 eigenbasis via LAPACK), the z/vz
// coordinate+velocity plumbing is wired into the non-RUNSPH path (Particles.h declaration,
// Particles.cpp DIM==3 alloc/free), and the MFM I/O + initial distribution round-trip z/vz
// through the /x and /v datasets. (Built and linked in 3D; not yet runtime-validated.)
// GIZMO_ELASTIC_FLUX zeroes the deviatoric stress passed to the Riemann solver and
// applies its own per-eigenvalue tensile damping in the separate stress flux, so the
// in-solver stress-damping branches TC2/TC3 would be no-ops: forbid the combination.
#if GIZMO_ELASTIC_FLUX && TENSILE_CORRECTION_2
#error "GIZMO_ELASTIC_FLUX and TENSILE_CORRECTION_2 are mutually exclusive: stress damping is done in the elastic flux, not the Riemann solver."
#endif
#if GIZMO_ELASTIC_FLUX && TENSILE_CORRECTION_3
#error "GIZMO_ELASTIC_FLUX and TENSILE_CORRECTION_3 are mutually exclusive: stress damping is done in the elastic flux, not the Riemann solver."
#endif

// Use HLL solver
#define HLLC_general_EOS 1

// Riemann related: Use dummy pressure to avoid ill defined riemann problem
#define USE_DUMMY_PRESSURE 1
#if USE_DUMMY_PRESSURE == 0 && TENSILE_CORRECTION
#error "Dummy pressure needs to be activated if tensile correction is to be used."
#endif
#define USE_HLL 0

// Use Roe Average for HLL solver. Otherwise direct estimate is used.
#define USE_ROE 1

// Convert HLL(C) flux to primitive vars for debugging
//#define DEBUG_HLL 1

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
#define ARTVISC 0

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
#define EXPLICIT_VOL_INTEGRATION 1
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
/// Velocity-divergence source for the explicit-volume advection. GIZMO drives the
/// continuity integration with a robust SPH kernel-weighted estimator
/// (Particle_DivVel, density.c:334/514: -sum_j dW/dr (dp.dv)/r, normalized by the
/// kernel-sum neighbour number = omega), NOT the matrix-inverse MFM gradient trace.
/// The MFM gradient is ill-conditioned at the ring's free surface and injects a
/// spurious divergence that bloats/drifts the rings before contact (absent in GIZMO).
/// NOTE: ported faithfully, but in the demonstrator the SPH estimator is slightly
/// noisier at the free surface than the MFM trace and did not improve the rings, so
/// it is left OFF by default; the working seed-fix run uses the MFM trace.
/// 1 = GIZMO SPH estimator (faithful, noisier here); 0 = MFM gradient trace (DEFAULT).
#define EXPLICIT_VOL_SPH_DIVV 0
/// Number of initial steps to run on the kernel density before seeding the advected
/// density. The demonstrator's t=0 kernel density is truncated by the incomplete
/// startup NNS (bulk ~0.89) and self-heals by re-summation at step 1; seeding the
/// (re-summation-free) advected density from that poisoned value induces a spurious
/// breathing/expansion (ring bloat + pre-contact drift) that GIZMO does not have
/// (GIZMO has no t=0 truncation). 1 = seed from the healed step-1 density (DEFAULT);
/// 0 = legacy seed from the truncated t=0 density (bloats).
#define EXPLICIT_VOL_SEED_SKIP 1

/// Seed the explicitly-integrated density from the per-material reference density
/// rho0 (the IC density) instead of the kernel sum, and use it from step 0 (no
/// seed-skip). This mirrors miluphcuda's INTEGRATE_DENSITY, which advects density
/// from the IC value: the kernel sum under-counts at a free surface (~25% deficit
/// here), which a stiff EOS (iron Tillotson A=128 GPa) turns into multi-GPa
/// spurious tension that ejects the surface layer and drives the internal energy
/// negative. A static cold surface seeded at rho0 has divv=0, so the advected
/// density stays rho0 and the surface is stress-free, exactly as in miluphcuda.
/// Only meaningful for EOS 1/2 (needs a material rho0). 0 = seed from kernel sum
/// (DEFAULT); 1 = seed from rho0 (free-surface solids in physical units).
#define EXPLICIT_VOL_SEED_RHO0 0
#if EXPLICIT_VOL_SEED_RHO0 && EOS == 0
#error "EXPLICIT_VOL_SEED_RHO0 needs a material reference density (EOS 1 or 2)."
#endif

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
#define DIAG_COND_ENABLE     0
#define DIAG_COND_TRIGGER    1000.
#define DIAG_WINDOW_STEPS    100

// using pressure floor to avoid predicting negative pressures
//#define PRESSURE_FLOOR -1

// Use density floor
#define DENSITY_FLOOR .000001

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
