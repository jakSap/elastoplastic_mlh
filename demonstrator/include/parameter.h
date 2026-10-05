//
// Created by Johannes Martin on 21.09.22.
//

#ifndef DEMONSTRATOR_PARAMETER_H
#define DEMONSTRATOR_PARAMETER_H

// Defaults since 2026-10-03 = the mto configuration that passes the whole testplan (rings, blocks,
// Taylor rods, rotating disc with 2nd-order angular momentum, all waves 2nd order): candidate
// testcases/mto/candidates/fin1_g2_eta1_mn.flags. Runtime part, not set here: smoothing-length
// tolerance sml.tol 1e-10 (SML_TOL is only the fallback); the material viscosity eta defaults to
// MATERIAL_ETA_DEFAULT 0.1 (config key `eta` overrides, `eta 0` switches it off). Flags documented as "0 = baseline"
// below refer to the pre-10-03 baseline, which testcases/mto/flags/base.flags still pins.

/// possible values 2 or 3 for 2D or 3D simulations
#define DIM 2

// Order of the least-squares polynomial gradient estimator:
//   1 = linear-exact  (classic Hopkins-2015 MFM first moment matrix)
//   2 = quadratic-exact (augmented with the symmetric Hessian terms)
//   3 = cubic-exact     (also the third-derivative terms)
//   4 = quartic-exact   (also the fourth-derivative terms)
#define GRADIENT_ORDER 1

// The augmented moment-matrix machinery (separate psijTildeGrad array, PxP
// solve, per-particle fallback) is shared by order 2 and 3; only the basis
// size GRAD_MAT_DIM and the polynomial basis grow with the order.
#define SECOND_ORDER_GRADIENTS (GRADIENT_ORDER >= 2)

// Use the FULL recovered Taylor polynomial in the face reconstruction, not
// only its gradient. Without this, GRADIENT_ORDER >= 2 merely makes the
// gradient itself higher-order accurate while the face extrapolation stays
// linear (f_face = f_i + grad.dx); the Hessian / third-derivative
// coefficients determined by the same least-squares solve are discarded.
// With this flag, compPolyCoefficients() recovers all P Taylor coefficients
// per particle and field (rho, P, v, and -- on the solid-HLLC path -- S), and
// the reconstruction adds the quadratic-and-above terms on top of the linear
// extrapolation. The linear part remains today's slope-limited, SPH-blended
// gradient (bit-identical when the higher coefficients are zero); particles
// on the GRAD_MAT_COND_MAX kappa-fallback get zero higher coefficients, and
// the pairwise limiter clamps the full polynomial face value, so the
// stabilization strategy of the 20260720 study carries over unchanged.
// Only meaningful when GRADIENT_ORDER >= 2 (no-op otherwise).
#define POLY_RECONSTRUCTION 0
#define POLY_RECONSTRUCTION_ACTIVE (POLY_RECONSTRUCTION && SECOND_ORDER_GRADIENTS)
// Stress-field polynomial additionally requires ELASTIC (S exists at all)
// and !GIZMO_ELASTIC_FLUX (that path uses cell-centered stress -- the
// deviatoric flux is SPH, not a face reconstruction -- so the S*PolyH
// arrays are only allocated when this is true; see the matching guard in
// Particles.h/Particles.cpp). Evaluated lazily wherever used, so it is safe
// to reference ELASTIC/GIZMO_ELASTIC_FLUX here even though both are defined
// later in this file.
#define POLY_RECONSTRUCTION_STRESS_ACTIVE (POLY_RECONSTRUCTION_ACTIVE && ELASTIC && !GIZMO_ELASTIC_FLUX)

// Size of the least-squares moment matrix = number of polynomial basis terms
// (no constant term). GRAD_MAT_DIM also sizes the LAPACK scratch buffers in
// Helper, so it is always defined.
//   order 1: DIM                                                (2D: 2,  3D: 3)
//   order 2: + DIM*(DIM+1)/2             symmetric Hessian       (2D: 5,  3D: 9)
//   order 3: + DIM*(DIM+1)*(DIM+2)/6     symmetric 3rd deriv     (2D: 9,  3D: 19)
//   order 4: + DIM*(DIM+1)*(DIM+2)*(DIM+3)/24  symmetric 4th deriv (2D: 14, 3D: 34)
#if GRADIENT_ORDER == 1
#define GRAD_MAT_DIM (DIM)
#elif GRADIENT_ORDER == 2
#define GRAD_MAT_DIM (DIM + (DIM*(DIM+1))/2)
#elif GRADIENT_ORDER == 3
#define GRAD_MAT_DIM (DIM + (DIM*(DIM+1))/2 + (DIM*(DIM+1)*(DIM+2))/6)
#else
#define GRAD_MAT_DIM (DIM + (DIM*(DIM+1))/2 + (DIM*(DIM+1)*(DIM+2))/6 + (DIM*(DIM+1)*(DIM+2)*(DIM+3))/24)
#endif

// Number of higher-than-linear Taylor coefficients stored per particle and
// field by POLY_RECONSTRUCTION (Hessian terms and above).
#define GRAD_NUM_HIGHER (GRAD_MAT_DIM - DIM)

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
// With SURFACE_GRAD_BLEND on it degrades to a pure singularity net (the smooth
// surface fade below does the discrimination), so it can be left loose.
#define GRAD_MAT_COND_MAX 1.0e5

// --- Surface-aware high-order gradient blend -------------------------------
// A more elegant surface switch than thresholding GRAD_MAT_COND_MAX. kappa_F of
// the *augmented* (order-dependent) moment matrix is a poor surface proxy: its
// bulk value grows steeply with the order (kappa_F ~ 9 / 48 / 205 for order
// 2/3/4) so the single threshold separating healthy bulk from a diverging
// one-sided surface stencil must be retuned per order AND per resolution -- and
// for order 4 the bulk value overlaps the surface tail so no threshold works.
//
// Instead, detect the surface from the FIRST-ORDER moment matrix
//   B_i = sum_j (x_j - x_i)(x_j - x_i)^T psi_j        (already built every step)
// via its eigenvalue ratio
//   s_i = lambda_min(B_i) / lambda_max(B_i)   in (0, 1].
// s_i is (a) dimensionless -> RESOLUTION-independent, and (b) evaluated at first
// order -> IDENTICAL for every GRADIENT_ORDER, so ONE threshold works for all
// orders and resolutions. s ~ 1 in an isotropic bulk stencil; s -> 0 as the
// stencil goes one-sided at a free surface (the standard Marrone et al. 2010
// delta-SPH surface detector). The higher-than-linear reconstruction is then
// faded with a smoothstep ramp
//   s >= SURFACE_BLEND_S_HI : w = 1  full high-order (bulk)
//   s <= SURFACE_BLEND_S_LO : w = 0  pure first order (surface)
// The gradient weights blend  w*high-order + (1-w)*first-order and the poly
// coefficients scale by w, so the order transition is spatially CONTINUOUS (no
// hard order jump injecting noise at the surface) and high-order usage is
// MAXIMISED: w = 1 across the whole bulk at every order, only the thin surface
// shell fades. Set to 0 to reproduce the pure GRAD_MAT_COND_MAX hard-switch
// behaviour of the 20260720 study exactly.
#define SURFACE_GRAD_BLEND 1
// Ramp endpoints on s = lambda_min/lambda_max. s ~ 1 (isotropic bulk) -> full
// high-order; ~0.4-0.5 at a flat free surface, lower at corners. The window
// [LO, HI] is where the reconstruction fades (in smoothstep) from first to full
// high order. Order- and resolution-independent by construction.
#define SURFACE_BLEND_S_LO 0.35
#define SURFACE_BLEND_S_HI 0.75

// SECOND detector: first-moment (centroid) asymmetry of the stencil,
//   xi_i = |sum_j W_ij (x_i - x_j)| / (h_i Omega_i)          (Reinhardt & Stadel)
// xi = 0 only for a genuinely symmetric stencil, and grows with ANY one-sidedness.
//
// This is needed because s = lambda_min/lambda_max is a SECOND-moment quantity and
// is structurally BLIND TO CORNERS: a corner stencil is missing roughly a quadrant,
// so the remaining support still spans both directions and B stays near-isotropic
// (measured: s ~ 0.81 on the most-damaged corner particles, i.e. ABOVE S_HI -> full
// high order), while a flat surface is missing a half-plane, collapses one principal
// direction, and is caught. The first moment has no such blind spot -- the centroid
// offset is in fact LARGEST at corners. Measured on the 300 most-damaged particles of
// the blocks collision: s flags 40%, xi flags 77%.
//
// The combined weight is the more conservative of the two ramps,
//   w = min( smoothstep((s - S_LO)/(S_HI - S_LO)),
//            smoothstep((XI_HI - xi)/(XI_HI - XI_LO)) ),
// so xi covers corners and s covers degenerate-direction stencils. xi is dimensionless
// (normalised by h_i Omega_i), so it is order- and resolution-independent like s.
// Calibration: bulk xi ~ 0; xi = 0.0103 is where the SURFACE_VOLCORR fit saturates
// (fce = 1); a flat free surface is xi ~ 0.149 (fce ~ 0.65); corners run to xi ~ 0.27.
#define SURFACE_BLEND_XI_LO 0.02
#define SURFACE_BLEND_XI_HI 0.15

/// define if periodic boundaries should be employed
#define PERIODIC_BOUNDARIES 0

#if POLY_RECONSTRUCTION_ACTIVE && PERIODIC_BOUNDARIES
#error "POLY_RECONSTRUCTION only wires the non-periodic path (compPolyCoefficients ignores ghost neighbours). Set POLY_RECONSTRUCTION to 0 for periodic runs."
#endif

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
#define SML_TOL 1e-10

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
/// Yield return of the von Mises models, GIZMO-compatible variant.
/// GIZMO (solids/elastic_physics.c:106-112, applied after every kick AND drift of S)
/// scales the deviator by f_Y = Y0^2/(NDim J2) with J2 = 1/2 S:S over the NDim x NDim
/// components -- the ratio itself, not its square root. The 2D J2 here uses the same
/// in-plane components, so with this flag the return is GIZMO's.
///   0 = radial return f = sqrt(Y^2/(d J2)), lands on the yield surface (baseline)
///   1 = GIZMO factor f = Y^2/(d J2)
#define PLASTICITY_YIELD_GIZMO 0
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
// Scale-free tensile correction (TENSILE_CORRECTION_LOCAL 1).
//
// Monaghan's f_ab = W(q)/W(dp_ref) has three dependencies it should not have:
// the reference spacing dp_ref = sqrt(pi/N_NN) is the ASSUMED spacing of an
// undisturbed lattice rather than the measured local one, the profile across
// the stencil IS the kernel shape, and the whole thing is feedforward -- on a
// perfect lattice it still applies f ~ 0.33..0.56 of damping (kernel-dependent),
// which is why the prefactor has to stay small and the correction is then too
// weak when it is actually needed.
//
// This variant replaces the ratio by
//     s_ab = r_ab / (0.5*(dp_a + dp_b)),   dp_i = (m_i/rho_i)^(1/DIM)
// and stores max(0, 1 - s_ab), so the consumers' PREFAC*pow(fab, POWER) yields
//     f = PREFAC * max(0, 1 - s_ab)^POWER
// s_ab is 1 on ANY uniform lattice regardless of kernel, N_NN or resolution, so
// f vanishes identically there (measured: 0.00000 for every configuration in
// the 27-run baseline, against 0.33-0.56 for the Monaghan form). It also stays
// 1 under uniform stretching -- r and dp_i grow together -- so the correction
// no longer fires on mere expansion, only on genuine pair collapse.
#define TENSILE_CORRECTION_LOCAL 0

// Upper clamp on the correction factor (TENSILE_CORRECTION_LOCAL only).
// Without it the (1 - s) profile keeps growing as a pair collapses and turns
// into an impulsive kick on an already-collapsed pair -- measured as negative
// internal energy at vmax = 0.60. Physically: never convert more than this
// fraction of the tensile traction into repulsion.
#define TENSILE_CORRECTION_FMAX 0.5

// Approach damping (TENSILE_CORRECTION_LOCAL only). 0 disables it.
//
// The repulsion above can push a collapsed pair apart, but it does no work
// during the APPROACH, so a fast pair arrives already collapsed and the
// repulsion then acts as an impulsive kick -- observed as negative internal
// energy at vmax = 0.60, and as a correction whose authority stops scaling
// with the drive once the FMAX clamp saturates.
//
// This boosts the EXISTING radial HLL dissipation by (1 + DAMP*f) for pairs
// that are both closer than the local spacing (f > 0) and approaching
// (dv.r < 0). Because it only scales a dissipation coefficient it inherits
// that channel's properties exactly: the force stays central and pairwise
// antisymmetric, so linear and angular momentum remain conserved to machine
// precision at any strength, and it can only REMOVE kinetic energy -- it can
// never inject any, whatever DAMP is set to.
#define TENSILE_CORRECTION_DAMP 8.0

// Two Parameters for the tensile correction facor, as in Monaghan, 2000
// (TENSILE_CORRECTION_LOCAL reuses both: PREFAC scales f, POWER sharpens the
// (1 - s) profile. Because f == 0 on a good lattice, PREFAC may be raised well
// above Monaghan's 0.2 without paying a permanent tax.)
#define TENSILE_CORRECTION_PREFAC 0.3
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

// Elastic HLLC-EP (Liu, Cheng & Liu 2019, Comput. Fluids 192, "MHLLCEP",
// elastic 3-wave case): solid HLLC built on the total normal traction
// t = p - S_nn throughout, with the contact condition t*_L = t*_R -- p and
// S_nn jump independently across the interface, so USE_DUMMY_PRESSURE is
// bypassed on this path -- and Davis-type wave speeds clamped to bracket the
// interface. The star density / deviatoric stress split (their eq. 33 log
// closure) is the hook for the later plastic-wave cases.
// TENSILE_CORRECTION_1 is supported on this path: the Monaghan anti-clumping
// repulsion fTC*max(0,-min(t_L,t_R)) is added to the normal momentum (and
// Sstar-weighted energy) flux after the solve -- the net effect of the legacy
// dummy-pressure+TC1 combination, without perturbing the wave fan. It is
// purely normal, so variant 2 keeps machine-exact L_z. Without it the elastic
// EP solvers are tensile-unstable at ring resolution dp<=0.05 (20260719 study:
// variant 1 fractures, variant 2 fills the contact zone with pairing noise).
//   0 = off (existing solvers)
//   1 = face-normal solver, replaces the solid HLLC on the ELASTIC path (A')
//   2 = pair-radial (B'): additionally project the effective face onto the
//       line of centers and drop the tangential traction, making every pair
//       flux central -> L_z conserved to machine precision by construction
#define ELASTIC_HLLC_EP 0

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
// ELASTIC_HLLC_EP == 2 reuses the radial face projection for its geometry
#define AM_RADIAL_FACE  (AM_METHOD == 2 || AM_METHOD == 7 || ELASTIC_HLLC_EP == 2)
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

// Feed the HLL velocity dissipation in the GIZMO elastic stress flux with the
// RECONSTRUCTED face jump instead of the raw particle velocity difference.
//
// Why: `dv = v_i - v_j` is non-zero for a rigid-body rotation (dv = omega x r_ij,
// purely transverse) even though the physical shear rate is exactly zero, so the
// transverse term `wt_t*dv` brakes any rotating body. The 20260720 rotating-disc
// null test measured -89% of L_z in a single rotation with no contact anywhere in
// the problem. This is the classic SPH failure mode the Balsara switch cures; the
// transverse HLL term here has no such guard. AM_METHOD 3 fixes L_z by dropping
// wt_t entirely, but that also discards the physical shear-wave damping.
//
// Both sides reconstruct to the SAME face point xij (compRiemannStatesLR), so for
// any linear velocity field -- rigid rotation and resolved shear alike -- the
// reconstructed jump vanishes identically and no dissipation is applied, while
// genuine discontinuities are still damped. That is what a Godunov scheme should
// do; the raw form is effectively 1st-order dissipative, hence the O(1) L_z error.
//
//   0 = off, raw velocity difference (baseline, bit-identical to pre-fix code)
//   1 = transverse (S-wave) channel on the face jump; longitudinal (P-wave)
//       channel keeps the raw, more diffusive form. Minimal change: targets
//       exactly the spurious rotational braking, leaves shock capturing alone.
//       VALIDATED (20260720): rotating disc L_z err 6.2e-1 -> 2.5e-5, rotating
//       block -99.9% -> -0.4%, offset rings -10.5% -> -1.8%, all matching or
//       slightly better than AM_METHOD 3, with ring integrity and machine-exact
//       E/p unchanged. Recommended setting.
//   2 = both channels on the face jump. REJECTED (20260720): removing the
//       longitudinal dissipation as well is UNSTABLE -- the rotating block
//       overshoots and reverses (L_z drift -290%, omega goes negative), the disc
//       degrades to -80% spin-down. The raw-difference P-wave term is load-
//       bearing for stability. Kept only to document the negative result.
//   3 = transverse channel on the RESIDUAL jump dv - (grad v_i + grad v_j)/2 . r_ij
//       with the UNLIMITED lab-frame gradients (option O1, mto 2026-10-03): zero for
//       any linear field like 1, but no slope limiter and no frame error (the
//       reconstructed rows WijR/WijL are rotated into the face basis by the Riemann
//       constructor, which makes 1 chiral). Longitudinal channel raw.
//   4 = both channels on the residual (cf. 2: expect the stability problem of 2
//       unless WTR_LIMITER keeps some longitudinal dissipation).
//   5 = hybrid (2D): both channels on the residual (as 4) in smooth elastic material, on the RAW jump
//       for pairs of different bodies (contact) and where a particle yields (weight
//       min(1, (1 - gamma)/RECON_FACE_DV5_WIDTH), gamma = radial-return factor): GIZMO-like
//       first-order dissipation in impacts and plastic flow, second order in elastic waves.
//   6 = strain-rate switch (2D): as 5, but the raw weight comes from the pair's strain rate,
//       s = h_ij max(|D_i|, |D_j|) / c_ij with D = sym(grad v) (raw gradients, zero for rigid
//       rotation): weight min(1, max(0, s/RECON_FACE_DV6_S0 - 1)), and 1 between bodies. Small-
//       amplitude waves and rigid rotation see the residual only; impacts and plastic flow the raw jump.
#define RECON_FACE_DV 7
#define RECON_FACE_DV5_WIDTH 0.05
#define RECON_FACE_DV6_S0 5e-3
//   7 = as 6 for the longitudinal channel (threshold RECON_FACE_DV6_S0; central, so it exerts no torque),
//       and a separate threshold RECON_FACE_DV7_S0T for the transverse channel (the non-central raw
//       friction that brakes rotation): the disc's elastic breathing can trigger 6 at low thresholds.
#define RECON_FACE_DV7_S0T 8e-3
// minmod of the residual (RECON_FACE_DV >= 3) and the raw jump, normal and tangential components:
// keeps the residual dissipation from turning anti-dissipative where the gradients are inexact.
#define RECON_RESID_MINMOD 0

// Numerical viscosity of the HLL damping in the GIZMO elastic stress flux (2026-10-05, branch num-visc).
// The HLL damping acts on a velocity jump in two channels: longitudinal (along r_ij) and transverse.
// Each channel's jump blends the RECONSTRUCTED face jump v_R - v_L (the HLLC's own limited,
// time-predicted face states at the pair midpoint, taken in the lab frame before the Riemann solver
// rotates them) with the BARE particle difference v_i - v_j:
//     jump = (1 - w) (v_R - v_L) + w (v_i - v_j),   w = NUM_VISC_L (longitudinal), NUM_VISC_T (transverse).
// w is a numerical viscosity and is named as one: the bare difference is O(h) for any smooth field,
// so every w > 0 makes the damping first order with strength w. w = 0 is the Godunov/MUSCL damping:
// zero for any linear field (rigid rotation included), O(h^2) in smooth flow, full strength at
// discontinuities and on grid-scale (zigzag) modes, which the reconstruction cannot follow.
// Replaces RECON_FACE_DV (which must be 0): DV 0 ~ (1, 1); DV 1 ~ (1, 0) and DV 2 ~ (0, 0), both
// without their frame error (DV 1/2 read the face states after the solver rotated them).
//   NUM_VISC 0 = off, RECON_FACE_DV decides (baseline, bit-identical)
//   NUM_VISC 1 = on, weights NUM_VISC_L / NUM_VISC_T in [0, 1]
#define NUM_VISC 0
#define NUM_VISC_L 1.
#define NUM_VISC_T 1.

// Artificial viscosity (2026-10-05, branch num-visc): quadratic (von Neumann-Richtmyer type)
// coefficient in the objective Morris pair form of MATERIAL_VISCOSITY,
//     eta_AV,i = rho_i (ART_VISC_C dx_i)^2 |D_i|,   dx_i = (m_i/rho_i)^(1/DIM),   |D| = sqrt(D:D),
// D = sym(grad v) from the unlimited gradients; pair coefficient (eta_AV,i + eta_AV,j)/2.
// It scales with h^2 and with the strain rate: O(h^2) in smooth flow (and ~amplitude^2, so negligible
// for small-amplitude waves at every resolution), strong in impacts, zero for rigid rotation (D = 0,
// and the pair form removes the spin). One formula everywhere, no switch. Energy via the work term;
// timestep dt <= 0.25 rho dx^2 / eta_AV.
//   0 = off (baseline, bit-identical)
//   1 = on, |D| of each particle (gradients; blind to pairing / zigzag modes)
//   2 = on, the pair's own objective strain rate s_ij = |(v_i - v_j) - W r_ij| / r_ij in place of |D|:
//       eta_ij = rho_ij (ART_VISC_C dx_ij)^2 s_ij (pair means of rho and dx^2). Equals |D r_hat| for a
//       smooth field, so still O(h^2) and zero for rigid rotation, but it also sees a pair approaching
//       its partner (pairing, clumping). The timestep uses the largest pair coefficient of a particle.
//   3 = LINEAR: eta_AV,i = ART_VISC_C rho_i c_i dx_i, c = P-wave speed sqrt((K + PWAVE_MU_FAC mu)/rho):
//       O(h) like the numerical viscosity of the bare HLL jump (eta_num ~ rho c r_ij), first order, but in
//       the objective pair form (no braking of rigid rotation), acting on the full pair difference.
//       Timestep dt <= 0.25 dx / (C c).
#define ART_VISC 0
#define ART_VISC_C 1.

// Nonlinearity switch for the GIZMO-like stabilisers (2D). phi_i = min(1, max(0, s_i/RECON_FACE_DV6_S0 - 1))
// with s_i = h_i |D_i| / c_i (raw gradients; D = sym(grad v) vanishes for rigid rotation), the same
// indicator as RECON_FACE_DV 6. Small-amplitude waves and rigid rotation have phi = 0, impacts and
// plastic flow phi = 1.
//   NL_SWITCH_TC   1: tensile-correction factor f scaled by max(phi_i, phi_j) inside a body (full
//                     across bodies): the Monaghan form damps on a perfect lattice too (f ~ 0.2), which
//                     gives a wave-speed error that does not converge.
//   NL_SWITCH_DIVV 1: explicit volume rate = phi * SPH divergence + (1 - phi) * matrix divergence (raw
//                     gradients): the SPH divergence spoils sound/P waves (c_err 8e-2 vs 1.4e-3).
#define NL_SWITCH_TC 1
#define NL_SWITCH_DIVV 1
//   NL_SWITCH_FACE 1: deviatoric stress flux on GIZMO's radial SPH area where phi = 1 (and across the
//                     contact), on the MFM face A_ij (as ELASTIC_FACE_MODE 2) where phi = 0. Needs
//                     ELASTIC_FACE_MODE 0.
#define NL_SWITCH_FACE 1
// amplitude part of the switch (1 = on): phi = max(phi_rate, clamp(s_amp/NL_AMP_S0 - 1, 0, 1)) with
// s_amp = max(|P|/K, sqrt(J2)/mu), the strain amplitude. The rate part alone switched the stabilisers off
// in quiet but finitely strained material (rod end faces after separation, stretched rings).
#define NL_SWITCH_AMP 0
#define NL_AMP_S0 1e-4
// memory time of the switch (0 = none): phi = max(phi_now, phi_old exp(-dt/TAU)), so a region keeps the
// stabilisers for ~TAU after the nonlinear phase (rod end faces relax after separation, where |D| -> 0)
#define NL_SWITCH_TAU 2.
// periodic images in the SPH divergence of the explicit volume (bug fix; 0 = old, images missing)
#define SPH_DIVV_GHOSTS 1

// Feed the deviatoric TRACTION in the GIZMO elastic stress flux with the
// spatially reconstructed face stress (Sij extrapolated to the shared face
// point xij, per side) instead of the raw cell-centered Sij[p]. Companion to
// RECON_FACE_DV: with both on, every dynamical input to the stress flux (stress
// AND velocity) is taken from the reconstructed face states, so the flux stops
// violating the scheme's own MUSCL reconstruction standard.
//
// NOTE the velocity is space+time centered (WijR carries the half-step), but the
// stress reconstruction here is SPATIAL only -- there is no dS/dt half-step in
// the Riemann reconstruction anywhere in this code (stress time integration lives
// in the separate RK2 integrateStressTensor). So this matches how the elastic
// HLLC path already reconstructs stress; full time-centering of S is a further step.
//
//   0 = off, raw cell-centered Sij[p] (baseline)
//   1 = per-side reconstructed face stress
#define RECON_FACE_TRACTION 0

/// Physical (material) viscosity in the GIZMO elastic stress flux -- objective
/// pair form (2026-10-02, mto). Per material, optional config key `eta` (dynamic
/// shear viscosity, default 0 => off for that material; Tillotson only).
///
/// Pair force on i (Morris 1997 form, same FNormT area as the elastic traction):
///     F_ij = -eta_ij * FNormT / r * ( (v_i - v_j) - W_ij (x_i - x_j) ),
/// eta_ij = harmonic mean of the two materials (as GIZMO's viscosity.h), and
/// W_ij = (W_i + W_j)/2 the mean SPIN (antisymmetric part of grad v) from the
/// UNLIMITED velocity gradients. Removing W r makes the force vanish exactly for
/// rigid rotation (objective, no angular-momentum drain on a spinning body),
/// while 2dx pairing/zigzag modes -- invisible to a gradient-based viscosity,
/// GIZMO's included -- are damped through the direct pair difference. That is
/// what RECON_FACE_DV 0 provides in its transverse channel (with the impedance
/// rho c_s instead of eta/r), minus the rigid-rotation braking.
/// Energy: the force enters cmag, so the work term carries it (see
/// ELASTIC_WORK_SIGN_FIX). Timestep: dt <= 0.25 rho L^2 / eta (GIZMO timestep.c).
///   0 = off (baseline, bit-identical)
///   1 = on (materials with eta > 0)
#define MATERIAL_VISCOSITY 1

/// Shear viscosity of a Tillotson material whose config block has no `eta` key
/// (only used where MATERIAL_VISCOSITY is on). 0.1 since 2026-10-03: the plasticity of the Taylor rods needs
/// some viscosity to stay well behaved (no clumping). Elastic-only waves and viscosity-free studies set
/// `eta 0` explicitly in their material block.
#define MATERIAL_ETA_DEFAULT 0.1

/// Kernel-moment normalisation of the MATERIAL_VISCOSITY pair coefficient
/// (2026-10-03, mto W4p/W4s). The pair force is eta * <operator> with
///     <operator> = M_i (lap v - div W),  M_i = 1/DIM sum_j V_j |W'(r,h_i) + W'(r,h_j)|/2 r_ij,
/// and M = 1 only for the continuum integral. On the wave lattice (h/dx 3.18,
/// 21 neighbours) M = 0.99187 at every resolution, so the viscous damping of
/// both modes is 0.8 % too weak and does not converge (measured 0.7 % at n64
/// and n128, fin1_g2_eta1). Dividing the coefficient by the pair mean
/// (M_i + M_j)/2 keeps the pair symmetric (momentum exact) and removes it.
/// At a free surface M drops to ~0.5; VISC_MOMENT_FLOOR bounds the division
/// there so a truncated support does not double the surface viscosity.
///   0 = off (baseline, bit-identical)
///   1 = coefficient / max(VISC_MOMENT_FLOOR, (M_i + M_j)/2)
#define VISC_MOMENT_NORM 1
#define VISC_MOMENT_FLOOR 0.95

/// Sign of the work term of the GIZMO elastic stress flux. The pair force on i
/// is +cmag (Fout momentum -= cmag), so the energy exchange must be +cmag . v_bar;
/// the port writes -cmag . v_bar (v_int = -(v_i+v_j)/2). Total energy is still
/// conserved (antisymmetric), but the kinetic/internal split is wrong: u is not
/// Galilean invariant and goes negative at moving surfaces (2026-09-26 report:
/// the flipped sign makes u invariant to round-off, mechanics unchanged).
///   0 = as ported (baseline)
///   1 = v_int = +(v_i+v_j)/2
#define ELASTIC_WORK_SIGN_FIX 0

// 2D elastic material model (mto P1 GIZMO comparison, 2026-10-02).
//   0 = plane strain without S_zz: dS = 2 mu (D - tr(D)/3 I) in the plane, P-wave
//       modulus K + 4/3 mu. The in-plane deviator is NOT traceless, so an in-plane
//       dilatation also loads S isotropically (tension under expansion, and it
//       enters J2, i.e. yields under pure dilatation).
//   1 = GIZMO's 2D model (solids/elastic_physics.c:131-132): dS = mu (grad v +
//       grad v^T) - mu tr(grad v) I, traceless in 2D; P-wave modulus K + mu
//       (eos.c:76 adds mu/rho to c^2), used for the HLL weights of the stress
//       flux, the HLLC wave speeds and the timestep as well.
#define STRESS_2D_GIZMO 1
#if STRESS_2D_GIZMO && DIM == 2
#define PWAVE_MU_FAC 1.0
#else
#define PWAVE_MU_FAC (4.0/3.0)
#endif

// GIZMO's star-pressure floor (hydro/reimann.h:547, get_wavespeeds_and_pressure_star):
// in the dummy-shifted frame the HLLC star pressure is clipped at 0 (and S* set to
// 0), so the Riemann channel never pulls harder than the dummy shift, i.e. than the
// smaller of the two particle pressures: separating faces get no velocity-driven
// tension. 0 = unclipped (default), 1 = GIZMO floor.
#define RIEMANN_PSTAR_FLOOR 1

// Velocity gradient for the stress rate in the first stress half-kick.
//   0 = the slope-LIMITED gradient (it is computed after slopeLimiter)
//   1 = the unlimited gradient, as GIZMO (gradients.c:1214: "needs to be done
//       before slope-limiting to use full velocity gradient information")
// The second half-kick uses freshly computed, unlimited gradients either way.
#define STRESS_RATE_RAW_GRAD 1

// Sign of the Jaumann rotation term in the stress rate (2D).
//   0 = objective (co-rotational) rate dS/dt += W S - S W, W = (grad v - grad v^T)/2
//   1 = GIZMO's sign (elastic_physics.c:133: dS += S R - R S with the same R),
//       which rotates the stress against the material rotation. Comparison only.
#define JAUMANN_SIGN_GIZMO 0

// HLLC wave-speed estimate of the x-split elastic solver.
//   0 = per side: S_L = v_L - a_L, S_R = v_R + a_R
//   1 = GIZMO (reimann.h:541-542, Davis): S_L = min(v_L,v_R) - max(a_L,a_R),
//       S_R = max(v_L,v_R) + max(a_L,a_R)
#define HLLC_WAVESPEED_GIZMO 1

#define STORE_RAW_VGRAD (MATERIAL_VISCOSITY || STRESS_RATE_RAW_GRAD || RECON_FACE_DV >= 3 || GNOG_PENALTY || NL_SWITCH_TC || NL_SWITCH_DIVV || NL_SWITCH_FACE || ART_VISC)

// Limited longitudinal dissipation of the GIZMO stress flux (option O5, Zhang-Hu-
// Adams / Rezavand form): the impedance speeds of wt_r are replaced by
// beta_K = min(3 max(-dv.rhat, 0), c_K), i.e. full P-wave dissipation only for
// pairs approaching at >= c/3, none for receding pairs (no viscous cohesion across
// an opening contact) and ~none for linear waves of small amplitude (removes the
// O(h) wave damping). 0 = off (wt_r with c_K), 1 = on.
#define WTR_LIMITER 0

// Contact without adhesion between different bodies (objective alternative to
// JAUMANN_SIGN_GIZMO for the P1 rebound). Bodies are the connected components of
// the first force pass (link: r_ij < 1.5 x the smaller nearest-neighbour distance),
// so T1/T2/P1 have two bodies and the disc and the wave boxes one. For a pair of
// different bodies an ATTRACTIVE normal component of the total pair force
// (Riemann + stress flux) is removed; compression and the tangential part are kept.
// The energy flux is corrected by the removed momentum flux times the mean pair
// velocity. 0 = off, 1 = on.
#define CONTACT_NO_TENSION 1

// Split-face form of the deviatoric stress divergence in the GIZMO stress flux
// (option O3, Mars ELASTIC_FACE_MODE 3): pair force S~_i a_ij - S~_j a_ji with the
// MFM half-faces a_ij = V_i psiTilde_j(x_i) (A_ij = a_ij - a_ji). Exactly torque-free
// for symmetric S, because sum_j (x_j - x_i) (x) psiTilde_j(x_i) = I, PROVIDED the
// tensile correction is per particle (option O4): S~_i = S_i with its tensile
// principal components damped once by f_i = max_j f_ij (TC_TOTAL_STRESS_GATE
// applies). The HLL dissipation terms are unchanged. 2D, non-periodic.
//   0 = off (eigen-projection on the radial SPH area, per-pair TC), 1 = on
#define SPLIT_FACE_STRESS 0

// GNOG hourglass penalty (Zhang, Wu, Lourenco, Hu, arXiv 2409.11474), as implemented in
// SPHinXsys (ShearStressRelaxationHourglassControl2ndHalf, with the GNOG-A rotation): a
// per-particle force F^h_i, integrated in time, against the residual of the pair velocity
// difference with respect to the linear prediction (raw, unlimited gradients):
//   v^_ij  = v_ij - 1/2 (grad v_i + grad v_j) r_ij,   r_ij = x_i - x_j, e_ij = r_ij/|r_ij|
//   dF^h_i = sum_j xi_ij G (v^_ij . e_ij) e_ij dW/dr V_i V_j dt / |r_ij|   (same body only)
// F^h_i is rotated each step by the particle's spin (Rodrigues). v^ vanishes for linear
// velocity fields (rigid rotation); each increment is a central pair force (momentum
// exact, torque-free). xi_ij = (xi_i + xi_j)/2 with xi_i = GNOG_XI (elastic) or
// GNOG_XI_PLASTIC times the radial-return factor min(1, Y/sigma_vM,trial) (paper: 4 and
// 0.2). The work F^h.v enters the total energy, the internal energy is not touched.
//   0 = off
//   1 = as above: one accumulated force per particle (SPHinXsys). NOT momentum-conserving per body
//       (the per-particle rotation and the persistent sums; report 2026-10-03-tipsntricks-abgleich:
//       free blocks keep accelerating, T2 energy +1.8 %). Kept for reproducibility. 2D, non-periodic.
//   2 = conservative: one history per unordered pair, u^_ij = int gamma_ij v^_ij dt (rotated by the
//       pair-mean spin), force F_i = xi G (u^_ij . e_ij) e_ij dW/dr V_i V_j / |r_ij| = -F_j, recomputed
//       every step (exact momentum, central -> exact AM), energy flux F.(v_i+v_j)/2 like every pair
//       flux (exact total energy). Histories of pairs that left the stencil are dropped (the force
//       vanishes there anyway with dW). 2D, periodic allowed.
#define GNOG_PENALTY 2
#define GNOG_XI 4.0
#define GNOG_XI_PLASTIC 0.2

// Gate of the deviatoric tensile correction (2D GIZMO stress flux in Particles.cpp; the embedded-stress
// Riemann solver, TENSILE_CORRECTION_3, in Riemann.cpp, 2D and 3D). The pressure part of the TC
// (TENSILE_CORRECTION_1, both Riemann solvers) acts on tension only by construction.
//   0 = GIZMO: every positive eigenvalue l_k of the DEVIATORIC stress S is damped
//       by (1 - f), also where the total principal stress l_k - P is compressive
//       (measured on P1: 6-16 % of the correction acts under compression)
//   1 = damp only the part of l_k that is tensile in the TOTAL stress (P = pressure, > 0 compressive):
//       l_k -> l_k - f min(l_k, max(0, l_k - P))       STANDARD since 2026-10-03 (user)
//   2 = 1 for pairs of different bodies (contact: no particle crossing), 0 inside a body (lets
//       the TC separate pairs collapsed under compression; the no-adhesion contact leaves the
//       rod end faces in residual compression, where 1 never fires). GIZMO stress flux only;
//       in Riemann.cpp it acts as 1.
#define TC_TOTAL_STRESS_GATE 2

// Angular-momentum TORQUE BUDGET instrumentation (WP0 of the 20260727 plan).
//
// Splits the spurious pair torque that drives the L_z drain into the four
// channels that can produce it, so the loss can be attributed to a mechanism
// instead of inferred from a parameter sweep:
//
//   tqHLLC : pressure / Riemann flux through the MFM face A_ij
//   tqDev  : deviatoric traction (GIZMO eigenvalue loop)
//   tqWtT  : transverse (shear-wave) HLL dissipation
//   tqWtR  : longitudinal HLL dissipation
//
// tqWtR is the SELF-CHECK: that term is built as a multiple of (dx,dy), i.e.
// exactly parallel to r_ij, so its pair torque r_ij x F_ij vanishes identically.
// If the reported tqWtR is not at round-off, the instrumentation is wrong.
// Likewise tqHLLC+tqDev+tqWtT+tqWtR must reproduce tqF particle by particle
// (up to the round-off of recomputing the channels separately).
//
// Also computes consResid[i], the anisotropy of the elastic stress area
// (|| Ghat_i - I ||_F). The torque-free condition for a pair force
// F_ij = sigma . A_ij with symmetric sigma is M = sum r_ij (x) A_ij ~ I; the MFM
// face satisfies it exactly, the raw SPH area of the elastic flux does not, and
// consResid measures by how much. Doubles as the free-surface classifier.
//
// Diagnostic only: the flux math is left byte-for-byte untouched, the channels
// are recomputed alongside it. 0 = off (production), 1 = on.
#define TORQUE_DIAG 0

// tqF (per-particle residual torque accumulator) is needed by the AM_METHOD 4-6
// repair strategies AND by the torque budget above.
#define TORQUE_ACCUM (AM_TORQUE_TRACK || TORQUE_DIAG)

// Contact-force budget (diagnostic, mto P1 welding study 2026-10-02).
//
// Splits the pair force that crosses between the two bodies of a two-body test
// (body = sign of x at the first force evaluation; T1, T2 and P1 are all split
// by x = 0) into channels and writes one line per step to contact_forces.csv in
// the working directory, as the x-force on the RIGHT body:
//   riem       isotropic HLLC momentum flux (incl. the dummy-pressure shift)
//   riem_tc    of which the tensile-correction repulsion f*pDummy
//   riem_floor what GIZMO's star-pressure floor (P_M >= 0 in the shifted frame,
//              reimann.h get_wavespeeds_and_pressure_star) would ADD
//   dev        deviatoric traction (with the tensile correction)
//   tc         of which the tensile correction
//   hll_r      longitudinal HLL dissipation (incl. the approach damping)
//   hll_t      transverse HLL dissipation
//   damp       of which the approach damping (TENSILE_CORRECTION_DAMP)
//   visc       MATERIAL_VISCOSITY
// plus, over ALL pairs, the summed magnitude of the tensile correction and of
// the approach damping, split by whether the total stress is tensile or
// compressive (does the tensile correction act only under tension?).
// Read-only: the flux math is untouched. 0 = off.
#define CONTACT_FORCE_DIAG 0

// Variationally consistent face AREA for the GIZMO elastic stress flux
// (WP2 of the 20260727 plan).
//
// The deviatoric flux integrates over the RAW SPH area FVec = FNormT * rhat
// (addGizmoElasticStressFlux), which carries no consistency guarantee. T1
// measures the consequence exactly: the net torque under a uniform stress S is
//     sum tqDev = S_xy (M_xx - M_yy) - (S_xx - S_yy) M_xy,
//     M = sum_pairs FNormT r (rhat (x) rhat) = 0.5 * sum_i G_i,
// reproduced to ratio 1.000 by torqueFromG(). It vanishes only when the GLOBAL
// moment M is isotropic -- true by symmetry on a pristine lattice, false as soon
// as the configuration is disordered, which is exactly the regime where L_z is
// actually lost. The MFM face A_ij has no such problem: B = E^-1 makes
// M = -2 V_tot I identically, and T1 measures its channel at 1e-17 relative even
// under 20% position jitter.
//
//   0 = raw SPH radial area (baseline, bit-identical to pre-WP2 code)
//
//   1 = pair-symmetric Bonet-Lok correction. Build the per-particle consistency
//       matrix G_i = sum_j FNormT_ij r_ij (rhat (x) rhat) and correct the area
//       with B_i = V_i G_i^-1, applied pair-symmetrically as 0.5*(B_i + B_j).
//       The symmetrisation is the point: the standard objection to gradient
//       correction is that B_i != B_j destroys pairwise antisymmetry and with it
//       exact linear momentum, and 0.5*(B_i+B_j) is invariant under i<->j, so
//       A_ij = -A_ji survives and p stays machine-exact (checkFluxSymmetry
//       confirms). Only the TRACTION is corrected -- wt_r/wt_t keep the scalar
//       FNormT, because making the acoustic impedances anisotropic is a separate
//       question and because a radial wt_r is torque-free by construction.
//
//   2 = use the MFM face A_ij itself. Isotropy is then exact rather than
//       approximate. Riskier: GIZMO chose the radial area specifically to avoid
//       the tensile instability, so the ring is the regression gate.
#define ELASTIC_FACE_MODE 0
// Reject the correction and fall back to B_i = I when G_i is this badly
// conditioned (one-sided surface stencils). Counted and reported.
#define ELASTIC_FACE_COND_MAX 1.0e3

// Angular-momentum SELF-TESTS (WP-T of the 20260727 plan). These overwrite the
// particle state with an analytically known field on step 0 and measure a single
// quantity, so a hypothesis that would otherwise cost a resolution campaign is
// settled in seconds. Run them on any 2D elastic IC that has a free surface --
// the standard offset-blocks IC is a good default (bulk + flat faces + corners).
//
//   0 = off, normal run.
//
//   1 = T1, UNIFORM-STRESS TORQUE. v := 0 and S := const, then one flux pass.
//       Reports the per-channel torque and exits. The torque-free condition
//       M ~ I predicts: the HLLC/MFM-face channel is machine zero even on a
//       strongly anisotropic surface stencil, while the elastic channel -- which
//       uses the raw SPH area, not the MFM face -- is O(1) * anisotropy.
//       THIS IS THE GATE FOR WP2: the corrected area must drive the elastic
//       channel to machine zero too. Requires TORQUE_DIAG 1.
//
//   2 = T2, RIGID-ROTATION RECONSTRUCTION. v := Omega x r, then one
//       reconstruction pass. Reports the reconstructed face velocity jump
//       against the raw particle difference and exits. RECON_FACE_DV's whole
//       justification is that a linear velocity field reconstructs to zero jump;
//       that holds with the limiter off and FAILS once the limiter clips vx and
//       vy by different alphas. THIS IS THE GATE FOR WP1: the jump must be at
//       round-off with SLOPE_LIMITING on as well as off.
//
//   3 = T3, PRE-STRESSED ROTATION. v := Omega x r AND S := const, then a NORMAL
//       run (no early exit). The existing rotating-disc null test is stress-FREE
//       and therefore blind to an inconsistency between the velocity gradient
//       the Jaumann terms in integrateStressTensor use and the one the flux pass
//       uses; such an inconsistency would be a torque source upstream of
//       everything else. Here S must simply rotate rigidly and L_z must hold.
#define AM_SELFTEST 0
#define AM_SELFTEST_OMEGA 0.1
// Deterministic position jitter for T1, as a fraction of the mean particle
// spacing (0 = off). ESSENTIAL for a meaningful T1: the net deviatoric torque is
// controlled by the GLOBAL moment M = sum_pairs FNormT r rhat (x) rhat, and on a
// symmetric configuration (two square blocks on a square lattice) M is isotropic
// BY SYMMETRY -- the O(1) per-pair torques cancel exactly and the test reports a
// false pass, even though every surface particle has a strongly anisotropic
// stencil. Jitter breaks that cancellation, which is also what the real problem
// does: L_z is lost during DISORDERED contact, not on a pristine lattice.
#define AM_SELFTEST_JITTER 0.0
// Uniform pressure forced in T1 (see compPressure). Without this the jitter
// perturbs rho and hence P, and the "uniform stress" premise the test rests on
// is violated -- the HLLC channel then reads a real torque that has nothing to do
// with the face geometry under test.
#define AM_SELFTEST_P 1.0
// Constant deviatoric stress for T1/T3 (traceless in 2D).
#define AM_SELFTEST_SXX  1.0
#define AM_SELFTEST_SXY  0.5
#define AM_SELFTEST_SYY -1.0

/** maximum interactions with ghost particles
 *  ignored when `PERIODIC_BOUNDARIES` is not set
**/
#define MAX_NUM_GHOST_INTERACTIONS 150

/// flag for slope limiting, 0: no slope limiting
#define SLOPE_LIMITING 1

/// slope limiting parameter, ignored when `SLOPE_LIMITING` is false
#define BETA 4.             // TODO: move to config
/// Relative round-off tolerance of the slope limiter's denominators (0 = absolute 1e-30, the old
/// behaviour). With 0 the limiter breaks mirror symmetry in the first step of P1 (a surface particle
/// keeps its density gradient, its mirror partner's is zeroed: exact 0 vs 1e-17 in the denominator).
#define SLOPE_LIMITER_EPS 1e-10

/// use pairwise limiter
#define PAIRWISE_LIMITER 1
#define PSI_1 .5            // TODO: move to config
#define PSI_2 .25           // TODO: move to config

/// meshless finite mass method instead of meshless finite volume
#define MESHLESS_FINITE_MASS 1

// --- Adiabatic (entropic) face correction for MFM -------------------------
//
// Hopkins (2015), MNRAS 450, 53 -- "GIZMO: a new class of accurate, mesh-free
// hydrodynamic simulation methods". Port of the MFM "face correction for
// adiabatic flows", GIZMO master, hydro/hydro_core_meshless.h:335-375.
//
// WHAT IT DOES. In MFM there is no mass flux, so the whole energy exchange
// across a face is the pressure work. The Godunov form takes it from the
// Riemann contact wave,
//
//     du_old = A_ij * P_M * (S_M + v_frame . nhat)                        (1)
//
// which is correct at a shock but is NOT entropy-conserving in smooth flow:
// S_M carries the solver's dissipation, so every face leaks a little entropy
// even where the flow is perfectly adiabatic. The alternative is the
// thermodynamic PdV work, evaluated from the actual rate of change of the
// cell volume (an SPH-style estimator), which is adiabatic by construction:
//
//     du_new = 1/2 * ( PdV_i - PdV_j + A_ij * P_M * (v_i.nhat + v_j.nhat) ) (2)
//     PdV_k  = (dW/dr)(r, h_k) * V_k^2 * P_M * (dx . dv)/r
//
// The correction swaps (1) for (2) -- `Fij[1] += du_new - du_old` -- but only
// where the face is not doing shock work, i.e. where the contact speed is
// small compared to the sound speed:
//
//     |S_M| / min(c_i, c_j) < ADIABATIC_FACE_EPS_BIG
//
// Both (1) and (2) are antisymmetric under i <-> j, so total energy stays
// conserved to machine precision whichever is used; only the entropy
// (i.e. the split between internal and kinetic) changes.
//
// INTERMEDIATE BAND. Between ADIABATIC_FACE_EPS_SMALL and ..._EPS_BIG, GIZMO
// keeps the entropic form only if it does not move energy the "wrong way"
// across the face relative to the entropic ordering P_i/rho_i vs P_j/rho_j
// (hydro_core_meshless.h:354-370); below EPS_SMALL it is used unconditionally.
// That guard is reproduced here.
//
// NOTE -- PAPER vs CODE. This correction is NOT in the Hopkins (2015) paper.
// The paper's Appendix D describes a *different*, energy-ratio-based
// energy-entropy switch for extremely high Mach numbers (the Bryan et al. 1995
// dual-energy formalism, with alpha_kin / alpha_grav ~ 1e-3), which in GIZMO is
// the separate, OFF-by-default ENERGY_ENTROPY_SWITCH_IS_ACTIVE. The correction
// implemented here is switched on |S_M|/c_s, is ON by default in GIZMO's
// MFM/MFV, and is documented only in the code and in the GIZMO user guide
// (under FORCE_ENTROPIC_EOS_BELOW, which retunes ..._EPS_SMALL).
//
// DEVIATIONS FROM GIZMO, all recorded because they matter for interpreting a
// convergence test:
//   * GIZMO weights PdV_k by the grad-h factor DhsmlNgbFactor; the demonstrator
//     has no such array (its smoothing length solve carries no Omega_grad-h
//     correction), so it is taken as 1. It is O(1) and -> 1 for a smooth,
//     well-sampled field, which is the regime this correction targets.
//   * GIZMO's literal PdV_k carries no m_j/m_i ratio, which makes (2) exact
//     only for equal particle masses. The wave tests here are equal-mass.
//   * GIZMO additionally disables the correction on a bad face-closure error
//     (leak_vs_tol) and forces it on at a dangerous condition number; the
//     demonstrator's condition-number gate is applied, the closure gate is not.
//
// GIZMO restricts this to ideal gases -- "below is defined for adiabatic ideal
// fluids, don't use for materials" -- and compiles it out for EOS_TILLOTSON and
// EOS_ELASTIC. The same restriction is enforced below: the correction refuses to
// build with a material EOS or with live solid dynamics.
//
//   0 = off (baseline, bit-identical to pre-correction code)
//   1 = on
#define ADIABATIC_FACE_CORRECTION 0

// ---------------------------------------------------------------------------
// Time-order study 2026-09-24: three single-purpose switches, all default 0/1 =
// baseline, added to find the O(dt) term that survives ADIABATIC_FACE_CORRECTION
// (2026-09-24-fixed-dt-leiter-adiabatische-korrektur: p stays 1.0 over six
// configurations, only the prefactor drops 5.9 -> 2.1).
// ---------------------------------------------------------------------------

/// Drop the advective energy remainder from the MFM star-state energy flux.
///
/// xSplitElasticHLLC builds the MFM flux by zeroing BOTH advective pieces of the
/// mass line -- rho_K v_K (never written) and rho_K S_K (starfac - 1) -- but the
/// energy line drops only rho_K v_K e_K and keeps rho_K S_K (starfac - 1) e_K,
/// because that statement sits behind #endif // MESHLESS_FINITE_MASS.  Mass, and
/// the three momentum components, each got an MFM variant; the energy line alone
/// did not.
///
/// The remainder is  rho_K S_K (S* - v~_K)/(S_K - S*) e_K, which vanishes exactly
/// under the zero-mass-flux condition S* = v~_K.  Decisive check: set S* = 0.
/// The face energy flux must then be analytically 0 (star state,
/// operator_splitting.tex Sec. 4.4); the code gives -rho_K v_K e_K, and with the
/// term removed it is exactly 0.  GIZMO has no such term at all: in the MFM frame
/// it sets Fluxes.p = P_M * S_M outright (reimann.h:411-412).
///
///   0 = baseline (remainder kept, as of 2026-09-24)
///   1 = remainder dropped under MESHLESS_FINITE_MASS
#define MFM_ENERGY_STARFIX 0

/// Second-order (trapezoidal) drift of the particle positions.
///
/// The step is kick-then-drift: updateState() applies the full-dt kick from
/// fluxes that ARE time-centred (the dt/2 MUSCL predictor in
/// compRiemannStatesLR), so v^{n+1} is second order -- but the drift then uses a
/// single endpoint velocity, x^{n+1} = x^n + dt v^{n+1} in moveParticles() (the
/// ELASTIC path) and x^{n+1} = x^n + dt v^n in updateStateAndPosition() (the
/// non-ELASTIC path).  Either way the position is first order: the local error
/// dt (v^{n+1} - v^{n+1/2}) = dt^2 a / 2 accumulates over T/dt steps to O(dt).
///
/// In a Lagrangian scheme the positions ARE the mesh, so an O(dt) position error
/// is an O(dt) volume -- hence density -- error, while v and u stay second order.
/// That asymmetry is the measured signature: L_1(v_x) fits p = 1.995 while the
/// density damping 1 - a_rho fits p = 1.0, and the coefficient is independent of
/// both the limiter and the resolution.
///
/// With x^{n+1} = x^n + dt (v^n + v^{n+1})/2 the update is velocity-Verlet to
/// O(dt^3).  The trapezoidal form is already in updateStateAndPosition(), commented
/// out with "stay consistent with choice of effective frame" -- see
/// DRIFT_KDK_FRAME for that consistency question, which this flag leaves alone.
///
///   0 = baseline (single-endpoint drift)
///   1 = trapezoidal drift, x += dt (v^n + v^{n+1})/2
#define DRIFT_KDK 1

/// MUSCL half-step (time) predictor in the face reconstruction.
///
/// Diagnostic only: 0 passes dt = 0 into compRiemannStatesLR, which is the single
/// place dt enters the reconstruction, so the states are reconstructed in space
/// but not centred in time.  That must be first order in dt.  It is here to
/// bracket the ladder: if the baseline coefficient a is unchanged with the
/// predictor off, the predictor is not what sets the O(dt) term.
///
///   1 = baseline (predict to t^n + dt/2)
///   0 = no time centring (deliberately first order)
#define RECON_TIME_PRED 1

/// Centre the kick interval on the flux evaluation (drift-kick-drift).
///
/// The loop builds geometry, states and faces at t^n and then kicks over
/// [t^n, t^n + dt].  t^n is the START of that interval, not its middle, so the
/// kick is off-centre by dt/2 -- an O(dt) error in the velocity even once the
/// drift itself is Verlet-correct (DRIFT_KDK).  That is the observed split after
/// the 2026-09-24 wave-1 sweep: with DRIFT_KDK the density amplitude is
/// dt-independent to 0.8 % over a factor 85 in dt, while the PHASE error stays
/// first order.  On the harmonic oscillator the same scheme gives exactly
/// velocity-Verlet for x and plain explicit Euler for v.
///
/// GIZMO has no such off-centring and needs no predictor for it: it drifts to t,
/// builds everything there, and applies +-dt/2 kicks around it
/// (DO_HALFSTEP_FOR_MESHLESS_METHODS is NOT enabled by default,
/// Template_Config.sh:530).  TIME_DKD does the same here -- drift dt/2, rebuild
/// the neighbour state, evaluate, kick the full dt, drift dt/2 -- so the states
/// are AT the midpoint rather than extrapolated to it and RECON_TIME_PRED
/// becomes unnecessary.
///
/// Cost: a second neighbour + density pass per step (~1.6x).  The VARIABLE_SML
/// retry loop is not repeated in the rebuild; nnsRadius was sized for this step
/// and half a drift does not move h_max appreciably in a smooth flow.
///
/// Combine with DRIFT_KDK 0: the two half drifts already sum to the Verlet
/// position update, and DRIFT_KDK would average the velocities a second time.
///
///   0 = baseline (evaluate at t^n, kick forward a full dt)
///   1 = drift dt/2, evaluate, kick dt, drift dt/2
#define TIME_DKD 0

/// Move the face frame velocity to the half step along with the states.
///
/// compRiemannStatesLR boosts the reconstructed states into the face frame
/// (W := v - vFrame) and THEN predicts them half a step forward, but never moves
/// vFrame itself -- the frame stays at t^n while the states are at t^n + dt/2.
/// The code has carried the TODO for this since the predictor was written:
/// "center vL and vR and update vFrame (compare to GIZMO code
/// hydro_core_meshless.h:178ff)".
///
/// GIZMO's answer (hydro_core_meshless.h:163-177, DO_HALFSTEP_FOR_MESHLESS_METHODS)
/// splits the predicted velocity change dv_k in two: the ANTIsymmetric half
/// (dv_i - dv_j)/2 stays on the states, the SYMMETRIC half (dv_i + dv_j)/2 moves
/// into v_frame.  Since our states are boosted before the prediction, that is the
/// same as subtracting (dv_i + dv_j)/2 from both and adding it to vFrame -- which
/// is what the implementation does, so the lab-frame prediction is unchanged and
/// only the frame moves.  Only the acceleration part of dv enters: the advective
/// terms transport the state, they do not move the frame.
///
/// Motivation: with DRIFT_KDK the O(dt) DISSIPATION is gone (2026-09-24 sweep,
/// a: 5.885 -> -0.198) but the residual O(dt) error is DISPERSIVE -- the phase
/// error at the coarsest rung grows 2.40e-3 -> 5.87e-3.  A face that moves at the
/// t^n velocity while the particles move at the centred one is an O(dt) advection
/// error, i.e. exactly a phase error.
///
/// Note GIZMO does NOT enable its half-step by default (Template_Config.sh:530
/// lists DO_HALFSTEP_FOR_MESHLESS_METHODS commented out, and calls the effect
/// "very, very weak") -- so GIZMO's second-order-in-time behaviour comes from the
/// leapfrog, not from the reconstruction.  RECON_TIME_PRED 0 together with
/// DRIFT_KDK 1 is that configuration.
///
/// The DIM == 3 branch takes the pressure acceleration only; the demonstrator's
/// own 3D velocity prediction is incomplete there ("TODO: update below for 3D").
///
///   0 = baseline (frame at t^n)
///   1 = frame advanced to t^n + dt/2
#define FRAME_HALFSTEP 0

/// |S_M|/c_s above which the face is treated as doing shock work and the
/// Godunov energy flux is kept unchanged. GIZMO: epsilon_entropic_eos_big,
/// "can be anything from (small number = more diffusive, less accurate entropy
/// conservation) to ~1.1-1.3 (least diffusive, most noisy)".
#define ADIABATIC_FACE_EPS_BIG 0.5

/// |S_M|/c_s below which the entropic form is used unconditionally. GIZMO:
/// epsilon_entropic_eos_small; should be << ADIABATIC_FACE_EPS_BIG.
#define ADIABATIC_FACE_EPS_SMALL 1.e-3

/// Condition number above which the entropic form is forced on regardless of
/// the intermediate-band guard. GIZMO: CONDITION_NUMBER_DANGER.
#define ADIABATIC_FACE_COND_DANGER 1.0e3

/// Report, once per step, the fraction of faces the correction fired on and the
/// largest relative gap |du_new - du_old| / |du_old| among them. In smooth flow
/// the two forms estimate the same physical work, so the gap must be small and
/// must shrink under resolution; a gap of order 2 means the sign convention of
/// (2) is inverted relative to this code's "flux leaves i" convention. Used to
/// validate the port. 0 = off.
#define ADIABATIC_FACE_DIAG 0

// Consistency guards live further down, once USE_HLLC is defined.
//
// ELASTIC is deliberately NOT forbidden: the ideal-gas soundwave build runs the
// ELASTIC branch with shear modulus mu = 0 and S = 0 (see the testcase's
// makeBuilds.sh), which is pure hydro. With a non-zero mu the deviatoric
// traction carries work that Eq. (2) does not model, so the correction must not
// be used there -- there is no compile-time check for that, since mu lives in
// the runtime material table.

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

#if TIME_DKD && DRIFT_KDK
#error "TIME_DKD already yields the Verlet position update through its two half drifts: set DRIFT_KDK 0."
#endif

// ADIABATIC_FACE_CORRECTION consistency guards (same reason: needs USE_HLLC)
#if ADIABATIC_FACE_CORRECTION && !MESHLESS_FINITE_MASS
#error "ADIABATIC_FACE_CORRECTION is an MFM correction: it needs MESHLESS_FINITE_MASS == 1."
#endif
#if ADIABATIC_FACE_CORRECTION && EOS != 0
#error "ADIABATIC_FACE_CORRECTION is for ideal gases only (GIZMO excludes EOS_TILLOTSON/EOS_ELASTIC): set EOS 0."
#endif
#if ADIABATIC_FACE_CORRECTION && !(USE_HLLC && ELASTIC)
// P_M and S_M are read off Riemann::HLLCFlux's x-split solver output, which only
// exists on the ELASTIC HLLC branch. That is the branch the ideal-gas soundwave
// builds actually use (ELASTIC 1, mu = 0); the plain solveHLLC1 path would need
// its own S_M hand-back and is deliberately not wired up rather than guessed at.
#error "ADIABATIC_FACE_CORRECTION needs USE_HLLC == 1 && ELASTIC == 1 (x-split HLLC path)."
#endif
#if RECON_FACE_DV && !GIZMO_ELASTIC_FLUX
#error "RECON_FACE_DV only affects the GIZMO elastic stress flux; enable GIZMO_ELASTIC_FLUX."
#endif
#if RECON_FACE_DV && AM_NO_TDISS
#error "RECON_FACE_DV and AM_METHOD 3/7 (wt_t=0) both target the transverse dissipation; pick one."
#endif
#if RECON_FACE_TRACTION && !GIZMO_ELASTIC_FLUX
#error "RECON_FACE_TRACTION only affects the GIZMO elastic stress flux; enable GIZMO_ELASTIC_FLUX."
#endif
// TORQUE_DIAG consistency guards
#if TORQUE_DIAG && DIM != 2
#error "TORQUE_DIAG reports the z-torque only; implemented for DIM == 2."
#endif
#if TORQUE_DIAG && !(USE_HLLC && GIZMO_ELASTIC_FLUX)
#error "TORQUE_DIAG splits the HLLC + GIZMO elastic channels; both must be active."
#endif
// (the ENFORCE_FLUX_SYM guard for TORQUE_DIAG lives below, after that flag is defined)
#if SPLIT_FACE_STRESS && (DIM != 2 || PERIODIC_BOUNDARIES || !GIZMO_ELASTIC_FLUX || !TENSILE_CORRECTION)
#error "SPLIT_FACE_STRESS: DIM 2, non-periodic, GIZMO elastic flux with TENSILE_CORRECTION only."
#endif
#if NL_SWITCH_FACE && (DIM != 2 || ELASTIC_FACE_MODE != 0)
#error "NL_SWITCH_FACE: DIM 2 with ELASTIC_FACE_MODE 0."
#endif
#if GNOG_PENALTY && (DIM != 2 || !ELASTIC || (GNOG_PENALTY == 1 && PERIODIC_BOUNDARIES))
#error "GNOG_PENALTY: DIM 2 and ELASTIC only; mode 1 also non-periodic."
#endif
#if NUM_VISC && (RECON_FACE_DV || !GIZMO_ELASTIC_FLUX || !USE_HLLC || DIM != 2 || AM_NO_TDISS || TORQUE_DIAG || CONTACT_FORCE_DIAG)
#error "NUM_VISC replaces RECON_FACE_DV (set it to 0); DIM 2, HLLC + GIZMO elastic flux, no AM_NO_TDISS / TORQUE_DIAG / CONTACT_FORCE_DIAG."
#endif
#if ART_VISC && !(GIZMO_ELASTIC_FLUX && DIM == 2 && (EOS == 1 || EOS == 2))
#error "ART_VISC lives in the 2D GIZMO elastic stress flux (EOS 1/2)."
#endif
#if RECON_FACE_DV >= 3 && DIM != 2
#error "RECON_FACE_DV 3/4 (residual jump) is implemented for DIM == 2 only."
#endif
#if (WTR_LIMITER || CONTACT_NO_TENSION) && !(DIM == 2 && USE_HLLC && GIZMO_ELASTIC_FLUX)
#error "WTR_LIMITER / CONTACT_NO_TENSION: DIM 2 with the HLLC + GIZMO elastic flux only."
#endif
#if CONTACT_FORCE_DIAG && !(DIM == 2 && USE_HLLC && GIZMO_ELASTIC_FLUX && !PERIODIC_BOUNDARIES)
#error "CONTACT_FORCE_DIAG: DIM 2, HLLC + GIZMO elastic flux, non-periodic only."
#endif
// ELASTIC_FACE_MODE consistency guards
#if ELASTIC_FACE_MODE && !GIZMO_ELASTIC_FLUX
#error "ELASTIC_FACE_MODE only affects the GIZMO elastic stress flux; enable GIZMO_ELASTIC_FLUX."
#endif
#if ELASTIC_FACE_MODE && DIM != 2
#error "ELASTIC_FACE_MODE is implemented for DIM == 2 only."
#endif
#if ELASTIC_FACE_MODE > 2
#error "ELASTIC_FACE_MODE must be 0 (raw), 1 (Bonet-Lok) or 2 (MFM face)."
#endif
#if ELASTIC_FACE_MODE == 2 && RECON_FACE_TRACTION
#error "ELASTIC_FACE_MODE 2 replaces the area, RECON_FACE_TRACTION the stress; combination untested."
#endif
// AM_SELFTEST consistency guards
#if AM_SELFTEST && !(DIM == 2 && ELASTIC)
#error "AM_SELFTEST seeds a 2D deviatoric stress state; needs DIM == 2 and ELASTIC 1."
#endif
#if AM_SELFTEST == 1 && !TORQUE_DIAG
#error "AM_SELFTEST 1 (T1) reports the per-channel torque; enable TORQUE_DIAG."
#endif
#if RECON_FACE_TRACTION && DIM != 2
#error "RECON_FACE_TRACTION is currently implemented for DIM == 2 only."
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

// ELASTIC_HLLC_EP consistency guards
#if ELASTIC_HLLC_EP && !(USE_HLLC && ELASTIC)
#error "ELASTIC_HLLC_EP requires USE_HLLC == 1 and ELASTIC == 1."
#endif
#if ELASTIC_HLLC_EP && GIZMO_ELASTIC_FLUX
#error "ELASTIC_HLLC_EP carries the deviatoric stress inside the Riemann solver; disable GIZMO_ELASTIC_FLUX."
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

#if TORQUE_DIAG && ENFORCE_FLUX_SYM != 1
#error "TORQUE_DIAG accumulates each pair once at compute time; needs ENFORCE_FLUX_SYM 1."
#endif
#if CONTACT_FORCE_DIAG && ENFORCE_FLUX_SYM != 1
#error "CONTACT_FORCE_DIAG accumulates each pair once at compute time; needs ENFORCE_FLUX_SYM 1."
#endif
#if CONTACT_NO_TENSION && ENFORCE_FLUX_SYM != 1
#error "CONTACT_NO_TENSION edits each pair flux once at compute time; needs ENFORCE_FLUX_SYM 1."
#endif

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

/// GIZMO's internal-energy limiter (kicks.c:333, predict.c:226): a kick may at most
/// halve u, `if (u_new < 0.5 u_old) u = 0.5 u_old`. Keeps u > 0 without a floor, at
/// the price of exact energy conservation: the energy added is accumulated and logged
/// once per step (WARN, only on steps where the limiter fired), so a run's energy
/// error can be split into scheme error and limiter injection.
///   0 = off (baseline, bit-identical)
///   1 = on
#define ENERGY_HALF_LIMIT 1

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

// PERIODIC_BOUNDARIES note: the closure asymmetry sum is now continued over the
// periodic images (Particles::compOmega(i, ghostParticles) resumes the partial
// first moment stashed in firstMoment[i]), and compDensity(ghostParticles)
// applies the same 1/FCE lift as the non-periodic path. A particle near a box
// face therefore sees a symmetric stencil and gets FCE = 1, as it should.

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
#define EXPLICIT_VOL_SPH_DIVV 1
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

// PERIODIC_BOUNDARIES note: ghosts hold no rhoExplicit/rhoKernel of their own -
// they do not need to. applyExplicitVolumeOverride() runs after BOTH density
// passes, so `rho` already carries the advected value by the time
// updateGhostState() copies it out to the images, and the div-v that drives the
// kick comes from the ghost-aware velocity gradients.

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
