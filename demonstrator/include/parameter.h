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
#define POLY_RECONSTRUCTION 1
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
/// 2D solid model (acts only in DIM 2; 3D is unchanged).
///   0 = plane strain: dS/dt = 2 mu (D - tr(D)/3 I) + Jaumann, i.e. the 3D rate
///       with dvz/dz = 0; longitudinal modulus K + 4/3 mu. S_zz is not stored,
///       so the in-plane S is not traceless and J2 over the 2x2 tensor is not
///       the plane-strain J2 (the pre-2026-09-30 behaviour).
///   1 = true 2D, GIZMO's model (solids/elastic_physics.c, trace_vel/NDim):
///       dS/dt = 2 mu (D - tr(D)/2 I) + Jaumann, S traceless in 2D, longitudinal
///       modulus K + mu in every wave speed (Riemann solvers, HLL weights, CFL,
///       damping budget, crack speed); matches GIZMO's sqrt(c_Till^2 + mu/rho).
/// Plasticity test plan (MA-Obsidian/2026-09-30-testplan-plastisch.md): 1.
#define SOLID_2D_MODEL 0
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

/// Yield return map applied after each stress update (include/Plasticity.h).
/// The yield surface is J2 = Y^2/DIM in both.
///   0 = radial return onto the surface, S <- S sqrt(Y^2/(DIM J2)).
///   1 = GIZMO (elastic_physics.c), S <- S Y^2/(DIM J2): no square root, so the
///       state lands inside the surface by the fraction it overshot (the Benz &
///       Asphaug / miluphcuda form); the difference vanishes as dt -> 0.
#define YIELD_RETURN 0
/// Accumulate the plastic shear strain gamma_p = sum (1 - f) sqrt(J2)/mu per
/// particle and write it as /plasticStrain (diagnostic only).
#define OUTPUT_PLASTIC_STRAIN 0

/// Longitudinal (P-wave) modulus K + ELASTIC_LONG_MU_FAC mu of the solid model.
#if DIM == 2 && SOLID_2D_MODEL == 1
#define ELASTIC_LONG_MU_FAC 1.0
#else
#define ELASTIC_LONG_MU_FAC (4./3.)
#endif

#if SOLID_2D_MODEL != 0 && SOLID_2D_MODEL != 1
#error "SOLID_2D_MODEL must be 0 (plane strain) or 1 (true 2D, GIZMO)."
#endif
#if YIELD_RETURN != 0 && YIELD_RETURN != 1
#error "YIELD_RETURN must be 0 (radial return) or 1 (GIZMO)."
#endif
#if OUTPUT_PLASTIC_STRAIN && !PLASTICITY_ANY
#error "OUTPUT_PLASTIC_STRAIN needs a plasticity model (PLASTICITY or *_PLASTICITY)."
#endif

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
#define TENSILE_CORRECTION_LOCAL 1

// Upper clamp on the correction factor (TENSILE_CORRECTION_LOCAL only).
// Without it the (1 - s) profile keeps growing as a pair collapses and turns
// into an impulsive kick on an already-collapsed pair -- measured as negative
// internal energy at vmax = 0.60. Physically: never convert more than this
// fraction of the tensile traction into repulsion.
#define TENSILE_CORRECTION_FMAX 1.0
// Separate clamp for the factor inside the Riemann solver (TENSILE_CORRECTION_LOCAL only),
// i.e. for TENSILE_CORRECTION_1's anti-clumping repulsion: the solve removes only
// (1 - f) of the dummy pressure, so the net repulsion on a tensile face is f*|p_min|,
// and f > 1 pushes a collapsing pair apart harder than the tension pulls it together.
// FMAX above keeps clamping every other consumer (the stress-flux damping (1 - f_i),
// where f > 1 would flip the sign of the tensile stress, and the approach damping).
// Under GIZMO_ELASTIC_FLUX 1 the deviator is zeroed before the solve, so this cap
// touches only the repulsion.
// Default 3.0 since 2026-09-30 (agR30 / g2R30): the only variant that removed the tensile
// clumping of the fine rings (d0.035) and it halves the block AM error at d0.00625
// (testcases/hires_20260929, tgate2_20260930). FMAX_HLLC = FMAX reproduces the old default.
#define TENSILE_CORRECTION_FMAX_HLLC 3.0

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
// precision at any strength. Semi-discretely it only REMOVES kinetic energy;
// the explicit kick does so only below lambda dt ~ 1 (see the limiter below:
// beyond it the overshoot cools single particles, the "cold streaks").
#define TENSILE_CORRECTION_DAMP 18.0

// Stability limit for the approach damping above (2D GIZMO stress flux only;
// the 3D branch has no DAMP boost).
//
// Why: the boosted radial dissipation drives the pair's relative normal
// velocity to zero at the rate lambda_ij = w_r (1 + DAMP f) / mu_ij,
// mu_ij = m_i m_j / (m_i + m_j). Unboosted, lambda dt is O(CFL) (~0.09 per
// nearest-neighbour pair); with DAMP 18 and f = FMAX = 1 it is ~1.7. The kick
// is explicit and, under TIME_KDK, driven by the PREDICTED velocity, whose
// linear stability limit is lambda dt = 1 (not 2). Beyond it the relative
// velocity overshoots, and because u = E/m - v^2/2 the overshoot shows up as
// cooling of single particles: cold streaks along the lattice rows of the
// offset blocks (MA-Obsidian/2026-09-29-damp-stability-limiter-idea.md).
// Both limited forms scale only the coefficient, so the force stays central
// and pairwise antisymmetric (L and p exact); the unboosted w_r is never
// reduced, so a pair that would not overshoot is untouched.
// dt is the velocity change one flux pass produces: under TIME_KDK the closing
// half of the previous step plus the opening half of this one,
// (dt_prev + dt)/2; otherwise dt.
// Tested 2026-09-29 (testcases/damp_limit_20260929, MA-Obsidian/2026-09-29-damp-limit):
// every mode removes the streaks at C <= 0.7 (blocks to d0.00625: 164 -> 0), but
// the T.1 rings d0.05 then clump at the contact tips -- p1d18f10's ring margin
// came partly from the overshoot. No C passes both. Default 1 / C 0.5 since
// 2026-09-29 evening (user priority: L_z convergence > no negE > clumping; test
// plan reference, study testcases/hires_20260929).
// NOTE lambda ~ 1/h and dt ~ h, so lambda dt is resolution-independent at fixed
// CFL: an active limiter means DAMP acts as "at most the fraction C of the
// relative normal velocity per step" -- a numerical, not a physical, strength.
//
//   0 = off (bit-identical to 5865de1)
//   1 = cap:  w = min(w_r (1 + DAMP f), max(w_r, C mu_ij / dt))
//   2 = exponential: w = max(w_r, mu_ij (1 - exp(-lambda_ij dt)) / dt), the
//       exact decay of the pair's relative normal velocity over dt;
//       = w_r (1 + DAMP f) for small lambda dt, never overshoots
//   3 = per-particle budget: w = max(w_r, w_r (1 + DAMP f) min(s_i, s_j)),
//       s_i = min(1, C / (dt Lambda_i)), Lambda_i = sum of lambda_ij over i's
//       boosted pairs (Gershgorin-type). Modes 1/2 bound each pair alone, but a
//       contact-face particle is hit by several closing partners at once
//       (measured with OUTPUT_DAMP_DIAG, blocks d0.05: sum lambda dt 1.46 at a
//       per-pair max of 0.58).
#define TENSILE_CORRECTION_DAMP_LIMIT 1
// C of modes 1 and 3 (fraction of the relative normal velocity removed per step at most)
#define TENSILE_CORRECTION_DAMP_LIMIT_C 0.5

// Traction gate for EVERY tensile-correction channel (user decision 2026-09-30,
// MA-Obsidian/2026-09-30.md): a correction may act only where the material is in
// effective tension, i.e. where the traction is negative on SOME plane: the largest
// principal value sigma_1 of the total stress sigma = -P I + S is > 0.
// (Version 1, commit 1c05cca, used the traction along the pair / face normal instead.
// Wrong: in the rings' outer surface layer under hoop tension the particles pair up
// RADIALLY, where the traction is ~0, so the gate switched the correction off exactly
// where the material was in tension -- found by the user 2026-09-30.)
// Clumping under compression (e.g. the contact edge of the colliding rings) is a
// different defect and must not be hidden by a tensile correction; with the gate
// every correction is visibly confined to tension, so a later tear can be
// attributed to physics (damage, yield) or to the correction.
// f itself stays geometric (pair collapse); the gate decides where it may act:
//   TENSILE_CORRECTION_1 (Riemann repulsion): unchanged. Under GIZMO_ELASTIC_FLUX it acts
//       only at P < 0 (hydrostatic tension), so it is already confined to tension.
//   stress-flux damping (1 - f_i): only principal components whose TOTAL principal
//       stress lambda_k - P is tensile (lambda_k > max(0, P)); still one tensor per
//       particle, so ELASTIC_FACE_MODE 3 stays torque-free.
//   approach damping (DAMP): only for pairs whose mean sigma_1 is > 0; symmetric in
//       (i, j), so the boosted force stays central and antisymmetric.
//   0 = off (the ungated scheme up to 2026-09-30), 1 = on (default since 2026-09-30).
// Builds without TENSILE_CORRECTION or without the GIZMO elastic flux (e.g. GEF 0
// comparisons, tc0) must set it to 0, like TENSILE_CORRECTION_DAMP_LIMIT.
#define TENSILE_CORRECTION_TRACTION_GATE 1
// Per-particle approach-damping diagnostic in the h5 dump (TENSILE_CORRECTION_LOCAL,
// 2D): dampLdt = sum over the particle's boosted pairs of lambda_ij dt (UNLIMITED,
// Gershgorin-type load), dampLdtMax = max over them, dampNb = number of boosted
// pairs, dampNcap = number of those the limiter reduced (this flux pass), and
// dampLdtRun = running max of dampLdt over the run. Numerics unchanged.
#define OUTPUT_DAMP_DIAG 0
// Per-particle internal-energy budget in the h5 dump (2026-09-30, to find what drives u < 0).
// Numerics unchanged (one extra Riemann solve per pair without the repulsion, used only for
// the split). With finite masses one kick changes u exactly by
//     du = -dt/m * sum_c ( E_c - v_mid . p_c ),   v_mid = (v_old + v_new)/2,
// where (E_c, p_c) are the energy and momentum fluxes of channel c; the running integral of
// each channel's du over the real kicks (the TIME_KDK predictor is not counted) is dumped:
//   ebH   hydro Riemann flux without the tensile repulsion (pressure work)
//   ebR   TENSILE_CORRECTION_1 repulsion (full Riemann flux minus ebH)
//   ebS   GIZMO stress flux, deviator, undamped
//   ebSD  its tensile damping (1 - f_i): damped minus undamped deviator flux
//   ebV   stress-flux dissipation (longitudinal wt_r without the DAMP boost, transverse wt_t)
//   ebD   approach-damping boost (DAMP): the extra longitudinal dissipation
// Check: u(t) - u(0) = ebH + ebR + ebS + ebSD + ebV + ebD to round-off.
#define OUTPUT_ENERGY_BUDGET 0

// Two Parameters for the tensile correction facor, as in Monaghan, 2000
// (TENSILE_CORRECTION_LOCAL reuses both: PREFAC scales f, POWER sharpens the
// (1 - s) profile. Because f == 0 on a good lattice, PREFAC may be raised well
// above Monaghan's 0.2 without paying a permanent tax.)
#define TENSILE_CORRECTION_PREFAC 4.0
#define TENSILE_CORRECTION_POWER 1.

// Per-PARTICLE instead of per-PAIR damping of the tensile principal stresses in
// the GIZMO elastic stress flux (GIZMO_ELASTIC_FLUX only).
//
// Why: the deviatoric traction of a pair, -0.5 * S . FVec, is not central for
// any single pair; the total torque vanishes only because each particle applies
// ONE symmetric S_i on all of its pairs (ELASTIC_FACE_MODE 3: exactly, through
// sum_j a_ij (x) x_ij = -V_i I). The per-pair (1 - f_ij) on the positive
// eigenvalues makes the effective stress of particle i differ from pair to pair
// and breaks that cancellation (disc plateau ~2e-5, 2026-09-26; suspected main
// source of the A.2 block L_z error, MA-Obsidian/2026-09-28-kontakt-ablation.md).
// Here each side uses its OWN particle's f_i, so S~_i is again one symmetric
// tensor per particle and FACE 3 stays torque-free at any strength. The stored S
// is not modified (reversible, flux-only -- a geometric softening of tension,
// not plasticity). f_i is built from the same clamped pair values
// PREFAC*pow(f_ab, POWER) (<= FMAX for _LOCAL 1), over normal and periodic pairs.
// The approach damping (DAMP) and TENSILE_CORRECTION_1's HLLC repulsion keep the
// pair value.
//
//   0 = per pair (GIZMO default, bit-identical to pre-2026-09-28 code)
//   1 = per particle, f_i = max_j f_ij
//   2 = per particle, f_i = mean_j f_ij over all neighbours (zeros included)
//   3 = ablation: f_i = 0, i.e. NO tensile damping in the stress flux at all
//       (DAMP and the HLLC repulsion still act); separates what the damping
//       itself does to clumping from what removing the per-pair form does
// Default since 2026-09-29 (test plan §1): 1 with DAMP 18, FMAX 1.0 (tuning MA-Obsidian/2026-09-29-tc-per-particle-tuning).
#define TENSILE_CORRECTION_PARTICLE 1
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
#define RECON_FACE_DV 2

// Take the RECON_FACE_DV jump from the SPATIAL reconstruction only, i.e. the
// limited face states before the MUSCL half-step predictor adds dt/2 * dW/dt.
//
// Why: the predictor advances each side with its own Taylor step,
//     v_f,i(t + dt/2) = v_i + (x_f - x_i).grad v_i + dt/2 [a_i - (v_i - v_F).grad v_i],
// which drops the mixed term dt/2 (x_f - x_i) . d_t grad v. For a linear field
// v = G(t) x the two sides then differ by dt/2 * dG/dt * x_ij instead of by zero,
// so a body whose rotation rate or strain rate changes (a spinning disc that
// breathes) feels a spurious transverse dissipation of order dt*h, and its
// antisymmetric part (the spin acceleration) is a torque. Measured 2026-09-26
// (testcases/am_timecentering/20260926): rotating disc d0.05 omega 0.2, the
// transverse channel carries 86 % of the residual L_z torque and halves with dt
// (3.2e-4 -> 1.7e-4 of L_z at t = 10, CFL 0.2 -> 0.1).
// The jump only has to vanish for linear fields at every instant; it does not
// have to be time-centred, since it is itself O(h^2) on smooth fields.
//
// FRAME BUG of mode 0 (found 2026-09-26 evening, testcases/am_timecentering/
// 20260926_kdkdrift): mode 0 reads the jump as WijR - WijL inside the stress
// flux, i.e. AFTER the Riemann constructor has rotated the face velocities in
// place into the face frame (Riemann.cpp, "Rotate velocities (xSplitElasticHLLC
// expects vx = normal velocity)"), but projects it on the LAB-frame r_ij. The
// jump therefore points the wrong way on every face whose normal is not +x;
// on faces with normal ~ -x it is sign-flipped (anti-dissipation). Modes 1 and
// 2 store a lab-frame copy during the reconstruction, before the solve.
// Evidence: TIME_KDK (whose jump is spatial by construction) converges to a
// different solution than the baseline/TIME_DKD as dt -> 0 (max|dv| 2.6e-3 at
// t = 0.02, dt = 1e-4, disc d0.05) with mode 0, and to the same one (2.5e-5)
// with mode 1; with RECON_FACE_DV 2 mode 0 disrupts the disc at t ~ 1.5 for
// any CFL. Part of what was attributed to the predictor above is this bug.
//
//   0 = jump read from the rotated face states (predictor included; FRAME BUG)
//   1 = spatial-only jump, lab frame (the predicted states still feed the HLLC solve)
//   2 = predicted jump, lab frame (mode 0 without the frame bug; diagnostic)
#define RECON_FACE_DV_SPATIAL 1

// Minmod-limit the RADIAL part of the RECON_FACE_DV jump against the raw
// difference v_i - v_j: same sign -> the smaller of the two, opposite sign -> 0.
// (The transverse part is not limited: against the one-signed raw difference of a
// rotating body a minmod brakes the spin, -2.5e-3 of L_z per rotation measured.)
//
// Why: the HLL term of a pair changes the kinetic energy at the rate
// -wt (dv_F . dv), which is a loss only while the reconstructed jump dv_F points
// the same way as the raw one. Where the two sides' extrapolations overshoot each
// other the term pumps kinetic energy out of u (the work term books it there), u
// goes negative and the body heats its own motion. Measured 2026-09-26: TIME_KDK +
// RECON_FACE_DV 2 disrupts the rotating disc at t ~ 10 on d0.025 and d0.0177
// (u -> -96, 4.3e6 negative-u lines). Rigid rotation keeps dv_F = 0; on smooth
// fields the limited jump stays O(h^2) (or 0), so the order is not affected.
//
//   0 = off, 1 = minmod against the raw difference (2D only)
#define RECON_FACE_DV_MINMOD 0

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
//
//   3 = split MFM half faces: each side's stress on its OWN half face,
//           F_ij = S_i a_ij - S_j a_ji,   a_ij = V_i psiTilde_j(x_i),
//       where A_ij = a_ij - a_ji. Mode 2 takes 0.5 (S_i + S_j) A_ij, and its
//       torque is exact only for a UNIFORM stress: expanding the pair sum, the
//       S_i a_ij part cancels through sum_j a_ij (x) x_ij = -V_i I (linear
//       consistency of psiTilde, B = E^-1), but the cross term
//       0.5 sum x_ij x S_j a_ij leaves a residual wherever S varies -- measured
//       2026-09-26 as the whole L_z residual of the rotating disc once the time
//       errors are gone (surface band, not convergent in d). The split form is
//       still pairwise antisymmetric (momentum exact) and its total torque is
//       sum_i eps : (S_i sum_j a_ij (x) x_ij) = -sum_i V_i eps : S_i = 0 for every
//       symmetric S_i. The tensile correction's per-pair (1 - f) and the
//       dissipation terms are the only non-central parts left.
#define ELASTIC_FACE_MODE 3
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
// Default since 2026-09-29 (test plan §1). On the plan's Tillotson (P independent of u) it changes only u;
// x, v, rho stay bitwise identical (blocks d0.05).
#define MFM_ENERGY_STARFIX 1

/// "Option B" (2026-09-24 TODO 3, operator_splitting.tex Sec. 4.4): evaluate the
/// MFM flux in the frame of the contact wave, S*' = 0, instead of zeroing the
/// advective pieces by hand in the solver frame.
///
/// In that frame the FULL HLLC star flux (MFV form) is analytically
///     (mass, energy, v_n, v_t) = (0, -W_t S_xy, t*, -S_xy),
///     t* = t_K + rho_K (S_K - v_K)(S* - v_K)   (same for K = L, R by the S* formula),
/// i.e. zero mass flux by construction, and transformed back to the solver frame
/// (surface moving with S*) it is  (0, S* t* - W_t S_xy, t*, -S_xy)  -- GIZMO's
/// MFM flux Fluxes.v = P_M, Fluxes.p = P_M S_M (reimann.h:411-412). The caller's
/// de-boost to the lab is unchanged. The tangential side (W_t, S_xy) is picked by
/// the sign of S* exactly as in the baseline.
///
/// The hand-zeroed baseline branch agrees only at S* = 0: its normal momentum is
/// t* + rho*_K S* (S* - v_K) and its energy (even with MFM_ENERGY_STARFIX 1) keeps
/// t_K (S* - v_K) S* / (S_K - S*) -- both second order in the velocity jump, so
/// invisible in linear waves, O(S*/c) at contacts. With this flag on,
/// MFM_ENERGY_STARFIX has no effect (the branch it patches is not reached).
///
///   0 = baseline (hand-zeroed MFM branch in the solver frame)
///   1 = exact contact-frame flux
#define MFM_CONTACT_FRAME 0

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
/// DEFAULT ON since 2026-09-25 (MA-Obsidian/2026-09-25-erste-ordnung-drift.md):
/// the endpoint drift is the demonstrator's O(dt) soundwave error.  Faces use
/// V = m/rho_explicit with position-based psi-tilde weights, so the background
/// pressure P0*sum_j A_ij couples the kernel density (the positions) into the
/// momentum equation; an off-centred drift then damps the wave at the rate
/// (P0/rho0) k^2 dt / 4.  GIZMO drifts KDK with v^{n+1/2}.  With 1 the ideal-gas
/// soundwave converges at p ~ 2 (dt ∝ h).  Absent for P0 = 0 (solid ICs).
///
/// CAVEAT (2026-09-26, MA-Obsidian/2026-09-26-am-zeitzentrierung.md): this
/// centres the DRIFT but not the KICK, which still uses forces built at x^n.
/// Then L^{n+1} - L^n = dt sum x^n x F + dt^2/2 sum v^{n+1/2} x F exactly, and the
/// second term is +omega^2 dt/2 per unit time on a rotating body although every
/// pair force is central (the endpoint drift conserves L exactly for central
/// forces). The Strang half steps of S and rhoExplicit after the kick also read
/// grad v^{n+1} on x^n, so rhoExplicit/rhoKernel drifts by +omega^2 dt per unit
/// time. For rotation use TIME_DKD (or TIME_KDK); see there.
///
///   0 = single-endpoint drift (pre-2026-09-25 behaviour, first order)
///   1 = trapezoidal drift, x += dt (v^n + v^{n+1})/2
#define DRIFT_KDK 0

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
#define RECON_TIME_PRED 0

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

/// GIZMO's kick-drift-kick: evaluate everything at the sync time t^n and
/// apply the result as two half kicks around it.
///
/// The baseline evaluates at t^n and kicks over [t^n, t^n + dt]; DRIFT_KDK
/// centres the drift but not the kick, and the Strang half steps of S and
/// rhoExplicit that bracket the flux pass read grad v^{n+1} on the geometry
/// x^n.  Measured 2026-09-26 (testcases/am_timecentering/20260926) on the
/// rotating disc:
///   * L_z grows by dt^2/2 sum_i v_i x F_i per step = +omega^2 dt/2 per unit
///     time (1.37 % measured / 1.42 % predicted, d0.025 omega 0.2, t = 100),
///     first order in dt although every pair force is central;
///   * rhoExplicit/rhoKernel drifts by +omega^2 dt per unit time and the disc
///     inflates (1.68 % / 1.72 % over one rotation, d0.05).
/// Both vanish when every rate is taken at one geometry and one time level.
///
/// Step n under TIME_KDK (one neighbour/flux pass per step, like the baseline):
///   1. positions x^n; the working arrays hold the drift-PREDICTED sync state
///      (v, u, S, rhoExplicit advanced from the half-step store by dt/2 with
///      the previous rates, GIZMO's VelPred).  Density, gradients, faces and
///      fluxes are built from it; the reconstruction has NO half-step predictor.
///   2. closing half kick of step n-1 from the half-step store with the NEW
///      rates -> synchronised state at t^n; the snapshot is written here.
///   3. opening half kick with the same rates -> half-step store (n + 1/2).
///   4. drift x^{n+1} = x^n + dt v^{n+1/2}; predict the working arrays to t^{n+1}.
/// Each half kick sits at one fixed geometry, so a central pair force changes
/// L_z by exactly zero, and S and rhoExplicit see grad v at the time and place
/// it was measured.  The dump moves from the top of the step to point 2 (the
/// state at the top is a prediction, not the solution).
///
/// Needs DRIFT_KDK 0, TIME_DKD 0 (both are other answers to the same question)
/// and RECON_TIME_PRED 0: the flux is centred on t^n by construction, and a
/// half-step predictor would push it to t^n + dt/2.
///
///   0 = baseline (evaluate at t^n, kick forward a full dt)
///   1 = kick-drift-kick around the sync time
// Default since 2026-09-29 (test plan §1): the TIME_KDK candidate = TIME_KDK 1, TIME_KDK_FIELD_DRIFT 2,
// DRIFT_KDK 0, RECON_TIME_PRED 0, RECON_FACE_DV 2, RECON_FACE_DV_SPATIAL 1, ELASTIC_FACE_MODE 3.
#define TIME_KDK 1

/// Under TIME_KDK, advance S and rhoExplicit in the DRIFT, like the positions,
/// instead of kicking them like v.
///
/// TIME_KDK 1 kicks S and rhoExplicit in half steps with rates taken from the
/// drift-PREDICTED state (v_pred = v^{n+1/2} + dt/2 a^{n-1}). Written in the
/// half-step variables that is Adams-Bashforth 2 for the v-S pair:
///   v^{n+1/2} = v^{n-1/2} + dt F((3 S^{n-1/2} - S^{n-3/2})/2),
///   S^{n+1/2} = S^{n-1/2} + dt G((3 v^{n-1/2} - v^{n-3/2})/2),
/// and AB2 amplifies every undamped oscillation by ~1 + (w dt)^4/4 per step
/// (toy oscillator: 1.00043 per step at w dt = 0.2, 1.027 at 0.5). Resolved waves
/// (w dt ~ 0.03) never notice; grid-scale modes (w dt ~ 0.5-1) are held only
/// by the HLL dissipation, which is weak at a free surface and O(h^2) with
/// RECON_FACE_DV 2. Measured 2026-09-26 on the concentric-ring disc d0.025:
/// the outer ring roughens (r scatter 3e-5 -> 6e-4 by t = 10, u 0.78..1.08),
/// the HLLC channel torques the disc up, and without the radial minmod
/// RECON_FACE_DV 2 disrupts it.
///
/// Here S and rhoExplicit are position-like: after the opening half kick the
/// velocity gradient of v^{n+1/2} is rebuilt on the flux-pass geometry x^n and
/// both are advanced by the full dt with it, before the particles drift with the
/// same v^{n+1/2}. v and u keep the KDK half kicks. The v-S pair is then a
/// leapfrog (stable for w dt < 2), S^n and rhoExplicit^n at the flux pass are
/// the solution, not a prediction, and the half-step store holds v and u only.
///
/// The gradient is taken on x^n but belongs to the drift x^n -> x^{n+1}.
/// Mode 2 maps it to the midpoint geometry: with G = grad_{x^n} v^{n+1/2} and
/// x^{n+1/2} = x^n + dt/2 v^{n+1/2}, grad_{x^{n+1/2}} v^{n+1/2} = G (I + dt/2 G)^-1,
/// exact for any linear velocity field. For rigid rotation G carries the spurious
/// divergence -omega^2 dt (the half kick adds the centripetal -omega^2 x/2 dt);
/// mode 1 keeps it (rhoExplicit drifts by +omega^2 dt per unit time, as in
/// the baseline), mode 2 cancels it.
///
///   0 = S and rhoExplicit kicked with the predicted rates (TIME_KDK 1 as is)
///   1 = S and rhoExplicit drifted with grad_{x^n} v^{n+1/2}
///   2 = as 1, gradient mapped to the midpoint geometry
#define TIME_KDK_FIELD_DRIFT 2

// ---------------------------------------------------------------------------
// Third-order study 2026-10-01 (branch thirdorder, MA-Obsidian/2026-10-01-dritte-ordnung.md).
// Hopkins (2015) sec. 2.8: higher order "by increasing the order of the gradient" -- plus
// higher-order reconstruction, time integration and quadrature/effective faces. Measured on
// the ideal-gas soundwave (truncation prototype testcases/thirdorder_20261001/proto): raising
// GRADIENT_ORDER alone leaves every equation second order; the momentum equation is held at
// O(h^2) by the face closure sum_j A_ij (the background pressure P0 sum_j A_ij), which is
// O(h^2) with the first-order psiTilde even for exact volumes, and O(h^4) with the
// cubic-exact weights. All three flags default to 0 (bit-identical).
// ---------------------------------------------------------------------------

/// Build the effective faces A_ij = V_i psiTilde_j(x_i) - V_j psiTilde_i(x_j) from the
/// higher-order least-squares gradient weights psijTildeGrad (GRADIENT_ORDER >= 2) instead of
/// the first-order psijTilde. A_ij stays antisymmetric for any weights, so conservation is
/// unchanged; any consistent gradient weights also keep sum_j psiTilde_j (x_j - x_i)^T = I.
/// With the cubic-exact weights (GRADIENT_ORDER 3) the face closure sum_j A_ij of a smoothly
/// deformed lattice drops from O(h^2) to O(h^4), and on a point-symmetric stencil the central
/// flux operator becomes fourth order even with LINEAR face reconstruction (the cubic weights
/// annihilate the O(D^3) midpoint term). Particles on the GRAD_MAT_COND_MAX fallback carry the
/// first-order weights there, so the faces fall back with them.
///   0 = first-order weights (Hopkins 2015 / GIZMO)
///   1 = higher-order weights psijTildeGrad
#define FACE_GRAD_WEIGHTS 0

/// With FACE_GRAD_WEIGHTS, 78 % of the faces of the square lattice point against x_j - x_i (the
/// negative outer coefficients of the cubic-exact weights, like a fourth-order difference), and
/// the Riemann solver upwinds along A_ij. 1 = hand the solver the two sides in geometric order
/// for those faces (A_ij . (x_j - x_i) < 0). Measured 2026-10-01: does NOT stabilise (the growth
/// comes from the background-pressure closure, see MOMENTUM_PREF; linearised growth 3.6 -> 5.0)
/// and, being a switch on the sign of A.dx, it injects a non-smooth error where faces nearly
/// perpendicular to dx flip with the wave phase. Kept for the record; leave at 0.
#define FACE_GRAD_WEIGHTS_SWAP 0

/// Seed rhoExplicit at t = 0 from an IC dataset /rho instead of the kernel-sum density (and
/// skip EXPLICIT_VOL_SEED_SKIP). The kernel sum is a smoothed density, its error O(h^2); a
/// third-order run needs the density to third order from the start, which the IC generator
/// knows exactly (point values rho(x_i)). The IC must carry /rho.
///   0 = seed from the kernel density (default)
///   1 = seed from /rho
#define EXPLICIT_VOL_SEED_IC 0

/// Remove a global reference pressure from the MOMENTUM flux: F_mom,ij = (P*_ij - P_ref) A_ij,
/// P_ref = volume-weighted mean pressure of the pass. A_ij = -A_ji, so sum_i sum_j A_ij P_ref = 0
/// and momentum stays exactly conserved; the energy flux keeps the full P* (compression work of
/// the background pressure is physical), so total energy stays conserved too.
/// Why: a uniform pressure exerts no force in the continuum, but P0 sum_j A_ij does in the scheme.
/// With the first-order faces that closure is O(h^2) (it limits the ideal-gas soundwave to second
/// order) and it is also the weak restoring force of the particle lattice; with FACE_GRAD_WEIGHTS it
/// is O(h^4) for smooth deformations but ANTI-restoring at mid wavenumbers: linearised MFM at rest
/// (proto/mfm_stab.py, n = 16) has real eigenvalues +3.6 (growth ~ 1/dx) with P0 = 0.6, none
/// without it. Removing P_ref leaves the lattice modes neutral (max Re lambda ~ 1e-9).
/// Only the constant part of the pressure is removed, so this cures the closure error exactly for
/// small perturbations of a uniform state, not for O(1) pressure variations.
///   0 = off
///   1 = P_ref = volume-weighted mean pressure
#define MOMENTUM_PREF 0

/// Third-order strong-stability-preserving Runge-Kutta (Shu & Osher 1988) instead of the
/// kick-drift schemes above:
///   U1 = U0 + dt L(U0),  U2 = 3/4 U0 + 1/4 (U1 + dt L(U1)),  U3 = 1/3 U0 + 2/3 (U2 + dt L(U2)).
/// Every stage rebuilds the full state at its own geometry (neighbour search, h, density,
/// gradients, faces, fluxes); the stage update is the Euler step x += dt v, (m v, E) -= dt F,
/// ln rhoExplicit -= dt div v, S += dt dS/dt, and the stages are combined linearly in x
/// (minimum image), v, the specific total energy (so E stays conserved), ln rhoExplicit and S.
/// No half-step predictor (RECON_TIME_PRED 0). Three flux passes per step. Fixed timestep only.
///   0 = off
///   1 = SSP-RK3
#define TIME_RK3 0

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
#if TIME_KDK && (DRIFT_KDK || TIME_DKD)
#error "TIME_KDK is its own time-centring scheme: set DRIFT_KDK 0 and TIME_DKD 0."
#endif
#if TIME_KDK && RECON_TIME_PRED
#error "TIME_KDK evaluates the fluxes at the sync time; the half-step predictor would move them to t + dt/2: set RECON_TIME_PRED 0."
#endif
#if FACE_GRAD_WEIGHTS && (!USE_HLLC || FRAGMENTATION || OUTPUT_ENERGY_BUDGET)
#error "FACE_GRAD_WEIGHTS: the backward-face swap is wired into the HLLC solves only (no damage, no energy budget)."
#endif
#if FACE_GRAD_WEIGHTS && !SECOND_ORDER_GRADIENTS
#error "FACE_GRAD_WEIGHTS builds the faces from the higher-order gradient weights: set GRADIENT_ORDER >= 2."
#endif
#if TIME_RK3 && (TIME_KDK || DRIFT_KDK || TIME_DKD || RECON_TIME_PRED || FRAME_HALFSTEP)
#error "TIME_RK3 is its own time integrator: set TIME_KDK, DRIFT_KDK, TIME_DKD, RECON_TIME_PRED and FRAME_HALFSTEP to 0."
#endif
#if TIME_RK3 && (ADAPTIVE_TIMESTEP || !ELASTIC)
#error "TIME_RK3 is implemented for a fixed timestep on the ELASTIC updateState/moveParticles path."
#endif
#if TIME_KDK_FIELD_DRIFT && !TIME_KDK
#error "TIME_KDK_FIELD_DRIFT changes how TIME_KDK advances S and rhoExplicit: enable TIME_KDK."
#endif
#if TIME_KDK_FIELD_DRIFT > 2
#error "TIME_KDK_FIELD_DRIFT must be 0 (kicked), 1 (drifted, grad on x^n) or 2 (drifted, midpoint grad)."
#endif
#if TENSILE_CORRECTION_PARTICLE && !(TENSILE_CORRECTION && GIZMO_ELASTIC_FLUX)
#error "TENSILE_CORRECTION_PARTICLE replaces the per-pair damping of the GIZMO stress flux: needs TENSILE_CORRECTION 1 and GIZMO_ELASTIC_FLUX 1."
#endif
#if TENSILE_CORRECTION_PARTICLE > 3
#error "TENSILE_CORRECTION_PARTICLE must be 0 (per pair), 1 (per particle, max), 2 (per particle, mean) or 3 (none)."
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
#if (TENSILE_CORRECTION_DAMP_LIMIT || OUTPUT_DAMP_DIAG) && !(TENSILE_CORRECTION && TENSILE_CORRECTION_LOCAL && GIZMO_ELASTIC_FLUX && DIM == 2)
#error "TENSILE_CORRECTION_DAMP_LIMIT / OUTPUT_DAMP_DIAG act on the approach damping, which exists only in the 2D GIZMO stress flux with TENSILE_CORRECTION_LOCAL."
#endif
#if TENSILE_CORRECTION_TRACTION_GATE && !(TENSILE_CORRECTION && GIZMO_ELASTIC_FLUX && USE_HLLC && !ELASTIC_HLLC_EP)
#error "TENSILE_CORRECTION_TRACTION_GATE is implemented for the GIZMO elastic flux on the x-split HLLC path."
#endif
#if TENSILE_CORRECTION_TRACTION_GATE && TENSILE_CORRECTION_DAMP_LIMIT == 3
#error "TENSILE_CORRECTION_TRACTION_GATE: the per-particle damping budget (DAMP_LIMIT 3) is not gated yet."
#endif
#if TENSILE_CORRECTION_DAMP_LIMIT > 3
#error "TENSILE_CORRECTION_DAMP_LIMIT must be 0 (off), 1 (cap), 2 (exponential) or 3 (per-particle budget)."
#endif
#if RECON_FACE_DV_MINMOD && !(RECON_FACE_DV && DIM == 2)
#error "RECON_FACE_DV_MINMOD limits the RECON_FACE_DV jump and is implemented for DIM == 2."
#endif
#if RECON_FACE_DV_SPATIAL > 2
#error "RECON_FACE_DV_SPATIAL must be 0 (rotated, legacy), 1 (spatial, lab frame) or 2 (predicted, lab frame)."
#endif
#if RECON_FACE_DV_SPATIAL && !RECON_FACE_DV
#error "RECON_FACE_DV_SPATIAL selects the source of the RECON_FACE_DV jump; enable RECON_FACE_DV."
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
// ELASTIC_FACE_MODE consistency guards
#if ELASTIC_FACE_MODE && !GIZMO_ELASTIC_FLUX
#error "ELASTIC_FACE_MODE only affects the GIZMO elastic stress flux; enable GIZMO_ELASTIC_FLUX."
#endif
#if ELASTIC_FACE_MODE && DIM != 2
#error "ELASTIC_FACE_MODE is implemented for DIM == 2 only."
#endif
#if ELASTIC_FACE_MODE > 3
#error "ELASTIC_FACE_MODE must be 0 (raw), 1 (Bonet-Lok), 2 (MFM face) or 3 (split MFM half faces)."
#endif
#if ELASTIC_FACE_MODE >= 2 && RECON_FACE_TRACTION
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

// floor for internal energy u, set to negative value to disable; may be non-integer (e.g. 0.01).
// Not conservative: the energy it adds is logged per kick ("Energy floor: clamped ..."), and u <= 0
// still logs the "Negative internal energy" ERROR before the clamp.
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

// TIME_KDK guards that read flags defined late in this file (MOVE_PARTICLES, RUNSPH)
#if TIME_KDK && !(ELASTIC && MOVE_PARTICLES)
#error "TIME_KDK is implemented on the ELASTIC updateState/moveParticles path only."
#endif
#if TIME_KDK && (FRAGMENTATION || AM_TORQUE_TRACK || RUNSPH)
#error "TIME_KDK: damage, the AM_METHOD 4-6 repairs and SPH are not wired into the half-kick store."
#endif
#if TIME_RK3 && (!MOVE_PARTICLES || FRAGMENTATION || AM_TORQUE_TRACK || RUNSPH || TORQUE_DIAG || AM_SELFTEST \
                 || TENSILE_CORRECTION_DAMP_LIMIT || OUTPUT_DAMP_DIAG || OUTPUT_ENERGY_BUDGET)
#error "TIME_RK3: moving particles only; damage, AM repairs/diagnostics, SPH, the damping limiter and the energy budget are not wired into the stages."
#endif
#if EXPLICIT_VOL_SEED_IC && (!EXPLICIT_VOL_INTEGRATION || EXPLICIT_VOL_SEED_RHO0)
#error "EXPLICIT_VOL_SEED_IC seeds rhoExplicit from the IC: needs EXPLICIT_VOL_INTEGRATION 1 and EXPLICIT_VOL_SEED_RHO0 0."
#endif

#if OUTPUT_ENERGY_BUDGET && !(DIM == 2 && ELASTIC && USE_HLLC && GIZMO_ELASTIC_FLUX && TENSILE_CORRECTION && ENFORCE_FLUX_SYM && MESHLESS_FINITE_MASS && !PERIODIC_BOUNDARIES && !FRAGMENTATION && !ADIABATIC_FACE_CORRECTION && AM_METHOD == 0)
#error "OUTPUT_ENERGY_BUDGET: 2D GIZMO elastic flux with TENSILE_CORRECTION, symmetric fluxes, finite mass, no periodic boundaries / fragmentation / AFC / AM corrections."
#endif
#endif //DEMONSTRATOR_PARAMETER_H
