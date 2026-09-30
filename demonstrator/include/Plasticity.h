#ifndef DEMONSTRATOR_PLASTICITY_H
#define DEMONSTRATOR_PLASTICITY_H

// Pointwise pieces of the solid model that Particles::integrateStressTensor
// uses and tests/test_plasticity.cpp checks directly:
//   * the deviatoric strain-rate part of dS/dt in 2D (SOLID_2D_MODEL),
//   * the second invariant J2 and the von Mises yield return (YIELD_RETURN).
// Header-only so the unit test exercises the same code as the solver. The
// model switches are arguments defaulting to the parameter.h flags, so the
// test can check every mode from one build; the solver always passes the
// defaults, which the compiler folds.

#include <cmath>
#include "parameter.h"

namespace Plasticity {

/// Deviatoric strain-rate term 2 mu (D - c tr(D) I) of the 2D stress rate,
/// without the Jaumann rotation terms.
///   model 0: plane strain, c = 1/3 (the 3D formula with dvz/dz = 0); the
///            in-plane S is not traceless and S_zz is not stored.
///   model 1: true 2D as in GIZMO (solids/elastic_physics.c, trace_vel/NDim),
///            c = 1/2; S is traceless in 2D.
inline void devRate2D(const double mu,
                      const double dvx_dx, const double dvx_dy,
                      const double dvy_dx, const double dvy_dy,
                      double &rxx, double &ryy, double &rxy,
                      const int model = SOLID_2D_MODEL){
    if (model == 1){
        rxx =  mu*(dvx_dx - dvy_dy);
        ryy = -rxx;
    } else {
        rxx =  (4./3.)*mu*dvx_dx - (2./3.)*mu*dvy_dy;
        ryy = -(2./3.)*mu*dvx_dx + (4./3.)*mu*dvy_dy;
    }
    rxy = mu * (dvy_dx + dvx_dy);
}

/// J2 = 1/2 S:S over the stored components (2D: the 2x2 tensor).
inline double J2_2D(const double sxx, const double syy, const double sxy){
    return 0.5 * (sxx*sxx + syy*syy + 2.0*sxy*sxy);
}
inline double J2_3D(const double sxx, const double syy, const double szz,
                    const double sxy, const double sxz, const double syz){
    return 0.5 * (sxx*sxx + syy*syy + szz*szz
                  + 2.0*(sxy*sxy + sxz*sxz + syz*syz));
}

/// Factor f by which S is scaled when J2 exceeds the yield surface
/// J2 = Y^2/dim; 1 when no return is needed.
///   mode 0: radial return onto the surface, f = sqrt(Y^2/(dim J2)).
///   mode 1: GIZMO (elastic_physics.c), f = Y^2/(dim J2), no square root:
///           lands inside the surface by the fraction it overshot; the
///           difference to mode 0 vanishes as dt -> 0.
inline double yieldFactor(const double J2, const double Y,
                          const int mode = YIELD_RETURN, const int dim = DIM){
    const double Y2_over_d = Y * Y / dim;
    if (!(J2 > Y2_over_d && J2 > 0.0)) return 1.0;
    if (mode == 1) return Y2_over_d / J2;
    return std::sqrt(Y2_over_d / J2);
}

/// Plastic shear strain carried off by a return with factor f:
/// d gamma_p = (sqrt(J2) - sqrt(J2_after))/mu = (1 - f) sqrt(J2)/mu.
/// For pure shear (sqrt(J2) = |S_xy| = mu gamma_el) it is the relaxed shear strain.
inline double plasticShearIncrement(const double J2, const double f, const double mu){
    return (mu > 0.0) ? (1.0 - f) * std::sqrt(J2) / mu : 0.0;
}

} // namespace Plasticity

#endif // DEMONSTRATOR_PLASTICITY_H
