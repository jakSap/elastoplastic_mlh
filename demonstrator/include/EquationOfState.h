//
// Created by Jakob Sappler on 18.09.25
//

#ifndef MESHLESSHYDRO_EQUATIONOFSTATE_H
#define MESHLESSHYDRO_EQUATIONOFSTATE_H

#include <cmath>
#include <cassert>
#include <vector>

#include "parameter.h"
#include "Logger.h"

#if EOS == 1
/// Per-material parameters for the Murnaghan EOS path. The shear modulus
/// `mu` is bundled here even though it is a constitutive (not strictly EOS)
/// quantity, so a single matId lookup yields all per-material parameters.
/// Future plasticity work can extend this struct with Y0.
struct MurnaghanMaterial {
    double K0;       ///< bulk-modulus reference
    double n;        ///< Murnaghan exponent
    double rho0;     ///< reference density
    double mu;       ///< shear modulus (constitutive)
};
#elif EOS == 2
/// Per-material parameters for the Tillotson EOS path. The shear modulus
/// `mu` is bundled here for the same reason as in MurnaghanMaterial.
/// References: Tillotson (1962); Melosh, "Impact Cratering", Appendix B;
/// Benz & Asphaug (1995).
struct TillotsonMaterial {
    double rho0;    ///< reference density
    double A;       ///< bulk modulus at rho0
    double B;       ///< nonlinear compression coefficient
    double u0;      ///< reference specific internal energy
    double a;       ///< Tillotson coefficient (compressed/cold)
    double b;       ///< Tillotson coefficient (compressed/cold)
    double alpha;   ///< exponential decay parameter (expanded regime)
    double beta;    ///< exponential decay parameter (expanded regime)
    double u_iv;    ///< specific energy of incipient vaporization
    double u_cv;    ///< specific energy of complete vaporization
    double mu;      ///< shear modulus (constitutive)
    // --- Solid-failure parameters (all default 0 => model inactive) ---
    double weibull_k;     ///< Weibull coefficient k [length^-DIM] (FRAGMENTATION)
    double weibull_m;     ///< Weibull exponent m (dimensionless)
    double Y0;            ///< cohesion / yield strength at zero pressure
    double mu_i;          ///< coefficient of internal friction (intact, Collins)
    double mu_d;          ///< coefficient of friction (damaged, Collins)
    double Y_M;           ///< von Mises limit / yield at high pressure
    double frictionAngle; ///< friction angle [rad] (Mohr-Coulomb / Drucker-Prager)
    double u_melt;        ///< melt specific energy (Collins melt degradation)
};
#endif

class EquationOfState {

public:
    EquationOfState(
#if EOS == 0
                const double hydro_gamma
#elif EOS == 1
                std::vector<MurnaghanMaterial> mats
#elif EOS == 2
                std::vector<TillotsonMaterial> mats
#endif // EOS
    );

#if EOS == 1 || EOS == 2
    double EOSPressure(int matId, const double &rho, const double &u);
    double EOSSoundSpeed(int matId, const double &rho, const double &u, const double &p);
    double EOSInternalEnergy(int matId, const double &rho, const double &p);
    double EOSEnergyFluxGamma(int matId, const double &rho, const double &p, const double &u);
    double EOSAdiabaticSoundSpeed(int matId, const double &rho, const double &p);
    double EOSGeneralGamma(int matId, const double &rho, const double &p);
    double EOSBulkModulus(int matId, const double &rho, const double &p);
    double EOSShearModulus(int matId);
#if EOS == 1
    const MurnaghanMaterial& EOSGetMaterial(int matId);
#else
    const TillotsonMaterial& EOSGetMaterial(int matId);
    /// Static (reference) Young's modulus E = 9Kmu/(3K+mu) with K = A (the
    /// reference bulk modulus). Used by the Grady-Kipp damage model.
    double EOSYoungModulus(int matId);
    /// Per-material yield strength Y(P, damage, u) selected at compile time by
    /// the *_PLASTICITY flags. `damage` is the full damage D in [0,1].
    /// Returns a large value when no per-material model is active.
    double EOSYieldStrength(int matId, const double &P, const double &damage,
                            const double &u);
#endif
#else
    double EOSPressure(const double &rho, const double &u);
    double EOSSoundSpeed(const double &rho, const double &u,
                        const double &p);
    double EOSInternalEnergy(const double &rho, const double &p);
    double EOSEnergyFluxGamma(const double &rho, const double &p, const double &u);

// For general EOS HLLC solver:
    double EOSAdiabaticSoundSpeed(const double &rho, const double &p);
    double EOSGeneralGamma(const double &rho, const double &p);
    double EOSBulkModulus(const double &rho, const double &p);
#endif // EOS

#if EOS == 0
    double EOSGetHydroGammaParam();
#endif

private:
#if EOS == 0
    const double hydro_gamma;
#elif EOS == 1
    std::vector<MurnaghanMaterial> materials;
#elif EOS == 2
    std::vector<TillotsonMaterial> materials;
#endif
};
#endif // MESHLESSHYDRO_EQUATIONOFSTATE_H
