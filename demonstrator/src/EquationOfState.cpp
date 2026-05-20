//
// Created by Jakob Sappler on 18.09.25
//

#include "../include/EquationOfState.h"

EquationOfState::EquationOfState(
#if EOS == 0
    const double hydro_gamma
#elif EOS == 1
    std::vector<MurnaghanMaterial> mats
#elif EOS == 2
    std::vector<TillotsonMaterial> mats
#endif // EOS
    ) :
# if EOS == 0
    hydro_gamma {hydro_gamma}
#elif EOS == 1
    materials { std::move(mats) }
#elif EOS == 2
    materials { std::move(mats) }
#endif //EOS
    {
#if EOS == 0
        Logger(INFO) << "Using ideal gas polytropic EOS with gamma = " << hydro_gamma;
#elif EOS == 1
        Logger(INFO) << "Using Murnaghan EOS with " << materials.size() << " material(s):";
        for (size_t k = 0; k < materials.size(); ++k){
            Logger(INFO) << "    > Material " << k
                         << ": K0=" << materials[k].K0
                         << ", n=" << materials[k].n
                         << ", rho0=" << materials[k].rho0
                         << ", mu=" << materials[k].mu;
        }
#elif EOS == 2
        Logger(INFO) << "Using Tillotson EOS with " << materials.size() << " material(s):";
        for (size_t k = 0; k < materials.size(); ++k){
            const TillotsonMaterial &m = materials[k];
            Logger(INFO) << "    > Material " << k
                         << ": rho0=" << m.rho0 << ", A=" << m.A << ", B=" << m.B
                         << ", u0=" << m.u0 << ", a=" << m.a << ", b=" << m.b
                         << ", alpha=" << m.alpha << ", beta=" << m.beta
                         << ", u_iv=" << m.u_iv << ", u_cv=" << m.u_cv
                         << ", mu=" << m.mu;
        }
#endif //EOS
}

#if EOS == 1
double EquationOfState::EOSPressure(int matId, const double &rho, const double &u){
#if DEBUG_LVL >= 1
    assert(matId >= 0 && matId < (int)materials.size());
#endif
    const MurnaghanMaterial &m = materials[matId];
    return m.K0 / m.n * (pow(rho / m.rho0, m.n) - 1);
}

double EquationOfState::EOSSoundSpeed(int matId, const double &rho, const double &u, const double &p){
#if DEBUG_LVL >= 1
    assert(matId >= 0 && matId < (int)materials.size());
#endif
    const MurnaghanMaterial &m = materials[matId];
    return m.K0 / m.rho0 * pow(rho / m.rho0, m.n - 1);
}

double EquationOfState::EOSInternalEnergy(int matId, const double &rho, const double &p){
#if DEBUG_LVL >= 1
    assert(matId >= 0 && matId < (int)materials.size());
#endif
    const MurnaghanMaterial &m = materials[matId];
    // u(rho) = integral_{rho0}^{rho} P(rho')/rho'^2 drho'
    if (std::abs(m.n - 1.0) < 1e-14) {
        // n = 1: u = K0/rho0 * ln(rho/rho0) + K0 * (1/rho - 1/rho0)
        return m.K0 / m.rho0 * log(rho / m.rho0) + m.K0 * (1.0 / rho - 1.0 / m.rho0);
    } else {
        // General n: u = K0/(n(n-1)rho0^n) * (rho^(n-1) - rho0^(n-1)) + K0/n * (1/rho - 1/rho0)
        return m.K0 / (m.n * (m.n - 1.0) * pow(m.rho0, m.n)) * (pow(rho, m.n - 1.0) - pow(m.rho0, m.n - 1.0))
             + m.K0 / m.n * (1.0 / rho - 1.0 / m.rho0);
    }
}

double EquationOfState::EOSEnergyFluxGamma(int matId, const double &rho, const double &p, const double &u){
#if DEBUG_LVL >= 1
    assert(matId >= 0 && matId < (int)materials.size());
#endif
    const MurnaghanMaterial &m = materials[matId];
    return m.K0 * pow(rho / m.rho0, m.n) / p;
}

double EquationOfState::EOSAdiabaticSoundSpeed(int matId, const double &rho, const double &p){
#if DEBUG_LVL >= 1
    assert(matId >= 0 && matId < (int)materials.size());
#endif
    const MurnaghanMaterial &m = materials[matId];
    const double eta = rho / m.rho0;
    double cs = m.n * pow(eta, m.n) / (pow(eta, m.n) - 1);
    assert(cs >= 0 && "Negative sound speed encountered");
    assert(cs > 0 && "Zero sound speed encountered");
    return cs;
}

double EquationOfState::EOSGeneralGamma(int matId, const double &rho, const double &p){
#if DEBUG_LVL >= 1
    assert(matId >= 0 && matId < (int)materials.size());
#endif
    const MurnaghanMaterial &m = materials[matId];
    return m.n * pow(rho / m.rho0, m.n) / (pow(rho / m.rho0, m.n) - 1);
}

double EquationOfState::EOSBulkModulus(int matId, const double &rho, const double &p){
#if DEBUG_LVL >= 1
    assert(matId >= 0 && matId < (int)materials.size());
#endif
    const MurnaghanMaterial &m = materials[matId];
    return m.K0 * pow(rho / m.rho0, m.n);
}

double EquationOfState::EOSShearModulus(int matId){
#if DEBUG_LVL >= 1
    assert(matId >= 0 && matId < (int)materials.size());
#endif
    return materials[matId].mu;
}

const MurnaghanMaterial& EquationOfState::EOSGetMaterial(int matId){
#if DEBUG_LVL >= 1
    assert(matId >= 0 && matId < (int)materials.size());
#endif
    return materials[matId];
}

#elif EOS == 2
// ============================================================================
// Tillotson EOS implementation. References: Tillotson (1962); Melosh,
// "Impact Cratering", Appendix B; Benz & Asphaug (1995). Piecewise across
// three regimes (compressed/cold, hot expanded, mixed). The shared helpers
// below return P together with dP/drho|_u and dP/du|_rho so the analytic
// sound-speed expression and the Newton inverter all use one code path.
// ============================================================================

namespace {

inline void tillotsonCompressed(const TillotsonMaterial &m,
                                 const double rho, const double u,
                                 double &P, double &dPdrho, double &dPdu) {
    const double eta = rho / m.rho0;
    const double mu  = eta - 1.0;
    const double inv_u0_eta2 = 1.0 / (m.u0 * eta * eta);
    const double omega = u * inv_u0_eta2 + 1.0;
    const double inv_omega  = 1.0 / omega;
    const double inv_omega2 = inv_omega * inv_omega;
    const double om_m1 = omega - 1.0;

    P      = (m.a + m.b * inv_omega) * rho * u + m.A * mu + m.B * mu * mu;
    dPdu   = (m.a + m.b * inv_omega2) * rho;
    dPdrho = (m.a + m.b * inv_omega) * u
             + 2.0 * m.b * u * om_m1 * inv_omega2
             + (m.A + 2.0 * m.B * mu) / m.rho0;
}

inline void tillotsonExpanded(const TillotsonMaterial &m,
                               const double rho, const double u,
                               double &P, double &dPdrho, double &dPdu) {
    const double eta = rho / m.rho0;
    const double mu  = eta - 1.0;
    const double inv_u0_eta2 = 1.0 / (m.u0 * eta * eta);
    const double omega = u * inv_u0_eta2 + 1.0;
    const double inv_omega  = 1.0 / omega;
    const double inv_omega2 = inv_omega * inv_omega;
    const double om_m1 = omega - 1.0;

    const double z    = 1.0 / eta - 1.0;          // rho0/rho - 1
    const double zp1  = z + 1.0;                  // 1/eta
    const double exp1 = std::exp(-m.beta * z);
    const double exp2 = std::exp(-m.alpha * z * z);

    const double Q = m.b * rho * u * inv_omega + m.A * mu * exp1;
    P = m.a * rho * u + Q * exp2;

    dPdu = m.a * rho + (m.b * rho * inv_omega2) * exp2;

    const double dQdrho = m.b * u * inv_omega
                          + 2.0 * m.b * u * om_m1 * inv_omega2
                          + m.A * exp1 / m.rho0
                          + m.A * mu * m.beta * exp1 * zp1 / rho;
    const double dexp2_drho = 2.0 * m.alpha * z * zp1 / rho * exp2;
    dPdrho = m.a * u + dQdrho * exp2 + Q * dexp2_drho;
}

inline void tillotsonPressureAndDerivs(const TillotsonMaterial &m,
                                        const double rho, const double u,
                                        double &P, double &dPdrho, double &dPdu) {
    const double eta = rho / m.rho0;
    if (eta >= 1.0 || u <= m.u_iv) {
        tillotsonCompressed(m, rho, u, P, dPdrho, dPdu);
    } else if (u >= m.u_cv) {
        tillotsonExpanded(m, rho, u, P, dPdrho, dPdu);
    } else {
        double Pc, dPcdrho, dPcdu;
        double Pe, dPedrho, dPedu;
        tillotsonCompressed(m, rho, u, Pc, dPcdrho, dPcdu);
        tillotsonExpanded(m, rho, u, Pe, dPedrho, dPedu);
        const double inv_du = 1.0 / (m.u_cv - m.u_iv);
        const double w  = (u - m.u_iv) * inv_du;
        P      = w * Pe + (1.0 - w) * Pc;
        dPdrho = w * dPedrho + (1.0 - w) * dPcdrho;
        dPdu   = w * dPedu   + (1.0 - w) * dPcdu + (Pe - Pc) * inv_du;
    }
}

/// Newton-Raphson solve P_tillotson(rho, u) = p_target for u.
/// Seeded from a compressed-regime linearization. Used by every method that
/// receives (rho, p) but needs u — namely the post-HLLC reconstruction path
/// and the IC-init pressure-to-energy step.
double tillotsonInternalEnergy(const TillotsonMaterial &m,
                                const double rho, const double p_target) {
    const double eta = rho / m.rho0;
    const double mu  = eta - 1.0;
    // Cold-region linearization: ignore the b/omega term to get an initial
    // guess. P_c ~ (a + b) * rho * u + A*mu + B*mu^2  when omega ~ 1.
    double u = (p_target - m.A * mu - m.B * mu * mu) / ((m.a + m.b) * rho);
    if (!std::isfinite(u) || u < 0.0) u = 0.0;

    const double tol = 1e-10 * std::max(std::abs(p_target), 1.0);
    for (int iter = 0; iter < 25; ++iter) {
        double P, dPdrho, dPdu;
        tillotsonPressureAndDerivs(m, rho, u, P, dPdrho, dPdu);
        const double resid = P - p_target;
        if (std::abs(resid) < tol) return u;
        if (std::abs(dPdu) < 1e-300) break;
        double du = -resid / dPdu;
        // Damped step: keep u >= 0 (Tillotson is undefined for u < 0).
        if (u + du < 0.0) du = -0.5 * u;
        u += du;
    }
    return u;
}

inline double tillotsonSoundSpeedSquared(const TillotsonMaterial &m,
                                          const double rho, const double u) {
    double P, dPdrho, dPdu;
    tillotsonPressureAndDerivs(m, rho, u, P, dPdrho, dPdu);
    double cs2 = dPdrho + P / (rho * rho) * dPdu;
    if (cs2 < 0.0) cs2 = 0.0;
    return cs2;
}

} // anonymous namespace

double EquationOfState::EOSPressure(int matId, const double &rho, const double &u){
#if DEBUG_LVL >= 1
    assert(matId >= 0 && matId < (int)materials.size());
#endif
    const TillotsonMaterial &m = materials[matId];
    double P, dPdrho, dPdu;
    tillotsonPressureAndDerivs(m, rho, u, P, dPdrho, dPdu);
    return P;
}

double EquationOfState::EOSSoundSpeed(int matId, const double &rho, const double &u, const double &p){
#if DEBUG_LVL >= 1
    assert(matId >= 0 && matId < (int)materials.size());
#endif
    const TillotsonMaterial &m = materials[matId];
    // u < 0 is the caller's sentinel for "use p instead" (see Particles.cpp:2148).
    const double u_use = (u >= 0.0) ? u : tillotsonInternalEnergy(m, rho, p);
    return tillotsonSoundSpeedSquared(m, rho, u_use);
}

double EquationOfState::EOSInternalEnergy(int matId, const double &rho, const double &p){
#if DEBUG_LVL >= 1
    assert(matId >= 0 && matId < (int)materials.size());
#endif
    return tillotsonInternalEnergy(materials[matId], rho, p);
}

double EquationOfState::EOSEnergyFluxGamma(int matId, const double &rho, const double &p, const double &u){
#if DEBUG_LVL >= 1
    assert(matId >= 0 && matId < (int)materials.size());
#endif
    const TillotsonMaterial &m = materials[matId];
    const double cs2 = tillotsonSoundSpeedSquared(m, rho, u);
    return cs2 * rho / p;
}

double EquationOfState::EOSAdiabaticSoundSpeed(int matId, const double &rho, const double &p){
#if DEBUG_LVL >= 1
    assert(matId >= 0 && matId < (int)materials.size());
#endif
    const TillotsonMaterial &m = materials[matId];
    const double u   = tillotsonInternalEnergy(m, rho, p);
    const double cs2 = tillotsonSoundSpeedSquared(m, rho, u);
    return cs2 * rho / p;
}

double EquationOfState::EOSGeneralGamma(int matId, const double &rho, const double &p){
#if DEBUG_LVL >= 1
    assert(matId >= 0 && matId < (int)materials.size());
#endif
    const TillotsonMaterial &m = materials[matId];
    const double u   = tillotsonInternalEnergy(m, rho, p);
    const double cs2 = tillotsonSoundSpeedSquared(m, rho, u);
    return cs2 * rho / p;
}

double EquationOfState::EOSBulkModulus(int matId, const double &rho, const double &p){
#if DEBUG_LVL >= 1
    assert(matId >= 0 && matId < (int)materials.size());
#endif
    const TillotsonMaterial &m = materials[matId];
    const double u   = tillotsonInternalEnergy(m, rho, p);
    const double cs2 = tillotsonSoundSpeedSquared(m, rho, u);
    return rho * cs2;
}

double EquationOfState::EOSShearModulus(int matId){
#if DEBUG_LVL >= 1
    assert(matId >= 0 && matId < (int)materials.size());
#endif
    return materials[matId].mu;
}

const TillotsonMaterial& EquationOfState::EOSGetMaterial(int matId){
#if DEBUG_LVL >= 1
    assert(matId >= 0 && matId < (int)materials.size());
#endif
    return materials[matId];
}

#else // EOS == 0

double EquationOfState::EOSPressure(const double &rho, const double &u){
    return (hydro_gamma - 1) * u * rho;
}

double EquationOfState::EOSSoundSpeed(const double &rho, const double &u,
                        const double &p){
    return sqrt(hydro_gamma * p / rho);
}

double EquationOfState::EOSInternalEnergy(const double &rho, const double &p){
    return p /((hydro_gamma - 1) * rho);
}

double EquationOfState::EOSEnergyFluxGamma(const double &rho, const double &p, const double &u){
    return hydro_gamma;
}

double EquationOfState::EOSAdiabaticSoundSpeed(const double &rho, const double &p){
    double cs = sqrt(hydro_gamma * p / rho);
    assert(cs >= 0 && "Negative sound speed encountered");
    assert(cs > 0 && "Zero sound speed encountered");
    return cs;
}

double EquationOfState::EOSBulkModulus(const double &rho, const double &p){
    return hydro_gamma * p;
}

double EquationOfState::EOSGeneralGamma(const double &rho, const double &p){
    return hydro_gamma;
}

double EquationOfState::EOSGetHydroGammaParam(){
    return hydro_gamma;
}
#endif // EOS
