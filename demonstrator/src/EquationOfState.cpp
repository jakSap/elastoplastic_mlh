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
    const double TIL_A, const double TIL_B, const double TIL_u0, const double TIL_a, const double TIL_b,
        const double TIL_alpha, const double TIL_beta, const double u_iv, const double TIL_u_cv
#endif // EOS
    ) :
# if EOS == 0
    hydro_gamma {hydro_gamma}
#elif EOS == 1
    materials { std::move(mats) }
#elif EOS == 2
    TIL_A {TIL_A}, TIL_B {TIL_B}, TIL_u0 {TIL_u0}, TIL_a {TIL_a}, TIL_b {TIL_b}, TIL_alpha {TIL_alpha}, TIL_beta {TIL_beta},
        TIL_u_iv {TIL_u_iv}, TIL_u_cv {TIL_u_cv}
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
        Logger(INFO) << "Using Tillotson EOS, parameter: ";
        Logger(INFO) << "TIL_A = " << TIL_A << "TIL_B =" << TIL_B << "TIL_u0 = " << TIL_u0;
        Logger(INFO) << "TIL_a = " << TIL_a << ", TIL_b = "<< TIL_b;
        Logger(INFO) << "TIL_alpha = " << TIL_alpha} << ", TIL_beta = " << TIL_beta ", TIL_u_iv =" <<TIL_u_iv << ", TIL_u_cv = " <<TIL_u_cv;
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
#else // EOS != 1

double EquationOfState::EOSPressure(const double &rho, const double &u){
#if EOS == 0 // Ideal gas
    return (hydro_gamma - 1) * u * rho;
#elif EOS == 2 // Tillotson
    return -1
#endif // EOS
}

double EquationOfState::EOSSoundSpeed(const double &rho, const double &u,
                        const double &p){
#if EOS == 0 // Ideal gas
    return sqrt(hydro_gamma * p / rho);
#elif EOS == 2 // Tillotson
    return -1; // TODO
#endif // EOS
}

double EquationOfState::EOSInternalEnergy(const double &rho, const double &p){
#if EOS == 0 // Ideal gas
    return p /((hydro_gamma - 1) * rho);
#endif // EOS
}

double EquationOfState::EOSEnergyFluxGamma(const double &rho, const double &p, const double &u){
#if EOS == 0
    return hydro_gamma;
#endif
}

double EquationOfState::EOSAdiabaticSoundSpeed(const double &rho, const double &p){
    double cs = 0;
#if EOS == 0
    cs = sqrt(hydro_gamma * p / rho);
#endif
    assert(cs >= 0 && "Negative sound speed encountered");
    assert(cs > 0 && "Zero sound speed encountered");
    return cs;
}

double EquationOfState::EOSBulkModulus(const double &rho, const double &p){
#if EOS == 0
    return hydro_gamma * p;
#endif
}

double EquationOfState::EOSGeneralGamma(const double &rho, const double &p){
#if EOS == 0
    return hydro_gamma;
#endif
}
#endif // EOS == 1

#if EOS == 0
double EquationOfState::EOSGetHydroGammaParam(){
    return hydro_gamma;
}
#endif
