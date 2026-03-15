//
// Created by Jakob Sappler on 18.09.25
//

#include "../include/EquationOfState.h"
EquationOfState::EquationOfState(
#if EOS == 0
    const double hydro_gamma
#elif EOS == 1
    const double K0, const double murn_n, const double rho0
#elif EOS == 2
    const double TIL_A, const double TIL_B, const double TIL_u0, const double TIL_a, const double TIL_b,
        const double TIL_alpha, const double TIL_beta, const double u_iv, const double TIL_u_cv
#endif // EOS
    ) :
# if EOS == 0
    hydro_gamma {hydro_gamma}
#elif EOS == 1
    K0 {K0}, murn_n {murn_n}, rho0 {rho0}
#elif EOS == 2
    TIL_A {TIL_A}, TIL_B {TIL_B}, TIL_u0 {TIL_u0}, TIL_a {TIL_a}, TIL_b {TIL_b}, TIL_alpha {TIL_alpha}, TIL_beta {TIL_beta},
        TIL_u_iv {TIL_u_iv}, TIL_u_cv {TIL_u_cv}
#endif //EOS
    {
#if EOS == 0
        Logger(INFO) << "Using ideal gas polytropic EOS with gamma = " << hydro_gamma;
#elif EOS == 1
        Logger(INFO) << "Using Murnaghan EOS, K0 = " << K0 <<
                            ", n = " << murn_n << ", rho0 = " << rho0;
#elif EOS == 2
        Logger(INFO) << "Using Tillotson EOS, parameter: ";
        Logger(INFO) << "TIL_A = " << TIL_A << "TIL_B =" << TIL_B << "TIL_u0 = " << TIL_u0;
        Logger(INFO) << "TIL_a = " << TIL_a << ", TIL_b = "<< TIL_b;
        Logger(INFO) << "TIL_alpha = " << TIL_alpha} << ", TIL_beta = " << TIL_beta ", TIL_u_iv =" <<TIL_u_iv << ", TIL_u_cv = " <<TIL_u_cv;
#endif //EOS
}

double EquationOfState::EOSPressure(const double &rho, const double &u){
#if EOS == 0 // Ideal gas
    return (hydro_gamma - 1) * u * rho;
#elif EOS == 1 // Murnaghan
    return K0 / murn_n * (pow(rho / rho0, murn_n) - 1);
#elif EOS == 2 // Tillotson
    return -1
#endif // EOS
}

double EquationOfState::EOSSoundSpeed(const double &rho, const double &u,
                        const double &p){
#if EOS == 0 // Ideal gas
    return sqrt(hydro_gamma * p / rho);
#elif EOS == 1 // Murnaghan
    return K0 / rho0 * pow(rho / rho0, murn_n - 1);
#elif EOS == 2 // Tillotson
    return -1; // TODO
#endif // EOS
}

double EquationOfState::EOSInternalEnergy(const double &rho, const double &p){
#if EOS == 0 // Ideal gas
    return p /((hydro_gamma - 1) * rho);
#elif EOS == 1 // Murnaghan
    return 1;
#endif // EOS
}

double EquationOfState::EOSEnergyFluxGamma(const double &rho, const double &p, const double &u){
#if EOS == 0
    return hydro_gamma;
#elif EOS == 1
    return K0 * pow(rho / rho0, murn_n) / p;
#endif
}

double EquationOfState::EOSAdiabaticSoundSpeed(const double &rho, const double &p){
    double cs = 0;
#if EOS == 0
    cs = sqrt(hydro_gamma * p / rho);
#elif EOS == 1
    const double eta = rho / rho0;
    cs = murn_n * pow(eta, murn_n) / (pow(eta, murn_n) - 1);
#endif
    assert(cs >= 0 && "Negative sound speed encountered");
    assert(cs > 0 && "Zero sound speed encountered");
    return cs;
}

double EquationOfState::EOSBulkModulus(const double &rho, const double &p){
#if EOS == 0
    return hydro_gamma * p;
#elif EOS == 1
    return K0 * pow(rho / rho0, murn_n);
#endif
}

// For general EOS: Compute Gamma = \partial ln(P) / \partial ln(rho)
double EquationOfState::EOSGeneralGamma(const double &rho, const double &p){
#if EOS == 0
    return hydro_gamma;
#elif EOS == 1
    return murn_n * pow(rho / rho0, murn_n) / (pow(rho / rho0, murn_n) - 1);
#endif
}

#if EOS == 0
double EquationOfState::EOSGetHydroGammaParam(){
    return hydro_gamma;
}

#elif EOS == 1
double EquationOfState::EOSGetMurnaghan_K0(){
    return K0;
}

double EquationOfState::EOSGetMurnaghan_n(){
    return murn_n;
}

double EquationOfState::EOSGetMurnaghan_rho0(){
    return rho0;
}
#endif
