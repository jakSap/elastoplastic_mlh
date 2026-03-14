//
// Created by Jakob Sappler on 18.09.25
//

#ifndef MESHLESSHYDRO_EQUATIONOFSTATE_H
#define MESHLESSHYDRO_EQUATIONOFSTATE_H

#include <cmath>
#include <cassert>

#include "parameter.h"
#include "Logger.h"

class EquationOfState {

public:
    EquationOfState(
#if EOS == 0
                const double hydro_gamma
#elif EOS == 1
                const double K0, const double murn_n, const double rho0
#endif // EOS
    );

    double EOSPressure(const double &rho, const double &u);
    double EOSSoundSpeed(const double &rho, const double &u,
                        const double &p);
    double EOSInternalEnergy(const double &rho, const double &p);
    double EOSEnergyFluxGamma(const double &rho, const double &p, const double &u);

// For general EOS HLLC solver:
    double EOSAdiabaticSoundSpeed(const double &rho, const double &p);
    double EOSGeneralGamma(const double &rho, const double &p);
#if EOS == 0
    double EOSGetHydroGammaParam();
#elif EOS == 1
    double EOSGetMurnaghan_K0();
    double EOSGetMurnaghan_n();
    double EOSGetMurnaghan_rho0();
#endif

private:
#if EOS == 0
    const double hydro_gamma;
#elif EOS == 1
    const double K0, murn_n, rho0;
#endif
};
#endif // MESHLESSHYDRO_EQUATIONOFSTATE_H