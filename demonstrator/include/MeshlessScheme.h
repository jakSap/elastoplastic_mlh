//
// Created by Johannes Martin on 21.09.22.
//

#ifndef DEMONSTRATOR_MESHLESSSCHEME_H
#define DEMONSTRATOR_MESHLESSSCHEME_H

#include <algorithm>
#include <iomanip>

#include "parameter.h"
#include "InitialDistribution.h"
#include "Logger.h"
#include "Domain.h"
#include "Helper.h"


class MeshlessScheme {

public:
    struct Configuration {
        std::string initFile;
        std::string outDir;
        double timeStep;
        double timeEnd;
        int h5DumpInterval;
        double periodicBoxLimits[2 * DIM];
        double kernelSize;
#if VARIABLE_SML
        // Per-particle smoothing length iteration parameters.
        // Defaults come from parameter.h; the config file may override them.
        double smlNNNTarget { (double)NNN_TARGET };
        double smlTol       { SML_TOL };
        int    smlMaxIter   { SML_MAX_ITER };
        double smlMaxFactor { SML_MAX_FACTOR };
        // Safeguards: hard bounds as factors of kernelSize, and the fractions
        // of bad-state particles that trigger WARN / ERROR+exit.
        double smlHMinFactor   { SML_H_MIN_FACTOR };
        double smlHMaxFactor   { SML_H_MAX_FACTOR };
        double smlWarnFraction { SML_WARN_FRACTION };
        double smlPanicFraction{ SML_PANIC_FRACTION };
#endif
#if EOS == 0
        double hydro_gamma; // adiabatic index
#elif EOS == 1
        double K0;
        double murn_n;
        double rho0;
#endif
    };

    MeshlessScheme(Configuration config, Particles *particles, Domain::Cell domain);
    ~MeshlessScheme();

    void run();

private:
    Configuration config;
    double timeStep;
    Particles *particles;
#if PERIODIC_BOUNDARIES
    Particles ghostParticles;
#endif
    Domain domain;
    Helper helper {};

};


#endif //DEMONSTRATOR_MESHLESSSCHEME_H
