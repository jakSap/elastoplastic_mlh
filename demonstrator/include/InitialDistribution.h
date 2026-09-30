//
// Created by Johannes Martin on 23.09.21.
//

#ifndef MESHLESSHYDRO_INITIALDISTRIBUTION_H
#define MESHLESSHYDRO_INITIALDISTRIBUTION_H

#include <highfive/H5File.hpp>

#include "Particles.h"

class InitialDistribution {
public:
    InitialDistribution(const std::string &file);

    int getNumberOfParticles() const { return numberOfParticles; };
    /// Copy IC data into `particles`. The smoothing length array is filled
    /// from the IC file if a `/sml` dataset is present; otherwise every
    /// particle gets `defaultSml` (typically `config.kernelSize`).
    void getAllParticles(Particles &particles, double defaultSml);

private:
    // containers to be filled from hdf5 file
    std::vector<double> m {}, u {}, sml {};
    std::vector<std::vector<double>> x {}, v {};
    std::vector<int> matId {};
    bool hasSml { false };
#if EXPLICIT_VOL_SEED_IC
    // /rho: density at t = 0 that seeds rhoExplicit (required under EXPLICIT_VOL_SEED_IC)
    std::vector<double> rhoIC {};
#endif
#if ELASTIC
    // Optional initial deviatoric stress, same dataset names as the snapshot
    // output (/Sxx, /Sxy, /Syy [, /Sxz, /Syz, /Szz]), so a snapshot can be fed
    // back in as an IC. All-or-nothing: absent -> S = 0, which is the right
    // start only for a state that is stress-free at t = 0. An elastic
    // travelling wave is not (S_xx = -4/3 mu ln(rho/rho0)), hence this hook.
    std::vector<double> Sxx {}, Sxy {}, Syy {};
#if DIM == 3
    std::vector<double> Sxz {}, Syz {}, Szz {};
#endif
    bool hasStress { false };
#endif
#if FRAGMENTATION
    // Per-particle Weibull activation strains (ragged, padded to MAX_NUM_FLAWS
    // on copy) and how many slots are populated per particle.
    std::vector<std::vector<double>> flaws {};
    std::vector<int> numFlaws {};
    bool hasFlaws { false };
#endif
    int numberOfParticles { 0 };
};


#endif //MESHLESSHYDRO_INITIALDISTRIBUTION_H
