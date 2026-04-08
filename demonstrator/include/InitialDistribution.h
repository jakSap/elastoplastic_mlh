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
    int numberOfParticles { 0 };
};


#endif //MESHLESSHYDRO_INITIALDISTRIBUTION_H
