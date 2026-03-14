//
// Created by Johannes Martin on 23.09.21.
//

#include "../include/InitialDistribution.h"

InitialDistribution::InitialDistribution(const std::string &file){
    HighFive::File h5file(file, HighFive::File::ReadOnly);

    // read datasets from file
    HighFive::DataSet mass = h5file.getDataSet("/m");
    HighFive::DataSet pos = h5file.getDataSet("/x");
    HighFive::DataSet vel = h5file.getDataSet("/v");
    HighFive::DataSet energy = h5file.getDataSet("/u");
#if USE_MATID
    HighFive::DataSet materialId = h5file.getDataSet("/materialId");
#endif

    // read data into containers
    mass.read(m);
    pos.read(x);
    vel.read(v);
    energy.read(u);
#if USE_MATID
    materialId.read(matId);

#endif
    // sanity check
    if (x.size() == v.size() && x.size() == m.size() && x.size() == u.size()
#if USE_MATID
    && x.size() == matId.size()
#endif
    ){
        numberOfParticles = x.size();
    } else {
        throw std::length_error("Length mismatch between mass, position and/or velocity vectors.");
    }
}

void InitialDistribution::getAllParticles(Particles &particles){
    std::vector<std::vector<double>>::iterator xit = x.begin();
    std::vector<std::vector<double>>::iterator vit = v.begin();
    std::vector<double>::iterator mit = m.begin();
    std::vector<double>::iterator uit = u.begin();
#if USE_MATID
    std::vector<int>::iterator matIdIt = matId.begin();

#endif
    int pCounter = 0;

    while (xit != x.end()){
        particles.m[pCounter] = *mit;
        particles.u[pCounter] = *uit;
        particles.x[pCounter] = (*xit)[0];
        particles.vx[pCounter] = (*vit)[0];
        particles.y[pCounter] = (*xit)[1];
        particles.vy[pCounter] = (*vit)[1];
#if USE_MATID
        particles.matId[pCounter] = *matIdIt;
#endif
#if DIM == 3
        particles.z[pCounter] = (*xit)[2];
        particles.vz[pCounter] = (*vit)[2];
#endif
#if USE_MATID
        ++matIdIt;
#endif
        ++uit;
        ++xit;
        ++vit;
        ++mit;
        ++pCounter;
    }
}
