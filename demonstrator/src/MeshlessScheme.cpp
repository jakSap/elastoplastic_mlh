//
// Created by Johannes Martin on 21.09.22.
//

#include <chrono>

#include "../include/MeshlessScheme.h"


MeshlessScheme::MeshlessScheme(Configuration config, Particles *particles,
                               Domain::Cell bounds) : config { config }, particles { particles },
                                                      domain(bounds)
#if PERIODIC_BOUNDARIES
                                                      , ghostParticles(DIM*(particles->N),
                                                                        particles->MeshlessEOS,
                                                                        true) // TODO: memory optimization
#endif
                                                      {

    Logger(INFO) << "    > Creating grid ... ";
    domain.createGrid(config.kernelSize);
    Logger(INFO) << "    > ... got " << domain.numGridCells << " cells";
#if VARIABLE_SML
    particles->smlNNNTarget = config.smlNNNTarget;
    particles->smlTol       = config.smlTol;
    particles->smlMaxIter   = config.smlMaxIter;
    particles->smlMaxFactor = config.smlMaxFactor;
    // Resolve the hMin/hMax factors against the configured kernel size so
    // Particles stores absolute bounds and does not need to know kernelSize.
    particles->smlHMin          = config.smlHMinFactor * config.kernelSize;
    particles->smlHMax          = config.smlHMaxFactor * config.kernelSize;
    particles->smlH0            = config.kernelSize;
    particles->smlWarnFraction  = config.smlWarnFraction;
    particles->smlPanicFraction = config.smlPanicFraction;
#endif
}

void MeshlessScheme::run(){
    auto tStart = std::chrono::steady_clock::now();
    double t = 0;
    int step = 0;

#if ADAPTIVE_TIMESTEP
    int numDumpTimes = (int)(config.timeEnd/config.timeStep)/config.h5DumpInterval+1;
    Logger(DEBUG) << "      > Times for file dump: " << numDumpTimes;
    double dumpTimes[numDumpTimes];
    for(int iDump=0; iDump<numDumpTimes; ++iDump){
        dumpTimes[iDump] = iDump*config.timeStep*config.h5DumpInterval;
        // Logger(DEBUG) << "        dumpTimes[" << iDump << "] = " << dumpTimes[iDump];
    }
    int dumpStep = 0;
    bool dump = true;
    bool dumpNext = false;
#endif // ADAPTIVE_TIMESTEP


    do {
        Logger(INFO) << "  > TIME: " << t << ", STEP: " << step;

        // Search radius used for the grid, NNS and ghost-particle generation.
        // With VARIABLE_SML, h_i may grow beyond config.kernelSize. To keep
        // neighbor relations symmetric (i in nnl[j] iff j in nnl[i]), we
        // search with the largest h currently in use; flux symmetry is then
        // preserved by the existing ENFORCE_FLUX_SYM machinery.
#if VARIABLE_SML
        const double nnsRadius = std::max(config.kernelSize, particles->hMax());
#else
        const double nnsRadius = config.kernelSize;
#endif

#if !PERIODIC_BOUNDARIES
        Logger(INFO) << "    > Computing domain limits ...";
        double domainLimits[DIM*2];
        particles->getDomainLimits(domainLimits);
        Domain::Cell boundingBox { domainLimits };
        domain.bounds = boundingBox;
        domain.printout();
        Logger(DEBUG) << "      > ... creating grid ...";
        domain.createGrid(nnsRadius);
        Logger(INFO) << "    > ... done.";
#else // PERIODIC_BOUNDARIES
#if VARIABLE_SML
        // The periodic grid is built once in the constructor with cell size
        // based on config.kernelSize. If h_max grows past the current cell
        // size, we must rebuild the grid so that gridNNS' 3^DIM cell scan
        // still covers the search radius.
        if (nnsRadius > domain.cellSizeX){
            Logger(DEBUG) << "      > h_max=" << nnsRadius
                          << " exceeds cellSizeX=" << domain.cellSizeX
                          << ", rebuilding periodic grid";
            domain.createGrid(nnsRadius);
        }
#endif // VARIABLE_SML
#endif
        Logger(INFO) << "    > Assigning particles ...";
        particles->assignParticlesAndCells(domain);
        Logger(INFO) << "    > ... done.";
#if PERIODIC_BOUNDARIES
        Logger(INFO) << "    > Creating ghost particles ...";
        //Logger(DEBUG) << "      > Creating ghost grid";
        //domain.createGhostGrid();
        Logger(DEBUG) << "      > Creating ghost particles ... ";
        particles->createGhostParticles(domain, ghostParticles, nnsRadius);
        Logger(DEBUG) << "      > ... found " << ghostParticles.N << " ghosts";
        Logger(INFO) << "    > ... done.";

#endif
        Logger(INFO) << "    > Nearest neighbor search";
        particles->gridNNS(domain, nnsRadius);
#if PERIODIC_BOUNDARIES
        Logger(DEBUG) << "      > Ghosts NNS";
        particles->ghostNNS(domain, ghostParticles, nnsRadius);
        //particles->printNoi();
#endif
#if VARIABLE_SML
        Logger(INFO) << "    > Updating smoothing lengths";
#if PERIODIC_BOUNDARIES
        particles->updateAllSmoothingLengths(ghostParticles);
#else
        particles->updateAllSmoothingLengths();
#endif
#endif // VARIABLE_SML
        Logger(INFO) << "    > Computing density";
        particles->compDensity();
#if PERIODIC_BOUNDARIES
        particles->compDensity(ghostParticles);
#endif

#if EXPLICIT_VOL_INTEGRATION
        // GIZMO HYDRO_EXPLICITLY_INTEGRATE_VOLUME: from step 1 onward, rho
        // downstream is the explicitly-integrated value, not the kernel sum.
        Logger(DEBUG) << "      > Explicit-volume override (substitute rho)";
        particles->applyExplicitVolumeOverride();
#endif

        Logger(INFO) << "    > Computing pressure";
        particles->compPressure();
        //particles->printDensity(config.gamma);

        Logger(DEBUG) << "      SANITY CHECK > V_tot = " << particles->sumVolume();
        Logger(DEBUG) << "      SANITY CHECK > M_tot = " << particles->sumMass();
        Logger(DEBUG) << "      SANITY CHECK > E_tot = " << particles->sumEnergy();
        Logger(DEBUG) << "      SANITY CHECK > px_tot = " << particles->sumMomentumX();
        Logger(DEBUG) << "      SANITY CHECK > py_tot = " << particles->sumMomentumY();
#if DIM == 3
        Logger(DEBUG) << "      SANITY CHECK > pz_tot = " << particles->sumMomentumZ();
#endif

#if ADAPTIVE_TIMESTEP
        Logger(INFO) << "    > Selecting global timestep ... ";
        timeStep = particles->compGlobalTimestep();
        //Logger(INFO) << "Time  > dt = " << timeStep << " selected.";
        if(dumpStep >= numDumpTimes){
            Logger(ERROR) << "Simulation did not abort after reaching timeEnd. Exiting.";
            exit(9);
        } else if((t+timeStep>=dumpTimes[dumpStep+1])&(t<dumpTimes[dumpStep+1]))
	{

            dumpNext = true;
            timeStep = dumpTimes[dumpStep+1]-t;
            Logger(INFO) << "Shorter timestep for punctual dumping";
        }
        Logger(INFO) << "Time  > dt = " << timeStep << " selected.";
#else // ADAPTIVE_TIMESTEP
        timeStep = config.timeStep;
#endif // ADAPTIVE_TIMESTEP

        Logger(INFO) << "    > Computing gradients";
#if PERIODIC_BOUNDARIES
        particles->updateGhostState(ghostParticles);
        particles->compPsijTilde(helper, ghostParticles);
        //Logger(DEBUG) << "      > Update ghost psij_xiTilde";
        //particles->updateGhostPsijTilde(ghostParticles);

        particles->gradient(particles->rho, particles->rhoGrad, ghostParticles.rho, ghostParticles);
        particles->gradient(particles->vx, particles->vxGrad, ghostParticles.vx, ghostParticles);
        particles->gradient(particles->vy, particles->vyGrad, ghostParticles.vy, ghostParticles);
#if DIM == 3
        particles->gradient(particles->vz, particles->vzGrad, ghostParticles.vz, ghostParticles);
#endif
        particles->gradient(particles->P, particles->PGrad, ghostParticles.P, ghostParticles);
#if ELASTIC
        // Stress tensor gradients
        particles->gradient(particles->Sxx, particles->SxxGrad, ghostParticles.Sxx, ghostParticles);
        particles->gradient(particles->Sxy, particles->SxyGrad, ghostParticles.Sxy, ghostParticles);
        particles->gradient(particles->Syy, particles->SyyGrad, ghostParticles.Syy, ghostParticles);
#if DIM == 3
        particles->gradient(particles->Sxz, particles->SxzGrad, ghostParticles.Sxz, ghostParticles);
        particles->gradient(particles->Syz, particles->SyzGrad, ghostParticles.Syz, ghostParticles);
        particles->gradient(particles->Szz, particles->SzzGrad, ghostParticles.Szz, ghostParticles);
#endif // DIM == 3
#endif // ELASTIC
        Logger(DEBUG) << "      > Update ghost gradients";
        particles->updateGhostGradients(ghostParticles);

#if SLOPE_LIMITING
        // TODO: Check slope limiter
        Logger(DEBUG) << "      > Limiting slopes";
        particles->slopeLimiter(&ghostParticles);
#if ELASTIC
        particles->slopeLimiterStress(&ghostParticles);
#endif
        Logger(DEBUG) << "      > Update limited ghost gradients";
        particles->updateGhostGradients(ghostParticles);
#endif // SLOPE_LIMITING
#else // PERIODIC_BOUNDARIES
        particles->compPsijTilde(helper);

#if OUTPUT_CONDITION_NUMBER
        // Diagnostic: identify the particles with the worst (largest)
        // condition number. A high conditionNumber means the E matrix is
        // near-singular and the gradient estimator is unreliable. We want
        // to correlate bad gradients with the EXTREME reconstruction events.
        {
            const int TOP = 5;
            int    worst[TOP];
            double worstVal[TOP];
            for (int k = 0; k < TOP; ++k) { worst[k] = -1; worstVal[k] = -1.; }
            for (int i = 0; i < particles->N; ++i){
                double c = particles->conditionNumber[i];
                if (!(c > worstVal[TOP-1])) continue;
                int k = TOP - 1;
                while (k > 0 && c > worstVal[k-1]){
                    worst[k] = worst[k-1];
                    worstVal[k] = worstVal[k-1];
                    --k;
                }
                worst[k] = i;
                worstVal[k] = c;
            }
            for (int k = 0; k < TOP; ++k){
                int i = worst[k];
                if (i < 0) break;
                Logger(DEBUG) << "      > worst cond[" << k << "] @ i=" << i
                              << " cond=" << worstVal[k]
                              << " x=[" << particles->x[i] << "," << particles->y[i] << "]"
                              << " v=[" << particles->vx[i] << "," << particles->vy[i] << "]"
                              << " rho=" << particles->rho[i]
                              << " sml=" << particles->sml[i];
            }
        }
#endif

        particles->gradient(particles->rho, particles->rhoGrad);
        particles->gradient(particles->vx, particles->vxGrad);
        particles->gradient(particles->vy, particles->vyGrad);
#if DIM == 3
        particles->gradient(particles->vz, particles->vzGrad);
#endif
        particles->gradient(particles->P, particles->PGrad);
#if ELASTIC
        particles->gradient(particles->Sxx, particles->SxxGrad);
        particles->gradient(particles->Sxy, particles->SxyGrad);
        particles->gradient(particles->Syy, particles->SyyGrad);
#if DIM == 3
        particles->gradient(particles->Sxz, particles->SxzGrad);
        particles->gradient(particles->Syz, particles->SyzGrad);
        particles->gradient(particles->Szz, particles->SzzGrad);
#endif
#endif // ELASTIC
#if SLOPE_LIMITING
        // TODO: check how to properly limit gradiens
        Logger(DEBUG) << "      > Limiting slopes";
        particles->slopeLimiter();
#if ELASTIC
        particles->slopeLimiterStress();
#endif // ELASTIC
#endif // SLOPE_LIMITING
#endif
#if ELASTIC
	Logger(INFO) << "    > Performing stress integration 1 / 2";
        particles->integrateStressTensor(timeStep / 2.0);
#if FRAGMENTATION
        particles->integrateDamage(timeStep / 2.0);
#endif
#if PERIODIC_BOUNDARIES
        Logger(DEBUG) << "      > Update ghost state after stress integration";
        particles->updateGhostState(ghostParticles);
#endif
#endif // ELASTIC
#if EXPLICIT_VOL_INTEGRATION
        // Half-step kick A: advance rhoExplicit by dt/2 using pre-flux gradients.
        // finalize=false: under EXPLICIT_VOL_FREEZE_RHO the working rho is left
        // frozen for the flux pass (GIZMO mode==0).
        Logger(DEBUG) << "      > Explicit-volume half-step kick A";
        particles->integrateExplicitVolumeHalfStep(timeStep / 2.0, /*finalize=*/false);
#endif
        Logger(INFO) << "    > Preparing Riemann solver";
        Logger(DEBUG) << "      > Computing effective faces";
        particles->compEffectiveFace();
#if PERIODIC_BOUNDARIES
        particles->compEffectiveFace(ghostParticles);
#endif // PERIODIC_BOUNDARIES

#if TENSILE_CORRECTION
        Logger(INFO) << "       > Computing fabMonaghan";
        particles->computeFabMonaghan();
#endif

        Logger(DEBUG) << "      > Computing fluxes";
        particles->compRiemannStatesLR(timeStep);

#if PERIODIC_BOUNDARIES
        Logger(DEBUG) << "      > Computing ghost fluxes";
        particles->compRiemannStatesLR(timeStep,
                                     ghostParticles);
        //Logger(DEBUG) << "Aborting for debugging.";
        //exit(6);

#endif// PERIODIC_BOUNDARIES

#if ADAPTIVE_TIMESTEP
        if (dump){
            dump = false;

#else // ADAPTIVE_TIMESTEP
        if (step % config.h5DumpInterval == 0) {
#endif // ADAPTIVE_TIMESTEP

            std::stringstream stepss;
            Logger(INFO) << "   > Dump particle distribution";

            stepss << std::setw(6) << std::setfill('0')
#if ADAPTIVE_TIMESTEP
            << dumpStep;
#else
            << step;
#endif // ADAPTIVE_TIMESTEP

            Logger(INFO) << "      > Dump particles to file";
            particles->dump2file(config.outDir + "/" + stepss.str() + std::string(".h5"), t);
// ???
#if ADAPTIVE_TIMESTEP
            ++dumpStep;
#else
            ++ step;
#endif // ADAPTIVE_TIMESTEP

#if DEBUG_LVL > 1
#if PERIODIC_BOUNDARIES
            Logger(INFO) << "      > Dump ghosts to file";
            ghostParticles.dump2file(config.outDir + "/" + stepss.str() + std::string("Ghosts.h5"), t);
            Logger(INFO) << "      > Dump NNL to file";
            particles->dumpNNL(config.outDir + "/" + stepss.str() + std::string("NNL.h5"), ghostParticles);
#else
            Logger(INFO) << "      > Dump NNL to file";
            particles->dumpNNL(config.outDir + "/" + stepss.str() + std::string("NNL.h5"));
#endif // PERIODIC_BOUNDARIES
#endif // DEBUG_LVL
        }
        if (t>=config.timeEnd){
            Logger(INFO) << "    > t = " << t << " -> FINISHED!";
            break;
        }


        Logger(INFO) << "    > Solving Riemann problems";
#if PERIODIC_BOUNDARIES
        particles->solveRiemannProblems(ghostParticles);
#else
        Particles ghostParticlesDummy { 0, particles->MeshlessEOS, true };
        particles->solveRiemannProblems(ghostParticlesDummy);
#endif

#if DEBUG_LVL
        Logger(DEBUG) << "    > Checking flux symmetry";
#if PERIODIC_BOUNDARIES
        particles->checkFluxSymmetry(&ghostParticles);
#else
        particles->checkFluxSymmetry();
#endif
#endif

        Logger(INFO) << "    > Collecting fluxes";

        particles->collectFluxes(helper);
#if !ELASTIC
	Logger(INFO) << "    > Updating state and position";
        particles->updateStateAndPosition(timeStep, domain);
#else // !ELASTIC

        Logger(INFO) << "    > Updating state";
        particles->updateState(timeStep);
        Logger(INFO) << "    > Recomputing velocity and stress gradients";
#if PERIODIC_BOUNDARIES
        particles->updateGhostState(ghostParticles);

        particles->gradient(particles->vx, particles->vxGrad, ghostParticles.vx, ghostParticles);
        particles->gradient(particles->vy, particles->vyGrad, ghostParticles.vy, ghostParticles);
#if DIM == 3
        particles->gradient(particles->vz, particles->vzGrad, ghostParticles.vz, ghostParticles);
#endif
        particles->gradient(particles->Sxx, particles->SxxGrad, ghostParticles.Sxx, ghostParticles);
        particles->gradient(particles->Sxy, particles->SxyGrad, ghostParticles.Sxy, ghostParticles);
        particles->gradient(particles->Syy, particles->SyyGrad, ghostParticles.Syy, ghostParticles);
#if DIM == 3
        particles->gradient(particles->Sxz, particles->SxzGrad, ghostParticles.Sxz, ghostParticles);
        particles->gradient(particles->Syz, particles->SyzGrad, ghostParticles.Syz, ghostParticles);
        particles->gradient(particles->Szz, particles->SzzGrad, ghostParticles.Szz, ghostParticles);
#endif

        particles->updateGhostGradients(ghostParticles);

#if SLOPE_LIMITING
        Logger(DEBUG) << "      > Limiting new gradients";
        particles->slopeLimiterStress(&ghostParticles);
        particles->updateGhostGradients(ghostParticles);
#endif // SLOPE_LIMITING

#else // !PERIODIC_BOUNDARIES
        particles->gradient(particles->vx, particles->vxGrad);
        particles->gradient(particles->vy, particles->vyGrad);
#if DIM == 3
        particles->gradient(particles->vz, particles->vzGrad);
#endif
        particles->gradient(particles->Sxx, particles->SxxGrad);
        particles->gradient(particles->Sxy, particles->SxyGrad);
        particles->gradient(particles->Syy, particles->SyyGrad);
#if DIM == 3
        particles->gradient(particles->Sxz, particles->SxzGrad);
        particles->gradient(particles->Syz, particles->SyzGrad);
        particles->gradient(particles->Szz, particles->SzzGrad);
#endif
#if SLOPE_LIMITING
        particles->slopeLimiterStress();
#endif
#endif // PERIODIC_BOUNDARIES

        Logger(INFO) << "    > Performing stress integration 2 / 2";
        particles->integrateStressTensor(timeStep / 2.0);
#if FRAGMENTATION
        particles->integrateDamage(timeStep / 2.0);
#endif

#if EXPLICIT_VOL_INTEGRATION
        // Half-step kick B: advance rhoExplicit by dt/2 using post-update gradients.
        // finalize=true: end-of-step swap rho <- rhoExplicit (GIZMO mode==1).
        Logger(DEBUG) << "      > Explicit-volume half-step kick B";
        particles->integrateExplicitVolumeHalfStep(timeStep / 2.0, /*finalize=*/true);
#endif

        //Logger(INFO) << "    > Moving particles";
        Logger(INFO) << "    > Moving particles";
        particles->moveParticles(timeStep, domain);
#endif // !ELASTIC

        t += timeStep;
        ++step;

#if ADAPTIVE_TIMESTEP
        if (dumpNext){
            dump = true;
            dumpNext = false;
        }
#endif // ADAPTIVE_TIMESTEP

        //Logger(DEBUG) << "    > t = " << t << ", step =  " << step
        //          << ", t_end = " << config.timeEnd;

        // DEBUGGING
        // TODO: remove
        //Logger(DEBUG) << "      SANITY CHECK > M_tot = " << particles->sumMass();

        //stepss = std::stringstream();;
        //Logger(INFO) << "   > Dump particle distribution";
        //stepss << std::setw(6) << std::setfill('0') << step;
        //Logger(INFO) << "      > Dump particles to file";
        //particles->dump2file(config.outDir + "/" + stepss.str() + std::string(".h5"));
        // END DEBUGGING

    } while(t<config.timeEnd+timeStep);

    double elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - tStart).count() / 60.0;
    Logger(INFO) << "Total wall-clock time: " << elapsed << " min";
}

MeshlessScheme::~MeshlessScheme(){}
