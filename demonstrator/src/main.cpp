//
// Created by Johannes Martin on 17.09.21.
//

// header only libraries
#include <cxxopts.hpp>

#include "../include/Logger.h"
#include "../include/ConfigParser.h"
#include "../include/MeshlessScheme.h"
#include "../include/SPH.h"
#include "../include/EquationOfState.h"


structlog LOGCFG = {};

int main(int argc, char *argv[]){

    cxxopts::Options cmdLineOptions { "mlh",
                                      "Demonstrator for the meshless hydrodynamic simulation methods MFV and MFM. Also does SPH with smoothed gradient" };
    cmdLineOptions.add_options()
            ("c,config", "Path to config file", cxxopts::value<std::string>()->default_value("config.info"))
            ("v,verbose", "More printouts for debugging")
            ("s,silent", "Suppress normal printouts")
            ("h,help", "Show this help");

    auto cmdLineOpts = cmdLineOptions.parse(argc, argv);

    if (cmdLineOpts.count("help")) {
        std::cout << cmdLineOptions.help() << std::endl;
        exit(0);
    }

    ConfigParser confP { cmdLineOpts["config"].as<std::string>() };

    // initialize Logger
    LOGCFG.headers = true;
    LOGCFG.level = cmdLineOpts.count("verbose") ? DEBUG : INFO;

    if (cmdLineOpts.count("silent")){
        if(cmdLineOpts.count("verbose")){
            throw std::invalid_argument("Command line options -s and -v are incompatible");
        } else {
            LOGCFG.level = WARN;
        }
    }


    Logger(INFO) << "Reading configuration ... ";
#if RUNSPH
    SPH::Configuration config;
#else
    MeshlessScheme::Configuration config;
#endif

    config.initFile = confP.getVal<std::string>("initFile");
    Logger(INFO) << "    > Initial distribution: " << config.initFile;
    config.outDir = confP.getVal<std::string>("outDir");
    Logger(INFO) << "    > Output directory: " << config.outDir;
    config.timeStep = confP.getVal<double>("timeStep");
    Logger(INFO) << "    > Time step: " << config.timeStep;
    config.timeEnd = confP.getVal<double>("timeEnd");
    Logger(INFO) << "    > End of simulation: " << config.timeEnd;
    config.h5DumpInterval = confP.getVal<int>("h5DumpInterval");
    Logger(INFO) << "    > Dump data to h5 file every " << config.h5DumpInterval << " steps";
    config.kernelSize = confP.getVal<double>("kernelSize");
    Logger(INFO) << "    > Using global kernel size h = " << config.kernelSize;
#if EOS == 0
    config.hydro_gamma = confP.getVal<double>("hydro_gamma");
    Logger(INFO) << "    > Adiabatic index for ideal gas EOS hydro_gamma = " << config.hydro_gamma;
    EquationOfState MeshlessEOS(config.hydro_gamma);
#elif EOS == 1
    config.K0 = confP.getVal<double>("MURN_K0");
    Logger(INFO) << "    > Bulk modulus for Murnaghan EOS K0 = " << config.K0;
    config.murn_n = confP.getVal<double>("MURN_n");
    Logger(INFO) << "    > Murnaghan exponent for Murnaghan EOS n = " << config.murn_n;
    config.rho0 = confP.getVal<double>("MURN_rho0");
    Logger(INFO) << "    > Relaxed density for Murnaghan EOS rho0 = " << config.rho0;
    EquationOfState MeshlessEOS(config.K0, config.murn_n, config.rho0);
#elif EOS == 2
    config.TIL_A = confP.getVal<double>("TIL_A");
    Logger(INFO) << "   > Tillotson EOS A_T = " << config.TIL_A;
    config.TIL_B = confP.getVal<double>("TIL_B");
    Logger(INFO) << "   > Tillotson EOS B_T = " << config.TIL_B;
    config.TIL_u0 = confP.getVal<double>("TIL_u0");
    Logger(INFO) << "   > Tillotson EOS TIL_u0 = " << config.TIL_u0;
    config.TIL_a = confP.getVal<double>("TIL_a");
    Logger(INFO) << "   > Tillotson EOS TIL_a = " << config.TIL_a;
    config.TIL_b = confP.getVal<double>("TIL_b");
    Logger(INFO) << "   > Tillotson EOS TIL_b = " << config.TIL_b;
    config.TIL_alpha = confP.getVal<double>("TIL_alpha");
    Logger(INFO) << "   > Tillotson EOS TIL_alpha = " << config.TIL_alpha;
    config.TIL_beta = confP.getVal<double>("TIL_beta");
    Logger(INFO) << "   > Tillotson EOS TIL_beta = " << config.TIL_beta;
    config.TIL_u_iv = confP.getVal<double>("TIL_u_iv");
    Logger(INFO) << "   > Energy of incipient vaporization Tillotson EOS TIL_u_iv = " << config.TIL_u_iv;
    config.TIL_u_cv = confP.getVal<double>("TIL_u_cv");
    Logger(INFO) << "   > Energy of complete vaporization Tillotson EOS TIL_u_cv = " << config.TIL_u_cv;
    EquationOfState MeshlessEOS(config.TIL_A, config.TIL_B, config.TIL_u0, config.TIL_a,
                    config.TIL_b, config.TIL_alpha, config.TIL_beta, config.TIL_u_iv, config.TIL_u_cv);
#endif // EOS
    Logger(INFO) << "    > Initialized Equation Of State";
#if PERIODIC_BOUNDARIES
    auto periodicBoxLimits = confP.getObj("periodicBoxLimits");
    config.periodicBoxLimits[0] = periodicBoxLimits.getVal<double>("lowerX");
    config.periodicBoxLimits[DIM] = periodicBoxLimits.getVal<double>("upperX");
    config.periodicBoxLimits[1] = periodicBoxLimits.getVal<double>("lowerY");
    config.periodicBoxLimits[DIM+1] = periodicBoxLimits.getVal<double>("upperY");
#if DIM == 3
    config.periodicBoxLimits[2] = periodicBoxLimits.getVal<double>("lowerZ");
    config.periodicBoxLimits[DIM+2] = periodicBoxLimits.getVal<double>("upperZ");
#endif
    std::string periodicBoxStr = "[";
    for (int i=0; i<2*DIM; i++){
        periodicBoxStr.append(std::to_string(config.periodicBoxLimits[i]));
        if(i<2*DIM-1) periodicBoxStr.append(", ");
    }
    Logger(INFO) << "    > Periodic boundaries within box: " << periodicBoxStr << "]";
#endif

    Logger(INFO) << "    > Reading initial distribution ...";

    InitialDistribution initDist { config.initFile };
    Particles particles { initDist.getNumberOfParticles(), &MeshlessEOS};
    initDist.getAllParticles(particles);

    Logger(INFO) << "    > N = " << particles.N;
    Logger(INFO) << "... done. Initializing simulation ...";

#if PERIODIC_BOUNDARIES
    double *domainLimits = config.periodicBoxLimits;
#else
    double domainLimits[DIM*2];
    particles.getDomainLimits(domainLimits);
#endif
    Domain::Cell boundingBox { domainLimits };


#if RUNSPH
    SPH algorithm {config, &particles, boundingBox};
#else // RUNSPH
    MeshlessScheme algorithm { config, &particles, boundingBox };
#endif // RUNSPH

    Logger(INFO) << "... done.";

    Logger(INFO) << "Starting time integration ...";
    algorithm.run();
    Logger(INFO) << "... done.";
    return 0;
}
