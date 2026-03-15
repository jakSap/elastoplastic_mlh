//
// Created by Johannes Martin on 28.09.22.
//

#include "../include/Riemann.h"
//#include "../include/riemannhelper.h"
Riemann::Riemann(double *WR, double *WL, double *vFrame, double *Aij, int i, EquationOfState &MeshlessEOS
#if USE_HLLC && ELASTIC
                 , double SxxR, double SxyR, double SyyR
                 , double SxxL, double SxyL, double SyyL
#if DIM == 3
                 , double SxzR, double SyzR, double SzzR
                 , double SxzL, double SyzL, double SzzL
#endif
#endif
                 ) :
                    WR { WR }, WL { WL }, vFrame { vFrame },  Aij { Aij }, i { i }, MeshlessEOS{MeshlessEOS}{

    // compute norm of effective face
    // double AijNorm = sqrt(Helper::dotProduct(Aij, Aij));

    // compute norm of effective face
    AijNorm = sqrt(Helper::dotProduct(Aij, Aij));
    
    // Logger(INFO) << "Gamma by EOS: " << MeshlessEOS.EOSPressure(WR[0], 0);

    if(AijNorm > 0){
        hatAij[0] = 1./AijNorm*Aij[0];
        hatAij[1] = 1./AijNorm*Aij[1];
#if DIM == 3
        hatAij[2] = 1./AijNorm*Aij[2];
#endif

    }
    else{
        Logger(WARN) << "Warning, AijNorm = 0 encountered, i = " << i;
        Logger(WARN) << "Aborting for debugging";
        // exit(6);
	
        hatAij[0] = 0;
        hatAij[1] = 0;
#if DIM == 3
        hatAij[2] = 0;
#endif
    }
    // if (i == 6){
    if (i == 12157){
       Logger(DEBUG) << "Aij = [" << Aij[0] << ", " << Aij[1] << "], AijNorm = " << AijNorm;
    }

#if USE_HLLC && ELASTIC
    //    Logger(DEBUG) << "Aij = [" << Aij[0] << ", " << Aij[1] << "], AijNorm = " << AijNorm;
    double SijR[DIM*DIM], SijL[DIM*DIM];
#if DIM == 2
    SijR[0] = SxxR; SijR[1] = SxyR;
    SijR[2] = SxyR; SijR[3] = SyyR;
    SijL[0] = SxxL; SijL[1] = SxyL;
    SijL[2] = SxyL; SijL[3] = SyyL;
#else
    SijR[0] = SxxR; SijR[1] = SxyR; SijR[2] = SxzR;
    SijR[3] = SxyR; SijR[4] = SyyR; SijR[5] = SyzR;
    SijR[6] = SxzR; SijR[7] = SyzR; SijR[8] = SzzR;
    SijL[0] = SxxL; SijL[1] = SxyL; SijL[2] = SxzL;
    SijL[3] = SxyL; SijL[4] = SyyL; SijL[5] = SyzL;
    SijL[6] = SxzL; SijL[7] = SyzL; SijL[8] = SzzL;
#endif
    // }
    double Lambda[DIM*DIM];
#if DIM == 2
    Helper::rotationMatrix2D(hatAij, unitX, Lambda);
#else
    Helper::rotationMatrix3D(hatAij, unitX, Lambda);
#endif

    // Compute Lambda^T
    double LambdaT[DIM*DIM];
    for (int a = 0; a < DIM; ++a)
        for (int b = 0; b < DIM; ++b)
            LambdaT[b + DIM*a] = Lambda[a + DIM*b];

    // S' = Lambda * S * Lambda^T
    double temp[DIM*DIM];
    Helper::matMul(Lambda, SijR, temp);
    Helper::matMul(temp, LambdaT, SijRotR);
    Helper::matMul(Lambda, SijL, temp);
    Helper::matMul(temp, LambdaT, SijRotL);

    // Store Lambda for flux back-rotation in HLLCFlux
    for (int a = 0; a < DIM*DIM; ++a) Lambda_[a] = Lambda[a];

    // Rotate velocities (xSplitElasticHLLC expects vx = normal velocity)
#if DIM == 2
    double vBufR[DIM] = { WR[2], WR[3] };
    double vBufL[DIM] = { WL[2], WL[3] };
    WR[2] = Lambda[0]*vBufR[0] + Lambda[1]*vBufR[1];
    WR[3] = Lambda[2]*vBufR[0] + Lambda[3]*vBufR[1];
    WL[2] = Lambda[0]*vBufL[0] + Lambda[1]*vBufL[1];
    WL[3] = Lambda[2]*vBufL[0] + Lambda[3]*vBufL[1];
#else
    double vBufR[DIM] = { WR[2], WR[3], WR[4] };
    double vBufL[DIM] = { WL[2], WL[3], WL[4] };
    WR[2] = Lambda[0]*vBufR[0] + Lambda[1]*vBufR[1] + Lambda[2]*vBufR[2];
    WR[3] = Lambda[3]*vBufR[0] + Lambda[4]*vBufR[1] + Lambda[5]*vBufR[2];
    WR[4] = Lambda[6]*vBufR[0] + Lambda[7]*vBufR[1] + Lambda[8]*vBufR[2];
    WL[2] = Lambda[0]*vBufL[0] + Lambda[1]*vBufL[1] + Lambda[2]*vBufL[2];
    WL[3] = Lambda[3]*vBufL[0] + Lambda[4]*vBufL[1] + Lambda[5]*vBufL[2];
    WL[4] = Lambda[6]*vBufL[0] + Lambda[7]*vBufL[1] + Lambda[8]*vBufL[2];
#endif

    // Rotate vFrame for consistent de-boosting in x-split solver
#if DIM == 2
    double vFrameBuf[DIM] = { vFrame[0], vFrame[1] };
    vFrame[0] = Lambda[0]*vFrameBuf[0] + Lambda[1]*vFrameBuf[1];
    vFrame[1] = Lambda[2]*vFrameBuf[0] + Lambda[3]*vFrameBuf[1];
#else
    double vFrameBuf[DIM] = { vFrame[0], vFrame[1], vFrame[2] };
    vFrame[0] = Lambda[0]*vFrameBuf[0] + Lambda[1]*vFrameBuf[1] + Lambda[2]*vFrameBuf[2];
    vFrame[1] = Lambda[3]*vFrameBuf[0] + Lambda[4]*vFrameBuf[1] + Lambda[5]*vFrameBuf[2];
    vFrame[2] = Lambda[6]*vFrameBuf[0] + Lambda[7]*vFrameBuf[1] + Lambda[8]*vFrameBuf[2];
#endif
#endif // USE_HLLC && ELASTIC

#if !USE_HLLC
//#if DIM == 3
//    Logger(ERROR) << "Rotation of states is not yet implemented in 3D. - Aborting.";
//    exit(4);
//#endif
    // rotate states to be aligned with hatAij
    double Lambda[DIM*DIM];
    //Helper::rotationMatrix2D(unitX, hatAij, Lambda);

#if DIM == 2
    Helper::rotationMatrix2D(hatAij, unitX, Lambda);

    double vBufR[DIM] = { WR[2], WR[3] };
    double vBufL[DIM] = { WL[2], WL[3] };

    // Implementing eq. A5 in GIZMO
    WR[2] = Lambda[0]*vBufR[0]+Lambda[1]*vBufR[1];
    WR[3] = Lambda[2]*vBufR[0]+Lambda[3]*vBufR[1];
    WL[2] = Lambda[0]*vBufL[0]+Lambda[1]*vBufL[1];
    WL[3] = Lambda[2]*vBufL[0]+Lambda[3]*vBufL[1];

    //Logger(DEBUG) << "hatAij = [" << hatAij[0] << ", " << hatAij[1] << "], |hatAij| = "
    //          << sqrt(Helper::dotProduct(hatAij, hatAij));
    //Logger(DEBUG) << "Lambda = [" << Lambda[0] << ", " << Lambda[1];
    //Logger(DEBUG) << "          " << Lambda[2] << ", " << Lambda[3] << "]";

    //Logger(DEBUG) << "vR_unrot = [" << WR[2] << ", " << WR[3] << "], vL_unrot = ["
    //          <<  WL[2] << ", " << WL[3] << "]";

    //if(i == 204){
    //    Logger(DEBUG) <<  "rhoR = " << WR[0] << ", rhoL = " << WL[0]
    //                  << ", vR = [" << WR[2] << ", " << WR[3]
    //                  << "], vL = [" <<  WL[2] << ", " << WL[3]
    //                  << "], PR = " << WR[1] << ", PL = " << WL[1];
    //}

    //if(i == 46){
    //    Logger(DEBUG) <<  "rhoR = " << WR[0] << ", rhoL = " << WL[0]
    //                  << ", vR = [" << WR[2] << ", " << WR[3]
    //                  << "], vL = [" <<  WL[2] << ", " << WL[3]
    //                  << "], PR = " << WR[1] << ", PL = " << WL[1];
    //}
#else
    Helper::rotationMatrix3D(hatAij, unitX, Lambda);

    //Logger(DEBUG) << "Lambda = [" << Lambda[0] << ", " << Lambda[1] << ", " << Lambda[2];
    //Logger(DEBUG) << "          " << Lambda[3] << ", " << Lambda[4] << ", " << Lambda[5];
    //Logger(DEBUG) << "          " << Lambda[6] << ", " << Lambda[7] << ", " << Lambda[8] << "]";

    double vBufR[DIM] = { WR[2], WR[3], WR[4] };
    double vBufL[DIM] = { WL[2], WL[3], WL[4] };

    WR[2] = Lambda[0]*vBufR[0]+Lambda[1]*vBufR[1]+Lambda[2]*vBufR[2];
    WR[3] = Lambda[3]*vBufR[0]+Lambda[4]*vBufR[1]+Lambda[5]*vBufR[2];
    WR[4] = Lambda[6]*vBufR[0]+Lambda[7]*vBufR[1]+Lambda[8]*vBufR[2];

    WL[2] = Lambda[0]*vBufL[0]+Lambda[1]*vBufL[1]+Lambda[2]*vBufL[2];
    WL[3] = Lambda[3]*vBufL[0]+Lambda[4]*vBufL[1]+Lambda[5]*vBufL[2];
    WL[4] = Lambda[6]*vBufL[0]+Lambda[7]*vBufL[1]+Lambda[8]*vBufL[2];
#endif
#endif // !USE_HLLC
}

#if USE_HLLC
// Add HLLC function for approximate Riemann Solver
void Riemann::HLLCFlux(double *Fij){
#if ELASTIC
    // x-split elastic HLLC: states and vFrame are already rotated in constructor
#if DEBUG_LVL
    if(i == 2){
        Logger(DEBUG) << "HLLCFlux i=" << i << ":"
                    << " WL=[rho=" << WL[0] << ", P=" << WL[1]
                    << ", vx=" << WL[2] << ", vy=" << WL[3]
#if DIM == 3
                    << ", vz=" << WL[4]
#endif
                    << "]";
        Logger(DEBUG) << "  i=" << i
                    << " WR=[rho=" << WR[0] << ", P=" << WR[1]
                    << ", vx=" << WR[2] << ", vy=" << WR[3]
#if DIM == 3
                    << ", vz=" << WR[4]
#endif
                    << "]";
        Logger(DEBUG) << "  i=" << i
                    << " SijRotL=[Sxx=" << SijRotL[0] << ", Sxy=" << SijRotL[1] << "]"
                    << " SijRotR=[Sxx=" << SijRotR[0] << ", Sxy=" << SijRotR[1] << "]";
    }
#endif
    HLLC::xSplitElasticHLLC(WR, WL, SijRotR, SijRotL, Fij, MeshlessEOS);

    // Back-rotate momentum fluxes to lab frame (Lambda^{-1} = Lambda^T)
#if DIM == 2
    double FijBuf[DIM] = { Fij[2], Fij[3] };
    Fij[2] = Lambda_[0]*FijBuf[0] + Lambda_[2]*FijBuf[1];
    Fij[3] = Lambda_[1]*FijBuf[0] + Lambda_[3]*FijBuf[1];
#else
    double FijBuf[DIM] = { Fij[2], Fij[3], Fij[4] };
    Fij[2] = Lambda_[0]*FijBuf[0] + Lambda_[3]*FijBuf[1] + Lambda_[6]*FijBuf[2];
    Fij[3] = Lambda_[1]*FijBuf[0] + Lambda_[4]*FijBuf[1] + Lambda_[7]*FijBuf[2];
    Fij[4] = Lambda_[2]*FijBuf[0] + Lambda_[5]*FijBuf[1] + Lambda_[8]*FijBuf[2];
#endif

    // De-boost to lab frame: back-rotate vFrame via Lambda^T, then Galilean correction
#if DIM == 2
    double vLab[DIM];
    vLab[0] = Lambda_[0]*vFrame[0] + Lambda_[2]*vFrame[1];
    vLab[1] = Lambda_[1]*vFrame[0] + Lambda_[3]*vFrame[1];
    const double v2 = vLab[0]*vLab[0] + vLab[1]*vLab[1];
    Fij[1] += vLab[0]*Fij[2] + vLab[1]*Fij[3] + 0.5 * v2 * Fij[0];
    Fij[2] += vLab[0] * Fij[0];
    Fij[3] += vLab[1] * Fij[0];
#else
    double vLab[DIM];
    vLab[0] = Lambda_[0]*vFrame[0] + Lambda_[3]*vFrame[1] + Lambda_[6]*vFrame[2];
    vLab[1] = Lambda_[1]*vFrame[0] + Lambda_[4]*vFrame[1] + Lambda_[7]*vFrame[2];
    vLab[2] = Lambda_[2]*vFrame[0] + Lambda_[5]*vFrame[1] + Lambda_[8]*vFrame[2];
    const double v2 = vLab[0]*vLab[0] + vLab[1]*vLab[1] + vLab[2]*vLab[2];
    Fij[1] += vLab[0]*Fij[2] + vLab[1]*Fij[3] + vLab[2]*Fij[4] + 0.5 * v2 * Fij[0];
    Fij[2] += vLab[0] * Fij[0];
    Fij[3] += vLab[1] * Fij[0];
    Fij[4] += vLab[2] * Fij[0];
#endif

    // Project onto effective faces
    Fij[0] *= AijNorm;
    Fij[1] *= AijNorm;
    Fij[2] *= AijNorm;
    Fij[3] *= AijNorm;
#if DIM == 3
    Fij[4] *= AijNorm;
#endif

#else // !ELASTIC
#if DIM==3
        HLLC::solveHLLC(WR, WL, hatAij, Fij, vFrame, gamma); // ToDo: 3d is not implemented

        // Rotate and project fluxes onto Aij
#else
        // Logger(DEBUG) << " WL: " << WL[0] << " " << WL[1] << " " << WL[2] << " " << WL[3];

#if USE_HLL
        const double gamma = MeshlessEOS.EOSGetHydroGammaParam();
        HLLC::HLL(WR, WL, Fij, gamma);
#else
        HLLC::solveHLLC1(WR, WL, hatAij, Fij, vFrame, MeshlessEOS);
#endif  // USE_HLL

        // Logger(DEBUG) << "i = " << i << " mF = " << Fij[0];
        // Rotate fluxes back into the Aij-direction:
        // double LambdaInv[DIM*DIM];

        // // Build rotation matrix to rotate from unitX to hatAij
        // Helper::rotationMatrix2D(unitX, hatAij, LambdaInv);
        // // Helper::rotationMatrix2D(unitX, hatAij, LambdaInv);

        // // Buffer velocity fluxes
        // double FijBuf[DIM] = { Fij[2], Fij[3] };

        // // Rotate fluxes: Eq A7
        // Fij[2] = LambdaInv[0]*FijBuf[0]+LambdaInv[1]*FijBuf[1];
        // Fij[3] = LambdaInv[2]*FijBuf[0]+LambdaInv[3]*FijBuf[1];

        double vFrame2 = pow(vFrame[0], 2) + pow(vFrame[1], 2);


        /// Compute fluxes projected onto Aij
        // Fij[0] *= (Aij[0] + Aij[1]); // mass flux


        // Project onto effective faces:
        Fij[0] *= AijNorm;
        Fij[1] *= AijNorm;
        //
        // Option 1: Hadamarf product
        // double v = sqrtf(pow(Fij[2], 2) + pow(Fij[3], 2));
        Fij[2] = Fij[2] * AijNorm;
        Fij[3] = Fij[3] * AijNorm;

#endif // DIM
#endif // !ELASTIC
}
#endif // USE_HLLC
// Exact Riemann solver
void Riemann::exact(double *Fij){
#if HLLC_general_EOS
    Logger(INFO) << "For general EOS, only HLLC solver is available, exiting.";
    exit(6);
#endif
#if EOS == 0
    const double gamma = MeshlessEOS.EOSGetHydroGammaParam();

    RiemannSolver solver { gamma };

    //if(WR[1] < 0. || WL[1] < 0.){
    //    Logger(ERROR) << "Negative pressure encountered. Very bad :( !!";
    //    Logger(DEBUG) << "rhoL = " << WL[0] << ", rhoR = " << WR[0]
    //                  << ", uL = " << WL[2] << ", uR = " << WR[2]
    //                  << ", PL = " << WL[1] << ", PR = " << WR[1];
    //}

    int flagLR = solver.solve(WL[0], WL[2], WL[1], WR[0], WR[2], WR[1],
                              rhoSol, vSol[0], PSol);

    // Fij[3] = 0.; // explicitly setting vy flux to zero
    //if (Fij[2] < 0.){
    //    Logger(DEBUG) << "Fij = [" << Fij[0] << " (mass), " << Fij[2] << " (vx), "
    //              << Fij[3] << " (vy), " << Fij[1] << " (energy)]";
    //}

    /// CONVERT SOLUTION TO FLUXES
#if DIM==3
    if (flagLR == 1){
        // right state sampled
        vSol[1] = WR[3];
        vSol[2] = WR[4];
    } else if (flagLR == -1){ // flagLR == -1
        // left state sampled
        vSol[1] = WL[3];
        vSol[2] = WL[4];
    } else { // flagLR == 0
        Logger(WARN) << "  > Vacuum state sampled. This is not expected.";
    }

    rotateAndProjectFluxes3D(Fij, gamma);

#else // DIM == 3
    if (flagLR == 1){
        // right state sampled
        vSol[1] = WR[3];
    } else if (flagLR == -1){ // flagLR == -1
        // left state sampled
        vSol[1] = WL[3];
    } else { // flagLR == 0
        Logger(WARN) << "  > Vacuum state sampled. This is not expected.";
    }

    //if (i == 46){
    //    Logger(DEBUG) << "riemann solution = [" << rhoSol << ", " << vSol[0] << ", " << vSol[1] << ", " << PSol << "]";
    //}

    rotateAndProjectFluxes2D(Fij);

    //if (i == 46){
    //    Logger(DEBUG) << "Fij = [" << Fij[0] << " (mass), " << Fij[2] << " (vx), "
    //                  << Fij[3] << " (vy), " << Fij[1] << " (energy)]";
    //}

#endif // DIM == 3
#else // EOS == 0
    Logger(INFO) << "For EOS other than the ideal gas EOS, only HLLC solver is available, exiting.";
    exit(6);
#endif // EOS == 0
}

#if DIM==2
void Riemann::rotateAndProjectFluxes2D(double *Fij){
#if HLLC_general_EOS
    Logger(INFO) << "Exact solver with general EOS not implemented yet, exiting";
    exit(6);
#endif
#if EOS == 0
    const double gamma = MeshlessEOS.EOSGetHydroGammaParam();
#else
    const double gamma = -1; // Placeholder
#endif

    // rotate back velocities to simulation coordinate frame
    double LambdaInv[DIM*DIM];
    //Helper::rotationMatrix2D(hatAij, unitX, LambdaInv);
    Helper::rotationMatrix2D(unitX, hatAij, LambdaInv);

    //Logger(DEBUG) << "vSol_unrot = [" << vSol[0] << ", " << vSol[1] << "]";
    //Logger(DEBUG) << "LambdaInv = [" << LambdaInv[0] << ", " << LambdaInv[1]
    //          << ", " << LambdaInv[2] << ", " << LambdaInv[3] << "]";

    double vSolBuf[DIM] = { vSol[0], vSol[1] };
    vSol[0] = LambdaInv[0]*vSolBuf[0]+LambdaInv[1]*vSolBuf[1];
    vSol[1] = LambdaInv[2]*vSolBuf[0]+LambdaInv[3]*vSolBuf[1];

    //if (i == 204){
    //    Logger(DEBUG) << "rhoSol = " << rhoSol << ", vSol_rot = [" << vSol[0] << ", " << vSol[1] << "], PSol = " << PSol;
    //}

    /// Compute fluxes projected onto Aij
    Fij[0] = Aij[0]*rhoSol*vSol[0] + Aij[1]*rhoSol*vSol[1]; // mass flux



    // De-boost: Velocity solution in the lab frame
    double vLab[DIM];
    vLab[0] = vSol[0] + vFrame[0];
    vLab[1] = vSol[1] + vFrame[1];

#if MESHLESS_FINITE_MASS
    // this has no direct physical meaning but effectively suppresses all advection terms
    vSol[0] = 0.;
    vSol[1] = 0.;
#endif // MESHLESS_FINITE_MASS

        
    Fij[2] = Aij[0]*(rhoSol*vLab[0]*vSol[0]+PSol) + Aij[1]*rhoSol*vLab[0]*vSol[1]; // vx flux
    Fij[3] = Aij[0]*rhoSol*vLab[1]*vSol[0] + Aij[1]*(rhoSol*vLab[1]*vSol[1]+PSol); // vy flux

    Fij[1] = Aij[0]*(vSol[0]*(PSol/(gamma-1.)+rhoSol*.5*Helper::dotProduct(vLab, vLab)) + PSol*vLab[0])
            + Aij[1]*(vSol[1]*(PSol/(gamma-1.)+rhoSol*.5*Helper::dotProduct(vLab, vLab)) + PSol*vLab[1]); // energy flux

    //double vFluxBuf[DIM] = { Fij[2], Fij[3] };
    //Fij[2] = LambdaInv[0]*vFluxBuf[0]+LambdaInv[1]*vFluxBuf[1];
    //Fij[3] = LambdaInv[2]*vFluxBuf[0]+LambdaInv[3]*vFluxBuf[1];
}

#else
void Riemann::rotateAndProjectFluxes3D(double *Fij){
#if HLLC_general_EOS
    Logger(Info << "Exact solver with general EOS not implemented yet, exiting");
    exit(6);
#endif
#if EOS == 0
    const double gamma = MeshlessEOS.EOSGetHydroGammaParam();
#else
    const double gamma = -1; // PLaceholder
#endif

    // rotate back velocities to simulation coordinate frame
    double LambdaInv[DIM*DIM];
    Helper::rotationMatrix3D(unitX, hatAij, LambdaInv);

    double vSolBuf[DIM] = { vSol[0], vSol[1], vSol[2] };
    vSol[0] = LambdaInv[0]*vSolBuf[0]+LambdaInv[1]*vSolBuf[1]+LambdaInv[2]*vSolBuf[2];
    vSol[1] = LambdaInv[3]*vSolBuf[0]+LambdaInv[4]*vSolBuf[1]+LambdaInv[5]*vSolBuf[2];
    vSol[2] = LambdaInv[6]*vSolBuf[0]+LambdaInv[7]*vSolBuf[1]+LambdaInv[8]*vSolBuf[2];

    //if (i == 204){
    //    Logger(DEBUG) << "rhoSol = " << rhoSol << ", vSol_rot = [" << vSol[0] << ", " << vSol[1] << "], PSol = " << PSol;
    //}

    /// Compute fluxes projected onto Aij
    Fij[0] = Aij[0]*rhoSol*vSol[0] + Aij[1]*rhoSol*vSol[1] + Aij[2]*rhoSol*vSol[2]; // mass flux

    // velocity solution in the lab frame
    double vLab[DIM];
    vLab[0] = vSol[0] + vFrame[0];
    vLab[1] = vSol[1] + vFrame[1];
    vLab[2] = vSol[2] + vFrame[2];

#if MESHLESS_FINITE_MASS
    // this has no direct physical meaning but effectively suppresses all advection terms
    vSol[0] = 0.;
    vSol[1] = 0.;
    vSol[2] = 0.;
#endif // MESHLESS_FINITE_MASS

    Fij[2] = Aij[0]*(rhoSol*vLab[0]*vSol[0]+PSol) + Aij[1]*rhoSol*vLab[0]*vSol[1] + Aij[2]*rhoSol*vLab[0]*vSol[2]; // vx flux
    Fij[3] = Aij[0]*rhoSol*vLab[1]*vSol[0] + Aij[1]*(rhoSol*vLab[1]*vSol[1]+PSol) + Aij[2]*rhoSol*vLab[1]*vSol[2]; // vy flux
    Fij[4] = Aij[0]*rhoSol*vLab[2]*vSol[0] + Aij[1]*rhoSol*vLab[2]*vSol[1] + Aij[2]*(rhoSol*vLab[2]*vSol[2]+PSol); // vz flux

    Fij[1] = Aij[0]*(vSol[0]*(PSol/(gamma-1.)+rhoSol*.5*Helper::dotProduct(vLab, vLab)) + PSol*vLab[0])
             + Aij[1]*(vSol[1]*(PSol/(gamma-1.)+rhoSol*.5*Helper::dotProduct(vLab, vLab)) + PSol*vLab[1])
             + Aij[2]*(vSol[2]*(PSol/(gamma-1.)+rhoSol*.5*Helper::dotProduct(vLab, vLab)) + PSol*vLab[2]); // energy flux

}
#endif
