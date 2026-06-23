//
// Created by Johannes Martin on 05.09.22.
//

#include "../include/Particles.h"

// Kernel functions live in include/Kernel.h (selected via KERNEL_FUNCTION).

Particles::Particles(int numParticles, EquationOfState *MeshlessEOS
            , bool ghosts
            ) : N { numParticles }, MeshlessEOS(MeshlessEOS),
            ghosts { ghosts }
{
    // allocate memory
#if USE_MATID
    matId = new int[numParticles];
#endif
    cell = new int[numParticles];
    m = new double[numParticles];
    u = new double[numParticles];
    rho = new double[numParticles];
    P = new double[numParticles];
    sml = new double[numParticles];
    x = new double[numParticles];
    y = new double[numParticles];
    vx = new double[numParticles];
    vy = new double[numParticles];
    rhoGrad = new double[numParticles][DIM];
    vxGrad = new double[numParticles][DIM];
    vyGrad = new double[numParticles][DIM];
    PGrad = new double[numParticles][DIM];

#if ELASTIC
    // value-initialize (zero) so stress starts at the Murnaghan equilibrium
    // for u = 0. Uninitialized memory here breaks stationary ICs where no
    // velocity gradient exists to wash out the garbage.
#if TENSILE_CORRECTION
    fabMonaghan = new double[numParticles*MAX_NUM_INTERACTIONS]();
#endif
    Sxx = new double[numParticles]();
    Sxy = new double[numParticles]();
    Syy = new double[numParticles]();

    SxxGrad = new double[numParticles][DIM]();
    SxyGrad = new double[numParticles][DIM]();
    SyyGrad = new double[numParticles][DIM]();

#if DIM == 3
    Sxz = new double[numParticles]();
    Syz = new double[numParticles]();
    Szz = new double[numParticles]();

    SxzGrad = new double[numParticles][DIM]();
    SyzGrad = new double[numParticles][DIM]();
    SzzGrad = new double[numParticles][DIM]();
#endif // DIM == 3
#if FRAGMENTATION
    damage         = new double[numParticles]();
    dddt           = new double[numParticles]();
    damageTotal    = new double[numParticles]();
    flaws          = new double[numParticles*MAX_NUM_FLAWS]();
    numFlaws       = new int[numParticles]();
    numActiveFlaws = new int[numParticles]();
#endif // FRAGMENTATION
#endif // ELASTIC

#if RUNSPH
    ax = new double[numParticles];
    ay = new double[numParticles];

    // ArtVisc thingy
    cs = new double[numParticles];
    dudtArtVisc = new double[numParticles];

    axArtVisc = new double[numParticles];
    ayArtVisc = new double[numParticles];
#if DIM == 3
    az = new double[numParticles];
    azArtVisc = new double[numParticles];
#endif


	//unimportant = new double[numParticles];
	dEdt = new double[numParticles];
	dn = new double[numParticles];
	drho = new double[numParticles];
#endif // RUNSPH

    //TODO: check if this is needed as array
    //B = new double[numParticles][DIM*DIM];
#if DIM == 3
    z = new double[numParticles];
    vz = new double[numParticles];
    vzGrad = new double[numParticles][DIM];
#endif
    omega = new double[numParticles];
#if SURFACE_VOLCORR
    fce = new double[numParticles];
    for (int i=0; i<numParticles; ++i) fce[i] = 1.;
#endif
#if EXPLICIT_VOL_INTEGRATION
    rhoExplicit = new double[numParticles]();
    rhoKernel   = new double[numParticles]();
#endif
#if OUTPUT_CONDITION_NUMBER
    conditionNumber = new double[numParticles];
#if DIM == 2
    lambdaMax = new double[numParticles];
    lambdaMin = new double[numParticles];
    eigenvecMin = new double[numParticles][DIM];
#endif
#endif
    if (!ghosts){
        nnl = new int[numParticles*MAX_NUM_INTERACTIONS];
        noi = new int[numParticles];
        psijTilde_xi = new double[numParticles*MAX_NUM_INTERACTIONS][DIM];
        Aij = new double[numParticles*MAX_NUM_INTERACTIONS][DIM];
        WijL = new double[numParticles*MAX_NUM_INTERACTIONS][DIM+2];
        WijR = new double[numParticles*MAX_NUM_INTERACTIONS][DIM+2];
        Fij = new double[numParticles*MAX_NUM_INTERACTIONS][DIM+2]; // TODO: this buffer should not be needed
        vFrame = new double[numParticles*MAX_NUM_INTERACTIONS][DIM];

        mF = new double[numParticles];
        vF = new double[numParticles][DIM];
        eF = new double[numParticles];

#if PERIODIC_BOUNDARIES
        // estimated memory allocation
        nnlGhosts = new int[numParticles*MAX_NUM_GHOST_INTERACTIONS];
        noiGhosts = new int[numParticles];
        Logger(DEBUG) << "Declared array size: " << numParticles*(DIM+1);
        ghostMap = new int[numParticles*(DIM+1)]; // TODO: this is only applicable for DIM==2
        psijTilde_xiGhosts = new double[numParticles*MAX_NUM_GHOST_INTERACTIONS][DIM];
        AijGhosts = new double[numParticles*MAX_NUM_GHOST_INTERACTIONS][DIM];
        WijLGhosts = new double[numParticles*MAX_NUM_GHOST_INTERACTIONS][DIM+2];
        WijRGhosts = new double[numParticles*MAX_NUM_GHOST_INTERACTIONS][DIM+2];
        FijGhosts = new double[numParticles*MAX_NUM_GHOST_INTERACTIONS][DIM+2]; // TODO: this buffer should not be needed
        vFrameGhosts = new double[numParticles*MAX_NUM_GHOST_INTERACTIONS][DIM];
    } else {
        parent = new int[numParticles]; // store index of original node
#endif
    }
#if DEBUG_LVL > 1
    if (ghosts){
        // This is necessary for dumping ghosts to file
        noi = new int[numParticles];
    }
#endif
}

Particles::~Particles() {
    // delete &MeshlessEOS;
#if USE_MATID
    delete[] matId;
#endif
    delete[] cell;
    delete[] m;
    delete[] u;
    delete[] P;
    delete[] rho;
    delete[] sml;
    delete[] x;
    delete[] y;
    delete[] vx;
    delete[] vy;
    delete[] rhoGrad;
    delete[] vxGrad;
    delete[] vyGrad;
    delete[] PGrad;
    delete[] omega;
#if SURFACE_VOLCORR
    delete[] fce;
#endif
#if EXPLICIT_VOL_INTEGRATION
    delete[] rhoExplicit;
    delete[] rhoKernel;
#endif
#if OUTPUT_CONDITION_NUMBER
    delete[] conditionNumber;
#if DIM == 2
    delete[] lambdaMax;
    delete[] lambdaMin;
    delete[] eigenvecMin;
#endif
#endif

#if ELASTIC
#if TENSILE_CORRECTION
    delete[] fabMonaghan;
#endif
    delete[] Sxx;
    delete[] Sxy;
    delete[] Syy;

    delete[] SxxGrad;
    delete[] SxyGrad;
    delete[] SyyGrad;
#if DIM == 3
    delete[] Sxz;
    delete[] Syz;
    delete[] Szz;

    delete[] SxzGrad;
    delete[] SyzGrad;
    delete[] SzzGrad;
#endif // DIM == 3
#if FRAGMENTATION
    delete[] damage;
    delete[] dddt;
    delete[] damageTotal;
    delete[] flaws;
    delete[] numFlaws;
    delete[] numActiveFlaws;
#endif // FRAGMENTATION
#endif // ELASTIC
#if RUNSPH
	//delete[] unimportant;
	delete[] dEdt;
	delete[] dn;
	delete[] drho;

    delete[] ax;
    delete[] ay;

    delete[] cs;
    delete[] dudtArtVisc;

    delete[] axArtVisc;
    delete[] ayArtVisc;
#if DIM == 3
    delete[] az;
    delete[] azArtVisc;
#endif
#endif // RUNSPH

    if (!ghosts) {
        delete[] psijTilde_xi;
        delete[] WijL;
        delete[] WijR;
        delete[] Fij;
        delete[] vFrame;
        delete[] mF;
        delete[] vF;
        delete[] eF;
#if DIM == 3
        delete[] z;
        delete[] vz;
        delete[] vzGrad;
#endif
        delete[] nnl;
        delete[] noi;
        delete[] Aij;
#if PERIODIC_BOUNDARIES
        delete[] nnlGhosts;
        delete[] noiGhosts;
        delete[] psijTilde_xiGhosts;
        delete[] AijGhosts;
        delete[] WijLGhosts;
        delete[] WijRGhosts;
        delete[] ghostMap;
        delete[] FijGhosts;
        delete[] vFrameGhosts;
    } else {
        delete[] parent;
#endif
    }
#if DEBUG_LVL > 1
    if (ghosts){
        // This is necessary for dumping ghosts to file
        delete[] noi;
    }
#endif
}

#if ELASTIC
#if TENSILE_CORRECTION
void Particles::computeFabMonaghan(){
    // Get normalized mean particle distance within kernel length
    // Eqivalent to \Delta p in Monaghan paper
    double deltaP = sqrt(3.14159 / smlNNNTarget);
    // And take the Kernel (follows KERNEL_FUNCTION, as GIZMO's kernel_main does)
    double kernelDeltaP = Kernel::W(deltaP, 1.);

    // Now loop through the interaction pairs:
    for(int i = 0; i < N; ++i){
        for(int jn = 0; jn < noi[i]; ++jn){
            // This way of looping means recomputing...
            int j = nnl[i*MAX_NUM_INTERACTIONS+jn];
            // Now work with i and j.
            double r2 = pow(x[j] - x[i], 2) + pow(y[j] - y[i], 2);
#if DIM == 3
            r2 += pow(z[j] - z[i], 2);
#endif
            double r = sqrt(r2);
            // q = r_ab / max(h_a, h_b) following Monaghan and GIZMO
            double q = r / std::max(sml[i], sml[j]);
            // Take the kernel, again..
            double kernelQ = Kernel::W(q, 1.);
            // To be consistent with Monaghan, dont apply the power of n here.
            fabMonaghan[i*MAX_NUM_INTERACTIONS+jn] = kernelQ / kernelDeltaP;
        }
    }
}
#endif // TENSILE_CORRECTION

#if GIZMO_ELASTIC_FLUX
// Port of GIZMO solids/elastic_stress_tensor_force.h (eigenvalue branch), in the
// demonstrator's lab frame. Computed once per pair (i<jj); the flux-symmetry copy
// negates it for jj, matching GIZMO's antisymmetric +i/-j application.
void Particles::addGizmoElasticStressFlux(int i, int jj, const double &f, double *Fout){
    const double dx = x[i] - x[jj], dy = y[i] - y[jj];
    const double r2 = dx*dx + dy*dy;
    if (r2 <= 0.) return;
    const double r = sqrt(r2), rinv = 1./r;

    // SPH 'effective area' along the line of centres (chosen to avoid tensile instability)
    const double dwk = std::fabs(Kernel::WDr(r, sml[i]) + Kernel::WDr(r, sml[jj]));
    const double rhoi = rho[i], rhoj = rho[jj];
    const double FNormT = m[i]*m[jj]*dwk / (rhoi*rhoj);
    const double FVec[2] = { FNormT*dx*rinv, FNormT*dy*rinv };

    // mean-velocity interface (the -0.5 carries GIZMO's required sign)
    const double v_int[2] = { -0.5*(vx[i]+vx[jj]), -0.5*(vy[i]+vy[jj]) };
    const double dv[2] = { vx[i]-vx[jj], vy[i]-vy[jj] };
    const double vdotr = dv[0]*dx + dv[1]*dy;

    // longitudinal (P-wave) and transverse (S-wave) acoustic impedances for HLL diffusion
#if EOS == 1 || EOS == 2
    const double Ki = MeshlessEOS->EOSBulkModulus(matId[i], rhoi, P[i]);
    const double Kj = MeshlessEOS->EOSBulkModulus(matId[jj], rhoj, P[jj]);
    const double mui = MeshlessEOS->EOSShearModulus(matId[i]);
    const double muj = MeshlessEOS->EOSShearModulus(matId[jj]);
#else
    const double Ki = MeshlessEOS->EOSBulkModulus(rhoi, P[i]);
    const double Kj = MeshlessEOS->EOSBulkModulus(rhoj, P[jj]);
    const double mui = 0., muj = 0.;
#endif
    const double ci  = sqrt((Ki + 4./3.*mui)/rhoi), cj  = sqrt((Kj + 4./3.*muj)/rhoj);
    const double cTi = sqrt(mui/rhoi),               cTj = sqrt(muj/rhoj);
    const double wt_r = rhoi*ci*rhoj*cj / (rhoi*ci + rhoj*cj) * FNormT;
    const double denT = rhoi*cTi + rhoj*cTj;
    const double wt_t = (denT > 0.) ? rhoi*cTi*rhoj*cTj / denT * FNormT : 0.;
    const double wt_rt = (wt_r - wt_t) * vdotr * rinv * rinv;

    // deviatoric stress force: project face onto stress eigenvectors, damp tensile
    // (positive) principal components by (1 - f); summed over both sides at wt=-0.5
    double cmag[2] = {0., 0.};
    for (int side = 0; side < 2; ++side){
        const int p = (side == 0) ? i : jj;
        double a = Sxx[p], b = Sxy[p], c = Syy[p];
#if FRAGMENTATION && DAMAGE_ACTS_ON_S
        // Grady-Kipp: damaged material carries less deviatoric stress (matches the
        // solid-HLLC path in Riemann.cpp). Pressure damage flows through the solver.
        const double dmg = 1.0 - damageTotal[p];
        a *= dmg; b *= dmg; c *= dmg;
#endif
        const double mean = 0.5*(a + c);
        const double dev  = sqrt(0.25*(a - c)*(a - c) + b*b);
        double l1 = mean + dev, l2 = mean - dev;
        double e1[2], e2[2];
        if (dev <= 1e-12*std::fabs(mean)){ e1[0]=1.; e1[1]=0.; e2[0]=0.; e2[1]=1.; }
        else {
            double v1x = b, v1y = l1 - a; const double wx = l1 - c, wy = b;
            if (wx*wx + wy*wy > v1x*v1x + v1y*v1y){ v1x = wx; v1y = wy; }
            const double nv = sqrt(v1x*v1x + v1y*v1y);
            e1[0] = v1x/nv; e1[1] = v1y/nv; e2[0] = -e1[1]; e2[1] = e1[0];
        }
        const double l[2] = { l1, l2 };
        const double *e[2] = { e1, e2 };
        for (int kk = 0; kk < 2; ++kk){
            double prefac = -0.5 * l[kk] * (FVec[0]*e[kk][0] + FVec[1]*e[kk][1]);
            if (l[kk] > 0.) prefac *= 1. - f;
            cmag[0] += prefac * e[kk][0];
            cmag[1] += prefac * e[kk][1];
        }
    }
    // HLL-type dissipation of velocity differences (longitudinal + transverse shear wave)
    cmag[0] -= wt_rt*dx + wt_t*dv[0];
    cmag[1] -= wt_rt*dy + wt_t*dv[1];
    const double cmag_E = cmag[0]*v_int[0] + cmag[1]*v_int[1];

    // demonstrator convention: dp_i/dt = -sum_j Fij, so add the negative of GIZMO's force
    Fout[1] -= cmag_E;
    Fout[2] -= cmag[0];
    Fout[3] -= cmag[1];
}
#endif // GIZMO_ELASTIC_FLUX

void Particles::integrateStressTensor(const double &dt) {
    for (int i = 0; i < N; ++i) {
#if EOS == 1 || EOS == 2
        const double mu = MeshlessEOS->EOSShearModulus(matId[i]);
#else
        const double mu = 0.;
#endif
        // Velocity gradients
        const double dvx_dx = vxGrad[i][0];
        const double dvx_dy = vxGrad[i][1];
        const double dvy_dx = vyGrad[i][0];
        const double dvy_dy = vyGrad[i][1];
#if DIM == 3
        const double dvx_dz = vxGrad[i][2];
        const double dvy_dz = vyGrad[i][2];
        const double dvz_dx = vzGrad[i][0];
        const double dvz_dy = vzGrad[i][1];
        const double dvz_dz = vzGrad[i][2];
#endif

        // Velocities
        const double vxi = vx[i];
        const double vyi = vy[i];
#if DIM == 3
        const double vzi = vz[i];
#endif

        // Current stress tensor values
        double sxx = Sxx[i];
        double syy = Syy[i];
        double sxy = Sxy[i];
#if DIM == 3
        double szz = Szz[i];
        double sxz = Sxz[i];
        double syz = Syz[i];
#endif

        // --- Constant terms (independent of S) ---

        // Traceless viscous strain rate (eq. 36-41, first line each)
#if DIM == 2
        const double visc_xx =  (4./3.)*mu*dvx_dx - (2./3.)*mu*dvy_dy;
        const double visc_yy = -(2./3.)*mu*dvx_dx + (4./3.)*mu*dvy_dy;
        const double visc_xy = mu * (dvy_dx + dvx_dy);
#elif DIM == 3
        const double visc_xx =  (4./3.)*mu*dvx_dx - (2./3.)*mu*dvy_dy - (2./3.)*mu*dvz_dz;
        const double visc_yy = -(2./3.)*mu*dvx_dx + (4./3.)*mu*dvy_dy - (2./3.)*mu*dvz_dz;
        const double visc_zz = -(2./3.)*mu*dvx_dx - (2./3.)*mu*dvy_dy + (4./3.)*mu*dvz_dz;
        const double visc_xy = mu * (dvy_dx + dvx_dy);
        const double visc_xz = mu * (dvz_dx + dvx_dz);
        const double visc_yz = mu * (dvz_dy + dvy_dz);
#endif

        // Advection: -v_gamma * d_gamma S^{alpha beta} (eq. 42)
#if DIM == 2
        const double adv_xx = -vxi*SxxGrad[i][0] - vyi*SxxGrad[i][1];
        const double adv_yy = -vxi*SyyGrad[i][0] - vyi*SyyGrad[i][1];
        const double adv_xy = -vxi*SxyGrad[i][0] - vyi*SxyGrad[i][1];
#elif DIM == 3
        const double adv_xx = -vxi*SxxGrad[i][0] - vyi*SxxGrad[i][1] - vzi*SxxGrad[i][2];
        const double adv_yy = -vxi*SyyGrad[i][0] - vyi*SyyGrad[i][1] - vzi*SyyGrad[i][2];
        const double adv_zz = -vxi*SzzGrad[i][0] - vyi*SzzGrad[i][1] - vzi*SzzGrad[i][2];
        const double adv_xy = -vxi*SxyGrad[i][0] - vyi*SxyGrad[i][1] - vzi*SxyGrad[i][2];
        const double adv_xz = -vxi*SxzGrad[i][0] - vyi*SxzGrad[i][1] - vzi*SxzGrad[i][2];
        const double adv_yz = -vxi*SyzGrad[i][0] - vyi*SyzGrad[i][1] - vzi*SyzGrad[i][2];
#endif

        // Rotation rate components (fixed during integration)
        const double R_xy = dvy_dx - dvx_dy; // d_x v_y - d_y v_x
#if DIM == 3
        const double R_xz = dvz_dx - dvx_dz; // d_x v_z - d_z v_x
        const double R_yz = dvz_dy - dvy_dz; // d_y v_z - d_z v_y
#endif

        // --- RK2 midpoint method ---
        // k1 = f(S^n): Jaumann rotation terms depend on current S
#if DIM == 2
        double k1_xx = visc_xx - sxy*R_xy;
        double k1_yy = visc_yy + sxy*R_xy;
        double k1_xy = visc_xy + 0.5*(sxx-syy)*R_xy;
#elif DIM == 3
        double k1_xx = visc_xx - sxy*R_xy - sxz*R_xz;
        double k1_yy = visc_yy + sxy*R_xy - syz*R_yz;
        double k1_zz = visc_zz + sxz*R_xz + syz*R_yz;
        double k1_xy = visc_xy + 0.5*(sxx-syy)*R_xy - 0.5*syz*R_xz - 0.5*sxz*R_yz;
        double k1_xz = visc_xz + 0.5*(sxx-szz)*R_xz - 0.5*syz*R_xy + 0.5*sxy*R_yz;
        double k1_yz = visc_yz + 0.5*(syy-szz)*R_yz + 0.5*sxz*R_xy + 0.5*sxy*R_xz;
#endif

        // S_mid = S^n + (dt/2) * k1
        double sxx_m = sxx + 0.5*dt*k1_xx;
        double syy_m = syy + 0.5*dt*k1_yy;
        double sxy_m = sxy + 0.5*dt*k1_xy;
#if DIM == 3
        double szz_m = szz + 0.5*dt*k1_zz;
        double sxz_m = sxz + 0.5*dt*k1_xz;
        double syz_m = syz + 0.5*dt*k1_yz;
#endif

        // k2 = f(S_mid): re-evaluate Jaumann terms at midpoint
#if DIM == 2
        double k2_xx = visc_xx - sxy_m*R_xy;
        double k2_yy = visc_yy + sxy_m*R_xy;
        double k2_xy = visc_xy + 0.5*(sxx_m-syy_m)*R_xy;
#elif DIM == 3
        double k2_xx = visc_xx - sxy_m*R_xy - sxz_m*R_xz;
        double k2_yy = visc_yy + sxy_m*R_xy - syz_m*R_yz;
        double k2_zz = visc_zz + sxz_m*R_xz + syz_m*R_yz;
        double k2_xy = visc_xy + 0.5*(sxx_m-syy_m)*R_xy - 0.5*syz_m*R_xz - 0.5*sxz_m*R_yz;
        double k2_xz = visc_xz + 0.5*(sxx_m-szz_m)*R_xz - 0.5*syz_m*R_xy + 0.5*sxy_m*R_yz;
        double k2_yz = visc_yz + 0.5*(syy_m-szz_m)*R_yz + 0.5*sxz_m*R_xy + 0.5*sxy_m*R_xz;
#endif

        // S^{n+1} = S^n + dt * k2
        Sxx[i] = sxx + dt*k2_xx;
        Syy[i] = syy + dt*k2_yy;
        Sxy[i] = sxy + dt*k2_xy;
#if DIM == 3
        Szz[i] = szz + dt*k2_zz;
        Sxz[i] = sxz + dt*k2_xz;
        Syz[i] = syz + dt*k2_yz;
#endif

#if PLASTICITY_ANY
        // --- Yield criterion: radial return to the von-Mises-equivalent Y ---
        {
#if PLASTICITY_MODEL_COUNT
            // Per-material yield strength Y(P, D, u). For Collins, D blends the
            // intact and damaged strength curves.
#if FRAGMENTATION
            const double Dtot = damageTotal[i];
#else
            const double Dtot = 0.;
#endif
            const double Y = MeshlessEOS->EOSYieldStrength(matId[i], P[i], Dtot, u[i]);
#else
            const double Y = YIELD_STRESS;   // legacy constant von Mises
#endif
#if DIM == 2
            const double J2 = 0.5 * (Sxx[i]*Sxx[i] + Syy[i]*Syy[i]
                                      + 2.0*Sxy[i]*Sxy[i]);
            const double Y2_over_d = Y * Y / 2.0;
#else
            const double J2 = 0.5 * (Sxx[i]*Sxx[i] + Syy[i]*Syy[i] + Szz[i]*Szz[i]
                                      + 2.0*(Sxy[i]*Sxy[i] + Sxz[i]*Sxz[i] + Syz[i]*Syz[i]));
            const double Y2_over_d = Y * Y / 3.0;
#endif
            // Radial return: scale S so that J2 -> Y2_over_d. Since J2 is
            // quadratic in S, the factor is the square root of the ratio.
            if (J2 > Y2_over_d && J2 > 0.0) {
                const double f = std::sqrt(Y2_over_d / J2);
                Sxx[i] *= f;  Syy[i] *= f;  Sxy[i] *= f;
#if DIM == 3
                Szz[i] *= f;  Sxz[i] *= f;  Syz[i] *= f;
#endif
            }
        }
#endif // PLASTICITY_ANY

    }
}

#if FRAGMENTATION
// Grady-Kipp damage evolution (Benz & Asphaug 1995). Advances the DIM-root
// damage variable `damage`; full damage is damageTotal = clamp(damage^DIM,0,1).
void Particles::integrateDamage(const double &dt){
    for (int i = 0; i < N; ++i){
        if (numFlaws[i] <= 0) { dddt[i] = 0.0; continue; }
        const double E   = MeshlessEOS->EOSYoungModulus(matId[i]);
        const double K   = MeshlessEOS->EOSGetMaterial(matId[i]).A;
        const double muS = MeshlessEOS->EOSShearModulus(matId[i]);
        const double Di  = damageTotal[i];

        // Maximum tensile principal stress of sigma = -P I + S.
        double sigma[DIM*DIM];
#if DIM == 2
        sigma[0] = -P[i] + Sxx[i]; sigma[1] = Sxy[i];
        sigma[2] = Sxy[i];         sigma[3] = -P[i] + Syy[i];
#else
        sigma[0]=-P[i]+Sxx[i]; sigma[1]=Sxy[i];       sigma[2]=Sxz[i];
        sigma[3]=Sxy[i];       sigma[4]=-P[i]+Syy[i]; sigma[5]=Syz[i];
        sigma[6]=Sxz[i];       sigma[7]=Syz[i];       sigma[8]=-P[i]+Szz[i];
#endif
        const double sigMax = Helper::maxEigenvalueSym(sigma);

        // Only tension (positive max principal stress) grows cracks.
        if (sigMax <= 0.0 || E <= 0.0) { dddt[i] = 0.0; continue; }

        // Local scalar strain, softened by the already-accumulated damage.
        const double local_strain = sigMax / ((1.0 - Di) * E);

        // Count flaws whose activation strain is exceeded (monotone in time).
        int nActive = 0;
        for (int k = 0; k < numFlaws[i]; ++k){
            if (flaws[i*MAX_NUM_FLAWS + k] < local_strain) ++nActive;
        }
        if (nActive > numActiveFlaws[i]) numActiveFlaws[i] = nActive;
        nActive = numActiveFlaws[i];

        if (nActive > 0){
            // Crack-growth speed: 0.4 x longitudinal elastic wave speed.
            const double cl2 = (K + 4.0/3.0 * muS * (1.0 - Di)) / rho[i];
            const double cg  = 0.4 * std::sqrt(cl2 > 0.0 ? cl2 : 0.0);
            dddt[i] = nActive * cg / sml[i];   // d(damage)/dt; damage is DIM-root
            damage[i] += dt * dddt[i];
        } else {
            dddt[i] = 0.0;
        }

        // Cap by the active-flaw fraction, then recompute the total damage.
        const double cap = std::pow((double)numActiveFlaws[i] / (double)numFlaws[i],
                                    1.0 / (double)DIM);
        if (damage[i] > cap) damage[i] = cap;
        if (damage[i] < 0.0) damage[i] = 0.0;
        double Dt = std::pow(damage[i], (double)DIM);
        if (Dt > 1.0) Dt = 1.0;
        damageTotal[i] = Dt;
    }
}
#endif // FRAGMENTATION
#endif // ELASTIC

#if !PERIODIC_BOUNDARIES
void Particles::getDomainLimits(double *domainLimits){

    double minX { std::numeric_limits<double>::max() };
    double maxX { std::numeric_limits<double>::min() };
    double minY { std::numeric_limits<double>::max() };
    double maxY { std::numeric_limits<double>::min() };
#if DIM == 3
    double minZ { std::numeric_limits<double>::max() };
    double maxZ { std::numeric_limits<double>::min() };
#endif

    for(int i=0; i<N; ++i){
        if (x[i] < minX){
            minX = x[i];
        } else if (x[i] > maxX){
            maxX = x[i];
        }
        if (y[i] < minY){
            minY = y[i];
        } else if (y[i] > maxY){
            maxY = y[i];
        }
#if DIM == 3
        if (z[i] < minZ){
            minZ = z[i];
        } else if (z[i] > maxZ){
            maxZ = z[i];
        }
#endif
    }

    domainLimits[0] = minX;
    domainLimits[DIM] = maxX;
    domainLimits[1] = minY;
    domainLimits[DIM+1] = maxY;
#if DIM == 3
    domainLimits[2] = minZ;
    domainLimits[DIM+2] = maxZ;
#endif
}
#endif // !PERIODIC_BOUNDARIES

void Particles::assignParticlesAndCells(Domain &domain){

    // reset particles assigned to grid cells
    for(int iGrid=0; iGrid<domain.numGridCells; ++iGrid){
        domain.grid[iGrid].prtcls = std::vector<int>();
    }

    for(int i=0; i<N; ++i){

        int floorX = floor((x[i]-domain.bounds.minX)/domain.cellSizeX);
        int floorY = floor((y[i]-domain.bounds.minY)/domain.cellSizeY);

        // This rarely happens when x or y is really close to the domain bounds
        if (floorX < 0){
            floorX = 0;
        } else if (floorX >= domain.cellsX){
            floorX = domain.cellsX - 1;
        }
        if (floorY < 0){
            floorY = 0;
        } else if (floorY >= domain.cellsY){
            floorY = domain.cellsY - 1;
        }

#if DIM==3
        int floorZ = floor((z[i]-domain.bounds.minZ)/domain.cellSizeZ);
        if (floorZ < 0){
            floorZ = 0;
        } else if (floorZ >= domain.cellsZ){
            floorZ = domain.cellsZ - 1;
        }
#endif


        int iGrid = floorX + floorY * domain.cellsX
#if DIM == 3
        + floorZ * domain.cellsX * domain.cellsY
#endif
        ;

        if (floorX < 0 || floorX >= domain.cellsX || floorY < 0 || floorY >= domain.cellsY){
           Logger(ERROR) << "  > Particle i = " << i << " cannot be properly assigned to search grid.";
           Logger(ERROR) << "    > x = " << x[i] << " -> floorX = " << floorX
                         << ", y = " << y[i] << " -> floorY = " << floorY;
        }
        //          << ", floor y = " << floorY;

        //Logger(DEBUG) << "      > Assigning particle@" << i << " = [" << x[i] << ", " << y[i]
//#if DIM == 3
//                      << ", " << z[i]
//#endif
//                      << "] to cell " << iGrid;

        domain.grid[iGrid].prtcls.push_back(i);
        cell[i] = iGrid; // assign cells to particles
    }
}

void Particles::gridNNS(Domain &domain, const double &kernelSize){
    int maxNoi = 0;
    long noiSum = 0;
    int nNearLimit = 0;
    // Warn threshold: log once per step how many particles are filling more
    // than ~80% of the nnl buffer, so the approach to MAX_NUM_INTERACTIONS
    // is visible well before the hard abort.
    const int nearLimit = (MAX_NUM_INTERACTIONS * 4) / 5;
    // loop over particles
    for(int i=0; i<N; ++i){
        int numSearchCells = pow(3, DIM);
        // search for nearest neighbors in the particle cell and neighbor cells
        int cells[numSearchCells];
        // neighboring cells (including particles cell)
        domain.getNeighborCells(cell[i], cells);
        // do nearest neighbor search
        int noiBuf = 0;
#if VARIABLE_SML
        // With variable smoothing lengths the grid is sized by the global
        // h_max, but a pair (i,j) only needs to be neighbours when the
        // distance is within max(sml[i], sml[j]).  Using the tighter
        // per-pair criterion keeps the neighbour lists compact and avoids
        // blowing MAX_NUM_INTERACTIONS for interior particles just because
        // a few boundary particles have large h.
        const double hi = sml[i];
#else
        const double hSqr = kernelSize * kernelSize;
#endif
        for (int iNeighbor=0; iNeighbor<numSearchCells; ++iNeighbor){
            // loop over particle indices in all
            if (cells[iNeighbor] < 0){
                // handle ghost cells in external function
            } else {
                for(auto const &iPrtcl : domain.grid[cells[iNeighbor]].prtcls){
                    if (iPrtcl != i){
                        double dSqr = pow(x[iPrtcl] - x[i], 2)
                                      + pow(y[iPrtcl] - y[i], 2);
#if DIM == 3
                        dSqr += pow(z[iPrtcl] - z[i], 2);
#endif
#if VARIABLE_SML
                        const double hPair = std::max(hi, sml[iPrtcl]);
                        const double hSqr = hPair * hPair;
#endif
                        if (dSqr < hSqr) {
                            if (noiBuf >= MAX_NUM_INTERACTIONS) {
                                Logger(ERROR) << "MAX_NUM_INTERACTIONS exceeded for particle "
                                              << i << " - Aborting.";
                                exit(1);
                            }
                            nnl[noiBuf + i * MAX_NUM_INTERACTIONS] = iPrtcl;
                            ++noiBuf;
                        }
                    }
                }
            }
        }
        if (noiBuf == 0){
            Logger(WARN) << "No neighbors for particle " << i << ". Caution.";
        }
        noi[i] = noiBuf;
        if (noiBuf > maxNoi) maxNoi = noiBuf;
        noiSum += noiBuf;
        if (noiBuf >= nearLimit) ++nNearLimit;
    }
    Logger(DEBUG) << "      > NNS noi[max/mean] = " << maxNoi << " / "
                  << ((double)noiSum / (double)N)
                  << " (MAX_NUM_INTERACTIONS = " << MAX_NUM_INTERACTIONS << ")";
    if (nNearLimit > 0){
        Logger(WARN) << "gridNNS: " << nNearLimit << " particles with noi >= "
                     << nearLimit << "/" << MAX_NUM_INTERACTIONS
                     << " (max noi = " << maxNoi
                     << "). Approaching the hard limit.";
    }
}

#if RUNSPH
// For comparable ICs: This sets the internal energies so that P = 2.5 everywhere
void Particles::setInternalEnergy(double Pressure, const double hydro_gamma){
    for (int i = 0; i < N; i++){
#if EOS == 1 || EOS == 2
        u[i] = MeshlessEOS->EOSInternalEnergy(matId[i], rho[i], Pressure);
#else
        u[i] = MeshlessEOS->EOSInternalEnergy(rho[i], Pressure);
#endif
    }
}

// Computes the density via kernel smoothing w/o ghost particles
void Particles::compDensitySPH(const double &kernelSize){
    double dnst;
    double dSqr;
    double r;
    int iP;
    for (int i = 0; i < N; i++){
        dnst = 0;
        dSqr = 0;
        r = 0.;
        for(int j = 0; j < noi[i]; j++){
            iP = nnl[j+i*MAX_NUM_INTERACTIONS];
            dSqr = pow(x[i] - x[iP], 2)
                        + pow(y[i] - y[iP], 2);
#if DIM == 3
            dSqr += pow(z[i] - z[iP], 2);
#endif
            r = sqrt(dSqr);
            dnst += m[i] * kernel(r, kernelSize);
            // For normalization, as in compOmega: Add self interaction
            rho[i] = dnst + m[j] * kernel(0., kernelSize);
        }
    }
}


void Particles::compAccSPH(const double &kernelSize){
    int iP;
    double dSqr;
    double r;
    double PRhoHost;
    double PRhoTarget;
    double PIij;
    for (int i = 0; i < N; i++){
        ax[i] = 0;
        ay[i] = 0;
        #if DIM == 3
        az[i] = 0;
        #endif
        PRhoHost = P[i] / pow(rho[i], 2);
        for (int j = 0; j < noi[i]; j++){
            iP = nnl[j+i*MAX_NUM_INTERACTIONS];

            dSqr = pow(x[i] - x[iP], 2)
            + pow(y[i] - y[iP], 2);
            #if DIM == 3
            dSqr += pow(z[i] - z[iP], 2);
            #endif
            r = sqrt(dSqr);
            PRhoTarget =  P[j] / pow(rho[j], 2);
            #if ARTVISC
            PIij = compPIij(i, iP, kernelSize);
            ax[i] += m[iP] * (PRhoHost + PRhoTarget + PIij)  * (x[iP] - x[i])/r * Kernel::dWdr(r, kernelSize);
            ay[i] += m[iP] * (PRhoHost + PRhoTarget + PIij)  * (y[iP] - y[i])/r * Kernel::dWdr(r, kernelSize);

            //ax[i] -= m[j] * (PRhoHost + PRhoTarget) * (x[i] - x[j])/r * Kernel::dWdr(r, kernelSize);
            //ay[i] -= m[j] * (PRhoHost + PRhoTarget) * (y[i] - y[j])/r * Kernel::dWdr(r, kernelSize);
            // Logger(DEBUG) << i << " " << j << " " << m[j]*(PRhoHost+PRhoTarget)*(x[i] - x[j])/r * Kernel::dWdr(r, kernelSize);
            #if DIM == 3
            az[i] -= m[j]*(PRhoHost+PRhoTarget+PIij)*(z[i] - z[j])/r * Kernel::dWdr(r, kernelSize);
            #endif
            #else
            ax[i] += m[iP] * (PRhoHost+PRhoTarget)  * (x[iP] - x[i])/r * Kernel::dWdr(r, kernelSize);
            ay[i] += m[iP] * (PRhoHost+PRhoTarget)  * (y[iP] - y[i])/r * Kernel::dWdr(r, kernelSize);

            //ax[i] -= m[j] * (PRhoHost + PRhoTarget) * (x[i] - x[j])/r * Kernel::dWdr(r, kernelSize);
            //ay[i] -= m[j] * (PRhoHost + PRhoTarget) * (y[i] - y[j])/r * Kernel::dWdr(r, kernelSize);
            // Logger(DEBUG) << i << " " << j << " " << m[j]*(PRhoHost+PRhoTarget)*(x[i] - x[j])/r * Kernel::dWdr(r, kernelSize);
            #if DIM == 3
            az[i] -= m[j]*(PRhoHost+PRhoTarget)*(z[i] - z[j])/r * Kernel::dWdr(r, kernelSize);
            #endif
            ax[i] -= m[i]*(PRhoHost+PRhoTarget)*(x[i] - x[j])/r * Kernel::dWdr(r, kernelSize);
            ay[i] -= m[i]*(PRhoHost+PRhoTarget)*(y[i] - y[j])/r * Kernel::dWdr(r, kernelSize);
            #if DIM == 3
            az[i] -= m[i]*(PRhoHost+PRhoTarget)*(z[i] - z[j])/r * Kernel::dWdr(r, kernelSize);
            #endif
            #endif
        }
    }
}


// Euler intgration for SPH
void Particles::eulerSPH(const double &dt, const Domain &domain){
    for (int i = 0; i < N; i++){
        vx[i] += ax[i]*dt;
        vy[i] += ay[i]*dt;
#if DIM == 3
        vz[i] += dt * az[i];
#endif

        x[i] += vx[i]*dt;
        y[i] += vy[i]*dt;
#if DIM == 3
        z[i] += vz[i]*dt;
#endif

#if PERIODIC_BOUNDARIES
        if (x[i] < domain.bounds.minX) {
            x[i] = domain.bounds.maxX - (domain.bounds.minX - x[i]);
        }
        else if (domain.bounds.maxX <= x[i]) {
            //Logger(DEBUG) << "X is " << x[i];
            x[i] = domain.bounds.minX + (x[i] - domain.bounds.maxX);
        }
        if (y[i] < domain.bounds.minY) {
            y[i] = domain.bounds.maxY - (domain.bounds.minY - y[i]);
        }
        else if (domain.bounds.maxY <= y[i]) {
            y[i] = domain.bounds.minY + (y[i] - domain.bounds.maxY);
        }
#if DIM ==3
        if (z[i] < domain.bounds.minZ) {
            z[i] = domain.bounds.maxZ - (domain.bounds.minZ - z[i]);
        }
        else if (domain.bounds.maxZ <= z[i]) {
            z[i] = domain.bounds.minZ + (z[i] - domain.bounds.maxZ);
        }
#endif
#endif

        // Check if particles are out of bounds
        // if (x[i] > domain.bounds.maxX){
            //     Logger(DEBUG) << "X is too big, out of bounds";
            // }
    }
}


// SPH energy computation functions:
void Particles::compuis(const double &dt, const double &kernelSize){
    //double u_i_old;
    for (int i = 0; i < N; i++){
        //u_i_old = u[i];
        u[i] += dEdt[i] * dt;
#if ARTVISC
        u[i] += dudtArtVisc[i] * dt;
#endif
        //std::cout << "u_" << i << ": " << u[i] << ", was " << u_i_old << std::endl;;
        if (u[i] < 0.){
            Logger(WARN) << "+++DANGER+++ negative specific internal energy encountered.";
        }
#if ENERGY_FLOOR >= 0
        if (u[i] < ENERGY_FLOOR){
            u[i] = ENERGY_FLOOR;
        }
#endif
    }
}


// Computing all omegas:
void Particles::compOmegas(const double &kernelSize){
    (void)kernelSize; // sml[i] is used inside compOmega
    for (int i = 0; i < N; i++){
        compOmega(i);
    }
}

// Implementing equations F5 and F6 in Hopkins` GIZMO Paper
void Particles::calcdndrho(const double &kernelSize){
    double r, dSqr, tmp;
    int iP;
    for (int i = 0; i < N; i++){
        dn[i] = 0;
        drho[i] = 0;
        for (int j = 0; j < noi[i]; j++){
            iP = nnl[j+i*MAX_NUM_INTERACTIONS];
            dSqr = pow(x[i] - x[iP], 2)
                        + pow(y[i] - y[iP], 2);
#if DIM == 3
            dSqr += pow(z[i] - z[iP], 2);
#endif
            r = sqrt(dSqr);
            dn[i] -= 1 / kernelSize * (DIM * Kernel::cubicSpline(r, kernelSize) + r*Kernel::dWdr(r, kernelSize));
            drho[i] -= m[iP]*1 / kernelSize * (DIM * Kernel::cubicSpline(r, kernelSize) + r*Kernel::dWdr(r, kernelSize));
        }
    }
}

// Implementing the sum in equation F3 in Hopkins` GIZMO Paper
void Particles::calcdE(const double &kernelSize){
    double tmp, dSqr, r, fij;
    int iP;
    for (int i = 0; i < N; i++){
        dEdt[i] = 0;
        for (int j = 0; j < noi[i]; j++){
            iP = nnl[j+i*MAX_NUM_INTERACTIONS];
            dSqr = pow(x[i] - x[iP], 2)
                        + pow(y[i] - y[iP], 2);
#if DIM == 3
            dSqr += pow(z[i] - z[iP], 2);
#endif
            r = sqrt(dSqr);
            fij = 1 - 1 / m[iP]
                    * kernelSize / (omega[i] * DIM) * drho[i]
                        / (1 + kernelSize /(omega[i] * DIM) * dn[i]);
            //fij = 1;
            dEdt[i]  += m[iP]
                        * ((vx[i]-vx[iP])
                        * (x[i]-x[iP])
                        + (vy[i]-vy[iP])
                        * (y[i]-y[iP])
#if DIM == 3
                        + (vz[i]-vz[iP])*(z[i]-z[iP])
#endif
                        ) / r * P[i]/pow(rho[i],2)*fij*Kernel::dWdr(r, kernelSize)
                        ;
        }
    }
}



#if PERIODIC_BOUNDARIES

// Computes the density via kernel smoothing w/ ghost particles
void Particles::compDensitySPH(const Particles &ghostParticles, const double &kernelSize){
    int iP;
    double dnst;
    double dSqr;
    double r;
    double ghostMass;
    for (int i = 0; i < N; ++i){
        dnst = 0;
        for(int j = 0; j < noi[i]; j++){
            iP = nnl[j+i*MAX_NUM_INTERACTIONS];
            dSqr = pow((x[i] - x[iP]), 2)
                        + pow((y[i] - y[iP]), 2);
#if DIM == 3
            dSqr += pow(z[i] - z[iP], 2);
#endif
            r = sqrt(dSqr);

            //Logger(DEBUG) << "density: i = " << i << ", j = " << iP << ", r = " << r;

            dnst += m[iP] * kernel(r, kernelSize);
        }
        // For normalization: Self interaction
        dnst += m[i] * kernel(0., kernelSize);

        // Ghosts
        //Logger(DEBUG) << "i: " << i << " noiGhosts: " << noiGhosts[i];
        //Logger(DEBUG) << "Random ghost mass: " << ghostParticles.m[0];
#if PERIODIC_BOUNDARIES
        for(int k = 0; k < noiGhosts[i]; ++k){
            iP = nnlGhosts[k+i*MAX_NUM_GHOST_INTERACTIONS];
            dSqr = pow(x[i] - ghostParticles.x[iP], 2)
                          + pow(y[i] - ghostParticles.y[iP], 2);
#if DIM == 3
            dSqr += pow(z[i] - ghostParticles.z[iP], 2);
#endif
            r = sqrt(dSqr);
            // TODO: include mass in update ghost state
            dnst += m[ghostParticles.parent[iP]] * kernel(r, kernelSize);
            //Logger(DEBUG) << "k = " << k << " mass:  " << m[ghostParticles.parent[iP]] << " W: " << kernel(r,kernelSize);
        }
#endif
        // if ((i - 0) % 30 == 0){
        //     Logger(DEBUG) << "density from ghosts: i = " << i << ", dnst = " << dnst;
        // }
        rho[i] = dnst;
    }
}


void Particles::compAccSPH(const Particles &ghostParticles, const double &kernelSize){
    int iP;
    // Calculate acceleration for each particle, w/o Ghosts
    // c.f. eq 8 in Monaghan: SPH and its diverse Applications, Annu.Rev. Fluid mechanics, 2012
    double dSqr;
    double r;
    double PRhoHost;
    double PRhoTarget;
    double PIij;
    for (int i = 0; i < N; i++){
        ax[i] = 0;
        ay[i] = 0;
#if DIM == 3
        az[i] = 0;
#endif
        PRhoHost = P[i] / pow(rho[i], 2);
        for (int j = 0; j < noi[i]; j++){
            iP = nnl[j+i*MAX_NUM_INTERACTIONS];

            dSqr = pow(x[i] - x[iP], 2)
                        + pow(y[i] - y[iP], 2);
#if DIM == 3
            dSqr += pow(z[i] - z[iP], 2);
#endif
            r = sqrt(dSqr);
            PRhoTarget =  P[iP] / pow(rho[iP], 2);
#if ARTVISC
            PIij = compPIij(i, iP, kernelSize); // TODO: what is 1.5 and 3. ???
            ax[i] += m[iP] * (PRhoHost + PRhoTarget + PIij)  * (x[iP] - x[i])/r * Kernel::dWdr(r, kernelSize);
            ay[i] += m[iP] * (PRhoHost + PRhoTarget + PIij)  * (y[iP] - y[i])/r * Kernel::dWdr(r, kernelSize);

            //ax[i] -= m[j] * (PRhoHost + PRhoTarget) * (x[i] - x[j])/r * Kernel::dWdr(r, kernelSize);
            //ay[i] -= m[j] * (PRhoHost + PRhoTarget) * (y[i] - y[j])/r * Kernel::dWdr(r, kernelSize);
            // Logger(DEBUG) << i << " " << j << " " << m[j]*(PRhoHost+PRhoTarget)*(x[i] - x[j])/r * Kernel::dWdr(r, kernelSize);
#if DIM == 3
            az[i] -= m[j]*(PRhoHost+PRhoTarget+PIij)*(z[i] - z[j])/r * Kernel::dWdr(r, kernelSize);
#endif
#else
            ax[i] += m[iP] * (PRhoHost+PRhoTarget)  * (x[iP] - x[i])/r * Kernel::dWdr(r, kernelSize);
            ay[i] += m[iP] * (PRhoHost+PRhoTarget)  * (y[iP] - y[i])/r * Kernel::dWdr(r, kernelSize);

//ax[i] -= m[j] * (PRhoHost + PRhoTarget) * (x[i] - x[j])/r * Kernel::dWdr(r, kernelSize);
//ay[i] -= m[j] * (PRhoHost + PRhoTarget) * (y[i] - y[j])/r * Kernel::dWdr(r, kernelSize);
// Logger(DEBUG) << i << " " << j << " " << m[j]*(PRhoHost+PRhoTarget)*(x[i] - x[j])/r * Kernel::dWdr(r, kernelSize);
#if DIM == 3
            az[i] -= m[j]*(PRhoHost+PRhoTarget)*(z[i] - z[j])/r * Kernel::dWdr(r, kernelSize);
#endif
#endif
        }

        // Ghosts
        for (int k = 0; k < noiGhosts[i]; k++){
            iP = nnlGhosts[k+i*MAX_NUM_GHOST_INTERACTIONS];
            dSqr = pow(x[i] - ghostParticles.x[iP], 2)
                          + pow(y[i] - ghostParticles.y[iP], 2);
#if DIM == 3
            dSqr += pow(z[i] - ghostParticles.z[iP], 2);
#endif
            r = sqrt(dSqr);
            PRhoTarget = ghostParticles.P[iP]
                / pow(ghostParticles.rho[iP], 2);
#if ArtVisc
            PIij = compPIij(ghostParticles, i, iP, kernelSize);
            ax[i] -= m[ghostParticles.parent[iP]]*(PRhoHost+PRhoTarget+PIij)
                * (x[i] - ghostParticles.x[iP])
                    / r * Kernel::dWdr(r, kernelSize);
            ay[i] -= m[ghostParticles.parent[iP]]*(PRhoHost+PRhoTarget+PIij)
                * (y[i] - ghostParticles.y[iP])
                    / r * Kernel::dWdr(r, kernelSize);
#if DIM == 3
            az[i] -= m[ghostParticles.parent[iP]]*(PRhoHost+PRhoTarget+PIij)
                * (z[i] - ghostParticles.z[iP]) / r * Kernel::dWdr(r, kernelSize);
#endif
#else
            ax[i] -= m[ghostParticles.parent[iP]]*(PRhoHost+PRhoTarget)
                * (x[i] - ghostParticles.x[iP])
                    / r * Kernel::dWdr(r, kernelSize);
            ay[i] -= m[ghostParticles.parent[iP]]*(PRhoHost+PRhoTarget)
                * (y[i] - ghostParticles.y[iP])
                    / r * Kernel::dWdr(r, kernelSize);
#if DIM == 3
            az[i] -= m[ghostParticles.parent[iP]]*(PRhoHost+PRhoTarget)
                * (z[i] - ghostParticles.z[iP]) / r * Kernel::dWdr(r, kernelSize);
#endif
#endif
        }
        // Logger(DEBUG) << "For i = " << i << ", ax and ay are " << ax[i] << " and " << ay[i];
    }
}

// Compute all omegas
void Particles::compOmegas(const Particles &ghostParticles, const double &kernelSize){
    (void)kernelSize; // sml[i] is used inside compOmega
    for (int i = 0; i < N; i++){
        compOmega(i, ghostParticles);
        //std::cout << i << " " << omega[i] << std::endl;
    }
}

// Implementing equations F5 and F6 in Hopkins` GIZMO Paper
void Particles::calcdndrho(const Particles &ghostParticles, const double &kernelSize){
    int iP;
    double r, dSqr, tmp;
    for (int i = 0; i < N; i++){
        dn[i] = 0;
        drho[i] = 0;
        for (int j = 0; j < noi[i]; j++){
            iP = nnl[j+i*MAX_NUM_INTERACTIONS];
            dSqr = pow(x[i] - x[iP], 2)
                        + pow(y[i] - y[iP], 2);
#if DIM == 3
            dSqr += pow(z[i] - z[iP], 2);
#endif
            r = sqrt(dSqr);

            dn[i] -= 1 / kernelSize
                * (DIM *Kernel::cubicSpline(r, kernelSize) + r / kernelSize * Kernel::dWdh(r, kernelSize));
            drho[i] -= m[iP] / kernelSize
                * (DIM *Kernel::cubicSpline(r, kernelSize) + r / kernelSize * Kernel::dWdh(r, kernelSize));
        }
        for (int k = 0; k < noiGhosts[i]; k++){
            iP = nnlGhosts[k+i*MAX_NUM_GHOST_INTERACTIONS];
            dSqr = pow(x[i] - ghostParticles.x[iP], 2)
            + pow(y[i] - ghostParticles.y[iP], 2);
#if DIM == 3
            dSqr += pow(z[i] - ghostParticles.z[iP], 2);
#endif
            r = sqrt(dSqr);
            dn[i] -= 1 / kernelSize
                * (DIM * Kernel::cubicSpline(r, kernelSize) + r/kernelSize * Kernel::dWdh(r, kernelSize));
            drho[i] -= m[ghostParticles.parent[iP]] / kernelSize
                * (DIM * Kernel::cubicSpline(r, kernelSize) + r/kernelSize * Kernel::dWdh(r, kernelSize));
        }
        //std::cout << "For i = " << i <<", drho is " << drho[i] << " dn " << dn[i] << " m " << m[i] << std::endl;
    }
}


// Implementing equation F3 in Hopkins` GIZMO Paper
void Particles::calcdE(const Particles &ghostParticles, const double &kernelSize){
    int iP;
    double dSqr, r, fij;
    for (int i = 0; i < N; i++){
        dEdt[i] = 0;
        for (int j = 0; j < noi[i]; j++){
            iP = nnl[j+i*MAX_NUM_INTERACTIONS];
            dSqr = pow(x[i] - x[iP], 2)
                        + pow(y[i] - y[iP], 2);
#if DIM == 3
            dSqr += pow(z[i] - z[iP], 2);
#endif
            r = sqrt(dSqr);
            fij = 1 - 1 / m[iP]
                    * kernelSize / (omega[i] * DIM) * drho[i]
                        / (1 + kernelSize /(omega[i] * DIM) * dn[i]);
            //fij = 1;
            dEdt[i]  += m[iP]
                        * ((vx[i]-vx[iP])
                        * (x[i]-x[iP])
                        + (vy[i]-vy[iP])
                        * (y[i]-y[iP])
#if DIM == 3
                        + (vz[i]-vz[iP])*(z[i]-z[iP])
#endif
                        ) / r * P[i]/pow(rho[i],2)*fij*Kernel::dWdr(r, kernelSize)
                        ;
        }
        // std::cout << "For i = " << i <<", fij is " << fij << std::endl;
        for (int k = 0; k < noiGhosts[i]; k++){
            iP = nnlGhosts[k+i*MAX_NUM_GHOST_INTERACTIONS];
            dSqr = pow(x[i] - ghostParticles.x[iP], 2)
            + pow(y[i] - ghostParticles.y[iP], 2);
#if DIM == 3
            dSqr += pow(z[i] - ghostParticles.z[iP], 2);
#endif
            r = sqrt(dSqr);

            fij = 1 - 1 / m[ghostParticles.parent[iP]]
                    * kernelSize / omega[i] / DIM * drho[i] / (1 + kernelSize /omega[i] / DIM * dn[i]);
            //fij = 1;
            dEdt[i]  += m[ghostParticles.parent[iP]]
                        * ((vx[i]-ghostParticles.vx[iP])* (x[i]-ghostParticles.x[iP])
                        + (vy[i]-ghostParticles.vy[iP]) * (y[i]-ghostParticles.y[iP])
#if DIM == 3
                        + (vz[i]-vz[nnl[k+i*MAX_NUM_GHOST_INTERACTIONS]])*(z[i]-z[nnl[k+i*MAX_NUM_GHOST_INTERACTIONS]])
#endif
                        ) / r * P[i]/pow(rho[i],2)*fij*Kernel::dWdr(r, kernelSize)
                    ;
        }
        //std::cout << "dE for i=" << i << " is " << dEdt[i] << std::endl;
    }
}

#endif

#if ARTVISC
// Artificial Viscocity!
// Compute cs
void Particles::compCs(const double hydro_gamma){
    for (int i = 0; i < N; i++){
        cs[i] = sqrt((hydro_gamma - 1) * u[i]);
    }
}

// Compute Mu_ij
double Particles::compMuij(int i, int j, const double &kernelSize){
    double numerator;
    numerator = ((vx[i] - vx[j]) * (x[i] - x[j])
            + (vy[i] - vy[j]) * (y[i] - y[j]));
#if DIM == 3
    numerator += (vz[i] - vz[j]) * (z[i] - z[j]);
#endif

    double rSqr;
    rSqr = pow(x[i] - x[j], 2) + pow(y[i] - y[j], 2);
#if DIM == 3
    rSqr += pow(z[i] - z[j], 2);
#endif

    return kernelSize * numerator / (sqrt(rSqr) + EPSMU * pow(kernelSize, 2));
}

// Compute PI_ij:
double Particles::compPIij(int i, int j, const double &kernelSize){
    double numerator;
    numerator = ((vx[i] - vx[j]) * (x[i] - x[j])
    + (vy[i] - vy[j]) * (y[i] - y[j]));
    #if DIM == 3
    numerator += (vz[i] - vz[j]) * (z[i] - z[j]);
    #endif

    if (numerator < 0){
        return 0;
    }
    else{
        double Muij = compMuij(i, j, kernelSize);

        double cij = (cs[i] + cs[j])/2;
        double rhoij = (rho[i] + rho[j]) / 2;

        return (- ALPHA_VISC * cij * Muij + BETA_VISC * pow(Muij, 2)) / rhoij;
    }
}

// Compute additional acceleration term for artificial viscocity
void Particles::compAccArtVisc(const double &kernelSize){
    double PIij, dSqr, r;
    // Interaction partner:
    int iP;
    for (int i = 0; i < N; i++){
        axArtVisc[i] = 0;
        ayArtVisc[i] = 0;
        #if DIM == 3
        azArtVisc[i] = 0;
        #endif
        for (int j = 0; j < noi[i]; j++){
            iP = nnl[j+i*MAX_NUM_INTERACTIONS];

            dSqr = pow(x[i] - x[iP], 2)
            + pow(y[i] - y[iP], 2);
            #if DIM == 3
            dSqr += pow(z[i] - z[iP], 2);
            #endif
            r = sqrt(dSqr);
            PIij = compPIij(i, iP, kernelSize);
            //Logger(DEBUG) << "PIij is " << PIij;
            axArtVisc[i] -= m[iP] * PIij * (x[i] - x[iP]) / r *  Kernel::dWdr(r, kernelSize);
            axArtVisc[i] -= m[iP] * PIij * (y[i] - y[iP]) / r *  Kernel::dWdr(r, kernelSize);
            #if DIM == 3
            azArtVisc[i] -= m[iP] * PIij * (z[i] - z[iP]) / r *  Kernel::dWdr(r, kernelSize);
            #endif
        }
    }
}

// Compute additional energy terms for artificial viscocity
void Particles::compUiArtVisc(const double &kernelSize){
    double PIij, dSqr, r;
    // Interaction partner:
    int iP;
    for (int i = 0; i < N; i++){
        dudtArtVisc[i] = 0;
        for (int j = 0; j < noi[i]; j++){
            iP = nnl[j+i*MAX_NUM_INTERACTIONS];

            dSqr = pow(x[i] - x[iP], 2)
            + pow(y[i] - y[iP], 2);
            #if DIM == 3
            dSqr += pow(z[i] - z[iP], 2);
            #endif
            r = sqrt(dSqr);
            // This is not really needed
            if (r <= 0){
                Logger(DEBUG) << "DANGER, r < 0";
            }
            PIij = compPIij(i, iP, kernelSize);
            dudtArtVisc[i] += .5 * m[iP] * PIij *  Kernel::dWdr(r, kernelSize)
            * ((vx[i] - vx[iP]) * (x[i] - x[iP])
            + (vy[i] - vy[iP]) * (y[i] - y[iP])
            #if DIM == 3
            + (vz[i] - vy[iP]) * (z[i] - z[iP])
            #endif
        ) / r;
    }
}
}

#if PERIODIC_BOUNDARIES
// Compute Mu_ij
double Particles::compMuij(const Particles &ghostParticles, int i, int j, const double &kernelSize){
    double numerator;
    numerator = ((vx[i] - ghostParticles.vx[j]) * (x[i] - ghostParticles.x[j])
            + (vy[i] - ghostParticles.vy[j]) * (y[i] - ghostParticles.y[j]));
#if DIM == 3
    numerator += (vz[i] - ghostParticles.vz[j]) * (z[i] - ghostParticles.z[j]);
#endif

    double rSqr;
    rSqr = pow(x[i] - ghostParticles.x[j], 2) + pow(y[i] - ghostParticles.y[j], 2);
#if DIM == 3
    rSqr += pow(z[i] - ghostParticles.z[j], 2);
#endif
    return (kernelSize * numerator / (sqrt(rSqr) + 0.000025 * pow(kernelSize, 2)));
}


// Compute PI_ij:
double Particles::compPIij(const Particles &ghostParticles, int i, int j, const double &kernelSize){
    double numerator;
    numerator = ((vx[i] - ghostParticles.vx[j]) * (x[i] - ghostParticles.x[j])
            + (vy[i] - ghostParticles.vy[j]) * (y[i] - ghostParticles.y[j]));
#if DIM == 3
    numerator += (vz[i] - ghostParticles.vz[j]) * (z[i] - ghostParticles.z[j]);
#endif

    if (numerator < 0){
        return 0;
    }
    else{
        double Muij = compMuij(ghostParticles, i, j, kernelSize);

        double cij = (cs[i] + cs[ghostParticles.parent[j]])/2;
        double rhoij = (rho[i] + ghostParticles.rho[j]) / 2;

        return (- ALPHA_VISC * cij * Muij + BETA_VISC * pow(Muij, 2)) / rhoij;
    }
}

// Compute additional acceleration term for artificial viscocity w/ periodic boundaries
void Particles::compAccArtVisc(const Particles &ghostParticles, const double &kernelSize){
    double PIij, dSqr, r;
    // Interaction partner:
    int iP;
    for (int i = 0; i < N; i++){
        for (int j = 0; j < noiGhosts[i]; j++){
            iP = nnlGhosts[j+i*MAX_NUM_GHOST_INTERACTIONS];

            dSqr = pow(x[i] - ghostParticles.x[iP], 2)
                        + pow(y[i] - ghostParticles.y[iP], 2);
#if DIM == 3
            dSqr += pow(z[i] - ghostParticles.z[iP], 2);
#endif
            r = sqrt(dSqr);
            PIij = compPIij(ghostParticles, i, iP, kernelSize);
            //Logger(DEBUG) << "	> PIij ghost is " << PIij;
            axArtVisc[i] -= m[ghostParticles.parent[iP]] * PIij * (x[i] - ghostParticles.x[iP]) / r *  Kernel::dWdr(r, kernelSize);
            axArtVisc[i] -= m[ghostParticles.parent[iP]] * PIij * (y[i] - ghostParticles.y[iP]) / r *  Kernel::dWdr(r, kernelSize);
#if DIM == 3
            azArtVisc[i] -= m[ghostParticles.parent[iP]] * PIij * (z[i] - ghostParticles.z[iP]) / r *  Kernel::dWdr(r, kernelSize);
#endif
        }
    }
}



// Compute additional energy terms for artificial viscocity
void Particles::compUiArtVisc(const Particles &ghostParticles, const double &kernelSize){
    double PIij, dSqr, r;
    // Interaction partner:
    int iP;
    for (int i = 0; i < N; i++){
        for (int j = 0; j < noiGhosts[i]; j++){
            iP = nnlGhosts[j+i*MAX_NUM_GHOST_INTERACTIONS];

            dSqr = pow(x[i] - ghostParticles.x[iP], 2)
                        + pow(y[i] - ghostParticles.y[iP], 2);
#if DIM == 3
            dSqr += pow(z[i] - ghostParticles.z[iP], 2);
#endif
            r = sqrt(dSqr);
            PIij = compPIij(ghostParticles, i, iP, kernelSize);
            dudtArtVisc[i] += .5 * m[ghostParticles.parent[iP]] * PIij *  Kernel::dWdr(r, kernelSize)
                * ((vx[i] - ghostParticles.vx[iP]) * (x[i] - ghostParticles.x[iP])
                    + (vy[i] - ghostParticles.vy[iP]) * (y[i] - ghostParticles.y[iP])
#if DIM == 3
                    + (vz[i] - ghostParticles.vy[iP]) * (z[i] - ghostParticles.z[iP])
#endif
                ) / r;
        }
    }
}

#endif

#endif
#endif // RUNSPH
// Calculate acceleration for each particle, w/o Ghosts
// c.f. eq 8 in Monaghan: SPH and its diverse Applications, Annu.Rev. Fluid mechanics, 2012


// void Particles::printDensity(const double &gamma){
//     for (int i = 0; i < N; i += 25){
//         Logger(DEBUG) << "i=" << i << ":rho= " << rho[i] << ", u= " << u[i] << "< p= " << P[i] << "test=" << rho[i]*u[i]*(1 - gamma);
//     }
// }




// Todo: Remove
// Implementing equation F2 in Hopkins` GIZMO Paper
// Storing v_i*dP_i/dt in the variable unimportant, so that only a scalar is needed
// void Particles::calcdP(const Particles &ghostParticles, const double &kernelSize){
//     int iP;
//     double dSqr, r, fij, fji, summand;
//     for (int i = 0; i < N; i++){
//         unimportant[i] = 0;
//         for (int j = 0; j < noi[i]; j++){
//             dSqr = pow((x[i] - x[iP]), 2)
//                         + pow((y[i] - y[iP]), 2);
// #if DIM == 3
//             dSqr += pow(z[i] - z[iP], 2);
// #endif
//             r = sqrt(dSqr);
//
//             fij = 1 - 1 / m[iP]
//                     * kernelSize /omega[i] / DIM * drho[i]
//                     / (1 + kernelSize /omega[i] / DIM * dn[i]);
//
//             summand = m[i]*m[iP] * Kernel::dWdr(r, kernelSize)
//                     * (P[i]/pow(rho[i], 2)*fij
//                         +  P[iP]/pow(rho[iP], 2)*fji);
//
//
//             unimportant[i] -= summand * (vx[i]*(x[i]- x[iP]) + vy[i]*(y[i]-y[iP])
// #if DIM == 3
//                             + vz[i]*(z[i]- z[iP])
// #endif
//                             ) / r;
//         }
//         // Ghosts!
//         for (int k = 0; k < noiGhosts[i]; k++){
//             dSqr = pow((x[i] - x[nnl[k+i*MAX_NUM_INTERACTIONS]]), 2)
//                         + pow((y[i] - y[iP]), 2);
// #if DIM == 3
//             dSqr += pow(z[i] - z[iP], 2);
// #endif
//             r = sqrt(dSqr);
//
//             fij = 1 - 1 / m[ghostParticles.parent[iP]]
//                     * kernelSize /omega[i] / DIM * drho[i]
//                     / (1 + kernelSize /omega[i] / DIM * dn[i]);
//
//             fji = 1 - 1 / m[i]
//                     * kernelSize / omega[iP] / DIM * drho[iP]
//                     / ( 1 + kernelSize / omega[iP] / DIM * dn[iP]);
//
//             summand = m[i]*m[ghostParticles.parent[iP]] * Kernel::dWdr(r, kernelSize)
//                     * (P[i]/pow(rho[i], 2)*fij
//                         +  P[iP]/pow(rho[iP], 2)*fji);
//
//
//             unimportant[i] -= summand * (vx[i]*(x[i]- x[iP]) + vy[i]*(y[i]-y[iP])
// #if DIM == 3
//                             + vz[i]*(z[i]- z[iP])
// #endif
//                             ) / r;
//         }
//     }
// }
//
//


// Implementing equation F2 in Hopkins` GIZMO Paper
// Storing v_i*dP_i/dt in the variable unimportant, so that only a scalar is needed
// void Particles::calcdP(const double &kernelSize){
//     double fij, fji, summand;
//     for (int i = 0; i < N; i++){
//         unimportant[i] = 0;
//         for (int j = 0; j < noi[i]; j++){
//             dSqr = pow((x[i] - x[iP]), 2)
//                         + pow((y[i] - y[iP]), 2);
// #if DIM == 3
//             dSqr += pow(z[i] - z[iP], 2);
// #endif
//             r = sqrt(dSqr);
//
//             fij = 1 - 1 / m[iP]
//                     * kernelSize /omega[i] / DIM * drho[i] / (1 + kernelSize /omega[i] / DIM * dn[i]);
//
//             fji = 1 - 1 / m[i] * kernelSize / omega[iP]/DIM*drho[iP]
//                     / ( 1 + kernelSize / omega[iP] / DIM * dn[iP]);
//
//             summand = m[i]*m[nnList[j+i*MAX_NUM_INTERACTIONS]] * Kernel::dWdr(kernelSize)
//                     * (P[i]/pow(rho[i], 2)*fij
//                         +  P[iP]/pow(rho[iP], 2)*fji);
//
//
//             unimportant[i] += summand * (vx[i]*(x[i]- x[iP]) + vy[i]*(y[i]-y[iP]));
// #if DIM == 3
//             unimportant[i] += summand * vz[i]*(z[i]- z[iP]);
// #endif
//         }
//     }
// }

// For HLLC solver: Calculate norm vector between two particles:
// n_unit needs to be pre-allocated
void Particles::calcNunit(const int i, const int j, double* n_unit){
    double dSqr = pow(x[j] - x[i], 2)
                        + pow(y[j] - y[i], 2);
#if DIM == 3
    dSqr += pow(z[j] - z[i], 2);
#endif
    double r = sqrt(dSqr);
    n_unit[0] = (x[j] - x[i]) / r;
    n_unit[1] = (y[j] - y[i]) / r;
#if DIM == 3
    n_unit[2] = (z[j] - z[i]) / r;
#endif
}

#if PERIODIC_BOUNDARIES

// For HLLC solver: Calculate norm vector between two particles, one of which is a ghost:
// n_unit needs to be pre-allocated
void Particles::calcNunit(const Particles &ghostParticles, const int i, const int j, double* n_unit){
    double dSqr = pow(ghostParticles.x[j] - x[i], 2)
                        + pow(ghostParticles.y[j] - y[i], 2);
#if DIM == 3
    dSqr += pow(ghostParticles.z[j] - z[i], 2);
#endif
    double r = sqrt(dSqr);
    n_unit[0] = (ghostParticles.x[j] - x[i]) / r;
    n_unit[1] = (ghostParticles.y[j] - y[i]) / r;
#if DIM == 3
    n_unit[2] = (ghostParticles.z[j] - z[i]) / r;
#endif
}
#endif // PERIODIC_BOUNDARIES

void Particles::compDensity(){
    for(int i=0; i<N; ++i){
//#if PERIODIC_BOUNDARIES
        compOmega(i);
#if SURFACE_VOLCORR
        // Hopkins/GIZMO surface volume closure: lift kernel-sum density at
        // free surfaces by 1/FCE_i so the boundary is single-cell-wide in rho.
        rho[i] = m[i]*omega[i] / fce[i];
#else
        rho[i] = m[i]*omega[i];
#endif
        if(rho[i] <= 0.){
            if (DENSITY_FLOOR < 0){
                Logger(DEBUG) << "Negative density encountered, i = " << i << ". Aborting for debugging.";
                exit(6);
            }
            else{
                Logger(INFO) << "Negative density encountered, i = " << i << ". Using density floor!";
                rho[i] = DENSITY_FLOOR;
            }
        }

        //if ((i - 0) % 30 == 0){
        //    Logger(DEBUG) << "density from ghosts: i = " << i << ", dnst = " << rho[i];
        //}
    }
}

#if EXPLICIT_VOL_INTEGRATION
void Particles::applyExplicitVolumeOverride(){
    // Mirror of GIZMO master/hydro/density.c:982-987. The kernel-sum density
    // is the relaxation target; downstream physics sees the integrated value.
#if EXPLICIT_VOL_SEED_RHO0
    // miluphcuda INTEGRATE_DENSITY-faithful path: seed the advected density from
    // the material reference rho0 (the IC density) and use it from step 0. The
    // kernel sum is fine in the bulk but under-counts at free surfaces; for a
    // stiff EOS that deficit is multi-GPa spurious tension. See parameter.h.
    if (!explicitVolInitialized){
        for (int i = 0; i < N; ++i){
            const double r0 = MeshlessEOS->EOSGetMaterial(matId[i]).rho0;
            rhoKernel[i]   = rho[i];
            rhoExplicit[i] = r0;
            rho[i]         = r0;
        }
        explicitVolInitialized = true;
        return;
    }
    for (int i = 0; i < N; ++i){
        rhoKernel[i] = rho[i];
        if (rhoExplicit[i] > 0.) rho[i] = rhoExplicit[i];
        else                     rhoExplicit[i] = rho[i];
    }
    return;
#endif
    // Skip the first EXPLICIT_VOL_SEED_SKIP step(s): keep the working rho at the
    // kernel density so the t=0 startup-NNS truncation (~0.89) self-heals by
    // re-summation before it seeds the (re-summation-free) advected density.
    if (explicitVolSeedSkip > 0){
        for (int i = 0; i < N; ++i) rhoKernel[i] = rho[i];
        --explicitVolSeedSkip;
        return;
    }
    if (!explicitVolInitialized){
        for (int i = 0; i < N; ++i){
            rhoKernel[i]   = rho[i];
            rhoExplicit[i] = rho[i];
        }
        explicitVolInitialized = true;
    } else {
        for (int i = 0; i < N; ++i){
            rhoKernel[i] = rho[i];
            // floor-protect: if a particle got density-floored above, do not
            // poison rhoExplicit. Keep evolving the integrated value instead.
            if (rhoExplicit[i] > 0.) rho[i] = rhoExplicit[i];
            else                     rhoExplicit[i] = rho[i];
        }
    }
}

#if EXPLICIT_VOL_SPH_DIVV
double Particles::sphDivV(int i){
    // GIZMO Particle_DivVel (density.c:334,514): -sum_j dW/dr (dp.dv)/r normalized
    // by the kernel-sum neighbour number (= omega, which includes the self term).
    // The kernel normalization cancels between numerator and omega. No ghosts:
    // EXPLICIT_VOL_INTEGRATION is guarded to the non-periodic path.
    double acc = 0.;
    for (int jn = 0; jn < noi[i]; ++jn){
        const int j = nnl[i*MAX_NUM_INTERACTIONS + jn];
        const double dpx = x[i]-x[j], dpy = y[i]-y[j];
        const double dvx = vx[i]-vx[j], dvy = vy[i]-vy[j];
        double r2 = dpx*dpx + dpy*dpy;
        double dpdv = dpx*dvx + dpy*dvy;
#if DIM == 3
        const double dpz = z[i]-z[j], dvz = vz[i]-vz[j];
        r2 += dpz*dpz; dpdv += dpz*dvz;
#endif
        if (r2 <= 0.) continue;
        const double r = std::sqrt(r2);
        acc -= Kernel::WDr(r, sml[i]) * dpdv / r;
    }
    return (omega[i] > 0.) ? acc / omega[i] : 0.;
}
#endif

void Particles::integrateExplicitVolumeHalfStep(const double &dt, bool finalize){
    // Mirror of GIZMO master/kicks.c:308-318. Called twice per outer step
    // with dt = full_step / 2 to bracket the flux/updateState block.
    // No-op until the advected density has been seeded (see applyExplicitVolumeOverride):
    // during the seed-skip step(s) rho is the kernel density and must not be advected.
    if (!explicitVolInitialized) return;
    for (int i = 0; i < N; ++i){
        // Advection step: d ln rho / dt = -div v.
#if EXPLICIT_VOL_SPH_DIVV
        double divV = sphDivV(i);
#else
        double divV = vxGrad[i][0] + vyGrad[i][1];
#if DIM == 3
        divV += vzGrad[i][2];
#endif
#endif
        double arg = divV * dt;
        if (arg >  EXPLICIT_VOL_DIVV_CLAMP) arg =  EXPLICIT_VOL_DIVV_CLAMP;
        if (arg < -EXPLICIT_VOL_DIVV_CLAMP) arg = -EXPLICIT_VOL_DIVV_CLAMP;
        rhoExplicit[i] *= std::exp(-arg);

        // Relaxation toward the kernel-sum density in log-space.
        double drho2 = rhoGrad[i][0]*rhoGrad[i][0]
                     + rhoGrad[i][1]*rhoGrad[i][1];
#if DIM == 3
        drho2 += rhoGrad[i][2]*rhoGrad[i][2];
#endif
        if (drho2 > 0. && rhoExplicit[i] > 0. && rhoKernel[i] > 0.){
            double Lgrad = rho[i] / std::sqrt(drho2);
            if (Lgrad < sml[i]) Lgrad = sml[i];

#if EOS == 1 || EOS == 2
            double cEff = MeshlessEOS->EOSSoundSpeed(matId[i], rhoExplicit[i], u[i], P[i]);
#else
            double cEff = MeshlessEOS->EOSSoundSpeed(rhoExplicit[i], u[i], P[i]);
#endif
#if ELASTIC && (EOS == 1 || EOS == 2)
            // Deviatoric-wave speed sqrt(mu/rho) may set a stricter
            // restoring-force timescale (GIZMO kicks.c:313-315).
            const double mu = MeshlessEOS->EOSShearModulus(matId[i]);
            if (mu > 0.){
                const double csDev = std::sqrt(mu / rhoExplicit[i]);
                if (csDev < cEff) cEff = csDev;
            }
#endif
            const double delta = EXPLICIT_VOL_RELAX_COEF * dt * cEff / Lgrad;
            const double q0 = std::log(rhoExplicit[i]);
            const double q1 = std::log(rhoKernel[i]);
            double qn;
            if (delta > 0.005){
                const double e = std::exp(-delta);
                qn = q0 * e + q1 * (1. - e);
            } else {
                qn = q0 + (q1 - q0) * delta * (1. - 0.5 * delta);
            }
#if EXPLICIT_VOL_RELAX
            // Relaxation toward the kernel-sum density actually fires.
            rhoExplicit[i] = std::exp(qn);
#else
            // GIZMO kicks.c:317 behaviour: discard the relaxation and keep the
            // pure-advection (Monaghan continuity) density. Matches the GIZMO
            // baseline that bounces; relaxing toward the surface-noisy kernel
            // density destabilises the free surface. (qn is computed above but
            // intentionally unused, mirroring GIZMO's dead-code block.)
            (void)qn;
            rhoExplicit[i] = std::exp(q0);
#endif
        }

#if EXPLICIT_VOL_FREEZE_RHO
        // GIZMO freezes the working density through the whole hydro pass and only
        // finalises Density <- Density_ExplicitInt at end-of-step (kicks.c:395,
        // mode==1). Mid-pass (kick A) we leave rho untouched so pressure,
        // gradients, faces and Riemann reconstruction all see the single override
        // density set by applyExplicitVolumeOverride(); rhoExplicit still
        // accumulates both half-kick advections and is swapped in next step.
        if (finalize) rho[i] = rhoExplicit[i];
#else
        // Legacy: keep the working rho synchronised with the just-updated
        // integrated value after every half-kick (desyncs P and rho for the flux
        // pass, since pressure was computed before kick A).
        (void)finalize;
        rho[i] = rhoExplicit[i];
#endif
    }
}
#endif // EXPLICIT_VOL_INTEGRATION

#if VARIABLE_SML
double Particles::hMax() const {
    double hm = 0.;
    for (int i = 0; i < N; ++i){
        if (sml[i] > hm) hm = sml[i];
    }
    return hm;
}

// Helper macro: kernel volume V(h) = c_V * h^DIM and dV/dh = c_V * DIM * h^(DIM-1).
// In 2D: V = pi h^2; in 3D: V = (4/3) pi h^3.
#if DIM == 2
#  define SML_V_COEF (M_PI)
#else
#  define SML_V_COEF ((4./3.)*M_PI)
#endif

// Single Newton iteration body. Sums over the (already-built) NNL of particle
// i, and over its ghost NNL when periodic. Returns true once converged.
// real-only variant
void Particles::updateAllSmoothingLengths(){
    int totalIters = 0;
    int maxIters = 0;
    int unconverged = 0;
    int clamped = 0;
    int stuck = 0;
    double hMinObs = std::numeric_limits<double>::max();
    double hMaxObs = 0.;
    double hSum = 0.;
    for (int i = 0; i < N; ++i){
        double h = sml[i];
        bool converged = false;
        bool rescued  = false;
        int it = 0;
        for (it = 0; it < smlMaxIter; ++it){
            // self contribution at r = 0
            double n_h   = Kernel::W(0., h);
            double dn_dh = Kernel::WDh(0., h);
            for (int j = 0; j < noi[i]; ++j){
                int iP = nnl[j + i*MAX_NUM_INTERACTIONS];
                double dSqr = pow(x[i] - x[iP], 2)
                            + pow(y[i] - y[iP], 2);
#if DIM == 3
                dSqr += pow(z[i] - z[iP], 2);
#endif
                double r = sqrt(dSqr);
                n_h   += Kernel::W(r, h);
                dn_dh += Kernel::WDh(r, h);
            }
#if DIM == 2
            double V  = SML_V_COEF * h * h;
            double dV = SML_V_COEF * 2. * h;
#else
            double V  = SML_V_COEF * h * h * h;
            double dV = SML_V_COEF * 3. * h * h;
#endif
            double f      = V * n_h - smlNNNTarget;
            double fPrime = dV * n_h + V * dn_dh;
            if (fPrime == 0.){
                // Self-only case: V(h)*W(0,h) is h-invariant, so no h solves
                // f = 0. Reset to the reference kernelSize so sml is not
                // frozen at a stale value (often a previous hMax clamp);
                // the next step's NNS may pick up neighbors and a normal
                // Newton iteration can resume. The particle is in a safe,
                // known state, so it does not contribute to the bad-state
                // fraction used by the warn / panic thresholds below.
                h = smlH0;
                ++stuck;
                rescued = true;
                break;
            }
            double dh = -f / fPrime;
            // clamp the update so |dh|/h <= (SML_MAX_FACTOR - 1)
            const double maxAbs = h * (smlMaxFactor - 1.);
            if (dh >  maxAbs) dh =  maxAbs;
            if (dh < -maxAbs) dh = -maxAbs;
            h += dh;
            // hard absolute bounds: never let a single particle escape the
            // band set by MeshlessScheme (resolved from config.kernelSize).
            if (h < smlHMin) h = smlHMin;
            if (h > smlHMax) h = smlHMax;
            if (std::fabs(dh) / h < smlTol){
                converged = true;
                ++it;
                break;
            }
        }
        sml[i] = h;
        totalIters += it;
        if (it > maxIters) maxIters = it;
        if (!converged && !rescued) ++unconverged;
        const bool atBound = (h <= smlHMin) || (h >= smlHMax);
        if (atBound) ++clamped;
        if (h < hMinObs) hMinObs = h;
        if (h > hMaxObs) hMaxObs = h;
        hSum += h;
    }
    Logger(DEBUG) << "      > sml iter: avg = " << ((double)totalIters/(double)N)
                  << ", max = " << maxIters
                  << ", unconverged = " << unconverged
                  << ", clamped = " << clamped
                  << ", stuck = " << stuck
                  << ", h[min/mean/max] = " << hMinObs << " / "
                  << (hSum/(double)N) << " / " << hMaxObs;
    // Union of the two bad-state sets is bounded above by their sum, which
    // is all we need to compare against the warn / panic fractions.
    const double badFrac = (double)(unconverged + clamped) / (double)N;
    if (badFrac > smlPanicFraction){
        Logger(ERROR) << "updateAllSmoothingLengths: "
                      << (unconverged + clamped) << "/" << N
                      << " particles in a bad state (unconverged or clamped), "
                      << "fraction " << badFrac
                      << " exceeds panicFraction " << smlPanicFraction
                      << " - Aborting.";
        exit(8);
    } else if (badFrac > smlWarnFraction){
        Logger(WARN) << "updateAllSmoothingLengths: "
                     << (unconverged + clamped) << "/" << N
                     << " particles in a bad state (unconverged or clamped), "
                     << "fraction " << badFrac
                     << " exceeds warnFraction " << smlWarnFraction;
    }
}

#if PERIODIC_BOUNDARIES
// ghost-aware variant: identical to the above but also sums the ghost NNL
void Particles::updateAllSmoothingLengths(const Particles &ghostParticles){
    int totalIters = 0;
    int maxIters = 0;
    int unconverged = 0;
    int clamped = 0;
    int stuck = 0;
    double hMinObs = std::numeric_limits<double>::max();
    double hMaxObs = 0.;
    double hSum = 0.;
    for (int i = 0; i < N; ++i){
        double h = sml[i];
        bool converged = false;
        bool rescued  = false;
        int it = 0;
        for (it = 0; it < smlMaxIter; ++it){
            double n_h   = Kernel::W(0., h);
            double dn_dh = Kernel::WDh(0., h);
            for (int j = 0; j < noi[i]; ++j){
                int iP = nnl[j + i*MAX_NUM_INTERACTIONS];
                double dSqr = pow(x[i] - x[iP], 2)
                            + pow(y[i] - y[iP], 2);
#if DIM == 3
                dSqr += pow(z[i] - z[iP], 2);
#endif
                double r = sqrt(dSqr);
                n_h   += Kernel::W(r, h);
                dn_dh += Kernel::WDh(r, h);
            }
            for (int j = 0; j < noiGhosts[i]; ++j){
                int iP = nnlGhosts[j + i*MAX_NUM_GHOST_INTERACTIONS];
                double dSqr = pow(x[i] - ghostParticles.x[iP], 2)
                            + pow(y[i] - ghostParticles.y[iP], 2);
#if DIM == 3
                dSqr += pow(z[i] - ghostParticles.z[iP], 2);
#endif
                double r = sqrt(dSqr);
                n_h   += Kernel::W(r, h);
                dn_dh += Kernel::WDh(r, h);
            }
#if DIM == 2
            double V  = SML_V_COEF * h * h;
            double dV = SML_V_COEF * 2. * h;
#else
            double V  = SML_V_COEF * h * h * h;
            double dV = SML_V_COEF * 3. * h * h;
#endif
            double f      = V * n_h - smlNNNTarget;
            double fPrime = dV * n_h + V * dn_dh;
            if (fPrime == 0.){
                // See the non-ghost variant for rationale.
                h = smlH0;
                ++stuck;
                rescued = true;
                break;
            }
            double dh = -f / fPrime;
            const double maxAbs = h * (smlMaxFactor - 1.);
            if (dh >  maxAbs) dh =  maxAbs;
            if (dh < -maxAbs) dh = -maxAbs;
            h += dh;
            // hard absolute bounds: see the non-ghost variant for rationale.
            if (h < smlHMin) h = smlHMin;
            if (h > smlHMax) h = smlHMax;
            if (std::fabs(dh) / h < smlTol){
                converged = true;
                ++it;
                break;
            }
        }
        sml[i] = h;
        totalIters += it;
        if (it > maxIters) maxIters = it;
        if (!converged && !rescued) ++unconverged;
        const bool atBound = (h <= smlHMin) || (h >= smlHMax);
        if (atBound) ++clamped;
        if (h < hMinObs) hMinObs = h;
        if (h > hMaxObs) hMaxObs = h;
        hSum += h;
    }
    Logger(DEBUG) << "      > sml iter (ghosts): avg = " << ((double)totalIters/(double)N)
                  << ", max = " << maxIters
                  << ", unconverged = " << unconverged
                  << ", clamped = " << clamped
                  << ", stuck = " << stuck
                  << ", h[min/mean/max] = " << hMinObs << " / "
                  << (hSum/(double)N) << " / " << hMaxObs;
    const double badFrac = (double)(unconverged + clamped) / (double)N;
    if (badFrac > smlPanicFraction){
        Logger(ERROR) << "updateAllSmoothingLengths: "
                      << (unconverged + clamped) << "/" << N
                      << " particles in a bad state (unconverged or clamped), "
                      << "fraction " << badFrac
                      << " exceeds panicFraction " << smlPanicFraction
                      << " - Aborting.";
        exit(8);
    } else if (badFrac > smlWarnFraction){
        Logger(WARN) << "updateAllSmoothingLengths: "
                     << (unconverged + clamped) << "/" << N
                     << " particles in a bad state (unconverged or clamped), "
                     << "fraction " << badFrac
                     << " exceeds warnFraction " << smlWarnFraction;
    }
}
#endif // PERIODIC_BOUNDARIES

#undef SML_V_COEF
#endif // VARIABLE_SML

void Particles::compOmega(int i){
    const double hi = sml[i];
    double omg = 0.;
#if SURFACE_VOLCORR
    // Asymmetry of the neighbour kernel sum: S_i = sum_j W_ij * (x_i - x_j).
    // Zero in a symmetric (bulk) stencil; nonzero at a free surface where
    // half the support is empty. xi_i = |S_i| / (h_i Omega_i) is the
    // dimensionless closure asymmetry from Reinhardt & Stadel 2017.
    double sx = 0., sy = 0.;
#if DIM == 3
    double sz = 0.;
#endif
#endif
    int iP;
    for (int j=0; j<noi[i]; ++j){
        iP = nnl[j+i*MAX_NUM_INTERACTIONS];
        double dx = x[i] - x[iP];
        double dy = y[i] - y[iP];
        double dSqr = dx*dx + dy*dy;
#if DIM == 3
        double dz = z[i] - z[iP];
        dSqr += dz*dz;
#endif
        double r = sqrt(dSqr);
        double wij = kernel(r, hi);
        omg += wij;
#if SURFACE_VOLCORR
        sx += wij * dx;
        sy += wij * dy;
#if DIM == 3
        sz += wij * dz;
#endif
#endif
    }
    omega[i] = omg + kernel(0., hi); // add self interaction to normalization factor
    if (omega[i] < 0){
        Logger(WARN) << "Negative Omega encountered, i = " << i;
    }
#if SURFACE_VOLCORR
    // Self-contribution to S_i is zero (r_ii = 0), so no diagonal term.
    double S2 = sx*sx + sy*sy;
#if DIM == 3
    S2 += sz*sz;
#endif
    double xi_asym = sqrt(S2) / (hi * omega[i]);
    double fce_raw = SURFACE_VOLCORR_A - SURFACE_VOLCORR_B * xi_asym;
    if (fce_raw > 1.0) fce_raw = 1.0;
    if (fce_raw < SURFACE_VOLCORR_FLOOR) fce_raw = SURFACE_VOLCORR_FLOOR;
    fce[i] = fce_raw;
#endif
}

void Particles::compPsijTilde(Helper &helper){
    int iP;
    for (int i=0; i<N; ++i){
        const double hi = sml[i];

        // reset buffer
        for (int k=0; k<DIM*DIM; ++k){
            B[k] = 0.;
        }

        //for (int alpha = 0; alpha < DIM; ++alpha) {
        //    psijTilde_xi[i][alpha] = 0.;
        //}

        xi[0] = x[i];
        xi[1] = y[i];
#if DIM==3
        xi[2] = z[i];
#endif

        for (int j=0; j<noi[i]; ++j){
            iP = nnl[j+i*MAX_NUM_INTERACTIONS];
            double dSqr = pow(x[i] - x[iP], 2)
                          + pow(y[i] - y[iP], 2);
#if DIM == 3
            dSqr += pow(z[i] - z[iP], 2);
#endif
            double r = sqrt(dSqr);
            double psij_xi = kernel(r, hi)/omega[i];

            xj[0] = x[iP];
            xj[1] = y[iP];
#if DIM==3
            xj[2] = z[iP];
#endif

            for (int alpha=0; alpha<DIM; ++alpha){
                for(int beta=0; beta<DIM; ++beta){
                    B[DIM*alpha+beta] += (xj[alpha] - xi[alpha])*(xj[beta] - xi[beta]) * psij_xi;
                }
            }
        }

#if OUTPUT_CONDITION_NUMBER
        double normE = 0;
        for (int alpha=0; alpha<DIM; ++alpha){
           for (int beta=0; beta<DIM; ++beta){
               normE += B[alpha*DIM+beta]*B[alpha*DIM+beta];
           }
        }
#if DIM == 2
        {
            double a = B[0], b = B[1], c = B[3];
            double trace = a + c;
            double disc = sqrt((a - c) * (a - c) + 4.0 * b * b);
            lambdaMax[i] = 0.5 * (trace + disc);
            lambdaMin[i] = 0.5 * (trace - disc);
            double ev_x = b;
            double ev_y = lambdaMin[i] - a;
            double norm = sqrt(ev_x * ev_x + ev_y * ev_y);
            if (norm > 0.) {
                eigenvecMin[i][0] = ev_x / norm;
                eigenvecMin[i][1] = ev_y / norm;
            } else {
                eigenvecMin[i][0] = 1.;
                eigenvecMin[i][1] = 0.;
            }
        }
#endif
#endif

        helper.inverseMatrix(B, DIM);

#if OUTPUT_CONDITION_NUMBER
        double normB = 0;
        for (int alpha=0; alpha<DIM; ++alpha){
           for (int beta=0; beta<DIM; ++beta){
               normB += B[alpha*DIM+beta]*B[alpha*DIM+beta];
           }
        }
        conditionNumber[i] = 1./DIM * sqrt(normE*normB);
#endif

        for (int j=0; j<noi[i]; ++j) {
            iP = nnl[j+i*MAX_NUM_INTERACTIONS];
            double dSqr = pow(x[i] - x[iP], 2)
                          + pow(y[i] - y[iP], 2);
#if DIM == 3
            dSqr += pow(z[i] - z[iP], 2);
#endif
            double r = sqrt(dSqr);
            double psij_xi = kernel(r, hi) / omega[i];

            xj[0] = x[iP];
            xj[1] = y[iP];

#if DIM == 3
            xj[2] = z[iP];
#endif
            for (int alpha = 0; alpha < DIM; ++alpha) {
                psijTilde_xi[j + i * MAX_NUM_INTERACTIONS][alpha] = 0.;
                for (int beta = 0; beta < DIM; ++beta) {
                    psijTilde_xi[j + i * MAX_NUM_INTERACTIONS][alpha] += B[DIM * alpha + beta] * (xj[beta] - xi[beta]) * psij_xi;
                }
            }
        }
    }
}


void Particles::gradient(double *f, double (*grad)[DIM]){
    for (int i=0; i<N; ++i) {
        for (int alpha = 0; alpha < DIM; ++alpha) {
            grad[i][alpha] = 0;
        }

        // Hopkins 2015 MFM <-> SPH blend weight from conditionNumber[i].
        double wSPH = 0.;
        if (COND_MAX_FOR_GRADIENT > 0.){
            const double c = conditionNumber[i];
            if (COND_BLEND_HI > COND_MAX_FOR_GRADIENT){
                if      (c <= COND_MAX_FOR_GRADIENT) wSPH = 0.;
                else if (c >= COND_BLEND_HI)         wSPH = 1.;
                else wSPH = (c - COND_MAX_FOR_GRADIENT)
                          / (COND_BLEND_HI - COND_MAX_FOR_GRADIENT);
            } else {
                wSPH = (c > COND_MAX_FOR_GRADIENT) ? 1. : 0.;
            }
        }
        const double wMFM = 1. - wSPH;

        for (int j = 0; j < noi[i]; ++j) {
            int jIdx = nnl[j + i * MAX_NUM_INTERACTIONS];
            if (COND_MAX_NEIGHBOR_FOR_GRADIENT > 0. && conditionNumber[jIdx] > COND_MAX_NEIGHBOR_FOR_GRADIENT) {
                continue;
            }
            if (noi[jIdx] < MIN_NOI_FOR_NEIGHBOR_USE) {
                continue;
            }
            if (wMFM > 0.){
                for (int alpha = 0; alpha < DIM; ++alpha) {
                    grad[i][alpha] += wMFM * (f[jIdx] - f[i])
                                      * psijTilde_xi[j + i * MAX_NUM_INTERACTIONS][alpha];
                }
            }
            if (wSPH > 0.){
                const double dx = x[i] - x[jIdx];
                const double dy = y[i] - y[jIdx];
                double r2 = dx*dx + dy*dy;
#if DIM == 3
                const double dz = z[i] - z[jIdx];
                r2 += dz*dz;
#endif
                const double r = sqrt(r2);
                if (r <= 0.) continue;
                const double dW = Kernel::WDr(r, sml[i]);
                const double w  = wSPH * m[jIdx] * (f[jIdx] - f[i]) * dW / (r * rho[i]);
                grad[i][0] += w * dx;
                grad[i][1] += w * dy;
#if DIM == 3
                grad[i][2] += w * dz;
#endif
            }
        }
    }
}

void Particles::compPressure(){
    for (int i=0; i<N; ++i){
#if EOS == 1 || EOS == 2
        P[i] = MeshlessEOS->EOSPressure(matId[i], rho[i], u[i]);
#else
        P[i] = MeshlessEOS->EOSPressure(rho[i], u[i]);
#endif
//#if DIM == 3
//        P[i] = (hydro_gamma-1.)*rho[i]*(u[i]+.5*(vx[i]*vx[i]+vy[i]*vy[i]+vz[i]*vz[i]));
//#else
//        P[i] = (hydro_gamma-1.)*rho[i]*(u[i]+.5*(vx[i]*vx[i]+vy[i]*vy[i]));
//#endif
#if EOS == 0
        if(P[i] <= 0.){
            Logger(WARN) << "Zero or negative pressure @" << i;
            Logger(WARN) << "Paricle " << i << " has density " << rho[i] << " and u[i] " << u[i] << " and pressure " << P[i];
        }
#endif
    }
}

void Particles::compEffectiveFace(){
    for (int i=0; i<N; ++i){
        if(i < 10){
            Logger(DEBUG) << "Particle i = " << i << ", noi[i] = " << noi[i];
        }
        for (int j=0; j<noi[i]; ++j){
            int ji = nnl[i*MAX_NUM_INTERACTIONS+j]; // index i of particle j
            // search neighbor i in nnl[] of j
            int ij;
            for(ij=0; ij<noi[ji]; ++ij){
                if (nnl[ij+ji*MAX_NUM_INTERACTIONS] == i) break;
            }
            // V_i = m_i/rho_i, exactly GIZMO's Mass/Density in face assembly
            // (hydro_core_meshless.h:46, hydro_evaluate.h:80,311). Deriving the
            // face volume from rho -- rather than the kernel value fce/omega --
            // keeps the geometry consistent with whatever density the rest of
            // the hydro pass sees:
            //   SURFACE_VOLCORR on, explicit off : rho = m*omega/fce
            //                                       -> m/rho = fce/omega (unchanged)
            //   SURFACE_VOLCORR off              : rho = m*omega
            //                                       -> m/rho = 1/omega   (unchanged)
            //   EXPLICIT_VOL_INTEGRATION on      : rho = rhoExplicit
            //                                       -> m/rho = m/rhoExplicit
            // The last case is the fix: with the explicit override the faces
            // must track rhoExplicit too, or the cell carries one density for
            // pressure and another for its volume, which injects a spurious
            // net outward flux (the "bloat"). GIZMO avoids this by overwriting
            // Density with Density_ExplicitInt (density.c:985) before building
            // the faces from Mass/Density.
            const double Vi  = m[i]  / rho[i];
            const double Vji = m[ji] / rho[ji];
            for (int alpha=0; alpha<DIM; ++alpha){
                Aij[i*MAX_NUM_INTERACTIONS+j][alpha] = Vi*psijTilde_xi[i*MAX_NUM_INTERACTIONS+j][alpha]
                        - Vji*psijTilde_xi[ij+ji*MAX_NUM_INTERACTIONS][alpha];
            }

            //if(i < 10){
            //    Logger(DEBUG) << "A[" << i << " -> " << ji << "] = [" << Aij[i*MAX_NUM_INTERACTIONS+j][0] << ", "
            //                  << Aij[i*MAX_NUM_INTERACTIONS+j][1] << ", " << Aij[i*MAX_NUM_INTERACTIONS+j][2] << "], xi = ["
            //                  << x[i] << ", " << y[i] << ", " << z[i] << "], xj = " << x[ji] << ", " << y[ji] << ", " << z[ji] << "]";
            //}
        }
    }
}

void Particles::slopeLimiter(Particles *ghostParticles){
    double *rhoGhost = nullptr, *vxGhost = nullptr, *vyGhost = nullptr,
#if DIM==3
        *vzGhost = nullptr,
#endif
    *PGhost = nullptr;
#if PERIODIC_BOUNDARIES
    rhoGhost = ghostParticles->rho;
    vxGhost = ghostParticles->vx;
    vyGhost = ghostParticles->vy;
#if DIM==3
    vzGhost = ghostParticles->vz;
#endif
    PGhost = ghostParticles->P;
#endif
    slopeLimiter(rho, rhoGrad, ghostParticles, rhoGhost);
    slopeLimiter(vx, vxGrad, ghostParticles, vxGhost);
    slopeLimiter(vy, vyGrad, ghostParticles, vyGhost);
#if DIM==3
    slopeLimiter(vz, vzGrad, ghostParticles, vzGhost);
#endif
    slopeLimiter(P, PGrad, ghostParticles, PGhost);
}

#if ELASTIC
void Particles::slopeLimiterStress(Particles *ghostParticles){
    double *SxxGhost = nullptr, *SxyGhost = nullptr, *SyyGhost = nullptr;
#if DIM==3
    double *SxzGhost = nullptr, *SyzGhost = nullptr, *SzzGhost = nullptr;
#endif
#if PERIODIC_BOUNDARIES
    SxxGhost = ghostParticles->Sxx;
    SxyGhost = ghostParticles->Sxy;
    SyyGhost = ghostParticles->Syy;
#if DIM==3
    SxzGhost = ghostParticles->Sxz;
    SyzGhost = ghostParticles->Syz;
    SzzGhost = ghostParticles->Szz;
#endif
#endif
    slopeLimiter(Sxx, SxxGrad, ghostParticles, SxxGhost);
    slopeLimiter(Sxy, SxyGrad, ghostParticles, SxyGhost);
    slopeLimiter(Syy, SyyGrad, ghostParticles, SyyGhost);
#if DIM==3
    slopeLimiter(Sxz, SxzGrad, ghostParticles, SxzGhost);
    slopeLimiter(Syz, SyzGrad, ghostParticles, SyzGhost);
    slopeLimiter(Szz, SzzGrad, ghostParticles, SzzGhost);
#endif
}
#endif // ELASTIC
void Particles::slopeLimiter(double *f, double (*grad)[DIM],
                             Particles *ghostParticles, double *fGhost){
    double xij[DIM], xijxi[DIM];
    //Logger(DEBUG) << "#################### NEXT GRADIENT ####################";

    for (int i=0; i<N; ++i){
        double psiMaxNgb { std::numeric_limits<double>::lowest() };
        double psiMinNgb { std::numeric_limits<double>::max() };
        double psiMaxMid { std::numeric_limits<double>::lowest() };
        double psiMinMid { std::numeric_limits<double>::max() };

        for (int jn=0; jn<noi[i]; ++jn) {
            int j = nnl[i * MAX_NUM_INTERACTIONS + jn];

            //xij[0] = x[i] + sml[i] / 2. * (x[j] - x[i]);
            //xij[1] = y[i] + sml[i] / 2. * (y[j] - y[i]);
#if FIRST_ORDER_QUAD_POINT
            xij[0] = (x[i] + x[j])/2.;
            xij[1] = (y[i] + y[j])/2.;
#else
            xij[0] = x[i] + sml[i] / 4. * (x[j] - x[i]);
            xij[1] = y[i] + sml[i] / 4. * (y[j] - y[i]);
#endif

            xijxi[0] = xij[0] - x[i];
            xijxi[1] = xij[1] - y[i];
#if DIM == 3
#if FIRST_ORDER_QUAD_POINT
            xij[2] = (z[i] + z[j])/2.;
#else
            //xij[2] = z[i] + sml[i]/2. * (z[j] - z[i]);
            xij[2] = z[i] + sml[i]/4. * (z[j] - z[i]);
#endif
            xijxi[2] = xij[2] - z[i];
#endif
            if (psiMaxNgb < f[j]) psiMaxNgb = f[j];
            if (psiMinNgb > f[j]) psiMinNgb = f[j];
            // reconstruct states at effective face
            double fij = f[i] + Helper::dotProduct(grad[i], xijxi);
            if (psiMaxMid < fij) psiMaxMid = fij;
            if (psiMinMid > fij) psiMinMid = fij;
        }
#if PERIODIC_BOUNDARIES
        for (int jn=0; jn<noiGhosts[i]; ++jn) {
            int j = nnlGhosts[i * MAX_NUM_GHOST_INTERACTIONS + jn];

            //xij[0] = x[i] + sml[i] / 2. * (ghostParticles->x[j] - x[i]);
            //xij[1] = y[i] + sml[i] / 2. * (ghostParticles->y[j] - y[i]);

#if FIRST_ORDER_QUAD_POINT
            xij[0] = (x[i] + ghostParticles->x[j])/2.;
            xij[1] = (y[i] + ghostParticles->y[j])/2.;
#else
            xij[0] = x[i] + sml[i] / 4. * (ghostParticles->x[j] - x[i]);
            xij[1] = y[i] + sml[i] / 4. * (ghostParticles->y[j] - y[i]);
#endif

            xijxi[0] = xij[0] - x[i];
            xijxi[1] = xij[1] - y[i];
#if DIM == 3
#if FIRST_ORDER_QUAD_POINT
            xij[2] = (z[i] + ghostParticles->z[j])/2.;
#else
            //xij[2] = z[i] + sml[i]/2. * (ghostParticles->z[j] - z[i]);
            xij[2] = z[i] + sml[i]/4. * (ghostParticles->z[j] - z[i]);
#endif

            xijxi[2] = xij[2] - z[i];
#endif
            if (psiMaxNgb < fGhost[j]) psiMaxNgb = fGhost[j];
            if (psiMinNgb > fGhost[j]) psiMinNgb = fGhost[j];
            // reconstruct states at effective face
            double fij = f[i] + Helper::dotProduct(grad[i], xijxi);
            if (psiMaxMid < fij) psiMaxMid = fij;
            if (psiMinMid > fij) psiMinMid = fij;
        }
#endif
        // update gradients if necessary
        double denomMax = psiMaxMid - f[i];
        double denomMin = f[i] - psiMinMid;
        double alphaMax = (std::fabs(denomMax) > 1e-30) ? std::fabs((psiMaxNgb - f[i]) / denomMax) : 1.0;
        double alphaMin = (std::fabs(denomMin) > 1e-30) ? std::fabs((f[i] - psiMinNgb) / denomMin) : 1.0;

        //if(i==6){
        //    Logger(DEBUG) << "alphaMin = " << alphaMin << ", alphaMax = " << alphaMax
        //                << ", psiMinNgb = " << psiMinNgb << ", psiMaxNgb = " << psiMaxNgb
        //                << ", psiMinNgb = " << psiMinMid << ", psiMaxNgb = " << psiMaxMid
        //                << ", f[i] = " << f[i];
        //}

#if DEBUG_LVL
        if(i == 184){
            Logger(DEBUG) << "slopeLimiter i=2: f[i]=" << f[i]
                          << " grad=[" << grad[i][0] << ", " << grad[i][1] << "]"
                          << " alphaMin=" << alphaMin << " alphaMax=" << alphaMax
                          << " psiMinNgb=" << psiMinNgb << " psiMaxNgb=" << psiMaxNgb
                          << " psiMinMid=" << psiMinMid << " psiMaxMid=" << psiMaxMid;
        }
#endif
        if (alphaMin <= alphaMax && BETA*alphaMin < 1.){
            //Logger(DEBUG) << "        > Limiting gradient with alphaMin@" << i << ", alpha = "
            //              << alphaMin << ", grad[0] = " << grad[i][0] << ", grad[1] = " << grad[i][1];
            grad[i][0] *= alphaMin;
            grad[i][1] *= alphaMin;
#if DIM==3
            grad[i][2] *= alphaMin;
#endif
        } else if (alphaMax <= alphaMin && BETA*alphaMax < 1.){
            //Logger(DEBUG) << "        > Limiting gradient with alphaMax@" << i << ", alpha = "
            //            << alphaMax << ", grad[0] = " << grad[i][0] << ", grad[1] = " << grad[i][1];
            grad[i][0] *= alphaMax;
            grad[i][1] *= alphaMax;
#if DIM==3
            grad[i][2] *= alphaMax;
#endif
        }
    }
}

double Particles::compGlobalTimestep(){
    double dt_ = std::numeric_limits<double>::max();
    for (int i=0; i<N; ++i){
        double vSig = std::numeric_limits<double>::min();
#if EOS == 1 || EOS == 2
        double ci = MeshlessEOS->EOSSoundSpeed(matId[i], rho[i], -1, P[i]);
#else
        double ci = MeshlessEOS->EOSSoundSpeed(rho[i], -1, P[i]);
#endif
        // searching for maximum signal speed
        for (int jn=0; jn<noi[i]; ++jn){
            int j = nnl[i*MAX_NUM_INTERACTIONS+jn];

#if EOS == 1 || EOS == 2
            double cj = MeshlessEOS->EOSSoundSpeed(matId[j], rho[j], -1, P[j]);
#else
            double cj = MeshlessEOS->EOSSoundSpeed(rho[j], -1, P[j]);
#endif

            double xij[DIM], vij[DIM];

            xij[0] = x[i] - x[j];
            xij[1] = y[i] - y[j];

            vij[0] = vx[i] - vx[j];
            vij[1] = vy[i] - vy[j];

#if DIM == 3
            xij[2] = z[i] - z[j];
            vij[2] = vz[i] - vz[j];
#endif
            double vijxij = Helper::dotProduct(vij, xij)/sqrt(Helper::dotProduct(xij, xij));
            vijxij = vijxij < 0. ? vijxij : 0.;

            double vSig_i = ci+cj-vijxij;
            vSig = vSig_i > vSig ? vSig_i : vSig;

        }

        // TODO: Note: Kernel size is double the actual kernel size
        double dt = CFL*sml[i]/vSig;
        dt_ = dt < dt_ ? dt : dt_;
    }

    return dt_;
}

#if ELASTIC
double Particles::compElasticTimestep(){
    double dt_ = std::numeric_limits<double>::max();
    for (int i=0; i<N; ++i){
        double vSig = std::numeric_limits<double>::min();

        // Elastic longitudinal wave speed (eq. 92)
#if EOS == 1 || EOS == 2
        double Ki  = MeshlessEOS->EOSBulkModulus(matId[i], rho[i], P[i]);
        double mui = MeshlessEOS->EOSShearModulus(matId[i]);
#else
        double Ki  = MeshlessEOS->EOSBulkModulus(rho[i], P[i]);
        double mui = 0.;
#endif
        double ci = sqrt((Ki + 4.0/3.0 * mui) / rho[i]);

        for (int jn=0; jn<noi[i]; ++jn){
            int j = nnl[i*MAX_NUM_INTERACTIONS+jn];

#if EOS == 1 || EOS == 2
            double Kj  = MeshlessEOS->EOSBulkModulus(matId[j], rho[j], P[j]);
            double muj = MeshlessEOS->EOSShearModulus(matId[j]);
#else
            double Kj  = MeshlessEOS->EOSBulkModulus(rho[j], P[j]);
            double muj = 0.;
#endif
            double cj = sqrt((Kj + 4.0/3.0 * muj) / rho[j]);

            double xij[DIM], vij[DIM];
            xij[0] = x[i] - x[j];
            xij[1] = y[i] - y[j];
            vij[0] = vx[i] - vx[j];
            vij[1] = vy[i] - vy[j];
#if DIM == 3
            xij[2] = z[i] - z[j];
            vij[2] = vz[i] - vz[j];
#endif
            double vijxij = Helper::dotProduct(vij, xij)/sqrt(Helper::dotProduct(xij, xij));
            vijxij = vijxij < 0. ? vijxij : 0.;

            double vSig_i = ci + cj - vijxij;
            vSig = vSig_i > vSig ? vSig_i : vSig;
        }

        double dt = CFL*sml[i]/vSig;
        dt_ = dt < dt_ ? dt : dt_;
    }
    return dt_;
}
#endif // ELASTIC

void Particles::compRiemannStatesLR(const double &dt){
    // Rate-limited diagnostics for free-surface reconstruction blow-ups.
    // Triggers when the final reconstructed rho falls below DENSITY_FLOOR
    // or the reconstructed velocity magnitude exceeds EXTREME_V.
    int nExtremeLogged = 0;
    const int MAX_EXTREME_LOGS = 5;
    const double EXTREME_V = 10.;

    for (int i=0; i<N; ++i){
        double xij[DIM];
        //double vFrame[DIM];
        // helper vectors
        double xijxi[DIM], xjxi[DIM], xijxj[DIM];

        for (int jn=0; jn<noi[i]; ++jn){
            int j = nnl[i*MAX_NUM_INTERACTIONS+jn];

            xjxi[0] = x[j] - x[i];
            xjxi[1] = y[j] - y[i];

            //xij[0] = x[i] + sml[i]/2. * xjxi[0];
            //xij[1] = y[i] + sml[i]/2. * xjxi[1];

#if FIRST_ORDER_QUAD_POINT
            xij[0] = (x[i] + x[j])/2.;
            xij[1] = (y[i] + y[j])/2.;

            xijxj[0] = .5*(x[i] - x[j]);
            xijxj[1] = .5*(y[i] - y[j]);

            xijxi[0] = .5*(x[j] - x[i]);
            xijxi[1] = .5*(y[j] - y[i]);

#else // !FIRST_ORDER_QUAD_POINT
            xij[0] = x[i] + sml[i]/4. * xjxi[0];
            xij[1] = y[i] + sml[i]/4. * xjxi[1];

            xijxi[0] = xij[0] - x[i];
            xijxi[1] = xij[1] - y[i];

            xijxj[0] = xij[0] - x[j];
            xijxj[1] = xij[1] - y[j];
#endif // FIRST_ORDER_QUAD_POINT

#if DIM==3
#if FIRST_ORDER_QUAD_POINT
            xij[2] = (z[i] + z[j])/2.;
            xijxj[2] = .5*(z[i] - z[j]);
            xijxi[2] = .5*(z[j] - z[i]);
#else // !FIRST_ORDER_QUAD_POINT
            xjxi[2] = z[j] - z[i];
            //xij[2] = z[i] + sml[i]/2. * xjxi[2];

            xij[2] = z[i] + sml[i]/4. * xjxi[2];

            xijxi[2] = xij[2] - z[i];
            xijxj[2] = xij[2] - z[j];
#endif // FIRST_ORDER_QUAD_POINT
#endif // DIM == 3

            int iW = i*MAX_NUM_INTERACTIONS+jn;
#if !MOVE_PARTICLES
            vFrame[iW][0] = 0.;
            vFrame[iW][1] = 0.;
#if DIM==3
            vFrame[iW][2] = 0.;
#endif // DIM == 3
#else // MOVE_PARTICLES
#if FIRST_ORDER_QUAD_POINT
            vFrame[iW][0] = (vx[i] + vx[j])/2.;
            vFrame[iW][1] = (vy[i] + vy[j])/2.;
#if DIM==3
            vFrame[iW][2] = (vz[i] + vz[j])/2.;
#endif // DIM == 3
#else // !FIRST_ORDER_QUAD_POINT
            double dotProd = Helper::dotProduct(xijxi, xjxi);
            double dSqr = Helper::dotProduct(xjxi, xjxi);

            vFrame[iW][0] = vx[i] + (vx[j]-vx[i]) * dotProd/dSqr;
            vFrame[iW][1] = vy[i] + (vy[j]-vy[i]) * dotProd/dSqr;
#if DIM==3
            vFrame[iW][2] = vz[i] + (vz[j]-vz[i]) * dotProd/dSqr;
#endif // DIM == 3
#endif // FIRST_ORDER_QUAD_POINT
#endif // MOVE_PARTICLES
            // boost frame to effective face
#if DEBUG_LVL
            if(i == 184 && jn == 0){
                Logger(DEBUG) << "compRiemannStatesLR i=2, j=" << j
                              << " RAW: vx[i]=" << vx[i] << ", vy[i]=" << vy[i]
                              << ", vFrame=[" << vFrame[iW][0] << ", " << vFrame[iW][1] << "]";
                Logger(DEBUG) << "  vxGrad[i]=[" << vxGrad[i][0] << ", " << vxGrad[i][1]
                              << "], vyGrad[i]=[" << vyGrad[i][0] << ", " << vyGrad[i][1] << "]";
                Logger(DEBUG) << "  PGrad[i]=[" << PGrad[i][0] << ", " << PGrad[i][1]
                              << "], rho[i]=" << rho[i];
                Logger(DEBUG) << "  xijxi=[" << xijxi[0] << ", " << xijxi[1] << "]";
#if ELASTIC
                Logger(DEBUG) << "  SxxGrad[i]=[" << SxxGrad[i][0] << ", " << SxxGrad[i][1]
                              << "], SxyGrad[i]=[" << SxyGrad[i][0] << ", " << SxyGrad[i][1]
                              << "], SyyGrad[i]=[" << SyyGrad[i][0] << ", " << SyyGrad[i][1] << "]";
#endif
            }
#endif
            WijR[iW][0] = rho[i];
            WijL[iW][0] = rho[j];
            WijR[iW][1] = P[i];
            WijL[iW][1] = P[j];

            WijR[iW][2] = vx[i] - vFrame[iW][0];
            WijL[iW][2] = vx[j] - vFrame[iW][0];
            WijR[iW][3] = vy[i] - vFrame[iW][1];
            WijL[iW][3] = vy[j] - vFrame[iW][1];

#if DIM == 3
            WijR[iW][4] = vz[i] - vFrame[iW][2];
            WijL[iW][4] = vz[j] - vFrame[iW][2];
#endif // DIM == 3

#if PAIRWISE_LIMITER
            double WijR_buf[DIM+2], WijL_buf[DIM+2];
            for(int nu=0; nu<DIM+2; ++nu){
                WijR_buf[nu] = WijR[iW][nu];
                WijL_buf[nu] = WijL[iW][nu];
            }
#endif

            //if (i == 46){// && jn == 28){
            //    Logger(DEBUG) << "        j = " << j
            //                  << ", rhoL = " << WijL[iW][0] << ", rhoR = " << WijR[iW][0]
            //                  << ", uL = " << WijL[iW][2] << ", uR = " << WijR[iW][2]
            //                  << ", PL = " << WijL[iW][1] << ", PR = " << WijR[iW][1];
            //}

            //if(i == 46){// && jn == 28){
            //    Logger(DEBUG) << "vFrame[iW] = [" << vFrame[iW][0]
            //                << ", " << vFrame[iW][1] << "]"
            //                << ", rhoGrad[i] = [" << rhoGrad[i][0]
            //                << ", " << rhoGrad[i][1] << "]";
            //    Logger(DEBUG) << "PGrad[i] = [" << PGrad[i][0] << ", " << PGrad[i][1]
            //                << "], PGrad[j] = [" << PGrad[j][0] << ", " << PGrad[j][1]
            //                << "], rho[i] = " << rho[i]
                //            << "], xijxi = [" << xijxi[0] << ", " << xijxi[1]
                //            << "], xijxj = [" << xijxj[0] << ", " << xijxj[1] << "]";
            //                << ", xj = [" << x[j] << ", " << y[j] << "] @" << j;
                //exit(5);
            //}
            // Raw (pre-extrap) snapshot for diagnostic logging.
            double WijR_raw[DIM+2], WijL_raw[DIM+2];
            for (int nu = 0; nu < DIM+2; ++nu){
                WijR_raw[nu] = WijR[iW][nu];
                WijL_raw[nu] = WijL[iW][nu];
            }

            // reconstruction at effective face
            WijR[iW][0] += Helper::dotProduct(rhoGrad[i], xijxi);
            WijL[iW][0] += Helper::dotProduct(rhoGrad[j], xijxj);
            WijR[iW][1] += Helper::dotProduct(PGrad[i], xijxi);
            WijL[iW][1] += Helper::dotProduct(PGrad[j], xijxj);
            WijR[iW][2] += Helper::dotProduct(vxGrad[i], xijxi);
            WijL[iW][2] += Helper::dotProduct(vxGrad[j], xijxj);
            WijR[iW][3] += Helper::dotProduct(vyGrad[i], xijxi);
            WijL[iW][3] += Helper::dotProduct(vyGrad[j], xijxj);
#if DIM==3
            WijR[iW][4] += Helper::dotProduct(vzGrad[i], xijxi);
            WijL[iW][4] += Helper::dotProduct(vzGrad[j], xijxj);
#endif // DIM == 3

            // Post-grad-extrap snapshot for diagnostics.
            double WijR_extrap[DIM+2], WijL_extrap[DIM+2];
            for (int nu = 0; nu < DIM+2; ++nu){
                WijR_extrap[nu] = WijR[iW][nu];
                WijL_extrap[nu] = WijL[iW][nu];
            }

#if DEBUG_LVL
            if(i == 184 && jn == 0){
                Logger(DEBUG) << "  AFTER GRAD EXTRAP: vxR=" << WijR[iW][2]
                              << ", vyR=" << WijR[iW][3];
            }
#endif

#if PAIRWISE_LIMITER
            double xijxi_abs = 0., xijxj_abs = 0., xjxi_abs = 0.;
            for(int alpha=0; alpha<DIM; ++alpha){
                xijxi_abs += xijxi[alpha]*xijxi[alpha];
                xijxj_abs += xijxj[alpha]*xijxj[alpha];
                xjxi_abs += xjxi[alpha]*xjxi[alpha];
            }
            xijxi_abs = sqrt(xijxi_abs);
            xijxj_abs = sqrt(xijxj_abs);
            xjxi_abs = sqrt(xjxi_abs);

            for(int nu=0; nu<DIM+2; ++nu){
                WijR[iW][nu] = pairwiseLimiter(WijR[iW][nu], WijR_buf[nu], WijL_buf[nu], xijxi_abs, xjxi_abs);
                WijL[iW][nu] = pairwiseLimiter(WijL[iW][nu], WijL_buf[nu], WijR_buf[nu], xijxj_abs, xjxi_abs);
            }
#endif // PAIRWISE_LIMITER

            // Post-pairwise-limiter snapshot for diagnostics.
            double WijR_pwl[DIM+2], WijL_pwl[DIM+2];
            for (int nu = 0; nu < DIM+2; ++nu){
                WijR_pwl[nu] = WijR[iW][nu];
                WijL_pwl[nu] = WijL[iW][nu];
            }

#if DEBUG_LVL
            if(i == 184 && jn == 0){
                Logger(DEBUG) << "  AFTER PAIRWISE LIM: vxR=" << WijR[iW][2]
                              << ", vyR=" << WijR[iW][3];
            }
#if EOS == 0
            if(WijR[iW][1] < 0. || WijL[iW][1] < 0.) {
                Logger(WARN) << "FACE RECONSTRUCTION > Negative pressure encountered@(i = " << i << ", j = " << j << ")";
                Logger(DEBUG) << "rhoL = " << WijL[iW][0] << ", rhoR = " << WijR[iW][0]
                              << ", uL = " << WijL[iW][2] << ", uR = " << WijR[iW][2]
                              << ", PL = " << WijL[iW][1] << ", PR = " << WijR[iW][1]
                              << ", PGrad[i] = [" << PGrad[i][0] << ", " << PGrad[i][1]
#if DIM == 3
                              << ", " << PGrad[i][2]
#endif // DIM == 3
                              << "] , PGrad[j] = [" << PGrad[j][0] << ", " << PGrad[j][1]
#if DIM == 3
                              << ", " << PGrad[j][2]
#endif // DIM == 3
                              << "]";
            }
#endif // EOS == 0
#endif // DEBUG_LVL
            //if (i == 46){
            //    Logger(DEBUG) << "        j = " << j
            //                  << ", rhoL = " << WijL[iW][0] << ", rhoR = " << WijR[iW][0]
            //                  << ", uL = " << WijL[iW][2] << ", uR = " << WijR[iW][2]
            //                  << ", PL = " << WijL[iW][1] << ", PR = " << WijR[iW][1];
            //}

            // predict half a timestep
            double viDiv = vxGrad[i][0] + vyGrad[i][1];
            double vjDiv = vxGrad[j][0] + vyGrad[j][1];
#if DIM==3
            viDiv += vzGrad[i][2];
            vjDiv += vzGrad[j][2];
#endif // DIM == 3
            // density
            WijR[iW][0] -= dt/2. * (rho[i] * viDiv + (vx[i]-vFrame[iW][0])*rhoGrad[i][0] + (vy[i]-vFrame[iW][1])*rhoGrad[i][1]);
            WijL[iW][0] -= dt/2. * (rho[j] * vjDiv + (vx[j]-vFrame[iW][0])*rhoGrad[j][0] + (vy[j]-vFrame[iW][1])*rhoGrad[j][1]);

            // energy
            // WijR[iW][1] -= dt/2. * (hydro_gamma*P[i] * viDiv + (vx[i]-vFrame[iW][0])*PGrad[i][0] + (vy[i]-vFrame[iW][1])*PGrad[i][1]);
            // WijL[iW][1] -= dt/2. * (hydro_gamma*P[j] * vjDiv + (vx[j]-vFrame[iW][0])*PGrad[j][0] + (vy[j]-vFrame[iW][1])*PGrad[j][1]);
#if EOS == 1 || EOS == 2
            double Ki = MeshlessEOS->EOSBulkModulus(matId[i], rho[i], P[i]);
            double Kj = MeshlessEOS->EOSBulkModulus(matId[j], rho[j], P[j]);
#else
            double Ki = MeshlessEOS->EOSBulkModulus(rho[i], P[i]);
            double Kj = MeshlessEOS->EOSBulkModulus(rho[j], P[j]);
#endif
            WijR[iW][1] -= dt/2. * (Ki * viDiv
                                        + (vx[i]-vFrame[iW][0])*PGrad[i][0] + (vy[i]-vFrame[iW][1])*PGrad[i][1]);
            WijL[iW][1] -= dt/2. * (Kj * vjDiv
                                        + (vx[j]-vFrame[iW][0])*PGrad[j][0] + (vy[j]-vFrame[iW][1])*PGrad[j][1]);

            // velocities
            // TODO: center vL and vR and update vFrame (compare to GIZMO code hydro_core_meshless.h:178ff) ??
            WijR[iW][2] -= dt/2. * (PGrad[i][0]/rho[i] + (vx[i]-vFrame[iW][0])*vxGrad[i][0] + (vy[i]-vFrame[iW][1])*vxGrad[i][1]);
            WijL[iW][2] -= dt/2. * (PGrad[j][0]/rho[j] + (vx[j]-vFrame[iW][0])*vxGrad[j][0] + (vy[j]-vFrame[iW][1])*vxGrad[j][1]);
            WijR[iW][3] -= dt/2. * (PGrad[i][1]/rho[i] + (vx[i]-vFrame[iW][0])*vyGrad[i][0] + (vy[i]-vFrame[iW][1])*vyGrad[i][1]);
            WijL[iW][3] -= dt/2. * (PGrad[j][1]/rho[j] + (vx[j]-vFrame[iW][0])*vyGrad[j][0] + (vy[j]-vFrame[iW][1])*vyGrad[j][1]);
#if ELASTIC
            // Elastic source: dv/dt += (1/ρ) ∇·S
            WijR[iW][2] += dt/2. * (SxxGrad[i][0] + SxyGrad[i][1]) / rho[i];
            WijL[iW][2] += dt/2. * (SxxGrad[j][0] + SxyGrad[j][1]) / rho[j];
            WijR[iW][3] += dt/2. * (SxyGrad[i][0] + SyyGrad[i][1]) / rho[i];
            WijL[iW][3] += dt/2. * (SxyGrad[j][0] + SyyGrad[j][1]) / rho[j];
#endif // ELASTIC
#if DEBUG_LVL
            if(i == 184 && jn == 0){
                Logger(DEBUG) << "  AFTER TIME PRED+ELASTIC: vxR=" << WijR[iW][2]
                              << ", vyR=" << WijR[iW][3];
            }
#endif
#if DIM==3
            // density
            WijR[iW][0] -= dt/2. * (vz[i]-vFrame[iW][2])*rhoGrad[i][2];
            WijL[iW][0] -= dt/2. * (vz[j]-vFrame[iW][2])*rhoGrad[j][2];

            // energy
            WijR[iW][1] -= dt/2. * (vz[i]-vFrame[iW][2])*PGrad[i][2];
            WijL[iW][1] -= dt/2. * (vz[j]-vFrame[iW][2])*PGrad[j][2];

#if DEBUG_LVL
#if EOS == 0
            if(WijR[iW][1] < 0. || WijL[iW][1] < 0.) {
                Logger(WARN) << "TIME PREDICTION > Negative pressure encountered@(i = " << i << ", j = " << j << ")";

                double timePredP_i = dt / 2. * (MeshlessEOS->EOSEnergyFluxGamma(rho[i], P[i], u[i]) * viDiv + (vx[i] - vFrame[iW][0]) * PGrad[i][0] +
                                                (vy[i] - vFrame[iW][1]) * PGrad[i][1]
#if DIM == 3
                                                + (vz[i] - vFrame[iW][2]) * PGrad[i][2]
#endif // DIM == 3
                );
                Logger(DEBUG) << "Pressure timestep prediction term @i: " << timePredP_i;

                double timePredP_j = dt / 2. * (MeshlessEOS->EOSEnergyFluxGamma(rho[j], P[j], u[j]) * viDiv + (vx[j] - vFrame[iW][0]) * PGrad[j][0] +
                                                (vy[j] - vFrame[iW][1]) * PGrad[j][1]
#if DIM == 3
                                                + (vz[j] - vFrame[iW][2]) * PGrad[j][2]
#endif // DIM == 3
                );
                Logger(DEBUG) << "Pressure timestep prediction term @j: " << timePredP_j;
            }
#endif // EOS == 0
#endif // DEBUG_LVL

            // velocities
            WijR[iW][2] -= dt/2. * (vz[i]-vFrame[iW][2])*vxGrad[i][2];
            WijL[iW][2] -= dt/2. * (vz[i]-vFrame[iW][2])*vxGrad[j][2];
            WijR[iW][3] -= dt/2. * (vz[i]-vFrame[iW][2])*vyGrad[i][2];
            WijL[iW][3] -= dt/2. * (vz[i]-vFrame[iW][2])*vyGrad[j][2];
#if DIM == 3
            WijR[iW][4] -= dt/2. * (PGrad[i][2]/rho[i] + (vx[i]-vFrame[iW][0])*vzGrad[i][0] + (vy[i]-vFrame[iW][1])*vzGrad[i][1] + (vz[i]-vFrame[iW][2])*vzGrad[i][2]);
            WijL[iW][4] -= dt/2. * (PGrad[j][2]/rho[j] + (vx[j]-vFrame[iW][0])*vzGrad[j][0] + (vy[j]-vFrame[iW][1])*vzGrad[j][1] + (vz[j]-vFrame[iW][2])*vzGrad[j][2]);
#endif // DIM == 3
#if ELASTIC
            // Elastic 3D source terms for velocity prediction
            WijR[iW][2] += dt/2. * SxzGrad[i][2] / rho[i];
            WijL[iW][2] += dt/2. * SxzGrad[j][2] / rho[j];
            WijR[iW][3] += dt/2. * SyzGrad[i][2] / rho[i];
            WijL[iW][3] += dt/2. * SyzGrad[j][2] / rho[j];
            WijR[iW][4] += dt/2. * (SxzGrad[i][0] + SyzGrad[i][1] + SzzGrad[i][2]) / rho[i];
            WijL[iW][4] += dt/2. * (SxzGrad[j][0] + SyzGrad[j][1] + SzzGrad[j][2]) / rho[j];
#endif // ELASTIC
#endif // DIM == 3

            // Free-surface diagnostic: dump the per-stage breakdown whenever
            // the final reconstructed rho falls below DENSITY_FLOOR or a
            // reconstructed velocity magnitude exceeds EXTREME_V. Rate-limited
            // to MAX_EXTREME_LOGS pairs per call so the log stays manageable.
            {
                bool pathologicalR = (WijR[iW][0] < DENSITY_FLOOR) ||
                                     (fabs(WijR[iW][2]) > EXTREME_V) ||
                                     (fabs(WijR[iW][3]) > EXTREME_V);
                bool pathologicalL = (WijL[iW][0] < DENSITY_FLOOR) ||
                                     (fabs(WijL[iW][2]) > EXTREME_V) ||
                                     (fabs(WijL[iW][3]) > EXTREME_V);
                if ((pathologicalR || pathologicalL) && nExtremeLogged < MAX_EXTREME_LOGS){
                    ++nExtremeLogged;
                    Logger(DEBUG) << "EXTREME reconstruction @ pair (i=" << i
                                  << ", j=" << j << "): xi=[" << x[i] << "," << y[i]
                                  << "] xj=[" << x[j] << "," << y[j]
                                  << "] noi[i]=" << noi[i] << " noi[j]=" << noi[j]
#if OUTPUT_CONDITION_NUMBER
                                  << " cond[i]=" << conditionNumber[i]
                                  << " cond[j]=" << conditionNumber[j]
#endif
                                  ;
                    Logger(DEBUG) << "  |xijxi|=" << sqrt(xijxi[0]*xijxi[0]+xijxi[1]*xijxi[1])
                                  << "  |xijxj|=" << sqrt(xijxj[0]*xijxj[0]+xijxj[1]*xijxj[1])
                                  << "  sml[i]=" << sml[i] << "  sml[j]=" << sml[j];
                    Logger(DEBUG) << "  rho[i]=" << rho[i] << "  rho[j]=" << rho[j]
                                  << "  P[i]="  << P[i]   << "  P[j]="  << P[j];
                    Logger(DEBUG) << "  rhoGrad[i]=[" << rhoGrad[i][0] << "," << rhoGrad[i][1]
                                  << "]  rhoGrad[j]=[" << rhoGrad[j][0] << "," << rhoGrad[j][1] << "]";
                    Logger(DEBUG) << "  PGrad[i]=[" << PGrad[i][0] << "," << PGrad[i][1]
                                  << "]  PGrad[j]=[" << PGrad[j][0] << "," << PGrad[j][1] << "]";
                    Logger(DEBUG) << "  vxGrad[i]=[" << vxGrad[i][0] << "," << vxGrad[i][1]
                                  << "]  vxGrad[j]=[" << vxGrad[j][0] << "," << vxGrad[j][1] << "]";
                    Logger(DEBUG) << "  vyGrad[i]=[" << vyGrad[i][0] << "," << vyGrad[i][1]
                                  << "]  vyGrad[j]=[" << vyGrad[j][0] << "," << vyGrad[j][1] << "]";
                    Logger(DEBUG) << "  RAW:    rhoR=" << WijR_raw[0]    << " PR=" << WijR_raw[1]
                                  << " vxR=" << WijR_raw[2] << " vyR=" << WijR_raw[3]
                                  << " | rhoL=" << WijL_raw[0]    << " PL=" << WijL_raw[1]
                                  << " vxL=" << WijL_raw[2] << " vyL=" << WijL_raw[3];
                    Logger(DEBUG) << "  EXTRAP: rhoR=" << WijR_extrap[0] << " PR=" << WijR_extrap[1]
                                  << " vxR=" << WijR_extrap[2] << " vyR=" << WijR_extrap[3]
                                  << " | rhoL=" << WijL_extrap[0] << " PL=" << WijL_extrap[1]
                                  << " vxL=" << WijL_extrap[2] << " vyL=" << WijL_extrap[3];
                    Logger(DEBUG) << "  PWL:    rhoR=" << WijR_pwl[0]    << " PR=" << WijR_pwl[1]
                                  << " vxR=" << WijR_pwl[2] << " vyR=" << WijR_pwl[3]
                                  << " | rhoL=" << WijL_pwl[0]    << " PL=" << WijL_pwl[1]
                                  << " vxL=" << WijL_pwl[2] << " vyL=" << WijL_pwl[3];
                    Logger(DEBUG) << "  FINAL:  rhoR=" << WijR[iW][0]    << " PR=" << WijR[iW][1]
                                  << " vxR=" << WijR[iW][2] << " vyR=" << WijR[iW][3]
                                  << " | rhoL=" << WijL[iW][0]    << " PL=" << WijL[iW][1]
                                  << " vxL=" << WijL[iW][2] << " vyL=" << WijL[iW][3];
                }
            }

            // Clamp reconstructed density to the floor. Free-surface particles
            // with steep rhoGrad can have extrapolation drive rho below zero at
            // the quadrature point; without this the corrupted state feeds HLLC
            // and NaN propagates through the subsequent flux update.
            if (DENSITY_FLOOR > 0.){
                if (WijR[iW][0] < DENSITY_FLOOR) WijR[iW][0] = DENSITY_FLOOR;
                if (WijL[iW][0] < DENSITY_FLOOR) WijL[iW][0] = DENSITY_FLOOR;
            }
        }
    }
}

double Particles::pairwiseLimiter(double phi0, double phi_i, double phi_j, double xijxi_abs, double xjxi_abs) {
    double phi_ = phi_i;

    /// calculate helper values
    double phi_ij = phi_i + xijxi_abs / xjxi_abs * (phi_j - phi_i);
    double phiMin, phiMax;
    if (phi_i < phi_j) {
        phiMin = phi_i;
        phiMax = phi_j;
    } else {
        phiMin = phi_j;
        phiMax = phi_i;
    }
    double delta1 = PSI_1 * abs(phi_i - phi_j);
    double delta2 = PSI_2 * abs(phi_i - phi_j);
    double phiMinus, phiPlus;
    if ((phiMax + delta1 >= 0. && phiMax >= 0.) || (phiMax + delta1 < 0. && phiMax < 0.)) {
        phiPlus = phiMax + delta1;
    } else {
        phiPlus = phiMax / (1. + delta1 / abs(phiMax));
    }
    if ((phiMin - delta1 >= 0. && phiMin >= 0.) || (phiMin - delta1 < 0. && phiMin < 0.)) {
        phiMinus = phiMin - delta1;
    } else {
        phiMinus = phiMin / (1. + delta1 / abs(phiMin));
    }

    /// actually compute the effective face limited value
    if (phi_i < phi_j) {
        double minPhiD2;
        if (phi_ij + delta2 < phi0){
            minPhiD2 = phi_ij + delta2;
        } else {
            minPhiD2 = phi0;
        }

        phi_ = phiMinus > minPhiD2 ? phiMinus : minPhiD2;

    } else if (phi_i > phi_j){
        double maxPhiD2;
        if (phi_ij - delta2 > phi0){
            maxPhiD2 = phi_ij - delta2;
        } else {
            maxPhiD2 = phi0;
        }

        phi_ = phiPlus < maxPhiD2 ? phiPlus : maxPhiD2;

    }
    return phi_;
}

void Particles::solveRiemannProblems(const Particles &ghostParticles){
#if USE_HLLC
    double n_unit[DIM];
#endif
    for (int i=0; i<N; ++i){

        // if (!(i % (N/VERBOSITY_PARTICLES))){
        //    Logger(DEBUG) << "        > i = " << i;
        // }
        // Logger(DEBUG) << "        > i = " << i << ", V = " << 1./omega[i];
        
        for (int j=0; j<noi[i]; ++j){
            // Logger(DEBUG) << "Ok " << j << " "<< i;
            int ii = i*MAX_NUM_INTERACTIONS+j; // interaction index
            // if (i == 9 && j == 11){
            //     Logger(DEBUG) << "i = " << i << ", ii = " << ii << ", j = " << nnl[ii];
            //     Logger(DEBUG) << "xi = [" << x[i] << ", " << y[i] << "], xj = ["
            //                  << x[nnl[ii]] << ", " << y[nnl[ii]] << "] , vi = ["
            //                  << vx[i] << ", " << vy[i] << "], vj = ["
            //                  << vx[nnl[ii]] << ", " << vy[nnl[ii]] << "]";
            //     Logger(DEBUG) << "vFrame = [" << vFrame[ii][0] << ", " << vFrame[ii][1] << "]";
            //     Logger(DEBUG) << "rhoL = " << WijL[ii][0] << ", rhoR = " << WijR[ii][0]
            //                  << ", uL = " << WijL[ii][2] << ", uR = " << WijR[ii][2]
            //                  << ", PL = " << WijL[ii][1] << ", PR = " << WijR[ii][1]
            //                  << ", Aij = [" << Aij[ii][0] << ", " << Aij[ii][1] << "]";
            // }

            // Logger(DEBUG) << "*WR = " << WijR[ii] << ", *WL = " << WijL[ii];
            // Logger(DEBUG) << "i = " << i << ", j = " << nnl[ii] << ", PR = " << WijR[ii][1] << ", PL = " << WijL[ii][1];
#if EOS == 0
            if(WijR[ii][1] < 0. || WijL[ii][1] < 0.){
                Logger(WARN) << "Negative pressure encountered@(i = " << i << ", j = " << j << ") Very bad :( !!";
//                 Logger(DEBUG) << "    > rhoL = " << WijL[ii][0] << ", rhoR = " << WijR[ii][0]
//                               << ", uL = " << WijL[ii][2] << ", uR = " << WijR[ii][2]
//                               << ", PL = " << WijL[ii][1] << ", PR = " << WijR[ii][1]
//                               << ", Aij = [" << Aij[ii][0] << ", " << Aij[ii][1]
// #if DIM ==3
//                               << ", " << Aij[ii][2]
// #endif
//                               << "]";

#if DEBUG_LVL > 1
                Logger(DEBUG) << "Aborting for debugging.";
                exit(6);
#endif
                //if(WijR[ii][1] < 0.){
                //    WijR[ii][1] = PRESSURE_FLOOR;
                //}
                //if(WijL[ii][1] < 0.){
                //    WijL[ii][1] = PRESSURE_FLOOR;
                //}
            }
#endif // EOS == 0
            bool compute = true;
            int iij;
#if ENFORCE_FLUX_SYM
            //int ii = i*MAX_NUM_INTERACTIONS+j; // interaction index i->j
            int ji = nnl[ii]; // index i of particle j
            if(ji<i) {
                compute = false;
                // search neighbor i in nnl[] of j
                int ij;
                for (ij = 0; ij < noi[ji]; ++ij) {
                    if (nnl[ij + ji * MAX_NUM_INTERACTIONS] == i) break;
                }
                iij = ji * MAX_NUM_INTERACTIONS + ij; // interaction index j->i
            }
#endif
            if (compute){
#if USE_HLLC
                // Logger(DEBUG) << " i = " << i << " j = " << nnl[ii] << " mF = " << Fij[ii][0];

            //    calcNunit(i, ii, n_unit);
                // Logger(DEBUG) << "Ok " << j << " "<< i;
#if ELASTIC
                // Extrapolate stress to face midpoints
                int jj = nnl[ii];
                double xijxi_loc[DIM], xijxj_loc[DIM];
                xijxi_loc[0] = .5*(x[jj] - x[i]);
                xijxi_loc[1] = .5*(y[jj] - y[i]);
                xijxj_loc[0] = -xijxi_loc[0];
                xijxj_loc[1] = -xijxi_loc[1];
#if DIM == 3
                xijxi_loc[2] = .5*(z[jj] - z[i]);
                xijxj_loc[2] = -xijxi_loc[2];
#endif
                double fSxxR = Sxx[i]  + Helper::dotProduct(SxxGrad[i], xijxi_loc);
                double fSxyR = Sxy[i]  + Helper::dotProduct(SxyGrad[i], xijxi_loc);
                double fSyyR = Syy[i]  + Helper::dotProduct(SyyGrad[i], xijxi_loc);
                double fSxxL = Sxx[jj] + Helper::dotProduct(SxxGrad[jj], xijxj_loc);
                double fSxyL = Sxy[jj] + Helper::dotProduct(SxyGrad[jj], xijxj_loc);
                double fSyyL = Syy[jj] + Helper::dotProduct(SyyGrad[jj], xijxj_loc);
#if GIZMO_ELASTIC_FLUX
                // GIZMO keeps the Riemann problem isotropic; the deviatoric stress
                // is added as a separate SPH flux after the solve (see below).
                fSxxR = fSxyR = fSyyR = fSxxL = fSxyL = fSyyL = 0.;
#endif
#if DIM == 3
                double fSxzR = Sxz[i]  + Helper::dotProduct(SxzGrad[i], xijxi_loc);
                double fSyzR = Syz[i]  + Helper::dotProduct(SyzGrad[i], xijxi_loc);
                double fSzzR = Szz[i]  + Helper::dotProduct(SzzGrad[i], xijxi_loc);
                double fSxzL = Sxz[jj] + Helper::dotProduct(SxzGrad[jj], xijxj_loc);
                double fSyzL = Syz[jj] + Helper::dotProduct(SyzGrad[jj], xijxj_loc);
                double fSzzL = Szz[jj] + Helper::dotProduct(SzzGrad[jj], xijxj_loc);
#endif
#endif
#if USE_MATID
                int matIdR_ii = matId[i];
                int matIdL_ii = matId[nnl[ii]];
#else
                int matIdR_ii = 0;
                int matIdL_ii = 0;
#endif
                Riemann solver { WijL[ii], WijR[ii], vFrame[ii], Aij[ii], i,
                                matIdL_ii, matIdR_ii, *MeshlessEOS
#if ELASTIC
#if TENSILE_CORRECTION
                                , fabMonaghan[ii]
#endif
                                , fSxxL, fSxyL, fSyyL
                                , fSxxR, fSxyR, fSyyR
#if DIM == 3
                                , fSxzL, fSyzL, fSzzL
                                , fSxzR, fSyzR, fSzzR
#endif
#if FRAGMENTATION
                                // R-slot = j-side, L-slot = i-side (matches Sij order)
                                , damageTotal[jj], damageTotal[i]
#endif
#endif
                };
                solver.HLLCFlux(Fij[ii]);
#if GIZMO_ELASTIC_FLUX
                double fTC = 0.;
#if TENSILE_CORRECTION
                fTC = TENSILE_CORRECTION_PREFAC * pow(fabMonaghan[ii], TENSILE_CORRECTION_POWER);
#endif
                addGizmoElasticStressFlux(i, jj, fTC, Fij[ii]);
#endif
#else
#if USE_MATID
                int matIdR_ii = matId[i];
                int matIdL_ii = matId[nnl[ii]];
#else
                int matIdR_ii = 0;
                int matIdL_ii = 0;
#endif
                Riemann solver { WijL[ii], WijR[ii], vFrame[ii], Aij[ii], i,
                                matIdL_ii, matIdR_ii, *MeshlessEOS};
                solver.exact(Fij[ii]);
#endif
            } else {
                // Logger(DEBUG) << " i = " << i << " j = " << nnl[ii] << " No compute";
                for(int d=0; d<DIM+2; ++d){
                    Fij[ii][d] = -Fij[iij][d];
                }
            }

            // if(i == 6){//&& j==11){
            //    Logger(DEBUG) << "Fluxes = [" << Fij[ii][0] << ", " << Fij[ii][1] << ", " << Fij[ii][2] << ", " << Fij[ii][3] << "]";
            // }
        }

#if PERIODIC_BOUNDARIES
        for (int j=0; j<noiGhosts[i]; ++j){
            int ii = i*MAX_NUM_GHOST_INTERACTIONS+j; // interaction index
            //Logger(DEBUG) << "i = " << i << ", ii = " << ii << ", jGhost = " << nnlGhosts[ii];
            //Logger(DEBUG) << "xi = [" << x[i] << ", " << y[i] << "], xj = ["
            //              << ghostParticles.x[nnlGhosts[ii]] << ", " << ghostParticles.y[nnlGhosts[ii]] << "]";
            //Logger(DEBUG) << "vFrameGhosts = [" << vFrameGhosts[ii][0] << ", " << vFrameGhosts[ii][1] << "]";

#if EOS == 0
            if(WijRGhosts[ii][1] < 0. || WijLGhosts[ii][1] < 0.){
                Logger(WARN) << "Negative pressure encountered@(i = " << i << ", jGhost = " << j << ") Very bad :( !!";
                Logger(DEBUG) << "rhoL = " << WijLGhosts[ii][0] << ", rhoR = " << WijRGhosts[ii][0]
                              << ", uL = " << WijLGhosts[ii][2] << ", uR = " << WijRGhosts[ii][2]
                              << ", PL = " << WijLGhosts[ii][1] << ", PR = " << WijRGhosts[ii][1];
#if DEBUG_LVL
                Logger(DEBUG) << "Aborting for debugging.";
                exit(6);
#endif
            }

#endif // EOS == 0
            bool compute = true;
            int iij;
#if ENFORCE_FLUX_SYM
            //int ii = i*MAX_NUM_GHOST_INTERACTIONS+j; // interaction index i->j
            int ji = nnlGhosts[ii]; // index i of particle j
            if (ghostParticles.parent[ji]<i){
                const int p = ghostParticles.parent[ji];
                // Image-aware reverse match: the forward ghost ji sits at
                // x[p]+delta; the reciprocal interaction is p with the image of
                // i at x[i]-delta. Matching by parent id alone is ambiguous when
                // i is reachable through more than one periodic image, so pick
                // the ghost-of-i closest to that expected position.
                const double tgtX = x[i] - (ghostParticles.x[ji] - x[p]);
                const double tgtY = y[i] - (ghostParticles.y[ji] - y[p]);
                double best = std::numeric_limits<double>::max();
                int ijBest = -1;
                for (int ij=0; ij<noiGhosts[p]; ++ij){
                    const int gg = nnlGhosts[ij+p*MAX_NUM_GHOST_INTERACTIONS];
                    if (ghostParticles.parent[gg] != i) continue;
                    double d = pow(ghostParticles.x[gg]-tgtX, 2)
                             + pow(ghostParticles.y[gg]-tgtY, 2);
                    if (d < best){ best = d; ijBest = ij; }
                }
                if (ijBest < 0){
                    // Non-reciprocal: p has no ghost of i. Compute directly
                    // rather than copy from a stale/garbage slot.
                    compute = true;
                } else {
                    compute = false;
                    iij = ijBest+p*MAX_NUM_GHOST_INTERACTIONS;
                }
            }
#endif

            if (compute) {
#if USE_HLLC
                // Logger(DEBUG) << " i = " << i << " j = " << nnl[ii] << " mF = " << Fij[ii][0];

                // calcNunit(i, ii, n_unit);

#if ELASTIC
                // Extrapolate stress to face midpoints (ghost interactions)
                int jj = nnlGhosts[ii];
                double xijxi_loc[DIM], xijxj_loc[DIM];
                xijxi_loc[0] = .5*(ghostParticles.x[jj] - x[i]);
                xijxi_loc[1] = .5*(ghostParticles.y[jj] - y[i]);
                xijxj_loc[0] = -xijxi_loc[0];
                xijxj_loc[1] = -xijxi_loc[1];
#if DIM == 3
                xijxi_loc[2] = .5*(ghostParticles.z[jj] - z[i]);
                xijxj_loc[2] = -xijxi_loc[2];
#endif
                double fSxxR = Sxx[i] + Helper::dotProduct(SxxGrad[i], xijxi_loc);
                double fSxyR = Sxy[i] + Helper::dotProduct(SxyGrad[i], xijxi_loc);
                double fSyyR = Syy[i] + Helper::dotProduct(SyyGrad[i], xijxi_loc);
                double fSxxL = ghostParticles.Sxx[jj] + Helper::dotProduct(ghostParticles.SxxGrad[jj], xijxj_loc);
                double fSxyL = ghostParticles.Sxy[jj] + Helper::dotProduct(ghostParticles.SxyGrad[jj], xijxj_loc);
                double fSyyL = ghostParticles.Syy[jj] + Helper::dotProduct(ghostParticles.SyyGrad[jj], xijxj_loc);
#if DIM == 3
                double fSxzR = Sxz[i] + Helper::dotProduct(SxzGrad[i], xijxi_loc);
                double fSyzR = Syz[i] + Helper::dotProduct(SyzGrad[i], xijxi_loc);
                double fSzzR = Szz[i] + Helper::dotProduct(SzzGrad[i], xijxi_loc);
                double fSxzL = ghostParticles.Sxz[jj] + Helper::dotProduct(ghostParticles.SxzGrad[jj], xijxj_loc);
                double fSyzL = ghostParticles.Syz[jj] + Helper::dotProduct(ghostParticles.SyzGrad[jj], xijxj_loc);
                double fSzzL = ghostParticles.Szz[jj] + Helper::dotProduct(ghostParticles.SzzGrad[jj], xijxj_loc);
#endif
#endif
#if USE_MATID
                int matIdR_g = matId[i];
                int matIdL_g = ghostParticles.matId[nnlGhosts[ii]];
#else
                int matIdR_g = 0;
                int matIdL_g = 0;
#endif
                Riemann solver{WijLGhosts[ii], WijRGhosts[ii], vFrameGhosts[ii], AijGhosts[ii], i,
                               matIdL_g, matIdR_g, *MeshlessEOS
#if ELASTIC
                               , fSxxR, fSxyR, fSyyR
                               , fSxxL, fSxyL, fSyyL
#if DIM == 3
                               , fSxzR, fSyzR, fSzzR
                               , fSxzL, fSyzL, fSzzL
#endif
#if FRAGMENTATION
                               // R-slot = i-side here; ghost-side damage is not tracked.
                               , damageTotal[i], 0.0
#endif
#endif
                };
                solver.HLLCFlux(FijGhosts[ii]);
#else
#if USE_MATID
                int matIdR_g = matId[i];
                int matIdL_g = ghostParticles.matId[nnlGhosts[ii]];
#else
                int matIdR_g = 0;
                int matIdL_g = 0;
#endif
                Riemann solver{WijLGhosts[ii], WijRGhosts[ii], vFrameGhosts[ii], AijGhosts[ii], i,
                               matIdL_g, matIdR_g, *MeshlessEOS};
                solver.exact(FijGhosts[ii]);
#endif
            } else {
                // Logger(DEBUG) << " i = " << i << " j = " << nnl[ii] << " No compute";

                for(int d=0; d<DIM+2; ++d){
                    FijGhosts[ii][d] = -FijGhosts[iij][d];
                }
            }
        }
        
#endif
    }
}

void Particles::collectFluxes(Helper &helper){
#if DIAG_COND_ENABLE
    static int diagStepCounter = 0;
    static int diagActivatedAt = -1;
    static int diagTarget = -1;
    ++diagStepCounter;
    bool diagActive = false;

    if (diagTarget < 0) {
        double worstCond = DIAG_COND_TRIGGER;
        int worstIdx = -1;
        for (int k = 0; k < N; ++k) {
            if (conditionNumber[k] > worstCond) {
                worstCond = conditionNumber[k];
                worstIdx = k;
            }
        }
        if (worstIdx >= 0) {
            diagTarget = worstIdx;
            diagActivatedAt = diagStepCounter;
            Logger(DEBUG) << "TARGET DIAG LOCKED i=" << diagTarget
                          << " cond=" << worstCond
                          << " diagStep=" << diagStepCounter
                          << " xi=[" << x[diagTarget] << "," << y[diagTarget] << "]"
                          << " vi=[" << vx[diagTarget] << "," << vy[diagTarget] << "]";
        }
    }
    if (diagTarget >= 0
        && diagStepCounter - diagActivatedAt < DIAG_WINDOW_STEPS) {
        diagActive = true;
    }
#endif
    for (int i=0; i<N; ++i){
        mF[i] = 0.;

        vF[i][0] = 0.;
        vF[i][1] = 0.;
#if DIM == 3
        vF[i][2] = 0.;
#endif
        eF[i] = 0.;

        //Logger(DEBUG) << "      > i = " << i;

        for(int j=0; j<noi[i]; ++j){
            int ii = j+i*MAX_NUM_INTERACTIONS;
            double AijNorm = sqrt(Helper::dotProduct(Aij[ii], Aij[ii]));

            /// MASS FLUXES
            //mF[i] += AijNorm*Fij[ii][0];
            mF[i] += Fij[ii][0];

            // Logger(DEBUG) << " i: " << i << " xi = [" << x[i] << ", " << y[i] << "]"
            //          << ", xj = [" << x[nnl[ii]] << ", " << y[nnl[ii]] << "]"
            //          << ", mF[ii] = " << Fij[ii][0] << ", AijNorm = " << AijNorm;

            // Debug mF = nan
            // Logger(DEBUG) << " i: " << i << " j: " << nnl[ii] << " mF[ii] = " << Fij[ii][0];

            /// VELOCITY FLUXES
            // add de-boosted velocities
            //vF[i][0] += AijNorm * (Fij[ii][2] + Fij[ii][0]*vFrame[ii][0]);
            //vF[i][1] += AijNorm * (Fij[ii][3] + Fij[ii][0]*vFrame[ii][1]);

            vF[i][0] += Fij[ii][2];
            vF[i][1] += Fij[ii][3];
#if DIM==3
            vF[i][2] += Fij[ii][4];
#endif

            /// ENERGY FLUXES
            // allocate buffer for energy update
            //double Fv[DIM];
            //Fv[0] = Fij[ii][2];
            //Fv[1] = Fij[ii][3];
            //eF[i] += AijNorm*(Fij[ii][1] + .5*Helper::dotProduct(vFrame[ii], vFrame[ii])*Fij[ii][0]
            //                  + Helper::dotProduct(vFrame[ii], Fv));

            eF[i] += Fij[ii][1];

#if DIAG_COND_ENABLE
            if (diagActive && i == diagTarget){
                int jIdx = nnl[ii];
                double dx = x[jIdx] - x[i];
                double dy = y[jIdx] - y[i];
                double rij = sqrt(dx*dx + dy*dy);
                Logger(DEBUG) << "TARGET DIAG PAIR i=" << i
                              << " j=" << jIdx
                              << " dxij=[" << dx << "," << dy << "]"
                              << " rij=" << rij
                              << " Fxy_pair=[" << Fij[ii][2] << "," << Fij[ii][3] << "]"
                              << " Aij=[" << Aij[ii][0] << "," << Aij[ii][1] << "]"
                              << " cond_j=" << conditionNumber[jIdx]
                              << " noi_j=" << noi[jIdx];
            }
#endif

        }

#if PERIODIC_BOUNDARIES
        for(int j=0; j<noiGhosts[i]; ++j){
            int ii = j+i*MAX_NUM_GHOST_INTERACTIONS;
            // double AijNorm = sqrt(Helper::dotProduct(AijGhosts[ii], AijGhosts[ii]));

            /// MASS FLUXES
            //mF[i] += AijNorm*FijGhosts[ii][0];
            mF[i] += FijGhosts[ii][0];

            //Logger(DEBUG) << "xi = [" << x[i] << ", " << y[i] << "]"
            //              << ", xjGhost = [" << ghostParticles.x[nnlGhosts[ii]] << ", " << ghostParticles.y[nnlGhosts[ii]] << "]"
            //              << ", mF[ii] = " << FijGhosts[ii][0] << ", AijNorm = " << AijNorm;

            /// VELOCITY FLUXES
            // add de-boosted velocities
            //vF[i][0] += AijNorm * (FijGhosts[ii][2] + FijGhosts[ii][0]*vFrameGhosts[ii][0]);
            //vF[i][1] += AijNorm * (FijGhosts[ii][3] + FijGhosts[ii][0]*vFrameGhosts[ii][1]);

            vF[i][0] += FijGhosts[ii][2];
            vF[i][1] += FijGhosts[ii][3];


            // TODO: implement z-component to work for 3D

            /// ENERGY FLUXES
            // allocate buffer for energy update
            //double Fv[DIM];
            //Fv[0] = FijGhosts[ii][2];
            //Fv[1] = FijGhosts[ii][3];
            //eF[i] += AijNorm*(FijGhosts[ii][1] + .5*Helper::dotProduct(vFrameGhosts[ii], vFrameGhosts[ii])*FijGhosts[ii][0]
            //                  + Helper::dotProduct(vFrameGhosts[ii], Fv));

            eF[i] += FijGhosts[ii][1];

            //if (i == 46){
            //    Logger(DEBUG) << "  > jGhost = " << nnlGhosts[ii] << ", AijNorm = " << AijNorm
            //                  << ", vFrame = [" << vFrameGhosts[ii][0] << ", " << vFrameGhosts[ii][1] << "]";
            //    Logger(DEBUG) << "  > Fm = " << mF[i] << ", Fv = [" << vF[i][0] << ", " << vF[i][1]
            //                  << "], Fe = " << eF[i];
            //}

            if (i == 184){
               // Logger(DEBUG) << "  > j = " << nnl[ii] << ", AijNorm = " << AijNorm
               //           << ", vFrame = [" << vFrame[ii][0] << ", " << vFrame[ii][1] << "]";
               Logger(DEBUG) << "  > j = " << nnl[ii] <<  " > Fm = " << FijGhosts[ii][0] << ", Fv = [" << FijGhosts[ii][2] << ", " << FijGhosts[ii][3]
                             << "], dFe = " << FijGhosts[ii][1]
                            //  << " v_ij = " << vFrame[ii][0] << " " << vFrame[i][1]
                             ;
            }
        }
#endif

#if DIAG_COND_ENABLE
        if (diagActive && i == diagTarget){
            // Per-step NET dump AFTER all neighbours summed. F_net = -vF
            // because updateState does Q[1] -= dt*vF[i][0].
            Logger(DEBUG) << "TARGET DIAG NET diagStep=" << diagStepCounter
                          << " activatedAt=" << diagActivatedAt
                          << " i=" << i
                          << " xi=[" << x[i] << "," << y[i] << "]"
                          << " vi=[" << vx[i] << "," << vy[i] << "]"
                          << " rho=" << rho[i] << " P=" << P[i]
                          << " sml=" << sml[i] << " noi=" << noi[i]
                          << " cond=" << conditionNumber[i]
                          << " vF=[" << vF[i][0] << "," << vF[i][1] << "]"
                          << " Fnet=[" << -vF[i][0] << "," << -vF[i][1] << "]";
        }
#endif
    }
}

void Particles::updateStateAndPosition(const double &dt, const Domain &domain){

    for(int i=0; i<N; ++i){

        // store velocity for position update
        double vxi = vx[i], vyi = vy[i];
#if DIM==3
        double vzi = vz[i];
#endif

        // create state vector without mass
        double Q[DIM+1];
#if DIM==3
        Q[0] = m[i]*(u[i] + .5*(vxi*vxi+vyi*vyi+vzi*vzi));
#else
        Q[0] = m[i]*(u[i] + .5*(vxi*vxi+vyi*vyi));
#endif
        Q[1] = m[i]*vxi;
        Q[2] = m[i]*vyi;
#if DIM==3
        Q[3] = m[i]*vzi;
#endif
        //if(i == 46){
        //    Logger(DEBUG) << "i = " << i << ", mass flux = " << mF[i]
        //                  << ", momentum flux = [" << vF[i][0] << ", " << vF[i][1] << "], energy flux = " << eF[i];
        //}

        // UPDATE MASS
#if !MESHLESS_FINITE_MASS
        m[i] -= dt*mF[i];
#endif
        if (m[i] <= 0.){
            Logger(ERROR) << "Negative mass. m[" << i << "] =" << m[i] << ", mF = " << mF[i];
            //m[i] = MASS_FLOOR;
        }

        // UPDATE VELOCITY
        Q[1] -= dt*vF[i][0];
        Q[2] -= dt*vF[i][1];
#if DIM == 3
        Q[3] -= dt*vF[i][2];
#endif
        vx[i] = Q[1]/m[i];
        vy[i] = Q[2]/m[i];
#if DIM==3
        vz[i] = Q[3]/m[i];
        //Logger(DEBUG) << "New velocity v = [" << vx[i] << ", " << vy[i] << ", " << vz[i] << "]";
#endif
        // UPDATE INTERNAL ENERGY
        Q[0] -= dt*eF[i];
#if DIM==3
        u[i] = Q[0]/m[i]-.5*(vx[i]*vx[i]+vy[i]*vy[i]+vz[i]*vz[i]);

        //Logger(DEBUG) << "Total energy Q_E = " << Q[0] << ", kinetic energy E_kin = "
        //          << .5*(vx[i]*vx[i]+vy[i]*vy[i]+vz[i]*vz[i]);

#else
        u[i] = Q[0]/m[i]-.5*(vx[i]*vx[i]+vy[i]*vy[i]);
#endif
#if EOS != 1 // For murnaghan this is ok
        if (u[i] <= 0.){
            Logger(ERROR) << "Negative internal energy. u[" << i << "] =" << u[i] << ", eF = " << eF[i];
        }
#endif

#if ENERGY_FLOOR >= 0
        if (u[i] < ENERGY_FLOOR){
            u[i] = ENERGY_FLOOR;
            Logger(ERROR) << "Clamped internal energy at " << ENERGY_FLOOR;
        }
#endif

#if MOVE_PARTICLES
        // MOVE PARTICLES
        //x[i] += .5*(vx[i]+vxi)*dt;
        //y[i] += .5*(vy[i]+vyi)*dt;
        // stay consistent with choice of effective frame
        x[i] += vxi*dt;
        y[i] += vyi*dt;

#if DIM==3
        //z[i] += .5*(vz[i]+vzi)*dt;
        z[i] += vzi*dt;
#endif
#if PERIODIC_BOUNDARIES
        if (x[i] < domain.bounds.minX) {
            x[i] = domain.bounds.maxX - (domain.bounds.minX - x[i]);
        } else if (domain.bounds.maxX <= x[i]) {
            x[i] = domain.bounds.minX + (x[i] - domain.bounds.maxX);
        }
        if (y[i] < domain.bounds.minY) {
            y[i] = domain.bounds.maxY - (domain.bounds.minY - y[i]);
        } else if (domain.bounds.maxY <= y[i]) {
            y[i] = domain.bounds.minY + (y[i] - domain.bounds.maxY);
        }
#if DIM ==3
        if (z[i] < domain.bounds.minZ) {
            z[i] = domain.bounds.maxZ - (domain.bounds.minZ - z[i]);
        } else if (domain.bounds.maxZ <= z[i]) {
            z[i] = domain.bounds.minZ + (z[i] - domain.bounds.maxZ);
        }
#endif
#endif // PERIODIC_BOUNDARIES
#endif // MOVE_PARTICLES
    }
}

void Particles::updateState(const double &dt){
    for(int i=0; i<N; ++i){
        double vxi = vx[i], vyi = vy[i];
#if DIM==3
        double vzi = vz[i];
#endif
        double Q[DIM+1];
#if DIM==3
        Q[0] = m[i]*(u[i] + .5*(vxi*vxi+vyi*vyi+vzi*vzi));
#else
        Q[0] = m[i]*(u[i] + .5*(vxi*vxi+vyi*vyi));
#endif
        Q[1] = m[i]*vxi;
        Q[2] = m[i]*vyi;
#if DIM==3
        Q[3] = m[i]*vzi;
#endif

#if !MESHLESS_FINITE_MASS
        m[i] -= dt*mF[i];
#endif
        if (m[i] <= 0.){
            Logger(ERROR) << "Negative mass. m[" << i << "] =" << m[i] << ", mF = " << mF[i];
        }

        Q[1] -= dt*vF[i][0];
        Q[2] -= dt*vF[i][1];
#if DIM == 3
        Q[3] -= dt*vF[i][2];
#endif
        vx[i] = Q[1]/m[i];
        vy[i] = Q[2]/m[i];
#if DIM==3
        vz[i] = Q[3]/m[i];
#endif

        Q[0] -= dt*eF[i];
#if DIM==3
        u[i] = Q[0]/m[i]-.5*(vx[i]*vx[i]+vy[i]*vy[i]+vz[i]*vz[i]);
#else
        u[i] = Q[0]/m[i]-.5*(vx[i]*vx[i]+vy[i]*vy[i]);
#endif
#if EOS != 1 // Ok for murnaghan
        if (u[i] <= 0.){
            Logger(ERROR) << "Negative internal energy. u[" << i << "] =" << u[i] << ", eF = " << eF[i];
        }
#endif
#if ENERGY_FLOOR >= 0
        if (u[i] < ENERGY_FLOOR){
            u[i] = ENERGY_FLOOR;
        }
#endif
    }
}

void Particles::moveParticles(const double &dt, const Domain &domain){
#if MOVE_PARTICLES
    for(int i=0; i<N; ++i){
        x[i] += vx[i]*dt;
        y[i] += vy[i]*dt;
#if DIM==3
        z[i] += vz[i]*dt;
#endif
#if PERIODIC_BOUNDARIES
        if (x[i] < domain.bounds.minX) {
            x[i] = domain.bounds.maxX - (domain.bounds.minX - x[i]);
        } else if (domain.bounds.maxX <= x[i]) {
            x[i] = domain.bounds.minX + (x[i] - domain.bounds.maxX);
        }
        if (y[i] < domain.bounds.minY) {
            y[i] = domain.bounds.maxY - (domain.bounds.minY - y[i]);
        } else if (domain.bounds.maxY <= y[i]) {
            y[i] = domain.bounds.minY + (y[i] - domain.bounds.maxY);
        }
#if DIM==3
        if (z[i] < domain.bounds.minZ) {
            z[i] = domain.bounds.maxZ - (domain.bounds.minZ - z[i]);
        } else if (domain.bounds.maxZ <= z[i]) {
            z[i] = domain.bounds.minZ + (z[i] - domain.bounds.maxZ);
        }
#endif
#endif // PERIODIC_BOUNDARIES
    }
#endif // MOVE_PARTICLES
}
#if PERIODIC_BOUNDARIES
void Particles::createGhostParticles(Domain &domain, Particles &ghostParticles,
                                     const double &kernelSize){
    int iGhost = 0;
    Logger(DEBUG) << "N = " << N;
    for(int i=0; i<N; ++i) {
        // Logger(DEBUG) << "i = " << i << ", N = "<< N;
        bool foundGhostX = false, foundGhostY = false;
        ghostMap[i*(DIM+1)] = -1;
        ghostMap[i*(DIM+1)+1] = -1;
        ghostMap[i*(DIM+1)+2] = -1;

        // x-direction
        if (x[i] <= domain.bounds.minX + kernelSize){ // && x[i] > domain.bounds.minX) {
            ghostParticles.x[iGhost] = domain.bounds.maxX + (x[i] - domain.bounds.minX);
            foundGhostX = true;
        } else if (domain.bounds.maxX - kernelSize < x[i]){ // && x[i] < domain.bounds.maxX) {
            ghostParticles.x[iGhost] = domain.bounds.minX - (domain.bounds.maxX - x[i]);
            foundGhostX = true;
        } else {
            ghostParticles.x[iGhost] = x[i];
        }

        // y-direction
        if (y[i] <= domain.bounds.minY + kernelSize){ // && y[i] > domain.bounds.minY) {
            ghostParticles.y[iGhost] = domain.bounds.maxY + (y[i] - domain.bounds.minY);
            foundGhostY = true;
        } else if (domain.bounds.maxY - kernelSize < y[i]){ // && y[i] < domain.bounds.maxY) {
            ghostParticles.y[iGhost] = domain.bounds.minY - (domain.bounds.maxY - y[i]);
            foundGhostY = true;
        } else {
            ghostParticles.y[iGhost] = y[i];
        }
        // Logger(DEBUG) << "corner-direction:";

        // 'corner' particle first if both are true
        if (foundGhostX || foundGhostY) {
            ghostMap[i*(DIM+1)] = iGhost;
            ghostParticles.parent[iGhost] = i;
            // Keep ghost sml current with the parent at creation time. ghostNNS
            // runs before updateGhostState, so without this the symmetric
            // max(h_i,h_j) neighbour criterion would read stale, mis-indexed
            // smoothing lengths and produce non-reciprocal ghost neighbour lists.
            ghostParticles.sml[iGhost] = sml[i];
            //Logger(DEBUG) << "particle@" << i << " = [" << x[i] << ", " << y[i] << "] makes "
            //          << "ghost@" << iGhost << " = [" << ghostParticles.x[iGhost] << ", "
            //          << ghostParticles.y[iGhost] << "]";
            ++iGhost;
        }
        // Logger(DEBUG) << "corner ok";

        // create DIM extra normal particles
        if (foundGhostX && foundGhostY){
            ghostParticles.x[iGhost] = x[i];
            if (y[i] <= domain.bounds.minY + kernelSize){ // && y[i] > domain.bounds.minY) {
                ghostParticles.y[iGhost] = domain.bounds.maxY + (y[i] - domain.bounds.minY);
            } else if (domain.bounds.maxY - kernelSize < y[i]){ // && y[i] < domain.bounds.maxY) {
                ghostParticles.y[iGhost] = domain.bounds.minY - (domain.bounds.maxY - y[i]);
            }
            ghostMap[i*(DIM+1)+1] = iGhost;
            ghostParticles.parent[iGhost] = i;
            ghostParticles.sml[iGhost] = sml[i];
            //Logger(DEBUG) << "particle@" << i << " = [" << x[i] << ", " << y[i] << "] makes "
            //              <<"ghost@" << iGhost << " = [" << ghostParticles.x[iGhost] << ", "
            //              << ghostParticles.y[iGhost] << "]";
            ++iGhost;
            if (x[i] <= domain.bounds.minX + kernelSize){ // && x[i] > domain.bounds.minX) {
                ghostParticles.x[iGhost] = domain.bounds.maxX + (x[i] - domain.bounds.minX);
            } else if (domain.bounds.maxX - kernelSize < x[i]){ // && x[i] < domain.bounds.maxX) {
                ghostParticles.x[iGhost] = domain.bounds.minX - (domain.bounds.maxX - x[i]);
            }
            ghostParticles.y[iGhost] = y[i];
            ghostMap[i*(DIM+1)+2] = iGhost;
            ghostParticles.parent[iGhost] = i;
            ghostParticles.sml[iGhost] = sml[i];
            //Logger(DEBUG) << "particle@" << i << " = [" << x[i] << ", " << y[i] << "] makes "
            //              << "ghost@" << iGhost << " = [" << ghostParticles.x[iGhost] << ", "
            //              << ghostParticles.y[iGhost] << "]";
            ++iGhost;
        }

        //foundGhostX = false;
        //foundGhostY = false;

#if DIM == 3
        Logger(ERROR) << "Ghost cells not implemented for 3D simulations. - Aborting.";
        exit(2);
#endif
    }
    ghostParticles.N = iGhost;
}

void Particles::updateGhostState(Particles &ghostParticles){
    for (int i=0; i<N*(DIM+1); ++i){
        if (ghostMap[i] >= 0){
#if USE_MATID
            ghostParticles.matId[ghostMap[i]] = matId[i/(DIM+1)];
#endif
            ghostParticles.rho[ghostMap[i]] = rho[i/(DIM+1)];
            ghostParticles.P[ghostMap[i]] = P[i/(DIM+1)];
            ghostParticles.omega[ghostMap[i]] = omega[i/(DIM+1)];
            ghostParticles.vx[ghostMap[i]] = vx[i/(DIM+1)];
            ghostParticles.vy[ghostMap[i]] = vy[i/(DIM+1)];
#if DIM==3
            ghostParticles.vz[ghostMap[i]] = vz[i/(DIM+1)];
#endif
            ghostParticles.sml[ghostMap[i]] = sml[i/(DIM+1)];
#if ELASTIC
            ghostParticles.Sxx[ghostMap[i]] = Sxx[i/(DIM+1)];
            ghostParticles.Sxy[ghostMap[i]] = Sxy[i/(DIM+1)];
            ghostParticles.Syy[ghostMap[i]] = Syy[i/(DIM+1)];
#if DIM == 3
            ghostParticles.Sxz[ghostMap[i]] = Sxz[i/(DIM+1)];
            ghostParticles.Syz[ghostMap[i]] = Syz[i/(DIM+1)];
            ghostParticles.Szz[ghostMap[i]] = Szz[i/(DIM+1)];
#endif // DIM == 3
#endif // ELASTIC
        }
    }
}

void Particles::updateGhostGradients(Particles &ghostParticles){
    for (int i=0; i<N*(DIM+1); ++i){
        if (ghostMap[i] >= 0){
            for (int alpha=0; alpha<DIM; ++alpha){
                ghostParticles.rhoGrad[ghostMap[i]][alpha] = rhoGrad[i/(DIM+1)][alpha];
                ghostParticles.vxGrad[ghostMap[i]][alpha] = vxGrad[i/(DIM+1)][alpha];
                ghostParticles.vyGrad[ghostMap[i]][alpha] = vyGrad[i/(DIM+1)][alpha];
#if DIM == 3
                ghostParticles.vzGrad[ghostMap[i]][alpha] = vzGrad[i/(DIM+1)][alpha];
#endif
                ghostParticles.PGrad[ghostMap[i]][alpha] = PGrad[i/(DIM+1)][alpha];
#if ELASTIC
                ghostParticles.SxxGrad[ghostMap[i]][alpha] = SxxGrad[i/(DIM+1)][alpha];
                ghostParticles.SxyGrad[ghostMap[i]][alpha] = SxyGrad[i/(DIM+1)][alpha];
                ghostParticles.SyyGrad[ghostMap[i]][alpha] = SyyGrad[i/(DIM+1)][alpha];
#if DIM == 3
                ghostParticles.SxzGrad[ghostMap[i]][alpha] = SxzGrad[i/(DIM+1)][alpha];
                ghostParticles.SyzGrad[ghostMap[i]][alpha] = SyzGrad[i/(DIM+1)][alpha];
                ghostParticles.SzzGrad[ghostMap[i]][alpha] = SzzGrad[i/(DIM+1)][alpha];
#endif // DIM == 3
#endif // ELASTIC
            }
        }
    }
}

/*void Particles::updateGhostPsijTilde(Particles &ghostParticles){
    for (int i=0; i<N*(DIM+1); ++i){
        if (ghostMap[i] >= 0){
            for (int j=0; j<noiGhosts[i/(DIM+1)]; ++j){
                for (int alpha=0; alpha<DIM; ++alpha){
                    ghostParticles.psijTilde_xiGhosts[j+ghostMap[i]*MAX_NUM_GHOST_INTERACTIONS][alpha]
                                    = psijTilde_xiGhosts[j+i/(DIM+1)*MAX_NUM_GHOST_INTERACTIONS][alpha];
                }
            }
        }
    }
}*/

void Particles::ghostNNS(Domain &domain, const Particles &ghostParticles, const double &kernelSize){
#if !VARIABLE_SML
    const double hSqr = kernelSize * kernelSize;
#endif
    for(int i=0; i<N; ++i){
        int noiBuf = 0;
#if VARIABLE_SML
        const double hi = sml[i];
#endif
        for(int iGhost=0; iGhost<ghostParticles.N; ++iGhost){
            double dSqr = pow(ghostParticles.x[iGhost] - x[i], 2)
                          + pow(ghostParticles.y[iGhost] - y[i], 2);
#if DIM == 3
            Logger(ERROR) << "Ghost cells not implemented for 3D simulations. - Aborting.";
        exit(2);
#endif
#if VARIABLE_SML
            // Use the parent's live sml directly rather than the ghost-array
            // copy: createGhostParticles repacks the ghost array every step and
            // ghostNNS runs before updateGhostState, so ghostParticles.sml[iGhost]
            // can be stale/mis-indexed. The parent index is always current.
            const double hPair = std::max(hi, sml[ghostParticles.parent[iGhost]]);
            const double hSqr = hPair * hPair;
#endif
            if (dSqr < hSqr){
                if(noiBuf >= MAX_NUM_GHOST_INTERACTIONS){
                    Logger(ERROR) << "MAX_NUM_GHOST_INTERACTIONS exceeded for particle "
                                  << i << " - Aborting.";
                    exit(3);
                }
                nnlGhosts[noiBuf+i*MAX_NUM_GHOST_INTERACTIONS] = iGhost;
                ++noiBuf;
            }
        }
        noiGhosts[i] = noiBuf;
    }
}

void Particles::compDensity(const Particles &ghostParticles){
    for(int i=0; i<N; ++i){
        compOmega(i, ghostParticles);
        rho[i] = m[i]*omega[i];
        if(rho[i] <= 0.){
            Logger(WARN) << "Zero or negative density @" << i;
        }
    }
}

void Particles::compOmega(int i, const Particles &ghostParticles){
    const double hi = sml[i];
    double omg = omega[i];
    for (int j=0; j<noiGhosts[i]; ++j){
        double dSqr = pow(x[i] - ghostParticles.x[nnlGhosts[j+i*MAX_NUM_GHOST_INTERACTIONS]], 2)
                      + pow(y[i] - ghostParticles.y[nnlGhosts[j+i*MAX_NUM_GHOST_INTERACTIONS]], 2);
#if DIM == 3
        dSqr += pow(z[i] - ghostParticles.z[nnlGhosts[j+i*MAX_NUM_GHOST_INTERACTIONS]], 2);
#endif
        double r = sqrt(dSqr);
        omg += kernel(r, hi);

        //Logger(DEBUG) << "x[" << i << "] = [" << x[i] << ", " <<  y[i] << "], x["
        //              << nnlGhosts[j+i*MAX_NUM_GHOST_INTERACTIONS] << "] = [" << ghostParticles.x[nnlGhosts[j+i*MAX_NUM_GHOST_INTERACTIONS]] << ", "
        //              << ghostParticles.y[nnlGhosts[j+i*MAX_NUM_GHOST_INTERACTIONS]] << "]";
    }
    omega[i] = omg;
    //Logger(DEBUG) << "V[" << i << "] = " << 1./omega[i] << ", noiTot = " << noi[i] + noiGhosts[i]
    //          << " noiGhosts = " << noiGhosts[i];
}

void Particles::compPsijTilde(Helper &helper, const Particles &ghostParticles){
    int iP;
    for (int i=0; i<N; ++i){
        const double hi = sml[i];

        // reset buffer
        for (int k=0; k<DIM*DIM; ++k){
            B[k] = 0.;
        }

        //for (int alpha = 0; alpha < DIM; ++alpha) {
        //    psijTilde_xi[i][alpha] = 0.;
        //}

        xi[0] = x[i];
        xi[1] = y[i];
#if DIM==3
        xi[2] = z[i];
#endif

        //Logger(DEBUG) << "In compPsijTilde(): noi[" << i << "] = " << noi[i]
        //        << ", noiGhosts[" << i << "] = " << noiGhosts[i];

        for (int j=0; j<noi[i]; ++j){
            iP = nnl[j+i*MAX_NUM_INTERACTIONS];
            double dSqr = pow(x[i] - x[iP], 2)
                          + pow(y[i] - y[iP], 2);
#if DIM == 3
            dSqr += pow(z[i] - z[iP], 2);
#endif
            double r = sqrt(dSqr);
            //double psij_xi = kernel(r, hi)/omega[nnl[j + i * MAX_NUM_INTERACTIONS]];
            double psij_xi = kernel(r, hi)/omega[i];

            xj[0] = x[iP];
            xj[1] = y[iP];
#if DIM==3
            xj[2] = z[iP];
#endif

            //if (i == 7){
            //    Logger(DEBUG) << "psij_xi = " << psij_xi << ", xj = [" << xj[0] << ", " << xj[1] << "]"
            //                  << "; xi = [" << xi[0] << ", " << xi[1] << "]";
            //}

            for (int alpha=0; alpha<DIM; ++alpha){
                for(int beta=0; beta<DIM; ++beta){
                    B[DIM*alpha+beta] += (xj[alpha] - xi[alpha])*(xj[beta] - xi[beta]) * psij_xi;
                }
            }
        }

        for (int j=0; j<noiGhosts[i]; ++j){

            double dSqr = pow(x[i] - ghostParticles.x[nnlGhosts[j+i*MAX_NUM_GHOST_INTERACTIONS]], 2)
                          + pow(y[i] - ghostParticles.y[nnlGhosts[j+i*MAX_NUM_GHOST_INTERACTIONS]], 2);
#if DIM == 3
            dSqr += pow(z[i] - ghostParticles.z[nnlGhosts[j+i*MAX_NUM_GHOST_INTERACTIONS]], 2);
#endif
            double r = sqrt(dSqr);
            //double psij_xi = kernel(r, hi)/ghostParticles.omega[nnlGhosts[j + i * MAX_NUM_GHOST_INTERACTIONS]];
            double psij_xi = kernel(r, hi)/omega[i];

            xjGhost[0] = ghostParticles.x[nnlGhosts[j+i*MAX_NUM_GHOST_INTERACTIONS]];
            xjGhost[1] = ghostParticles.y[nnlGhosts[j+i*MAX_NUM_GHOST_INTERACTIONS]];
#if DIM==3
            xj[2] = ghostParticles.z[nnlGhosts[j+i*MAX_NUM_GHOST_INTERACTIONS]];
#endif
            //if (i == 7){
            //    Logger(DEBUG) << "psij_xi = " << psij_xi << ", xjGhost = [" << xjGhost[0] << ", " << xjGhost[1] << "]"
            //                << "; xi = [" << xi[0] << ", " << xi[1] << "]";
            //}

            for (int alpha=0; alpha<DIM; ++alpha){
                for(int beta=0; beta<DIM; ++beta){
                    B[DIM*alpha+beta] += (xjGhost[alpha] - xi[alpha])*(xjGhost[beta] - xi[beta]) * psij_xi;
                }
            }
        }

        //Logger(DEBUG) << "E = " << "[" << B[0] << ", " << B[1] << ", " << B[2] << ", " << B[3] << "]";

        // needed for sanity check of matrix E
        double normE = 0;
        for (int alpha=0; alpha<DIM; ++alpha){
           for (int beta=0; beta<DIM; ++beta){
               normE += B[alpha*DIM+beta]*B[alpha*DIM+beta];
           }
        }

#if OUTPUT_CONDITION_NUMBER && DIM == 2
        {
            double a = B[0], b = B[1], c = B[3];
            double trace = a + c;
            double disc = sqrt((a - c) * (a - c) + 4.0 * b * b);
            lambdaMax[i] = 0.5 * (trace + disc);
            lambdaMin[i] = 0.5 * (trace - disc);
            double ev_x = b;
            double ev_y = lambdaMin[i] - a;
            double norm = sqrt(ev_x * ev_x + ev_y * ev_y);
            if (norm > 0.) {
                eigenvecMin[i][0] = ev_x / norm;
                eigenvecMin[i][1] = ev_y / norm;
            } else {
                eigenvecMin[i][0] = 1.;
                eigenvecMin[i][1] = 0.;
            }
        }
#endif

        helper.inverseMatrix(B, DIM);

        double normB = 0;
        for (int alpha=0; alpha<DIM; ++alpha){
           for (int beta=0; beta<DIM; ++beta){
               normB += B[alpha*DIM+beta]*B[alpha*DIM+beta];
           }
        }

        // Check whether Matrix E is ill-conditioned
        double Ncond = 1./DIM * sqrt(normE*normB);
#if OUTPUT_CONDITION_NUMBER
        conditionNumber[i] = Ncond;
#endif
        // if (Ncond > 1){
        //     Logger(DEBUG) << "Ncond@" << i << " = " << Ncond;
        // }
        //if (i == 7) {
        //    Logger(DEBUG) << "B = " << "[" << B[0] << ", " << B[1] << ", " << B[2] << ", " << B[3] << "]";
        //}
        //Logger(DEBUG) << "noi[" << i << "] = " << noi[i] << ", noiGhosts[" << i << "] = " << noiGhosts[i];
        //exit(7);

        for (int j=0; j<noi[i]; ++j) {
            iP = nnl[j+i*MAX_NUM_INTERACTIONS];
            double dSqr = pow(x[i] - x[iP], 2)
                          + pow(y[i] - y[iP], 2);
#if DIM == 3
            dSqr += pow(z[i] - z[iP], 2);
#endif
            double r = sqrt(dSqr);
            //double psij_xi = kernel(r, hi) / omega[nnl[j + i * MAX_NUM_INTERACTIONS]];
            double psij_xi = kernel(r, hi) / omega[i];

            xj[0] = x[iP];
            xj[1] = y[iP];

            for (int alpha = 0; alpha < DIM; ++alpha) {
                //psijTilde_xi[nnl[j + i * MAX_NUM_INTERACTIONS]+i*MAX_NUM_INTERACTIONS][alpha] = 0.;
                psijTilde_xi[j + i * MAX_NUM_INTERACTIONS][alpha] = 0.;
                for (int beta = 0; beta < DIM; ++beta) {
                    //psijTilde_xi[nnl[j + i * MAX_NUM_INTERACTIONS]+i*MAX_NUM_INTERACTIONS][alpha] += B[alpha * DIM + beta] * (xj[beta] - xi[beta]) * psij_xi;
                    psijTilde_xi[j + i * MAX_NUM_INTERACTIONS][alpha] += B[alpha * DIM + beta] * (xj[beta] - xi[beta]) * psij_xi;
                }

                //if(i == 86){
                //    Logger(DEBUG) << "psijTilde_xi[" << alpha << "]@" << nnl[j + i * MAX_NUM_INTERACTIONS] << " = " << psijTilde_xi[nnl[j + i * MAX_NUM_INTERACTIONS]][alpha];
                //}
            }
        }

        for (int j=0; j<noiGhosts[i]; ++j) {

            double dSqr = pow(x[i] - ghostParticles.x[nnlGhosts[j + i * MAX_NUM_GHOST_INTERACTIONS]], 2)
                          + pow(y[i] - ghostParticles.y[nnlGhosts[j + i * MAX_NUM_GHOST_INTERACTIONS]], 2);
#if DIM == 3
            dSqr += pow(z[i] - ghostParticles.z[nnlGhosts[j+i*MAX_NUM_GHOST_INTERACTIONS]], 2);
#endif
            double r = sqrt(dSqr);
            //double psij_xi = kernel(r, hi) / ghostParticles.omega[nnlGhosts[j+i*MAX_NUM_GHOST_INTERACTIONS]];
            double psij_xi = kernel(r, hi) / omega[i];

            xjGhost[0] = ghostParticles.x[nnlGhosts[j+i*MAX_NUM_GHOST_INTERACTIONS]];
            xjGhost[1] = ghostParticles.y[nnlGhosts[j+i*MAX_NUM_GHOST_INTERACTIONS]];

            for (int alpha = 0; alpha < DIM; ++alpha) {
                //psijTilde_xiGhosts[nnlGhosts[j + i * MAX_NUM_GHOST_INTERACTIONS]+i*MAX_NUM_GHOST_INTERACTIONS][alpha] = 0.;
                psijTilde_xiGhosts[j + i * MAX_NUM_GHOST_INTERACTIONS][alpha] = 0.;
                for (int beta = 0; beta < DIM; ++beta) {
                    //psijTilde_xiGhosts[nnlGhosts[j + i * MAX_NUM_GHOST_INTERACTIONS]+i*MAX_NUM_GHOST_INTERACTIONS][alpha] += B[alpha * DIM + beta] * (xjGhost[beta] - xi[beta]) * psij_xi;
                    psijTilde_xiGhosts[j + i * MAX_NUM_GHOST_INTERACTIONS][alpha] += B[alpha * DIM + beta] * (xjGhost[beta] - xi[beta]) * psij_xi;
                }
                //if(i == 86){
                //    Logger(DEBUG) << "psijTildeGhost_xi[" << alpha << "]@" << nnlGhosts[j + i * MAX_NUM_GHOST_INTERACTIONS] << " = " << ghostParticles.psijTilde_xi[nnlGhosts[j + i * MAX_NUM_GHOST_INTERACTIONS]][alpha];
                //}
            }
        }
    }
}

void Particles::gradient(double *f, double (*grad)[DIM], double *fGhost, const Particles &ghostParticles){
    for (int i=0; i<N; ++i) {
        for (int alpha = 0; alpha < DIM; ++alpha) {
            grad[i][alpha] = 0;
        }

        // Hopkins 2015 MFM <-> SPH blend weight from conditionNumber[i].
        double wSPH = 0.;
        if (COND_MAX_FOR_GRADIENT > 0.){
            const double c = conditionNumber[i];
            if (COND_BLEND_HI > COND_MAX_FOR_GRADIENT){
                if      (c <= COND_MAX_FOR_GRADIENT) wSPH = 0.;
                else if (c >= COND_BLEND_HI)         wSPH = 1.;
                else wSPH = (c - COND_MAX_FOR_GRADIENT)
                          / (COND_BLEND_HI - COND_MAX_FOR_GRADIENT);
            } else {
                wSPH = (c > COND_MAX_FOR_GRADIENT) ? 1. : 0.;
            }
        }
        const double wMFM = 1. - wSPH;

        for (int j = 0; j < noi[i]; ++j) {
            int jIdx = nnl[j + i * MAX_NUM_INTERACTIONS];
            if (COND_MAX_NEIGHBOR_FOR_GRADIENT > 0. && conditionNumber[jIdx] > COND_MAX_NEIGHBOR_FOR_GRADIENT) {
                continue;
            }
            if (noi[jIdx] < MIN_NOI_FOR_NEIGHBOR_USE) {
                continue;
            }
            if (wMFM > 0.){
                for (int alpha = 0; alpha < DIM; ++alpha) {
                    grad[i][alpha] += wMFM * (f[jIdx] - f[i])
                                      * psijTilde_xi[j + i * MAX_NUM_INTERACTIONS][alpha];
                }
            }
            if (wSPH > 0.){
                const double dx = x[i] - x[jIdx];
                const double dy = y[i] - y[jIdx];
                double r2 = dx*dx + dy*dy;
#if DIM == 3
                const double dz = z[i] - z[jIdx];
                r2 += dz*dz;
#endif
                const double r = sqrt(r2);
                if (r <= 0.) continue;
                const double dW = Kernel::WDr(r, sml[i]);
                const double w  = wSPH * m[jIdx] * (f[jIdx] - f[i]) * dW / (r * rho[i]);
                grad[i][0] += w * dx;
                grad[i][1] += w * dy;
#if DIM == 3
                grad[i][2] += w * dz;
#endif
            }
        }

        for (int j = 0; j < noiGhosts[i]; ++j) {
            int gIdx = nnlGhosts[j + i * MAX_NUM_GHOST_INTERACTIONS];
            if (wMFM > 0.){
                for (int alpha = 0; alpha < DIM; ++alpha) {
                    grad[i][alpha] += wMFM * (fGhost[gIdx] - f[i])
                                      * psijTilde_xiGhosts[j + i * MAX_NUM_GHOST_INTERACTIONS][alpha];
                }
            }
            if (wSPH > 0.){
                const double dx = x[i] - ghostParticles.x[gIdx];
                const double dy = y[i] - ghostParticles.y[gIdx];
                double r2 = dx*dx + dy*dy;
#if DIM == 3
                const double dz = z[i] - ghostParticles.z[gIdx];
                r2 += dz*dz;
#endif
                const double r = sqrt(r2);
                if (r <= 0.) continue;
                const double dW = Kernel::WDr(r, sml[i]);
                const double w  = wSPH * ghostParticles.m[gIdx] * (fGhost[gIdx] - f[i])
                                  * dW / (r * rho[i]);
                grad[i][0] += w * dx;
                grad[i][1] += w * dy;
#if DIM == 3
                grad[i][2] += w * dz;
#endif
            }
        }

        //Logger(DEBUG) << "        > grad[" << i << "] = [" << grad[i][0] << ", "
        //          << grad[i][1] << "]";

        //for (int alpha=0; alpha < DIM; ++alpha ){
        //    if(i == 86) Logger(DEBUG) << "grad[" << i << "][" << alpha << "] = " << grad[i][alpha];
        //}
    }
}

void Particles::compEffectiveFace(const Particles &ghostParticles){
    for (int i=0; i<N; ++i){

        //Logger(DEBUG) << "      > i = " << i << ", x = [" << x[i] << ", " << y[i] <<
        //          "], omega = " << omega[i];

        for (int j=0; j<noiGhosts[i]; ++j){
            int ii = i*MAX_NUM_GHOST_INTERACTIONS+j;
            int ji = nnlGhosts[ii]; // index i of particle j
            const int p = ghostParticles.parent[ji];
            // Image-aware reverse match (see solveRiemannProblems): pick the
            // ghost-of-i in p's list closest to x[i]-delta, not the first parent
            // match, so the face is antisymmetric across multi-image seams.
            const double tgtX = x[i] - (ghostParticles.x[ji] - x[p]);
            const double tgtY = y[i] - (ghostParticles.y[ji] - y[p]);
            int ij = -1;
            double best = std::numeric_limits<double>::max();
            for (int ijc=0; ijc<noiGhosts[p]; ++ijc){
                const int gg = nnlGhosts[ijc+p*MAX_NUM_GHOST_INTERACTIONS];
                if (ghostParticles.parent[gg] != i) continue;
                double d = pow(ghostParticles.x[gg]-tgtX, 2)
                         + pow(ghostParticles.y[gg]-tgtY, 2);
                if (d < best){ best = d; ij = ijc; }
            }

            //Logger(DEBUG) << "j = " << j << ", ii = " << ii << ", ji = " << ji
            //          << ", ij = " << ij;
            //Logger(DEBUG) << "omegaGhost[ji] = " << ghostParticles.omega[ji]
            //          << ", psijTilde_xiGhost[ii] = [" << psijTilde_xiGhosts[ii][0] << ", " << psijTilde_xiGhosts[ii][1]
            //          << "], psijTilde_xiGhost[ji] = [" << psijTilde_xiGhosts[ij+ghostParticles.parent[ji]*MAX_NUM_GHOST_INTERACTIONS][0]
            //          << ", " << psijTilde_xiGhosts[ij+ghostParticles.parent[ji]*MAX_NUM_GHOST_INTERACTIONS][1] << "]";

            for (int alpha=0; alpha<DIM; ++alpha){
                // Drop the reverse term if p has no ghost of i (non-reciprocal):
                // indexing a stale slot would corrupt the face.
                const double rev = (ij < 0) ? 0.
                    : 1./ghostParticles.omega[ji]
                      * psijTilde_xiGhosts[ij+p*MAX_NUM_GHOST_INTERACTIONS][alpha];
                AijGhosts[ii][alpha] = 1./omega[i]*psijTilde_xiGhosts[ii][alpha] - rev;
            }
        }
    }
}

void Particles::compRiemannStatesLR(const double &dt,
                                  const Particles &ghostParticles){

    for (int i=0; i<N; ++i){
        double xij[DIM];
        //double vFrame[DIM];
        // helper vectors
        double xijxi[DIM], xjxi[DIM], xijxj[DIM];

        for (int jn=0; jn<noiGhosts[i]; ++jn){

            int j = nnlGhosts[i*MAX_NUM_GHOST_INTERACTIONS+jn];

            xjxi[0] = ghostParticles.x[j] - x[i];
            xjxi[1] = ghostParticles.y[j] - y[i];

            //xij[0] = x[i] + sml[i]/2. * xjxi[0];
            //xij[1] = y[i] + sml[i]/2. * xjxi[1];


#if FIRST_ORDER_QUAD_POINT
            xij[0] = (x[i] + ghostParticles.x[j])/2.;
            xij[1] = (y[i] + ghostParticles.y[j])/2.;

            xijxj[0] = .5*(x[i] - ghostParticles.x[j]);
            xijxj[1] = .5*(y[i] - ghostParticles.y[j]);

            xijxi[0] = .5*(ghostParticles.x[j] - x[i]);
            xijxi[1] = .5*(ghostParticles.y[j] - y[i]);

#else // !FIRST_ORDER_QUAD_POINT
            xij[0] = x[i] + sml[i]/4. * xjxi[0];
            xij[1] = y[i] + sml[i]/4. * xjxi[1];
            xijxi[0] = xij[0] - x[i];
            xijxi[1] = xij[1] - y[i];

            xijxj[0] = xij[0] - ghostParticles.x[j];
            xijxj[1] = xij[1] - ghostParticles.y[j];
#endif // FIRST_ORDER_QUAD_POINT

#if DIM==3
            xjxi[2] = ghostParticles.z[j] - z[i];
            //xij[2] = z[i] + sml[i]/2. * xjxi[2];
            xij[2] = z[i] + sml[i]/4. * xjxi[2];
            xijxi[2] = xij[2] - z[i];
            xijxj[2] = xij[2] - ghostParticles.z[j];
            // TODO: add FIRST_ORDER_QUAD_POINT
#endif // DIM == 3

            int iW = i*MAX_NUM_GHOST_INTERACTIONS+jn;
#if !MOVE_PARTICLES
            vFrameGhosts[iW][0] = 0.;
            vFrameGhosts[iW][1] = 0.;
#if DIM==3
            vFrameGhosts[iW][2] = 0.;
#endif // DIM == 3
#else // !MOVE_PARTICLES
#if FIRST_ORDER_QUAD_POINT
            vFrameGhosts[iW][0] = (vx[i] + ghostParticles.vx[j])/2.;
            vFrameGhosts[iW][1] = (vy[i] + ghostParticles.vy[j])/2.;
#if DIM==3
            vFrameGhosts[iW][2] = (vz[i] + ghostParticles.vz[j])/2.;
#endif // DIM == 3
#else // !FIRST_ORDER_QUAD_POINT
            double dotProd = Helper::dotProduct(xijxi, xjxi);
            double dSqr = Helper::dotProduct(xjxi, xjxi);

            vFrameGhosts[iW][0] = vx[i] + (ghostParticles.vx[j]-vx[i]) * dotProd/dSqr;
            vFrameGhosts[iW][1] = vy[i] + (ghostParticles.vy[j]-vy[i]) * dotProd/dSqr;
#if DIM==3
            vFrameGhosts[iW][2] = vz[i] + (ghostParticles.vz[j]-vz[i]) * dotProd/dSqr;
#endif // DIM == 3
#endif // FIRST_ORDER_QUAD_POINT
#endif // !MOVE_PARTICLES
            /*if(i == 394){
                Logger(DEBUG) << "vFrame[0] = " << vFrameGhosts[iW][0]
                              << ", vFrame[1] = " << vFrameGhosts[iW][1]
                              << ", rhoGrad[i][0] = " << rhoGrad[i][0]
                              << ", rhoGrad[i][1] = " << rhoGrad[i][1]
                              << ", rho[i] = " << rho[i]
                              << ", xijxi = [" << xijxi[0] << ", " << xijxi[1]
                              << "], xij = [" << xij[0] << ", " << xij[1]
                              << "], xj = [" << ghostParticles.x[j] << ", " << ghostParticles.y[j] << "] @" << j;
                //exit(5);
            }*/
            // boost frame to effective face
            WijRGhosts[iW][0] = rho[i];
            WijLGhosts[iW][0] = ghostParticles.rho[j];
            WijRGhosts[iW][1] = P[i];
            WijLGhosts[iW][1] = ghostParticles.P[j];
            WijRGhosts[iW][2] = vx[i] - vFrameGhosts[iW][0];
            WijLGhosts[iW][2] = ghostParticles.vx[j] - vFrameGhosts[iW][0];
            WijRGhosts[iW][3] = vy[i] - vFrameGhosts[iW][1];
            WijLGhosts[iW][3] = ghostParticles.vy[j] - vFrameGhosts[iW][1];
#if DIM == 3
            WijRGhosts[iW][4] = vz[i] - vFrameGhosts[iW][2];
            WijLGhosts[iW][4] = ghostParticles.vz[j] - vFrameGhosts[iW][2];
#endif // DIM == 3

            // reconstruction at effective face
            WijRGhosts[iW][0] += Helper::dotProduct(rhoGrad[i], xijxi);
            WijLGhosts[iW][0] += Helper::dotProduct(ghostParticles.rhoGrad[j], xijxj);
            WijRGhosts[iW][1] += Helper::dotProduct(PGrad[i], xijxi);
            WijLGhosts[iW][1] += Helper::dotProduct(ghostParticles.PGrad[j], xijxj);
            WijRGhosts[iW][2] += Helper::dotProduct(vxGrad[i], xijxi);
            WijLGhosts[iW][2] += Helper::dotProduct(ghostParticles.vxGrad[j], xijxj);
            WijRGhosts[iW][3] += Helper::dotProduct(vyGrad[i], xijxi);
            WijLGhosts[iW][3] += Helper::dotProduct(ghostParticles.vyGrad[j], xijxj);
#if DIM==3
            WijRGhosts[iW][4] += Helper::dotProduct(vzGrad[i], xijxi);
            WijLGhosts[iW][4] += Helper::dotProduct(ghostParticles.vzGrad[j], xijxj);
#endif

            //if (i == 1843 && iW == 28){
            //    Logger(DEBUG) << "        j = " << j
            //                  << ", rhoL = " << WijLGhosts[iW][0] << ", rhoR = " << WijRGhosts[iW][0]
            //                  << ", uL = " << WijLGhosts[iW][2] << ", uR = " << WijRGhosts[iW][2]
            //                  << ", PL = " << WijLGhosts[iW][1] << ", PR = " << WijRGhosts[iW][1];
            //}

            // predict half a timestep
            double viDiv = vxGrad[i][0] + vyGrad[i][1];
            double vjDiv = ghostParticles.vxGrad[j][0] + ghostParticles.vyGrad[j][1];
#if DIM==3
            viDiv += vzGrad[i][2];
            vjDiv += ghostParticles.vzGrad[j][2];
#endif // DIM == 3

#if EOS == 1 || EOS == 2
            double Ki = MeshlessEOS->EOSBulkModulus(matId[i], rho[i], P[i]);
            double Kj = MeshlessEOS->EOSBulkModulus(ghostParticles.matId[j], ghostParticles.rho[j], ghostParticles.P[j]);
#else
            double Ki = MeshlessEOS->EOSBulkModulus(rho[i], P[i]);
            double Kj = MeshlessEOS->EOSBulkModulus(ghostParticles.rho[j], ghostParticles.P[j]);
#endif

            WijRGhosts[iW][0] -= dt/2. * (rho[i] * viDiv + (vx[i]-vFrameGhosts[iW][0])*rhoGrad[i][0] + (vy[i]-vFrameGhosts[iW][1])*rhoGrad[i][1]);
            WijLGhosts[iW][0] -= dt/2. * (ghostParticles.rho[j] * vjDiv + (ghostParticles.vx[j]-vFrameGhosts[iW][0])*ghostParticles.rhoGrad[j][0] + (ghostParticles.vy[j]-vFrameGhosts[iW][1])*ghostParticles.rhoGrad[j][1]);
            WijRGhosts[iW][1] -= dt/2. * (Ki * viDiv + (vx[i]-vFrameGhosts[iW][0])*PGrad[i][0] + (vy[i]-vFrameGhosts[iW][1])*PGrad[i][1]);
            WijLGhosts[iW][1] -= dt/2. * (Kj * vjDiv + (ghostParticles.vx[j]-vFrameGhosts[iW][0])*ghostParticles.PGrad[j][0] + (ghostParticles.vy[j]-vFrameGhosts[iW][1])*ghostParticles.PGrad[j][1]);
            WijRGhosts[iW][2] -= dt/2. * (PGrad[i][0]/rho[i] + (vx[i] - vFrameGhosts[iW][0])*vxGrad[i][0] + (vy[i] - vFrameGhosts[iW][1])*vxGrad[i][1]);
            WijLGhosts[iW][2] -= dt/2. * (ghostParticles.PGrad[j][0]/ghostParticles.rho[j] + (ghostParticles.vx[j]-vFrameGhosts[iW][0])*ghostParticles.vxGrad[j][0] + (ghostParticles.vy[j]-vFrameGhosts[iW][1])*ghostParticles.vxGrad[j][1]);
            WijRGhosts[iW][3] -= dt/2. * (PGrad[i][1]/rho[i] + (vx[i] - vFrameGhosts[iW][0])*vyGrad[i][0] + (vy[i] - vFrameGhosts[iW][1])*vyGrad[i][1]);
            WijLGhosts[iW][3] -= dt/2. * (ghostParticles.PGrad[j][1]/ghostParticles.rho[j] + (ghostParticles.vx[j]-vFrameGhosts[iW][0])*ghostParticles.vyGrad[j][0] + (ghostParticles.vy[j]-vFrameGhosts[iW][1])*ghostParticles.vyGrad[j][1]);
#if ELASTIC
            // Elastic source: dv/dt += (1/ρ) ∇·S
            WijRGhosts[iW][2] += dt/2. * (SxxGrad[i][0] + SxyGrad[i][1]) / rho[i];
            WijLGhosts[iW][2] += dt/2. * (ghostParticles.SxxGrad[j][0] + ghostParticles.SxyGrad[j][1]) / ghostParticles.rho[j];
            WijRGhosts[iW][3] += dt/2. * (SxyGrad[i][0] + SyyGrad[i][1]) / rho[i];
            WijLGhosts[iW][3] += dt/2. * (ghostParticles.SxyGrad[j][0] + ghostParticles.SyyGrad[j][1]) / ghostParticles.rho[j];
#endif // ELASTIC
#if DIM==3
            // TODO: update below for 3D
            WijRGhosts[iW][4] -= dt/2. * PGrad[i][2]/rho[i];
            WijLGhosts[iW][4] -= dt/2. * ghostParticles.PGrad[j][2]/ghostParticles.rho[j];
#if ELASTIC
            WijRGhosts[iW][2] += dt/2. * SxzGrad[i][2] / rho[i];
            WijLGhosts[iW][2] += dt/2. * ghostParticles.SxzGrad[j][2] / ghostParticles.rho[j];
            WijRGhosts[iW][3] += dt/2. * SyzGrad[i][2] / rho[i];
            WijLGhosts[iW][3] += dt/2. * ghostParticles.SyzGrad[j][2] / ghostParticles.rho[j];
            WijRGhosts[iW][4] += dt/2. * (SxzGrad[i][0] + SyzGrad[i][1] + SzzGrad[i][2]) / rho[i];
            WijLGhosts[iW][4] += dt/2. * (ghostParticles.SxzGrad[j][0] + ghostParticles.SyzGrad[j][1] + ghostParticles.SzzGrad[j][2]) / ghostParticles.rho[j];
#endif // ELASTIC
#endif // DIM == 3

            if (DENSITY_FLOOR > 0.){
                if (WijRGhosts[iW][0] < DENSITY_FLOOR) WijRGhosts[iW][0] = DENSITY_FLOOR;
                if (WijLGhosts[iW][0] < DENSITY_FLOOR) WijLGhosts[iW][0] = DENSITY_FLOOR;
            }

        }
    }
}

/// debugging function to printout the number of (ghost-)interactions
void Particles::printNoi(){
    for (int i=0; i<N; ++i) {
        Logger(DEBUG) << "          > i = " << i << ", x = [" << x[i] << ", " << y[i] <<
                  "], : noi = " << noi[i] <<
                  " + " << noiGhosts[i] << " = " << noi[i] + noiGhosts[i];
    }
}

/// debugging function dumping nearest neighbors to file
void Particles::dumpNNL(std::string filename, const Particles &ghostParticles){
    // open output file
    HighFive::File h5File { filename, HighFive::File::ReadWrite |
                                      HighFive::File::Create |
                                      HighFive::File::Truncate };

    for (int i=0; i<N; ++i){

        int noiTot = noi[i] + noiGhosts[i];

        std::vector<int> nnList(noiTot);

        std::vector<std::vector<double>> nnlPrtcls {};
        std::vector<std::vector<double>> nnlAij {};
        std::vector<std::vector<double>> nnl_vFrame {};
        std::vector<std::vector<double>> nnlWijR {};
        std::vector<std::vector<double>> nnlWijL {};

        std::vector<size_t> dataSpaceDims(2);
        dataSpaceDims[0] = std::size_t(noiTot);
        dataSpaceDims[1] = DIM;

        for (int j=0; j<noi[i]; ++j){

            int ii = j+i*MAX_NUM_INTERACTIONS;

            nnList[j] = nnl[ii];

            nnlPrtcls.push_back(std::vector<double>(DIM));
            nnlPrtcls[j][0] = x[nnl[ii]];
            nnlPrtcls[j][1] = y[nnl[ii]];
#if DIM == 3
            nnlPrtcls[j][2] = z[nnl[ii]];
#endif
            nnlAij.push_back(std::vector<double>(DIM));
            nnlAij[j][0] = Aij[ii][0];
            nnlAij[j][1] = Aij[ii][1];
#if DIM == 3
            nnlAij[j][2] = Aij[ii][2];
#endif

            nnl_vFrame.push_back(std::vector<double>(DIM));
            nnl_vFrame[j][0] = vFrame[ii][0];
            nnl_vFrame[j][1] = vFrame[ii][1];
#if DIM == 3
            nnl_vFrame[j][2] = vFrame[ii][2];
#endif
        }

        for (int j=0; j<noiGhosts[i]; ++j){

            int ii = j+i*MAX_NUM_GHOST_INTERACTIONS;

            nnList[j+noi[i]] = nnlGhosts[ii];

            nnlPrtcls.push_back(std::vector<double>(DIM));
            nnlPrtcls[j+noi[i]][0] = ghostParticles.x[nnlGhosts[ii]];
            nnlPrtcls[j+noi[i]][1] = ghostParticles.y[nnlGhosts[ii]];
#if DIM == 3
            nnlPrtcls[j+noi[i]][2] = ghostParticles.z[nnlGhosts[ii]];
#endif
            nnlAij.push_back(std::vector<double>(DIM));
            nnlAij[j+noi[i]][0] = AijGhosts[ii][0];
            nnlAij[j+noi[i]][1] = AijGhosts[ii][1];
#if DIM == 3
            nnlAij[j+noi[i]][2] = AijGhosts[ii][2];
#endif

            nnl_vFrame.push_back(std::vector<double>(DIM));
            nnl_vFrame[j+noi[i]][0] = vFrameGhosts[ii][0];
            nnl_vFrame[j+noi[i]][1] = vFrameGhosts[ii][1];
#if DIM == 3
            nnl_vFrame[j+noi[i]][2] = vFrameGhosts[ii][2];
#endif

        }

        HighFive::DataSet nnListDataSet = h5File.createDataSet<int>("/nnl" +std::to_string(i),
                                                                    HighFive::DataSpace(noiTot));
        nnListDataSet.write(nnList);

        HighFive::DataSet nnlDataSet = h5File.createDataSet<double>("/nnlPrtcls" + std::to_string(i),
                                                                    HighFive::DataSpace(dataSpaceDims));
        nnlDataSet.write(nnlPrtcls);

        HighFive::DataSet AijDataSet = h5File.createDataSet<double>("/Aij" + std::to_string(i),
                                                                    HighFive::DataSpace(dataSpaceDims));
        AijDataSet.write(nnlAij);

        HighFive::DataSet vFrameDataSet = h5File.createDataSet<double>("/vFrame" + std::to_string(i),
                                                                    HighFive::DataSpace(dataSpaceDims));
        vFrameDataSet.write(nnl_vFrame);

    }
}
#endif
void Particles::dumpNNL(std::string filename){
    HighFive::File h5File { filename, HighFive::File::ReadWrite |
                                      HighFive::File::Create |
                                      HighFive::File::Truncate };

    for (int i=0; i<N; ++i){

        std::vector<int> nnList(noi[i]);

        std::vector<std::vector<double>> nnlPrtcls {};
        std::vector<std::vector<double>> nnlAij {};
        std::vector<std::vector<double>> nnl_vFrame {};

        std::vector<size_t> dataSpaceDims(2);
        dataSpaceDims[0] = std::size_t(noi[i]);
        dataSpaceDims[1] = DIM;

        for (int j=0; j<noi[i]; ++j){

            int ii = j+i*MAX_NUM_INTERACTIONS;

            nnList[j] = nnl[ii];

            nnlPrtcls.push_back(std::vector<double>(DIM));
            nnlPrtcls[j][0] = x[nnl[ii]];
            nnlPrtcls[j][1] = y[nnl[ii]];
#if DIM == 3
            nnlPrtcls[j][2] = z[nnl[ii]];
#endif
            nnlAij.push_back(std::vector<double>(DIM));
            nnlAij[j][0] = Aij[ii][0];
            nnlAij[j][1] = Aij[ii][1];
#if DIM == 3
            nnlAij[j][2] = Aij[ii][2];
#endif

            nnl_vFrame.push_back(std::vector<double>(DIM));
            nnl_vFrame[j][0] = vFrame[ii][0];
            nnl_vFrame[j][1] = vFrame[ii][1];
#if DIM == 3
            nnl_vFrame[j][2] = vFrame[ii][2];
#endif
        }

        HighFive::DataSet nnListDataSet = h5File.createDataSet<int>("/nnl" +std::to_string(i),
                                                                    HighFive::DataSpace(noi[i]));
        nnListDataSet.write(nnList);

        HighFive::DataSet nnlDataSet = h5File.createDataSet<double>("/nnlPrtcls" + std::to_string(i),
                                                                    HighFive::DataSpace(dataSpaceDims));
        nnlDataSet.write(nnlPrtcls);

        HighFive::DataSet AijDataSet = h5File.createDataSet<double>("/Aij" + std::to_string(i),
                                                                    HighFive::DataSpace(dataSpaceDims));
        AijDataSet.write(nnlAij);

        HighFive::DataSet vFrameDataSet = h5File.createDataSet<double>("/vFrame" + std::to_string(i),
                                                                    HighFive::DataSpace(dataSpaceDims));
        vFrameDataSet.write(nnl_vFrame);

    }
}

// TODO: remove below
//void Particles::move(const double &dt, Domain &domain){
//
//    for(int i=0; i<N; ++i) {
//
//        //Logger(DEBUG) << "v@" << i <<  " = [" << vx[i] << ", " << vy[i] << "]"
//        //            << ", x[n] = [" << x[i] << ", " << y[i] << "]";
//
//        x[i] = x[i] + vx[i] * dt;
//        y[i] = y[i] + vy[i] * dt;
//#if DIM == 3
//        z[i] = z[i] +vz[i] * dt;
//#endif
//#if PERIODIC_BOUNDARIES
//        if (x[i] < domain.bounds.minX) {
//            x[i] = domain.bounds.maxX - (domain.bounds.minX - x[i]);
//        } else if (domain.bounds.maxX <= x[i]) {
//            x[i] = domain.bounds.minX + (x[i] - domain.bounds.maxX);
//        }
//        if (y[i] < domain.bounds.minY) {
//            y[i] = domain.bounds.maxY - (domain.bounds.minY - y[i]);
//        } else if (domain.bounds.maxY <= y[i]) {
//            y[i] = domain.bounds.minY + (y[i] - domain.bounds.maxY);
//        }
//#if DIM ==3
//        if (z[i] < domain.bounds.minZ) {
//            z[i] = domain.bounds.maxZ - (domain.bounds.minZ - z[i]);
//        } else if (domain.bounds.maxZ <= z[i]) {
//            z[i] = domain.bounds.minZ + (z[i] - domain.bounds.maxZ);
//        }
//#endif
//#endif
//        //Logger(DEBUG) << "                               x[n+1] = ["
//        //          << x[i] << ", " << y[i] << "]";
//    }
//}

/// Sanity check functions
double Particles::sumVolume(){
    double V = 0.;
    for (int i=0; i<N; ++i){
        // V_i = m_i/rho_i, mirroring compEffectiveFace and GIZMO's Mass/Density.
        // Tracks rhoExplicit under EXPLICIT_VOL_INTEGRATION and reduces to
        // fce/omega (SURFACE_VOLCORR) or 1/omega otherwise -- see the comment in
        // compEffectiveFace for the equivalences.
        V += m[i] / rho[i];
    }
    return V;
}

double Particles::sumMass(){
    double M = 0.;
    for (int i=0; i<N; ++i){
        if (std::isnan(m[i])){
            Logger(WARN) << "!! m[" << i <<"] " << "is nan. !!";
#if DEBUG_LVL
            exit(6);
#endif
        }
        M += m[i];
    }
    return M;
}

double Particles::sumEnergy(){
    double E = 0.;
    for (int i=0; i<N; ++i){
#if DIM == 2
        E += m[i]*(u[i] + .5*(vx[i]*vx[i]+vy[i]*vy[i]));
#else
        E += m[i]*(u[i] + .5*(vx[i]*vx[i]+vy[i]*vy[i]+vz[i]*vz[i]));
#endif
    }
    //std::cout << "E is " <<  E << std::endl;
    return E;
}

double Particles::sumMomentumX(){
    double momX = 0.;
    for (int i=0; i<N; ++i){
        momX += m[i]*vx[i];
    }
    return momX;
}

double Particles::sumMomentumY(){
    double momY = 0.;
    for (int i=0; i<N; ++i){
        momY += m[i]*vy[i];
    }
    return momY;
}

#if DIM == 3
double Particles::sumMomentumZ(){
    double momZ = 0.;
    for (int i=0; i<N; ++i){
        momZ += m[i]*vz[i];
    }
    return momZ;
}
#endif

void Particles::checkFluxSymmetry(Particles *ghostParticles){
    for (int i=0; i<N; ++i){
        for (int j=0; j<noi[i]; ++j){
            int ii = i*MAX_NUM_INTERACTIONS+j; // interaction index i->j

            int ji = nnl[ii]; // index i of particle j
            // search neighbor i in nnl[] of j
            int ij;
            for(ij=0; ij<noi[ji]; ++ij){
                if (nnl[ij+ji*MAX_NUM_INTERACTIONS] == i) break;
            }
            int iij = ji*MAX_NUM_INTERACTIONS+ij; // interaction index j->i

            bool notSym = false;
            if (Fij[iij][0] + Fij[ii][0] > FLUX_SYM_TOL){
                Logger(WARN) << "  > Mass fluxes are NOT symmetric for " << i << " -> " << j
                             << " => Fij = " << Fij[ii][0] << ", Fji = " << Fij[iij][0];
                notSym = true;
            }
            if (Fij[iij][1] + Fij[ii][1] > FLUX_SYM_TOL){
                Logger(WARN) << "  > Energy fluxes are NOT symmetric for " << i << " -> " << j
                             << " => Fij = " << Fij[ii][1] << ", Fji = " << Fij[iij][1];
                notSym = true;
            }
            if (Fij[iij][2] + Fij[ii][2] > FLUX_SYM_TOL){
                Logger(WARN) << "  > Momentum fluxes (x) are NOT symmetric for " << i << " -> " << j
                             << " => Fij = " << Fij[ii][2] << ", Fji = " << Fij[iij][2];
                notSym = true;
            }
            if (Fij[iij][3] + Fij[ii][3] > FLUX_SYM_TOL){
                Logger(WARN) << "  > Momentum fluxes (y) are NOT symmetric for " << i << " -> " << j
                             << " => Fij = " << Fij[ii][3] << ", Fji = " << Fij[iij][3];
                notSym = true;
            }
#if DIM == 3
            if (Fij[iij][4] + Fij[ii][4] > FLUX_SYM_TOL){
                Logger(WARN) << "  > Momentum fluxes (z) are NOT symmetric for " << i << " -> " << j
                             << " => Fij = " << Fij[ii][4] << ", Fji = " << Fij[iij][4];
                notSym = true;
            }
#endif
            if (notSym){
                Logger(INFO) << "  > Aij = [" << Aij[ii][0] << ", " << Aij[ii][1] << "], Aji = ["
                             << Aij[iij][0] << ", " << Aij[iij][1] << "]";
            }
        }
#if PERIODIC_BOUNDARIES
        for (int j=0; j<noiGhosts[i]; ++j) {
            int ii = i * MAX_NUM_GHOST_INTERACTIONS + j; // interaction index i->j

            int ji = nnlGhosts[ii]; // index i of particle j
            const int p = ghostParticles->parent[ji];
            // Image-aware reverse match (see solveRiemannProblems): the
            // reciprocal interaction is p with the image of i at x[i]-delta.
            const double tgtX = x[i] - (ghostParticles->x[ji] - x[p]);
            const double tgtY = y[i] - (ghostParticles->y[ji] - y[p]);
            int ij = -1;
            double best = std::numeric_limits<double>::max();
            for (int ijc = 0; ijc < noiGhosts[p]; ++ijc) {
                const int gg = nnlGhosts[ijc + p * MAX_NUM_GHOST_INTERACTIONS];
                if (ghostParticles->parent[gg] != i) continue;
                double d = pow(ghostParticles->x[gg]-tgtX, 2)
                         + pow(ghostParticles->y[gg]-tgtY, 2);
                if (d < best){ best = d; ij = ijc; }
            }
            if (ij < 0) continue; // non-reciprocal: no partner to compare against
            int iij = ij + p * MAX_NUM_GHOST_INTERACTIONS; // interaction index j->i

            bool notSym = false;
            if (FijGhosts[iij][0] + FijGhosts[ii][0] > FLUX_SYM_TOL) {
                Logger(WARN) << "  > Ghosts mass fluxes are NOT symmetric for " << i << " -> " << j
                             << " => Fij = " << FijGhosts[ii][0] << ", Fji = " << FijGhosts[iij][0];
                notSym = true;
            }
            if (FijGhosts[iij][1] + FijGhosts[ii][1] > FLUX_SYM_TOL) {
                Logger(WARN) << "  > Ghosts energy fluxes are NOT symmetric for " << i << " -> " << j
                             << " => Fij = " << FijGhosts[ii][1] << ", Fji = " << FijGhosts[iij][1];
                notSym = true;
            }
            if (FijGhosts[iij][2] + FijGhosts[ii][2] > FLUX_SYM_TOL) {
                Logger(WARN) << "  > Ghosts momentum fluxes (x) are NOT symmetric for " << i << " -> " << j
                             << " => Fij = " << FijGhosts[ii][2] << ", Fji = " << FijGhosts[iij][2];
                notSym = true;
            }
            if (FijGhosts[iij][3] + FijGhosts[ii][3] > FLUX_SYM_TOL) {
                Logger(WARN) << "  > Ghosts momentum fluxes (y) are NOT symmetric for " << i << " -> " << j
                             << " => Fij = " << FijGhosts[ii][3] << ", Fji = " << FijGhosts[iij][3];
                notSym = true;
            }
            if (notSym) {
                Logger(INFO) << "  > AijGhost = [" << AijGhosts[ii][0] << ", " << AijGhosts[ii][1] << "], Aji = ["
                             << AijGhosts[iij][0] << ", " << AijGhosts[iij][1] << "]";
            }
        }
#endif
    }
}

void Particles::dump2file(std::string filename, double simTime){
    // open output file
    HighFive::File h5File { filename, HighFive::File::ReadWrite |
                                      HighFive::File::Create |
                                      HighFive::File::Truncate };

    // dimensions for datasets containing vectors
    std::vector<size_t> dataSpaceDims(2);
    dataSpaceDims[0] = std::size_t(N); // number of particles
    dataSpaceDims[1] = DIM;

    // create datasets
    // TODO: Create a h5 object holding all meta data
    HighFive::DataSet timeDataSet = h5File.createDataSet<double>("/time", HighFive::DataSpace(1));
    HighFive::DataSet totalMassDataSet = h5File.createDataSet<double>("/totalMass", HighFive::DataSpace(1));
    HighFive::DataSet energyDataSet = h5File.createDataSet<double>("/energy", HighFive::DataSpace(1));
    HighFive::DataSet xMomDataSet = h5File.createDataSet<double>("/xMomentum", HighFive::DataSpace(1));
    HighFive::DataSet yMomDataSet = h5File.createDataSet<double>("/yMomentum", HighFive::DataSpace(1));
#if DIM == 3
    HighFive::DataSet zMomDataSet = h5File.createDataSet<double>("/zMomentum", HighFive::DataSpace(1));
#endif
    HighFive::DataSet rhoDataSet = h5File.createDataSet<double>("/rho", HighFive::DataSpace(N));
    HighFive::DataSet mDataSet = h5File.createDataSet<double>("/m", HighFive::DataSpace(N));
    HighFive::DataSet uDataSet = h5File.createDataSet<double>("/u", HighFive::DataSpace(N));
    HighFive::DataSet posDataSet = h5File.createDataSet<double>("/x", HighFive::DataSpace(dataSpaceDims));
    HighFive::DataSet velDataSet = h5File.createDataSet<double>("/v", HighFive::DataSpace(dataSpaceDims));
    HighFive::DataSet rhoGradDataSet = h5File.createDataSet<double>("/rhoGrad", HighFive::DataSpace(dataSpaceDims));
    HighFive::DataSet PDataSet = h5File.createDataSet<double>("/P", HighFive::DataSpace(N));
    HighFive::DataSet smlDataSet = h5File.createDataSet<double>("/sml", HighFive::DataSpace(N));
    HighFive::DataSet noiDataSet = h5File.createDataSet<int>("/noi", HighFive::DataSpace(N));
#if OUTPUT_CONDITION_NUMBER
    HighFive::DataSet condDataSet = h5File.createDataSet<double>("/conditionNumber", HighFive::DataSpace(N));
#if DIM == 2
    HighFive::DataSet lambdaMaxDataSet = h5File.createDataSet<double>("/lambdaMax", HighFive::DataSpace(N));
    HighFive::DataSet lambdaMinDataSet = h5File.createDataSet<double>("/lambdaMin", HighFive::DataSpace(N));
    HighFive::DataSet eigenvecMinDataSet = h5File.createDataSet<double>("/eigenvecMin", HighFive::DataSpace(dataSpaceDims));
#endif
#endif
#if SURFACE_VOLCORR
    HighFive::DataSet fceDataSet = h5File.createDataSet<double>("/fce", HighFive::DataSpace(N));
#endif
#if EXPLICIT_VOL_INTEGRATION
    HighFive::DataSet rhoExplicitDataSet = h5File.createDataSet<double>("/rhoExplicit", HighFive::DataSpace(N));
    HighFive::DataSet rhoKernelDataSet   = h5File.createDataSet<double>("/rhoKernel",   HighFive::DataSpace(N));
#endif
#if ELASTIC
    HighFive::DataSet SxxDataSet = h5File.createDataSet<double>("/Sxx", HighFive::DataSpace(N));
    HighFive::DataSet SxyDataSet = h5File.createDataSet<double>("/Sxy", HighFive::DataSpace(N));
    HighFive::DataSet SyyDataSet = h5File.createDataSet<double>("/Syy", HighFive::DataSpace(N));
#if DIM == 3
    HighFive::DataSet SxzDataSet = h5File.createDataSet<double>("/Sxz", HighFive::DataSpace(N));
    HighFive::DataSet SyzDataSet = h5File.createDataSet<double>("/Syz", HighFive::DataSpace(N));
    HighFive::DataSet SzzDataSet = h5File.createDataSet<double>("/Szz", HighFive::DataSpace(N));
#endif
#if FRAGMENTATION
    HighFive::DataSet damageDataSet      = h5File.createDataSet<double>("/damage", HighFive::DataSpace(N));
    HighFive::DataSet damageTotalDataSet = h5File.createDataSet<double>("/damageTotal", HighFive::DataSpace(N));
    HighFive::DataSet numActiveFlawsDataSet = h5File.createDataSet<int>("/numActiveFlaws", HighFive::DataSpace(N));
#endif
#endif

    // containers for particle data
    std::vector<double> timeVec({ simTime });
    std::vector<double> totalMassVec({ sumMass() });
    std::vector<double> energyVec({ sumEnergy() });
    std::vector<double> xMomVec({ sumMomentumX() });
    std::vector<double> yMomVec({ sumMomentumY() });
#if DIM ==3
    std::vector<double> zMomVec({ sumMomentumZ() });
#endif
    std::vector<double> rhoVec(rho, rho+N);
    std::vector<double> mVec(m, m+N);
    std::vector<double> uVec(u, u+N);
    std::vector<double> PVec(P, P+N);
    std::vector<int> noiVec(noi, noi+N);

    std::vector<std::vector<double>> posVec(N);
    std::vector<std::vector<double>> velVec(N);

    std::vector<std::vector<double>> rhoGradVec(N);
#if ELASTIC
    std::vector<double> SxxVec(Sxx, Sxx+N);
    std::vector<double> SxyVec(Sxy, Sxy+N);
    std::vector<double> SyyVec(Syy, Syy+N);
#if DIM == 3
    std::vector<double> SxzVec(Sxz, Sxz+N);
    std::vector<double> SyzVec(Syz, Syz+N);
    std::vector<double> SzzVec(Szz, Szz+N);
#endif
#if FRAGMENTATION
    std::vector<double> damageVec(damage, damage+N);
    std::vector<double> damageTotalVec(damageTotal, damageTotal+N);
    std::vector<int> numActiveFlawsVec(numActiveFlaws, numActiveFlaws+N);
#endif
#endif // ELASTIC

    // fill containers with data
    std::vector<double> posBuf(DIM);
    std::vector<double> velBuf(DIM);
    std::vector<double> rhoGradBuf(DIM);
    for(int i=0; i<N; ++i){
        // Logger(DEBUG) << " Pressure of particle i = " << i << " is P = " << P[i] << ", vector at i is " << PVec.at(i);
        //Logger(DEBUG) << "      > Dumping particle @"  << i;
        // position
        posBuf[0] = x[i];
        posBuf[1] = y[i];
#if DIM == 3
        posBuf[2] = z[i];
#endif
        posVec[i] = posBuf;

        // velocity
        velBuf[0] = vx[i];
        velBuf[1] = vy[i];
#if DIM == 3
        velBuf[2] = vz[i];
#endif
        velVec[i] = velBuf;

        // density gradient
        rhoGradBuf[0] = rhoGrad[i][0];
        rhoGradBuf[1] = rhoGrad[i][1];
#if DIM == 3
        rhoGradBuf[2] = rhoGrad[i][2];
#endif
        rhoGradVec[i] = rhoGradBuf;
    }
    // write data
    timeDataSet.write(timeVec); // dummy vec containing one element
    totalMassDataSet.write(totalMassVec);
    energyDataSet.write(energyVec);
    xMomDataSet.write(xMomVec);
    yMomDataSet.write(yMomVec);
#if DIM == 3
    zMomDataSet.write(zMomVec);
#endif
    rhoDataSet.write(rhoVec);
    mDataSet.write(mVec);
    uDataSet.write(uVec);
    PDataSet.write(PVec);
    {
        std::vector<double> smlVec(sml, sml+N);
        smlDataSet.write(smlVec);
    }
    noiDataSet.write(noi);
#if SURFACE_VOLCORR
    {
        std::vector<double> fceVec(fce, fce+N);
        fceDataSet.write(fceVec);
    }
#endif
#if EXPLICIT_VOL_INTEGRATION
    {
        std::vector<double> rhoExplicitVec(rhoExplicit, rhoExplicit+N);
        std::vector<double> rhoKernelVec(rhoKernel,     rhoKernel+N);
        rhoExplicitDataSet.write(rhoExplicitVec);
        rhoKernelDataSet.write(rhoKernelVec);
    }
#endif
#if OUTPUT_CONDITION_NUMBER
    std::vector<double> condVec(conditionNumber, conditionNumber+N);
    condDataSet.write(condVec);
#if DIM == 2
    std::vector<double> lambdaMaxVec(lambdaMax, lambdaMax+N);
    std::vector<double> lambdaMinVec(lambdaMin, lambdaMin+N);
    lambdaMaxDataSet.write(lambdaMaxVec);
    lambdaMinDataSet.write(lambdaMinVec);
    {
        std::vector<std::vector<double>> eigenvecMinVec(N);
        std::vector<double> evBuf(DIM);
        for (int i = 0; i < N; ++i) {
            evBuf[0] = eigenvecMin[i][0];
            evBuf[1] = eigenvecMin[i][1];
            eigenvecMinVec[i] = evBuf;
        }
        eigenvecMinDataSet.write(eigenvecMinVec);
    }
#endif
#endif
    posDataSet.write(posVec);
    velDataSet.write(velVec);
    rhoGradDataSet.write(rhoGradVec);

#if ELASTIC
    SxxDataSet.write(SxxVec);
    SxyDataSet.write(SxyVec);
    SyyDataSet.write(SyyVec);
#if DIM == 3
    SxzDataSet.write(SxzVec);
    SyzDataSet.write(SyzVec);
    SzzDataSet.write(SzzVec);
#endif
#if FRAGMENTATION
    damageDataSet.write(damageVec);
    damageTotalDataSet.write(damageTotalVec);
    numActiveFlawsDataSet.write(numActiveFlawsVec);
#endif
#endif
}
