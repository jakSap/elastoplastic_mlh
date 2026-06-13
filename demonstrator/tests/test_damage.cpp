// Unit test for the solid-failure support code:
//   * Helper::maxEigenvalueSym   (max tensile principal stress, real code)
//   * EquationOfState::EOSYoungModulus / EOSYieldStrength (real code)
//
// Uses the project's parameter.h (DIM / EOS / *_PLASTICITY as configured), so
// to exercise a particular plasticity model build it with that flag enabled.
//
// Build (from demonstrator/):  make test-damage

#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <vector>

#include "Helper.h"
#include "EquationOfState.h"

// Logger needs this global (normally defined in main.cpp).
structlog LOGCFG = {};

static int failures = 0;
static void check(bool ok, const char *what){
    std::printf("  %-48s %s\n", what, ok ? "OK" : "FAIL");
    if (!ok) ++failures;
}

// Build a symmetric matrix S = Q diag(l) Q^T and verify maxEigenvalueSym.
static void testEigen(){
    std::srand(7);
    double worst = 0.0;
    for (int t = 0; t < 20000; ++t){
#if DIM == 2
        const double a = 2*((double)std::rand()/RAND_MAX)-1;
        const double l0 = 4*((double)std::rand()/RAND_MAX)-2;
        const double l1 = 4*((double)std::rand()/RAND_MAX)-2;
        const double c = std::cos(a), s = std::sin(a);
        double S[4];
        S[0] = c*c*l0 + s*s*l1; S[3] = s*s*l0 + c*c*l1;
        S[1] = c*s*(l0-l1);     S[2] = S[1];
        const double lm = std::max(l0, l1);
#else
        double v[3][3];
        for (int i=0;i<3;i++) for (int j=0;j<3;j++) v[i][j]=2*((double)std::rand()/RAND_MAX)-1;
        auto dot=[&](int a,int b){return v[a][0]*v[b][0]+v[a][1]*v[b][1]+v[a][2]*v[b][2];};
        auto nrm=[&](int a){double n=std::sqrt(dot(a,a));for(int j=0;j<3;j++)v[a][j]/=n;};
        nrm(0);
        {double d=dot(1,0);for(int j=0;j<3;j++)v[1][j]-=d*v[0][j];} nrm(1);
        {double d0=dot(2,0),d1=dot(2,1);for(int j=0;j<3;j++)v[2][j]-=d0*v[0][j]+d1*v[1][j];} nrm(2);
        double l[3]={4*((double)std::rand()/RAND_MAX)-2,4*((double)std::rand()/RAND_MAX)-2,
                     4*((double)std::rand()/RAND_MAX)-2};
        double S[9];
        for(int i=0;i<3;i++)for(int j=0;j<3;j++){double sij=0;for(int k=0;k<3;k++)sij+=v[k][i]*l[k]*v[k][j];S[i*3+j]=sij;}
        const double lm = std::max(l[0], std::max(l[1], l[2]));
#endif
        const double e = std::fabs(Helper::maxEigenvalueSym(S) - lm);
        if (e > worst) worst = e;
    }
    std::printf("  maxEigenvalueSym (DIM=%d) maxAbsErr = %.3e\n", DIM, worst);
    check(worst < 1e-9, "maxEigenvalueSym matches known eigenvalues");
}

#if EOS == 2
static EquationOfState makeEOS(){
    TillotsonMaterial m {};
    m.rho0 = 1.0; m.A = 1.0; m.B = 0.0; m.u0 = 1.0; m.mu = 0.22;
    m.u_iv = 0.05; m.u_cv = 0.18;
    m.Y0 = 0.1; m.Y_M = 0.5; m.mu_i = 1.5; m.mu_d = 0.8;
    m.frictionAngle = 0.6; m.u_melt = 1.0;
    std::vector<TillotsonMaterial> mats { m };
    return EquationOfState(std::move(mats));
}

static void testEOS(){
    EquationOfState eos = makeEOS();
    // E = 9 K mu / (3K + mu) with K = A = 1, mu = 0.22.
    const double E = eos.EOSYoungModulus(0);
    const double Eref = 9.0*1.0*0.22/(3.0*1.0+0.22);
    check(std::fabs(E - Eref) < 1e-12, "EOSYoungModulus = 9Kmu/(3K+mu)");

#if PLASTICITY_MODEL_COUNT
    // Generic invariants that hold for any selected pressure-dependent model.
    const double Y_lowP  = eos.EOSYieldStrength(0, 0.0, 0.0, 0.0);
    const double Y_highP = eos.EOSYieldStrength(0, 1.0, 0.0, 0.0);
    check(Y_lowP >= 0.0 && Y_highP >= 0.0, "yield strength non-negative");
#if VON_MISES_PLASTICITY
    check(std::fabs(Y_lowP - 0.1) < 1e-12 && std::fabs(Y_highP - 0.1) < 1e-12,
          "von Mises yield constant = Y0");
#else
    check(Y_highP >= Y_lowP, "pressure-dependent yield increases with P");
#endif
#if COLLINS_PLASTICITY && !COLLINS_PLASTICITY_SIMPLE
    // Damage blends intact (D=0) and damaged (D=1) curves: Y(D=1) <= Y(D=0).
    const double Yi = eos.EOSYieldStrength(0, 1.0, 0.0, 0.0);
    const double Yd = eos.EOSYieldStrength(0, 1.0, 1.0, 0.0);
    check(Yd <= Yi + 1e-12, "Collins: fully damaged yield <= intact yield");
#endif
#else
    // No per-material model selected: a huge sentinel (no yielding).
    check(eos.EOSYieldStrength(0, 0.0, 0.0, 0.0) > 1e29, "no model => huge yield sentinel");
#endif
}
#endif // EOS == 2

int main(){
    std::printf("DIM=%d EOS=%d PLASTICITY_MODEL_COUNT=%d\n",
                DIM, EOS, PLASTICITY_MODEL_COUNT);
    testEigen();
#if EOS == 2
    testEOS();
#endif
    if (failures){ std::printf("FAILED (%d)\n", failures); return 1; }
    std::printf("All damage/plasticity unit tests passed.\n");
    return 0;
}
