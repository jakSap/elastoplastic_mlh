// Unit test P0a of the plasticity test plan (MA-Obsidian/2026-09-30-testplan-plastisch.md):
// the pointwise solid model in include/Plasticity.h, which Particles::integrateStressTensor
// calls -- 2D deviatoric rate (SOLID_2D_MODEL 0/1), J2, the yield return (YIELD_RETURN 0/1)
// and the plastic-strain increment -- plus a material-point driver built from them.
//
// Every mode is passed explicitly, so one build checks all of them regardless of the
// flags in parameter.h. Exit code = number of failed checks.
//
// Build + run (from demonstrator/):  make test-plasticity

#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <algorithm>

#include "Plasticity.h"

static int failures = 0;
static void check(bool ok, const char *what, double val = NAN){
    if (std::isnan(val)) std::printf("  %-66s %s\n", what, ok ? "OK" : "FAIL");
    else                 std::printf("  %-66s %s  (%.3e)\n", what, ok ? "OK" : "FAIL", val);
    if (!ok) ++failures;
}
static double urand(){ return 2.0*((double)std::rand()/RAND_MAX) - 1.0; }

// Material of the test plan (P.1.2): mu = 6/13, Y = 0.01.
static const double MU = 6.0/13.0, Y = 0.01;

// 1. Below the yield surface the return is the identity (bitwise: f == 1 exactly).
static void testIdentityBelowYield(){
    std::srand(11);
    bool ok = true;
    for (int t = 0; t < 100000; ++t){
        for (int dim = 2; dim <= 3; ++dim){
            const double J2max = Y*Y/dim;
            const double J2 = J2max * (0.5 + 0.5*urand());   // in [0, J2max]
            for (int mode = 0; mode <= 1; ++mode)
                if (Plasticity::yieldFactor(J2, Y, mode, dim) != 1.0) ok = false;
        }
        if (Plasticity::yieldFactor(Y*Y/2., Y, 0, 2) != 1.0) ok = false;   // on the surface
    }
    check(ok, "return is the identity for J2 <= Y^2/d (both modes, 2D+3D)");
}

// 2./3. Return maps on random stress states above the surface.
static void testReturnMaps(){
    std::srand(12);
    double err0_2 = 0., err0_3 = 0., err1 = 0., dirErr = 0.;
    for (int t = 0; t < 100000; ++t){
        const double s = Y * (1. + 20.*std::fabs(urand()));   // up to 21 Y
        // 2D, traceless (true 2D) and non-traceless states alike
        double sxx = s*urand(), syy = s*urand(), sxy = s*urand();
        double J2 = Plasticity::J2_2D(sxx, syy, sxy);
        double f0 = Plasticity::yieldFactor(J2, Y, 0, 2);
        if (f0 < 1.0){
            const double J2a = Plasticity::J2_2D(f0*sxx, f0*syy, f0*sxy);
            err0_2 = std::max(err0_2, std::fabs(std::sqrt(2.*J2a)/Y - 1.));
            // direction preserved: the scaled tensor is parallel to the trial
            const double n = std::sqrt(sxx*sxx + syy*syy + sxy*sxy);
            const double na = std::sqrt(f0*f0*(sxx*sxx + syy*syy + sxy*sxy));
            dirErr = std::max(dirErr, std::fabs((f0*sxx)/na - sxx/n));
            const double f1 = Plasticity::yieldFactor(J2, Y, 1, 2);
            const double ref = (Y*Y/2.) / J2;     // GIZMO: f_Y = Y0^2/(NDim J2)
            err1 = std::max(err1, std::fabs(f1 - ref)/ref);
        }
        // 3D
        double s3[6]; for (double &c : s3) c = s*urand();
        J2 = Plasticity::J2_3D(s3[0], s3[1], s3[2], s3[3], s3[4], s3[5]);
        const double ref3 = 0.5*(s3[0]*s3[0] + s3[1]*s3[1] + s3[2]*s3[2])
                          + (s3[3]*s3[3] + s3[4]*s3[4] + s3[5]*s3[5]);
        err0_3 = std::max(err0_3, std::fabs(J2 - ref3)/ref3);
        f0 = Plasticity::yieldFactor(J2, Y, 0, 3);
        if (f0 < 1.0){
            const double J2a = Plasticity::J2_3D(f0*s3[0], f0*s3[1], f0*s3[2], f0*s3[3], f0*s3[4], f0*s3[5]);
            err0_3 = std::max(err0_3, std::fabs(std::sqrt(3.*J2a)/Y - 1.));
        }
    }
    check(err0_2 < 1e-13, "YIELD_RETURN 0, 2D: sqrt(2 J2)/Y = 1 after the return", err0_2);
    check(dirErr < 1e-13, "YIELD_RETURN 0: direction of S unchanged", dirErr);
    check(err0_3 < 1e-13, "YIELD_RETURN 0, 3D: J2 form and sqrt(3 J2)/Y = 1 after return", err0_3);
    check(err1 < 1e-15,   "YIELD_RETURN 1: f = Y^2/(2 J2) exactly (GIZMO formula)", err1);
}

// 4. 2D deviatoric rate: model 1 traceless and = mu(dvx_dx - dvy_dy); model 0 = the
//    pre-2026-09-30 plane-strain expression bit for bit.
static void testDevRate(){
    std::srand(13);
    double tr1 = 0., e1 = 0.; bool same0 = true;
    for (int t = 0; t < 100000; ++t){
        const double a = urand(), b = urand(), c = urand(), d = urand(), mu = std::fabs(urand());
        double rxx, ryy, rxy;
        Plasticity::devRate2D(mu, a, b, c, d, rxx, ryy, rxy, 1);
        tr1 = std::max(tr1, std::fabs(rxx + ryy));
        e1 = std::max(e1, std::fabs(rxx - mu*(a - d)) + std::fabs(rxy - mu*(c + b)));
        Plasticity::devRate2D(mu, a, b, c, d, rxx, ryy, rxy, 0);
        const double oxx =  (4./3.)*mu*a - (2./3.)*mu*d;   // old Particles.cpp lines
        const double oyy = -(2./3.)*mu*a + (4./3.)*mu*d;
        const double oxy = mu * (c + b);
        if (rxx != oxx || ryy != oyy || rxy != oxy) same0 = false;
    }
    check(tr1 == 0.0, "SOLID_2D_MODEL 1: rate traceless in 2D", tr1);
    check(e1 < 1e-15, "SOLID_2D_MODEL 1: rxx = mu(dvx_dx - dvy_dy), rxy = mu(dvy_dx + dvx_dy)", e1);
    check(same0, "SOLID_2D_MODEL 0: old plane-strain rate reproduced bitwise");
}

// Material point: S' = devRate(L) + Jaumann, RK2 midpoint as in integrateStressTensor,
// then the return. L = velocity gradient (dvx_dx, dvx_dy, dvy_dx, dvy_dy), constant.
struct Point { double sxx = 0., syy = 0., sxy = 0., gp = 0.; };
static void step(Point &p, const double L[4], double dt, int model, int mode, double Yl = Y){
    double vxx, vyy, vxy;
    Plasticity::devRate2D(MU, L[0], L[1], L[2], L[3], vxx, vyy, vxy, model);
    const double R = L[2] - L[1];                          // d_x v_y - d_y v_x
    auto rate = [&](double sxx, double syy, double sxy, double &kxx, double &kyy, double &kxy){
        kxx = vxx - sxy*R; kyy = vyy + sxy*R; kxy = vxy + 0.5*(sxx - syy)*R; };
    double k1xx, k1yy, k1xy, k2xx, k2yy, k2xy;
    rate(p.sxx, p.syy, p.sxy, k1xx, k1yy, k1xy);
    rate(p.sxx + 0.5*dt*k1xx, p.syy + 0.5*dt*k1yy, p.sxy + 0.5*dt*k1xy, k2xx, k2yy, k2xy);
    p.sxx += dt*k2xx; p.syy += dt*k2yy; p.sxy += dt*k2xy;
    const double J2 = Plasticity::J2_2D(p.sxx, p.syy, p.sxy);
    const double f = Plasticity::yieldFactor(J2, Yl, mode, 2);
    if (f < 1.0){
        p.sxx *= f; p.syy *= f; p.sxy *= f;
        p.gp += Plasticity::plasticShearIncrement(J2, f, MU);
    }
}
static void drive(Point &p, const double L[4], double T, int n, int model, int mode, double Yl = Y){
    for (int k = 0; k < n; ++k) step(p, L, T/n, model, mode, Yl);
}

// 5. Material-point paths with the true-2D model.
static void testMaterialPoint(){
    const double tauY = Y/std::sqrt(2.), gY = tauY/MU;    // 0.0070711, 0.0153206
    // (a) pure shear to 3 gamma_Y: plateau tau_Y, gamma_p = gamma - gamma_Y
    {
        const double g = 1.0, L[4] = {0., g/2., g/2., 0.};
        Point p; drive(p, L, 3.*gY, 3000, 1, 0);
        check(std::fabs(p.sxy - tauY) < 1e-15 && p.sxx == 0. && p.syy == 0.,
              "pure shear, mode 0: plateau S_xy = Y/sqrt(2), no normal stress", std::fabs(p.sxy - tauY));
        check(std::fabs(p.gp - 2.*gY) < 1e-12, "pure shear, mode 0: gamma_p = gamma - gamma_Y",
              std::fabs(p.gp - 2.*gY));
        Point q; drive(q, L, 3.*gY, 3000, 0, 0);
        check(std::fabs(q.sxy - tauY) < 1e-15, "pure shear, model 0 (plane strain): same plateau",
              std::fabs(q.sxy - tauY));
    }
    // (b) uniaxial compression: S_xx -> -Y/sqrt(2), S_yy -> +Y/sqrt(2) (model 1)
    {
        const double L[4] = {-1., 0., 0., 0.};
        Point p; drive(p, L, 3.*gY, 3000, 1, 0);
        check(std::fabs(p.sxx + tauY) < 1e-15 && std::fabs(p.syy - tauY) < 1e-15,
              "uniaxial strain, model 1: S_xx = -Y/sqrt(2) = -S_yy", std::fabs(p.sxx + tauY));
        // model 0 (plane-strain rate, 2x2 J2): S_yy = -S_xx/2, J2 = 5/8 S_xx^2 -> |S_xx| = 2Y/sqrt(5)
        Point q; drive(q, L, 3.*gY, 3000, 0, 0);
        const double hyb = 2.*Y/std::sqrt(5.);
        check(std::fabs(-q.sxx - hyb) < 1e-12, "uniaxial strain, model 0: old hybrid |S_xx| = 2Y/sqrt(5)",
              std::fabs(-q.sxx - hyb));
    }
    // (c) load - unload - reverse (pure shear): unloading slope mu, reverse yield at -tau_Y
    {
        const double L[4] = {0., 0.5, 0.5, 0.}, Lr[4] = {0., -0.5, -0.5, 0.};
        Point p; drive(p, L, 3.*gY, 3000, 1, 0);
        const double gp0 = p.gp;
        drive(p, Lr, 1.9*gY, 1900, 1, 0);                  // elastic unloading, 1.9 gamma_Y
        check(std::fabs(p.sxy - (tauY - MU*1.9*gY)) < 1e-15 && p.gp == gp0,
              "unloading is elastic with slope mu (gamma_p frozen)", std::fabs(p.sxy - (tauY - MU*1.9*gY)));
        drive(p, Lr, 2.*gY, 2000, 1, 0);                   // through reverse yield
        check(std::fabs(p.sxy + tauY) < 1e-15, "reverse yield at S_xy = -Y/sqrt(2)", std::fabs(p.sxy + tauY));
        check(std::fabs(p.gp - gp0 - 1.9*gY) < 1e-12, "reverse plastic strain = 3.9 - 2 gamma_Y",
              std::fabs(p.gp - gp0 - 1.9*gY));
    }
    // (d) GIZMO return (mode 1): plateau error vs step size is first order
    {
        const double L[4] = {0., 0.5, 0.5, 0.};
        double err[3]; const int ns[3] = {60, 240, 960};
        for (int k = 0; k < 3; ++k){ Point p; drive(p, L, 3.*gY, ns[k], 1, 1); err[k] = std::fabs(p.sxy - tauY); }
        const double pOrd = std::log(err[1]/err[2]) / std::log(4.);
        std::printf("  mode 1 plateau error: %.3e %.3e %.3e (steps 60/240/960)\n", err[0], err[1], err[2]);
        check(pOrd > 0.8 && pOrd < 1.2, "YIELD_RETURN 1: plateau error converges at order 1", pOrd);
        check(err[2] < 0.01*tauY, "YIELD_RETURN 1: plateau within 1% of tau_Y at 960 steps", err[2]/tauY);
    }
    // (e) simple shear v_x = y, elastic (Y huge): Jaumann (Dienes) solution, model-independent,
    //     S_xy = mu sin(g), S_xx = -S_yy = mu (1 - cos g); the RK2 midpoint converges at order 2.
    {
        const double L[4] = {0., 1., 0., 0.}, g = 1.0;
        double err[2]; const int ns[2] = {500, 2000};
        for (int k = 0; k < 2; ++k){
            Point p; drive(p, L, g, ns[k], 1, 0, 1e30);
            err[k] = std::fabs(p.sxy - MU*std::sin(g)) + std::fabs(p.sxx - MU*(1. - std::cos(g)))
                   + std::fabs(p.syy + MU*(1. - std::cos(g)));
        }
        const double pOrd = std::log(err[0]/err[1]) / std::log(4.);
        check(err[1] < 1e-7*MU, "simple shear to g = 1 rad: Jaumann solution (2000 steps)", err[1]/MU);
        check(pOrd > 1.8 && pOrd < 2.2, "simple shear: RK2 midpoint order 2", pOrd);
    }
}

int main(){
    std::printf("Plasticity unit test (P0a), parameter.h: DIM=%d SOLID_2D_MODEL=%d YIELD_RETURN=%d\n",
                DIM, SOLID_2D_MODEL, YIELD_RETURN);
    testIdentityBelowYield();
    testReturnMaps();
    testDevRate();
    testMaterialPoint();
    std::printf("%s: %d failed check(s)\n", failures ? "FAILED" : "PASSED", failures);
    return failures;
}
