// Standalone unit test: verify that Kernel::cubicSplineDh agrees with a
// central finite difference of Kernel::cubicSpline.
//
// Build (from demonstrator/):
//     make test-kernel
// or directly:
//     g++ -std=c++11 -I include -DDIM=2 tests/test_kernel.cpp \
//         src/Particles.cpp -o bin/test_kernel_2d
//
// Note: linking the full Particles.cpp pulls in many symbols. To keep this
// test self-contained we only need the two Kernel functions, so we re-declare
// and re-define a tiny copy here. The reference implementation lives in
// src/Particles.cpp; if it changes, update this file accordingly.

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>

#ifndef DIM
#define DIM 2
#endif

namespace Kernel {
    double cubicSpline(const double &r, const double &h);
    double cubicSplineDh(const double &r, const double &h);
}

// --- copied verbatim from src/Particles.cpp -------------------------------
double Kernel::cubicSpline(const double &r, const double &h) {
    double h2 = h/2.;
#if DIM == 2
    const double sigma = 10./(7.*M_PI*h2*h2);
#else
    const double sigma = 1./(M_PI*h2*h2*h2);
#endif
    const double q = r/h2;
    if (0. <= q && q <= 1.){
        return sigma*(1.-3./2.*q*q*(1.-q/2.));
    } else if (1. < q && q < 2.){
        return sigma/4.*pow(2.-q, 3.);
    } else {
        return 0.;
    }
}

double Kernel::cubicSplineDh(const double &r, const double &h){
    const double h2 = h/2.;
#if DIM == 2
    const double sigma = 10./(7.*M_PI*h2*h2);
#else
    const double sigma = 1./(M_PI*h2*h2*h2);
#endif
    const double q = r/h2;
    double f, fPrime;
    if (0. <= q && q <= 1.){
        f = 1. - 3./2.*q*q*(1. - q/2.);
        fPrime = -3.*q + 9./4.*q*q;
    } else if (1. < q && q < 2.){
        f = pow(2.-q, 3.)/4.;
        fPrime = -3./4.*pow(2.-q, 2.);
    } else {
        return 0.;
    }
    return -(sigma/h) * ((double)DIM * f + q * fPrime);
}
// --------------------------------------------------------------------------

int main(){
    const double h = 0.7;
    const double eps = 1e-6;
    // Sample radii across the support [0, h]. Skip the kink points r = h/2
    // and r = h where the analytic derivative is continuous but the
    // finite-difference is sensitive to neighbouring evaluation points.
    const double sample[] = {
        0.001*h, 0.05*h, 0.10*h, 0.20*h, 0.35*h,
        0.55*h, 0.65*h, 0.75*h, 0.85*h, 0.95*h
    };
    const int nSample = sizeof(sample)/sizeof(double);

    int failures = 0;
    double maxRelErr = 0.;
    for (int i = 0; i < nSample; ++i){
        const double r = sample[i];
        const double analytic = Kernel::cubicSplineDh(r, h);
        const double fd = (Kernel::cubicSpline(r, h+eps)
                          - Kernel::cubicSpline(r, h-eps)) / (2.*eps);
        const double absErr = std::fabs(analytic - fd);
        const double scale = std::max(1e-12, std::fabs(fd));
        const double relErr = absErr / scale;
        if (relErr > maxRelErr) maxRelErr = relErr;
        const bool ok = (absErr < 1e-7) || (relErr < 1e-5);
        std::printf("  r/h=%.3f  analytic=% .6e  FD=% .6e  relErr=%.2e  %s\n",
                    r/h, analytic, fd, relErr, ok ? "OK" : "FAIL");
        if (!ok) ++failures;
    }
    // Outside the support both must be exactly zero.
    if (Kernel::cubicSplineDh(1.5*h, h) != 0.0){
        std::printf("  out-of-support nonzero!\n");
        ++failures;
    }

    std::printf("DIM=%d  max relative error = %.2e\n", DIM, maxRelErr);
    if (failures){
        std::printf("FAILED (%d)\n", failures);
        return 1;
    }
    std::printf("PASSED\n");
    return 0;
}
