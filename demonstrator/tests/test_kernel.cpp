// Standalone unit test for include/Kernel.h: verifies the analytic dW/dh and
// dW/dr of every kernel family against central finite differences, and checks
// the 2D/3D normalization integral of W.
//
// Build (from demonstrator/):
//     make test-kernel
// or directly:
//     g++ -std=c++11 -O2 -I include -DDIM=2 tests/test_kernel.cpp -o bin/test_kernel_2d

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>

#ifndef DIM
#define DIM 2
#endif

#include "Kernel.h"

typedef double (*KernelFn)(const double&, const double&);

static int failures = 0;

// FD check of dF agreeing with d/dx F along the given argument (0: r, 1: h)
static void checkDerivative(const char *name, KernelFn F, KernelFn dF, int wrt){
    const double h = 0.7;
    const double eps = 1e-6;
    // skip kink points (r = h/2 for the cubic spline) and the support edge
    const double sample[] = {
        0.001, 0.05, 0.10, 0.20, 0.35,
        0.55, 0.65, 0.75, 0.85, 0.95
    };
    const int nSample = sizeof(sample)/sizeof(double);
    double maxRelErr = 0.;
    for (int i = 0; i < nSample; ++i){
        const double r = sample[i]*h;
        const double analytic = dF(r, h);
        const double fd = (wrt == 0)
            ? (F(r+eps, h) - F(r-eps, h)) / (2.*eps)
            : (F(r, h+eps) - F(r, h-eps)) / (2.*eps);
        const double absErr = std::fabs(analytic - fd);
        const double relErr = absErr / std::max(1e-12, std::fabs(fd));
        if (relErr > maxRelErr) maxRelErr = relErr;
        if (!((absErr < 1e-6) || (relErr < 1e-5))){
            std::printf("  %s r/h=%.3f analytic=% .6e FD=% .6e relErr=%.2e FAIL\n",
                        name, sample[i], analytic, fd, relErr);
            ++failures;
        }
    }
    if (dF(1.5*h, h) != 0.0){
        std::printf("  %s out-of-support nonzero!\n", name);
        ++failures;
    }
    std::printf("  %-16s max relative error = %.2e\n", name, maxRelErr);
}

// radial quadrature of W: must integrate to 1 over the kernel support
static void checkNorm(const char *name, KernelFn F){
    const double h = 0.7;
    const int n = 200000;
    double sum = 0.;
    for (int i = 0; i < n; ++i){
        const double r = (i + .5) * h / n;
#if DIM == 2
        sum += F(r, h) * 2.*M_PI*r * (h/n);
#else
        sum += F(r, h) * 4.*M_PI*r*r * (h/n);
#endif
    }
    const bool ok = std::fabs(sum - 1.) < 1e-4;
    std::printf("  %-16s integral = %.6f  %s\n", name, sum, ok ? "OK" : "FAIL");
    if (!ok) ++failures;
}

int main(){
    std::printf("DIM=%d, KERNEL_FUNCTION=%d\n", DIM, KERNEL_FUNCTION);

    checkNorm("cubicSpline", Kernel::cubicSpline);
    checkNorm("wendlandC2", Kernel::wendlandC2);
    checkNorm("wendlandC4", Kernel::wendlandC4);
    checkNorm("wendlandC6", Kernel::wendlandC6);

    checkDerivative("cubicSplineDh", Kernel::cubicSpline, Kernel::cubicSplineDh, 1);
    checkDerivative("wendlandC2Dh", Kernel::wendlandC2, Kernel::wendlandC2Dh, 1);
    checkDerivative("wendlandC4Dh", Kernel::wendlandC4, Kernel::wendlandC4Dh, 1);
    checkDerivative("wendlandC6Dh", Kernel::wendlandC6, Kernel::wendlandC6Dh, 1);

    checkDerivative("wendlandC2Dr", Kernel::wendlandC2, Kernel::wendlandC2Dr, 0);
    checkDerivative("wendlandC4Dr", Kernel::wendlandC4, Kernel::wendlandC4Dr, 0);
    checkDerivative("wendlandC6Dr", Kernel::wendlandC6, Kernel::wendlandC6Dr, 0);
#if DIM == 2
    // legacy dWdr carries the 2D norm only
    checkDerivative("dWdr(cubic)", Kernel::cubicSpline, Kernel::dWdr, 0);
#endif

    if (failures){
        std::printf("FAILED (%d)\n", failures);
        return 1;
    }
    std::printf("PASSED\n");
    return 0;
}
