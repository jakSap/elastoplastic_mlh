#ifndef MESHLESSHYDRO_KERNEL_H
#define MESHLESSHYDRO_KERNEL_H

#include <cmath>

// Allow standalone use (unit tests) with externally defined DIM/KERNEL_FUNCTION.
#ifndef DIM
#include "parameter.h"
#endif
#ifndef KERNEL_FUNCTION
#define KERNEL_FUNCTION 3
#endif

// All kernels use the GIZMO convention: compact support radius = h, u = r/h.
namespace Kernel {

    inline double cubicSpline(const double &r, const double &h) {

        // TODO: remove this
        double h2 = h/2.;
#if DIM == 2
        const double sigma = 10./(7.*M_PI*h2*h2);
#else // DIM == 3
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

    // dW/dh = -(sigma/h) * (DIM * f(q) + q * f'(q)) with q = 2r/h.
    inline double cubicSplineDh(const double &r, const double &h){
        const double h2 = h/2.;
#if DIM == 2
        const double sigma = 10./(7.*M_PI*h2*h2);
#else // DIM == 3
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

    // dW(r, h)/dr. Scalar. For needs to be multiplied with (x_i - x_j)/r
    inline double dWdr(const double &r, const double &h){
        const double sigma = 40./(7.*M_PI);
        const double q = r/h;
        if (0. <= q && q < 1./2.){
            return 6 * sigma / pow(h, DIM + 1) * (3*pow(q, 2) - 2 * q);
        } else if (1./2. <= q && q <= 1.){
            return 6 * sigma / pow(h, DIM + 1) * -1 * pow((1 - q), 2);
        } else {
            return 0.;
        }
    }

    // dW/dh (legacy SPH variant, 2D cubic spline)
    inline double dWdh(const double &r, const double &h){
        const double sigma = 40./(7.*M_PI*h*h*h);
        const double q = r/h;
        if (0. <= q && q <= 1./2.){
            return 2 * sigma / pow(h, 6) * (12 * h * r*r - pow(h, 3) - 15 * pow(r, 3));
        } else if (1./2. < q && q < 1.){
            return 2 * sigma / pow(h, 6) * pow((h - r), 2) * (5 * r - 2 * h);
        } else {
            return 0.;
        }
    }

// Wendland kernels as in GIZMO kernel.h (2D/3D forms and norms).
#if DIM == 2
#define WC2_NORM (7./M_PI)
#define WC4_NORM (9./M_PI)
#define WC6_NORM (78./(7.*M_PI))
#else // DIM == 3
#define WC2_NORM (21./(2.*M_PI))
#define WC4_NORM (495./(32.*M_PI))
#define WC6_NORM (1365./(64.*M_PI))
#endif

    inline double wendlandC2(const double &r, const double &h){
        const double u = r/h;
        if (u >= 1.) return 0.;
        const double t = 1.-u, t2 = t*t;
        return WC2_NORM/pow(h, DIM) * t2*t2*(1.+4.*u);
    }
    inline double wendlandC2Dr(const double &r, const double &h){
        const double u = r/h;
        if (u >= 1.) return 0.;
        const double t = 1.-u;
        return WC2_NORM/pow(h, DIM+1) * (-20.*u*t*t*t);
    }
    inline double wendlandC2Dh(const double &r, const double &h){
        const double u = r/h;
        if (u >= 1.) return 0.;
        const double t = 1.-u, t3 = t*t*t;
        return -WC2_NORM/pow(h, DIM+1) * ((double)DIM*t3*t*(1.+4.*u) - 20.*u*u*t3);
    }

    inline double wendlandC4(const double &r, const double &h){
        const double u = r/h;
        if (u >= 1.) return 0.;
        const double t = 1.-u, t3 = t*t*t;
        return WC4_NORM/pow(h, DIM) * t3*t3*(1.+6.*u+(35./3.)*u*u);
    }
    inline double wendlandC4Dr(const double &r, const double &h){
        const double u = r/h;
        if (u >= 1.) return 0.;
        const double t = 1.-u, t2 = t*t;
        return WC4_NORM/pow(h, DIM+1) * (-(56./3.)*u*t2*t2*t*(1.+5.*u));
    }
    inline double wendlandC4Dh(const double &r, const double &h){
        const double u = r/h;
        if (u >= 1.) return 0.;
        const double t = 1.-u, t2 = t*t, t5 = t2*t2*t;
        return -WC4_NORM/pow(h, DIM+1)
               * ((double)DIM*t5*t*(1.+6.*u+(35./3.)*u*u) - (56./3.)*u*u*t5*(1.+5.*u));
    }

    inline double wendlandC6(const double &r, const double &h){
        const double u = r/h;
        if (u >= 1.) return 0.;
        const double t = 1.-u, t2 = t*t, t8 = t2*t2*t2*t2;
        return WC6_NORM/pow(h, DIM) * t8*(1.+8.*u+25.*u*u+32.*u*u*u);
    }
    inline double wendlandC6Dr(const double &r, const double &h){
        const double u = r/h;
        if (u >= 1.) return 0.;
        const double t = 1.-u, t2 = t*t, t7 = t2*t2*t2*t;
        return WC6_NORM/pow(h, DIM+1) * (-22.*u*t7*(1.+7.*u+16.*u*u));
    }
    inline double wendlandC6Dh(const double &r, const double &h){
        const double u = r/h;
        if (u >= 1.) return 0.;
        const double t = 1.-u, t2 = t*t, t7 = t2*t2*t2*t;
        return -WC6_NORM/pow(h, DIM+1)
               * ((double)DIM*t7*t*(1.+8.*u+25.*u*u+32.*u*u*u) - 22.*u*u*t7*(1.+7.*u+16.*u*u));
    }

    // Generic kernel selected by KERNEL_FUNCTION (GIZMO numbering:
    // 3 = cubic spline, 6 = Wendland C2, 7 = Wendland C4, 9 = Wendland C6).
#if KERNEL_FUNCTION == 3
    inline double W(const double &r, const double &h){ return cubicSpline(r, h); }
    inline double WDh(const double &r, const double &h){ return cubicSplineDh(r, h); }
    inline double WDr(const double &r, const double &h){ return dWdr(r, h); }
#elif KERNEL_FUNCTION == 6
    inline double W(const double &r, const double &h){ return wendlandC2(r, h); }
    inline double WDh(const double &r, const double &h){ return wendlandC2Dh(r, h); }
    inline double WDr(const double &r, const double &h){ return wendlandC2Dr(r, h); }
#elif KERNEL_FUNCTION == 7
    inline double W(const double &r, const double &h){ return wendlandC4(r, h); }
    inline double WDh(const double &r, const double &h){ return wendlandC4Dh(r, h); }
    inline double WDr(const double &r, const double &h){ return wendlandC4Dr(r, h); }
#elif KERNEL_FUNCTION == 9
    inline double W(const double &r, const double &h){ return wendlandC6(r, h); }
    inline double WDh(const double &r, const double &h){ return wendlandC6Dh(r, h); }
    inline double WDr(const double &r, const double &h){ return wendlandC6Dr(r, h); }
#else
#error "Unsupported KERNEL_FUNCTION. Supported: 3 (cubic spline), 6/7/9 (Wendland C2/C4/C6)."
#endif

} // namespace Kernel

#endif // MESHLESSHYDRO_KERNEL_H
