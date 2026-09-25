#!/usr/bin/env python3
"""2D linear travelling ELASTIC wave -- Hopkins' soundwave (2015, sec. 4.2.1)
with the shear modulus switched on.

Same geometry as testcases/soundwave_2d (periodic [0,1) x [0,1/4), square
lattice of n x n/4 equal-mass particles, wave in +x, uniform in y, A = 1e-6),
Tillotson EOS in its Hopkins-canonical form (rho0 = A = 1, everything else 0,
so P = A (rho/rho0 - 1) and K = A = 1 at rho0), but with mu > 0.

What changes with mu > 0
------------------------
With v_y = 0 and no y-dependence, integrateStressTensor (Particles.cpp:1150ff)
reduces to

    dS_xx/dt =  (4/3) mu dv_x/dx,   dS_yy/dt = -(2/3) mu dv_x/dx,   dS_xy/dt = 0

(the Jaumann terms vanish, R_xy = 0). Along a particle path dv_x/dx = -Dln rho/Dt,
so the deviatoric stress is slaved to the density,

    S_xx = -(4/3) mu ln(rho/rho_ref),   S_yy = +(2/3) mu ln(rho/rho_ref),

with rho_ref the density at which that particle's stress was zero. The normal
traction is t = P - S_xx = P + (4/3) mu ln(rho/rho0), so the linearised wave
speed is the longitudinal (P-wave) speed

    c_l = sqrt((K + 4/3 mu) / rho0)

which is exactly what xSplitElasticHLLC uses for its wave fan (HLLC.cpp:640).

This is why the IC must carry a stress. If S = 0 at t = 0 the deviatoric part
references each particle's OWN initial density while the pressure references
rho0, the state is not an eigenmode, and the wave splits into left- and
right-going parts -- at mu = K about 57% of the traction perturbation is
mis-seeded. The right-travelling eigenmode is

    rho(x) = rho_bar (1 + A sin(k x))
    v_x(x) = c_l A sin(k x)
    S_xx   = -(4/3) mu ln(rho/rho0),  S_yy = +(2/3) mu ln(rho/rho0),  S_xy = 0

and it is written to /Sxx, /Sxy, /Syy, which InitialDistribution reads when
present (same names as the snapshot output, so a snapshot can be fed back in).

mu defaults to Poisson's ratio nu = 0.3, i.e. mu = 3K(1-2nu)/(2(1+nu)) = 6/13,
giving c_l = sqrt(21/13) = 1.27098. The wave period is lambda/c_l, NOT 1.

Kernel-sum density offset: as in testcases/soundwave_2d, rho_i = m_i sum_j W_ij
on this lattice sits a constant above m/dx^2 that depends only on N_NN
(+4.238e-3 at N_NN = 32). --correct-nnn (ON by default here) divides the masses
by 1 + offset so the realised mean density is rho_bar = 1; the offset is
computed from the solver's own N_NN definition, not hard-coded. With corrected
masses P_bar = 0, so the wave has genuine tensile half-cycles.

u does not enter this Tillotson (a = b = B = 0 -> P = A mu_till), but it must
stay positive (u <= 0 logs an ERROR every step), so it gets the isentropic
profile around u_bar = 0.01 as in the mu = 0 case.

Usage:
    ./generateIC.py -n 64
    ./generateIC.py --all                  # n = 64, 128, 256
    ./generateIC.py --all --mu 1.0         # nu = 0.125, c_l = sqrt(7/3)
"""

import argparse
import numpy as np
import h5py

RESOLUTIONS = (64, 128, 256)

RHO_BAR = 1.0
LAMBDA = 1.0
BOX_X = 1.0
BOX_Y = 0.25

# Tillotson, Hopkins-canonical (as testcases/soundwave_2d): P = A (rho/rho0 - 1),
# K = rho dP/drho = A at rho = rho0.
TILL_RHO0, TILL_A = 1.0, 1.0
U_BAR = 0.01
# Poisson's ratio 0.3 at K = 1: mu = 3 K (1 - 2 nu) / (2 (1 + nu))
MU_DEFAULT = 3.0 * TILL_A * (1.0 - 2.0 * 0.3) / (2.0 * (1.0 + 0.3))


def eos_state(rho):
    """Tillotson pressure and isentropic u(rho) (cold branch, a = b = B = 0)."""
    A, r0 = TILL_A, TILL_RHO0

    def e(r):  # int P / rho^2 drho
        return A / r0 * np.log(r) + A / r
    return A * (rho / r0 - 1.0), U_BAR + e(rho) - e(RHO_BAR)


def longitudinal_speed(mu):
    """c_l = sqrt((K + 4/3 mu) / rho0), K = TILL_A at rho0 (HLLC.cpp:640)."""
    return np.sqrt((TILL_A + 4.0 / 3.0 * mu) / TILL_RHO0)


def lattice_kernel_sum(h_over_dx, R=12):
    """dx^2 * sum_j W(r_j, h) over the infinite square lattice (Wendland C2, 2D)."""
    i = np.arange(-R, R + 1)
    X, Y = np.meshgrid(i, i)
    u = np.hypot(X, Y) / h_over_dx
    t = np.clip(1.0 - u, 0.0, None)
    return (7.0 / np.pi / h_over_dx**2 * t**4 * (1.0 + 4.0 * u)).sum()


def kernel_offset(nnn):
    """Kernel-sum density offset for target neighbour number nnn (bisection on h/dx)."""
    lo, hi = 1.0, 10.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if np.pi * mid**2 * lattice_kernel_sum(mid) < nnn:
            lo = mid
        else:
            hi = mid
    return lattice_kernel_sum(0.5 * (lo + hi)) - 1.0


def lagrangian_to_euler(q, A, k):
    """Solve x + (A/k)(1 - cos kx) = q for x."""
    x = q.copy()
    for _ in range(50):
        f = x + A / k * (1.0 - np.cos(k * x)) - q
        dx = f / (1.0 + A * np.sin(k * x))
        x -= dx
        if np.max(np.abs(dx)) < 1e-17:
            break
    return x


def make_ic(n, A, mu, correct_nnn=32, out=None):
    k = 2.0 * np.pi / LAMBDA
    dx = BOX_X / n
    ny = int(round(BOX_Y / dx))
    if abs(ny * dx - BOX_Y) > 1e-12:
        raise ValueError("box height {} is not a multiple of dx = {}".format(BOX_Y, dx))
    q = (np.arange(n) + 0.5) * dx
    xs = lagrangian_to_euler(q, A, k)
    ys = (np.arange(ny) + 0.5) * dx
    X, Y = np.meshgrid(xs, ys, indexing="ij")
    pos = np.column_stack([X.ravel(), Y.ravel()])
    N = len(pos)

    rho = RHO_BAR * (1.0 + A * np.sin(k * pos[:, 0]))
    P, u = eos_state(rho)
    c_l = longitudinal_speed(mu)

    vel = np.zeros_like(pos)
    vel[:, 0] = c_l * A * np.sin(k * pos[:, 0])

    # hypo-elastic stress of the right-travelling eigenmode
    lnr = np.log(rho / TILL_RHO0)
    Sxx = -(4.0 / 3.0) * mu * lnr
    Syy = +(2.0 / 3.0) * mu * lnr
    Sxy = np.zeros(N)

    offset = 0.0 if correct_nnn is None else kernel_offset(correct_nnn)
    m = np.full(N, RHO_BAR * BOX_X * BOX_Y / N / (1.0 + offset))
    matId = np.zeros(N, dtype=np.int8)

    fname = out or "soundwave_elastic_N{}{}.h5".format(
        n, "" if correct_nnn is None else "_corr")
    with h5py.File(fname, "w") as f:
        f.create_dataset("x", data=pos)
        f.create_dataset("v", data=vel)
        f.create_dataset("m", data=m)
        f.create_dataset("u", data=u)
        f.create_dataset("materialId", data=matId)
        f.create_dataset("Sxx", data=Sxx)
        f.create_dataset("Sxy", data=Sxy)
        f.create_dataset("Syy", data=Syy)
        # metadata for the analysis (not read by the demonstrator)
        for key, val in (("eos", "tillotson"), ("n", n), ("ny", ny), ("amplitude", A),
                         ("wavelength", LAMBDA), ("box_x", BOX_X), ("box_y", BOX_Y),
                         ("rho_bar", RHO_BAR), ("mu", mu), ("K0", TILL_A),
                         ("c_l", c_l), ("period", LAMBDA / c_l),
                         ("P_bar", 0.0), ("kernel_offset", offset),
                         ("correct_nnn", 0 if correct_nnn is None else correct_nnn)):
            f.attrs[key] = val

    print("Wrote {}: {} x {} = {} particles, dx = {:g}, mu = {:.9g}, "
          "c_l = {:.9f}, period = {:.9f}, S_xx in [{:+.4e}, {:+.4e}], "
          "mass factor 1/(1 {:+.8e})".format(
              fname, n, ny, N, dx, mu, c_l, LAMBDA / c_l,
              Sxx.min(), Sxx.max(), offset))


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-n", "--resolution", type=int, help="particles per wavelength")
    ap.add_argument("-A", "--amplitude", type=float, default=1e-6,
                    help="delta_rho/rho (default: 1e-6, Hopkins 2015)")
    ap.add_argument("--mu", type=float, default=MU_DEFAULT,
                    help="shear modulus (default: %(default).9g, i.e. Poisson 0.3 at K = 1)")
    ap.add_argument("--correct-nnn", type=int, default=32, metavar="N_NN",
                    help="divide masses by 1 + kernel-sum offset for this N_NN "
                         "(default: %(default)s; use --no-correct for nominal masses)")
    ap.add_argument("--no-correct", action="store_true",
                    help="nominal masses, <rho> = 1 + offset (P_bar > 0, no tensile half-cycles)")
    ap.add_argument("-o", "--out", help="output file name (single resolution only)")
    ap.add_argument("--all", action="store_true",
                    help="generate all resolutions {}".format(RESOLUTIONS))
    args = ap.parse_args()

    corr = None if args.no_correct else args.correct_nnn
    if args.all:
        for n in RESOLUTIONS:
            make_ic(n, args.amplitude, args.mu, corr)
    else:
        if args.resolution is None:
            ap.error("need -n, or --all")
        make_ic(args.resolution, args.amplitude, args.mu, corr, args.out)
