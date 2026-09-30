#!/usr/bin/env python3
"""2D linear travelling soundwave (Hopkins 2015, sec. 4.2.1; Stone et al. 2008).

Hopkins: periodic unit domain, gamma = 5/3 gas with unit mean density and unit
sound speed (P = 3/5), travelling wave with delta_rho/rho = 1e-6 and unit
wavelength, evolved for one period (t = 1), after which the solution equals
the initial condition. L1 = 1/N sum_i |rho_i - rho(x_i)|.

2D version here: periodic rectangle [0, 1) x [0, 1/4), square lattice of
n x n/4 equal-mass particles (n = particles per wavelength), wave travelling
in +x (uniform in y):

    rho(x) = rho_bar * (1 + A sin(k x)),   k = 2 pi / lambda
    v_x(x) = c * A sin(k x),               v_y = 0
    u(x)   = isentropic: u_bar + int_{rho_bar}^{rho} P/rho'^2 drho'

The density profile is imposed by displacing the lattice (equal masses
m = rho_bar dx^2), so the kernel-sum density reproduces rho(x): each particle's
Lagrangian coordinate q is mapped to x by solving the cumulative-mass relation

    x + (A/k) (1 - cos(k x)) = q     (Newton, to machine precision).

All three EOS are normalised to rho_bar = 1 and c = 1 at the nominal mean
state, with the deviatoric stress switched off (mu = 0 in the config):

    ideal      gamma = 5/3, P_bar = 3/5, u_bar = 0.9    (Hopkins)
    murnaghan  K0 = 1, n = 2, rho0 = 1 -> P_bar = 0, c^2 = K0/rho0 = 1
    tillotson  rho0 = A = 1, a = b = B = alpha = beta = 0
               -> P = A (rho/rho0 - 1), P_bar = 0, c^2 = A/rho0 = 1
               u_bar = 0.01 > 0 (u <= 0 is reported as an ERROR every step
               for EOS != 1), below u_iv = 0.05 so the cold branch is used.

For Murnaghan u does not enter the EOS; it gets the same u_bar for symmetry.

Kernel-sum density offset: the demonstrator's density rho_i = m_i sum_j W_ij
(Wendland C2, self term included) on this lattice sits above m/dx^2 by a
constant that depends only on h/dx, i.e. on N_NN: +4.238e-3 for N_NN = 32,
+7.30e-4 for N_NN = 64. With nominal masses (default) the realised mean state
is rho = 1 + offset, which shifts P_bar (ideal 0.6025, solids +4.2e-3) and,
for Murnaghan (c ~ sqrt(rho)), the sound speed to 1.00212.

--correct-nnn N_NN divides the masses by (1 + offset(N_NN)), so the solver's
kernel-sum mean density is rho_bar. offset(N_NN) is computed, not hard-coded:
h/dx is solved from N_NN = pi h^2 sum_j W(r_j, h) (the solver's definition,
Particles.cpp:2622) on the infinite lattice, then offset = dx^2 sum_j W - 1.
This reproduces the solver's measured offset to 3e-11 (N_NN = 32) and 3e-10
(N_NN = 64, limited by the smoothing-length tolerance). Output file:
soundwave_<eos>_N<n>[_nnn<N_NN>]_corr.h5 (the _nnn part only if N_NN != 32).

Usage:
    ./generateIC.py -n 64 -e ideal
    ./generateIC.py -n 64 -e murnaghan --correct-nnn 32
    ./generateIC.py --all            # 3 EOS x {64, 128, 256, 512}, nominal masses
"""

import argparse
import numpy as np
import h5py

RESOLUTIONS = (64, 128, 256, 512)
EOS_NAMES = ("ideal", "murnaghan", "tillotson")

RHO_BAR = 1.0
LAMBDA = 1.0
BOX_X = 1.0
BOX_Y = 0.25

# ideal gas
GAMMA = 5.0 / 3.0
P_BAR_IDEAL = 3.0 / 5.0
# Murnaghan (repo standard: testcases/rip_test)
MURN_K0, MURN_N, MURN_RHO0 = 1.0, 2.0, 1.0
# Tillotson, Hopkins-canonical (testcases/periodic_wall_bounce)
TILL_RHO0, TILL_A = 1.0, 1.0
# mean specific internal energy for the solid EOS
U_BAR_SOLID = 1.0   # e0 >= 1 directive; dP/du = 0 for a=b=B=alpha=beta=0, so dynamics unchanged


def eos_state(eos, rho):
    """Pressure, sound speed and isentropic u(rho) (relative to rho_bar)."""
    if eos == "ideal":
        u_bar = P_BAR_IDEAL / ((GAMMA - 1.0) * RHO_BAR)
        u = u_bar * (rho / RHO_BAR) ** (GAMMA - 1.0)
        P = (GAMMA - 1.0) * rho * u
        c = np.sqrt(GAMMA * P / rho)
    elif eos == "murnaghan":
        K0, n, r0 = MURN_K0, MURN_N, MURN_RHO0

        def e(r):  # int P/rho^2 drho, as EquationOfState::EOSInternalEnergy
            return (K0 / (n * (n - 1.0) * r0**n) * (r**(n - 1.0) - r0**(n - 1.0))
                    + K0 / n * (1.0 / r - 1.0 / r0))
        u = U_BAR_SOLID + e(rho) - e(RHO_BAR)
        P = K0 / n * ((rho / r0) ** n - 1.0)
        c = np.sqrt(K0 / r0 * (rho / r0) ** (n - 1.0))
    elif eos == "tillotson":
        A, r0 = TILL_A, TILL_RHO0

        def e(r):  # int A (r/r0 - 1) / r^2 dr
            return A / r0 * np.log(r) + A / r
        u = U_BAR_SOLID + e(rho) - e(RHO_BAR)
        P = A * (rho / r0 - 1.0)
        c = np.sqrt(A / r0) * np.ones_like(rho)
    else:
        raise ValueError(eos)
    return P, c, u


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


def make_ic(n, eos, A, correct_nnn=None):
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
    P_bar, c_bar, _ = eos_state(eos, np.array([RHO_BAR]))
    P_bar, c_bar = float(P_bar[0]), float(c_bar[0])
    P, _, u = eos_state(eos, rho)

    vel = np.zeros_like(pos)
    vel[:, 0] = c_bar * A * np.sin(k * pos[:, 0])
    offset = 0.0 if correct_nnn is None else kernel_offset(correct_nnn)
    m = np.full(N, RHO_BAR * BOX_X * BOX_Y / N / (1.0 + offset))
    matId = np.zeros(N, dtype=np.int8)

    fname = "soundwave_{}_N{}".format(eos, n)
    if correct_nnn is not None:
        fname += ("" if correct_nnn == 32 else "_nnn{}".format(correct_nnn)) + "_corr"
    fname += ".h5"
    with h5py.File(fname, "w") as f:
        f.create_dataset("x", data=pos)
        f.create_dataset("v", data=vel)
        f.create_dataset("m", data=m)
        f.create_dataset("u", data=u)
        f.create_dataset("materialId", data=matId)
        # exact point-value density rho(x_i): seeds rhoExplicit under EXPLICIT_VOL_SEED_IC
        # (ignored otherwise). thirdorder study 2026-10-01.
        f.create_dataset("rho", data=rho)
        # metadata for the analysis (not read by the demonstrator)
        for key, val in (("eos", eos), ("n", n), ("ny", ny), ("amplitude", A),
                         ("wavelength", LAMBDA), ("box_x", BOX_X), ("box_y", BOX_Y),
                         ("rho_bar", RHO_BAR), ("c_bar", c_bar), ("P_bar", P_bar),
                         ("kernel_offset", offset),
                         ("correct_nnn", 0 if correct_nnn is None else correct_nnn)):
            f.attrs[key] = val

    print("Wrote {}: {} x {} = {} particles, dx = {:g}, c = {:.6f}, P_bar = {:.6g}, "
          "u in [{:.12g}, {:.12g}], max|x-q| = {:.3g}, mass factor 1/(1 {:+.8e})".format(
              fname, n, ny, N, dx, c_bar, P_bar, u.min(), u.max(), np.max(np.abs(xs - q)), offset))


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-n", "--resolution", type=int, help="particles per unit length")
    ap.add_argument("-e", "--eos", choices=EOS_NAMES, help="equation of state")
    ap.add_argument("-A", "--amplitude", type=float, default=1e-6,
                    help="delta_rho/rho (default: 1e-6, Hopkins 2015)")
    ap.add_argument("--correct-nnn", type=int, default=None, metavar="N_NN",
                    help="divide masses by 1 + kernel-sum offset for this N_NN (writes *_corr.h5)")
    ap.add_argument("--all", action="store_true",
                    help="generate all EOS x resolutions {}".format(RESOLUTIONS))
    args = ap.parse_args()

    if args.all:
        for eos in EOS_NAMES:
            for n in RESOLUTIONS:
                make_ic(n, eos, args.amplitude, args.correct_nnn)
    else:
        if args.resolution is None or args.eos is None:
            ap.error("need -n and -e, or --all")
        make_ic(args.resolution, args.eos, args.amplitude, args.correct_nnn)
