#!/usr/bin/env python3
"""2D linear travelling ELASTIC wave (mu > 0) for the third-order study -- the 2026-09-25 elastic eigenmode
(testcases/soundwave_2d/20260925_xsplit_eos/generateIC_elastic.py) for Murnaghan AND Tillotson, plus the exact
point-value density /rho that seeds rhoExplicit under EXPLICIT_VOL_SEED_IC.

Geometry as generateIC.py: periodic [0,1) x [0,1/4), n x n/4 equal-mass square lattice displaced to the wave,
wave in +x, uniform in y. Both EOS have K = 1 at rho0 = 1 (identical linear problem, different nonlinearity):

    murnaghan  K0 = 1, n = 2, rho0 = 1:     P = (rho^2 - 1)/2,  dP/drho = rho
    tillotson  rho0 = A = 1, rest 0:        P = rho - 1,        dP/drho = 1

With v_y = 0 and no y-dependence the hypo-elastic law reduces to dS_xx/dt = (4/3) mu dv_x/dx (plane strain,
SOLID_2D_MODEL 0) and dS_yy/dt = -(2/3) mu dv_x/dx, so the stress is slaved to the density,
S_xx = -(4/3) mu ln(rho/rho0), S_yy = +(2/3) mu ln(rho/rho0), and the normal stress acts like a barotropic
pressure Pi(rho) = P(rho) + (4/3) mu ln(rho/rho0) with c_l^2 = dPi/drho = P'(rho) + (4/3) mu/rho
(at rho0: (K + 4/3 mu)/rho0). The right-travelling linear eigenmode is

    rho = rho_bar (1 + A sin kx),  v_x = c_l A sin kx,  S as above,  u isentropic around u_bar = 1.

Masses are NOMINAL (m = rho_bar dx^2) by default: the density the scheme uses is rhoExplicit, seeded from /rho,
so the kernel-sum offset does not enter and V = m/rho_E tiles the box exactly. (--correct-nnn as in
generateIC.py would only matter without EXPLICIT_VOL_SEED_IC.) P_bar = 0: half of each cycle is tensile.

u does not enter either EOS here (a = b = B = alpha = beta = 0 for Tillotson; Murnaghan is barotropic); it gets
the isentropic profile around u_bar = 1 (e0 >= 1 directive) including the elastic work.

--mode shear: the transverse wave instead (exercises the tangential traction of xSplitElasticHLLC):
    v_y = c_s A sin kx,  S_xy = -mu A sin kx  (S_xy = mu du_y/dx, travelling +x with c_s = sqrt(mu/rho0)),
    rho = rho_bar, v_x = 0, lattice sheared by u_y = (A/k) cos kx (consistent with the stress). Linear to O(A^2):
    the Jaumann terms only create S_xx, S_yy ~ A^2, which feed back on v_y at O(A^3).

Usage:
    ./generateIC_elastic.py -n 64 -e murnaghan            # -> ics/elastic_murnaghan_N64.h5
    ./generateIC_elastic.py -n 64 -e tillotson --mode shear   # -> ics/shear_tillotson_N64.h5
    ./generateIC_elastic.py -n 64 -e tillotson -A 1e-7     # -> ics/elastic_tillotson_A1e-07_N64.h5
"""
import argparse, os
import numpy as np
import h5py
from generateIC import lagrangian_to_euler, kernel_offset

HERE = os.path.dirname(os.path.abspath(__file__))
RHO_BAR, LAMBDA, BOX_X, BOX_Y = 1.0, 1.0, 1.0, 0.25
MURN_K0, MURN_N, MURN_RHO0 = 1.0, 2.0, 1.0
TILL_RHO0, TILL_A = 1.0, 1.0
U_BAR = 1.0
MU_DEFAULT = 3.0 * 1.0 * (1.0 - 2.0 * 0.3) / (2.0 * (1.0 + 0.3))   # Poisson 0.3 at K = 1: 6/13


def pressure(eos, rho):
    if eos == "murnaghan":
        return MURN_K0 / MURN_N * ((rho / MURN_RHO0) ** MURN_N - 1.0)
    return TILL_A * (rho / TILL_RHO0 - 1.0)


def dPdrho(eos, rho):
    if eos == "murnaghan":
        return MURN_K0 / MURN_RHO0 * (rho / MURN_RHO0) ** (MURN_N - 1.0)
    return TILL_A / TILL_RHO0 * np.ones_like(rho)


def rho0_of(eos):
    return MURN_RHO0 if eos == "murnaghan" else TILL_RHO0


def wave_constants(eos, mu, rho_bar=RHO_BAR):
    """c0 = c_l(rho_bar) and the nonlinearity c1 = (rho/c) dc/drho of Pi = P + 4/3 mu ln rho (analysis.py)."""
    r = np.array([rho_bar])
    c2 = float(dPdrho(eos, r)[0]) + 4.0 / 3.0 * mu / rho_bar
    eps = 1e-6 * rho_bar
    d2 = (float(dPdrho(eos, r + eps)[0]) - float(dPdrho(eos, r - eps)[0])) / (2 * eps) - 4.0 / 3.0 * mu / rho_bar**2
    return np.sqrt(c2), rho_bar * d2 / (2.0 * c2)


def internal_energy(eos, rho, mu):
    """u_bar + int_{rho_bar}^{rho} Pi/rho'^2 drho' (isentropic, elastic work included), by quadrature."""
    out = np.empty_like(rho)
    g = np.linspace(0.0, 1.0, 65)
    for i, r in enumerate(rho):
        rr = RHO_BAR + (r - RHO_BAR) * g
        Pi = pressure(eos, rr) + 4.0 / 3.0 * mu * np.log(rr / rho0_of(eos))
        out[i] = U_BAR + np.trapz(Pi / rr**2, rr)
    return out


def make_ic(n, eos, A, mu, correct_nnn=None, mode="long"):
    k = 2.0 * np.pi / LAMBDA
    dx = BOX_X / n
    ny = int(round(BOX_Y / dx))
    if abs(ny * dx - BOX_Y) > 1e-12:
        raise ValueError("box height not a multiple of dx")
    q = (np.arange(n) + 0.5) * dx
    xs = lagrangian_to_euler(q, A, k)
    ys = (np.arange(ny) + 0.5) * dx
    X, Y = np.meshgrid(xs, ys, indexing="ij")
    if mode == "shear":
        X, Y = np.meshgrid(q, ys, indexing="ij")
        Y = Y + A / k * np.cos(k * X)
    pos = np.column_stack([X.ravel(), Y.ravel()])
    N = len(pos)
    if mode == "shear":
        return make_shear(n, ny, eos, A, mu, k, pos, correct_nnn)

    rho = RHO_BAR * (1.0 + A * np.sin(k * pos[:, 0]))
    c_l, c1 = wave_constants(eos, mu)
    vel = np.zeros_like(pos)
    vel[:, 0] = c_l * A * np.sin(k * pos[:, 0])
    lnr = np.log(rho / rho0_of(eos))
    Sxx, Syy, Sxy = -(4.0 / 3.0) * mu * lnr, (2.0 / 3.0) * mu * lnr, np.zeros(N)
    u = internal_energy(eos, rho, mu)
    offset = 0.0 if correct_nnn is None else kernel_offset(correct_nnn)
    m = np.full(N, RHO_BAR * BOX_X * BOX_Y / N / (1.0 + offset))

    name = "elastic_%s%s_N%d%s.h5" % (eos, "" if A == 1e-6 else "_A%.0e" % A, n,
                                      "" if correct_nnn is None else "_corr")
    os.makedirs(os.path.join(HERE, "ics"), exist_ok=True)
    fn = os.path.join(HERE, "ics", name)
    with h5py.File(fn, "w") as f:
        f.create_dataset("x", data=pos)
        f.create_dataset("v", data=vel)
        f.create_dataset("m", data=m)
        f.create_dataset("u", data=u)
        f.create_dataset("materialId", data=np.zeros(N, dtype=np.int8))
        f.create_dataset("Sxx", data=Sxx)
        f.create_dataset("Sxy", data=Sxy)
        f.create_dataset("Syy", data=Syy)
        f.create_dataset("rho", data=rho)
        for key, val in (("eos", eos), ("n", n), ("ny", ny), ("amplitude", A), ("wavelength", LAMBDA),
                         ("box_x", BOX_X), ("box_y", BOX_Y), ("rho_bar", RHO_BAR), ("mu", mu),
                         ("c_bar", c_l), ("c_l", c_l), ("c1", c1), ("period", LAMBDA / c_l),
                         ("P_bar", 0.0), ("kernel_offset", offset),
                         ("correct_nnn", 0 if correct_nnn is None else correct_nnn)):
            f.attrs[key] = val
    print("Wrote %s: %d particles, mu = %.9g, c_l = %.9f, c1 = %+.6f, period = %.9f, u in [%.6f, %.6f]"
          % (name, N, mu, c_l, c1, LAMBDA / c_l, u.min(), u.max()))
    return fn


def make_shear(n, ny, eos, A, mu, k, pos, correct_nnn):
    N = len(pos)
    c_s = np.sqrt(mu / RHO_BAR)
    c_l, _ = wave_constants(eos, mu)
    s = np.sin(k * pos[:, 0])
    vel = np.zeros_like(pos)
    vel[:, 1] = c_s * A * s
    offset = 0.0 if correct_nnn is None else kernel_offset(correct_nnn)
    m = np.full(N, RHO_BAR * BOX_X * BOX_Y / N / (1.0 + offset))
    name = "shear_%s%s_N%d%s.h5" % (eos, "" if A == 1e-6 else "_A%.0e" % A, n, "" if correct_nnn is None else "_corr")
    fn = os.path.join(HERE, "ics", name)
    with h5py.File(fn, "w") as f:
        f.create_dataset("x", data=pos)
        f.create_dataset("v", data=vel)
        f.create_dataset("m", data=m)
        f.create_dataset("u", data=U_BAR + 0.5 * mu * (A * s)**2 / RHO_BAR)   # + elastic energy density / rho
        f.create_dataset("materialId", data=np.zeros(N, dtype=np.int8))
        f.create_dataset("Sxx", data=np.zeros(N))
        f.create_dataset("Sxy", data=-mu * A * s)
        f.create_dataset("Syy", data=np.zeros(N))
        f.create_dataset("rho", data=np.full(N, RHO_BAR))
        for key, val in (("eos", eos), ("mode", "shear"), ("n", n), ("ny", ny), ("amplitude", A),
                         ("wavelength", LAMBDA), ("box_x", BOX_X), ("box_y", BOX_Y), ("rho_bar", RHO_BAR),
                         ("mu", mu), ("c_s", c_s), ("c_l", c_l), ("period", LAMBDA / c_s), ("P_bar", 0.0),
                         ("kernel_offset", offset), ("correct_nnn", 0 if correct_nnn is None else correct_nnn)):
            f.attrs[key] = val
    print("Wrote %s: %d particles, mu = %.9g, c_s = %.9f, c_l = %.9f, period = %.9f"
          % (name, N, mu, c_s, c_l, LAMBDA / c_s))
    return fn


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-n", type=int, nargs="+", required=True)
    ap.add_argument("-e", "--eos", choices=("murnaghan", "tillotson"), nargs="+", required=True)
    ap.add_argument("-A", "--amplitude", type=float, default=1e-6)
    ap.add_argument("--mu", type=float, default=MU_DEFAULT)
    ap.add_argument("--correct-nnn", type=int, default=None)
    ap.add_argument("--mode", choices=("long", "shear"), default="long")
    a = ap.parse_args()
    for e in a.eos:
        for n in a.n:
            make_ic(n, e, a.amplitude, a.mu, a.correct_nnn, a.mode)
