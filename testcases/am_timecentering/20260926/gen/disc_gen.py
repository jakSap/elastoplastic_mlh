#!/usr/bin/env python3
"""Rigidly rotating elastic disc -- the angular-momentum zero test (TEST_PLAN.md A.1).

A 2D elastic disc of radius R, uniform density, free surface, set spinning as a rigid body:

    v = omega x r = omega (-y, x),   rho = rho_0,   S = 0,   P = 0

Angular momentum L_z = sum_i m_i (x_i v_y - y_i v_x) is exactly conserved by the continuum
equations, and the discrete scheme should conserve it too. It does not: the meshless
momentum exchange is along the line connecting two particles only if the face geometry is
symmetric, and the error is what RECON_FACE_DV 1 exists to reduce (68-99.8 % loss without it
in the miluphcuda disc study, < 1 % with it).

The disc is not in elastic equilibrium -- rigid rotation needs a centrifugal stress that S = 0
does not provide -- so it breathes slightly. That is intended and harmless: the test is about
L_z, which is conserved whatever the disc does internally, and the breathing keeps the free
surface and the tensile branch active, which is what makes it a solid-dynamics test rather
than a kinematic one.

Material: Hopkins-canonical Tillotson (rho0 = A = 1, everything else 0 -> P = A (rho/rho0 - 1),
K = 1) with mu = 0.22, the same material the colliding-rings comparison uses, so the two
gate tests share a material and only the geometry differs.

Usage:
    ./generateDiscIC.py                 # default: R = 1, d = 0.025, omega = 0.2
    ./generateDiscIC.py --omega 0.05 -o disc_w005.h5
"""

import argparse

import h5py
import numpy as np

RHO0 = 1.0
MU = 0.22
# u must stay > 0 (u <= 0 logs an ERROR every step). With this Tillotson (a = b = B = 0)
# the pressure does not depend on u at all, but u still evolves with the energy equation and
# a free surface doing expansion work drives it down -- u = 0.01 produced tens of thousands
# of negative-energy errors. The colliding-rings IC uses u = 1.0; same value here, so both
# gate tests sit in the same state and neither spends the run logging errors.
U_BAR = 1.0


def make_disc(radius, spacing, omega, out):
    n = int(np.ceil(radius / spacing)) + 2
    i = np.arange(-n, n + 1) * spacing
    X, Y = np.meshgrid(i, i, indexing="ij")
    inside = X * X + Y * Y <= radius * radius
    pos = np.column_stack([X[inside], Y[inside]])
    N = len(pos)

    vel = np.column_stack([-omega * pos[:, 1], omega * pos[:, 0]])
    m = np.full(N, RHO0 * spacing * spacing)
    u = np.full(N, U_BAR)
    matId = np.zeros(N, dtype=np.int8)
    zero = np.zeros(N)

    Lz = float(np.sum(m * (pos[:, 0] * vel[:, 1] - pos[:, 1] * vel[:, 0])))
    # analytic rigid-disc value for a uniform disc of mass M: L_z = 1/2 M R^2 omega
    M = float(m.sum())
    Lz_exact = 0.5 * M * radius * radius * omega

    with h5py.File(out, "w") as f:
        f.create_dataset("x", data=pos)
        f.create_dataset("v", data=vel)
        f.create_dataset("m", data=m)
        f.create_dataset("u", data=u)
        f.create_dataset("materialId", data=matId)
        for name in ("Sxx", "Sxy", "Syy"):
            f.create_dataset(name, data=zero)
        for key, val in (("radius", radius), ("spacing", spacing), ("omega", omega),
                         ("rho0", RHO0), ("mu", MU), ("N", N), ("mass", M),
                         ("Lz", Lz), ("Lz_rigid", Lz_exact),
                         ("period", 2.0 * np.pi / omega)):
            f.attrs[key] = val

    print("Wrote {}: N = {}, R = {:g}, d = {:g}, omega = {:g}, M = {:.6g}, "
          "L_z = {:.9g} (rigid-disc value {:.9g}, ratio {:.6f}), T_rot = {:.4f}".format(
              out, N, radius, spacing, omega, M, Lz, Lz_exact, Lz / Lz_exact,
              2 * np.pi / omega))
    return out


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-R", "--radius", type=float, default=1.0)
    ap.add_argument("-d", "--spacing", type=float, default=0.025)
    ap.add_argument("-w", "--omega", type=float, default=0.2)
    ap.add_argument("-o", "--out", default=None)
    args = ap.parse_args()
    out = args.out or "disc_R{:g}_d{:g}_w{:g}.h5".format(args.radius, args.spacing, args.omega)
    make_disc(args.radius, args.spacing, args.omega, out)
