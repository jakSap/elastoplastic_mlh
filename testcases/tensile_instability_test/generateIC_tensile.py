#!/usr/bin/env python3
"""
Tensile-instability stretching test IC for the newMFM demonstrator.

A uniform lattice of single-material (matId=0) particles fills the periodic box
[-10, 10] x [-1, 1] at the Murnaghan reference density rho0. A step-function
stretching velocity field

    vx = -v  for x < 0,   vx = +v  for x >= 0,   vy = 0

is imposed: the two halves fly apart, the density drops below rho0, and the
Murnaghan EOS

    P = K0/n * ((rho/rho0)^n - 1)

returns P < 0 (uniform tension). With the tensile correction disabled, this is
the cleanest pure-fluid driver of the tensile instability: clumping should grow
with the cubic spline and be suppressed by the Wendland C2 kernel.

Single material, normalized units (rho0 = K0 = 1, c = 1). u is inert for
Murnaghan (P depends only on rho) and is set to 0.

Run:  python3 generateIC_tensile.py [--velocity 1.0] [--delta-p 0.05] [--plot]
"""

import os
import argparse
import numpy as np
import h5py

# ----------------------------------------------------------------------------
# Parameters (normalized units; must match config_*.info materials + box)
# ----------------------------------------------------------------------------
x_lo, x_hi = -10.0, 10.0  # periodic box in x
y_lo, y_hi = -1.0, 1.0    # periodic box in y
dp_default = 0.05          # default particle spacing
v_default = 1.0            # default bulk velocity magnitude
rho0 = 1.0                 # Murnaghan reference density (P = 0 here)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--velocity", type=float, default=v_default,
                    help="bulk velocity magnitude (default: %(default)s)")
    ap.add_argument("--delta-p", type=float, default=dp_default,
                    help="particle spacing (default: %(default)s)")
    ap.add_argument("--plot", action="store_true")
    args = ap.parse_args()
    v = args.velocity
    delta_p = args.delta_p

    out_file = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                            f"tensile_ic_v{v:g}_dp{delta_p:g}.h5")

    # Cell-centred lattice: positions sit at (lo + (i+0.5)*dp) so the wrap-around
    # spacing across each periodic seam is exactly delta_p (no double/empty row).
    nx = int(round((x_hi - x_lo) / delta_p))
    ny = int(round((y_hi - y_lo) / delta_p))
    xs = x_lo + (np.arange(nx) + 0.5) * delta_p
    ys = y_lo + (np.arange(ny) + 0.5) * delta_p
    gx, gy = np.meshgrid(xs, ys)
    gx = gx.ravel()
    gy = gy.ravel()
    N = gx.size

    x = np.zeros((N, 2))
    x[:, 0] = gx
    x[:, 1] = gy

    vel = np.zeros((N, 2))
    vel[:, 0] = np.where(gx < 0, -v, v)   # vx = -v for x<0, +v for x>=0

    mass = delta_p ** 2 * rho0  # uniform density rho0 on the lattice
    m = np.ones(N) * mass
    u = np.zeros(N)             # inert for Murnaghan
    materialId = np.zeros(N, dtype=np.int8)

    with h5py.File(out_file, "w") as h5f:
        h5f.create_dataset("x", data=x)
        h5f.create_dataset("v", data=vel)
        h5f.create_dataset("m", data=m)
        h5f.create_dataset("u", data=u)
        h5f.create_dataset("materialId", data=materialId)
        h5f.create_dataset("time", data=0.0 * m)

    print(f"N = {N}  ({nx} x {ny}),  delta_p = {delta_p}")
    print(f"box [{x_lo},{x_hi}] x [{y_lo},{y_hi}],  v = {v},  rho0 = {rho0}")
    print(f"mass/particle = {mass},  |vx| = {v}")
    print(f"Saved {out_file}")

    if args.plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(10, 2), dpi=130)
        sc = ax.scatter(x[:, 0], x[:, 1], s=1, c=vel[:, 0], cmap="coolwarm")
        ax.set_aspect("equal")
        fig.colorbar(sc, label="vx")
        ax.set_title(f"tensile stretching IC  N={N}")
        png = out_file.replace(".h5", ".png")
        plt.savefig(png, bbox_inches="tight")
        print("Saved", png)


if __name__ == "__main__":
    main()
