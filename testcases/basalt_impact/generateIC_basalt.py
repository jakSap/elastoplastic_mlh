#!/usr/bin/env python3
"""
2D basalt impact IC (miluphcuda's brittle-fragmentation benchmark, after
Nakamura & Fujiwara 1991 / Benz & Asphaug 1994), adapted to the newMFM
demonstrator.

A circular basalt TARGET disk at the origin is struck by a small, fast basalt
PROJECTILE disk coming in from the left (+x velocity). Single material
(matId = 0; the demonstrator has no multi-material ICs yet), normalized units
consistent with the colliding-rings testcase (rho0 = A = c = 1). The target
starts cold and at rest (u = 0, v = 0); the projectile carries v_imp.

Weibull activation flaws for the Grady-Kipp model are written into the IC
(/numFlaws, /flaws); they are ignored by the solver unless FRAGMENTATION is on.
max_flaws MUST match MAX_NUM_FLAWS in demonstrator/include/parameter.h.

Run:  python3 generateIC_basalt.py [--plot]
"""

import os
import sys
import argparse

import numpy as np

# shared Grady-Kipp Weibull flaw generator (testcases/weibull_flaws.py)
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from weibull_flaws import write_flaws

import h5py

# ----------------------------------------------------------------------------
# Parameters (normalized units)
# ----------------------------------------------------------------------------
delta_p   = 0.10          # particle spacing -> ~8k particles
R_target  = 5.0           # target disk radius
R_proj    = 0.8           # projectile disk radius
gap       = 0.1           # initial gap between projectile and target
v_imp     = 0.2           # impact speed; c = sqrt(A/rho) = 1, so Mach ~0.2.
                          # The scheme is tuned for the gentle (Mach<~0.3) regime
                          # of the colliding-rings benchmark; faster impacts
                          # collapse the CFL timestep.
rho0      = 1.0
u_const   = 0.0           # cold start; for the linear Tillotson (a=b=B=0) P
                          # depends only on density, so u is inert here.

# Weibull flaws -- calibrated (see run_sweep notes) so flaws activate in the
# high-tension region of this impact but not everywhere.
weibull_k = 5.0e12
weibull_m = 9.5           # basalt value (Nakamura)
max_flaws = 1             # MUST match MAX_NUM_FLAWS in parameter.h
flaw_seed = 0

out_file  = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                         "basalt_impact-2D.h5")


def disk(cx, cy, R, dp):
    """Grid points inside a disk of radius R centered at (cx, cy)."""
    n = int(np.ceil(2 * R / dp)) + 2
    xs = cx + (np.arange(-n, n + 1)) * dp
    ys = cy + (np.arange(-n, n + 1)) * dp
    gx, gy = np.meshgrid(xs, ys)
    gx = gx.ravel(); gy = gy.ravel()
    mask = (gx - cx) ** 2 + (gy - cy) ** 2 <= R * R
    return gx[mask], gy[mask]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--plot", action="store_true")
    args = ap.parse_args()

    # target at origin, at rest
    tx, ty = disk(0.0, 0.0, R_target, delta_p)
    # projectile to the left, moving +x
    cx = -(R_target + gap + R_proj)
    px, py = disk(cx, 0.0, R_proj, delta_p)

    N_t, N_p = tx.size, px.size
    N = N_t + N_p

    x = np.zeros((N, 2))
    x[:N_t, 0] = tx; x[:N_t, 1] = ty
    x[N_t:, 0] = px; x[N_t:, 1] = py

    v = np.zeros((N, 2))
    v[N_t:, 0] = v_imp                      # only the projectile moves

    mass = delta_p ** 2 * rho0
    m = np.ones(N) * mass
    u = np.ones(N) * u_const
    materialId = np.zeros(N, dtype=np.int8)  # single material

    h5f = h5py.File(out_file, "w")
    h5f.create_dataset("x", data=x)
    h5f.create_dataset("v", data=v)
    h5f.create_dataset("m", data=m)
    h5f.create_dataset("u", data=u)
    h5f.create_dataset("materialId", data=materialId)
    h5f.create_dataset("time", data=0.0 * m)
    write_flaws(h5f, m / rho0, weibull_k, weibull_m, max_flaws, flaw_seed)
    h5f.close()

    print(f"target N = {N_t}, projectile N = {N_p}, total N = {N}")
    print(f"v_imp = {v_imp}, projectile center x = {cx:.3f}")
    print(f"Saved {out_file}")

    if args.plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, axp = plt.subplots(figsize=(7, 6), dpi=130)
        axp.scatter(x[:N_t, 0], x[:N_t, 1], s=2, c="0.4", label="target")
        axp.scatter(x[N_t:, 0], x[N_t:, 1], s=2, c="C3", label="projectile")
        axp.set_aspect("equal"); axp.legend()
        axp.set_title(f"basalt impact IC  N={N}  v_imp={v_imp}")
        png = out_file.replace(".h5", ".png")
        plt.savefig(png); print("Saved", png)


if __name__ == "__main__":
    main()
