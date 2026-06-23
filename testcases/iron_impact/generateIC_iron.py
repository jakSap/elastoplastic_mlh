#!/usr/bin/env python3
"""
2D iron-impact IC for the newMFM demonstrator -- the same particle set as the 2D
miluphcuda plasticity benchmark (miluphcuda_bare/plasticity_impact_2d), just
written in the demonstrator's HDF5 format.

An iron bullet (disk, R = 1.25 cm) strikes an iron block (10 cm wide x 8 cm deep,
top at y = 0) straight down (-y) at 1000 m/s. Single material (matId = 0).
Elastic + von Mises plastic (Y = 2.5 GPa), Tillotson EOS -- all set in
config.info / parameter.h, not here.

SI units throughout (matches miluphcuda). dx = 2.5 mm, rho0 = 7800 kg/m^3.
The deviatoric stress starts at zero; the solver allocates/initializes S itself.

Run:  python3 generateIC_iron.py [--plot]
"""
import os
import argparse
import numpy as np
import h5py

# --------------------------------------------------------------------------
# geometry / material (identical to plasticity_impact_2d/create_ic.py)
# --------------------------------------------------------------------------
dx   = 0.0025                  # particle spacing [m] (2.5 mm)
sml0 = 2.4 * dx                # initial smoothing length: demonstrator support
                              # radius = h, matched to miluphcuda cubic-spline
                              # support 2*sml_milu = 2*1.2*dx = 2.4*dx.

block_x = (-0.05, 0.05)       # 10 cm wide
block_y = (-0.08, 0.0)        #  8 cm deep, top surface at y = 0

bullet_R   = 0.0125           # radius 1.25 cm
bullet_gap = 0.0025           # gap between bullet bottom and block top
bullet_cx  = 0.0
bullet_cy  = bullet_R + bullet_gap
impact_vel = 1000.0           # impact speed [m/s], directed -y

rho0 = 7.8e3                  # iron reference density
u0   = 0.0                   # cold start

out_file = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "iron_impact-2D.h5")


def lattice_box(xr, yr, dx):
    xs = np.arange(xr[0], xr[1] + 0.5 * dx, dx)
    ys = np.arange(yr[0], yr[1] + 0.5 * dx, dx)
    X, Y = np.meshgrid(xs, ys, indexing="ij")
    return X.ravel(), Y.ravel()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--plot", action="store_true")
    args = ap.parse_args()

    # ---- block (iron) ----
    bx, by = lattice_box(block_x, block_y, dx)
    # ---- bullet (iron disk) ----
    sx, sy = lattice_box((bullet_cx - bullet_R, bullet_cx + bullet_R),
                         (bullet_cy - bullet_R, bullet_cy + bullet_R), dx)
    inside = (sx - bullet_cx)**2 + (sy - bullet_cy)**2 <= bullet_R**2
    sx, sy = sx[inside], sy[inside]

    N_b, N_p = bx.size, sx.size
    N = N_b + N_p

    x = np.zeros((N, 2))
    x[:N_b, 0] = bx; x[:N_b, 1] = by
    x[N_b:, 0] = sx; x[N_b:, 1] = sy

    v = np.zeros((N, 2))
    v[N_b:, 1] = -impact_vel               # only the bullet moves, in -y

    m   = np.full(N, rho0 * dx**2)         # 2D mass per unit thickness
    u   = np.full(N, u0)
    sml = np.full(N, sml0)
    matId = np.zeros(N, dtype=np.int8)

    with h5py.File(out_file, "w") as f:
        f.create_dataset("x", data=x)
        f.create_dataset("v", data=v)
        f.create_dataset("m", data=m)
        f.create_dataset("u", data=u)
        f.create_dataset("sml", data=sml)
        f.create_dataset("materialId", data=matId)
        f.create_dataset("time", data=0.0 * m)

    print(f"block particles : {N_b}")
    print(f"bullet particles: {N_p}")
    print(f"total particles : {N}")
    print(f"mass/particle   : {rho0*dx**2:.6e} kg, sml0 = {sml0:.4e} m")
    print(f"Saved {out_file}")

    if args.plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(6, 6), dpi=130)
        ax.scatter(x[:N_b, 0], x[:N_b, 1], s=3, c="0.4", label="block")
        ax.scatter(x[N_b:, 0], x[N_b:, 1], s=3, c="C3", label="bullet")
        ax.set_aspect("equal"); ax.legend()
        ax.set_title(f"iron impact IC  N={N}")
        png = out_file.replace(".h5", ".png")
        plt.savefig(png, bbox_inches="tight"); print("Saved", png)


if __name__ == "__main__":
    main()
