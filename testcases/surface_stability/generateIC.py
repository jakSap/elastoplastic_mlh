#!/usr/bin/env python3
"""
Initial conditions for the surface-stability test case.

Generates a stationary 2D elastic offset-lune: an outer disk with a smaller,
off-centre disk removed. The inner disk is required to lie fully inside the
outer one (offset + R_in < R_out), so the lune has finite thickness
everywhere -- no cusping horns that would starve particles of neighbors.

The shape still carries both surface types in a single connected body:

  - a convex outer circle
  - a concave inner circle, closer to one side (thin vs thick lobes)

Particles carry zero velocity and zero internal energy (stress-free Murnaghan
equilibrium), so any surface motion observed in the run is a pure numerical
artifact of the scheme at free surfaces.

Also writes a per-particle 'region' label (0 = interior, 1 = convex outer,
2 = concave inner, 3 = thin-lobe near-surface) so the analysis script can
split diagnostics by local geometry.
"""

import argparse
import numpy as np
import h5py as h5
import matplotlib.pyplot as plt


def build_crescent(delta_p, R_out, R_in, offset):
    """Return (pos, region_label) for particles on a regular lattice that lie
    inside the outer disk but outside the inner (offset) disk.

    region labels:
        0 = interior
        1 = on convex outer arc (within ~2 dp)
        2 = on concave inner arc (within ~2 dp)
        3 = near a horn (both outer and inner arc within tolerance)
    """
    # bounding box of the outer disk
    N_side = int(np.ceil(2 * R_out / delta_p)) + 2
    coords = (np.arange(N_side) - (N_side - 1) / 2.0) * delta_p
    xx, yy = np.meshgrid(coords, coords, indexing='ij')
    x = xx.ravel()
    y = yy.ravel()

    r_outer = np.sqrt(x * x + y * y)
    r_inner = np.sqrt((x - offset) ** 2 + y * y)

    inside = (r_outer <= R_out) & (r_inner >= R_in)

    x = x[inside]
    y = y[inside]
    r_outer = r_outer[inside]
    r_inner = r_inner[inside]

    tol = 2.0 * delta_p
    on_outer = (R_out - r_outer) < tol
    on_inner = (r_inner - R_in) < tol

    region = np.zeros(len(x), dtype=np.int32)
    region[on_outer] = 1
    region[on_inner] = 2
    # thin lobe: the side opposite to the offset direction is the thinnest;
    # flag outer-boundary particles there so we can see whether thin-lobe
    # convex surface frays differently from the thick-lobe one.
    thin_side = on_outer & (x > 0)
    region[thin_side] = 3

    return np.column_stack([x, y]), region


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generate a crescent IC for the surface stability test.")
    parser.add_argument("--delta_p", "-d", type=float, default=0.05,
                        help="particle spacing (default: 0.05)")
    parser.add_argument("--R_out", type=float, default=3.0,
                        help="outer disk radius (default: 3.0)")
    parser.add_argument("--R_in", type=float, default=2.0,
                        help="radius of the subtracted inner disk (default: 2.0)")
    parser.add_argument("--offset", type=float, default=0.5,
                        help="x-offset of the subtracted disk (default: 0.5, "
                             "must satisfy offset + R_in < R_out)")
    parser.add_argument("--density", type=float, default=1.0,
                        help="reference density (default: 1.0)")
    parser.add_argument("--plot", action="store_true",
                        help="save a preview plot")
    args = parser.parse_args()

    if args.offset + args.R_in >= args.R_out:
        raise SystemExit(
            f"offset + R_in = {args.offset + args.R_in} must be < R_out = {args.R_out}; "
            "otherwise the shape cusps at horn tips and particles lose all neighbors.")
    pos, region = build_crescent(args.delta_p, args.R_out, args.R_in, args.offset)
    N = len(pos)

    mass = args.density * args.delta_p ** 2
    m = np.full(N, mass)
    u = np.zeros(N)              # Murnaghan stress-free equilibrium
    vel = np.zeros((N, 2))       # stationary
    materialId = np.zeros(N, dtype=np.int8)

    n_conv = int((region == 1).sum())
    n_conc = int((region == 2).sum())
    n_thin = int((region == 3).sum())
    n_int  = int((region == 0).sum())

    print(f"Offset lune: R_out={args.R_out}, R_in={args.R_in}, offset={args.offset}")
    print(f"  delta_p = {args.delta_p}")
    print(f"  N = {N}  (interior {n_int}, convex outer (thick lobe) {n_conv}, "
          f"concave inner {n_conc}, convex outer (thin lobe) {n_thin})")
    print(f"  mass per particle = {mass}")

    filename = f"lune_deltap{args.delta_p}-2D.h5"
    with h5.File(filename, "w") as h5f:
        h5f.create_dataset("x",           data=pos)
        h5f.create_dataset("v",           data=vel)
        h5f.create_dataset("m",           data=m)
        h5f.create_dataset("u",           data=u)
        h5f.create_dataset("materialId",  data=materialId)
        h5f.create_dataset("time",        data=0.0 * m)
        h5f.create_dataset("region",      data=region)

    print(f"Wrote {filename}")

    if args.plot:
        fig, ax = plt.subplots(figsize=(6, 6), dpi=150)
        colors = np.array(['#888', '#1f77b4', '#d62728', '#2ca02c'])
        ax.scatter(pos[:, 0], pos[:, 1], c=colors[region], s=3)
        ax.set_aspect('equal')
        ax.set_title("Offset lune IC: interior / convex(thick) / concave / convex(thin)")
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        plt.tight_layout()
        plotname = filename.replace('.h5', '.png')
        fig.savefig(plotname)
        print(f"Wrote {plotname}")
