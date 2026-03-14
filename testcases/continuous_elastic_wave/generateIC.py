#!/usr/bin/env python3

import argparse
import numpy as np
import h5py as h5
import matplotlib.pyplot as plt

def cubic_spline(q):
    """Cubic spline kernel W(q) for 0 <= q <= 2."""
    w = np.zeros_like(q)
    mask1 = (q >= 0) & (q <= 1)
    mask2 = (q > 1) & (q <= 2)
    w[mask1] = 1 - 1.5 * q[mask1]**2 + 0.75 * q[mask1]**3
    w[mask2] = 0.25 * (2 - q[mask2])**3
    return w

if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Create an initial condition HDF5 file for the continuous elastic wave test case.")
    parser.add_argument("--a", "-a", metavar="float", type=float, required=True,
                        help="box width (x-extent)")
    parser.add_argument("--b", "-b", metavar="float", type=float, required=True,
                        help="box height (y-extent)")
    parser.add_argument("--h", metavar="float", type=float, default=0.0,
                        help="strip width in x-direction (0 = single-particle strip)")
    parser.add_argument("--vx", metavar="float", type=float, required=True,
                        help="maximum x-velocity in the strip")
    parser.add_argument("--d", "-d", metavar="float", type=float, required=True,
                        help="inter-particle spacing (square lattice)")
    parser.add_argument("--feather", action="store_true",
                        help="use cubic spline kernel to interpolate vx across strip")
    parser.add_argument("--output", "-o", metavar="string", type=str, default=None,
                        help="output filename (auto-generated if omitted)")
    parser.add_argument("--plot", action="store_true",
                        help="show matplotlib plot of IC")
    parser.add_argument("--vbulk", metavar="float", type = float, default = 0.0,
                        help = "bulk velocity, default 0")
    args = parser.parse_args()

    a = getattr(args, 'a')
    b = getattr(args, 'b')
    h = getattr(args, 'h')
    vx_max = args.vx
    d = getattr(args, 'd')

    # Generate 2D square lattice filling [0, a) x [0, b)
    xs = np.arange(0, a, d)
    ys = np.arange(0, b, d)
    xx, yy = np.meshgrid(xs, ys, indexing='ij')
    pos = np.column_stack([xx.ravel(), yy.ravel()])
    N = len(pos)

    # Uniform properties
    rho = 1.0
    u_val = 1.0

    # Velocity assignment
    x_center = a / 2.0
    vx_arr = np.zeros(N)
    dx = np.abs(pos[:, 0] - x_center)

    if h == 0:
        # Single-particle strip: column closest to center
        unique_x = np.unique(pos[:, 0])
        closest_x = unique_x[np.argmin(np.abs(unique_x - x_center))]
        mask = pos[:, 0] == closest_x
        vx_arr[mask] = vx_max
    elif not args.feather:
        # Constant strip
        mask = dx <= h / 2.0
        vx_arr[mask] = vx_max
    else:
        # Feathered strip with cubic spline
        q = 2.0 * dx / h
        vx_arr = vx_max * cubic_spline(q)

    vx_arr += args.vbulk
    vel = np.column_stack([vx_arr, np.zeros(N)])

    m_arr = np.full(N, rho * d**2)
    u_arr = np.full(N, u_val)
    mat = np.zeros(N, dtype=np.int8)
    t_arr = np.zeros(N)

    # Optional plot
    if args.plot:
        fig, axes = plt.subplots(1, 2, figsize=(12, 5), dpi=150)

        # 2D scatter plot colored by vx
        sc = axes[0].scatter(pos[:, 0], pos[:, 1], c=vx_arr, s=4., cmap='coolwarm')
        fig.colorbar(sc, ax=axes[0], label="$v_x$")
        axes[0].set_aspect('equal')
        axes[0].set_xlabel("$x$")
        axes[0].set_ylabel("$y$")
        axes[0].set_title("Elastic wave IC")

        # 1D velocity distribution: vx vs x (one row of particles)
        unique_y = np.unique(pos[:, 1])
        mid_y = unique_y[len(unique_y) // 2]
        row_mask = pos[:, 1] == mid_y
        axes[1].plot(pos[row_mask, 0], vx_arr[row_mask], 'o-', ms=2)
        axes[1].set_xlabel("$x$")
        axes[1].set_ylabel("$v_x$")
        axes[1].set_title("Velocity distribution (midline)")
        plt.tight_layout()
        plt.show()

    # Write HDF5
    if args.output:
        filename = args.output
    else:
        filename = f"elastic_wave_d{d}.h5"

    with h5.File(filename, "w") as f:
        f.create_dataset("x",          data=pos)
        f.create_dataset("v",          data=vel)
        f.create_dataset("m",          data=m_arr)
        f.create_dataset("u",          data=u_arr)
        f.create_dataset("materialId", data=mat)
        f.create_dataset("time",       data=t_arr)

    print(f"Generated {N} particles → {filename}")
