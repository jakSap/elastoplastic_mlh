#!/usr/bin/env python3

import argparse
import numpy as np
import matplotlib.pyplot as plt
import h5py as h5

"""
Generates initial conditions for a single 2D elastic block of size 0.5x0.5,
centered in a 1x1 domain (coordinates in [-0.5, 0.5] x [-0.5, 0.5]).
Particles are placed on a regular lattice with spacing delta_p and
given a uniform initial velocity in the x-direction.

Intended for use with the Murnaghan EOS (EOS=1). The default internal energy
u=0 corresponds to the stress-free equilibrium state at rho = rho0.
"""

DIM = 2

if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Create an initial condition HDF5 file for the single elastic block test case.")
    parser.add_argument("--delta_p", "-d", metavar="float", type=float,
                        default=0.025, help="particle spacing (default: 0.025)")
    parser.add_argument("--velocity", "-v", metavar="float", type=float,
                        default=0.059, help="x-velocity of the block (default: 0.059)")
    parser.add_argument("--vy", metavar="float", type=float,
                        default=0.0, help="y-velocity of the block (default: 0.0)")
    parser.add_argument("--radius", "-r", metavar="float", type=float,
                        default=0.25, help="half-size of the block (default: 0.25)")
    parser.add_argument("--density", metavar="float", type=float,
                        default=1.0, help="reference particle density (default: 1.0)")
    parser.add_argument("--u", metavar="float", type=float,
                        default=0.0, help="specific internal energy (default: 0.0, equilibrium for Murnaghan at rho=rho0)")
    parser.add_argument("--vx-gradient", action="store_true",
                        help="apply a linear vx gradient: vx=0 at x=0, ramping to full velocity at x=+-0.25")
    parser.add_argument("--fillUp", "-f", action="store_true",
                        help="fill up coordinates to 3D with z=0")
    parser.add_argument("--plot", action="store_true",
                        help="plot initial configuration")

    args = parser.parse_args()

    delta_p = args.delta_p
    v_p     = args.velocity
    vy_p    = args.vy
    density = args.density
    u_const = args.u
    fillUp  = args.fillUp

    # Block centered at origin within the 1x1 domain
    half = args.radius
    block_size = 2.0 * half

    # 2D particle mass from area per particle
    mass = delta_p**2 * density

    print("Generating single elastic block initial conditions ...")
    print(f"  delta_p = {delta_p}, block = {block_size} x {block_size}, centered at (0, 0)")
    print(f"  velocity = {v_p}, density = {density}, u = {u_const}")
    print(f"  mass per particle = {mass}")

    # Regular lattice: particle centres at -half + delta_p/2, ..., half - delta_p/2
    N_per_side = int(round(block_size / delta_p))
    coords = np.arange(N_per_side) * delta_p - half + delta_p / 2.0

    xx, yy = np.meshgrid(coords, coords, indexing='ij')
    x_flat = xx.ravel()
    y_flat = yy.ravel()
    N = len(x_flat)

    if fillUp:
        pos = np.column_stack([x_flat, y_flat, np.zeros(N)])
        vel = np.zeros((N, DIM + 1))
    else:
        pos = np.column_stack([x_flat, y_flat])
        vel = np.zeros((N, DIM))

    if args.vx_gradient:
        # Antisymmetric triangle wave ("edgy sine"):
        # positive hump on [-0.25, 0], negative hump on [0, 0.25], sum = 0
        quarter = half / 2.0  # 0.125
        envelope = (quarter - np.abs(np.abs(x_flat) - quarter)) / quarter
        vel[:, 0] = v_p * envelope * np.sign(-x_flat)
    else:
        vel[:, 0] = v_p

    vel[:, 1] = vy_p

    mean_vx = np.mean(vel[:, 0])
    print('Mean vx: ', mean_vx)
    mean_vy = np.mean(vel[:, 1])
    print('Mean vy: ', mean_vy)

    m          = np.ones(N) * mass
    u          = np.ones(N) * u_const
    materialId = np.zeros(N, dtype=np.int8)

    suffix   = "3D" if fillUp else "2D"
    filename = f"single_block_deltap{delta_p}-{suffix}.h5"

    outH5 = h5.File(filename, "w")
    outH5.create_dataset("x",          data=pos)
    outH5.create_dataset("v",          data=vel)
    outH5.create_dataset("m",          data=m)
    outH5.create_dataset("u",          data=u)
    outH5.create_dataset("materialId", data=materialId)
    outH5.create_dataset("time",       data=0.0 * m)
    outH5.close()

    print(f"Number of particles: {N}  ({N_per_side} x {N_per_side})")
    print(f"... done. Output written to {filename}")

    if args.plot:
        print("Plotting initial configuration ...")
        fig, ax = plt.subplots(figsize=(6, 6), dpi=150)
        sc = ax.scatter(pos[:, 0], pos[:, 1], c=m, s=2.)
        fig.colorbar(sc, ax=ax, label="m[i]")
        # draw domain and block boundaries
        domain = plt.Rectangle((-0.5, -0.5), 1.0, 1.0,
                                linewidth=1, edgecolor="grey", facecolor="none", linestyle="--")
        block  = plt.Rectangle((-half, -half), block_size, block_size,
                                linewidth=1, edgecolor="red",  facecolor="none", linestyle="-")
        ax.add_patch(domain)
        ax.add_patch(block)
        ax.set_title("Single elastic block initial conditions")
        ax.set_xlabel("$x$")
        ax.set_ylabel("$y$")
        ax.set_aspect("equal")
        ax.set_xlim(-0.6, 0.6)
        ax.set_ylim(-0.6, 0.6)
        plt.tight_layout()
        plotname = f"single_block_deltap{delta_p}-{suffix}.png"
        plt.savefig(plotname)
        plt.show()
        print(f"... saved to {plotname}")
