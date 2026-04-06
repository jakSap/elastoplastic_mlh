#!/usr/bin/env python3

import argparse
import numpy as np
import matplotlib.pyplot as plt
import h5py as h5

"""
Generates initial conditions for two 2D elastic blocks of size 0.5x0.5,
placed symmetrically about x=0 in a 2x1 domain ([-1,1] x [-0.5,0.5]).
Left block is centered at (-0.5, 0), right block at (0.5, 0),
with a gap of 0.5 between them. Blocks move towards each other and
meet at x=0.

Intended for use with the Murnaghan EOS (EOS=1). The default internal energy
u=0 corresponds to the stress-free equilibrium state at rho = rho0.
"""

DIM = 2

if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Create an initial condition HDF5 file for the two colliding elastic blocks test case.")
    parser.add_argument("--delta_p", "-d", metavar="float", type=float,
                        default=0.025, help="particle spacing (default: 0.025)")
    parser.add_argument("--velocity", "-v", metavar="float", type=float,
                        default=0.059, help="x-velocity magnitude of each block (default: 0.059)")
    parser.add_argument("--closing_velocity", "-c", metavar="float", type=float,
                        default=None, help="total closing velocity of the two blocks (= 2 * velocity); overrides --velocity if given")
    parser.add_argument("--density", metavar="float", type=float,
                        default=1.0, help="reference particle density (default: 1.0)")
    parser.add_argument("--u", metavar="float", type=float,
                        default=0.0, help="specific internal energy (default: 0.0, equilibrium for Murnaghan at rho=rho0)")
    parser.add_argument("--fillUp", "-f", action="store_true",
                        help="fill up coordinates to 3D with z=0")
    parser.add_argument("--plot", action="store_true",
                        help="plot initial configuration")

    args = parser.parse_args()

    delta_p = args.delta_p
    density = args.density
    u_const = args.u
    fillUp  = args.fillUp

    if args.closing_velocity is not None:
        v_p = args.closing_velocity / 2.0
    else:
        v_p = args.velocity

    # Each block: 0.5 x 0.5
    block_size = 0.5
    half = block_size / 2.0  # 0.25

    # Block centers: left at (-0.5, 0), right at (0.5, 0)
    # Gap between the two blocks: 0.5 (from x=-0.25 to x=0.25)
    offset = 0.5

    # 2D particle mass from area per particle
    mass = delta_p**2 * density

    print("Generating two colliding elastic blocks initial conditions ...")
    print(f"  delta_p = {delta_p}, block = {block_size} x {block_size}")
    print(f"  left block center=(-{offset}, 0), right block center=({offset}, 0)")
    print(f"  individual velocity = {v_p}, closing velocity = {2*v_p}")
    print(f"  density = {density}, u = {u_const}")
    print(f"  mass per particle = {mass}")

    # Regular lattice for one block
    N_per_side = int(round(block_size / delta_p))
    coords = np.arange(N_per_side) * delta_p - half + delta_p / 2.0

    xx, yy = np.meshgrid(coords, coords, indexing='ij')
    x_flat = xx.ravel()
    y_flat = yy.ravel()
    N_block = len(x_flat)

    # Left block: centered at (-offset, 0), velocity +v_p in x
    x_left = x_flat - offset
    y_left = y_flat

    # Right block: centered at (+offset, 0), velocity -v_p in x
    x_right = x_flat + offset
    y_right = y_flat

    x_all = np.concatenate([x_left, x_right])
    y_all = np.concatenate([y_left, y_right])
    N = len(x_all)

    if fillUp:
        pos = np.column_stack([x_all, y_all, np.zeros(N)])
        vel = np.zeros((N, DIM + 1))
    else:
        pos = np.column_stack([x_all, y_all])
        vel = np.zeros((N, DIM))

    # Left block: +v_p, right block: -v_p
    vel[:N_block, 0] =  v_p
    vel[N_block:, 0] = -v_p

    m          = np.ones(N) * mass
    u          = np.ones(N) * u_const
    materialId = np.zeros(N, dtype=np.int8)
    materialId[N_block:] = 1  # right block gets material id 1

    suffix   = "3D" if fillUp else "2D"
    filename = f"two_blocks_deltap{delta_p}-{suffix}.h5"

    outH5 = h5.File(filename, "w")
    outH5.create_dataset("x",          data=pos)
    outH5.create_dataset("v",          data=vel)
    outH5.create_dataset("m",          data=m)
    outH5.create_dataset("u",          data=u)
    outH5.create_dataset("materialId", data=materialId)
    outH5.create_dataset("time",       data=0.0 * m)
    outH5.close()

    print(f"Number of particles: {N}  ({N_block} per block, {N_per_side} x {N_per_side} each)")
    print(f"... done. Output written to {filename}")

    if args.plot:
        print("Plotting initial configuration ...")
        fig, ax = plt.subplots(figsize=(10, 5), dpi=150)
        colors = np.where(materialId == 0, "blue", "red")
        ax.scatter(pos[:, 0], pos[:, 1], c=colors, s=2.)
        # draw domain and block boundaries
        domain = plt.Rectangle((-1.0, -0.5), 2.0, 1.0,
                                linewidth=1, edgecolor="grey", facecolor="none", linestyle="--")
        blk_l  = plt.Rectangle((-offset - half, -half), block_size, block_size,
                                linewidth=1, edgecolor="blue", facecolor="none", linestyle="-")
        blk_r  = plt.Rectangle(( offset - half, -half), block_size, block_size,
                                linewidth=1, edgecolor="red",  facecolor="none", linestyle="-")
        ax.add_patch(domain)
        ax.add_patch(blk_l)
        ax.add_patch(blk_r)
        ax.set_title("Two colliding elastic blocks initial conditions")
        ax.set_xlabel("$x$")
        ax.set_ylabel("$y$")
        ax.set_aspect("equal")
        ax.set_xlim(-1.1, 1.1)
        ax.set_ylim(-0.6, 0.6)
        plt.tight_layout()
        plotname = f"two_blocks_deltap{delta_p}-{suffix}.png"
        plt.savefig(plotname)
        print(f"... saved to {plotname}")
