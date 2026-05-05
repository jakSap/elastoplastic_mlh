#!/usr/bin/env python3

import argparse
import numpy as np
import matplotlib.pyplot as plt
import h5py as h5

"""
Creates two 2-D rings (z = 0) for the colliding rings testcase using concentric
circles of particles.  Each circle at radius r carries N = round(2*pi*r/delta_p)
evenly-spaced particles; alternate circles are offset by half an arc-step to
approximate hexagonal close-packing.

Reference: Gray, Monaghan, Swift – SPH elastic dynamics, CMAME (2001).
"""

DIM = 2


def make_ring(r_inner, r_outer, delta_p):
    """Return (x, y) arrays for one ring built from concentric circles."""
    xs, ys = [], []
    radii = np.arange(r_inner, r_outer + 1e-9, delta_p)
    for idx, r in enumerate(radii):
        n = max(1, round(2 * np.pi * r / delta_p))
        offset = (0.5 * delta_p / r) if (idx % 2 == 1) else 0.0
        angles = np.linspace(0, 2 * np.pi, n, endpoint=False) + offset
        xs.extend(r * np.cos(angles))
        ys.extend(r * np.sin(angles))
    return np.array(xs), np.array(ys)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Create an initial condition HDF5 file for the colliding rings "
                    "test case using concentric-circle particle placement.")
    parser.add_argument("--delta_p", "-d", metavar="float", type=float,
                        default=0.2, help="particle spacing (default: 0.2)")
    parser.add_argument("--r_inner", metavar="float", type=float,
                        default=3.0, help="inner ring radius (default: 3.0)")
    parser.add_argument("--r_outer", metavar="float", type=float,
                        default=4.0, help="outer ring radius (default: 4.0)")
    parser.add_argument("--distance", "-D", metavar="float", type=float,
                        default=2.0,
                        help="gap between the two rings' outer surfaces (default: 2.0)")
    parser.add_argument("--velocity", "-v", metavar="float", type=float,
                        default=0.059, help="projected speed of each ring (default: 0.059)")
    parser.add_argument("--density", metavar="float", type=float,
                        default=1.0, help="particle density (default: 1.0)")
    parser.add_argument("--pressure", "-P", metavar="float", type=float,
                        default=2.5, help="pressure constant (default: 2.5)")
    parser.add_argument("--gamma", metavar="float", type=float,
                        default=5.0/3.0, help="adiabatic index (default: 5/3)")
    parser.add_argument("--fillUp", "-f", action="store_true",
                        help="fill up coordinates to 3D with z=0")
    parser.add_argument("--plot", action="store_true",
                        help="plot initial configuration")

    args = parser.parse_args()

    delta_p  = args.delta_p
    r_inner  = args.r_inner
    r_outer  = args.r_outer
    distance = args.distance
    v_p      = args.velocity
    density  = args.density
    P        = args.pressure
    gamma    = args.gamma
    fillUp   = args.fillUp

    # centre-to-centre separation so outer surfaces are 'distance' apart
    x_shift = r_outer + distance / 2.0

    mass    = delta_p**2 * density
    u_const = P / ((gamma - 1.0) * density)

    print("Generating colliding rings initial conditions (concentric method) ...")
    print(f"  delta_p = {delta_p}, r_inner = {r_inner}, r_outer = {r_outer}")
    print(f"  distance = {distance}  →  x_shift = {x_shift}")
    print(f"  v_p = {v_p}, density = {density}")
    print(f"  mass = {mass}, u = {u_const}")

    xs, ys = make_ring(r_inner, r_outer, delta_p)
    N = len(xs)

    dim_out = DIM + 1 if fillUp else DIM

    r_ring1 = np.zeros((N, dim_out))
    r_ring2 = np.zeros((N, dim_out))
    v1      = np.zeros((N, dim_out))
    v2      = np.zeros((N, dim_out))

    r_ring1[:, 0] = xs - x_shift
    r_ring1[:, 1] = ys
    v1[:, 0]      = v_p

    r_ring2[:, 0] = xs + x_shift
    r_ring2[:, 1] = ys
    v2[:, 0]      = -v_p

    r_final = np.concatenate((r_ring1, r_ring2))
    v_final = np.concatenate((v1, v2))

    m          = np.ones(2 * N) * mass
    materialId = np.zeros(2 * N, dtype=np.int8)
    u          = np.ones(2 * N) * u_const

    suffix   = "3D" if fillUp else "2D"
    filename = f"rings_concentric_deltap{delta_p}-{suffix}.h5"

    outH5 = h5.File(filename, "w")
    outH5.create_dataset("x",          data=r_final)
    outH5.create_dataset("v",          data=v_final)
    outH5.create_dataset("m",          data=m)
    outH5.create_dataset("u",          data=u)
    outH5.create_dataset("materialId", data=materialId)
    outH5.create_dataset("time",       data=0.0 * m)
    outH5.close()

    print(f"Number of particles per ring: {N}  (total: {2 * N})")
    print(f"... done. Output written to {filename}")

    if args.plot:
        print("Plotting initial configuration ...")
        fig, ax = plt.subplots(figsize=(10, 5), dpi=200)
        ax.scatter(r_final[:, 0], r_final[:, 1], s=1.5, linewidths=0)
        ax.set_aspect("equal")
        ax.set_title("Colliding rings IC – concentric method")
        ax.set_xlabel("$x$")
        ax.set_ylabel("$y$")
        plt.tight_layout()
        plotname = f"rings_concentric_deltap{delta_p}-{suffix}.png"
        plt.savefig(plotname)
        print(f"... saved to {plotname}")
