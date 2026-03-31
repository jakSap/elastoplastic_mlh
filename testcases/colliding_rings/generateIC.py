#!/usr/bin/env python3

import argparse
import numpy as np
import matplotlib.pyplot as plt
import h5py as h5

"""
This program creates two 2-D rings (z = 0) around the origin which are then shifted to their final positions.
Used for colliding rings testcase, see Gray, Monaghan, Swift SPH elastic dynamics, journal of Computer methods
in applied mechanics and engineering (2001).
"""

DIM = 2

if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Create an initial condition HDF5 file for the colliding rings test case.")
    parser.add_argument("--delta_p", "-d", metavar="float", type=float,
                        default=0.2, help="particle spacing (default: 0.2)")
    parser.add_argument("--r_inner", metavar="float", type=float,
                        default=3.0, help="inner ring radius (default: 3.0)")
    parser.add_argument("--r_outer", metavar="float", type=float,
                        default=4.0, help="outer ring radius (default: 4.0)")
    parser.add_argument("--shift", "-s", metavar="float", type=float,
                        default=5, help="x-axis shift of rings from origin (default: 5)")
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

    delta_p = args.delta_p
    r_inner = args.r_inner
    r_outer = args.r_outer
    shift = args.shift
    v_p = args.velocity
    density = args.density
    P = args.pressure
    gamma = args.gamma
    fillUp = args.fillUp

    mass = delta_p**2 * density
    u_const = P / ((gamma - 1.0) * density)

    print("Generating colliding rings initial conditions ...")
    print(f"  delta_p = {delta_p}, r_inner = {r_inner}, r_outer = {r_outer}")
    print(f"  shift = {shift}, v_p = {v_p}, density = {density}")
    print(f"  mass = {mass}, u = {u_const}")

    # create initial distribution through creating a 2d grid with dims 2*r_outer x 2*r_outer
    # then delete particles which are not on the ring
    N_length = int(2 * r_outer / delta_p)
    N_square = N_length**2

    # coordinates of particles in square
    r = np.zeros((N_square, DIM))

    # 2D meshgrid
    a = np.mgrid[0:N_length, 0:N_length]

    # create square
    for i in range(DIM):
        k = 0
        for j in range(N_square):
            if j % N_length == 0 and j > 0:
                k += 1
                k = k % N_length
            r[j, i] = (a[i, k, j % N_length] - (N_length - 1) / 2) * delta_p

    # count particles in one ring
    N = 0
    arr = np.zeros(N_square)
    for i in range(N_square):
        radius = np.sqrt(r[i, 0]**2 + r[i, 1]**2)
        if r_outer >= radius >= r_inner:
            N += 1
            arr[i] = 1

    # construct two rings with N particles which then are shifted along the x-axis
    if fillUp:
        r_ring = np.zeros((N, DIM + 1))
        r_ring2 = np.zeros((N, DIM + 1))
        v = np.zeros((N, DIM + 1))
        v2 = np.zeros((N, DIM + 1))
    else:
        r_ring = np.zeros((N, DIM))
        r_ring2 = np.zeros((N, DIM))
        v = np.zeros((N, DIM))
        v2 = np.zeros((N, DIM))

    m = np.ones(2 * N) * mass
    materialId = np.zeros(2 * N, dtype=np.int8)
    u = np.ones(2 * N) * u_const

    # create ring 1
    counter = 0
    for i in range(N_square):
        if arr[i] == 1:
            r_ring[counter, 0] = r[i, 0] - shift
            r_ring[counter, 1] = r[i, 1]
            if fillUp:
                r_ring[counter, 2] = 0
            v[counter, 0] = v_p
            counter += 1

    # create ring 2
    counter = 0
    for i in range(N_square):
        if arr[i] == 1:
            r_ring2[counter, 0] = r[i, 0] + shift
            r_ring2[counter, 1] = r[i, 1]
            if fillUp:
                r_ring2[counter, 2] = 0
            v2[counter, 0] = -v_p
            counter += 1

    # put two rings in one array
    r_final = np.concatenate((r_ring, r_ring2))
    v_final = np.concatenate((v, v2))

    suffix = "3D" if fillUp else "2D"
    filename = f"rings_deltap{delta_p}-{suffix}.h5"

    outH5 = h5.File(filename, "w")
    outH5.create_dataset("x", data=r_final)
    outH5.create_dataset("v", data=v_final)
    outH5.create_dataset("m", data=m)
    outH5.create_dataset("u", data=u)
    outH5.create_dataset("materialId", data=materialId)
    outH5.create_dataset("time", data=0.0 * m)
    outH5.close()

    print(f"Number of particles: {2 * N}")
    print(f"... done. Output written to {filename}")

    if args.plot:
        print("Plotting initial configuration ...")
        fig, ax = plt.subplots(figsize=(8, 6), dpi=200)
        # sc = ax.scatter(r_final[:, 0], r_final[:, 1], c=np.linalg.norm(m, axis=1), s=1.)
        sc = ax.scatter(r_final[:, 0], r_final[:, 1], c=m, s=1.)
        fig.colorbar(sc, ax=ax, label=r"m[i]")
        ax.set_title("Colliding rings initial conditions")
        ax.set_xlabel("$x$")
        ax.set_ylabel("$y$")
        ax.set_aspect("equal")
        plt.tight_layout()
        plotname = f"rings_deltap{delta_p}-{suffix}.png"
        plt.savefig(plotname)
        print(f"... saved to {plotname}")
