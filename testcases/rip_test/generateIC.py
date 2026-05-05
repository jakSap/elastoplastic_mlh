#!/usr/bin/env python3

import argparse
import numpy as np
import matplotlib.pyplot as plt
import h5py as h5

"""
Creates an IC for the rip test: two horizontal bars of the same length stacked
above each other separated by a gap.  The upper bar has width w; the lower bar
has width w + 2h with symmetric V-shaped cuts of depth h on both its top and
bottom faces, centred in x.  At the centre the two cuts together remove 2h,
leaving exactly thickness w.  An optional rigid-body velocity and/or velocity
gradient in y can be applied to pull the bars apart.
"""

DIM = 2


def make_rectangle(x_min, x_max, y_min, y_max, delta_p):
    """Return (x, y) arrays for particles on a regular grid within a rectangle."""
    xs_1d = np.arange(x_min, x_max + 0.5 * delta_p, delta_p)
    ys_1d = np.arange(y_min, y_max + 0.5 * delta_p, delta_p)
    xx, yy = np.meshgrid(xs_1d, ys_1d)
    return xx.ravel(), yy.ravel()


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Create an IC HDF5 file for the rip test: a narrow bar (width w) "
                    "above a wide bar (width w+h) with a V-shaped centre cut of depth h.")
    parser.add_argument("--delta_p", "-d", metavar="float", type=float,
                        default=0.1,
                        help="particle spacing (default: 0.1)")
    parser.add_argument("--length", "-l", metavar="float", type=float,
                        default=10.0,
                        help="bar length in x, same for both bars (default: 10.0)")
    parser.add_argument("--width", "-w", metavar="float", type=float,
                        default=1.0,
                        help="narrow bar width in y (default: 1.0)")
    parser.add_argument("--cut_depth", "-H", metavar="float", type=float,
                        default=1.0,
                        help="V-cut depth h; wide bar has width w+h, "
                             "leaving thickness w at the cut tip (default: 1.0)")
    parser.add_argument("--cut_width", "-W", metavar="float", type=float,
                        default=2.0,
                        help="full width cw of the V-cut in x (default: 2.0)")
    parser.add_argument("--gap", "-g", metavar="float", type=float,
                        default=1.0,
                        help="gap between the two bars in y (default: 1.0)")
    parser.add_argument("--velocity", "-v", metavar="float", type=float,
                        default=0.0,
                        help="rigid-body separation velocity: upper bar gets +v, "
                             "lower bar gets -v in y (default: 0.0)")
    parser.add_argument("--gradient", metavar="float", type=float,
                        default=0.0,
                        help="initial velocity gradient in y: v_y = gradient * y "
                             "applied to all particles; positive pulls bars apart "
                             "(default: 0.0)")
    parser.add_argument("--density", metavar="float", type=float,
                        default=1.0, help="particle density (default: 1.0)")
    parser.add_argument("--pressure", "-P", metavar="float", type=float,
                        default=2.5, help="pressure constant (default: 2.5)")
    parser.add_argument("--gamma", metavar="float", type=float,
                        default=5.0 / 3.0, help="adiabatic index (default: 5/3)")
    parser.add_argument("--fillUp", "-f", action="store_true",
                        help="fill up coordinates to 3D with z=0")
    parser.add_argument("--plot", action="store_true",
                        help="plot initial configuration")

    args = parser.parse_args()

    delta_p     = args.delta_p
    L           = args.length
    w           = args.width
    h           = args.cut_depth   # also the extra width of the wide bar
    cw          = args.cut_width
    gap         = args.gap
    v_rigid     = args.velocity
    strain_rate = args.gradient
    density     = args.density
    P           = args.pressure
    gamma       = args.gamma
    fillUp      = args.fillUp

    mass    = delta_p**2 * density
    u_const = P / ((gamma - 1.0) * density)

    print("Generating rip test initial conditions ...")
    print(f"  delta_p = {delta_p}")
    print(f"  Narrow bar:  length={L},  width={w}")
    print(f"  Wide bar:    length={L},  width={w + h}  (= w + h)")
    print(f"  V-cut:       depth={h},   width={cw}  →  tip thickness = {w}")
    print(f"  Gap = {gap}")
    print(f"  v_rigid = {v_rigid},  strain_rate = {strain_rate}")
    print(f"  mass = {mass:.4g},  u = {u_const:.4g}")

    # ---- Narrow bar (upper): y in [gap/2, gap/2 + w], x in [-L/2, L/2] ----
    xs1, ys1 = make_rectangle(-L / 2, L / 2, gap / 2, gap / 2 + w, delta_p)

    # ---- Wide bar (lower): y in [-(gap/2 + w + h), -gap/2], x in [-L/2, L/2] ----
    y_top_wide = -gap / 2
    y_bot_wide = -(gap / 2 + w + h)
    xs2, ys2 = make_rectangle(-L / 2, L / 2, y_bot_wide, y_top_wide, delta_p)

    # V-shaped cut on the top face of the wide bar, centred at x=0.
    # Depth tapers linearly from h at x=0 to 0 at x=±cw/2.
    cut_depth_at_x = np.where(np.abs(xs2) <= cw / 2,
                              h * (1.0 - np.abs(xs2) / (cw / 2)),
                              0.0)
    cut_mask = ys2 >= y_top_wide - cut_depth_at_x
    xs2 = xs2[~cut_mask]
    ys2 = ys2[~cut_mask]

    N1      = len(xs1)
    N2      = len(xs2)
    N_total = N1 + N2

    dim_out = DIM + 1 if fillUp else DIM

    r1 = np.zeros((N1, dim_out))
    r1[:, 0] = xs1
    r1[:, 1] = ys1

    r2 = np.zeros((N2, dim_out))
    r2[:, 0] = xs2
    r2[:, 1] = ys2

    v1 = np.zeros((N1, dim_out))
    v2 = np.zeros((N2, dim_out))
    v1[:, 1] = +v_rigid + strain_rate * ys1
    v2[:, 1] = -v_rigid + strain_rate * ys2

    r_final = np.concatenate((r1, r2))
    v_final = np.concatenate((v1, v2))

    m          = np.ones(N_total) * mass
    materialId = np.zeros(N_total, dtype=np.int8)
    u          = np.ones(N_total) * u_const

    suffix   = "3D" if fillUp else "2D"
    filename = f"rip_deltap{delta_p}-{suffix}.h5"

    outH5 = h5.File(filename, "w")
    outH5.create_dataset("x",          data=r_final)
    outH5.create_dataset("v",          data=v_final)
    outH5.create_dataset("m",          data=m)
    outH5.create_dataset("u",          data=u)
    outH5.create_dataset("materialId", data=materialId)
    outH5.create_dataset("time",       data=0.0 * m)
    outH5.close()

    print(f"  N_narrow = {N1},  N_wide = {N2},  N_total = {N_total}")
    print(f"... done. Output written to {filename}")

    if args.plot:
        print("Plotting initial configuration ...")
        fig, ax = plt.subplots(figsize=(12, 5), dpi=200)
        ax.scatter(xs1, ys1, s=1.5, linewidths=0, color="tab:blue",   label="narrow bar")
        ax.scatter(xs2, ys2, s=1.5, linewidths=0, color="tab:orange", label="wide bar")
        ax.set_aspect("equal")
        ax.set_title("Rip test initial conditions")
        ax.set_xlabel("$x$")
        ax.set_ylabel("$y$")
        ax.legend(markerscale=5)
        plt.tight_layout()
        plotname = f"rip_deltap{delta_p}-{suffix}.png"
        plt.savefig(plotname)
        print(f"... saved to {plotname}")
