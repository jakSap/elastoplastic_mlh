#!/usr/bin/env python3

import argparse
import os
import sys

import numpy as np
import matplotlib.pyplot as plt
import h5py as h5

# shared Grady-Kipp Weibull flaw generator (testcases/weibull_flaws.py)
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from weibull_flaws import write_flaws

"""
Creates an IC for the rip test: two horizontal bars of the same length stacked
above each other separated by a gap.  The upper bar has width w; the lower bar
has width w + 2h with symmetric V-shaped cuts of depth h on both its top and
bottom faces, centred in x.  At the centre the two cuts together remove 2h,
leaving exactly thickness w.  An initial x-velocity profile splits each bar
into three equal parts: the inner third is at rest, while the outer thirds ramp
linearly from zero at the inner boundary to ±v_max at the bar ends (left end
moves left, right end moves right).
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
                        help="max pulling velocity: the outer thirds of each bar "
                             "ramp linearly from 0 at x=±L/6 to ±v at x=±L/2 "
                             "(left end −v, right end +v) (default: 0.0)")
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
    v_max       = args.velocity
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
    print(f"  v_max = {v_max}  (outer-third ramp in x)")
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

    # x-velocity profile: three equal segments of length L/3.
    # Inner third  |x| < L/6          : v_x = 0
    # Left outer   x in [-L/2, -L/6]  : v_x ramps from 0 to -v_max
    # Right outer  x in [ L/6,  L/2]  : v_x ramps from 0 to +v_max
    x_inner = L / 6.0

    def vx_profile(xs):
        v = np.zeros(len(xs))
        left  = xs < -x_inner
        right = xs >  x_inner
        v[left]  = v_max * (xs[left]  + x_inner) * 3.0 / L
        v[right] = v_max * (xs[right] - x_inner) * 3.0 / L
        return v

    v1 = np.zeros((N1, dim_out))
    v2 = np.zeros((N2, dim_out))
    v1[:, 0] = vx_profile(xs1)
    v2[:, 0] = vx_profile(xs2)

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
    # Grady-Kipp Weibull flaws (read by the solver only when FRAGMENTATION=1;
    # harmless otherwise). max_flaws MUST match MAX_NUM_FLAWS in parameter.h.
    write_flaws(outH5, m / density, k=1.0e4, m=8.0, max_flaws=1, seed=0)
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
