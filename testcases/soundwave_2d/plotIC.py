#!/usr/bin/env python3
"""Velocity plots of the soundwave initial conditions.

One column per resolution (64 x 16, 128 x 32, 256 x 64 on [0,1) x [0,1/4)):
  top    - v_x over the periodic box (diverging map, zero = neutral gray)
  bottom - v_x(x) of every particle against the analytic c A sin(k x)

The EOS only changes u, not v (c = 1 for all three), so the velocity is taken
from --eos and checked bit-for-bit against the other two ICs.

Usage: ./plotIC.py [--eos ideal] [-o plots/ic_velocity.png]
"""

import argparse
import os
import numpy as np
import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

RESOLUTIONS = (64, 128, 256)
EOS_NAMES = ("ideal", "murnaghan", "tillotson")

SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_2 = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"
SERIES_1 = "#2a78d6"
# diverging blue <-> red, neutral gray midpoint, equal steps per arm
DIVERGING = LinearSegmentedColormap.from_list(
    "bluered", ["#184f95", "#6da7ec", "#f0efec", "#ef8e8d", "#b52f2f"])

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.size": 9,
    "text.color": INK,
    "axes.labelcolor": INK_2,
    "axes.edgecolor": AXIS,
    "axes.linewidth": 0.8,
    "xtick.color": MUTED,
    "ytick.color": MUTED,
    "xtick.labelcolor": INK_2,
    "ytick.labelcolor": INK_2,
    "figure.facecolor": SURFACE,
    "axes.facecolor": SURFACE,
    "savefig.facecolor": SURFACE,
})


def load(eos, n):
    with h5py.File("soundwave_{}_N{}.h5".format(eos, n), "r") as f:
        return f["x"][:], f["v"][:], dict(f.attrs)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--eos", choices=EOS_NAMES, default="ideal")
    ap.add_argument("-o", "--output", default="plots/ic_velocity.png")
    args = ap.parse_args()
    os.chdir(os.path.dirname(os.path.abspath(__file__)))

    fig, axes = plt.subplots(2, 3, figsize=(11, 4.6),
                             gridspec_kw={"height_ratios": [0.42, 1], "hspace": 0.45,
                                          "wspace": 0.28})
    scale = 1e6  # plot in units of 1e-6

    for col, n in enumerate(RESOLUTIONS):
        x, v, attrs = load(args.eos, n)
        for other in EOS_NAMES:
            _, v_o, _ = load(other, n)
            if not np.array_equal(v, v_o):
                raise SystemExit("v differs between {} and {} at N{}".format(args.eos, other, n))
        if np.any(v[:, 1] != 0.0):
            raise SystemExit("v_y != 0 in N{}".format(n))

        A, c, lam = attrs["amplitude"], attrs["c_bar"], attrs["wavelength"]
        vmax = c * A * scale

        # IC lattice is written in meshgrid "ij" order: reshape to (nx, ny)
        ny, box_y = int(attrs["ny"]), float(attrs["box_y"])
        X = x[:, 0].reshape(n, ny)
        Y = x[:, 1].reshape(n, ny)
        VX = v[:, 0].reshape(n, ny) * scale

        ax = axes[0, col]
        mesh = ax.pcolormesh(X, Y, VX, cmap=DIVERGING, vmin=-vmax, vmax=vmax,
                             shading="nearest", rasterized=True)
        ax.set_aspect("equal")
        ax.set_xlim(0, 1)
        ax.set_ylim(0, box_y)
        ax.set_yticks([0, box_y])
        ax.set_title(r"${} \times {}$".format(n, ny), color=INK, fontsize=10, pad=6)
        ax.set_xlabel("x")
        if col == 0:
            ax.set_ylabel("y")
        ax.tick_params(length=2.5, width=0.6)

        ax = axes[1, col]
        xa = np.linspace(0, 1, 1000)
        ax.axhline(0, color=AXIS, lw=0.8, zorder=1)
        ax.plot(xa, c * A * np.sin(2 * np.pi * xa / lam) * scale, color=INK_2, lw=1.0,
                zorder=2, label=r"analytic $cA\sin(kx)$")
        cols = slice(None, None, n)  # one particle per lattice column; all rows equal
        ax.plot(x[cols, 0], v[cols, 0] * scale, "o", ms=3.2 if n == 64 else 2.2,
                color=SERIES_1, mec=SURFACE, mew=0.4, zorder=3,
                label="particles (all {} rows identical)".format(ny))
        ax.set_xlim(0, 1)
        ax.set_ylim(-1.25 * vmax, 1.25 * vmax)
        ax.grid(True, color=GRID, lw=0.6)
        ax.set_axisbelow(True)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        ax.set_xlabel("x")
        if col == 0:
            ax.set_ylabel(r"$v_x\ [10^{-6}]$")
            ax.legend(loc="lower left", frameon=False, fontsize=8, labelcolor=INK_2)
        ax.tick_params(length=2.5, width=0.6)

    cbar = fig.colorbar(mesh, ax=axes[0, :].tolist(), fraction=0.025, pad=0.02)
    cbar.set_label(r"$v_x\ [10^{-6}]$", color=INK_2)
    cbar.outline.set_edgecolor(AXIS)
    cbar.ax.tick_params(length=2.5, width=0.6)

    fig.suptitle("Soundwave IC velocity: identical for ideal gas, Murnaghan and Tillotson "
                 r"($A = 10^{-6}$, $c = 1$, $\lambda = 1$, $v_y \equiv 0$)",
                 color=INK, fontsize=11, x=0.06, ha="left", y=0.985)

    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
    fig.savefig(args.output, dpi=200, bbox_inches="tight")
    print("Wrote", args.output)
