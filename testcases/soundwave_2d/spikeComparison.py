#!/usr/bin/env python3
"""Demonstrator vs GIZMO on the soundwave: is the error spike at the wave extrema
(where dv/dx changes sign) a property of the method?

Runs compared (all on the same uncorrected ICs, N_NN = 32, Wendland C2):
  demonstrator ideal gas  output_ideal_N<n>
  GIZMO ideal gas         output_ideal_N<n>_gizmo  (see ~/Documents/GIZMO/soundwave_2d)
  demonstrator Tillotson  output_tillotson_N<n>    (barotropic reference)
  GIZMO Tillotson         output_tillotson_N<n>_gizmo  (table and plots/spike_spacetime_tillotson.png)

Shape residual: the y-averaged perturbation (density: (<rho>_y - <rho>)/(<rho> A),
velocity: <v_x>_y / A) minus its own best-fit sine a sin(k(x - t) - phi). This
removes the smooth amplitude and phase error, so what remains is the non-sinusoidal
part of the error, i.e. the spikes. "Near an extremum": within 2.5 dx of a zero of
cos(k(x - t) - phi).

Writes plots/spike_profiles_t1.png, plots/spike_spacetime.png, plots/spike_spacetime_tillotson.png,
plots/spike_vs_t.png, plots/gizmo_comparison.txt.
"""

import os
import numpy as np
import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D

from makeMovie import snapshots, sine_fit
from convergence import analyse, order

A, K = 1e-6, 2.0 * np.pi
SURFACE, INK, INK_2, MUTED, GRID, AXIS = "#fcfcfb", "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7"
SERIES = [("ideal", "", "demonstrator, ideal gas", "#2a78d6"),
          ("ideal", "_gizmo", "GIZMO, ideal gas", "#eb6834"),
          ("tillotson", "", "demonstrator, Tillotson", "#1baf7a")]
GIZMO_TILLOTSON = ("tillotson", "_gizmo", "GIZMO, Tillotson", "#eb6834")
ALL_SERIES = SERIES + [GIZMO_TILLOTSON]                      # table
TILLOTSON_SERIES = [SERIES[2], GIZMO_TILLOTSON]               # plots/spike_spacetime_tillotson.png
DIVERGING = LinearSegmentedColormap.from_list("bluered", ["#184f95", "#6da7ec", "#f0efec", "#ef8e8d", "#b52f2f"])

plt.rcParams.update({
    "font.family": "sans-serif", "font.size": 9, "text.color": INK,
    "axes.labelcolor": INK_2, "axes.edgecolor": AXIS, "axes.linewidth": 0.8,
    "xtick.color": MUTED, "ytick.color": MUTED, "xtick.labelcolor": INK_2, "ytick.labelcolor": INK_2,
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
})


def profiles(path, n):
    with h5py.File(path, "r") as d:
        x, rho, vx = d["x"][:, 0], d["rho"][:], d["v"][:, 0]
    col = np.floor(x * n).astype(int) % n
    cnt = np.bincount(col, minlength=n)
    prho = np.bincount(col, weights=rho, minlength=n) / cnt
    pvx = np.bincount(col, weights=vx, minlength=n) / cnt
    xc = (np.arange(n) + 0.5) / n
    return xc, (prho - rho.mean()) / (rho.mean() * A), pvx / A


def shape_residual(xc, y, t):
    a, phi = sine_fit(xc, y, K, t)
    r = y - a * np.sin(K * (xc - t) - phi)
    s = (xc - t - phi / K - 0.25) % 0.5          # distance to nearest extremum of the fitted wave
    dist = np.minimum(s, 0.5 - s)
    return r, dist, a, phi


def series(stem, n):
    files, times = snapshots("output_" + stem)
    out = {"t": np.array(times), "r_rho": [], "r_vx": [], "near_rho": [], "far_rho": [],
           "near_vx": [], "far_vx": [], "phi": []}
    for f, t in zip(files, times):
        xc, yrho, yvx = profiles(f, n)
        for key, y in (("rho", yrho), ("vx", yvx)):
            r, dist, a, phi = shape_residual(xc, y, t)
            near = dist <= 2.5 / n
            out["r_" + key].append(r)
            out["near_" + key].append(np.abs(r[near]).max())
            out["far_" + key].append(np.sqrt(np.mean(r[~near] ** 2)))
            if key == "rho":
                out["phi"].append(phi)
    out["xc"] = xc
    for key in list(out):
        if key not in ("t", "xc"):
            out[key] = np.array(out[key])
    return out


def style(ax):
    ax.grid(True, color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.tick_params(length=2.5, width=0.6)


def spacetime(data, panels, out, n=64):
    """Space-time maps of the density shape residual, one panel per run, shared colour scale."""
    fig, axes = plt.subplots(1, len(panels), figsize=(4 * len(panels), 4.2), sharey=True,
                             gridspec_kw={"wspace": 0.08})
    vmax = max(np.abs(data[(e, s, n)]["r_rho"]).max() for e, s, _, _ in panels)
    for ax, (eos, sfx, label, _) in zip(axes, panels):
        d = data[(eos, sfx, n)]
        mesh = ax.pcolormesh(d["xc"], d["t"], d["r_rho"], cmap=DIVERGING, vmin=-vmax, vmax=vmax,
                             shading="nearest", rasterized=True)
        for off in (0.25, 0.75):                     # extremum paths x = t + off (mod 1)
            tt = np.linspace(0, 1, 400)
            xx = (tt + off) % 1.0
            brk = np.where(np.diff(xx) < 0)[0] + 1
            for seg_t, seg_x in zip(np.split(tt, brk), np.split(xx, brk)):
                ax.plot(seg_x, seg_t, color=INK, lw=0.6, alpha=0.5)
        ax.set_title("{}, n = {}".format(label, n), fontsize=10, color=INK, loc="left")
        ax.set_xlabel("x")
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.tick_params(length=2.5, width=0.6)
    axes[0].set_ylabel("t")
    cb = fig.colorbar(mesh, ax=axes.tolist(), fraction=0.02 * 3 / len(panels), pad=0.015 * 3 / len(panels))
    cb.set_label("density shape residual [A]", color=INK_2)
    cb.outline.set_edgecolor(AXIS)
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    data = {(eos, sfx, n): series("{}_N{}{}".format(eos, n, sfx), n)
            for eos, sfx, _, _ in ALL_SERIES for n in (64, 128)}

    # ---- table -----------------------------------------------------------------
    lines = ["Soundwave, demonstrator vs GIZMO (same uncorrected ICs, N_NN = 32, Wendland C2), t = 1", ""]
    lines.append("{:24s} {:>10s} {:>10s} {:>9s} {:>10s} {:>10s} {:>10s} {:>9s}".format(
        "run", "<rho>-1", "L1_rho", "a_rho1", "phi1", "dc_disp", "L1_vx", "dE/E"))
    res = {}
    for eos, sfx, _, _ in ALL_SERIES:
        for n in (64, 128):
            r = analyse(eos, n, 32, False, sfx)
            res[(eos, sfx, n)] = r
            lines.append("{:24s} {:+10.3e} {:10.3e} {:9.5f} {:+10.3e} {:+10.3e} {:10.3e} {:+9.1e}".format(
                r["run"], r["rho_mean_m1_t1"], r["L1_rho"], r["a_rho_t1"], r["phi_rho_t1"], r["dc_disp"],
                r["L1_vx"], r["dE"]))
    lines.append("")
    lines.append("orders 64 -> 128:   L1_rho  |1-a_rho|  |dc_disp|   spike_rho(t=1)   spike_vx(t=1)")
    for eos, sfx, label, _ in ALL_SERIES:
        r1, r2 = res[(eos, sfx, 64)], res[(eos, sfx, 128)]
        d1, d2 = data[(eos, sfx, 64)], data[(eos, sfx, 128)]
        lines.append("  {:26s} {:6.2f} {:9.2f} {:11.2f} {:16.2f} {:15.2f}".format(
            label, order(r1["L1_rho"], r2["L1_rho"]), order(abs(1 - r1["a_rho_t1"]), abs(1 - r2["a_rho_t1"])),
            order(abs(r1["dc_disp"]), abs(r2["dc_disp"])),
            order(d1["near_rho"][-1], d2["near_rho"][-1]), order(d1["near_vx"][-1], d2["near_vx"][-1])))
    lines.append("")
    lines.append("shape residual [A] at t = 1: max |r| within 2.5 dx of an extremum / rms elsewhere / ratio")
    for eos, sfx, label, _ in ALL_SERIES:
        for n in (64, 128):
            d = data[(eos, sfx, n)]
            lines.append("  {:26s} n={:3d}  rho {:.2e} / {:.2e} / {:4.1f}   vx {:.2e} / {:.2e} / {:4.1f}".format(
                label, n, d["near_rho"][-1], d["far_rho"][-1], d["near_rho"][-1] / d["far_rho"][-1],
                d["near_vx"][-1], d["far_vx"][-1], d["near_vx"][-1] / d["far_vx"][-1]))
    text = "\n".join(lines)
    with open("plots/gizmo_comparison.txt", "w") as f:
        f.write(text + "\n")
    print(text)

    # ---- figure 1: shape residual at t = 1 ---------------------------------------
    fig, axes = plt.subplots(2, 2, figsize=(11, 5.6), sharex=True, gridspec_kw={"hspace": 0.15, "wspace": 0.18})
    for col, n in enumerate((64, 128)):
        for row, key in enumerate(("rho", "vx")):
            ax = axes[row, col]
            style(ax)
            ax.axhline(0, color=AXIS, lw=0.8, zorder=1)
            phi = data[("ideal", "_gizmo", n)]["phi"][-1]
            for xe in (0.25, 0.75):
                ax.axvline((xe + 1.0 + phi / K) % 1.0, color=MUTED, lw=0.8, zorder=1)
            for eos, sfx, label, c in SERIES:
                d = data[(eos, sfx, n)]
                ax.plot(d["xc"], d["r_" + key][-1], "-o", color=c, lw=1.5, ms=3 if n == 64 else 0,
                        mec=SURFACE, mew=0.5, zorder=3, label=label)
            ax.set_xlim(0, 1)
            if row == 0:
                ax.set_title("n = {} per wavelength, t = 1".format(n), fontsize=10, color=INK, loc="left")
            else:
                ax.set_xlabel("x")
            if col == 0:
                ax.set_ylabel((r"density" if key == "rho" else r"velocity $v_x$") + "\nshape residual [A]")
    axes[0, 0].legend(loc="upper left", frameon=False, fontsize=8, labelcolor=INK_2)
    axes[1, 1].text(0.25, axes[1, 1].get_ylim()[1], "  extremum of the wave\n  (dv/dx changes sign)",
                    color=MUTED, fontsize=8, va="top")
    fig.savefig("plots/spike_profiles_t1.png", dpi=200, bbox_inches="tight")
    plt.close(fig)

    # ---- figure 2: space-time maps of the density shape residual, n = 64 ---------
    for panels, out in ((SERIES, "plots/spike_spacetime.png"),
                        (TILLOTSON_SERIES, "plots/spike_spacetime_tillotson.png")):
        spacetime(data, panels, out)
    print("Wrote plots/spike_spacetime.png, plots/spike_spacetime_tillotson.png")

    # ---- figure 3: spike size vs time --------------------------------------------
    fig, axes = plt.subplots(1, 2, figsize=(11, 3.6), sharey=True, gridspec_kw={"wspace": 0.08})
    for ax, key, title in zip(axes, ("rho", "vx"), ("density", "velocity $v_x$")):
        style(ax)
        for eos, sfx, label, c in SERIES:
            for n, lw, mk in ((64, 1.2, dict(marker="o", markevery=10, ms=5, mfc=SURFACE, mec=c)), (128, 2.0, {})):
                d = data[(eos, sfx, n)]
                ax.plot(d["t"], d["near_" + key], color=c, lw=lw, **mk)
        ax.set_yscale("log")
        ax.set_xlim(0, 1)
        ax.set_xlabel("t")
        ax.set_title("max shape residual near an extremum, " + title, fontsize=9.5, color=INK, loc="left")
    axes[0].set_ylabel("[A]")
    handles = [Line2D([], [], color=c, lw=2, label=label) for _, _, label, c in SERIES]
    handles += [Line2D([], [], color=INK_2, lw=1.2, marker="o", ms=5, mfc=SURFACE, mec=INK_2, label="n = 64"),
                Line2D([], [], color=INK_2, lw=2.0, label="n = 128")]
    fig.legend(handles=handles, loc="upper center", ncol=5, frameon=False, fontsize=8.5,
               bbox_to_anchor=(0.5, 1.08), labelcolor=INK_2)
    fig.savefig("plots/spike_vs_t.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print("Wrote plots/spike_profiles_t1.png, plots/spike_spacetime.png, plots/spike_vs_t.png")
