#!/usr/bin/env python3
"""Limiter experiment: analysis (runs from limiterTest.sh here and in ~/Documents/GIZMO/soundwave_2d).

Hypothesis: the error spikes at the wave extrema are caused by the slope/pairwise
limiters, and the demonstrator's ~2x larger spikes (vs GIZMO) come from its limiter.

Tillotson, N_NN = 32, n = 64 and 128, same ICs. Variants per code:
  ref       both limiters (reference runs)
  limSL     slope limiter only
  limPW     pairwise limiter only
  limNone   no limiter
  limSLfix  (demonstrator only) both limiters, slope limiter as Hopkins (2015) eq. B2

Metrics at t = 1 (shape residual as in spikeComparison.py: y-averaged perturbation
minus its own best-fit sine): max |r| within 2.5 dx of an extremum ("spike"), rms
elsewhere ("background"); plus L1_rho, amplitude a_rho(1), dispersion from convergence.py.

Writes plots/limiter_summary.txt, plots/limiter_spacetime.png, plots/limiter_spikes.png.
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from convergence import analyse, order
from spikeComparison import series, DIVERGING, SURFACE, INK, INK_2, MUTED, GRID, AXIS, style

VARIANTS = [("ref", "both limiters"), ("limSL", "slope only"), ("limPW", "pairwise only"),
            ("limNone", "no limiter"), ("limSLfix", "both, slope\nlimiter as eq. B2")]
CODES = [("demo", "demonstrator", "#2a78d6"), ("gizmo", "GIZMO", "#eb6834")]
NS = (64, 128)


def suffix(code, var):
    s = "" if var == "ref" else "_" + var
    return s + ("_gizmo" if code == "gizmo" else "")


def exists(code, var, n):
    return os.path.isdir("output_tillotson_N{}{}".format(n, suffix(code, var)))


if __name__ == "__main__":
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    data, res = {}, {}
    for code, _, _ in CODES:
        for var, _ in VARIANTS:
            for n in NS:
                if exists(code, var, n):
                    stem = "tillotson_N{}{}".format(n, suffix(code, var))
                    data[(code, var, n)] = series(stem, n)
                    res[(code, var, n)] = analyse("tillotson", n, 32, False, suffix(code, var))

    # ---- table -------------------------------------------------------------------
    lines = ["Limiter experiment: Tillotson soundwave, N_NN = 32, same ICs, t = 1", "",
             "{:8s} {:9s} {:>4s} {:>10s} {:>10s} {:>6s} {:>10s} {:>9s} {:>10s}".format(
                 "code", "variant", "n", "spike_rho", "bg_rho", "ratio", "L1_rho", "a_rho1", "dc_disp")]
    for (code, var, n), d in data.items():
        r = res[(code, var, n)]
        lines.append("{:8s} {:9s} {:4d} {:10.3e} {:10.3e} {:6.2f} {:10.3e} {:9.5f} {:+10.3e}".format(
            code, var, n, d["near_rho"][-1], d["far_rho"][-1], d["near_rho"][-1] / d["far_rho"][-1],
            r["L1_rho"], r["a_rho_t1"], r["dc_disp"]))
    lines += ["", "orders 64 -> 128:  spike_rho  bg_rho  L1_rho  |1-a_rho|  |dc_disp|"]
    for code, _, _ in CODES:
        for var, _ in VARIANTS:
            if all((code, var, n) in data for n in NS):
                d1, d2 = data[(code, var, 64)], data[(code, var, 128)]
                r1, r2 = res[(code, var, 64)], res[(code, var, 128)]
                lines.append("  {:6s} {:9s} {:9.2f} {:7.2f} {:7.2f} {:10.2f} {:10.2f}".format(
                    code, var, order(d1["near_rho"][-1], d2["near_rho"][-1]),
                    order(d1["far_rho"][-1], d2["far_rho"][-1]), order(r1["L1_rho"], r2["L1_rho"]),
                    order(abs(1 - r1["a_rho_t1"]), abs(1 - r2["a_rho_t1"])),
                    order(abs(r1["dc_disp"]), abs(r2["dc_disp"]))))
    lines += ["", "demonstrator / GIZMO, same variant:  spike_rho   bg_rho   L1_rho"]
    for var, _ in VARIANTS:
        for n in NS:
            if ("demo", var, n) in data and ("gizmo", var, n) in data:
                a, b = data[("demo", var, n)], data[("gizmo", var, n)]
                lines.append("  {:9s} n={:3d} {:18.2f} {:8.2f} {:8.2f}".format(
                    var, n, a["near_rho"][-1] / b["near_rho"][-1], a["far_rho"][-1] / b["far_rho"][-1],
                    res[("demo", var, n)]["L1_rho"] / res[("gizmo", var, n)]["L1_rho"]))
    text = "\n".join(lines)
    with open("plots/limiter_summary.txt", "w") as f:
        f.write(text + "\n")
    print(text)

    # ---- figure 1: space-time maps, n = 64, rows = code, columns = variant ----------
    n = 64
    fig, axes = plt.subplots(2, len(VARIANTS), figsize=(3.3 * len(VARIANTS), 7.4), sharex=True, sharey=True,
                             gridspec_kw={"wspace": 0.08, "hspace": 0.22})
    vmax = max(np.abs(d["r_rho"]).max() for (c, v, nn), d in data.items() if nn == n)
    mesh = None
    for row, (code, clabel, _) in enumerate(CODES):
        for col, (var, vlabel) in enumerate(VARIANTS):
            ax = axes[row, col]
            if (code, var, n) not in data:
                ax.axis("off")
                continue
            d = data[(code, var, n)]
            mesh = ax.pcolormesh(d["xc"], d["t"], d["r_rho"], cmap=DIVERGING, vmin=-vmax, vmax=vmax,
                                 shading="nearest", rasterized=True)
            for off in (0.25, 0.75):
                tt = np.linspace(0, 1, 400)
                xx = (tt + off) % 1.0
                brk = np.where(np.diff(xx) < 0)[0] + 1
                for seg_t, seg_x in zip(np.split(tt, brk), np.split(xx, brk)):
                    ax.plot(seg_x, seg_t, color=INK, lw=0.5, alpha=0.4)
            ax.set_title("{}, {}".format(clabel, vlabel.replace("\n", " ")), fontsize=9, color=INK, loc="left")
            ax.text(0.02, 0.02, "max |r| = {:.1e}".format(np.abs(d["r_rho"]).max()), transform=ax.transAxes,
                    fontsize=7.5, color=INK_2, bbox=dict(facecolor=SURFACE, edgecolor="none", pad=1.5, alpha=0.85))
            ax.set_xlim(0, 1)
            ax.set_ylim(0, 1)
            ax.tick_params(length=2.5, width=0.6)
            if col == 0:
                ax.set_ylabel("t")
            if row == 1 or (code, var) == ("demo", "limSLfix"):
                ax.set_xlabel("x")
                ax.tick_params(labelbottom=True)
    cb = fig.colorbar(mesh, ax=axes.ravel().tolist(), fraction=0.015, pad=0.012)
    cb.set_label("density shape residual [A], Tillotson, n = 64", color=INK_2)
    cb.outline.set_edgecolor(AXIS)
    fig.savefig("plots/limiter_spacetime.png", dpi=200, bbox_inches="tight")
    plt.close(fig)

    # ---- figure 2: spike and background vs variant --------------------------------
    fig, axes = plt.subplots(1, 2, figsize=(11, 3.8), sharey=True, gridspec_kw={"wspace": 0.08})
    xs = np.arange(len(VARIANTS))
    for ax, key, title in zip(axes, ("near_rho", "far_rho"),
                              ("spike: max |shape residual| near an extremum, t = 1",
                               "background: rms shape residual elsewhere, t = 1")):
        style(ax)
        for k, (code, _, c) in enumerate(CODES):
            for nn, filled in ((64, False), (128, True)):
                pts = [(i, data[(code, var, nn)][key][-1]) for i, (var, _) in enumerate(VARIANTS)
                       if (code, var, nn) in data]
                if not pts:
                    continue
                xi, yi = zip(*pts)
                ax.plot(np.array(xi) + (k - 0.5) * 0.22, yi, "o", color=c, ms=7, mew=1.5,
                        mfc=c if filled else SURFACE, mec=c, zorder=3)
        ax.set_yscale("log")
        ax.set_xticks(xs)
        ax.set_xticklabels([lab for _, lab in VARIANTS], fontsize=8)
        ax.set_xlim(-0.5, len(VARIANTS) - 0.5)
        ax.set_title(title, fontsize=9.5, color=INK, loc="left")
    axes[0].set_ylabel("[A]")
    handles = [Line2D([], [], ls="", marker="o", ms=7, color=c, label=lab) for _, lab, c in CODES]
    handles += [Line2D([], [], ls="", marker="o", ms=7, mfc=SURFACE, mec=INK_2, mew=1.5, label="n = 64"),
                Line2D([], [], ls="", marker="o", ms=7, color=INK_2, label="n = 128")]
    fig.legend(handles=handles, loc="upper center", ncol=4, frameon=False, fontsize=8.5,
               bbox_to_anchor=(0.5, 1.07), labelcolor=INK_2)
    fig.savefig("plots/limiter_spikes.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print("Wrote plots/limiter_summary.txt, plots/limiter_spacetime.png, plots/limiter_spikes.png")
