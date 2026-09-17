#!/usr/bin/env python3
"""Murnaghan soundwave: does GIZMO reproduce the demonstrator's missing convergence?

Demonstrator: uncorrected Murnaghan (n = 2, K0 = rho0 = 1) does not converge in L1
(orders -0.17, -0.13): the kernel-sum density offset <rho> - 1 = +4.24e-3 gives
c = sqrt(<rho>) = 1.00212, i.e. a resolution-independent phase error of
2 pi (c - 1) t = 1.33e-2 rad per period; with corrected masses it converges at order 2.

GIZMO has no Murnaghan EOS; runs use its Tillotson with a = b = 0, A = K0, B = K0/2,
which is identical for n = 2 (see ~/Documents/GIZMO/soundwave_2d/makeRuns.py).
With corrected masses P crosses zero, so GIZMO's tensile correction acts
(solids/elastic_physics.c: f = 0.2 (W(r)/W(dp))^4, ~0.2 for nearest lattice
neighbours); one extra GIZMO set has it switched off (build_tillotson_notc, f = 0).

Runs (N_NN = 32, n = 64/128/256, same ICs in both codes):
  demonstrator  output_{murnaghan,tillotson}_N<n>, output_murnaghan_N<n>_corr
  GIZMO         the same with suffix _gizmo; output_murnaghan_N<n>_corr_notc_gizmo

Writes plots/murnaghan_gizmo.txt and plots/murnaghan_gizmo.png.
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from convergence import analyse, order, style_axes, EOS_COLOR, SURFACE, INK, INK_2, MUTED

NS = (64, 128, 256)
A = 1e-6
MURN, TILL, NOTC = EOS_COLOR["murnaghan"], EOS_COLOR["tillotson"], "#7a5ab8"
# (code, suffix, eos, corr, label, colour, marker, gizmo-style)
SERIES = [("demonstrator", "", "murnaghan", False, "Murnaghan, uncorrected", MURN, "o", False),
          ("demonstrator", "", "murnaghan", True, "Murnaghan, masses corrected", MURN, "s", False),
          ("demonstrator", "", "tillotson", False, "Tillotson", TILL, "o", False),
          ("GIZMO", "_gizmo", "murnaghan", False, "Murnaghan, uncorrected", MURN, "o", True),
          ("GIZMO", "_gizmo", "murnaghan", True, "Murnaghan, masses corrected", MURN, "s", True),
          ("GIZMO", "_notc_gizmo", "murnaghan", True, "Murnaghan, corrected, tensile corr. off", NOTC, "^", True),
          ("GIZMO", "_gizmo", "tillotson", False, "Tillotson", TILL, "o", True)]


def key(s):
    return s[0] + ", " + s[4]


if __name__ == "__main__":
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    res = {(key(s), n): analyse(s[2], n, 32, s[3], s[1]) for s in SERIES for n in NS}

    lines = ["Murnaghan soundwave (n = 2, K0 = rho0 = 1), N_NN = 32, t = 1: demonstrator vs GIZMO", "",
             "{:32s} {:>4s} {:>10s} {:>10s} {:>10s} {:>9s} {:>10s} {:>10s} {:>10s} {:>10s}".format(
                 "run", "n", "<rho>-1", "L1_nom", "L1_rho", "a_rho1", "phi1", "c_num-1", "c_exp-1", "dc_disp")]
    for s in SERIES:
        for n in NS:
            r = res[(key(s), n)]
            lines.append("{:32s} {:4d} {:+10.3e} {:10.3e} {:10.3e} {:9.5f} {:+10.3e} {:+10.3e} {:+10.3e} {:+10.3e}".format(
                r["run"], n, r["rho_mean_m1_t1"], r["L1_rho_nom"], r["L1_rho"], r["a_rho_t1"],
                r["phi_rho_t1"], r["c_num_m1"], r["c_exp_m1"], r["dc_disp"]))
    lines += ["", "orders [64->128, 128->256]:                  L1_rho        |1-a_rho|       |phi1|       |c_num-1|"]
    for s in SERIES:
        parts = []
        for f in (lambda r: r["L1_rho"], lambda r: abs(1 - r["a_rho_t1"]),
                  lambda r: abs(r["phi_rho_t1"]), lambda r: abs(r["c_num_m1"])):
            e = [f(res[(key(s), n)]) for n in NS]
            parts.append("{:5.2f}, {:5.2f}".format(order(e[0], e[1]), order(e[1], e[2])))
        lines.append("  {:52s} ".format(key(s)) + "   ".join(parts))
    lines += ["", "phase from the offset alone, 2 pi (c_exp - 1) t, vs phi(Murnaghan uncorr.) - phi(Tillotson), t = 1:"]
    for code in ("demonstrator", "GIZMO"):
        pred = [2 * np.pi * res[(code + ", Murnaghan, uncorrected", n)]["c_exp_m1"] for n in NS]
        meas = [res[(code + ", Murnaghan, uncorrected", n)]["phi_rho_t1"] - res[(code + ", Tillotson", n)]["phi_rho_t1"]
                for n in NS]
        lines.append("  {:12s} predicted {}   measured {}".format(
            code, " / ".join("{:.3e}".format(p) for p in pred), " / ".join("{:.3e}".format(m) for m in meas)))
    lines += ["", "error budget, uncorrected Murnaghan, t = 1: L1 of a pure phase error (2/pi) A |phi|,"
              " of a pure amplitude error (2/pi) A |1-a|, measured L1_rho"]
    for code in ("demonstrator", "GIZMO"):
        for n in NS:
            r = res[(code + ", Murnaghan, uncorrected", n)]
            lines.append("  {:12s} n={:3d}  phase {:.2e}   amplitude {:.2e}   measured {:.2e}".format(
                code, n, 2 / np.pi * A * abs(r["phi_rho_t1"]), 2 / np.pi * A * abs(1 - r["a_rho_t1"]), r["L1_rho"]))
    text = "\n".join(lines)
    with open("plots/murnaghan_gizmo.txt", "w") as f:
        f.write(text + "\n")
    print(text)

    # ---- figure ----------------------------------------------------------------------
    ns = np.array(NS)
    panels = [("L1", r"$L_1(\rho)$, offset removed", "log"),
              ("phi", r"phase $\varphi(t{=}1)$ [rad], + = ahead", "symlog"),
              ("amp", r"amplitude error $|1 - a_\rho(t{=}1)|$", "log")]
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.0), gridspec_kw={"wspace": 0.32})
    for ax, (pk, title, yscale) in zip(axes, panels):
        style_axes(ax, logx=True, logy=(yscale == "log"))
        if yscale == "symlog":
            ax.set_yscale("symlog", linthresh=1e-3)
        for s in SERIES:
            rs = [res[(key(s), n)] for n in NS]
            y = {"L1": [r["L1_rho"] for r in rs], "phi": [r["phi_rho_t1"] for r in rs],
                 "amp": [abs(1 - r["a_rho_t1"]) for r in rs]}[pk]
            gz, c = s[7], s[5]
            ax.plot(ns, y, "--" if gz else "-", marker=s[6], color=c, lw=1.3 if gz else 1.9, ms=7,
                    mfc=SURFACE if gz else c, mec=c, mew=1.4, zorder=3)
        if pk == "phi":
            ax.axhline(0, color=MUTED, lw=0.8, zorder=1)
            off = 2 * np.pi * res[("demonstrator, Murnaghan, uncorrected", 64)]["c_exp_m1"]
            ax.axhline(off, color=MUTED, lw=0.9, ls=":", zorder=1)
            ax.text(56, off * 1.15, r"offset: $2\pi(\sqrt{\langle\rho\rangle}-1)$", color=MUTED, fontsize=8, va="bottom")
        if pk == "L1":
            ref = res[("demonstrator, Tillotson", 128)]["L1_rho"] * 1.6
            ax.plot(ns, ref * (ns / 128.0) ** -2, ":", color=MUTED, lw=1.2, zorder=2)
            ax.text(256, ref * 0.25, r"  $\propto n^{-2}$", color=MUTED, fontsize=8, va="center")
        ax.set_xticks(ns)
        ax.set_xticklabels([str(n) for n in ns])
        ax.set_xlim(52, 330)
        ax.set_xlabel("particles per wavelength n")
        ax.set_title(title, fontsize=9.5, color=INK, loc="left")
    handles = [Line2D([], [], color=MURN, marker="o", ms=7, lw=1.8, label="Murnaghan, uncorrected"),
               Line2D([], [], color=MURN, marker="s", ms=7, lw=1.8, label="Murnaghan, masses corrected"),
               Line2D([], [], color=NOTC, marker="^", ms=7, lw=1.3, ls="--", mfc=SURFACE,
                      label="GIZMO Murnaghan corrected, tensile corr. off"),
               Line2D([], [], color=TILL, marker="o", ms=7, lw=1.8, label="Tillotson"),
               Line2D([], [], color=INK_2, lw=1.9, ls="-", marker="o", ms=7, label="demonstrator"),
               Line2D([], [], color=INK_2, lw=1.3, ls="--", marker="o", ms=7, mfc=SURFACE, mec=INK_2, label="GIZMO")]
    fig.legend(handles=handles, loc="upper center", ncol=6, frameon=False, fontsize=8.5,
               bbox_to_anchor=(0.5, 1.08), labelcolor=INK_2)
    fig.savefig("plots/murnaghan_gizmo.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print("Wrote plots/murnaghan_gizmo.txt, plots/murnaghan_gizmo.png")
