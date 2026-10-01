#!/usr/bin/env python3
"""Second figure of the thirdorder study: face reconstruction on top of P_ref + IC seed + RK3 + GO3 weights.
L1(rho)/A at t = 1 (A = 1e-7 ladders; qT at A = 1e-6) and the damping 1 - a (signed: negative = amplitude grows)."""
import os, sys, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from analysis import metrics
HERE = os.path.dirname(os.path.abspath(__file__))
SURF, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e4e3df"
C = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948"]
def series(tag, k=0):
    ns, e = [], []
    for n in (32, 64, 128, 256, 512):
        d = os.path.join(HERE, "run_%s_N%d/output" % (tag, n))
        if not os.path.isdir(d): continue
        r = metrics(d, n, 1.0)
        if r is None: continue
        ns.append(n); e.append(r[k])
    return np.array(ns), np.array(e)
S = [("t3p", "cubic faces, linear recon, RK3", C[0], "-"),
     ("qT", "cubic faces, Taylor quadratic recon, RK3", C[3], "-"),
     ("q23r4", "cubic faces, quadratic recon x 2/3, RK4", C[1], "-"),
     ("f1q23r4", "first-order faces, quadratic recon x 2/3, RK4", C[2], "-")]
plt.rcParams.update({"font.size": 10, "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2,
                     "ytick.color": INK2, "text.color": INK})
fig, axs = plt.subplots(1, 2, figsize=(11, 4.6), facecolor=SURF)
for ax, k, title, ylab in ((axs[0], 6, "(a) L1(rho) / A vs weakly nonlinear reference", "L1(rho) / A   at t = 1"),
                           (axs[1], 2, "(b) damping 1 - a (filled: > 0, open: < 0, amplitude grows)", "|1 - a|")):
    ax.set_facecolor(SURF); ax.grid(True, color=GRID, lw=0.8); ax.set_axisbelow(True)
    for s in ("top", "right"): ax.spines[s].set_visible(False)
    guides = ((3, 2e-3, 512), (4, 6e-5, 512)) if k == 6 else ((3, 2e-3, 512), (5, 3e-5, 256))
    for p, y0, nend in guides:
        n0 = np.array([32, nend])
        ax.plot(n0, y0*(n0/32.)**(-p), ls=(0, (4, 3)), color=INK2, lw=1)
        ax.text(n0[1]*1.05, y0*(n0[1]/32.)**(-p), "p = %d" % p, color=INK2, va="center", fontsize=9)
    for tag, lab, col, ls in S:
        ns, e = series(tag, k)
        if len(ns) == 0: continue
        ax.plot(ns, np.abs(e), ls, color=col, lw=2, label=lab)
        pos = e > 0
        ax.plot(ns[pos], np.abs(e[pos]), "o", color=col, ms=6, mec=SURF, mew=1.5)
        ax.plot(ns[~pos], np.abs(e[~pos]), "o", color=SURF, ms=6, mec=col, mew=1.8)
    ax.set_xscale("log", base=2); ax.set_yscale("log")
    ax.set_xticks([32, 64, 128, 256, 512]); ax.set_xticklabels(["32", "64", "128", "256", "512"])
    ax.set_xlim(26, 800); ax.set_ylim(3e-10, 1e-2)
    ax.set_xlabel("particles per wavelength n"); ax.set_ylabel(ylab)
    ax.set_title(title, loc="left", color=INK, fontsize=11)
axs[0].legend(frameon=False, fontsize=8, loc="lower left")
fig.suptitle("Face reconstruction and time integration (all: P_ref, exact IC density, GO 3, CFL 1/3, A = 1e-6)",
             x=0.01, ha="left", color=INK, fontsize=11)
fig.tight_layout()
out = os.path.join(HERE, "fig_recon.png"); fig.savefig(out, dpi=150, facecolor=SURF); print(out)
