#!/usr/bin/env python3
"""Convergence figure of the thirdorder study: L1(rho)/A at t = 1 vs n. (a) today's scheme, Hopkins' recipe
(gradient order 3 alone) and the full third-order scheme; (b) the full scheme with one ingredient removed."""
import os, sys, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from analysis import metrics
HERE = os.path.dirname(os.path.abspath(__file__))
SURF, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e4e3df"
C = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948"]
def series(tag):
    ns, e = [], []
    for n in (32, 64, 128, 256, 512):
        d = os.path.join(HERE, "run_%s_N%d/output" % (tag, n))
        if not os.path.isdir(d): continue
        r = metrics(d, n, 1.0)
        if r is None: continue
        ns.append(n); e.append(r[0])
    return np.array(ns), np.array(e)
A = [("t3p", "full third-order scheme", C[0]), ("b1", "today's scheme (GEF 0, KDK, limiters)", C[1]),
     ("hop", "gradient order 3 only", C[2])]
B = [("a_pref", "- P_ref (unstable)", C[7]), ("a_f1", "- higher-order faces", C[3]), ("a_seed", "- IC density seed", C[4]),
     ("a_kdk", "- RK3 (TIME_KDK)", C[5]), ("a_plim", "+ pairwise limiter", C[6]), ("t3p", "full", C[0])]
plt.rcParams.update({"font.size": 10, "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2,
                     "ytick.color": INK2, "text.color": INK})
fig, axs = plt.subplots(1, 2, figsize=(11, 4.6), facecolor=SURF, sharey=True)
for ax, S, title in ((axs[0], A, "(a) scheme"), (axs[1], B, "(b) full scheme, one ingredient changed")):
    ax.set_facecolor(SURF)
    ax.grid(True, which="major", color=GRID, lw=0.8); ax.set_axisbelow(True)
    for s in ("top", "right"): ax.spines[s].set_visible(False)
    n0 = np.array([32, 512])
    for p, y0, lab in ((2, 2e-2, "p = 2"), (3, 1.4e-3, "p = 3")):
        ax.plot(n0, y0*(n0/32.)**(-p), ls=(0, (4, 3)), color=INK2, lw=1)
        ax.text(n0[1]*1.05, y0*(n0[1]/32.)**(-p), lab, color=INK2, va="center", fontsize=9)
    for tag, lab, col in S:
        ns, e = series(tag)
        if len(ns) == 0: continue
        ok = e < 1
        ax.plot(ns[ok], e[ok], "-o", color=col, lw=2, ms=6, mec=SURF, mew=1.5, label=lab)
        if (~ok).any():
            ax.annotate("blows up (L1 = %.0f)" % e[~ok][0], (ns[ok][-1], e[ok][-1]), xytext=(30, 40),
                        textcoords="offset points", color=INK2, fontsize=8.5,
                        arrowprops=dict(arrowstyle="-", color=INK2, lw=0.8))
    ax.set_xscale("log", base=2); ax.set_yscale("log")
    ax.set_xticks([32, 64, 128, 256, 512]); ax.set_xticklabels(["32", "64", "128", "256", "512"])
    ax.set_xlim(26, 800); ax.set_ylim(1e-7, 0.2)
    ax.set_xlabel("particles per wavelength n")
    ax.set_title(title, loc="left", color=INK, fontsize=11)
    ax.legend(frameon=False, fontsize=8.5, loc="lower left")
axs[0].set_ylabel("L1(rho) / A   at t = 1")
fig.suptitle("Ideal-gas soundwave, xSplitElasticHLLC, mu = 0, fixed CFL (dt = dx/3c)", x=0.01, ha="left",
             color=INK, fontsize=11)
fig.tight_layout()
out = os.path.join(HERE, "fig_convergence.png"); fig.savefig(out, dpi=150, facecolor=SURF); print(out)
