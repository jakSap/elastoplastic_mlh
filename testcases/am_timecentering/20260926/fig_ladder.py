#!/usr/bin/env python3
"""Figure 2: L_z error (max over one rotation) vs particle spacing d, rotating disc omega 0.2.
Left: square-lattice ICs (disc_gen.py), right: concentric-ring ICs (disc_rings_gen.py). Slope-1 and slope-2 guides."""
import os, re, glob, numpy as np, h5py, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
os.chdir(os.path.dirname(os.path.abspath(__file__)))
C = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]; INK, INK2, GRID, SURF = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
def err(run):
    fs = sorted(glob.glob(f"runs/{run}/output/*.h5")); L = []
    for f in fs:
        with h5py.File(f) as h:
            x = h["x"][()]; v = h["v"][()]; m = h["m"][()]; L.append(np.sum(m*(x[:, 0]*v[:, 1]-x[:, 1]*v[:, 0])))
            t = h["time"][()][0]
    L = np.array(L); return np.abs((L-L[0])/L[0]).max(), t
def ladder(prefix, v):
    pts = []
    for p in glob.glob(f"runs/L_{prefix}{v}_d*"):
        if re.sub(r"_d[0-9.]+$", "", os.path.basename(p))[2:] != prefix + v: continue
        if not os.path.exists(p+"/exit_code") or open(p+"/exit_code").read().strip() != "0": continue
        e, t = err(os.path.basename(p))
        if t >= 30.99: pts.append((float(re.search(r"_d([0-9.]+)$", p).group(1)), e))
    return np.array(sorted(pts))
# series: (colour slot, label, square-lattice variant, concentric-ring variant); fixed colour per scheme
V = [(0, "Standard (DRIFT_KDK 1, RECON_FACE_DV 1)", "std", "stdx"),
     (1, "TIME_KDK + DV 2 + Minmod + Fläche 3", None, "gkdk_fm3_dv2_mmr"),
     (2, "TIME_DKD + räuml. Sprung + DV 2 + Fläche 3 (Quadrat: + Minmod)", "dkd_sp_fm3_dv2_mmr", "dkd_sp_fm3_dv2"),
     (3, "TIME_DKD + räuml. Sprung + Fläche 3, TC aus (DV 1)", "dkd_sp_fm3_tc0", "dkd_sp_fm3_tc0")]
plt.rcParams.update({"font.size": 9, "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2,
                     "axes.spines.top": False, "axes.spines.right": False})
fig, ax = plt.subplots(1, 2, figsize=(11, 4.9), facecolor=SURF, sharey=True)
handles = {}
for j, (a, pre, title) in enumerate(zip(ax, ("", "R"), ("Quadratgitter-Ausschnitt", "konzentrische Ringe"))):
    a.set_facecolor(SURF); a.grid(True, which="both", color=GRID, lw=0.5); a.set_xscale("log"); a.set_yscale("log")
    for k, lab, vq, vr in V:
        v = vq if pre == "" else vr
        if v is None: continue
        P = ladder(pre, v)
        if len(P):
            h, = a.plot(P[:, 0], P[:, 1], "o-", color=C[k], lw=2, ms=6, mec=SURF, mew=1.5, label=lab); handles[k] = h
    d = np.array([0.016, 0.11])
    for s_, y0, lab in ((1, 3e-3, "p = 1"), (2, 6e-5, "p = 2")):
        a.plot(d, y0*(d/0.025)**s_, "--", color=INK2, lw=1); a.text(d[-1], y0*(d[-1]/0.025)**s_, " "+lab, color=INK2, va="center", fontsize=8)
    a.set_xlabel("Teilchenabstand d"); a.set_title(title + ", CFL 0.2 (dt ∝ d), 1 Umdrehung", loc="left", color=INK)
ax[0].set_ylabel(r"max$_t\,|\Delta L_z/L_z|$")
fig.legend([handles[k] for k in sorted(handles)], [V[k][1] for k in sorted(handles)], frameon=False, fontsize=8,
           loc="lower center", ncol=2)
fig.tight_layout(rect=(0, 0.1, 1, 1))
fig.savefig("figs/am_ladder.png", dpi=130, facecolor=SURF); print("figs/am_ladder.png")
