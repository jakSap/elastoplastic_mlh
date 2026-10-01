#!/usr/bin/env python3
"""Rotating-disc metrics (thirdorder study, 2026-10-01) per run and snapshot:
  dLz    L_z(t)/L_z(0) - 1            (the A.1 metric; sign kept: < 0 = spin-down, > 0 = injection)
  dp     |P(t) - P(0)| / (M omega R)   (linear momentum, machine-exact in MFM -> anything else is a bug)
  dE     E(t)/E(0) - 1, E = sum m (u + v^2/2)
  rmax   max particle radius           (flag: > 1.2 R means material is leaving the body)
  vrig   <|v - omega_eff x r|>/(omega R), omega_eff = L_z/I   (breathing + lattice noise, rigid part removed)
  umin   min u                         (negative internal energy = instability signal; IC has u = 1)
  dw     differential rotation (<omega_loc>_{r<0.35} - <omega_loc>_{0.75<r<0.9})/omega_0, omega_loc = (x v_y - y v_x)/r^2
         (core spin-up against the outer ring; 0 for rigid rotation, surface rows excluded). The disc breathes
         (S = 0 start), which modulates dw by +-0.1 from snapshot to snapshot: the table gives the mean over the
         last rotation (32 snapshots) and the plot the running one-rotation mean.
Usage: disc/disc_analysis.py d tag [tag ...]   -> table at t_end + disc/fig_disc_d<d>.png"""
import os, sys, glob, re, numpy as np, h5py
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__))
SURF, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e4e3df"
C = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948"]
LABEL = {"dstd": "standard (GEF 1, KDK)", "dxs": "GEF 0, 2nd order (KDK)", "d3": "third order, cubic faces (unstable)",
         "d3fb": "third order, faded faces (FACE_GRAD_WEIGHTS 2)", "d3fbx0": "third order, faded faces, no XSPLIT_SHEAR_STAR",
         "d4fb": "fourth order (RK4 + quad recon), faded faces", "d4f1": "fourth order, first-order faces",
         "d3f1": "third order, first-order faces"}
COLOR = {"dstd": C[1], "dxs": C[3], "d3": C[6], "d3fb": C[0], "d3fbx0": C[4], "d4fb": C[2], "d4f1": C[5], "d3f1": C[7]}


def series(run):
    rows = []
    for fn in sorted(glob.glob(os.path.join(run, "output", "*.h5"))):
        with h5py.File(fn) as h:
            t = float(np.ravel(h["time"])[0]); x = h["x"][:]; v = h["v"][:]; m = h["m"][:]; u = h["u"][:]
        L = np.sum(m*(x[:, 0]*v[:, 1] - x[:, 1]*v[:, 0])); I = np.sum(m*(x**2).sum(1))
        om = L/I
        vr = v - om*np.column_stack([-x[:, 1], x[:, 0]])
        r = np.sqrt((x**2).sum(1)); wl = (x[:, 0]*v[:, 1] - x[:, 1]*v[:, 0])/np.maximum(r**2, 1e-12)
        core, rim = (r > 1e-9) & (r < 0.35), (r > 0.75) & (r < 0.9)
        rows.append((t, L, np.sum(m[:, None]*v, 0), np.sum(m*(u + .5*(v**2).sum(1))), r.max(),
                     np.mean(np.sqrt((vr**2).sum(1))), u.min(), m.sum(), (wl[core].mean() - wl[rim].mean())/0.2))
    return rows


def errors(run):
    try:
        log = open(os.path.join(run, "run.log"), errors="replace").read()
        ec = open(os.path.join(run, "exit_code")).read().strip()
    except Exception:
        return "running", 0
    return ec, log.count("ERROR")


if __name__ == "__main__":
    d = sys.argv[1]; tags = sys.argv[2:]
    om0, R = 0.2, 1.0
    out = ["d = %s; values at the last snapshot" % d,
           "%-6s %6s %7s %10s %9s %10s %6s %9s %7s %5s %8s %6s" % ("tag", "exit", "ERRORs", "t_end", "dLz", "dp", "rmax",
                                                                  "dE", "vrig", "umin", "dw", "wall")]
    fig, axs = plt.subplots(1, 3, figsize=(15, 5.6), facecolor=SURF)
    for ax in axs:
        ax.set_facecolor(SURF); ax.grid(True, color=GRID, lw=0.8); ax.set_axisbelow(True)
        for s in ("top", "right"): ax.spines[s].set_visible(False)
    for tag in tags:
        run = os.path.join(HERE, "run_%s_d%s" % (tag, d))
        if not os.path.isdir(run): continue
        S = series(run)
        if not S: continue
        t = np.array([r[0] for r in S]); L = np.array([r[1] for r in S]); P = np.array([r[2] for r in S])
        E = np.array([r[3] for r in S]); M = S[0][7]
        dLz = L/L[0] - 1; dp = np.abs(P - P[0]).max(1)/(M*om0*R); dE = E/E[0] - 1
        vrig = np.array([r[5] for r in S])/(om0*R)
        ec, ne = errors(run)
        try:
            wall = float(re.findall(r"Total wall-clock time: ([0-9.]+) min", open(os.path.join(run, "run.log")).read())[-1])
        except Exception:
            wall = np.nan
        dw = np.array([r[8] for r in S])
        nw = min(32, len(dw)); dwm = np.convolve(dw, np.ones(nw)/nw, mode="valid")
        out.append("%-6s %6s %7d %10.3f %+9.2e %10.2e %6.3f %+9.1e %7.1e %5.2f %+8.4f %6.1f" % (
            tag, ec, ne, t[-1], dLz[-1], dp.max(), S[-1][4], dE[-1], vrig[-1], min(r[6] for r in S), dwm[-1], wall))
        col = COLOR.get(tag, C[6]); lab = LABEL.get(tag, tag)
        trot = t/(2*np.pi/om0)
        axs[0].plot(trot, dLz, color=col, lw=2, label=lab)
        axs[1].plot(trot, np.abs(dE) + 1e-17, color=col, lw=2, label=lab)
        axs[2].plot(trot, [r[6] for r in S], color=col, lw=2, label=lab)
    axs[0].set_ylabel("L_z(t)/L_z(0) - 1  (symlog)"); axs[0].set_yscale("symlog", linthresh=1e-5)
    axs[0].set_title("(a) angular momentum", loc="left", fontsize=11)
    axs[1].set_ylabel("|E(t)/E(0) - 1|"); axs[1].set_yscale("log"); axs[1].set_title("(b) total energy", loc="left", fontsize=11)
    axs[2].set_ylabel("min u  (IC: u = 1)"); axs[2].set_ylim(0.6, 1.02)
    axs[2].set_title("(c) internal-energy drain at the surface", loc="left", fontsize=11)
    for ax in axs: ax.set_xlabel("rotations")
    h, l = axs[0].get_legend_handles_labels()
    fig.legend(h, l, frameon=False, fontsize=8.5, loc="lower center", ncol=3)
    fig.suptitle("Rotating disc (test plan A.1), d = %s, omega = 0.2, Tillotson mu = 0.22, fixed dt" % d,
                 x=0.01, ha="left", fontsize=11)
    fig.tight_layout(rect=(0, 0.13, 1, 1))
    fig.savefig(os.path.join(HERE, "fig_disc_d%s.png" % d), dpi=140, facecolor=SURF)
    print("\n".join(out))
    open(os.path.join(HERE, "disc_analysis_d%s.txt" % d), "w").write("\n".join(out) + "\n")
