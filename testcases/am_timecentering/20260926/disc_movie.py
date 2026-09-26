#!/usr/bin/env python3
"""Movie of one body run (rings / blocks / disc): particles coloured by rho/rho0 (turbo, 0.9..1.1, colorbar),
disc additionally shows the Lagrangian angle tag (initial polar angle) so rotation is visible.
Usage: body_movie.py runs/T1_rings/ref_d0.05 movies/T1_rings_ref_d0.05.mp4 [--every 1]"""
import sys, os, glob, argparse, numpy as np, h5py, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.animation import FFMpegWriter
ap = argparse.ArgumentParser(); ap.add_argument("run"); ap.add_argument("out"); ap.add_argument("--every", type=int, default=1)
a = ap.parse_args()
fs = sorted(glob.glob(os.path.join(a.run, "output", "*.h5")))[::a.every]
with h5py.File(fs[0]) as h: x0 = h["x"][:]
disc = True
tag = np.arctan2(x0[:, 1], x0[:, 0])
# fixed view: 1.3 x the initial body radius about the initial centre, so a disrupting run
# does not zoom out until the body is invisible; particles outside are counted in the title
c0 = x0.mean(0); R0 = np.max(np.linalg.norm(x0 - c0, axis=1)); lo, hi = c0 - 1.3*R0, c0 + 1.3*R0; pad = 0.
from scipy.spatial import cKDTree
d0 = np.median(cKDTree(x0).query(x0, 2)[0][:, 1])
fig, axs = plt.subplots(1, 2 if disc else 1, figsize=(11 if disc else 8, 5), facecolor="#fcfcfb", squeeze=False)
# marker = one particle spacing on screen (points^2), so the body reads as a filled solid at every resolution
axw_in = fig.get_figwidth()*(0.38 if disc else 0.75); span = (hi-lo).max() + 2*pad
ms = max(0.5, (0.9*d0/span*axw_in*72)**2)
w = FFMpegWriter(fps=20)
with w.saving(fig, a.out, dpi=110):
    for f in fs:
        with h5py.File(f) as h:
            t = float(np.ravel(h["time"])[0]); x = h["x"][:]; rho = h["rho"][:]
        for ax in axs[0]: ax.clear(); ax.set_facecolor("#fcfcfb"); ax.set_aspect("equal")
        ax = axs[0][0]
        sc = ax.scatter(x[:, 0], x[:, 1], c=rho, s=ms, cmap="turbo", vmin=0.9, vmax=1.1, lw=0)
        if f == fs[0]: fig.colorbar(sc, ax=ax, shrink=.8, pad=.02, label="ρ/ρ0")
        ax.set_xlim(lo[0]-pad, hi[0]+pad); ax.set_ylim(lo[1]-pad, hi[1]+pad)
        nout = int(np.sum(np.any((x < lo) | (x > hi), axis=1)))
        ax.set_title("%s   t = %.2f   rho %.3f..%.3f%s" % (os.path.basename(a.run.rstrip("/")), t, rho.min(), rho.max(),
                     "   %d outside view" % nout if nout else ""), fontsize=9, loc="left")
        if disc:
            ax2 = axs[0][1]; ax2.scatter(x[:, 0], x[:, 1], c=tag, s=ms, cmap="twilight", lw=0)
            ax2.set_xlim(lo[0]-pad, hi[0]+pad); ax2.set_ylim(lo[1]-pad, hi[1]+pad); ax2.set_title("Lagrange-Winkelmarke", fontsize=9, loc="left")
        w.grab_frame()
print("wrote", a.out, len(fs), "frames")
