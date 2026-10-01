#!/usr/bin/env python3
"""Side-by-side movie of rotating-disc variants (same dt and dump cadence, so snapshot k is the same time in every
run). Per column: top = Lagrangian material tag (initial polar angle, the "pizza slices" of discMovie.py; a density
plot of an axisymmetric disc shows nothing), bottom = local spin (x v_y - y v_x)/r^2 / omega_0 - 1 (diverging,
gray = rigid rotation at the initial rate). Column title: L_z(t)/L_z(0) - 1.
Usage: disc/disc_movie.py d tag [tag ...] [-j 6]  -> disc/movies/disc_d<d>.mp4"""
import os, sys, glob, subprocess, shutil, numpy as np, h5py
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from concurrent.futures import ProcessPoolExecutor
HERE = os.path.dirname(os.path.abspath(__file__))
SURF, INK2 = "#fcfcfb", "#52514e"
LABEL = {"dstd": "standard (GEF 1)", "dxs": "GEF 0, 2nd order", "d3": "third order (cubic faces)",
         "d3x0": "3rd, no shear star", "d4": "fourth order (cubic faces)", "d3fb": "third order (faded faces)",
         "d4fb": "fourth order (faded faces)", "d3fbx0": "3rd, faded, no shear star", "d4f1": "4th, first-order faces",
         "d3f1": "3rd, first-order faces"}
OM0, RV = 0.2, 1.3


def frame(job):
    k, files, tags, theta0, L0, out, ms = job
    fig, axs = plt.subplots(2, len(tags), figsize=(3.3*len(tags), 6.9), facecolor=SURF, squeeze=False)
    t = None
    for c, (fn, tag) in enumerate(zip(files, tags)):
        for ax in axs[:, c]:
            ax.set_facecolor(SURF); ax.set_aspect("equal"); ax.set_xlim(-RV, RV); ax.set_ylim(-RV, RV)
            ax.set_xticks([]); ax.set_yticks([])
        if fn is None:
            axs[0, c].set_title("%s\n(no snapshot)" % LABEL.get(tag, tag), fontsize=9); continue
        with h5py.File(fn) as h:
            x = h["x"][:]; v = h["v"][:]; m = h["m"][:]; t = float(np.ravel(h["time"])[0])
        L = np.sum(m*(x[:, 0]*v[:, 1] - x[:, 1]*v[:, 0]))
        r2 = (x**2).sum(1)
        spin = np.where(r2 > 1e-12, (x[:, 0]*v[:, 1] - x[:, 1]*v[:, 0])/np.maximum(r2, 1e-12), np.nan)/OM0 - 1
        axs[0, c].scatter(x[:, 0], x[:, 1], c=theta0[c], s=ms, cmap="hsv", vmin=-np.pi, vmax=np.pi, linewidths=0)
        sc = axs[1, c].scatter(x[:, 0], x[:, 1], c=spin, s=ms, cmap="RdBu_r", vmin=-0.1, vmax=0.1, linewidths=0)
        axs[0, c].set_title("%s\nL_z/L_z0 - 1 = %+.2e" % (LABEL.get(tag, tag), L/L0[c] - 1), fontsize=9)
    fig.subplots_adjust(left=0.02, right=0.9, top=0.86, bottom=0.03, wspace=0.08, hspace=0.12)
    cb = fig.colorbar(sc, cax=fig.add_axes([0.92, 0.06, 0.012, 0.36]))
    cb.set_label("local spin / omega_0 - 1")
    fig.suptitle("rotating disc, t = %.2f (%.3f rotations); top: material tag, bottom: local spin (oscillates with the breathing of the S = 0 disc)" % (
        t, t*OM0/(2*np.pi)), x=0.01, ha="left", fontsize=10)
    fig.savefig(os.path.join(out, "%05d.png" % k), dpi=90, facecolor=SURF)
    plt.close(fig)


if __name__ == "__main__":
    args = [a for a in sys.argv[1:]]
    nj = 6
    if "-j" in args:
        i = args.index("-j"); nj = int(args[i + 1]); del args[i:i + 2]
    d, tags = args[0], args[1:]
    runs = [os.path.join(HERE, "run_%s_d%s" % (t, d), "output") for t in tags]
    files = [sorted(glob.glob(os.path.join(r, "*.h5"))) for r in runs]
    nk = max(len(f) for f in files)
    theta0, L0 = [], []
    for f in files:
        with h5py.File(f[0]) as h:
            x = h["x"][:]; v = h["v"][:]; m = h["m"][:]
        theta0.append(np.arctan2(x[:, 1], x[:, 0])); L0.append(np.sum(m*(x[:, 0]*v[:, 1] - x[:, 1]*v[:, 0])))
    out = os.path.join(HERE, "frames_d%s" % d); shutil.rmtree(out, ignore_errors=True); os.makedirs(out)
    ms = 6.0 if float(d) >= 0.05 else 1.6
    jobs = [(k, [f[k] if k < len(f) else None for f in files], tags, theta0, L0, out, ms) for k in range(nk)]
    with ProcessPoolExecutor(nj) as ex: list(ex.map(frame, jobs))
    os.makedirs(os.path.join(HERE, "movies"), exist_ok=True)
    mp4 = os.path.join(HERE, "movies", "disc_d%s.mp4" % d)
    subprocess.run(["ffmpeg", "-y", "-loglevel", "error", "-framerate", "12", "-pattern_type", "glob", "-i",
                    os.path.join(out, "[0-9]*.png"), "-vf", "pad=ceil(iw/2)*2:ceil(ih/2)*2", "-c:v", "libx264",
                    "-pix_fmt", "yuv420p", mp4], check=True)
    print(mp4, nk, "frames (kept in", out, "for the frame check)")
