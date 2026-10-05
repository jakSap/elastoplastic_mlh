#!/usr/bin/env python3
"""Particle renders of the GIZMO rotating-disc series (runs/disc_<rung>_f0.5): density and |S| per particle.

    ./disc_renders.py OUTDIR [--snap -1] [--prefix 2026-10-05-gizmo-scheibe-render]
    ./disc_renders.py OUTDIR --movie [--prefix 2026-10-05-gizmo-scheibe-film]

One PNG per run at the chosen snapshot (default: last, t = T/2). With --movie: one mp4 per run over all snapshots.
Colour scales are common to all runs (1st..99.9th percentile of density, 0..max of |S|); for the movies
over all runs and all snapshots.
"""
import argparse, glob, os
import h5py
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import animation

HERE = os.path.dirname(os.path.abspath(__file__))
OFFSET = 5.0


def load(run, k):
    snaps = sorted(glob.glob(os.path.join(HERE, "runs", run, "output", "snapshot_*.hdf5")))
    with h5py.File(snaps[k], "r") as f:
        g = f["PartType0"]
        xy = g["Coordinates"][:][:, :2] - OFFSET
        return (float(f["Header"].attrs["Time"]), xy, g["Density"][:],
                np.linalg.norm(g["StressTensor"][:].reshape(-1, 9), axis=1))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("outdir")
    ap.add_argument("--snap", type=int, default=-1)
    ap.add_argument("--prefix", default="2026-10-05-gizmo-scheibe-render")
    ap.add_argument("--movie", action="store_true")
    a = ap.parse_args()
    if a.movie and a.prefix.endswith("render"):
        a.prefix = a.prefix[:-len("render")] + "film"

    runs = sorted((os.path.basename(d) for d in glob.glob(os.path.join(HERE, "runs", "disc_d0.0*_f0.5"))),
                  key=lambda r: -float(r.split("_")[1][1:]))
    lim, panel_in = 1.15, 4.6            # axes half-width, panel width [in]
    pt_per_unit = panel_in * 72 / (2 * lim)

    if a.movie:
        nsnap = {r: len(glob.glob(os.path.join(HERE, "runs", r, "output", "snapshot_*.hdf5"))) for r in runs}
        frames = {r: [load(r, k) for k in range(nsnap[r])] for r in runs}
        allf = [f for r in runs for f in frames[r]]
        rlo, rhi = np.percentile(np.concatenate([f[2] for f in allf]), [1, 99.9])
        smax = max(f[3].max() for f in allf)
        for r in runs:
            d0s = r.split("_")[1][1:]
            s = (float(d0s) * pt_per_unit * 1.15) ** 2
            fig, axs = plt.subplots(1, 2, figsize=(2 * panel_in + 2.2, panel_in + 0.9), dpi=120)
            t, xy, rho, S = frames[r][0]
            scs = []
            for ax, val, cmap, vmin, vmax, lab in ((axs[0], rho, "viridis", rlo, rhi, "Dichte $\\rho$"),
                                                   (axs[1], S, "magma", 0, smax, "$|S|$ (Frobenius-Norm)")):
                sc = ax.scatter(xy[:, 0], xy[:, 1], c=val, s=s, cmap=cmap, vmin=vmin, vmax=vmax, linewidths=0)
                ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.set_aspect("equal")
                ax.set_xlabel("$x$"); ax.set_ylabel("$y$")
                fig.colorbar(sc, ax=ax, shrink=0.8, label=lab)
                scs.append(sc)
            title = fig.suptitle("")
            fig.tight_layout()

            def frame(k, r=r, d0s=d0s):
                t, xy, rho, S = frames[r][k]
                for sc, val in zip(scs, (rho, S)):
                    sc.set_offsets(xy)
                    sc.set_array(val)
                title.set_text("GIZMO Scheibe d0=%s (N=%d)   t = %.2f  (T/2 = 15.71)" % (d0s, len(rho), t))

            out = os.path.join(a.outdir, "%s-d%s.mp4" % (a.prefix, d0s.replace(".", "")))
            ani = animation.FuncAnimation(fig, frame, frames=nsnap[r])
            ani.save(out, writer=animation.FFMpegWriter(fps=2, codec="libx264", extra_args=["-pix_fmt", "yuv420p", "-vf", "pad=ceil(iw/2)*2:ceil(ih/2)*2"]))
            plt.close(fig)
            print("wrote", out, "(%d frames)" % nsnap[r])
        return

    data = {r: load(r, a.snap) for r in runs}
    rho_all = np.concatenate([d[2] for d in data.values()])
    rlo, rhi = np.percentile(rho_all, [1, 99.9])
    smax = max(d[3].max() for d in data.values())
    for r in runs:
        d0s = r.split("_")[1][1:]
        t, xy, rho, S = data[r]
        s = (float(d0s) * pt_per_unit * 1.15) ** 2
        fig, axs = plt.subplots(1, 2, figsize=(2 * panel_in + 2.2, panel_in + 0.9), dpi=150)
        for ax, val, cmap, vmin, vmax, lab in ((axs[0], rho, "viridis", rlo, rhi, "Dichte $\\rho$"),
                                               (axs[1], S, "magma", 0, smax, "$|S|$ (Frobenius-Norm)")):
            o = np.argsort(val)
            sc = ax.scatter(xy[o, 0], xy[o, 1], c=val[o], s=s, cmap=cmap, vmin=vmin, vmax=vmax, linewidths=0)
            ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.set_aspect("equal")
            ax.set_xlabel("$x$"); ax.set_ylabel("$y$")
            fig.colorbar(sc, ax=ax, shrink=0.8, label=lab)
        fig.suptitle("GIZMO Scheibe d0=%s (N=%d)   t = %.2f" % (d0s, len(rho), t))
        fig.tight_layout()
        out = os.path.join(a.outdir, "%s-d%s.png" % (a.prefix, d0s.replace(".", "")))
        fig.savefig(out); plt.close(fig)
        print("wrote", out, " rho %.4f..%.4f  |S|max %.3e" % (rho.min(), rho.max(), S.max()))


if __name__ == "__main__":
    main()
