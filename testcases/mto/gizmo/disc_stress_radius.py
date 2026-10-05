#!/usr/bin/env python3
"""Pressure and polar stress components against radius for the GIZMO rotating-disc series (runs/disc_<rung>_f0.5).

    ./disc_stress_radius.py OUTDIR [--fps 1] [--quantities mag p xx yy xy] [--prefix 2026-10-05-gizmo-scheibe]

One movie per quantity, one panel per run (2 x 4 grid), one frame per snapshot (t = 0 ... T/2). Each particle is a
dot at (r, q), r = distance from the disc centre. Axes are common to all panels and frames of a movie.

Quantities (S = deviatoric StressTensor, in-plane 2x2 block; e_r, e_t = radial / tangential unit vector of the particle):
    mag  |S|, Frobenius norm of the full 3x3 tensor
    p    Pressure (snapshot field)
    xx   S_xx' = e_t.S.e_t  (tangential-tangential)
    yy   S_yy' = e_r.S.e_r  (radial-radial)
    xy   S_xy' = e_r.S.e_t  (polar shear S_r,theta; (e_t, e_r) is a left-handed pair, so this is the standard polar sign)
S_zz = 0 and tr S = 0 in the plane, hence S_yy' = -S_xx' exactly.
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

# key -> (file slug, y label, title text)
QUANTITIES = {
    "mag": ("stress-radius", "$|S|$", "Spannungsbetrag $|S|$"),
    "p": ("druck-radius", "$P$", "Druck $P$"),
    "xx": ("sxx-tangential-radius", "$S_{xx'}$ (tangential)", "$S_{xx'}$ (tangential)"),
    "yy": ("syy-radial-radius", "$S_{yy'}$ (radial)", "$S_{yy'}$ (radial)"),
    "xy": ("sxy-polar-radius", "$S_{xy'}$ ($=S_{r\\theta}$)", "$S_{xy'}$ ($=S_{r\\theta}$)"),
}


def load(run):
    out = []
    for s in sorted(glob.glob(os.path.join(HERE, "runs", run, "output", "snapshot_*.hdf5"))):
        with h5py.File(s, "r") as f:
            g = f["PartType0"]
            S = g["StressTensor"][:].reshape(-1, 3, 3)
            xy = g["Coordinates"][:][:, :2] - OFFSET
            phi = np.arctan2(xy[:, 1], xy[:, 0])
            er = np.stack([np.cos(phi), np.sin(phi)], 1)
            et = np.stack([-np.sin(phi), np.cos(phi)], 1)
            S2 = S[:, :2, :2]
            q = {"mag": np.linalg.norm(S.reshape(-1, 9), axis=1),
                 "p": g["Pressure"][:],
                 "xx": np.einsum("ni,nij,nj->n", et, S2, et),
                 "yy": np.einsum("ni,nij,nj->n", er, S2, er),
                 "xy": np.einsum("ni,nij,nj->n", er, S2, et)}
            out.append((float(f["Header"].attrs["Time"]), np.hypot(xy[:, 0], xy[:, 1]), q))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("outdir")
    ap.add_argument("--fps", type=float, default=1.0)
    ap.add_argument("--quantities", nargs="+", default=["p", "xx", "yy", "xy"], choices=list(QUANTITIES))
    ap.add_argument("--prefix", default="2026-10-05-gizmo-scheibe")
    a = ap.parse_args()

    runs = sorted((os.path.basename(d) for d in glob.glob(os.path.join(HERE, "runs", "disc_d0.0*_f0.5"))),
                  key=lambda r: -float(r.split("_")[1][1:]))
    data = {r: load(r) for r in runs}
    nfr = min(len(v) for v in data.values())
    rmax = 1.05 * max(f[1].max() for v in data.values() for f in v)

    for key in a.quantities:
        slug, ylab, ttext = QUANTITIES[key]
        lo = min(f[2][key].min() for v in data.values() for f in v)
        hi = max(f[2][key].max() for v in data.values() for f in v)
        pad = 0.05 * (hi - lo)
        fig, axs = plt.subplots(2, 4, figsize=(16, 7.6), dpi=110, sharex=True, sharey=True)
        pts = []
        for ax, r in zip(axs.flat, runs):
            t, rad, q = data[r][0]
            pts.append(ax.scatter(rad, q[key], s=1.5, c="#4c78a8", alpha=0.35, linewidths=0, rasterized=True))
            ax.set_xlim(0, rmax); ax.set_ylim(lo - pad, hi + pad); ax.grid(alpha=0.25)
            ax.axhline(0, color="k", lw=0.6)
            ax.set_title("d0=%s (N=%d)" % (r.split("_")[1][1:], len(rad)), fontsize=10)
        for ax in axs[1]:
            ax.set_xlabel("Radius $r$")
        for ax in axs[:, 0]:
            ax.set_ylabel(ylab)
        title = fig.suptitle("")
        fig.tight_layout(rect=(0, 0, 1, 0.96))

        def frame(k, key=key, ttext=ttext):
            for sc, r in zip(pts, runs):
                t, rad, q = data[r][k]
                sc.set_offsets(np.column_stack([rad, q[key]]))
            title.set_text("GIZMO Scheibe: %s gegen Radius   t = %.2f  (T/2 = 15.71)" % (ttext, t))

        out = os.path.join(a.outdir, "%s-%s.mp4" % (a.prefix, slug))
        ani = animation.FuncAnimation(fig, frame, frames=nfr)
        ani.save(out, writer=animation.FFMpegWriter(fps=a.fps, codec="libx264",
                                                    extra_args=["-pix_fmt", "yuv420p", "-vf", "pad=ceil(iw/2)*2:ceil(ih/2)*2"]))
        plt.close(fig)
        print("wrote", out, "(%d frames, %g fps, range %.3e .. %.3e)" % (nfr, a.fps, lo, hi))


if __name__ == "__main__":
    main()
