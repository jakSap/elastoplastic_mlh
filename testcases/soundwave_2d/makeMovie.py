#!/usr/bin/env python3
"""Movie of the y-averaged density rho(x, t) of one soundwave run.

Each frame, from one snapshot:
  top    - y-averaged density perturbation of the working density `rho`
           (the explicitly integrated one), per lattice column,
           (<rho>_y(x) - <rho>) / (<rho> A), against the analytic wave
           sin(k (x - c t)) with the nominal c = 1
  bottom - residual numerical - analytic, same units

<rho> is the domain mean at that time. It is subtracted because the kernel-sum
density carries a constant offset (+4.2e-3 at N_NN = 32) that is 4000x the
wave amplitude; its value is printed in the frame.

    ./makeMovie.py ideal_N64            # -> movies/ideal_N64.mp4
"""

import argparse
import glob
import os
import numpy as np
import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.animation import FFMpegWriter

SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_2 = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"
SERIES_1 = "#2a78d6"

plt.rcParams.update({
    "font.family": "sans-serif", "font.size": 9, "text.color": INK,
    "axes.labelcolor": INK_2, "axes.edgecolor": AXIS, "axes.linewidth": 0.8,
    "xtick.color": MUTED, "ytick.color": MUTED,
    "xtick.labelcolor": INK_2, "ytick.labelcolor": INK_2,
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
})


def parse_stem(stem):
    """<eos>_N<n>[_nnn<N_NN>][_corr] -> (eos, n, N_NN)."""
    parts = stem.split("_")
    eos, n, nnn = parts[0], int(parts[1][1:]), 32
    for part in parts[2:]:
        if part.startswith("nnn"):
            nnn = int(part[3:])
    return eos, n, nnn


def snapshots(outdir):
    files = glob.glob(os.path.join(outdir, "*.h5"))
    times = []
    for f in files:
        with h5py.File(f, "r") as d:
            times.append(float(np.ravel(d["time"])[0]))
    order = np.argsort(times)
    return [files[i] for i in order], [times[i] for i in order]


def column_profile(path, n):
    """Return column centres, y-averaged rho per column and the domain mean."""
    with h5py.File(path, "r") as d:
        x = d["x"][:, 0]
        rho = d["rho"][:]
    col = np.floor(x * n).astype(int) % n
    counts = np.bincount(col, minlength=n)
    prof = np.bincount(col, weights=rho, minlength=n) / np.maximum(counts, 1)
    return (np.arange(n) + 0.5) / n, prof, rho.mean()


def sine_fit(xc, y, k, t):
    """Amplitude and phase shift of y ~ a sin(k (x - t) - phi).

    phi > 0: the numerical wave is AHEAD of the c = 1 wave; a wave with speed c
    gives phi = k (c - 1) t.
    """
    s, c = np.sin(k * (xc - t)), np.cos(k * (xc - t))
    M = np.column_stack([s, c])
    (p, q), *_ = np.linalg.lstsq(M, y, rcond=None)
    return np.hypot(p, q), np.arctan2(-q, p)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("stem", help="run name, e.g. ideal_N64 or murnaghan_N128_nnn64")
    ap.add_argument("-d", "--dir", default=None, help="snapshot directory (default output_<stem>)")
    ap.add_argument("-o", "--out", default=None, help="output .mp4 (default movies/<stem>.mp4)")
    ap.add_argument("--fps", type=int, default=20)
    ap.add_argument("--dpi", type=int, default=130)
    args = ap.parse_args()
    os.chdir(os.path.dirname(os.path.abspath(__file__)))

    eos, n, nnn = parse_stem(args.stem)
    with h5py.File("soundwave_{}_N{}.h5".format(eos, n), "r") as f:
        A, lam, ny = f.attrs["amplitude"], f.attrs["wavelength"], int(f.attrs["ny"])
    k = 2.0 * np.pi / lam
    outdir = args.dir or "output_" + args.stem
    files, times = snapshots(outdir)
    if not files:
        raise SystemExit("no snapshots in " + outdir)

    frames = []
    for f, t in zip(files, times):
        xc, prof, mean = column_profile(f, n)
        pert = (prof - mean) / (mean * A)
        ana = np.sin(k * (xc - t))
        amp, phi = sine_fit(xc, pert, k, t)
        frames.append((t, xc, pert, pert - ana, mean, amp, phi))
    rmax = max(np.abs(fr[3]).max() for fr in frames)
    rlim = 1.15 * max(rmax, 1e-3)

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 5.4), sharex=True,
                                   gridspec_kw={"height_ratios": [1.6, 1], "hspace": 0.12})
    xa = np.linspace(0, 1, 800)
    for ax in (ax1, ax2):
        ax.set_xlim(0, 1)
        ax.grid(True, color=GRID, lw=0.6)
        ax.set_axisbelow(True)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        ax.axhline(0, color=AXIS, lw=0.8, zorder=1)
        ax.tick_params(length=2.5, width=0.6)
    ax1.set_ylim(-1.35, 1.35)
    ax1.set_ylabel(r"$(\langle\rho\rangle_y - \langle\rho\rangle)\,/\,(\langle\rho\rangle A)$")
    ax2.set_ylim(-rlim, rlim)
    ax2.set_ylabel("numerical − analytic  [A]")
    ax2.set_xlabel("x")
    ana_line, = ax1.plot([], [], color=INK_2, lw=1.0, zorder=2, label=r"analytic $\sin(k(x-t))$, $c=1$")
    num_line, = ax1.plot([], [], "o-", color=SERIES_1, lw=1.2, ms=2.5 if n <= 64 else 0,
                         zorder=3, label="y-averaged particles ({} per column)".format(ny))
    res_line, = ax2.plot([], [], color=SERIES_1, lw=1.2, zorder=3)
    ax1.legend(loc="upper right", frameon=False, fontsize=8, labelcolor=INK_2)
    title = fig.suptitle("", x=0.08, ha="left", y=0.97, fontsize=11, color=INK)
    info = ax1.text(0.01, 0.03, "", transform=ax1.transAxes, fontsize=8, color=INK_2,
                    family="monospace")

    os.makedirs("movies", exist_ok=True)
    out = args.out or os.path.join("movies", args.stem + ".mp4")
    writer = FFMpegWriter(fps=args.fps, codec="libx264", extra_args=["-pix_fmt", "yuv420p"])
    with writer.saving(fig, out, dpi=args.dpi):
        for t, xc, pert, res, mean, amp, phi in frames:
            ana_line.set_data(xa, np.sin(k * (xa - t)))
            num_line.set_data(xc, pert)
            res_line.set_data(xc, res)
            title.set_text("{}Soundwave {}{}  (n = {}, N_NN = {})    t = {:.2f}".format(
                "GIZMO  " if "_gizmo" in args.stem else "",
                eos, ", masses corrected" if "_corr" in args.stem else "", n, nnn, t))
            info.set_text("<rho>-1 = {:+.4e}   fitted amplitude = {:.5f} A   phase shift = {:+.2e} rad (+ = ahead)"
                          .format(mean - 1.0, amp, phi))
            writer.grab_frame()
    plt.close(fig)
    print("Wrote {} ({} frames, t = {:.3f} .. {:.3f})".format(out, len(frames), times[0], times[-1]))
