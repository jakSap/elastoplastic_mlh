#!/usr/bin/env python3
"""Movie of the y-averaged density and deviatoric stress of one elastic-wave run.

Each frame, from one snapshot:
  top    - y-averaged density perturbation of the working density `rho`
           (the explicitly integrated one), per lattice column,
           (<rho>_y(x) - <rho>) / (<rho> A), against the analytic wave
           sin(k (x - c_l t))
  middle - y-averaged S_xx in units of -(4/3) mu A, i.e. the same sine if the
           hypo-elastic stress stays slaved to the density
  bottom - density residual numerical - analytic, in units of A

<rho> is the domain mean at that time; with corrected masses it is 1 to ~1e-10.
Time is also shown in wave passings t/T, T = lambda/c_l.

    ./makeMovie.py elastic_N64_corr      # -> movies/elastic_N64_corr.mp4
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
SERIES_2 = "#eb6834"

plt.rcParams.update({
    "font.family": "sans-serif", "font.size": 9, "text.color": INK,
    "axes.labelcolor": INK_2, "axes.edgecolor": AXIS, "axes.linewidth": 0.8,
    "xtick.color": MUTED, "ytick.color": MUTED,
    "xtick.labelcolor": INK_2, "ytick.labelcolor": INK_2,
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
})


def ic_name(stem):
    """elastic_N<n>_corr -> soundwave_elastic_N<n>_corr.h5"""
    return "soundwave_{}.h5".format(stem)


def resolution(stem):
    return int([p for p in stem.split("_") if p.startswith("N")][0][1:])


def snapshots(outdir):
    files = glob.glob(os.path.join(outdir, "*.h5"))
    times = []
    for f in files:
        with h5py.File(f, "r") as d:
            times.append(float(np.ravel(d["time"])[0]))
    order = np.argsort(times)
    return [files[i] for i in order], [times[i] for i in order]


def column_profile(path, n, fields=("rho",)):
    """Column centres, y-averaged fields per column, and the domain mean of each."""
    with h5py.File(path, "r") as d:
        x = d["x"][:, 0]
        vals = [d[f][:] for f in fields]
    col = np.floor(x * n).astype(int) % n
    counts = np.maximum(np.bincount(col, minlength=n), 1)
    profs = [np.bincount(col, weights=v, minlength=n) / counts for v in vals]
    means = [v.mean() for v in vals]
    return (np.arange(n) + 0.5) / n, profs, means


def sine_fit(xc, y, k, phase):
    """Amplitude and phase shift of y ~ a sin(k xc - phase - phi).

    phi > 0: the numerical wave is AHEAD of the reference wave.
    """
    s, c = np.sin(k * xc - phase), np.cos(k * xc - phase)
    M = np.column_stack([s, c])
    (p, q), *_ = np.linalg.lstsq(M, y, rcond=None)
    return np.hypot(p, q), np.arctan2(-q, p)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("stem", help="run name, e.g. elastic_N64_corr")
    ap.add_argument("-d", "--dir", default=None, help="snapshot directory (default output_<stem>)")
    ap.add_argument("-o", "--out", default=None, help="output .mp4 (default movies/<stem>.mp4)")
    ap.add_argument("--fps", type=int, default=20)
    ap.add_argument("--dpi", type=int, default=130)
    args = ap.parse_args()
    os.chdir(os.path.dirname(os.path.abspath(__file__)))

    n = resolution(args.stem)
    with h5py.File(ic_name(args.stem), "r") as f:
        A, lam, ny = float(f.attrs["amplitude"]), float(f.attrs["wavelength"]), int(f.attrs["ny"])
        mu, c_l, T = float(f.attrs["mu"]), float(f.attrs["c_l"]), float(f.attrs["period"])
    k = 2.0 * np.pi / lam
    S_SCALE = (4.0 / 3.0) * mu * A      # |S_xx| of the analytic wave
    outdir = args.dir or "output_" + args.stem
    files, times = snapshots(outdir)
    if not files:
        raise SystemExit("no snapshots in " + outdir)

    frames = []
    for f, t in zip(files, times):
        xc, (rho_p, sxx_p), (rho_m, _) = column_profile(f, n, ("rho", "Sxx"))
        pert = (rho_p - rho_m) / (rho_m * A)
        sper = -sxx_p / S_SCALE
        phase = k * c_l * t
        ana = np.sin(k * xc - phase)
        amp, phi = sine_fit(xc, pert, k, phase)
        frames.append((t, xc, pert, sper, pert - ana, rho_m, amp, phi))
    rmax = max(np.abs(fr[4]).max() for fr in frames)
    rlim = 1.15 * max(rmax, 1e-3)

    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(8, 7.0), sharex=True,
                                        gridspec_kw={"height_ratios": [1.5, 1.1, 1.0],
                                                     "hspace": 0.12})
    xa = np.linspace(0, 1, 800)
    for ax in (ax1, ax2, ax3):
        ax.set_xlim(0, 1)
        ax.grid(True, color=GRID, lw=0.6)
        ax.set_axisbelow(True)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        ax.axhline(0, color=AXIS, lw=0.8, zorder=1)
        ax.tick_params(length=2.5, width=0.6)
    ax1.set_ylim(-1.35, 1.35)
    ax1.set_ylabel(r"$(\langle\rho\rangle_y - \langle\rho\rangle)\,/\,(\langle\rho\rangle A)$")
    ax2.set_ylim(-1.35, 1.35)
    ax2.set_ylabel(r"$-\langle S_{xx}\rangle_y\,/\,(\frac{4}{3}\mu A)$")
    ax3.set_ylim(-rlim, rlim)
    ax3.set_ylabel(r"$\rho$: numerical − analytic  [A]")
    ax3.set_xlabel("x")
    ana_line, = ax1.plot([], [], color=INK_2, lw=1.0, zorder=2,
                         label=r"analytic $\sin(k(x-c_l t))$, $c_l$ = {:.5f}".format(c_l))
    num_line, = ax1.plot([], [], "o-", color=SERIES_1, lw=1.2, ms=2.5 if n <= 64 else 0,
                         zorder=3, label="y-averaged particles ({} per column)".format(ny))
    ana2_line, = ax2.plot([], [], color=INK_2, lw=1.0, zorder=2)
    sxx_line, = ax2.plot([], [], "o-", color=SERIES_2, lw=1.2, ms=2.5 if n <= 64 else 0,
                         zorder=3, label=r"$S_{xx}$ (slaved to $\rho$ if it matches)")
    res_line, = ax3.plot([], [], color=SERIES_1, lw=1.2, zorder=3)
    ax1.legend(loc="upper right", frameon=False, fontsize=8, labelcolor=INK_2)
    ax2.legend(loc="upper right", frameon=False, fontsize=8, labelcolor=INK_2)
    title = fig.suptitle("", x=0.08, ha="left", y=0.975, fontsize=11, color=INK)
    info = ax1.text(0.01, 0.02, "", transform=ax1.transAxes, fontsize=7.5, color=INK_2,
                    family="monospace", va="bottom")

    os.makedirs("movies", exist_ok=True)
    out = args.out or os.path.join("movies", args.stem + ".mp4")
    writer = FFMpegWriter(fps=args.fps, codec="libx264", extra_args=["-pix_fmt", "yuv420p"])
    with writer.saving(fig, out, dpi=args.dpi):
        for t, xc, pert, sper, res, mean, amp, phi in frames:
            ana = np.sin(k * xa - k * c_l * t)
            ana_line.set_data(xa, ana)
            ana2_line.set_data(xa, ana)
            num_line.set_data(xc, pert)
            sxx_line.set_data(xc, sper)
            res_line.set_data(xc, res)
            title.set_text("Elastic wave, Tillotson, $\\mu$ = {:.4f}   "
                           "n = {}, $N_{{NN}}$ = 32, masses corrected".format(mu, n))
            info.set_text("t = {:.4f} = {:.2f} passings   <rho>-1 = {:+.3e}\n"
                          "amplitude = {:.5f} A   phase = {:+.2e} rad (+ = ahead)"
                          .format(t, t / T, mean - 1.0, amp, phi))
            writer.grab_frame()
    plt.close(fig)
    print("Wrote {} ({} frames, t = {:.3f} .. {:.3f} = {:.2f} passings)".format(
        out, len(frames), times[0], times[-1], times[-1] / T))
