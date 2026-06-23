#!/usr/bin/env python3
"""
4x4 panel movie of the tensile-instability sweep.
  Rows : periodic_cubic | periodic_wc2 | free_cubic | free_wc2
  Cols : strain rate 0.01 | 0.1 | 0.5 | 1
  Color: density rho, shared global scale clipped to [0, 2]

Usage:
  python3 make_sweep_movie.py [--sweepdir DIR] [--outdir DIR] [--fps N] [--workers N]
"""

import os, glob, argparse, subprocess
from multiprocessing import Pool

import numpy as np
import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

BASE = os.path.dirname(os.path.abspath(__file__))

COMBOS      = ["periodic_cubic", "periodic_wc2", "free_cubic", "free_wc2"]
STRAINRATES = ["0.01", "0.1", "0.5", "1"]

ROW_LABELS = [
    "cubic\nperiodic",
    "WC2\nperiodic",
    "cubic\nfree",
    "WC2\nfree",
]
COL_LABELS = [f"ṡ = {sr}" for sr in STRAINRATES]

EXIT_COLORS = {0: "#2ecc71", 6: "#e74c3c", 8: "#e67e22", 137: "#9b59b6"}


def get_snapshots(sweepdir, combo, sr):
    d = os.path.join(sweepdir, f"{combo}_sr{sr}", "output")
    return sorted(glob.glob(os.path.join(d, "*.h5")))


def get_exit_code(sweepdir, combo, sr):
    p = os.path.join(sweepdir, f"{combo}_sr{sr}", "exit_code")
    try:
        return int(open(p).read().strip())
    except Exception:
        return None


def read_rho_x_t(path):
    with h5py.File(path, 'r') as f:
        x   = f['x'][:]
        rho = f['rho'][:]
        t   = float(f['time'][0])
    return x, rho, t


def render_frame(task):
    frame_idx, snapmap, exit_codes, vmin, vmax, xlims, frames_dir = task

    fig, axes = plt.subplots(4, 4, figsize=(20, 7), dpi=110,
                             gridspec_kw=dict(hspace=0.08, wspace=0.06))

    current_t = None

    for row in range(4):
        for col in range(4):
            ax = axes[row, col]
            path, is_last = snapmap[(row, col)]

            if path is not None:
                x, rho, snap_t = read_rho_x_t(path)
                ax.scatter(x[:, 0], x[:, 1], c=rho, s=0.2,
                           vmin=vmin, vmax=vmax, cmap='plasma',
                           linewidths=0, rasterized=True)
                if current_t is None or not is_last:
                    current_t = snap_t
                if is_last:
                    ec = exit_codes[(row, col)]
                    ec_col = EXIT_COLORS.get(ec, 'white')
                    for spine in ax.spines.values():
                        spine.set_edgecolor(ec_col)
                        spine.set_linewidth(2.0)
                    ax.text(0.98, 0.97, f"exit={ec}",
                            transform=ax.transAxes, ha='right', va='top',
                            fontsize=5, color=ec_col,
                            bbox=dict(boxstyle='round,pad=0.1', fc='black', alpha=0.5))
            else:
                ax.set_facecolor('#0a0a0a')

            ax.set_xlim(xlims[col])
            ax.set_ylim(-1.35, 1.35)
            ax.tick_params(labelsize=5, length=2, pad=1)

            if row < 3:
                ax.set_xticklabels([])
            if col > 0:
                ax.set_yticklabels([])

            if row == 0:
                ax.set_title(COL_LABELS[col], fontsize=8, pad=3)
            if col == 0:
                ax.set_ylabel(ROW_LABELS[row], fontsize=7, labelpad=2)

    t_str = f"t = {current_t:.3f}" if current_t is not None else ""
    fig.suptitle(f"Tensile instability — density ρ,   {t_str}   "
                 f"(exit border: ✓ green | panic 6 red | FPE 8 orange | OOM 137 purple)",
                 fontsize=8, y=0.995)

    # Shared colorbar
    sm = plt.cm.ScalarMappable(cmap='plasma',
                                norm=plt.Normalize(vmin=vmin, vmax=vmax))
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=axes[:, -1], fraction=0.04, pad=0.01,
                        shrink=0.9, aspect=40)
    cbar.set_label(r'$\rho$', fontsize=8)
    cbar.ax.tick_params(labelsize=6)

    outpath = os.path.join(frames_dir, f"frame_{frame_idx:04d}.png")
    plt.savefig(outpath, dpi=110, bbox_inches='tight', facecolor='#1a1a2e')
    plt.close(fig)
    return frame_idx


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sweepdir", default=None,
                    help="sweep output dir (default: most recent YYYYMMDD under testcase dir)")
    ap.add_argument("--outdir", default=None,
                    help="output dir for frames + movie (default: <sweepdir>/movie)")
    ap.add_argument("--fps",     type=int, default=8)
    ap.add_argument("--workers", type=int, default=os.cpu_count())
    args = ap.parse_args()

    if args.sweepdir is None:
        # pick latest YYYYMMDD directory
        candidates = sorted(glob.glob(os.path.join(BASE, "2???????/")))
        if not candidates:
            print("No sweep date directory found."); return
        args.sweepdir = candidates[-1].rstrip("/")
    print(f"Sweep dir: {args.sweepdir}")

    if args.outdir is None:
        args.outdir = os.path.join(args.sweepdir, "movie")
    frames_dir = os.path.join(args.outdir, "frames")
    os.makedirs(frames_dir, exist_ok=True)

    # Collect snapshots and exit codes
    all_snaps  = {}   # (row, col) -> [path, ...]
    exit_codes = {}   # (row, col) -> int|None
    for row, combo in enumerate(COMBOS):
        for col, sr in enumerate(STRAINRATES):
            all_snaps[(row, col)]  = get_snapshots(args.sweepdir, combo, sr)
            exit_codes[(row, col)] = get_exit_code(args.sweepdir, combo, sr)

    n_frames = max((len(v) for v in all_snaps.values()), default=0)
    if n_frames == 0:
        print("No snapshots found."); return
    print(f"Total frames: {n_frames}")

    # Global rho range (sample first+last snap per run, clip to [0,2])
    print("Scanning rho range...")
    rho_max = 0.0
    for snaps in all_snaps.values():
        for path in ([snaps[0], snaps[-1]] if len(snaps) > 1 else snaps):
            with h5py.File(path, 'r') as f:
                rho_max = max(rho_max, float(np.nanpercentile(f['rho'][:], 99.5)))
    vmin, vmax = 0.0, min(rho_max, 2.5)
    print(f"Color scale: [{vmin:.3f}, {vmax:.3f}]")

    # Per-column x limits (from all snapshots in that column)
    print("Scanning x extents per column...")
    xlims = []
    for col, sr in enumerate(STRAINRATES):
        xlo, xhi = np.inf, -np.inf
        for row in range(4):
            for path in all_snaps[(row, col)]:
                with h5py.File(path, 'r') as f:
                    xi = f['x'][:, 0]
                xlo = min(xlo, float(xi.min()))
                xhi = max(xhi, float(xi.max()))
        margin = 0.04 * (xhi - xlo)
        xlims.append((xlo - margin, xhi + margin))
        print(f"  sr={sr}: x ∈ [{xlo:.2f}, {xhi:.2f}]")

    # Build per-frame snapmap: (row,col) -> (path|None, is_last_frame)
    tasks = []
    for fi in range(n_frames):
        snapmap = {}
        for row in range(4):
            for col in range(4):
                snaps = all_snaps[(row, col)]
                if not snaps:
                    snapmap[(row, col)] = (None, True)
                elif fi < len(snaps):
                    snapmap[(row, col)] = (snaps[fi], fi == len(snaps) - 1)
                else:
                    snapmap[(row, col)] = (snaps[-1], True)   # hold last frame
        tasks.append((fi, snapmap, exit_codes, vmin, vmax, xlims, frames_dir))

    print(f"Rendering {n_frames} frames with {args.workers} workers...")
    with Pool(args.workers) as pool:
        for i, _ in enumerate(pool.imap_unordered(render_frame, tasks), 1):
            print(f"  {i}/{n_frames}", end='\r', flush=True)
    print()

    movie_path = os.path.join(args.outdir, "tensile_sweep.mp4")
    cmd = ["ffmpeg", "-y",
           "-framerate", str(args.fps),
           "-i", os.path.join(frames_dir, "frame_%04d.png"),
           "-vf", "scale=trunc(iw/2)*2:trunc(ih/2)*2",
           "-c:v", "libx264", "-pix_fmt", "yuv420p", "-crf", "17",
           movie_path]
    print("Assembling:", " ".join(cmd))
    subprocess.run(cmd, check=True)
    print(f"Movie saved → {movie_path}")


if __name__ == "__main__":
    main()
