#!/usr/bin/env python3
"""
Diagnostic analysis for the surface-stability crescent test.

Reads all dumps in output/, pulls condition number, h_i, and noi, splits the
particles by their initial geometric region (convex outer / concave inner /
horn / interior) using the 'region' label stored in the IC file, and plots
how each diagnostic evolves per region.

Usage:
    python3 analyze.py                   # all standard plots
    python3 analyze.py --snapshots 0 5 10 15 20
"""

import argparse
import glob
import os
import numpy as np
import matplotlib.pyplot as plt
import h5py as h5


REGION_NAMES = {0: 'interior', 1: 'convex (thick lobe)',
                2: 'concave inner', 3: 'convex (thin lobe)'}
REGION_COLORS = {0: '#888888', 1: '#1f77b4', 2: '#d62728', 3: '#2ca02c'}


def load_region(ic_file):
    with h5.File(ic_file, 'r') as f:
        return f['region'][:]


def sorted_dumps(output_dir):
    files = sorted(glob.glob(os.path.join(output_dir, '??????.h5')))
    if not files:
        raise FileNotFoundError(f"No HDF5 outputs in {output_dir}")
    return files


def load_time(fname):
    with h5.File(fname, 'r') as f:
        return float(f['time'][0])


def time_series_by_region(files, region, kernel_size):
    """Collect max cond(E), max h_i/h_0, fraction clamped-high per region per time."""
    times = []
    stats = {r: {'cond_max': [], 'cond_p95': [],
                 'h_over_h0_max': [], 'h_over_h0_p95': [],
                 'noi_min': [], 'noi_mean': []} for r in REGION_NAMES}

    for f in files:
        with h5.File(f, 'r') as d:
            times.append(float(d['time'][0]))
            cond = d['conditionNumber'][:] if 'conditionNumber' in d else None
            sml = d['sml'][:] if 'sml' in d else None
            noi = d['noi'][:] if 'noi' in d else None

        for r in REGION_NAMES:
            mask = (region == r)
            if mask.sum() == 0:
                for k in stats[r]:
                    stats[r][k].append(np.nan)
                continue
            if cond is not None:
                c = cond[mask]
                stats[r]['cond_max'].append(float(np.max(c)))
                stats[r]['cond_p95'].append(float(np.percentile(c, 95)))
            else:
                stats[r]['cond_max'].append(np.nan)
                stats[r]['cond_p95'].append(np.nan)
            if sml is not None:
                h = sml[mask] / kernel_size
                stats[r]['h_over_h0_max'].append(float(np.max(h)))
                stats[r]['h_over_h0_p95'].append(float(np.percentile(h, 95)))
            else:
                stats[r]['h_over_h0_max'].append(np.nan)
                stats[r]['h_over_h0_p95'].append(np.nan)
            if noi is not None:
                n = noi[mask]
                stats[r]['noi_min'].append(int(np.min(n)))
                stats[r]['noi_mean'].append(float(np.mean(n)))
            else:
                stats[r]['noi_min'].append(np.nan)
                stats[r]['noi_mean'].append(np.nan)

    return np.array(times), stats


def plot_timeseries(times, stats, save_path):
    fig, axes = plt.subplots(3, 1, figsize=(9, 10), dpi=150, sharex=True)

    for r in sorted(REGION_NAMES):
        axes[0].plot(times, stats[r]['cond_p95'], color=REGION_COLORS[r],
                     label=f"{REGION_NAMES[r]} (p95)")
        axes[0].plot(times, stats[r]['cond_max'], color=REGION_COLORS[r],
                     ls='--', alpha=0.5)
        axes[1].plot(times, stats[r]['h_over_h0_p95'], color=REGION_COLORS[r],
                     label=f"{REGION_NAMES[r]} (p95)")
        axes[1].plot(times, stats[r]['h_over_h0_max'], color=REGION_COLORS[r],
                     ls='--', alpha=0.5)
        axes[2].plot(times, stats[r]['noi_mean'], color=REGION_COLORS[r],
                     label=f"{REGION_NAMES[r]} (mean)")
        axes[2].plot(times, stats[r]['noi_min'], color=REGION_COLORS[r],
                     ls='--', alpha=0.5)

    axes[0].set_ylabel(r'cond$(E)$')
    axes[0].set_yscale('log')
    axes[0].legend(ncol=2, fontsize=8)
    axes[0].grid(alpha=0.3)

    axes[1].set_ylabel(r'$h_i / h_0$')
    axes[1].legend(ncol=2, fontsize=8)
    axes[1].grid(alpha=0.3)

    axes[2].set_ylabel('noi')
    axes[2].set_xlabel('time')
    axes[2].legend(ncol=2, fontsize=8)
    axes[2].grid(alpha=0.3)

    fig.suptitle("Per-region diagnostics (dashed = max/min, solid = p95/mean)")
    plt.tight_layout()
    fig.savefig(save_path)
    print(f"Wrote {save_path}")


def plot_snapshots(files, snap_times, kernel_size, save_path):
    file_times = np.array([load_time(f) for f in files])

    def nearest(t):
        return files[int(np.argmin(np.abs(file_times - t)))]

    ncols = len(snap_times)
    fig, axes = plt.subplots(3, ncols, figsize=(4 * ncols, 11), dpi=150,
                             squeeze=False)

    for col, tw in enumerate(snap_times):
        fname = nearest(tw)
        with h5.File(fname, 'r') as d:
            t = float(d['time'][0])
            pos = d['x'][:]
            cond = d['conditionNumber'][:]
            sml = d['sml'][:]
            noi = d['noi'][:]

        for row, (field, label, cmap, vmin, vmax) in enumerate([
            (np.log10(np.clip(cond, 1, None)), r'log10 cond$(E)$', 'viridis', 0, 3),
            (sml / kernel_size,                r'$h_i / h_0$',     'plasma',  0.9, 2.0),
            (noi,                              'noi',              'cividis', 0, 80),
        ]):
            ax = axes[row, col]
            sc = ax.scatter(pos[:, 0], pos[:, 1], c=field, s=1.5,
                            cmap=cmap, vmin=vmin, vmax=vmax)
            ax.set_aspect('equal')
            if row == 0:
                ax.set_title(f"t = {t:.2f}")
            if col == 0:
                ax.set_ylabel(label)
            if col == ncols - 1:
                plt.colorbar(sc, ax=ax, fraction=0.04)

    plt.tight_layout()
    fig.savefig(save_path)
    print(f"Wrote {save_path}")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--dir', '-d', default='output',
                        help='output directory (default: output)')
    parser.add_argument('--ic', default='lune_deltap0.05-2D.h5',
                        help='IC file with region labels')
    parser.add_argument('--kernel-size', type=float, default=0.1,
                        help='global kernel size h_0 (default: 0.1)')
    parser.add_argument('--snapshots', nargs='+', type=float,
                        default=[0, 2, 5, 10, 20])
    args = parser.parse_args()

    script_dir = os.path.dirname(os.path.abspath(__file__))
    ic_path = os.path.join(script_dir, args.ic)
    output_dir = os.path.join(script_dir, args.dir)

    region = load_region(ic_path)
    files = sorted_dumps(output_dir)
    print(f"Loaded {len(files)} dumps from {output_dir}")

    times, stats = time_series_by_region(files, region, args.kernel_size)
    plot_timeseries(times, stats,
                    os.path.join(script_dir, 'diagnostics_timeseries.png'))
    plot_snapshots(files, args.snapshots, args.kernel_size,
                   os.path.join(script_dir, 'diagnostics_snapshots.png'))
