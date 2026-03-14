#!/usr/bin/env python3

import argparse
import os
import pathlib
import tempfile
import shutil
import numpy as np
import h5py as h5
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from multiprocessing import Pool

NUM_BINS = 20

def prescan_quantity_binned(files, key):
    """Pre-scan all files to find global min/max of bin-averaged values for a quantity."""
    lo, hi = None, None
    for f in files:
        with h5.File(f, 'r') as d:
            pos = d["x"][:]
            vals = d[key][()]
        x = pos[:, 0]
        x_min, x_max = float(x.min()), float(x.max())
        bin_edges = np.linspace(x_min, x_max, NUM_BINS + 1)
        indices = np.digitize(x, bin_edges) - 1
        indices = np.clip(indices, 0, NUM_BINS - 1)
        for b in range(NUM_BINS):
            mask = indices == b
            if mask.any():
                mean_val = float(vals[mask].mean())
                if lo is None or mean_val < lo:
                    lo = mean_val
                if hi is None or mean_val > hi:
                    hi = mean_val
    if lo is None or np.isnan(lo):
        lo = 0.
    if hi is None or np.isnan(hi):
        hi = 1.
    margin = 0.05 * (hi - lo) if hi > lo else 0.5
    return lo - margin, hi + margin


def create_wave_plot(h5File, outDir, ylimits, openBorders=False):
    data = h5.File(h5File, 'r')
    pos = data["x"][:]
    time = data["time"][0]

    quantities = [
        ("rho", data["rho"][()],  r"Density $\varrho$"),
        ("P",   data["P"][()],    r"Pressure $P$"),
        ("u",   data["u"][()],    r"Internal energy $u$"),
        ("Sxx", data["Sxx"][:],   r"$S_{xx}$"),
        ("Sxy", data["Sxy"][:],   r"$S_{xy}$"),
        ("Syy", data["Syy"][:],   r"$S_{yy}$"),
    ]

    x = pos[:, 0]
    x_min, x_max = float(x.min()), float(x.max())
    bin_edges = np.linspace(x_min, x_max, NUM_BINS + 1)
    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    indices = np.digitize(x, bin_edges) - 1
    indices = np.clip(indices, 0, NUM_BINS - 1)

    plt.rcParams.update({'font.size': 12})
    fig, axes = plt.subplots(2, 3, figsize=(10, 6), dpi=300)
    axes_flat = axes.flatten()

    for ax, (key, vals, label) in zip(axes_flat, quantities):
        bin_means = np.full(NUM_BINS, np.nan)
        for b in range(NUM_BINS):
            mask = indices == b
            if mask.any():
                bin_means[b] = vals[mask].mean()

        ax.plot(bin_centers, bin_means, '-o', markersize=3, linewidth=1.2)
        ax.set_title(label)
        ax.set_xlabel("$x$")
        ax.set_ylabel(label)
        if ylimits and key in ylimits:
            ax.set_ylim(ylimits[key])

    data.close()

    fig.suptitle(r"$t = " + f"{time:.4f}" + r"$", fontsize=16)
    plt.tight_layout()
    out = outDir + "/wave" + pathlib.Path(h5File).stem + ".png"
    print("Saving figure to", out)
    plt.savefig(out)
    plt.close()


def _worker_init(tmpdir):
    """Give each pool worker its own TeX cache directory."""
    from matplotlib.texmanager import TexManager
    cache_dir = os.path.join(tmpdir, str(os.getpid()), "tex.cache")
    os.makedirs(cache_dir, exist_ok=True)
    TexManager._texcache = cache_dir


def _worker(task):
    h5File, outDir, ylimits, openBorders = task
    create_wave_plot(h5File, outDir, ylimits, openBorders)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Plot y-averaged simulation quantities vs x-position (binned 1D profiles).")
    parser.add_argument("--simOutputDir", "-d", metavar="string", type=str,
                        help="output directory of simulation", required=True)
    parser.add_argument("--outDir", "-o", metavar="string", type=str,
                        help="output directory for generated plots", default="output")
    parser.add_argument("--openBorders", "-b", action="store_true",
                        help="Adjust plot domain to show all real particles")
    parser.add_argument("--continue", "-c", dest="continue_", action="store_true",
                        help="skip h5 files whose plots already exist in outDir")
    parser.add_argument("--workers", "-w", metavar="int", type=int, default=os.cpu_count(),
                        help="number of parallel worker processes (default: all CPUs)")

    args = parser.parse_args()

    plt.rc('text', usetex=True)

    os.makedirs(args.outDir, exist_ok=True)

    print("Examining files in", args.simOutputDir, "...")

    h5Files = sorted([f for f in pathlib.Path(args.simOutputDir).glob('*.h5')
                       if "NNL" not in str(f) and "Ghost" not in str(f)])

    if args.continue_:
        existing = sorted(pathlib.Path(args.outDir).glob('wave[0-9]*.png'))
        if existing:
            last_stem = existing[-1].stem[len("wave"):]
            h5Files = [f for f in h5Files if f.stem > last_stem]
            print(f"Continuing after {last_stem}, {len(h5Files)} file(s) remaining.")
        else:
            print("No existing plots found in outDir, starting from the beginning.")

    ylimits = None
    if len(h5Files) > 0:
        print(f"Pre-scanning {len(h5Files)} file(s) for y-axis limits (all 6 quantities)...")
        ylimits = {}
        for key in ["rho", "P", "u", "Sxx", "Sxy", "Syy"]:
            lo, hi = prescan_quantity_binned(h5Files, key)
            ylimits[key] = (lo, hi)
            print(f"  {key}: ymin={lo:.4g}, ymax={hi:.4g}")

    tasks = [(f, args.outDir, ylimits, args.openBorders) for f in h5Files]

    nworkers = min(args.workers, len(tasks)) if tasks else 1
    print(f"Rendering {len(tasks)} plot(s) with {nworkers} worker(s)...")
    tmpdir = tempfile.mkdtemp(prefix="mpl_workers_")
    try:
        with Pool(nworkers, initializer=_worker_init, initargs=(tmpdir,)) as pool:
            for i, _ in enumerate(pool.imap_unordered(_worker, tasks), 1):
                print(f"  {i}/{len(tasks)} done", end='\r', flush=True)
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)
    print()

    print("... done.")
