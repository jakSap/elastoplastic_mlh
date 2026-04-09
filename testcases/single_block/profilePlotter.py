#!/usr/bin/env python3

import argparse
import os
import pathlib
import numpy as np
import h5py as h5
import matplotlib
matplotlib.use('Agg')  # non-interactive backend, safe for multiprocessing
import matplotlib.pyplot as plt
from matplotlib import cm as mpl_cm
from multiprocessing import Pool

# Pretty LaTeX labels for known quantities. Falls back to the raw key.
QUANTITY_LABELS = {
    'rho': r'$\varrho$',
    'P': r'$P$',
    'u': r'$u$',
    'm': r'$m$',
    'noi': r'noi',
    'Sxx': r'$S_{xx}$',
    'Sxy': r'$S_{xy}$',
    'Syy': r'$S_{yy}$',
    'h': r'$h$',
    'conditionNumber': r'$N_\mathrm{cond}$',
    'lambdaMin': r'$\lambda_\mathrm{min}$',
    'lambdaMax': r'$\lambda_\mathrm{max}$',
    'v.x': r'$v_x$', 'v.y': r'$v_y$', 'v': r'$|v|$',
    'rhoGrad.x': r'$\nabla\varrho_x$', 'rhoGrad.y': r'$\nabla\varrho_y$',
    'rhoGrad': r'$|\nabla\varrho|$',
    'eigenvecMin.x': r'$e_{\min,x}$', 'eigenvecMin.y': r'$e_{\min,y}$',
    'eigenvecMin': r'$|e_{\min}|$',
    'x.x': r'$x$', 'x.y': r'$y$',
}


def quantityLabel(q):
    return QUANTITY_LABELS.get(q, q)


def resolveQuantity(data, quantity):
    """Return a 1D numpy array of length N for the requested quantity.

    Resolution rules:
      - "key.x" / "key.y": component of a 2D dataset.
      - bare key for a 2D dataset: vector magnitude.
      - bare key for a 1D dataset: the values themselves.
    """
    if '.' in quantity:
        key, comp = quantity.split('.', 1)
        if key not in data:
            raise KeyError(f"unknown h5 key '{key}' (in quantity '{quantity}')")
        arr = data[key]
        if arr.ndim != 2 or arr.shape[1] < 2:
            raise ValueError(f"'{key}' is not a 2D vector dataset, cannot take component '{comp}'")
        if comp == 'x':
            return arr[:, 0]
        if comp == 'y':
            return arr[:, 1]
        raise ValueError(f"unknown component '{comp}' (use 'x' or 'y')")
    else:
        if quantity not in data:
            raise KeyError(f"unknown h5 key '{quantity}'")
        arr = data[quantity]
        if arr.ndim == 1:
            return arr[:]
        if arr.ndim == 2 and arr.shape[1] >= 2:
            v = arr[:]
            return np.sqrt(np.sum(v * v, axis=1))
        raise ValueError(f"cannot reduce dataset '{quantity}' with shape {arr.shape} to a scalar per particle")


def validateQuantity(simOutputDir, quantity):
    """Open the first snapshot once to validate the quantity string up-front."""
    files = sorted(pathlib.Path(simOutputDir).glob('*.h5'))
    if not files:
        raise FileNotFoundError(f"no *.h5 files found in {simOutputDir}")
    with h5.File(files[0], 'r') as data:
        resolveQuantity(data, quantity)  # raises if invalid


def createProfilePlot(h5File, outDir, quantity, axis, index, useLog, plotMean):
    data = h5.File(h5File, 'r')
    pos = data["x"][:]
    time = data["time"][0]
    N = pos.shape[0]

    plot_axis = 0 if axis == 'x' else 1
    ortho_axis = 1 - plot_axis

    ortho = pos[:, ortho_axis]
    plot_coord = pos[:, plot_axis]

    values = resolveQuantity(data, quantity)

    plt.rcParams.update({'font.size': 18})
    fig, ax = plt.subplots(figsize=(8, 5), dpi=200)

    if plotMean:
        if index >= 0:
            print(f"Note: --mean is set, ignoring -i {index} for {h5File.name}.")

        # Bin the *plot* axis at width edge_plot / sqrt(N) and aggregate.
        plot_min = float(plot_coord.min())
        plot_max = float(plot_coord.max())
        n_bins_plot = max(1, int(round(np.sqrt(N))))
        plot_bin_edges = np.linspace(plot_min, plot_max, n_bins_plot + 1)
        plot_bin_idx = np.clip(np.digitize(plot_coord, plot_bin_edges) - 1,
                               0, n_bins_plot - 1)

        centers, means, lo, hi = [], [], [], []
        for b in range(n_bins_plot):
            mask = (plot_bin_idx == b)
            if not np.any(mask):
                continue
            v = values[mask]
            centers.append(0.5 * (plot_bin_edges[b] + plot_bin_edges[b + 1]))
            means.append(float(np.mean(v)))
            lo.append(float(np.percentile(v, 2.5)))
            hi.append(float(np.percentile(v, 97.5)))
        centers = np.array(centers)
        means = np.array(means)
        lo = np.array(lo)
        hi = np.array(hi)

        ax.fill_between(centers, lo, hi, color='tab:blue', alpha=0.25,
                        label=r'95\% envelope')
        ax.plot(centers, means, color='tab:blue', linewidth=2.0, label='mean')
        ax.legend(loc='best', fontsize=12)
    else:
        ortho_min = float(ortho.min())
        ortho_max = float(ortho.max())

        n_bins = max(1, int(round(np.sqrt(N))))
        bin_edges = np.linspace(ortho_min, ortho_max, n_bins + 1)
        bin_idx = np.clip(np.digitize(ortho, bin_edges) - 1, 0, n_bins - 1)

        highlight_bin = -1
        if index >= 0:
            if index < N:
                highlight_bin = int(bin_idx[index])
            else:
                print(f"Warning: index {index} out of range for {h5File.name} (N={N}), ignoring -i.")

        cmap = mpl_cm.get_cmap('viridis')

        # Draw non-highlighted bins first.
        for b in range(n_bins):
            if b == highlight_bin:
                continue
            mask = (bin_idx == b)
            if not np.any(mask):
                continue
            order = np.argsort(plot_coord[mask])
            xs = plot_coord[mask][order]
            ys = values[mask][order]
            color = cmap(b / max(1, n_bins - 1))
            if highlight_bin >= 0:
                ax.plot(xs, ys, color=color, linewidth=0.8, alpha=0.45)
            else:
                ax.plot(xs, ys, color=color, linewidth=1.0, alpha=0.9)

        # Draw the highlighted bin last so it sits on top.
        if highlight_bin >= 0:
            mask = (bin_idx == highlight_bin)
            if np.any(mask):
                order = np.argsort(plot_coord[mask])
                xs = plot_coord[mask][order]
                ys = values[mask][order]
                ax.plot(xs, ys, color='tab:red', linewidth=2.5,
                        label=f'bin of particle {index}')
                ax.legend(loc='best', fontsize=12)

    if useLog:
        ax.set_yscale('log')

    qlabel = quantityLabel(quantity)
    title_kind = "Mean profile" if plotMean else "Profile"
    ax.set_title(rf"{title_kind} of {qlabel} vs ${axis}$ at $t = {time:.4f}$")
    ax.set_xlabel(rf"${axis}$")
    ax.set_ylabel(qlabel)
    ax.grid()
    plt.tight_layout()

    q_safe = quantity.replace('.', '_')
    mean_tag = "mean_" if plotMean else ""
    outname = f"profile_{q_safe}_vs_{axis}_{mean_tag}{pathlib.Path(h5File).stem}.png"
    outpath = os.path.join(outDir, outname)
    print("Saving figure to", outpath)
    plt.savefig(outpath)
    # plt.show()
    plt.close(fig)
    data.close()


def _worker(args_tuple):
    createProfilePlot(*args_tuple)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Plot 1D profiles of an h5 quantity along x or y, "
                    "binning the orthogonal axis at width edge_length / sqrt(N).")
    parser.add_argument("--simOutputDir", "-d", metavar="string", type=str,
                        help="output directory of simulation", required=True)
    parser.add_argument("--outDir", "-o", metavar="string", type=str,
                        help="output directory for generated profile plots", default="output")
    parser.add_argument("--quantity", "-q", metavar="string", type=str, required=True,
                        help="quantity to plot, e.g. 'rho', 'P', 'v.x', 'rhoGrad.y', "
                             "or a bare vector key like 'v' for its magnitude")
    parser.add_argument("--axis", "-a", metavar="x|y", type=str, required=True,
                        choices=['x', 'y'],
                        help="coordinate axis to plot the profile against")
    parser.add_argument("--index", "-i", metavar="int", type=int, default=-1,
                        help="particle index whose containing bin is highlighted with a thicker line")
    parser.add_argument("--tstart", "-t", metavar="float", type=float, default=None,
                        help="only plot snapshots with simulation time >= tstart")
    parser.add_argument("--continue", "-c", dest="continue_", action="store_true",
                        help="skip snapshots whose profile plot already exists in outDir")
    parser.add_argument("--workers", "-w", metavar="int", type=int, default=os.cpu_count(),
                        help="number of parallel worker processes (default: all CPUs)")
    parser.add_argument("--log", action="store_true", help="use a logarithmic y-axis")
    parser.add_argument("--mean", action="store_true",
                        help="plot a single averaged profile (mean over particles in each plot-axis bin) "
                             "with a shaded 95%% envelope (2.5--97.5 percentile)")
    args = parser.parse_args()

    plt.rc('text', usetex=True)
    plt.rc('text.latex', preamble=r'\usepackage{amsmath}')

    print("Examining files in", args.simOutputDir, "...")

    # Validate the quantity once on the main process so we fail fast and not inside workers.
    validateQuantity(args.simOutputDir, args.quantity)

    h5Files = sorted([f for f in pathlib.Path(args.simOutputDir).glob('*.h5')
                      if "NNL" not in str(f) and "Ghost" not in str(f)])

    q_safe = args.quantity.replace('.', '_')
    mean_tag = "mean_" if args.mean else ""
    prefix = f"profile_{q_safe}_vs_{args.axis}_{mean_tag}"

    os.makedirs(args.outDir, exist_ok=True)

    if args.continue_:
        existing = sorted(pathlib.Path(args.outDir).glob(f'{prefix}[0-9]*.png'))
        if existing:
            last_stem = existing[-1].stem[len(prefix):]
            h5Files = [f for f in h5Files if f.stem > last_stem]
            print(f"Continuing after {last_stem}, {len(h5Files)} file(s) remaining.")
        else:
            print("No existing plots found in outDir, starting from the beginning.")

    if args.tstart is not None:
        before = len(h5Files)
        h5Files = [f for f in h5Files
                   if h5.File(f, 'r')["time"][0] >= args.tstart]
        print(f"tstart filter: {before} -> {len(h5Files)} files.")

    if not h5Files:
        print("No files to plot.")
        raise SystemExit(0)

    print(f"Plotting {len(h5Files)} snapshot(s) with {args.workers} worker(s) ...")

    tasks = [(f, args.outDir, args.quantity, args.axis, args.index, args.log, args.mean)
             for f in h5Files]

    if args.workers <= 1:
        for t in tasks:
            _worker(t)
    else:
        with Pool(args.workers) as pool:
            pool.map(_worker, tasks)

    print("... done.")
