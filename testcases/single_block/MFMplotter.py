#!/usr/bin/env python3
"""
MFMplotter.py — grid-based plotter using the meshless finite mass (MFM)
volume partition of unity from Hopkins (2015) / Martin's master thesis
(eqs. 2.25-2.27). For each grid cell center the chosen quantity is
reconstructed as

    Q(x) = sum_j Q_j * W(|x - x_j|, h) / sum_k W(|x - x_k|, h)

with the same 2D cubic spline kernel used by the demonstrator
(Kernel::cubicSpline in demonstrator/src/Particles.cpp).
"""

import argparse
import os
import pathlib
import numpy as np
import h5py as h5
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.axes_grid1 import make_axes_locatable
from multiprocessing import Pool
from scipy.spatial import cKDTree
from scipy.sparse import coo_matrix


# ----------------------------------------------------------------------
# Kernel — must match Kernel::cubicSpline in demonstrator/src/Particles.cpp
# ----------------------------------------------------------------------
def cubic_spline_2d(r, h):
    """2D cubic spline kernel with support radius h (i.e. q = r/(h/2) <= 2)."""
    h2 = h / 2.0
    sigma = 10.0 / (7.0 * np.pi * h2 * h2)
    q = r / h2
    W = np.zeros_like(r, dtype=np.float64)
    m1 = (q >= 0.0) & (q <= 1.0)
    m2 = (q > 1.0) & (q < 2.0)
    W[m1] = sigma * (1.0 - 1.5 * q[m1] ** 2 * (1.0 - 0.5 * q[m1]))
    W[m2] = sigma / 4.0 * (2.0 - q[m2]) ** 3
    return W


# ----------------------------------------------------------------------
# Per-particle quantity loader for the chosen mode
# ----------------------------------------------------------------------
def load_quantity(data, key):
    if key == "S":
        Sxx = data["Sxx"][:]
        Sxy = data["Sxy"][:]
        Syy = data["Syy"][:]
        return np.sqrt(Sxx**2 + Sxy**2 + Syy**2)
    return np.asarray(data[key][()], dtype=np.float64)


def quantity_label(key):
    return {
        "rho": (r"$\varrho$",       "rho"),
        "P":   (r"$P$",             "P"),
        "u":   (r"$u$",             "u"),
        "S":   (r"$|S|$",           "S"),
        "noi": (r"$N_\mathrm{int}$","noi"),
    }[key]


def get_quantity_key(args):
    if args.pressure: return "P"
    if args.energy:   return "u"
    if args.stress:   return "S"
    if args.noi:      return "noi"
    return "rho"


# ----------------------------------------------------------------------
# Reconstruction on a regular grid
# ----------------------------------------------------------------------
def reconstruct_on_grid(pos, vals, grid_centers, h):
    """
    pos          : (N, 2)  particle positions
    vals         : (N,)    per-particle quantity values
    grid_centers : (M, 2)  cell centers (flattened)
    h            : float   kernel support radius
    Returns      : (M,)    reconstructed values (NaN where no neighbours)

    Uses sparse_distance_matrix between the grid-centers tree and the
    particle tree to obtain all (cell, particle) pairs within distance h
    in one shot, then accumulates omega and sum(Q*W) via np.add.at.
    """
    M = grid_centers.shape[0]
    N = pos.shape[0]

    grid_tree = cKDTree(grid_centers)
    part_tree = cKDTree(pos)

    # COO sparse matrix of pairwise distances <= h
    dmat = grid_tree.sparse_distance_matrix(part_tree, max_distance=h,
                                            output_type='coo_matrix')
    if dmat.nnz == 0:
        return np.full(M, np.nan)

    rows = dmat.row          # cell index
    cols = dmat.col          # particle index
    dists = dmat.data        # distance

    W = cubic_spline_2d(dists, h)

    omega = np.zeros(M, dtype=np.float64)
    QW    = np.zeros(M, dtype=np.float64)
    np.add.at(omega, rows, W)
    np.add.at(QW,    rows, vals[cols] * W)

    out = np.full(M, np.nan)
    nz = omega > 0.0
    out[nz] = QW[nz] / omega[nz]
    return out


# ----------------------------------------------------------------------
# Workers — split into "compute the reconstructed grid" and "render it"
# so the colormap limits can be derived from the actual cell values
# (not the particle values, which contain outliers).
# ----------------------------------------------------------------------
def _compute_worker(task):
    (h5File, qkey, h, grid_centers, ny, nx) = task
    with h5.File(h5File, 'r') as data:
        pos  = data["x"][:]
        time = float(data["time"][0])
        vals = load_quantity(data, qkey)
    Q_flat = reconstruct_on_grid(pos, vals, grid_centers, h)
    return (str(h5File), time, Q_flat.reshape(ny, nx))


def _render_worker(task):
    (h5File, outDir, qkey, time, Q,
     xedges, yedges, vmin, vmax, prefix) = task

    plt.rcParams.update({'font.size': 18})
    fig, ax = plt.subplots(figsize=(7, 6), dpi=200)

    cmap = cm.get_cmap('viridis').copy()
    cmap.set_bad('white')

    pcm = ax.pcolormesh(xedges, yedges, np.ma.masked_invalid(Q),
                        cmap=cmap, vmin=vmin, vmax=vmax,
                        shading='flat')

    ax.set_aspect('equal')
    ax.set_xlim(xedges[0], xedges[-1])
    ax.set_ylim(yedges[0], yedges[-1])

    label, _ = quantity_label(qkey)
    ax.set_title(label + r" at $t = " + f"{time:.4f}" + r"$")
    ax.set_xlabel("$x$")
    ax.set_ylabel("$y$")

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    fig.colorbar(pcm, cax=cax)

    plt.tight_layout()
    out = os.path.join(outDir, "MFM" + prefix + pathlib.Path(h5File).stem + ".png")
    print("Saving figure to", out)
    plt.savefig(out)
    plt.close(fig)


# ----------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Plot particle quantities reconstructed on a regular "
                    "grid using the MFM volume partition of unity.")
    parser.add_argument("--simOutputDir", "-d", metavar="string", type=str,
                        required=True, help="output directory of simulation")
    parser.add_argument("--outDir", "-o", metavar="string", type=str,
                        default="output",
                        help="output directory for generated plots")
    parser.add_argument("--kernel", "-H", metavar="float", type=float,
                        required=True,
                        help="kernel support radius h used for reconstruction")
    parser.add_argument("--cellSize", "-c", metavar="float", type=float,
                        default=None,
                        help="grid cell size (default: equal to --kernel)")
    # Quantity selection (single, mutually exclusive)
    qgroup = parser.add_mutually_exclusive_group()
    qgroup.add_argument("--pressure", "-P", action="store_true",
                        help="plot pressure instead of density")
    qgroup.add_argument("--energy", "-u", action="store_true",
                        help="plot internal energy instead of density")
    qgroup.add_argument("--stress", "-S", action="store_true",
                        help="plot |S| = sqrt(Sxx^2 + Sxy^2 + Syy^2)")
    qgroup.add_argument("--noi", "-n", action="store_true",
                        help="plot number of interactions instead of density")
    # Bookkeeping
    parser.add_argument("--continue", "-C", dest="continue_", action="store_true",
                        help="skip h5 files whose plots already exist in outDir")
    parser.add_argument("--workers", "-w", metavar="int", type=int,
                        default=os.cpu_count(),
                        help="number of parallel worker processes")
    parser.add_argument("--tstart", "-t", metavar="float", type=float, default=None,
                        help="only plot files with simulation time >= tstart")

    args = parser.parse_args()

    h = args.kernel
    cellSize = args.cellSize if args.cellSize is not None else h
    if cellSize > h:
        print(f"WARNING: cellSize ({cellSize}) > kernel h ({h}); kernel is undersampled.")
    if cellSize < h / 8.0:
        print(f"WARNING: cellSize ({cellSize}) < h/8; reconstruction will be very slow.")

    plt.rc('text', usetex=True)

    qkey = get_quantity_key(args)
    _, prefix = quantity_label(qkey)

    print("Examining files in", args.simOutputDir, "...")
    h5Files = sorted([f for f in pathlib.Path(args.simOutputDir).glob('*.h5')
                      if "NNL" not in str(f) and "Ghost" not in str(f)])

    if args.continue_:
        existing = sorted(pathlib.Path(args.outDir).glob(f'MFM{prefix}[0-9]*.png'))
        if existing:
            last_stem = existing[-1].stem[len("MFM" + prefix):]
            h5Files = [f for f in h5Files if f.stem > last_stem]
            print(f"Continuing after {last_stem}, {len(h5Files)} file(s) remaining.")
        else:
            print("No existing plots found in outDir, starting from the beginning.")

    if args.tstart is not None:
        before = len(h5Files)
        h5Files = [f for f in h5Files
                   if h5.File(f, 'r')["time"][0] >= args.tstart]
        print(f"--tstart={args.tstart}: kept {len(h5Files)}/{before} file(s).")

    if not h5Files:
        print("No h5 files to process. Exiting.")
        raise SystemExit(0)

    # ------------------------------------------------------------------
    # Phase 1: domain pre-scan (positions only — cheap)
    # ------------------------------------------------------------------
    print(f"Pre-scanning {len(h5Files)} file(s) for domain bounds...")
    xmin = ymin = np.inf
    xmax = ymax = -np.inf
    for f in h5Files:
        with h5.File(f, 'r') as d:
            pos = d["x"][:]
        if pos.size:
            xmin = min(xmin, float(pos[:, 0].min()))
            xmax = max(xmax, float(pos[:, 0].max()))
            ymin = min(ymin, float(pos[:, 1].min()))
            ymax = max(ymax, float(pos[:, 1].max()))
    print(f"  domain: x in [{xmin:.4g}, {xmax:.4g}], y in [{ymin:.4g}, {ymax:.4g}]")

    # ------------------------------------------------------------------
    # Phase 2: build the grid (once, shared by all frames)
    # ------------------------------------------------------------------
    nx = max(1, int(np.ceil((xmax - xmin) / cellSize)))
    ny = max(1, int(np.ceil((ymax - ymin) / cellSize)))
    xedges = xmin + np.arange(nx + 1) * cellSize
    yedges = ymin + np.arange(ny + 1) * cellSize
    xc = 0.5 * (xedges[:-1] + xedges[1:])
    yc = 0.5 * (yedges[:-1] + yedges[1:])
    Xc, Yc = np.meshgrid(xc, yc)  # shape (ny, nx)
    grid_centers = np.column_stack([Xc.ravel(), Yc.ravel()])  # (ny*nx, 2)
    print(f"  grid: nx={nx}, ny={ny}, cellSize={cellSize}, total cells={nx*ny}")

    os.makedirs(args.outDir, exist_ok=True)

    # ------------------------------------------------------------------
    # Phase 3: reconstruct cell values for every frame (parallel)
    # We must do this BEFORE picking vmin/vmax so that the colormap
    # limits reflect the actual cell values rather than particle outliers.
    # ------------------------------------------------------------------
    nworkers = min(args.workers, len(h5Files))
    compute_tasks = [(f, qkey, h, grid_centers, ny, nx) for f in h5Files]
    print(f"Reconstructing {len(compute_tasks)} frame(s) with {nworkers} worker(s)...")
    results = []  # list of (h5File_str, time, Q_2d)
    with Pool(nworkers) as pool:
        for i, r in enumerate(pool.imap_unordered(_compute_worker, compute_tasks), 1):
            results.append(r)
            print(f"  reconstruct {i}/{len(compute_tasks)}", end='\r', flush=True)
    print()

    # ------------------------------------------------------------------
    # Phase 4: derive vmin/vmax from the reconstructed cell values
    # ------------------------------------------------------------------
    all_finite = np.concatenate(
        [Q[np.isfinite(Q)].ravel() for _, _, Q in results]
    )
    if all_finite.size == 0:
        vmin, vmax = 0.0, 1.0
    else:
        vmin = float(all_finite.min())
        vmax = float(all_finite.max())
    if vmax <= vmin:
        vmin, vmax = vmin - 0.5, vmin + 0.5
    print(f"  cell-value range: vmin={vmin:.4g}, vmax={vmax:.4g}")

    # ------------------------------------------------------------------
    # Phase 5: render in parallel using the cell-derived limits
    # ------------------------------------------------------------------
    render_tasks = [(h5File, args.outDir, qkey, time, Q,
                     xedges, yedges, vmin, vmax, prefix)
                    for (h5File, time, Q) in results]
    print(f"Rendering {len(render_tasks)} plot(s) with {nworkers} worker(s)...")
    with Pool(nworkers) as pool:
        for i, _ in enumerate(pool.imap_unordered(_render_worker, render_tasks), 1):
            print(f"  render {i}/{len(render_tasks)}", end='\r', flush=True)
    print()

    print("... done.")
