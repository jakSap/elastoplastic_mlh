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
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import LogNorm
from multiprocessing import Pool

QUANTITIES = ["rho", "P", "u", "Sxx", "Sxy", "Syy"]
LABELS = {
    "rho": r"$\varrho$",
    "P":   r"$P$",
    "u":   r"$u$",
    "Sxx": r"$S_{xx}$",
    "Sxy": r"$S_{xy}$",
    "Syy": r"$S_{yy}$",
}


def get_domain_bounds(pos, openBorders):
    if openBorders:
        margin = 0.05 * max(pos[:, 0].ptp(), pos[:, 1].ptp())
        return (pos[:, 0].min() - margin, pos[:, 0].max() + margin,
                pos[:, 1].min() - margin, pos[:, 1].max() + margin)
    return (0., 1., 0., 1.)


def bin_particles(pos, h_grid, openBorders):
    """Return grid info and cell index per particle."""
    xmin, xmax, ymin, ymax = get_domain_bounds(pos, openBorders)
    nx = max(1, int(np.ceil((xmax - xmin) / h_grid)))
    ny = max(1, int(np.ceil((ymax - ymin) / h_grid)))
    # actual cell sizes (may be slightly different from h_grid at edges)
    dx = (xmax - xmin) / nx
    dy = (ymax - ymin) / ny

    ix = np.clip(((pos[:, 0] - xmin) / dx).astype(int), 0, nx - 1)
    iy = np.clip(((pos[:, 1] - ymin) / dy).astype(int), 0, ny - 1)
    cell = iy * nx + ix

    grid_info = dict(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax,
                     nx=nx, ny=ny, dx=dx, dy=dy)
    return grid_info, cell


def compute_grid_averages(vals, mass, cell, ncells):
    """Return (number-weighted, mass-weighted) 1-D arrays of length ncells."""
    count = np.bincount(cell, minlength=ncells).astype(float)
    val_sum = np.bincount(cell, weights=vals, minlength=ncells)
    mass_sum = np.bincount(cell, weights=mass, minlength=ncells)
    mass_val_sum = np.bincount(cell, weights=mass * vals, minlength=ncells)

    nw = np.full(ncells, np.nan)
    mw = np.full(ncells, np.nan)
    occupied = count > 0
    nw[occupied] = val_sum[occupied] / count[occupied]
    mw_denom = mass_sum.copy()
    mw_denom[mw_denom == 0] = np.nan
    mw[occupied] = mass_val_sum[occupied] / mw_denom[occupied]
    return nw, mw


def make_grid_meshes(grid_info):
    """Return X, Y meshes for pcolormesh (cell edges)."""
    g = grid_info
    x_edges = np.linspace(g['xmin'], g['xmax'], g['nx'] + 1)
    y_edges = np.linspace(g['ymin'], g['ymax'], g['ny'] + 1)
    return np.meshgrid(x_edges, y_edges)


def plot_grid(ax, X, Y, Z, grid_info, vmin=None, vmax=None, norm=None):
    """Plot a single pcolormesh on ax, return the mappable."""
    kw = dict(cmap='viridis')
    if norm is not None:
        kw['norm'] = norm
    else:
        if vmin is not None:
            kw['vmin'] = vmin
        if vmax is not None:
            kw['vmax'] = vmax
    Z2d = Z.reshape(grid_info['ny'], grid_info['nx'])
    pcm = ax.pcolormesh(X, Y, Z2d, **kw)
    ax.set_xlim(grid_info['xmin'], grid_info['xmax'])
    ax.set_ylim(grid_info['ymin'], grid_info['ymax'])
    ax.set_aspect('equal')
    return pcm


# ---------------------------------------------------------------------------
# Combined plot (2x3) — one figure per weighting mode
# ---------------------------------------------------------------------------
def createCombinedGridPlot(h5File, outDir, h_grid, openBorders, vminmax_n, vminmax_m,
                           diff=False, first_frame_data=None, diff_vminmax=None):
    data = h5.File(h5File, 'r')
    pos = data["x"][:]
    mass = data["m"][()]
    time = data["time"][0]

    qdata = {k: data[k][()] for k in QUANTITIES}
    data.close()

    grid_info, cell = bin_particles(pos, h_grid, openBorders)
    ncells = grid_info['nx'] * grid_info['ny']
    X, Y = make_grid_meshes(grid_info)

    show_diff = diff and first_frame_data is not None

    for mode, prefix, vm_dict in [("number", "gridN_comb", vminmax_n),
                                   ("mass",   "gridM_comb", vminmax_m)]:
        plt.rcParams.update({'font.size': 12})
        fig, axes = plt.subplots(2, 3, figsize=(10, 6), dpi=300)
        axes_flat = axes.flatten()

        if show_diff:
            DIFF_FLOOR = 1e-12
            for idx, key in enumerate(QUANTITIES):
                ax = axes_flat[idx]
                diff_vals = np.abs(qdata[key] - first_frame_data[key])
                diff_vals = np.clip(diff_vals, DIFF_FLOOR, None)
                nw, mw = compute_grid_averages(diff_vals, mass, cell, ncells)
                grid_vals = nw if mode == "number" else mw
                dm = diff_vminmax[key] if diff_vminmax and key in diff_vminmax \
                     else (DIFF_FLOOR, DIFF_FLOOR * 10.)
                norm = LogNorm(vmin=max(dm[0], DIFF_FLOOR),
                               vmax=max(dm[1], max(dm[0], DIFF_FLOOR) * 10.))
                pcm = plot_grid(ax, X, Y, grid_vals, grid_info, norm=norm)
                ax.set_title(r"$|\Delta$" + LABELS[key] + r"$|$")
                ax.set_xlabel("$x$"); ax.set_ylabel("$y$")
                divider = make_axes_locatable(ax)
                cax = divider.append_axes("right", size="5%", pad=0.05)
                fig.colorbar(pcm, cax=cax, format="%.1e")
        else:
            for idx, key in enumerate(QUANTITIES):
                ax = axes_flat[idx]
                nw, mw = compute_grid_averages(qdata[key], mass, cell, ncells)
                grid_vals = nw if mode == "number" else mw
                vm = vm_dict[key] if vm_dict and key in vm_dict else (None, None)
                vlo, vhi = vm
                if vlo is not None and vhi is not None and vlo == vhi:
                    vlo, vhi = vlo - 0.5, vhi + 0.5
                pcm = plot_grid(ax, X, Y, grid_vals, grid_info, vmin=vlo, vmax=vhi)
                ax.set_title(LABELS[key])
                ax.set_xlabel("$x$"); ax.set_ylabel("$y$")
                divider = make_axes_locatable(ax)
                cax = divider.append_axes("right", size="5%", pad=0.05)
                fig.colorbar(pcm, cax=cax)

        weight_label = "number-weighted" if mode == "number" else "mass-weighted"
        fig.suptitle(r"$t = " + f"{time:.4f}" + r"$" + f"  ({weight_label})", fontsize=16)
        plt.tight_layout()
        out = os.path.join(outDir, prefix + pathlib.Path(h5File).stem + ".png")
        print("Saving figure to", out)
        plt.savefig(out)
        plt.close()


# ---------------------------------------------------------------------------
# Single-quantity plot (side-by-side: number-weighted | mass-weighted)
# ---------------------------------------------------------------------------
def createSingleGridPlot(h5File, outDir, h_grid, key, openBorders,
                         vmin_n=None, vmax_n=None, vmin_m=None, vmax_m=None):
    data = h5.File(h5File, 'r')
    pos = data["x"][:]
    mass = data["m"][()]
    time = data["time"][0]
    vals = data[key][()]
    data.close()

    grid_info, cell = bin_particles(pos, h_grid, openBorders)
    ncells = grid_info['nx'] * grid_info['ny']
    X, Y = make_grid_meshes(grid_info)
    nw, mw = compute_grid_averages(vals, mass, cell, ncells)

    plt.rcParams.update({'font.size': 14})
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5), dpi=300)

    for ax, gv, vlo, vhi, wlabel in [
        (ax1, nw, vmin_n, vmax_n, "number-weighted"),
        (ax2, mw, vmin_m, vmax_m, "mass-weighted"),
    ]:
        pcm = plot_grid(ax, X, Y, gv, grid_info, vmin=vlo, vmax=vhi)
        ax.set_title(f"{LABELS[key]}  ({wlabel})")
        ax.set_xlabel("$x$"); ax.set_ylabel("$y$")
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        fig.colorbar(pcm, cax=cax)

    fig.suptitle(r"$t = " + f"{time:.4f}" + r"$", fontsize=16)
    plt.tight_layout()
    prefix = "grid_" + key + "_"
    out = os.path.join(outDir, prefix + pathlib.Path(h5File).stem + ".png")
    print("Saving figure to", out)
    plt.savefig(out)
    plt.close()


# ---------------------------------------------------------------------------
# Worker
# ---------------------------------------------------------------------------
def _worker_init(tmpdir):
    from matplotlib.texmanager import TexManager
    cache_dir = os.path.join(tmpdir, str(os.getpid()), "tex.cache")
    os.makedirs(cache_dir, exist_ok=True)
    TexManager._texcache = cache_dir


def _worker(task):
    (h5File, outDir, h_grid, openBorders, combined,
     key, vmin_n, vmax_n, vmin_m, vmax_m,
     combined_vminmax_n, combined_vminmax_m,
     diff, first_frame_data, diff_vminmax) = task
    if combined:
        createCombinedGridPlot(h5File, outDir, h_grid, openBorders,
                               combined_vminmax_n, combined_vminmax_m,
                               diff, first_frame_data, diff_vminmax)
    else:
        createSingleGridPlot(h5File, outDir, h_grid, key, openBorders,
                             vmin_n, vmax_n, vmin_m, vmax_m)


# ---------------------------------------------------------------------------
# Pre-scan helpers
# ---------------------------------------------------------------------------
def prescan_quantity_grid(files, key, h_grid, openBorders, mode):
    """Scan all files for global vmin/vmax of grid-averaged quantity."""
    lo, hi = None, None
    for f in files:
        with h5.File(f, 'r') as d:
            pos = d["x"][:]
            mass = d["m"][()]
            vals = d[key][()]
        grid_info, cell = bin_particles(pos, h_grid, openBorders)
        ncells = grid_info['nx'] * grid_info['ny']
        nw, mw = compute_grid_averages(vals, mass, cell, ncells)
        gv = nw if mode == "number" else mw
        gv = gv[np.isfinite(gv)]
        if len(gv) == 0:
            continue
        flo, fhi = float(gv.min()), float(gv.max())
        if lo is None or flo < lo: lo = flo
        if hi is None or fhi > hi: hi = fhi
    if lo is None or np.isnan(lo): lo = 0.
    if hi is None or np.isnan(hi): hi = 1.
    margin = 0.25 * (hi - lo) if hi > lo else 0.5
    return lo - margin, hi + margin


def get_output_prefix(args):
    if args.combined:   return "gridN_comb"
    elif args.pressure: return "grid_P_"
    elif args.energy:   return "grid_u_"
    elif args.noi:      return "grid_noi_"
    elif args.stress:   return "grid_Syy_"
    else:               return "grid_rho_"


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Grid-averaged plotting: partition domain into squares and plot weighted quantities.")
    parser.add_argument("--simOutputDir", "-d", metavar="str", type=str, required=True,
                        help="output directory of simulation")
    parser.add_argument("--outDir", "-o", metavar="str", type=str, default="output",
                        help="output directory for generated plots")
    parser.add_argument("--gridSize", "-k", metavar="float", type=float, required=True,
                        help="edge length of grid squares")
    parser.add_argument("--pressure", "-P", action="store_true", help="plot pressure")
    parser.add_argument("--energy", "-u", action="store_true", help="plot internal energy")
    parser.add_argument("--stress", "-S", action="store_true", help="plot Syy stress component")
    parser.add_argument("--noi", "-n", action="store_true", help="plot number of interactions")
    parser.add_argument("--combined", "-C", action="store_true",
                        help="plot all 6 quantities in 2x3 combined figures")
    parser.add_argument("--openBorders", "-b", action="store_true",
                        help="adjust domain to particle extent")
    parser.add_argument("--continue", "-c", dest="continue_", action="store_true",
                        help="skip h5 files whose plots already exist")
    parser.add_argument("--workers", "-w", metavar="int", type=int, default=os.cpu_count(),
                        help="number of parallel workers (default: all CPUs)")
    parser.add_argument("--diff", "-D", action="store_true",
                        help="(with --combined) show log-scale |diff to first frame|")
    parser.add_argument("--plotGhosts", "-G", action="store_true",
                        help="also plot ghost cells")

    args = parser.parse_args()

    if args.diff and not args.combined:
        print("WARNING: --diff has no effect without --combined.")

    plt.rc('text', usetex=True)

    print("Examining files in", args.simOutputDir, "...")
    h5Files = sorted([f for f in pathlib.Path(args.simOutputDir).glob('*.h5')
                      if "NNL" not in str(f) and (args.plotGhosts or "Ghost" not in str(f))])
    all_h5Files = list(h5Files)

    if args.continue_:
        prefix = get_output_prefix(args)
        existing = sorted(pathlib.Path(args.outDir).glob(f'{prefix}[0-9]*.png'))
        if existing:
            last_stem = existing[-1].stem[len(prefix):]
            h5Files = [f for f in h5Files if f.stem > last_stem]
            print(f"Continuing after {last_stem}, {len(h5Files)} file(s) remaining.")
        else:
            print("No existing plots found, starting from beginning.")

    # Determine quantity key for single-quantity mode
    if args.pressure:     colorKey = "P"
    elif args.energy:     colorKey = "u"
    elif args.noi:        colorKey = "noi"
    elif args.stress:     colorKey = "Syy"
    else:                 colorKey = "rho"

    # Pre-scan for global color limits
    vmin_n, vmax_n, vmin_m, vmax_m = None, None, None, None
    combined_vminmax_n, combined_vminmax_m = None, None
    first_frame_data = None
    diff_vminmax = None

    if len(h5Files) > 0:
        if args.combined:
            print(f"Pre-scanning {len(h5Files)} file(s) for grid color limits (all 6 quantities, both weightings)...")
            combined_vminmax_n = {}
            combined_vminmax_m = {}
            for key in QUANTITIES:
                for mode, vm_dict in [("number", combined_vminmax_n), ("mass", combined_vminmax_m)]:
                    lo, hi = prescan_quantity_grid(h5Files, key, args.gridSize, args.openBorders, mode)
                    vm_dict[key] = (lo, hi)
                    print(f"  {key} ({mode}): vmin={lo:.4g}, vmax={hi:.4g}")

            if args.diff and all_h5Files:
                DIFF_FLOOR = 1e-12
                first_frame_data = {}
                with h5.File(all_h5Files[0], 'r') as d:
                    for key in QUANTITIES:
                        first_frame_data[key] = d[key][()]
                print(f"First frame loaded from {all_h5Files[0].name} for diff computation.")

                diff_vminmax = {}
                print(f"Pre-scanning {len(all_h5Files)} file(s) for diff color limits...")
                for key in QUANTITIES:
                    d_lo, d_hi = None, None
                    for f in all_h5Files:
                        with h5.File(f, 'r') as d:
                            diff_vals = np.abs(d[key][()] - first_frame_data[key])
                        flo = float(np.clip(diff_vals, DIFF_FLOOR, None).min())
                        fhi = float(diff_vals.max())
                        if d_lo is None or flo < d_lo: d_lo = flo
                        if d_hi is None or fhi > d_hi: d_hi = fhi
                    if not d_lo or d_lo <= 0.: d_lo = DIFF_FLOOR
                    if not d_hi or d_hi <= d_lo: d_hi = d_lo * 10.
                    diff_vminmax[key] = (d_lo, d_hi)
                    print(f"  diff {key}: vmin={d_lo:.4g}, vmax={d_hi:.4g}")
        else:
            print(f"Pre-scanning {len(h5Files)} file(s) for grid color limits ({colorKey})...")
            vmin_n, vmax_n = prescan_quantity_grid(h5Files, colorKey, args.gridSize, args.openBorders, "number")
            vmin_m, vmax_m = prescan_quantity_grid(h5Files, colorKey, args.gridSize, args.openBorders, "mass")
            print(f"  number-weighted: vmin={vmin_n:.4g}, vmax={vmax_n:.4g}")
            print(f"  mass-weighted:   vmin={vmin_m:.4g}, vmax={vmax_m:.4g}")

    tasks = [(str(f), args.outDir, args.gridSize, args.openBorders, args.combined,
              colorKey, vmin_n, vmax_n, vmin_m, vmax_m,
              combined_vminmax_n, combined_vminmax_m,
              args.diff, first_frame_data, diff_vminmax)
             for f in h5Files]

    nworkers = min(args.workers, len(tasks)) if tasks else 1
    print(f"Rendering {len(tasks)} plot(s) with {nworkers} worker(s)...")
    tmpdir = tempfile.mkdtemp(prefix="mpl_grid_workers_")
    try:
        with Pool(nworkers, initializer=_worker_init, initargs=(tmpdir,)) as pool:
            for i, _ in enumerate(pool.imap_unordered(_worker, tasks), 1):
                print(f"  {i}/{len(tasks)} done", end='\r', flush=True)
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)
    print()
    print("... done.")
