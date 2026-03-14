#!/usr/bin/env python3

import argparse
import os
import pathlib
import tempfile
import shutil
import numpy as np
import h5py as h5
import matplotlib
matplotlib.use('Agg')  # non-interactive backend, safe for multiprocessing
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import LogNorm
from multiprocessing import Pool

#MAX_NUM_INTERACTIONS = 1000
pref = 'S'
def setDomainLimits(ax, pos, h5File, openBorders):
    if openBorders:
        margin = 0.05 * max(pos[:,0].max() - pos[:,0].min(), pos[:,1].max() - pos[:,1].min())
        # if (np.isnan(pos[:,0].min()) || (pos[:,0].max()) np.isnan)
        ax.set_xlim((pos[:,0].min() - margin, pos[:,0].max() + margin))
        ax.set_ylim((pos[:,1].min() - margin, pos[:,1].max() + margin))
    elif "Ghosts" not in str(h5File):
        ax.set_xlim((0., 1.))
        ax.set_ylim((0., 1.))

def createPlot(h5File, outDir, plotGrad, plotVel, stress, iNNL, openBorders=False, vmin=None, vmax=None, markerSize=1.):
    data = h5.File(h5File, 'r')
    pos = data["x"][:]

    time = data["time"][0]
    if stress:
        Sxx = np.array(data["Sxx"][:])
        Sxy = np.array(data["Sxy"][:])
        Syy = np.array(data["Syy"][:])
        SAbs = np.sqrt(Sxx**2 + Sxy**2 + Syy**2)
        cm = Syy
    else:
        rho = data["rho"][()]
        cm = rho
    #P = data["P"][()]
    plt.rcParams.update({'font.size': 18})
    fig, ax = plt.subplots(figsize=(7,5), dpi=500)
    #rhoPlt = ax.scatter(pos[:,0], pos[:,1], c=rho, s=500.) # good for ~100 particles
    #rhoPlt = ax.scatter(pos[:,0], pos[:,1], c=rho, s=150.) # good for ~400 particles
    #rhoPlt = ax.scatter(pos[:,0], pos[:,1], c=rho, s=100.) # good for ~900 particles
    # rhoPlt = ax.scatter(pos[:,0], pos[:,1], c=rho, s=12.) # good for 10**4 particles
    # rhoPlt = ax.scatter(pos[:,0], pos[:,1], c=rho, s=5.) # good for 128**2 particles
    rhoPlt = ax.scatter(pos[:,0], pos[:,1], c=cm, s=markerSize)#, vmin=vmin, vmax=vmax) # good for 128**2 particles

    setDomainLimits(ax, pos, h5File, openBorders)

    # Plot gradient
    if plotGrad and not plotVel:
        plotGradient(data["rhoGrad"][:], pos, ax)
    elif not plotGrad and plotVel:
        plotVelocity(data["v"][:], pos, ax)
    elif plotGrad and plotVel:
        print("WARNING: command line arguments '--plotVelocity' and '--plotGradient' are incompatible. - Plotting neither.")

    # plot NNL for particle i
    if iNNL > -1 and "Ghosts" not in str(h5File):
        plotNNL(h5File, iNNL, pos, ax)

    if stress:
        plt.title(r"Color coded " + pref + r" at $t = " + f"{time:.4f}" + r"$")
    else:
        plt.title(r"Color coded density $\varrho$ at $t = " + f"{time:.4f}" + r"$")
    #plt.title(r"Color coded pressure $P$")
    plt.xlabel("$x$")
    plt.ylabel("$y$")

    # create an axes on the right side of ax. The width of cax will be 5%
    # of ax and the padding between cax and ax will be fixed at 0.05 inch.
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    fig.colorbar(rhoPlt, cax=cax) #, orientation='horizontal')
    #fig.colorbar(PPlt, ax=ax)
    
    ax.set_aspect('equal')
    
    plt.tight_layout()
    print("Saving figure to", outDir + "/" + pref + pathlib.Path(h5File).stem + ".png")
    if stress:
        plt.savefig(outDir +'/'+ pref + pathlib.Path(h5File).stem + ".png")
    else:
        plt.savefig(outDir + "/" + pathlib.Path(h5File).stem + ".png")    
    plt.close()
    #plt.show()

def plotGradient(grad, pos, ax):
    ax.quiver(pos[:,0], pos[:,1], grad[:,0], grad[:,1], angles='xy', scale_units='xy', scale=1.)
    #ax.quiver(pos[:,0], pos[:,1], grad[:,0], grad[:,1], angles='xy', scale=.01)
    
    #for i, rhoGrad in enumerate(grad):
    #    if np.linalg.norm(rhoGrad) > .05:
    #        print("rhoGrad @", i, "=", rhoGrad)

def plotVelocity(vel, pos, ax):
    ax.quiver(pos[:,0], pos[:,1], vel[:,0], vel[:,1], angles='xy', scale_units='xy', scale=.5)

def get_output_prefix(args):
    if args.combined:   return "comb"
    elif args.pressure: return "P"
    elif args.energy:   return "u"
    elif args.noi:      return "noi"
    elif args.stress:   return pref
    else:               return ""

def plotNNL(h5File, iNNL, pos, ax):
    data = h5.File(str(h5File).replace(".h5", "NNL.h5"), 'r')
    posNNL = data["nnlPrtcls"+str(iNNL)][:]
    ax.scatter(pos[iNNL,0], pos[iNNL,1], s=5., marker='x', color='r')
    ax.scatter(posNNL[:,0], posNNL[:,1], s=5, marker='x', color='m')

def createEnergyPlot(h5File, outDir, openBorders=False, vmin=None, vmax=None, markerSize=1.):
    data = h5.File(h5File, 'r')
    pos = data["x"][:]

    u = data["u"][()]
    fig, ax = plt.subplots(figsize=(8,6), dpi=200)
    #uPlt = ax.scatter(pos[:,0], pos[:,1], c=u, s=100.) # good for ~900 particles
    #uPlt = ax.scatter(pos[:,0], pos[:,1], c=u, s=200.) # good for ~400 particles
    uPlt = ax.scatter(pos[:,0], pos[:,1], c=u, s=markerSize, vmin=vmin, vmax=vmax)

    setDomainLimits(ax, pos, h5File, openBorders)

    fig.colorbar(uPlt, ax=ax)
    plt.title(r"Color coded internal energy $u$")
    plt.xlabel("$x$")
    plt.ylabel("$y$")
    plt.tight_layout()
    print("Saving figure to", outDir + "/u" + pathlib.Path(h5File).stem + ".png")
    plt.savefig(outDir + "/u" + pathlib.Path(h5File).stem + ".png")
    plt.close()

def createPressurePlot(h5File, outDir, openBorders=False, vmin=None, vmax=None, markerSize=1.):
    data = h5.File(h5File, 'r')
    pos = data["x"][:]

    time = data["time"][0]

    P = data["P"][()]

    plt.rcParams.update({'font.size': 18})
    fig, ax = plt.subplots(figsize=(7,6), dpi=200)
    #fig, ax = plt.subplots(figsize=(8,6), dpi=200)
    #PPlt = ax.scatter(pos[:,0], pos[:,1], c=P, s=200.) # good for ~400 particles
    #PPlt = ax.scatter(pos[:,0], pos[:,1], c=P, s=100.) # good for ~900 particles
    #PPlt = ax.scatter(pos[:,0], pos[:,1], c=P, s=10.) # good for 10**4 particles
    PPlt = ax.scatter(pos[:,0], pos[:,1], c=P, s=markerSize, vmin=vmin, vmax=vmax)

    setDomainLimits(ax, pos, h5File, openBorders)

    plt.title(r"Color coded pressure $P$ at $t = " + f"{time:.4f}" + r"$")
    plt.xlabel("$x$")
    plt.ylabel("$y$")    

    # create an axes on the right side of ax. The width of cax will be 5%
    # of ax and the padding between cax and ax will be fixed at 0.05 inch.
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)  

    fig.colorbar(PPlt, cax=cax)

    ax.set_aspect('equal')
    
    plt.tight_layout()
    print("Saving figure to", outDir + "/P" + pathlib.Path(h5File).stem + ".png")
    plt.savefig(outDir + "/P" + pathlib.Path(h5File).stem + ".png")
    plt.close()

def createNoiPlot(h5File, outDir, openBorders=False, vmin=None, vmax=None, markerSize=1.):
    data = h5.File(h5File, 'r')
    pos = data["x"][:]

    noi = data["noi"][()]
    fig, ax = plt.subplots(figsize=(8,6), dpi=200)
    noiPlt = ax.scatter(pos[:,0], pos[:,1], c=noi, s=markerSize, vmin=vmin, vmax=vmax)

    setDomainLimits(ax, pos, h5File, openBorders)
    fig.colorbar(noiPlt, ax=ax)
    plt.title(r"Color coded number of interactions")
    plt.xlabel("$x$")
    plt.ylabel("$y$")
    plt.tight_layout()
    print("Saving figure to", outDir + "/noi" + pathlib.Path(h5File).stem + ".png")
    plt.savefig(outDir + "/noi" + pathlib.Path(h5File).stem + ".png")
    plt.close()

def createCombinedPlot(h5File, outDir, openBorders=False, vminmax=None, markerSize=1.,
                       diff=False, first_frame_data=None, diff_vminmax=None):
    data = h5.File(h5File, 'r')
    pos  = data["x"][:]
    time = data["time"][0]

    quantities = [
        ("rho", data["rho"][()],  r"Density $\varrho$"),
        ("P",   data["P"][()],    r"Pressure $P$"),
        ("u",   data["u"][()],    r"Internal energy $u$"),
        ("Sxx", data["Sxx"][:],   r"$S_{xx}$"),
        ("Sxy", data["Sxy"][:],   r"$S_{xy}$"),
        ("Syy", data["Syy"][:],   r"$S_{yy}$"),
    ]
    quantities_dict = {key: vals for key, vals, _ in quantities}

    plt.rcParams.update({'font.size': 12})
    show_diff = diff and first_frame_data is not None
    fig, axes = plt.subplots(2, 3, figsize=(10, 6), dpi=300)
    axes_flat = axes.flatten()

    if show_diff:
        DIFF_FLOOR = 1e-12
        diff_meta = [
            (r"$|\Delta\varrho|$", "rho"),
            (r"$|\Delta P|$",      "P"),
            (r"$|\Delta u|$",      "u"),
            (r"$|\Delta S_{xx}|$", "Sxx"),
            (r"$|\Delta S_{xy}|$", "Sxy"),
            (r"$|\Delta S_{yy}|$", "Syy"),
        ]
        for idx, (dlabel, key) in enumerate(diff_meta):
            ax = axes_flat[idx]
            diff_vals = np.abs(quantities_dict[key] - first_frame_data[key])
            diff_vals = np.clip(diff_vals, DIFF_FLOOR, None)
            dm = diff_vminmax[key] if diff_vminmax and key in diff_vminmax \
                 else (DIFF_FLOOR, DIFF_FLOOR * 10.)
            norm = LogNorm(vmin=max(dm[0], DIFF_FLOOR),
                           vmax=max(dm[1], max(dm[0], DIFF_FLOOR) * 10.))
            sc = ax.scatter(pos[:,0], pos[:,1], c=diff_vals, s=markerSize,
                            norm=norm, cmap='viridis')
            setDomainLimits(ax, pos, h5File, openBorders)
            ax.set_aspect('equal')
            ax.set_title(dlabel)
            ax.set_xlabel("$x$")
            ax.set_ylabel("$y$")
            divider = make_axes_locatable(ax)
            cax = divider.append_axes("right", size="5%", pad=0.05)
            fig.colorbar(sc, cax=cax, format="%.1e")
    else:
        for ax, (key, vals, label) in zip(axes_flat, quantities):
            vm = vminmax[key] if vminmax and key in vminmax else (None, None)

            vlo, vhi = vm
            if vlo is not None and vhi is not None and vlo == vhi:
                vlo, vhi = vlo - 0.5, vhi + 0.5
            # Replace NaN values with 0 to avoid matplotlib StopIteration error
            vals = np.where(np.isfinite(vals), vals, 0.)
            sc = ax.scatter(pos[:,0], pos[:,1], c=vals, s=markerSize, vmin=vlo, vmax=vhi)

            setDomainLimits(ax, pos, h5File, openBorders)
            ax.set_aspect('equal')
            ax.set_title(label)
            ax.set_xlabel("$x$")
            ax.set_ylabel("$y$")
            divider = make_axes_locatable(ax)
            cax = divider.append_axes("right", size="5%", pad=0.05)
            fig.colorbar(sc, cax=cax)

    fig.suptitle(r"$t = " + f"{time:.4f}" + r"$", fontsize=16)
    plt.tight_layout()
    out = outDir + "/comb" + pathlib.Path(h5File).stem + ".png"
    print("Saving figure to", out)
    plt.savefig(out)
    plt.close()


def _worker_init(tmpdir):
    """Give each pool worker its own TeX cache directory so concurrent LaTeX
    compilations don't corrupt each other's cache files."""
    from matplotlib.texmanager import TexManager
    cache_dir = os.path.join(tmpdir, str(os.getpid()), "tex.cache")
    os.makedirs(cache_dir, exist_ok=True)
    TexManager._texcache = cache_dir


def _worker(task):
    (h5File, outDir, grad, vel, stress, iNNL, borders, vmin, vmax,
     pressure, energy, noi, combined, combined_vminmax,
     markerSize, diff, first_frame_data, diff_vminmax) = task
    if combined:
        createCombinedPlot(h5File, outDir, borders, combined_vminmax, markerSize,
                           diff, first_frame_data, diff_vminmax)
    elif pressure:
        createPressurePlot(h5File, outDir, borders, vmin, vmax, markerSize)
    elif energy:
        createEnergyPlot(h5File, outDir, borders, vmin, vmax, markerSize)
    elif noi:
        createNoiPlot(h5File, outDir, borders, vmin, vmax, markerSize)
    else:
        createPlot(h5File, outDir, grad, vel, stress, iNNL, borders, vmin, vmax, markerSize)


if __name__=="__main__":
    parser = argparse.ArgumentParser(description="Plot density of results from Kelvin-Helmholtz test case.")
    parser.add_argument("--simOutputDir", "-d", metavar="string", type=str, help="output directory of simulation", required=True)
    parser.add_argument("--outDir", "-o", metavar="string", type=str, help="output directory for generated plots", default="output")
    parser.add_argument("--plotGradient", "-g", action="store_true", help="plot density gradients")
    parser.add_argument("--plotGhosts", "-G", action="store_true", help="also plot ghost cells in an extra file")
    parser.add_argument("--pressure", "-P", action="store_true", help="plot pressure instead of density")
    parser.add_argument("--energy", "-u", action="store_true", help="plot internal energy instead of density")
    parser.add_argument("--stress", "-S", action="store_true", help="Plot max value of stress")
    parser.add_argument("--combined", "-C", action="store_true",
                        help="Plot rho, P, u, Sxx, Sxy, Syy in a 2x3 combined figure")
    parser.add_argument("--noi", "-n", action="store_true", help="plot number of interactions instead of density")
    parser.add_argument("--plotVelocity", "-v", action="store_true", help="plot velocity")
    parser.add_argument("--iNNL", "-i", metavar="int", type=int, help="plot NNL for particles i", default=-1)
    parser.add_argument("--openBorders", "-b", action="store_true", help="Adjust plot domain to show all real particles")
    parser.add_argument("--continue", "-c", dest="continue_", action="store_true",
                        help="skip h5 files whose plots already exist in outDir")
    parser.add_argument("--workers", "-w", metavar="int", type=int, default=os.cpu_count(),
                        help="number of parallel worker processes (default: all CPUs)")
    parser.add_argument("--markerSize", "-m", metavar="float", type=float, default=1.,
                        help="scatter marker size for all plots (default: 1.0)")
    parser.add_argument("--diff", "-D", action="store_true",
                        help="(with --combined) add two rows showing log-scale |diff to first frame|")

    args = parser.parse_args()

    if args.diff and not args.combined:
        print("WARNING: --diff has no effect without --combined.")

    plt.rc('text', usetex=True)

    print("Examining files in", args.simOutputDir, "...")

    h5Files = sorted([f for f in pathlib.Path(args.simOutputDir).glob('*.h5')
                       if "NNL" not in str(f) and (args.plotGhosts or "Ghost" not in str(f))])
    all_h5Files = list(h5Files)  # full sorted list, before --continue filtering

    if args.continue_:
        prefix = get_output_prefix(args)
        existing = sorted(pathlib.Path(args.outDir).glob(f'{prefix}[0-9]*.png'))
        if existing:
            last_stem = existing[-1].stem[len(prefix):]
            h5Files = [f for f in h5Files if f.stem > last_stem]
            print(f"Continuing after {last_stem}, {len(h5Files)} file(s) remaining.")
        else:
            print("No existing plots found in outDir, starting from the beginning.")

    def prescan_quantity(files, key):
        lo, hi = None, None
        for f in files:
            with h5.File(f, 'r') as d:
                vals = d[key][()]
            flo, fhi = float(vals.min()), float(vals.max())
            if lo is None or flo < lo: lo = flo
            if hi is None or fhi > hi: hi = fhi
        if lo is None or np.isnan(lo): lo = 0.
        if hi is None or np.isnan(hi): hi = 1.
        margin = 0.25 * (hi - lo) if hi > lo else 0.5
        return lo - margin, hi + margin

    # Pre-scan all files for global color limits (full pass, 25% margin)
    vmin, vmax = None, None
    combined_vminmax = None
    first_frame_data = None
    diff_vminmax = None
    if len(h5Files) > 0:
        if args.combined:
            print(f"Pre-scanning {len(h5Files)} file(s) for color limits (all 6 quantities)...")
            combined_vminmax = {}
            for key in ["rho", "P", "u", "Sxx", "Sxy", "Syy"]:
                lo, hi = prescan_quantity(h5Files, key)
                combined_vminmax[key] = (lo, hi)
                print(f"  {key}: vmin={lo:.4g}, vmax={hi:.4g}")

            if args.diff and all_h5Files:
                DIFF_FLOOR = 1e-12
                first_frame_data = {}
                with h5.File(all_h5Files[0], 'r') as d:
                    for key in ["rho", "P", "u", "Sxx", "Sxy", "Syy"]:
                        first_frame_data[key] = d[key][()]
                print(f"First frame loaded from {all_h5Files[0].name} for diff computation.")

                diff_vminmax = {}
                print(f"Pre-scanning {len(all_h5Files)} file(s) for diff color limits...")
                for key in ["rho", "P", "u", "Sxx", "Sxy", "Syy"]:
                    d_lo, d_hi = None, None
                    for f in all_h5Files:
                        with h5.File(f, 'r') as d:
                            diff = np.abs(d[key][()] - first_frame_data[key])
                        flo = float(np.clip(diff, DIFF_FLOOR, None).min())
                        fhi = float(diff.max())
                        if d_lo is None or flo < d_lo: d_lo = flo
                        if d_hi is None or fhi > d_hi: d_hi = fhi
                    if not d_lo or d_lo <= 0.: d_lo = DIFF_FLOOR
                    if not d_hi or d_hi <= d_lo: d_hi = d_lo * 10.
                    diff_vminmax[key] = (d_lo, d_hi)
                    print(f"  diff {key}: vmin={d_lo:.4g}, vmax={d_hi:.4g}")
        else:
            # Helper to get the plotted quantity key
            if args.pressure:
                colorKey = "P"
            elif args.energy:
                colorKey = "u"
            elif args.noi:
                colorKey = "noi"
            else:
                colorKey = "rho"
            print(f"Pre-scanning {len(h5Files)} file(s) for color limits...")
            vmin, vmax = prescan_quantity(h5Files, colorKey)
            print(f"Global color limits: vmin={vmin:.4g}, vmax={vmax:.4g} (25% margin)")

    if args.plotGradient or args.iNNL > -1:
        if args.pressure or args.energy or args.noi:
            print("WARNING: '--plotGradient' and '--iNNL' are ignored for this plot mode.")

    tasks = [(f, args.outDir, args.plotGradient, args.plotVelocity, args.stress,
              args.iNNL, args.openBorders, vmin, vmax,
              args.pressure, args.energy, args.noi, args.combined, combined_vminmax,
              args.markerSize, args.diff, first_frame_data, diff_vminmax)
             for f in h5Files]

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
    
