#!/usr/bin/env python3

import argparse
import os
import pathlib
import numpy as np
import h5py as h5
import matplotlib
matplotlib.use('Agg')  # non-interactive backend, safe for multiprocessing
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import LogNorm
from matplotlib.patches import Circle
from multiprocessing import Pool
import tempfile


def _init_tex_worker():
    """Give each worker process its own matplotlib TeX cache directory.

    With usetex=True, all workers otherwise share ~/.cache/matplotlib/tex.cache;
    concurrent processes writing the same cached .png clobber each other, and a
    reader then hits a half-written file ("SyntaxError: not a PNG file"). A
    per-process cache removes the contention entirely."""
    from matplotlib import texmanager
    d = tempfile.mkdtemp(prefix=f'mpltex_{os.getpid()}_')
    texmanager.TexManager.texcache = d

#MAX_NUM_INTERACTIONS = 1000
pref = 'S'

# Envelope (matId=2) particles are rendered at 1/3 the marker size of the main
# material so the solid ring stays visually dominant.
ENVELOPE_MATID = 2
ENVELOPE_SCALE = 1.0 / 3.0


def is_gizmo_file(path):
    """GADGET-style snapshot (PartType0/...) vs the demonstrator's flat layout.
    We use the .hdf5 extension as the cheap discriminator; demonstrator writes .h5."""
    return str(path).endswith('.hdf5')


def open_snapshot(path):
    """Open a snapshot in either format. For demonstrator files, returns the
    h5py.File (so all existing data['key'][...] calls work). For GIZMO files,
    returns a dict-of-ndarrays that supports the same access pattern under the
    demonstrator's key names. GIZMO snapshots lack stress / gradient / NNL
    fields, so plot modes that need those will fall back to a clear KeyError."""
    if is_gizmo_file(path):
        with h5.File(path, 'r') as f:
            coords = f['PartType0/Coordinates'][:]
            out = {
                'x'   : coords[:, :2],
                'rho' : f['PartType0/Density'][:],
                'm'   : f['PartType0/Masses'][:],
                'time': np.array([float(f['Header'].attrs['Time'])]),
            }
            if 'PartType0/Velocities'      in f: out['v']   = f['PartType0/Velocities'][:, :2]
            if 'PartType0/Pressure'        in f: out['P']   = f['PartType0/Pressure'][:]
            if 'PartType0/InternalEnergy'  in f: out['u']   = f['PartType0/InternalEnergy'][:]
            if 'PartType0/SmoothingLength' in f: out['sml'] = f['PartType0/SmoothingLength'][:]
        return out
    return h5.File(path, 'r')

def _marker_sizes(data, markerSize):
    """Return a per-particle size array if a materialId field is present,
    scaling envelope particles down. Otherwise return the scalar markerSize."""
    if "materialId" not in data:
        return markerSize
    matId = data["materialId"][:]
    sizes = np.full(matId.shape, float(markerSize))
    sizes[matId == ENVELOPE_MATID] = float(markerSize) * ENVELOPE_SCALE
    return sizes
def setDomainLimits(ax, pos, h5File, openBorders):
    if openBorders:
        margin = 0.05 * max(pos[:,0].max() - pos[:,0].min(), pos[:,1].max() - pos[:,1].min())
        if (np.isnan(pos[:,0].min())):
            # | (np.isnan(pos[:,0].max()))):
            ax.set_xlim((-8, 8))
            ax.set_ylim((-4.5, 4.5))
        else:
            ax.set_xlim((pos[:,0].min() - margin, pos[:,0].max() + margin))
            ax.set_ylim((pos[:,1].min() - margin, pos[:,1].max() + margin))
    elif "Ghosts" not in str(h5File):
        ax.set_xlim((0., 1.))
        ax.set_ylim((0., 1.))

def createPlot(h5File, outDir, plotGrad, plotVel, stress, iNNL, openBorders=False, vmin=None, vmax=None, markerSize=1., iHi=-1, dpi=200):
    data = open_snapshot(h5File)
    gizmo = is_gizmo_file(h5File)
    pos = data["x"][:]

    time = data["time"][0]
    if stress:
        if gizmo:
            raise KeyError("--stress is unavailable for GIZMO snapshots (no Sxx/Sxy/Syy fields)")
        Sxx = np.array(data["Sxx"][:])
        Sxy = np.array(data["Sxy"][:])
        Syy = np.array(data["Syy"][:])
        SAbs = np.sqrt(Sxx**2 + Sxy**2 + Syy**2)
        cm = SAbs
    else:
        rho = data["rho"][()]
        cm = rho
    #P = data["P"][()]
    plt.rcParams.update({'font.size': 18})
    fig, ax = plt.subplots(figsize=(7,5), dpi=dpi)
    #rhoPlt = ax.scatter(pos[:,0], pos[:,1], c=rho, s=500.) # good for ~100 particles
    #rhoPlt = ax.scatter(pos[:,0], pos[:,1], c=rho, s=150.) # good for ~400 particles
    #rhoPlt = ax.scatter(pos[:,0], pos[:,1], c=rho, s=100.) # good for ~900 particles
    # rhoPlt = ax.scatter(pos[:,0], pos[:,1], c=rho, s=12.) # good for 10**4 particles
    # rhoPlt = ax.scatter(pos[:,0], pos[:,1], c=rho, s=5.) # good for 128**2 particles
    rhoPlt = ax.scatter(pos[:,0], pos[:,1], c=cm, s=_marker_sizes(data, markerSize))#, vmin=vmin, vmax=vmax) # good for 128**2 particles

    setDomainLimits(ax, pos, h5File, openBorders)

    # Plot gradient
    if plotGrad and not plotVel:
        if gizmo:
            print(f"WARNING: --plotGradient unavailable for GIZMO snapshot {h5File}, skipping.")
        else:
            plotGradient(data["rhoGrad"][:], pos, ax)
    elif not plotGrad and plotVel:
        plotVelocity(data["v"][:], pos, ax)
    elif plotGrad and plotVel:
        print("WARNING: command line arguments '--plotVelocity' and '--plotGradient' are incompatible. - Plotting neither.")

    # plot NNL for particle i
    if iNNL > -1 and "Ghosts" not in str(h5File):
        if gizmo:
            print(f"WARNING: --iNNL unavailable for GIZMO snapshot {h5File}, skipping.")
        else:
            plotNNL(h5File, iNNL, pos, ax)

    # plot kernel circle for particle i
    if iHi > -1 and "Ghosts" not in str(h5File):
        plotKernelCircle(data, iHi, pos, ax)

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
        plt.savefig(outDir + "/" + pathlib.Path(h5File).stem + ".png", dpi=dpi)
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
    if args.combined:          return "comb"
    elif args.pressure:        return "P"
    elif args.energy:          return "u"
    elif args.noi:             return "noi"
    elif args.conditionNumber: return "ncond"
    elif args.stress:          return pref
    else:                      return ""

def plotKernelCircle(data, iHi, pos, ax):
    sml = data["sml"][:]
    h = float(sml[iHi])
    h_avg = float(np.mean(sml))
    print(f"  Particle {iHi}: h = {h:.6g}  |  avg h (all particles) = {h_avg:.6g}")
    ax.scatter(pos[iHi, 0], pos[iHi, 1], s=5., marker='+', color='red', zorder=5)
    circle = Circle((pos[iHi, 0], pos[iHi, 1]), h, fill=False,
                    edgecolor='red', linewidth=0.5, zorder=5)
    ax.add_patch(circle)
    ax.annotate(f"$h_i={h:.4g}$\n$\\bar{{h}}={h_avg:.4g}$",
                xy=(pos[iHi, 0], pos[iHi, 1]),
                xytext=(-40, 10), textcoords='offset points',
                fontsize=6, color='red', zorder=6)

def plotNNL(h5File, iNNL, pos, ax):
    data = h5.File(str(h5File).replace(".h5", "NNL.h5"), 'r')
    posNNL = data["nnlPrtcls"+str(iNNL)][:]
    ax.scatter(pos[iNNL,0], pos[iNNL,1], s=5., marker='x', color='r')
    ax.scatter(posNNL[:,0], posNNL[:,1], s=5, marker='x', color='m')

def createEnergyPlot(h5File, outDir, openBorders=False, vmin=None, vmax=None, markerSize=1., iHi=-1, dpi=200):
    data = open_snapshot(h5File)
    if "u" not in data:
        print(f"WARNING: internal energy unavailable for {h5File} (no u field), skipping.")
        return
    pos = data["x"][:]

    u = data["u"][()]
    fig, ax = plt.subplots(figsize=(8,6), dpi=dpi)
    #uPlt = ax.scatter(pos[:,0], pos[:,1], c=u, s=100.) # good for ~900 particles
    #uPlt = ax.scatter(pos[:,0], pos[:,1], c=u, s=200.) # good for ~400 particles
    uPlt = ax.scatter(pos[:,0], pos[:,1], c=u, s=_marker_sizes(data, markerSize), vmin=vmin, vmax=vmax)

    setDomainLimits(ax, pos, h5File, openBorders)

    if iHi > -1:
        plotKernelCircle(data, iHi, pos, ax)

    fig.colorbar(uPlt, ax=ax)
    plt.title(r"Color coded internal energy $u$")
    plt.xlabel("$x$")
    plt.ylabel("$y$")
    plt.tight_layout()
    print("Saving figure to", outDir + "/u" + pathlib.Path(h5File).stem + ".png")
    plt.savefig(outDir + "/u" + pathlib.Path(h5File).stem + ".png")
    plt.close()

def createPressurePlot(h5File, outDir, openBorders=False, vmin=None, vmax=None, markerSize=1., iHi=-1, dpi=200):
    data = open_snapshot(h5File)
    if "P" not in data:
        print(f"WARNING: pressure unavailable for {h5File} (GIZMO snapshots carry no Pressure field), skipping.")
        return
    pos = data["x"][:]

    time = data["time"][0]

    P = data["P"][()]

    plt.rcParams.update({'font.size': 18})
    fig, ax = plt.subplots(figsize=(7,6), dpi=dpi)
    #fig, ax = plt.subplots(figsize=(8,6), dpi=200)
    #PPlt = ax.scatter(pos[:,0], pos[:,1], c=P, s=200.) # good for ~400 particles
    #PPlt = ax.scatter(pos[:,0], pos[:,1], c=P, s=100.) # good for ~900 particles
    #PPlt = ax.scatter(pos[:,0], pos[:,1], c=P, s=10.) # good for 10**4 particles
    PPlt = ax.scatter(pos[:,0], pos[:,1], c=P, s=_marker_sizes(data, markerSize), vmin=vmin, vmax=vmax)

    setDomainLimits(ax, pos, h5File, openBorders)

    if iHi > -1:
        plotKernelCircle(data, iHi, pos, ax)

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

def createNoiPlot(h5File, outDir, openBorders=False, vmin=None, vmax=None, markerSize=1., iHi=-1, dpi=200):
    data = open_snapshot(h5File)
    if "noi" not in data:
        print(f"WARNING: number of interactions unavailable for {h5File} (no noi field), skipping.")
        return
    pos = data["x"][:]

    noi = data["noi"][()]
    fig, ax = plt.subplots(figsize=(8,6), dpi=dpi)
    noiPlt = ax.scatter(pos[:,0], pos[:,1], c=noi, s=_marker_sizes(data, markerSize), vmin=vmin, vmax=vmax)

    setDomainLimits(ax, pos, h5File, openBorders)

    if iHi > -1:
        plotKernelCircle(data, iHi, pos, ax)

    fig.colorbar(noiPlt, ax=ax)
    plt.title(r"Color coded number of interactions")
    plt.xlabel("$x$")
    plt.ylabel("$y$")
    plt.tight_layout()
    print("Saving figure to", outDir + "/noi" + pathlib.Path(h5File).stem + ".png")
    plt.savefig(outDir + "/noi" + pathlib.Path(h5File).stem + ".png")
    plt.close()

def createConditionNumberPlot(h5File, outDir, openBorders=False, vmin=None, vmax=None, markerSize=1., threshold=None, iHi=-1, dpi=200):
    data = open_snapshot(h5File)
    if "conditionNumber" not in data:
        print(f"WARNING: condition number unavailable for {h5File} (no conditionNumber field), skipping.")
        return
    pos = data["x"][:]
    time = data["time"][0]
    ncond = data["conditionNumber"][()]

    plt.rcParams.update({'font.size': 18})
    fig, ax = plt.subplots(figsize=(7,6), dpi=dpi)
    ncondPlt = ax.scatter(pos[:,0], pos[:,1], c=ncond, s=_marker_sizes(data, markerSize), vmin=vmin, vmax=vmax)

    if threshold is not None:
        mask = ncond > threshold[0]
        if "eigenvecMin" in data:
            eigvec_thresh = data["eigenvecMin"][:]
            mask &= np.abs(eigvec_thresh[:,0]) < threshold[1]
            mask &= np.abs(eigvec_thresh[:,1]) < threshold[2]
        if np.any(mask):
            ax.scatter(pos[mask,0], pos[mask,1], s=markerSize*4,
                       facecolors='none', edgecolors='r', linewidths=0.5)

    setDomainLimits(ax, pos, h5File, openBorders)

    if iHi > -1:
        plotKernelCircle(data, iHi, pos, ax)

    title = r"Condition number $N_{\mathrm{cond}}$ at $t = " + f"{time:.4f}" + r"$"
    if threshold is not None:
        title += (f"\n(red: $N_{{\\mathrm{{cond}}}} > {threshold[0]}$, "
                  f"$|e_x| < {threshold[1]}$, $|e_y| < {threshold[2]}$, "
                  f"count: {np.sum(mask)})")
    plt.title(title)
    plt.xlabel("$x$")
    plt.ylabel("$y$")

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    fig.colorbar(ncondPlt, cax=cax)

    ax.set_aspect('equal')

    plt.tight_layout()
    print("Saving figure to", outDir + "/ncond" + pathlib.Path(h5File).stem + ".png")
    plt.savefig(outDir + "/ncond" + pathlib.Path(h5File).stem + ".png")
    plt.close()

    # Histogram of condition numbers and eigenvector components
    has_eigvec = "eigenvecMin" in data
    nrows = 2 if has_eigvec else 1
    fig, axes = plt.subplots(nrows, 1, figsize=(7, 5 * nrows), dpi=dpi)
    if nrows == 1:
        axes = [axes]

    axes[0].hist(ncond, bins=100, edgecolor='black', linewidth=0.5)
    if threshold is not None:
        axes[0].axvline(threshold[0], color='r', linestyle='--', linewidth=1.5,
                        label=f'threshold = {threshold[0]}')
        axes[0].legend()
    axes[0].set_xlabel(r"$N_{\mathrm{cond}}$")
    axes[0].set_ylabel("Count")
    axes[0].set_title(r"Condition number distribution at $t = " + f"{time:.4f}" + r"$")

    if has_eigvec:
        eigvec = data["eigenvecMin"][:]
        axes[1].hist(np.abs(eigvec[:,0]), bins=100, edgecolor='black', linewidth=0.5, alpha=0.7, label=r"$|e_x|$")
        axes[1].hist(np.abs(eigvec[:,1]), bins=100, edgecolor='black', linewidth=0.5, alpha=0.7, label=r"$|e_y|$")
        if threshold is not None:
            axes[1].axvline(threshold[1], color='b', linestyle='--', linewidth=1.0,
                            label=rf'$|e_x| = {threshold[1]}$')
            axes[1].axvline(threshold[2], color='orange', linestyle='--', linewidth=1.0,
                            label=rf'$|e_y| = {threshold[2]}$')
        axes[1].legend()
        axes[1].set_xlabel(r"$|$Component$|$")
        axes[1].set_ylabel("Count")
        axes[1].set_title(r"Eigenvector $|\vec{e}_{\lambda_{\min}}|$ components at $t = " + f"{time:.4f}" + r"$")

    plt.tight_layout()
    histPath = outDir + "/ncond_hist" + pathlib.Path(h5File).stem + ".png"
    print("Saving figure to", histPath)
    plt.savefig(histPath)
    plt.close()

    # Eigenvalue scatter plots (side-by-side)
    if "lambdaMin" in data and "lambdaMax" in data:
        lambdaMin = data["lambdaMin"][()]
        lambdaMax = data["lambdaMax"][()]

        fig, axes = plt.subplots(1, 2, figsize=(14, 5), dpi=dpi)
        for ax, vals, label in zip(axes, [lambdaMin, lambdaMax],
                                   [r"$\lambda_{\min}$", r"$\lambda_{\max}$"]):
            sc = ax.scatter(pos[:,0], pos[:,1], c=vals, s=markerSize)
            setDomainLimits(ax, pos, h5File, openBorders)
            ax.set_aspect('equal')
            ax.set_title(label + r" at $t = " + f"{time:.4f}" + r"$")
            ax.set_xlabel("$x$")
            ax.set_ylabel("$y$")
            divider = make_axes_locatable(ax)
            cax = divider.append_axes("right", size="5%", pad=0.05)
            fig.colorbar(sc, cax=cax)
        plt.tight_layout()
        eigPath = outDir + "/eigenvalues" + pathlib.Path(h5File).stem + ".png"
        print("Saving figure to", eigPath)
        plt.savefig(eigPath)
        plt.close()

    # Eigenvector quiver plot
    if "eigenvecMin" in data:
        eigvec = data["eigenvecMin"][:]
        fig, ax = plt.subplots(figsize=(7, 6), dpi=dpi)
        sc = ax.scatter(pos[:,0], pos[:,1], c=ncond, s=markerSize, vmin=vmin, vmax=vmax)
        step = max(1, len(pos) // 1000)
        ax.quiver(pos[::step,0], pos[::step,1],
                  eigvec[::step,0], eigvec[::step,1],
                  angles='xy', scale_units='xy', scale=30., width=0.002, color='red', alpha=0.6)
        setDomainLimits(ax, pos, h5File, openBorders)
        ax.set_aspect('equal')
        ax.set_title(r"$\vec{e}_{\lambda_{\min}}$ at $t = " + f"{time:.4f}" + r"$")
        ax.set_xlabel("$x$")
        ax.set_ylabel("$y$")
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        fig.colorbar(sc, cax=cax)
        plt.tight_layout()
        evecPath = outDir + "/eigvecMin" + pathlib.Path(h5File).stem + ".png"
        print("Saving figure to", evecPath)
        plt.savefig(evecPath)
        plt.close()

def createCombinedPlot(h5File, outDir, openBorders=False, vminmax=None, markerSize=1.,
                       diff=False, first_frame_data=None, diff_vminmax=None, iHi=-1, dpi=200):
    data = open_snapshot(h5File)
    pos  = data["x"][:]
    time = data["time"][0]

    # Only include quantities the snapshot actually carries; GIZMO snapshots
    # provide rho and u but none of the stress fields, so the figure shrinks
    # gracefully instead of crashing on a missing key.
    candidate_quantities = [
        ("rho", r"Density $\varrho$"),
        ("P",   r"Pressure $P$"),
        ("u",   r"Internal energy $u$"),
        ("Sxx", r"$S_{xx}$"),
        ("Sxy", r"$S_{xy}$"),
        ("Syy", r"$S_{yy}$"),
    ]
    quantities = [(key, data[key][()], label)
                  for key, label in candidate_quantities if key in data]
    quantities_dict = {key: vals for key, vals, _ in quantities}

    plt.rcParams.update({'font.size': 12})
    show_diff = diff and first_frame_data is not None
    fig, axes = plt.subplots(2, 3, figsize=(10, 6), dpi=dpi)
    axes_flat = axes.flatten()

    if show_diff:
        DIFF_FLOOR = 1e-12
        diff_meta = [(dlabel, key) for dlabel, key in [
            (r"$|\Delta\varrho|$", "rho"),
            (r"$|\Delta P|$",      "P"),
            (r"$|\Delta u|$",      "u"),
            (r"$|\Delta S_{xx}|$", "Sxx"),
            (r"$|\Delta S_{xy}|$", "Sxy"),
            (r"$|\Delta S_{yy}|$", "Syy"),
        ] if key in quantities_dict and key in first_frame_data]
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
            sc = ax.scatter(pos[:,0], pos[:,1], c=vals, s=_marker_sizes(data, markerSize), vmin=vlo, vmax=vhi)

            setDomainLimits(ax, pos, h5File, openBorders)
            if iHi > -1:
                plotKernelCircle(data, iHi, pos, ax)
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


def _worker(task):
    (h5File, outDir, grad, vel, stress, iNNL, borders, vmin, vmax,
     pressure, energy, noi, conditionNumber, condThreshold,
     combined, combined_vminmax,
     markerSize, diff, first_frame_data, diff_vminmax, iHi, dpi) = task
    if combined:
        createCombinedPlot(h5File, outDir, borders, combined_vminmax, markerSize,
                           diff, first_frame_data, diff_vminmax, iHi, dpi)
    elif pressure:
        createPressurePlot(h5File, outDir, borders, vmin, vmax, markerSize, iHi, dpi)
    elif energy:
        createEnergyPlot(h5File, outDir, borders, vmin, vmax, markerSize, iHi, dpi)
    elif noi:
        createNoiPlot(h5File, outDir, borders, vmin, vmax, markerSize, iHi, dpi)
    elif conditionNumber:
        createConditionNumberPlot(h5File, outDir, borders, vmin, vmax, markerSize, condThreshold, iHi, dpi)
    else:
        createPlot(h5File, outDir, grad, vel, stress, iNNL, borders, vmin, vmax, markerSize, iHi, dpi)


if __name__=="__main__":
    parser = argparse.ArgumentParser(description="Plot density of results from Kelvin-Helmholtz test case.")
    parser.add_argument("--simOutputDir", "-d", metavar="string", type=str, help="output directory of simulation", required=True)
    parser.add_argument("--outDir", "-o", metavar="string", type=str, help="output directory for generated plots", default="output")
    parser.add_argument("--plotGradient", "-g", action="store_true", help="plot density gradients")
    parser.add_argument("--plotGhosts", "-G", action="store_true", help="also plot ghost cells in an extra file")
    parser.add_argument("--pressure", "-P", action="store_true", help="plot pressure instead of density")
    parser.add_argument("--energy", "-u", action="store_true", help="plot internal energy instead of density")
    parser.add_argument("--stress", "-S", action="store_true", help="Plot absolute value of stress (|S| = sqrt(Sxx^2 + Sxy^2 + Syy^2))")
    parser.add_argument("--combined", "-C", action="store_true",
                        help="Plot rho, P, u, Sxx, Sxy, Syy in a 2x3 combined figure")
    parser.add_argument("--noi", "-n", action="store_true", help="plot number of interactions instead of density")
    parser.add_argument("--conditionNumber", "-N", action="store_true", help="plot condition number of gradient estimation matrix")
    parser.add_argument("--threshold", "-T", nargs=3, metavar=("NCOND", "EX", "EY"), type=float, default=None,
                        help="highlight particles meeting all three: Ncond > NCOND, |e_x| < EX, |e_y| < EY (use with -N)")
    parser.add_argument("--plotVelocity", "-v", action="store_true", help="plot velocity")
    parser.add_argument("--iNNL", "-i", metavar="int", type=int, help="plot NNL for particles i", default=-1)
    parser.add_argument("--hi", metavar="int", type=int, help="plot kernel length h of particle i as a thin circle", default=-1)
    parser.add_argument("--openBorders", "-b", action="store_true", help="Adjust plot domain to show all real particles")
    parser.add_argument("--continue", "-c", dest="continue_", action="store_true",
                        help="skip h5 files whose plots already exist in outDir")
    parser.add_argument("--workers", "-w", metavar="int", type=int, default=os.cpu_count(),
                        help="number of parallel worker processes (default: all CPUs)")
    parser.add_argument("--markerSize", "-m", metavar="float", type=float, default=1.,
                        help="scatter marker size for all plots (default: 1.0)")
    parser.add_argument("--tstart", "-t", metavar="float", type=float, default=None,
                        help="only plot files with simulation time >= tstart")
    parser.add_argument("--diff", "-D", action="store_true",
                        help="(with --combined) add two rows showing log-scale |diff to first frame|")
    parser.add_argument("--dpi", metavar="int", type=int, default=200,
                        help="output resolution in dots per inch (default: 200)")

    args = parser.parse_args()

    if args.diff and not args.combined:
        print("WARNING: --diff has no effect without --combined.")
    if args.threshold is not None and not args.conditionNumber:
        print("WARNING: --threshold has no effect without --conditionNumber.")

    plt.rc('text', usetex=True)
    
    print("Examining files in", args.simOutputDir, "...")

    # Pick up both demonstrator (.h5) and GIZMO (.hdf5) snapshots.
    h5Files = sorted([f for f in list(pathlib.Path(args.simOutputDir).glob('*.h5'))
                                  + list(pathlib.Path(args.simOutputDir).glob('*.hdf5'))
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

    if args.tstart is not None:
        before = len(h5Files)
        h5Files = [f for f in h5Files
                    if open_snapshot(f)["time"][0] >= args.tstart]
        print(f"--tstart={args.tstart}: kept {len(h5Files)}/{before} file(s) with t >= {args.tstart}")

    def prescan_quantity(files, key):
        lo, hi = None, None
        for f in files:
            d = open_snapshot(f)
            if key not in d:
                # GIZMO snapshots lack the demonstrator-only fields (rhoGrad,
                # Sxx, conditionNumber, ...); skip rather than crash.
                continue
            vals = d[key][()] if hasattr(d[key], 'shape') else np.asarray(d[key])
            flo, fhi = float(np.nanmin(vals)), float(np.nanmax(vals))
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
                d = open_snapshot(all_h5Files[0])
                for key in ["rho", "P", "u", "Sxx", "Sxy", "Syy"]:
                    if key in d:
                        first_frame_data[key] = d[key][()]
                print(f"First frame loaded from {all_h5Files[0].name} for diff computation.")

                diff_vminmax = {}
                print(f"Pre-scanning {len(all_h5Files)} file(s) for diff color limits...")
                for key in first_frame_data:
                    d_lo, d_hi = None, None
                    for f in all_h5Files:
                        d = open_snapshot(f)
                        if key not in d:
                            continue
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
            elif args.conditionNumber:
                colorKey = "conditionNumber"
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
              args.pressure, args.energy, args.noi, args.conditionNumber,
              args.threshold, args.combined, combined_vminmax,
              args.markerSize, args.diff, first_frame_data, diff_vminmax, args.hi, args.dpi)
             for f in h5Files]

    nworkers = min(args.workers, len(tasks)) if tasks else 1
    print(f"Rendering {len(tasks)} plot(s) with {nworkers} worker(s)...")
    with Pool(nworkers, initializer=_init_tex_worker) as pool:
        for i, _ in enumerate(pool.imap_unordered(_worker, tasks), 1):
            print(f"  {i}/{len(tasks)} done", end='\r', flush=True)
    print()
            

    print("... done.")
    
