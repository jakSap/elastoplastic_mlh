#!/usr/bin/env python3

import argparse
import pathlib
import numpy as np
import h5py as h5
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.lines import Line2D
from mpl_toolkits.axes_grid1 import make_axes_locatable


def plot_interleaved(series, markerSize=4, nbands=24, seed=None):
    """Scatter every series with a randomized per-point draw order.

    `series` is a list of (x, y, color, marker, label) tuples. All points from
    all series are pooled and each is assigned a random integer "depth band";
    bands are drawn back-to-front. Because the band is chosen per point (not per
    series), no single quantity ends up uniformly on top of the others where
    they overlap. Returns proxy legend handles (one per series)."""
    rng = np.random.default_rng(seed)

    X, Y, C, M = [], [], [], []
    handles = []
    for x, y, color, marker, label in series:
        x = np.asarray(x, dtype=float)
        y = np.asarray(y, dtype=float)
        n = x.size
        X.append(x)
        Y.append(y)
        C.append(np.tile(mcolors.to_rgba(color), (n, 1)))
        M.append(np.full(n, marker))
        handles.append(Line2D([0], [0], linestyle='none', marker=marker,
                              color=color, markersize=markerSize, label=label))

    X = np.concatenate(X)
    Y = np.concatenate(Y)
    C = np.concatenate(C)
    M = np.concatenate(M)
    bands = rng.integers(0, nbands, size=X.size)

    # scatter needs a uniform marker per call, so group by (marker, band)
    for marker in np.unique(M):
        mMask = M == marker
        for b in range(nbands):
            sel = mMask & (bands == b)
            if np.any(sel):
                plt.scatter(X[sel], Y[sel], c=C[sel], marker=str(marker),
                            s=markerSize**2, zorder=b)
    return handles

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Plot conserved qunatities over time for the Kelvin-Helmholtz test case.")
    parser.add_argument("--simOutputDir", "-d", metavar="string", type=str, help="output directory of simulation", required=True)
    parser.add_argument("--label", "-l", metavar="string", type=str, help="label appended to the output file name", default="Unlabeled")
    parser.add_argument("--MFM", "-M", action="store_true", help="do not plot mass as the error is zero")
    parser.add_argument("--simOutputDir2", "-d2", metavar="string", type=str, help="output directory of second simulation (optional)", default=None)
    parser.add_argument("--legend1", "-L1", metavar="string", type=str, help="legend label for sim 1 (used only with -d2)", default="Sim 1")
    parser.add_argument("--legend2", "-L2", metavar="string", type=str, help="legend label for sim 2 (used only with -d2)", default="Sim 2")
    parser.add_argument("--markerSize", "-ms", metavar="float", type=float, help="marker size", default=4)
    parser.add_argument("--tstart", "-t", metavar="float", type=float, help="exclude data with time < tstart; reference is taken at the first frame >= tstart", default=0.0)
    parser.add_argument("--angularMomentum", "-A", action="store_true", help="also compute and plot angular momentum L_z = sum_i m_i (x_i v_{y,i} - y_i v_{x,i}) from per-particle data")
    parser.add_argument("--nbands", metavar="int", type=int, help="number of random depth bands used to interleave overlapping points", default=24)
    parser.add_argument("--seed", metavar="int", type=int, help="random seed for the point interleaving (omit for non-reproducible shuffling)", default=None)
    args = parser.parse_args()

    plt.rc('text', usetex=True)
    #plt.rc('text', usetex=True)
    plt.rc('text.latex', preamble=r'\usepackage{amsmath}')
    #\renewcommand\vec[1]{\bm{#1}}')

    def load_sim(directory, tstart=0.0, withAngMom=False):
        records = []
        for h5File in pathlib.Path(directory).glob('*.h5'):
            data = h5.File(h5File, 'r')
            rec = {
                "time":   data["time"][0],
                "mass":   data["totalMass"][0],
                "energy": data["energy"][0],
                "momX":   data["xMomentum"][0],
                "momY":   data["yMomentum"][0],
            }
            if withAngMom:
                m = data["m"][:]
                x = data["x"][:]
                v = data["v"][:]
                # z-component of angular momentum about the origin (2D)
                rec["angL"] = np.sum(m * (x[:, 0]*v[:, 1] - x[:, 1]*v[:, 0]))
            records.append(rec)

        # sort by time so the reference is genuinely the earliest kept frame
        # (Path.glob order is arbitrary), then drop frames before tstart
        records.sort(key=lambda r: r["time"])
        records = [r for r in records if r["time"] >= tstart]
        if not records:
            raise SystemExit(f"No frames with time >= {tstart} found in {directory}")

        ref = records[0]
        time   = [r["time"] for r in records]
        mass   = [abs(r["mass"]   - ref["mass"])   for r in records]
        energy = [abs(r["energy"] - ref["energy"]) for r in records]
        momX   = [abs(r["momX"]   - ref["momX"])   for r in records]
        momY   = [abs(r["momY"]   - ref["momY"])   for r in records]
        angL   = [abs(r["angL"]   - ref["angL"])   for r in records] if withAngMom else None
        return mass, energy, momX, momY, angL, time

    print("Examining files in", args.simOutputDir, "...")
    mass, energy, momX, momY, angL, time = load_sim(args.simOutputDir, args.tstart, args.angularMomentum)

    print("... plotting ... ")

    plt.rcParams.update({'font.size': 18})

    fig, ax = plt.subplots(figsize=(14,12), dpi=200)

    print("M_tot =",  mass)
    print("pX_tot =", momX)
    print("pY_tot =", momY)
    print("E_tot =", energy)
    if args.angularMomentum:
        print("L_tot =", angL)

    series = []  # (x, y, color, marker, label)

    if args.simOutputDir2 is None:
        # Single-sim mode
        if not args.MFM:
            series.append((time, mass,   'tab:red',    'o', r'$\Delta M_\text{tot}$'))
        series.append((time, momX,   'tab:blue',   'v', r'$\Delta p_{x, \text{tot}}$'))
        series.append((time, momY,   'tab:orange', '^', r'$\Delta p_{y, \text{tot}}$'))
        series.append((time, energy, 'tab:purple', 'D', r'$\Delta E_\text{tot}$'))
        if args.angularMomentum:
            series.append((time, angL, 'tab:green', 's', r'$\Delta L_{z, \text{tot}}$'))
    else:
        # Comparison mode
        print("Examining files in", args.simOutputDir2, "...")
        mass2, energy2, momX2, momY2, angL2, time2 = load_sim(args.simOutputDir2, args.tstart, args.angularMomentum)
        print("M_tot2 =",  mass2)
        print("pX_tot2 =", momX2)
        print("pY_tot2 =", momY2)
        print("E_tot2 =", energy2)
        if args.angularMomentum:
            print("L_tot2 =", angL2)

        L1 = args.legend1
        L2 = args.legend2

        # Marker shape identifies the quantity (same shape for both runs);
        # the two colors per quantity are high-contrast so the runs never
        # look similarly coloured.
        if not args.MFM:
            series.append((time,  mass,   '#d62728', 'o', L1 + r': $\Delta M_\text{tot}$'))      # red
            series.append((time2, mass2,  '#1f77b4', 'o', L2 + r': $\Delta M_\text{tot}$'))      # blue
        series.append((time,  momX,   '#ff7f0e', 'v', L1 + r': $\Delta p_{x, \text{tot}}$'))     # orange
        series.append((time2, momX2,  '#9467bd', 'v', L2 + r': $\Delta p_{x, \text{tot}}$'))     # purple
        series.append((time,  momY,   '#2ca02c', '^', L1 + r': $\Delta p_{y, \text{tot}}$'))     # green
        series.append((time2, momY2,  '#e377c2', '^', L2 + r': $\Delta p_{y, \text{tot}}$'))     # pink
        series.append((time,  energy, '#8c564b', 's', L1 + r': $\Delta E_\text{tot}$'))          # brown
        series.append((time2, energy2,'#17becf', 's', L2 + r': $\Delta E_\text{tot}$'))          # cyan
        if args.angularMomentum:
            series.append((time,  angL,  '#000000', 'D', L1 + r': $\Delta L_{z, \text{tot}}$'))  # black
            series.append((time2, angL2, '#bcbd22', 'D', L2 + r': $\Delta L_{z, \text{tot}}$'))  # olive

    handles = plot_interleaved(series, markerSize=args.markerSize, nbands=args.nbands, seed=args.seed)

    plt.yscale('log')

    plt.title(r"Absolute numerical error $\Delta_\text{num}$ over time $t$")
    plt.xlabel(r"Time $t$")
    plt.ylabel(r"$\Delta_\text{num}$")

    plt.legend(handles=handles, loc='lower left')
    plt.grid()

    plt.tight_layout()
    print("Saving figure to conservation" + args.label + ".png")
    plt.savefig("conservation" + args.label + ".png")
    plt.close()

    print("... done.")

