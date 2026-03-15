#!/usr/bin/env python3

import argparse
import pathlib
import numpy as np
import h5py as h5
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Plot conserved qunatities over time for the Kelvin-Helmholtz test case.")
    parser.add_argument("--simOutputDir", "-d", metavar="string", type=str, help="output directory of simulation", required=True)
    parser.add_argument("--label", "-l", metavar="string", type=str, help="label appended to the output file name", default="Unlabeled")
    parser.add_argument("--MFM", "-M", action="store_true", help="do not plot mass as the error is zero")
    parser.add_argument("--simOutputDir2", "-d2", metavar="string", type=str, help="output directory of second simulation (optional)", default=None)
    parser.add_argument("--legend1", "-L1", metavar="string", type=str, help="legend label for sim 1 (used only with -d2)", default="Sim 1")
    parser.add_argument("--legend2", "-L2", metavar="string", type=str, help="legend label for sim 2 (used only with -d2)", default="Sim 2")
    parser.add_argument("--markerSize", "-ms", metavar="float", type=float, help="marker size", default=4)
    args = parser.parse_args()

    plt.rc('text', usetex=True)
    #plt.rc('text', usetex=True)
    plt.rc('text.latex', preamble=r'\usepackage{amsmath}')
    #\renewcommand\vec[1]{\bm{#1}}')

    def load_sim(directory):
        mass = []
        energy = []
        momX = []
        momY = []
        time = []
        setReference = True
        for h5File in pathlib.Path(directory).glob('*.h5'):
            data = h5.File(h5File, 'r')
            if setReference:
                refMass   = data["totalMass"][0]
                refEnergy = data["energy"][0]
                refMomX   = data["xMomentum"][0]
                refMomY   = data["yMomentum"][0]
                setReference = False
            time.append(data["time"][0])
            mass.append(abs(data["totalMass"][0]  - refMass))
            energy.append(abs(data["energy"][0]   - refEnergy))
            momX.append(abs(data["xMomentum"][0]  - refMomX))
            momY.append(abs(data["yMomentum"][0]  - refMomY))
        return mass, energy, momX, momY, time

    print("Examining files in", args.simOutputDir, "...")
    mass, energy, momX, momY, time = load_sim(args.simOutputDir)

    print("... plotting ... ")

    plt.rcParams.update({'font.size': 18})

    fig, ax = plt.subplots(figsize=(14,12), dpi=200)

    print("M_tot =",  mass)
    print("pX_tot =", momX)
    print("pY_tot =", momY)
    print("E_tot =", energy)

    kw = dict(linestyle='none', markersize=args.markerSize)

    if args.simOutputDir2 is None:
        # Single-sim mode
        if not args.MFM:
            plt.plot(time, mass,   color='tab:red',    marker='o', label=r'$\Delta M_\text{tot}$',      **kw)
        plt.plot(time, momX,   color='tab:blue',   marker='v', label=r'$\Delta p_{x, \text{tot}}$', **kw)
        plt.plot(time, momY,   color='tab:orange', marker='^', label=r'$\Delta p_{y, \text{tot}}$', **kw)
        plt.plot(time, energy, color='tab:purple', marker='x', label=r'$\Delta E_\text{tot}$',      **kw)
    else:
        # Comparison mode
        print("Examining files in", args.simOutputDir2, "...")
        mass2, energy2, momX2, momY2, time2 = load_sim(args.simOutputDir2)
        print("M_tot2 =",  mass2)
        print("pX_tot2 =", momX2)
        print("pY_tot2 =", momY2)
        print("E_tot2 =", energy2)

        L1 = args.legend1
        L2 = args.legend2

        if not args.MFM:
            plt.plot(time,  mass,   color='tab:red',      marker='o', label=L1 + r': $\Delta M_\text{tot}$',      **kw)
            plt.plot(time2, mass2,  color='lightcoral',   marker='s', label=L2 + r': $\Delta M_\text{tot}$',      **kw)
        plt.plot(time,  momX,   color='tab:blue',     marker='v', label=L1 + r': $\Delta p_{x, \text{tot}}$', **kw)
        plt.plot(time2, momX2,  color='cornflowerblue',marker='^', label=L2 + r': $\Delta p_{x, \text{tot}}$', **kw)
        plt.plot(time,  momY,   color='tab:orange',   marker='^', label=L1 + r': $\Delta p_{y, \text{tot}}$', **kw)
        plt.plot(time2, momY2,  color='gold',         marker='v', label=L2 + r': $\Delta p_{y, \text{tot}}$', **kw)
        plt.plot(time,  energy, color='tab:purple',   marker='x', label=L1 + r': $\Delta E_\text{tot}$',      **kw)
        plt.plot(time2, energy2,color='mediumpurple', marker='+', label=L2 + r': $\Delta E_\text{tot}$',      **kw)

    plt.yscale('log')

    plt.title(r"Absolute numerical error $\Delta_\text{num}$ over time $t$")
    plt.xlabel(r"Time $t$")
    plt.ylabel(r"$\Delta_\text{num}$")

    plt.legend(loc='lower left')
    plt.grid()

    plt.tight_layout()
    print("Saving figure to conservation" + args.label + ".png")
    plt.savefig("conservation" + args.label + ".png")
    plt.close()

    print("... done.")

