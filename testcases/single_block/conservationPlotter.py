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
    parser.add_argument("--includeElastic", "-e", action="store_true", help="add elastic strain energy to total energy and plot stress sum deviation")
    parser.add_argument("--shearModulus", "-G", metavar="float", type=float, default=1.0, help="shear modulus G for elastic energy computation (default: 1.0)")
    args = parser.parse_args()
    
    plt.rc('text', usetex=True)
    #plt.rc('text', usetex=True)
    plt.rc('text.latex', preamble=r'\usepackage{amsmath}')
    #\renewcommand\vec[1]{\bm{#1}}')
    
    print("Examining files in", args.simOutputDir, "...")

    mass = []
    energy = []
    momX = []
    momY = []
    sumSAbs = []
    time = []

    setReference = True

    def elastic_energy(data, G):
        Sxx = data["Sxx"][:]
        Sxy = data["Sxy"][:]
        Syy = data["Syy"][:]
        V   = data["m"][:] / data["rho"][:]
        return np.sum(V * (Sxx**2 + 2.*Sxy**2 + Syy**2)) / (4. * G)

    def stress_sum(data):
        Sxx = data["Sxx"][:]
        Sxy = data["Sxy"][:]
        Syy = data["Syy"][:]
        return float(np.sum(np.sqrt(Sxx**2 + Sxy**2 + Syy**2)))

    for h5File in sorted(pathlib.Path(args.simOutputDir).glob('*.h5')):
        data = h5.File(h5File, 'r')
        E = data["energy"][0]
        if args.includeElastic:
            E += elastic_energy(data, args.shearModulus)
        if setReference:
            refMass   = data["totalMass"][0]
            refEnergy = E
            refMomX   = data["xMomentum"][0]
            refMomY   = data["yMomentum"][0]
            if args.includeElastic:
                refSumSAbs = stress_sum(data)
            setReference = False
        time.append(data["time"][0])
        mass.append(abs(data["totalMass"][0] - refMass))
        energy.append(abs(E - refEnergy))
        momX.append(abs(data["xMomentum"][0] - refMomX))
        momY.append(abs(data["yMomentum"][0] - refMomY))
        if args.includeElastic:
            sumSAbs.append(abs(stress_sum(data) - refSumSAbs))

    print("... plotting ... ")    

    plt.rcParams.update({'font.size': 18})
    
    fig, ax = plt.subplots(figsize=(14,12), dpi=200)

    print("M_tot =",  mass)
    print("pX_tot =", momX)
    print("pY_tot =", momY)
    print("E_tot =", energy)
    
    order = np.argsort(time)
    time   = np.array(time)[order]
    mass   = np.array(mass)[order]
    energy = np.array(energy)[order]
    momX   = np.array(momX)[order]
    momY   = np.array(momY)[order]
    if args.includeElastic:
        sumSAbs = np.array(sumSAbs)[order]

    if not args.MFM:
        plt.plot(time, mass, 'ro', label=r'$\Delta M_\text{tot}$')
    plt.plot(time, momX, 'bv', label=r'$\Delta p_{x, \text{tot}}$')
    plt.plot(time, momY, 'g^', label=r'$\Delta p_{y, \text{tot}}$')
    plt.plot(time, energy, 'kx', label=r'$\Delta E_\text{tot}$')
    if args.includeElastic:
        plt.plot(time, sumSAbs, 'ms', label=r'$\Delta \sum |S|$')
    plt.yscale('log')
    
    plt.title(r"Absolute numerical error $\Delta_\text{num}$ over time $t$")
    #plt.title(r"Color coded pressure $P$")
    plt.xlabel(r"Time $t$")
    plt.ylabel(r"$\Delta_\text{num}$")

    plt.legend(loc='lower left')
    plt.grid()
    
    plt.tight_layout()
    print("Saving figure to conservation" + args.label + ".png")
    plt.savefig("conservation" + args.label + ".png")
    plt.close()
    
    print("... done.")
    
