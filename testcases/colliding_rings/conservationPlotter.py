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
    parser.add_argument("--elasticEnergy", "-E", action="store_true",
                        help="plot the elastic strain energy instead of the kinetic+internal energy")
    parser.add_argument("--shearModulus", "-mu", metavar="float", type=float, default=0.22,
                        help="shear modulus used for the elastic strain energy (default 0.22, the rings test case)")
    args = parser.parse_args()

    plt.rc('text', usetex=True)
    #plt.rc('text', usetex=True)
    plt.rc('text.latex', preamble=r'\usepackage{amsmath}')
    #\renewcommand\vec[1]{\bm{#1}}')

    def detect_format(data):
        """Identify the output format from the datasets present in an open HDF5 file."""
        if "PartType0" in data:
            return "gizmo"
        if "totalMass" in data and "xMomentum" in data:
            return "demonstrator"
        if "e" in data and "deviatoric_stress" in data:
            return "miluphcuda"
        # Fall back on the demonstrator layout if particle fields look familiar.
        if {"m", "x", "v"}.issubset(set(data.keys())):
            return "miluphcuda"
        raise ValueError("Unrecognised output format; datasets: " + ", ".join(sorted(data.keys())))

    def read_particles(data, fmt):
        """Return (time, m, x, v, u, ss, rho) with x, v as (N, 2), u the specific
        internal energy, ss the deviatoric-stress double contraction S_ij S_ij and
        rho the density.

        All conserved quantities are derived uniformly from these per-particle fields,
        so the three formats are compared on an identical definition. For the
        demonstrator this reproduces its own stored totals exactly (its `energy`
        equals kinetic + internal energy)."""
        if fmt == "gizmo":
            p = data["PartType0"]
            m = p["Masses"][:].astype(np.float64)
            x = p["Coordinates"][:, :2].astype(np.float64)
            v = p["Velocities"][:, :2].astype(np.float64)
            u = p["InternalEnergy"][:].astype(np.float64)
            rho = p["Density"][:].astype(np.float64)
            # Full 3x3 traceless deviator (pressure stored separately): S:S = sum of squares.
            S = p["StressTensor"][:].astype(np.float64)
            ss = np.sum(S * S, axis=1)
            t = float(data["Header"].attrs["Time"])
        elif fmt == "demonstrator":
            m = data["m"][:].astype(np.float64)
            x = data["x"][:].astype(np.float64)
            v = data["v"][:].astype(np.float64)
            u = data["u"][:].astype(np.float64)
            rho = data["rho"][:].astype(np.float64)
            Sxx = data["Sxx"][:].astype(np.float64)
            Sxy = data["Sxy"][:].astype(np.float64)
            Syy = data["Syy"][:].astype(np.float64)
            # Symmetric 2D deviator: S:S = Sxx^2 + 2 Sxy^2 + Syy^2.
            ss = Sxx**2 + 2.0 * Sxy**2 + Syy**2
            t = float(data["time"][0])
        elif fmt == "miluphcuda":
            m = data["m"][:].astype(np.float64)
            x = data["x"][:].astype(np.float64)
            v = data["v"][:].astype(np.float64)
            u = data["e"][:].astype(np.float64)
            rho = data["rho"][:].astype(np.float64)
            # deviatoric_stress stores all four components [Sxx, Sxy, Syx, Syy].
            S = data["deviatoric_stress"][:].astype(np.float64)
            ss = np.sum(S * S, axis=1)
            t = float(data["time"][0])
        else:
            raise ValueError("Unknown format: " + fmt)
        return t, m, x[:, :2], v[:, :2], u, ss, rho

    def conserved(m, x, v, u, ss, rho, mu, elastic):
        """Total mass, x/y momentum, z angular momentum and energy.

        The energy slot holds the elastic strain energy E = sum_i V_i S:S_i / (4 mu)
        when `elastic` is set, otherwise the kinetic + internal energy."""
        mass   = np.sum(m)
        momX   = np.sum(m * v[:, 0])
        momY   = np.sum(m * v[:, 1])
        angMom = np.sum(m * (x[:, 0] * v[:, 1] - x[:, 1] * v[:, 0]))
        if elastic:
            energy = np.sum((m / rho) * ss / (4.0 * mu))
        else:
            energy = np.sum(m * (0.5 * (v[:, 0]**2 + v[:, 1]**2) + u))
        return mass, energy, momX, momY, angMom

    def load_sim(directory):
        # Gather every snapshot (.h5 demonstrator/miluphcuda, .hdf5 GIZMO).
        files = sorted(pathlib.Path(directory).glob('*.h5')) \
              + sorted(pathlib.Path(directory).glob('*.hdf5'))
        if not files:
            raise FileNotFoundError("No .h5 or .hdf5 snapshots found in " + str(directory))

        records = []
        fmt = None
        for h5File in files:
            with h5.File(h5File, 'r') as data:
                f = detect_format(data)
                if fmt is None:
                    fmt = f
                    print("  detected format:", fmt)
                t, m, x, v, u, ss, rho = read_particles(data, f)
            records.append((t,) + conserved(m, x, v, u, ss, rho,
                                            args.shearModulus, args.elasticEnergy))

        # Sort by time so the reference is the earliest snapshot, not glob order.
        records.sort(key=lambda r: r[0])
        time   = [r[0] for r in records]
        refMass, refEnergy, refMomX, refMomY, refAngMom = records[0][1:]

        mass   = [abs(r[1] - refMass)   for r in records]
        energy = [abs(r[2] - refEnergy) for r in records]
        momX   = [abs(r[3] - refMomX)   for r in records]
        momY   = [abs(r[4] - refMomY)   for r in records]
        angMom = [abs(r[5] - refAngMom) for r in records]
        return mass, energy, momX, momY, angMom, time

    print("Examining files in", args.simOutputDir, "...")
    mass, energy, momX, momY, angMom, time = load_sim(args.simOutputDir)

    print("... plotting ... ")

    plt.rcParams.update({'font.size': 18})

    fig, ax = plt.subplots(figsize=(14,12), dpi=200)

    print("M_tot =",  mass)
    print("pX_tot =", momX)
    print("pY_tot =", momY)
    print("L_tot =", angMom)
    print("E_tot =", energy)

    kw = dict(linestyle='none', markersize=args.markerSize)

    energyLabel = r'\Delta E_\text{elastic}' if args.elasticEnergy else r'\Delta E_\text{tot}'

    # One row per conserved quantity. In comparison mode the two simulations share a
    # quantity but are drawn in complementary (opposite-hue) colours so that the two
    # nearly-coincident curves stay cleanly distinguishable; markers differ as well.
    #              quantity label,                 sim1 colour,   sim2 colour,    m1,  m2
    quantities = [
        ('mass',   r'\Delta M_\text{tot}',         'crimson',     'teal',        'o', 's'),
        ('momX',   r'\Delta p_{x, \text{tot}}',    'darkorange',  'royalblue',   'v', '^'),
        ('momY',   r'\Delta p_{y, \text{tot}}',    'goldenrod',   'darkviolet',  '^', 'v'),
        ('angMom', r'\Delta L_{z, \text{tot}}',    'forestgreen', 'magenta',     'd', '*'),
        ('energy', energyLabel,                     'purple',      'yellowgreen', 'x', '+'),
    ]

    series1 = dict(mass=mass, momX=momX, momY=momY, angMom=angMom, energy=energy)

    if args.simOutputDir2 is None:
        # Single-sim mode
        for key, lbl, c1, _c2, m1, _m2 in quantities:
            if key == 'mass' and args.MFM:
                continue
            plt.plot(time, series1[key], color=c1, marker=m1, label=r'$' + lbl + r'$', **kw)
    else:
        # Comparison mode
        print("Examining files in", args.simOutputDir2, "...")
        mass2, energy2, momX2, momY2, angMom2, time2 = load_sim(args.simOutputDir2)
        print("M_tot2 =",  mass2)
        print("pX_tot2 =", momX2)
        print("pY_tot2 =", momY2)
        print("L_tot2 =", angMom2)
        print("E_tot2 =", energy2)

        L1 = args.legend1
        L2 = args.legend2
        series2 = dict(mass=mass2, momX=momX2, momY=momY2, angMom=angMom2, energy=energy2)

        for key, lbl, c1, c2, m1, m2 in quantities:
            if key == 'mass' and args.MFM:
                continue
            plt.plot(time,  series1[key], color=c1, marker=m1, label=L1 + r': $' + lbl + r'$', **kw)
            plt.plot(time2, series2[key], color=c2, marker=m2, label=L2 + r': $' + lbl + r'$', **kw)

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

