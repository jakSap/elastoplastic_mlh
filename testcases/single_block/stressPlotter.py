#!/usr/bin/env python3
"""
stressPlotter.py — plot stress tensor components (Sxx, Sxy, Syy) and the
absolute value |S| = sqrt(Sxx^2 + Sxy^2 + Syy^2) over time, each curve
normalized individually to its own maximum so they share one axis.

For each h5 snapshot the per-particle field is reduced to a scalar via
sum_i |.|, i.e. the L1 sum over all particles.
"""

import argparse
import pathlib
import numpy as np
import h5py as h5
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Plot normalized stress tensor components and |S| over time.")
    parser.add_argument("--simOutputDir", "-d", metavar="string", type=str,
                        required=True, help="output directory of simulation")
    parser.add_argument("--label", "-l", metavar="string", type=str,
                        default="Unlabeled",
                        help="label appended to the output file name")
    parser.add_argument("--tmin", metavar="float", type=float, default=None,
                        help="lower bound of the time axis")
    parser.add_argument("--tmax", metavar="float", type=float, default=None,
                        help="upper bound of the time axis")
    args = parser.parse_args()

    SCALE = 1.3  # uniform 30% size increase

    plt.rc('text', usetex=True)
    plt.rc('text.latex', preamble=r'\usepackage{amsmath}')

    print("Examining files in", args.simOutputDir, "...")

    time   = []
    sxxSum = []
    sxySum = []
    syySum = []
    sAbsSum = []
    eKin   = []
    eInt   = []

    for h5File in sorted(pathlib.Path(args.simOutputDir).glob('*.h5')):
        if "NNL" in str(h5File) or "Ghost" in str(h5File):
            continue
        with h5.File(h5File, 'r') as data:
            Sxx = data["Sxx"][:]
            Sxy = data["Sxy"][:]
            Syy = data["Syy"][:]
            m   = data["m"][:]
            v   = data["v"][:]
            u   = data["u"][:]
            time.append(float(data["time"][0]))
        sxxSum.append(float(np.sum(np.abs(Sxx))))
        sxySum.append(float(np.sum(np.abs(Sxy))))
        syySum.append(float(np.sum(np.abs(Syy))))
        sAbsSum.append(float(np.sum(np.sqrt(Sxx**2 + Sxy**2 + Syy**2))))
        eKin.append(float(0.5 * np.sum(m * (v[:, 0]**2 + v[:, 1]**2))))
        eInt.append(float(np.sum(m * u)))

    time    = np.array(time)
    sxxSum  = np.array(sxxSum)
    sxySum  = np.array(sxySum)
    syySum  = np.array(syySum)
    sAbsSum = np.array(sAbsSum)
    eKin    = np.array(eKin)
    eInt    = np.array(eInt)

    order = np.argsort(time)
    time    = time[order]
    sxxSum  = sxxSum[order]
    sxySum  = sxySum[order]
    syySum  = syySum[order]
    sAbsSum = sAbsSum[order]
    eKin    = eKin[order]
    eInt    = eInt[order]

    def normalize(a):
        m = np.max(np.abs(a))
        return a / m if m > 0 else a, m

    sxxN,  sxxMax  = normalize(sxxSum)
    sxyN,  sxyMax  = normalize(sxySum)
    syyN,  syyMax  = normalize(syySum)
    sAbsN, sAbsMax = normalize(sAbsSum)
    eKinN, eKinMax = normalize(eKin)
    eIntN, eIntMax = normalize(eInt)

    print(f"Sxx max:  {sxxMax:.4g}")
    print(f"Sxy max:  {sxyMax:.4g}")
    print(f"Syy max:  {syyMax:.4g}")
    print(f"|S| max:  {sAbsMax:.4g}")
    print(f"E_kin max: {eKinMax:.4g}")
    print(f"E_int max: {eIntMax:.4g}")

    plt.rcParams.update({
        'font.size':         18  * SCALE,
        'axes.titlesize':    18  * SCALE,
        'axes.labelsize':    18  * SCALE,
        'xtick.labelsize':   18  * SCALE,
        'ytick.labelsize':   18  * SCALE,
        'legend.fontsize':   16  * SCALE,
        'lines.linewidth':   1.5 * SCALE,
        'axes.linewidth':    0.8 * SCALE,
        'xtick.major.width': 0.8 * SCALE,
        'ytick.major.width': 0.8 * SCALE,
        'xtick.major.size':  3.5 * SCALE,
        'ytick.major.size':  3.5 * SCALE,
    })
    fig, ax = plt.subplots(figsize=(12 * SCALE, 8 * SCALE), dpi=200)

    ax.plot(time, sxxN,  'r-',  label=rf"$\sum_i |S_{{xx}}|$ (max $= {sxxMax:.3g}$)")
    ax.plot(time, sxyN,  'g--', label=rf"$\sum_i |S_{{xy}}|$ (max $= {sxyMax:.3g}$)")
    ax.plot(time, syyN,  'b-.', label=rf"$\sum_i |S_{{yy}}|$ (max $= {syyMax:.3g}$)")
    ax.plot(time, sAbsN, 'k:',  label=rf"$\sum_i |S|$ (max $= {sAbsMax:.3g}$)",
            linewidth=2.5 * SCALE)
    ax.plot(time, eKinN, 'm-',  label=rf"$E_\text{{kin}}$ (max $= {eKinMax:.3g}$)")
    ax.plot(time, eIntN, 'c--', label=rf"$E_\text{{int}}$ (max $= {eIntMax:.3g}$)")

    ax.set_title(r"Normalized stress tensor components over time")
    ax.set_xlabel(r"Time $t$")
    ax.set_ylabel(r"Normalized $\sum_i |\,\cdot\,| / \max$")
    ax.legend(loc='best')
    ax.grid(True)

    if args.tmin is not None or args.tmax is not None:
        cur_lo, cur_hi = ax.get_xlim()
        ax.set_xlim(args.tmin if args.tmin is not None else cur_lo,
                    args.tmax if args.tmax is not None else cur_hi)

    plt.tight_layout()
    out = "stress" + args.label + ".png"
    print("Saving figure to", out)
    plt.savefig(out)
    plt.close()

    print("... done.")
