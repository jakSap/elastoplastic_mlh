#!/usr/bin/env python3

import argparse
import pathlib
import numpy as np
import h5py as h5
import matplotlib.pyplot as plt

PARTICLE_KEYS = ['P', 'Sxx', 'Sxy', 'Syy', 'm', 'noi', 'rho', 'u']
VECTOR_KEYS = ['v', 'x', 'rhoGrad']

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Plot all attributes of a single particle over time.")
    parser.add_argument("--simOutputDir", "-d", metavar="string", type=str, help="output directory of simulation", required=True)
    parser.add_argument("--label", "-l", metavar="string", type=str, help="label appended to the output file name", default="Unlabeled")
    parser.add_argument("--index", "-i", metavar="int", type=int, help="particle index", required=True)
    parser.add_argument("--log", action="store_true", help="use logarithmic y-axis")
    parser.add_argument("--tstart", "-t", metavar="float", type=float, help="start time for plotting (earlier data is discarded)", default=None)
    parser.add_argument("--clamp", action="store_true", help="limit y-axis to [0, max(values)]")
    parser.add_argument("--vline", metavar="float", type=float, help="draw a vertical line at this time value", default=None)
    args = parser.parse_args()

    plt.rc('text', usetex=True)
    plt.rc('text.latex', preamble=r'\usepackage{amsmath}')

    print("Examining files in", args.simOutputDir, "...")

    h5Files = sorted(pathlib.Path(args.simOutputDir).glob('*.h5'))

    time = []
    scalars = {k: [] for k in PARTICLE_KEYS}
    vectors = {}
    for k in VECTOR_KEYS:
        vectors[k + '_0'] = []
        vectors[k + '_1'] = []

    for h5File in h5Files:
        data = h5.File(h5File, 'r')
        i = args.index
        if i >= data[PARTICLE_KEYS[0]].shape[0]:
            print(f"Warning: index {i} out of range in {h5File}, skipping.")
            data.close()
            continue
        t = data["time"][0]
        if args.tstart is not None and t < args.tstart:
            data.close()
            continue
        time.append(t)
        for k in PARTICLE_KEYS:
            scalars[k].append(data[k][i])
        for k in VECTOR_KEYS:
            vectors[k + '_0'].append(data[k][i, 0])
            vectors[k + '_1'].append(data[k][i, 1])
        data.close()

    all_series = {}
    all_series.update(scalars)
    all_series.update(vectors)

    labels = {
        'P': r'$P$', 'Sxx': r'$S_{xx}$', 'Sxy': r'$S_{xy}$', 'Syy': r'$S_{yy}$',
        'm': r'$m$', 'noi': r'noi', 'rho': r'$\rho$', 'u': r'$u$',
        'v_0': r'$v_x$', 'v_1': r'$v_y$',
        'x_0': r'$x$', 'x_1': r'$y$',
        'rhoGrad_0': r'$\nabla\rho_x$', 'rhoGrad_1': r'$\nabla\rho_y$',
    }

    print("... plotting ...")

    plt.rcParams.update({'font.size': 25})
    fig, ax = plt.subplots(figsize=(12,18), dpi=200)

    for key, values in all_series.items():
        arr = np.array(values)
        initial = arr[0] if arr[0] != 0.0 else 1.0
        ax.plot(time, arr / initial, label=labels.get(key, key))

    if args.log:
        ax.set_yscale('log')

    if args.clamp:
        all_vals = np.concatenate([np.array(v) / (np.array(v)[0] if np.array(v)[0] != 0.0 else 1.0) for v in all_series.values()])
        # ax.set_ylim(0, np.max(all_vals))
        ax.set_ylim(-.1, 10)

    if args.vline is not None:
        ax.axvline(x=args.vline, color='black', linestyle='--', linewidth=1)

    ax.set_title(f"Particle {args.index} attributes over time")
    ax.set_xlabel(r"Time $t$")
    ax.set_ylabel(r"Value / initial value")
    ax.legend(loc='best', ncol=3, fontsize=10)
    ax.grid()

    plt.tight_layout()
    outname = f"particle_{args.index}_{args.label}.png"
    print(f"Saving figure to {outname}")
    plt.savefig(outname)
    plt.show()
    plt.close()

    print("... done.")
