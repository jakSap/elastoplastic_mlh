#!/usr/bin/env python3
"""Compare the elastic-wave runs across CFL numbers at 3 wave passings.

Splits the stress-slaving defect into its time and space parts: with
defect(dt) = a dt + b and dt scaled by 1/4,

    a dt_full = (defect_0.2 - defect_0.05) / (1 - 1/4),   b = defect - a dt.

See cflTest.sh for what this is for.
"""

import glob
import sys
import numpy as np
import h5py

from makeMovie import ic_name


def at_3T(outdir, T):
    best = None
    for f in glob.glob(outdir + "/*.h5"):
        with h5py.File(f, "r") as d:
            t = float(np.ravel(d["time"])[0])
        if best is None or abs(t - 3 * T) < abs(best[1] - 3 * T):
            best = (f, t)
    return best


def diagnostics(path, t, md, n):
    A, mu, c_l = md["amplitude"], md["mu"], md["c_l"]
    S, k = (4.0 / 3.0) * mu * A, 2 * np.pi / md["wavelength"]
    with h5py.File(path, "r") as d:
        x, rho, sxx = d["x"][:, 0], d["rho"][:], d["Sxx"][:]
    phase = k * c_l * t
    col = np.floor(x * n).astype(int) % n
    cnt = np.maximum(np.bincount(col, minlength=n), 1)
    xc = (np.arange(n) + 0.5) / n
    defect = (sxx + (4.0 / 3.0) * mu * np.log(rho)) / S
    prof = np.bincount(col, weights=defect, minlength=n) / cnt
    pert = (np.bincount(col, weights=rho, minlength=n) / cnt - rho.mean()) / (rho.mean() * A)
    M = np.column_stack([np.sin(k * xc - phase), np.cos(k * xc - phase)])
    (p, q), *_ = np.linalg.lstsq(M, pert, rcond=None)
    return dict(defect=np.abs(prof).max(),
                L1_rho=np.mean(np.abs(rho - rho.mean() * (1 + A * np.sin(k * x - phase)))),
                a_rho=float(np.hypot(p, q)), phi=float(np.arctan2(-q, p)))


CFLS = [(0.20, ""), (0.05, "_cfl05"), (0.0125, "_cfl0125")]


def fit_power(E, dts):
    """err = a dt^p + b from >= 3 points on a constant-ratio dt ladder."""
    r = dts[0] / dts[1]
    num, den = E[0] - E[1], E[1] - E[2]
    if num <= 0 or den <= 0:
        return float("nan"), float("nan"), float("nan")
    p = np.log(num / den) / np.log(r)
    a = num / (1.0 - r ** -p)
    return p, a, E[0] - a


def modal_dt(outdir):
    import re
    dts = [float(m) for m in re.findall(r"dt = ([0-9.eE+-]+) selected", open(outdir + "/run.log").read())]
    v, c = np.unique(np.round(dts, 12), return_counts=True)
    return float(v[np.argmax(c)])


if __name__ == "__main__":
    ns = [int(a) for a in sys.argv[1:]] or [64]
    for n in ns:
        with h5py.File(ic_name("elastic_N{}_corr".format(n)), "r") as f:
            md = {k: float(f.attrs[k]) for k in ("amplitude", "wavelength", "mu", "c_l", "period")}
        vals, cfls, dts = [], [], []
        for cfl, sfx in CFLS:
            outdir = "output_elastic_N{}_corr{}".format(n, sfx)
            got = at_3T(outdir, md["period"])
            if not got:
                continue
            vals.append(diagnostics(got[0], got[1], md, n))
            cfls.append(cfl)
            dts.append(modal_dt(outdir))
        if len(vals) < 2:
            print("n = {}: need at least two CFL runs".format(n))
            continue
        print("\n=== n = {}   dt " .format(n)
              + " ".join("{:.4e}".format(d) for d in dts)
              + "   ratios " + " ".join("{:.3f}".format(a / b) for a, b in zip(dts, dts[1:])))
        print("  {:14s}".format("metric")
              + "".join("{:>13s}".format("CFL {:g}".format(c)) for c in cfls)
              + "{:>7s}{:>13s}{:>7s}{:>13s}".format("p", "time part", "%", "floor b"))
        for key, label in (("defect", "S-defect"), ("L1_rho", "L1(rho)")):
            E = [v[key] for v in vals]
            row = "  {:14s}".format(label) + "".join("{:13.4e}".format(e) for e in E)
            if len(E) >= 3:
                pw, a, b = fit_power(E, dts)
                row += "{:7.2f}{:13.4e}{:6.0f}%{:13.4e}".format(pw, a, 100 * a / E[0], b)
            else:
                r = dts[0] / dts[1]
                a = (E[0] - E[1]) / (1.0 - 1.0 / r)
                row += "{:>7s}{:13.4e}{:6.0f}%{:13.4e}".format("(1)", a, 100 * a / E[0], E[0] - a)
            print(row)
        print("  {:14s}".format("a_rho") + "".join("{:13.6f}".format(v["a_rho"]) for v in vals))
        if len(vals) < 3:
            print("  p assumed = 1; run ./cflTest.sh -c .0125 to measure it")
