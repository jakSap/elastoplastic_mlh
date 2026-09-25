#!/usr/bin/env python3
"""Compare the soundwave runs at t = 1 (one period) across CFL numbers.

A pure n-scan at fixed CFL cannot separate time from space: dt = CFL h / v_sig,
so dt ~ h and an O(dt) term looks exactly like an O(h) term. Holding n fixed and
cutting CFL breaks that link.

Two CFL points are still not enough to say WHICH power of dt: a dt and a dt^2
term both shrink when dt does. With three points spaced by a constant factor r
(here 0.2 / 0.05 / 0.0125, r = 4) and the model err(dt) = a dt^p + b,

    (E1 - E2) / (E2 - E3) = r^p        ->   p = log_r[(E1-E2)/(E2-E3)]
    a dt1^p = (E1 - E2) / (1 - r^-p),       b = E1 - a dt1^p

p ~ 1 is the dangerous case: in an n-scan it masquerades as first-order spatial
error. p ~ 2 is harmless -- dt^2 ~ h^2 is indistinguishable from, and no worse
than, the spatial second order. b is the floor that refining dt alone cannot pass.

Metrics: L1(rho), L1(vx), the amplitude deficit 1 - a_rho(t=1) (the damping that
converged at order 1 for the ideal gas on 17.09.) and |c_num - c_exp| (the
dispersion). Usage: ./cflAnalysis.py [STEM ...]
"""

import os
import re
import sys
import glob
import numpy as np
import h5py

from makeMovie import snapshots, column_profile, sine_fit

GAMMA = 5.0 / 3.0
DT_RE = re.compile(r"dt = ([0-9.eE+-]+) selected")


def parse_stem(stem):
    """<eos>_N<n>[_nnn<N_NN>][_corr][_limNone] -> (eos, n, corr)."""
    parts = stem.split("_")
    return parts[0], int(parts[1][1:]), "corr" in parts[2:]


def ic_name(eos, n, corr):
    cand = "soundwave_{}_N{}{}.h5".format(eos, n, "_corr" if corr else "")
    return cand if os.path.exists(cand) else "soundwave_{}_N{}.h5".format(eos, n)


def typical_dt(outdir):
    """Modal timestep of the run: the unshortened CFL step (dump steps are cut)."""
    dts = [float(m) for m in DT_RE.findall(open(outdir + "/run.log").read())]
    if not dts:
        return float("nan")
    vals, counts = np.unique(np.round(dts, 12), return_counts=True)
    return float(vals[np.argmax(counts)])


def diagnostics(outdir, eos, n, corr):
    with h5py.File(ic_name(eos, n, corr), "r") as f:
        A, lam = float(f.attrs["amplitude"]), float(f.attrs["wavelength"])
    k = 2.0 * np.pi / lam
    files, times = snapshots(outdir)
    i = int(np.argmin([abs(t - 1.0) for t in times]))
    path, t = files[i], times[i]
    with h5py.File(path, "r") as d:
        x, rho, vx, u = d["x"][:, 0], d["rho"][:], d["v"][:, 0], d["u"][:]
    mean = rho.mean()
    wave = np.sin(k * (x - t))
    xc, prof, m = column_profile(path, n)
    a_rho, phi = sine_fit(xc, (prof - m) / (m * A), k, t)
    c_num = 1.0 + phi / (k * t)
    if eos == "ideal":
        c_exp = np.sqrt(GAMMA * (GAMMA - 1.0) * u.mean())
    elif eos == "murnaghan":
        c_exp = np.sqrt(mean)
    else:
        c_exp = 1.0
    return dict(t=t, dt=typical_dt(outdir),
                L1_rho=float(np.mean(np.abs(rho - mean * (1.0 + A * wave)))),
                L1_vx=float(np.mean(np.abs(vx - A * wave))),
                amp_def=float(1.0 - a_rho),
                dc_disp=float(abs(c_num - c_exp)),
                a_rho=float(a_rho), phi=float(phi))


METRICS = [("L1_rho", "L1(rho)"), ("L1_vx", "L1(vx)"),
           ("amp_def", "1 - a_rho"), ("dc_disp", "|c_num-c_exp|")]

# (label, output suffix); the CFL 0.2 reference run has no suffix
CFLS = [(0.20, ""), (0.05, "_cfl05"), (0.0125, "_cfl0125")]


def fit_power(E, dts):
    """err = a dt^p + b from >= 3 points on a constant-ratio dt ladder."""
    r = dts[0] / dts[1]
    num, den = E[0] - E[1], E[1] - E[2]
    if num <= 0 or den <= 0:
        return float("nan"), float("nan"), float("nan")
    p = np.log(num / den) / np.log(r)
    a_dt1 = num / (1.0 - r ** -p)
    return p, a_dt1, E[0] - a_dt1


if __name__ == "__main__":
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    stems = sys.argv[1:] or ["ideal_N64", "murnaghan_N64_corr", "tillotson_N64"]
    for stem in stems:
        eos, n, corr = parse_stem(stem)
        vals, cfls = [], []
        for cfl, sfx in CFLS:
            outdir = "output_" + stem + sfx
            if not glob.glob(outdir + "/*.h5"):
                continue
            vals.append(diagnostics(outdir, eos, n, corr))
            cfls.append(cfl)
        if len(vals) < 2:
            print("{:22s} needs at least two CFL runs (have {})".format(stem, len(vals)))
            continue
        dts = [v["dt"] for v in vals]
        print("\n=== {}   t = {:.4f}".format(stem, vals[0]["t"]))
        print("    CFL " + "  ".join("{:.4f}".format(c) for c in cfls)
              + "   dt " + " ".join("{:.4e}".format(d) for d in dts)
              + "   ratios " + " ".join("{:.3f}".format(a / b) for a, b in zip(dts, dts[1:])))
        head = "  {:14s}".format("metric") + "".join("{:>12s}".format("CFL {:g}".format(c)) for c in cfls)
        print(head + "{:>7s}{:>13s}{:>7s}{:>13s}".format("p", "time part", "%", "floor b"))
        for key, label in METRICS:
            E = [v[key] for v in vals]
            row = "  {:14s}".format(label) + "".join("{:12.4e}".format(e) for e in E)
            if len(E) >= 3:
                p, a_dt1, b = fit_power(E, dts)
                row += "{:7.2f}{:13.4e}{:6.0f}%{:13.4e}".format(p, a_dt1, 100 * a_dt1 / E[0], b)
            else:
                # two points only: assume p = 1, and say so
                r = dts[0] / dts[1]
                a_dt1 = (E[0] - E[1]) / (1.0 - 1.0 / r)
                row += "{:>7s}{:13.4e}{:6.0f}%{:13.4e}".format("(1)", a_dt1, 100 * a_dt1 / E[0], E[0] - a_dt1)
            print(row)
        print("  {:14s}".format("a_rho") + "".join("{:12.6f}".format(v["a_rho"]) for v in vals))
        print("  {:14s}".format("phi [rad]") + "".join("{:+12.4e}".format(v["phi"]) for v in vals))
        if len(vals) < 3:
            print("  p assumed = 1 (only two CFL points); run ./cflTest.sh -c .0125 to measure it")
