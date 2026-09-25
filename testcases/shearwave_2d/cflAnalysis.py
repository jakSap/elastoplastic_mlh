#!/usr/bin/env python3
"""Compare the CFL 0.2 and CFL 0.05 shear-wave runs at 3 passings.

Same decomposition as testcases/soundwave_2d_elastic/cflAnalysis.py, applied to the mode
whose restoring force is entirely deviatoric. With defect(dt) = a dt + b and dt cut by 4,

    a dt_full = (defect_0.2 - defect_0.05)/(1 - 1/4),   b = defect - a dt_full,

so a ratio of 4 is pure O(dt) and a ratio of 1 is pure O(h).

The quantity that plays the role of the P-wave's stress-slaving defect is the
travelling-wave coupling

    max |S_xy + rho c_s v_y| / (rho_0 c_s^2 A),

which is the relation the IC is built on and which a phase offset between the stress
integrator and the momentum/density update breaks. See cflTest.sh.
"""

import glob
import sys

import h5py
import numpy as np


def ic_name(n):
    return "shearwave_N{}_corr.h5".format(n)


def at_3T(outdir, T):
    best = None
    for f in glob.glob(outdir + "/*.h5"):
        with h5py.File(f, "r") as d:
            t = float(np.ravel(d["time"])[0])
        if best is None or abs(t - 3 * T) < abs(best[1] - 3 * T):
            best = (f, t)
    return best


def diagnostics(path, t, md, n):
    A, c_s, rho0 = md["amplitude"], md["c_wave"], md["rho_bar"]
    S_SCALE, V_SCALE = rho0 * c_s * c_s * A, c_s * A
    k = 2 * np.pi / md["wavelength"]
    with h5py.File(path, "r") as d:
        x, vy = d["x"][:, 0], d["v"][:, 1]
        rho, sxy = d["rho"][:], d["Sxy"][:]
    phase = k * c_s * t
    wave = np.sin(k * x - phase)
    col = np.floor(x * n).astype(int) % n
    cnt = np.maximum(np.bincount(col, minlength=n), 1)
    xc = (np.arange(n) + 0.5) / n
    coupling = (sxy + rho * c_s * vy) / S_SCALE
    prof = np.bincount(col, weights=coupling, minlength=n) / cnt
    pert = np.bincount(col, weights=vy, minlength=n) / cnt / V_SCALE
    M = np.column_stack([np.sin(k * xc - phase), np.cos(k * xc - phase)])
    (p, q), *_ = np.linalg.lstsq(M, pert, rcond=None)
    return dict(coupling=np.abs(prof).max(),
                L1_vy=np.mean(np.abs(vy - V_SCALE * wave)),
                a_vy=float(np.hypot(p, q)), phi=float(np.arctan2(-q, p)))


if __name__ == "__main__":
    ns = [int(a) for a in sys.argv[1:]] or [64]
    print("{:>5s} {:>9s} {:>14s} {:>12s} {:>10s} {:>12s}".format(
        "n", "CFL", "max|coupling|", "L1_vy", "a_vy", "phi"))
    for n in ns:
        with h5py.File(ic_name(n), "r") as f:
            md = {k: float(f.attrs[k]) for k in
                  ("amplitude", "wavelength", "mu", "c_wave", "period", "rho_bar")}
        vals = {}
        for cfl, outdir in ((0.20, "output_shear_N{}_corr".format(n)),
                            (0.05, "output_shear_N{}_corr_cfl05".format(n))):
            got = at_3T(outdir, md["period"])
            if not got:
                print("  {:>3d} {:>9.2f}  (no snapshots in {})".format(n, cfl, outdir))
                continue
            vals[cfl] = diagnostics(got[0], got[1], md, n)
            d = vals[cfl]
            print("{:>5d} {:>9.2f} {:>14.4e} {:>12.4e} {:>10.6f} {:>+12.4e}".format(
                n, cfl, d["coupling"], d["L1_vy"], d["a_vy"], d["phi"]))
        if 0.20 in vals and 0.05 in vals:
            for key in ("coupling", "L1_vy"):
                hi, lo = vals[0.20][key], vals[0.05][key]
                time_part = (hi - lo) / 0.75
                print("      {:<13s} ratio {:.2f} (4 = pure O(dt));  "
                      "time part {:.4e} ({:.0f}%), space part {:.4e}".format(
                          key, hi / lo, time_part, 100 * time_part / hi, hi - time_part))
