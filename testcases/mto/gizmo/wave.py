#!/usr/bin/env python3
"""GIZMO control for the Tillotson wave ladders W1t / W2 / W3 (does GIZMO converge to the right speed?).

    ./wave.py W1t W2 W3 [-n 32 64 128] [-j 3]     run, convert to the demonstrator layout, analyse

Binary: ~/Documents/GIZMO/build_mto_wave/GIZMO = build_mto_elastic + BOX_PERIODIC, BOX_LONG_Y=0.25
(the [0,1) x [0,0.25) wave box) + MTO_READ_STRESS (read_ic.c/init.c: keep the IC's StressTensor;
stock GIZMO zeroes it, which would split the travelling wave into two). GIZMO's deviator is
trace-free in 2D, so the ICs are the _pm1 ones (P-modulus K + mu) that the demonstrator uses under
STRESS_2D_GIZMO 1. Same IC, same snapshot times, same analysis (lib/analysis.analyse_wave) as the
candidates' runs. Results: runs/<test>_n<n>/result.json.
"""
import argparse
import concurrent.futures as cf
import glob
import json
import os
import subprocess
import sys

import h5py
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))
sys.path.insert(0, HERE)
import tests as T          # noqa: E402
import analysis as A       # noqa: E402
from rod import PARAMS     # noqa: E402

GIZMO = os.path.expanduser("~/Documents/GIZMO/build_mto_wave/GIZMO")


def setup(test, n):
    r = T.rung(test, "n{}".format(n))
    src = os.path.join(HERE, "..", "ics", "{}_n{}_pm1.h5".format(test, n))
    if not os.path.exists(src):
        raise SystemExit("missing " + src)
    D = os.path.join(HERE, "runs", "{}_n{}".format(test, n))
    os.makedirs(os.path.join(D, "raw"), exist_ok=True)
    with h5py.File(src, "r") as f:
        x, v, m, u = f["x"][:], f["v"][:], f["m"][:], f["u"][:]
        Sxx, Sxy, Syy = f["Sxx"][:], f["Sxy"][:], f["Syy"][:]
    N = len(m)
    pos = np.zeros((N, 3)); pos[:, :2] = x
    vel = np.zeros((N, 3)); vel[:, :2] = v
    S9 = np.zeros((N, 9)); S9[:, 0], S9[:, 1], S9[:, 3], S9[:, 4] = Sxx, Sxy, Sxy, Syy
    with h5py.File(os.path.join(D, "ic.hdf5"), "w") as f:
        h = f.create_group("Header")
        npart = np.array([N, 0, 0, 0, 0, 0], dtype=np.int32)
        h.attrs["NumPart_ThisFile"] = npart
        h.attrs["NumPart_Total"] = npart.astype(np.uint32)
        h.attrs["NumPart_Total_HighWord"] = np.zeros(6, dtype=np.uint32)
        h.attrs["MassTable"] = np.zeros(6)
        for k, val in (("Time", 0.0), ("Redshift", 0.0), ("BoxSize", 1.0),
                       ("NumFilesPerSnapshot", 1), ("Omega0", 0.0), ("OmegaLambda", 0.0),
                       ("HubbleParam", 1.0), ("Flag_Sfr", 0), ("Flag_Cooling", 0),
                       ("Flag_StellarAge", 0), ("Flag_Metals", 0), ("Flag_Feedback", 0),
                       ("Flag_DoublePrecision", 1)):
            h.attrs[k] = val
        g = f.create_group("PartType0")
        g.create_dataset("Coordinates", data=pos)
        g.create_dataset("Velocities", data=vel)
        g.create_dataset("Masses", data=m)
        g.create_dataset("InternalEnergy", data=u)
        g.create_dataset("ParticleIDs", data=np.arange(1, N + 1, dtype=np.uint32))
        g.create_dataset("CompositionType", data=np.zeros(N, dtype=np.int32))
        g.create_dataset("StressTensor", data=S9)
    mat = r["materials"][0]
    text = PARAMS.format(tag="{}_n{}".format(test, n), dump=r["dump"], t_end=r["t_end"], box=1.0,
                         mu=mat["mu"], Y0=1e20)
    text = text.replace("OutputDir                   output", "OutputDir                   raw")
    text = text.replace("Taylor rods", "Tillotson wave").replace("rod.py", "wave.py")
    open(os.path.join(D, "wave.params"), "w").write(text)
    return D, r, src


def convert(D):
    out = os.path.join(D, "output")
    os.makedirs(out, exist_ok=True)
    for old in glob.glob(os.path.join(out, "*.h5")):
        os.remove(old)
    for k, snap in enumerate(sorted(glob.glob(os.path.join(D, "raw", "snapshot_*.hdf5")))):
        with h5py.File(snap, "r") as f:
            t = float(f["Header"].attrs["Time"])
            g = f["PartType0"]
            o = np.argsort(g["ParticleIDs"][:])
            x, v = g["Coordinates"][:][o, :2], g["Velocities"][:][o, :2]
            m, rho, u = g["Masses"][:][o], g["Density"][:][o], g["InternalEnergy"][:][o]
        with h5py.File(os.path.join(out, "{:06d}.h5".format(k)), "w") as d:
            d["x"], d["v"], d["m"], d["rho"], d["u"] = x, v, m, rho, u
            d["time"] = np.array([t])
    return out


def run(job):
    test, n = job
    D, r, ic = setup(test, n)
    with open(os.path.join(D, "run.log"), "w") as log:
        rc = subprocess.run([GIZMO, "wave.params"], cwd=D, stdout=log, stderr=subprocess.STDOUT).returncode
    convert(D)
    res = dict(test=test, tag="n{}".format(n), cand="GIZMO", exit=rc)
    try:
        res.update({k: v for k, v in A.analyse_wave(D, dict(ic=ic, tag="n{}".format(n), materials=r["materials"]),
                                                    test).items()})
    except Exception as e:   # noqa: BLE001
        res["error"] = repr(e)
    json.dump(res, open(os.path.join(D, "result.json"), "w"), indent=1, default=float)
    return res


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("tests", nargs="+")
    ap.add_argument("-n", type=int, nargs="+", default=[32, 64, 128])
    ap.add_argument("-j", type=int, default=3)
    a = ap.parse_args()
    jobs = [(t, n) for n in a.n for t in a.tests]
    with cf.ThreadPoolExecutor(a.j) as ex:
        for r in ex.map(run, jobs):
            print("{} {}: exit {} c_err {:+.3e} damp_v {:.3e} {}".format(
                r["test"], r["tag"], r["exit"], r.get("c_err", float("nan")), r.get("damp_v", float("nan")),
                r.get("error", "")), flush=True)


if __name__ == "__main__":
    main()
