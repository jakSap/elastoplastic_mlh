#!/usr/bin/env python3
"""Mirror-symmetry growth in the first steps of P1: rerun a candidate's P1 rung with dumps every
DT up to T, print the worst mirror mismatch of x and v per dump.   ./symtest.py CAND [w16] [T] [DT]"""
import os, re, sys, subprocess, numpy as np, h5py, glob, shutil
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lib.analysis import snapshots, load, sym_partners
cand = sys.argv[1]; rung = sys.argv[2] if len(sys.argv) > 2 else "w16"
T = float(sys.argv[3]) if len(sys.argv) > 3 else 0.04; DT = float(sys.argv[4]) if len(sys.argv) > 4 else 0.002
run = "runs/{}/P1/{}".format(cand, rung)
if not os.path.exists(run + "/config.info"):          # flag-only variant: borrow the reference config
    run = "runs/pm_obj0/P1/{}".format(rung)
    subprocess.run(["./mto.py", "build", cand, "-v", "plastic"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
D = "runs/symtest/{}_{}".format(cand, rung); shutil.rmtree(D, ignore_errors=True); os.makedirs(D + "/output")
cfg = open(run + "/config.info").read()
cfg = re.sub(r"(?m)^outDir .*$", "outDir   " + os.path.abspath(D + "/output"), cfg)
cfg = re.sub(r"(?m)^timeStep .*$", "timeStep {}".format(DT), cfg)
cfg = re.sub(r"(?m)^timeEnd .*$", "timeEnd  {}".format(T), cfg)
open(D + "/config.info", "w").write(cfg)
binary = os.path.abspath("builds/{}/plastic/bin/mlh".format(cand))
subprocess.run([binary, "-c", "config.info"], cwd=D, stdout=open(D + "/run.log", "w"), stderr=subprocess.STDOUT)
d0 = {"w10": 1/40., "w16": 1/64., "w24": 1/96.}[rung]
with h5py.File("ics/P1_{}.h5".format(rung)) as f: x0 = f["x"][:]
maps = sym_partners(x0, "P1")
for t, p in snapshots(D + "/output"):
    s = load(p, ("x", "v")); x, v = s["x"], s["v"]
    ex = max(np.abs(x - x[j]*[sx, sy]).max() for (sx, sy), j in maps) / d0
    ev = max(np.abs(v - v[j]*[sx, sy]).max() for (sx, sy), j in maps)
    print("{} t={:.4f}  dx/d0={:.2e}  dv={:.2e}".format(cand, t, ex, ev))
