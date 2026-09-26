#!/usr/bin/env python3
"""tqmap.py run [...]: where is the deviatoric-traction torque? Per snapshot: -sum tqDev (dL/dt from the dev channel,
relative to L0), split into surface band (r > R_edge - 3 d) and bulk; time-mean over the rotation; plus the same for tqWtT, tqHLLC."""
import sys, os, glob, re, numpy as np, h5py
os.chdir(os.path.dirname(os.path.abspath(__file__)))
for run in sys.argv[1:]:
    fs = sorted(glob.glob(f"runs/{run}/output/*.h5"))[2:]
    d = float(re.search(r"_d([0-9.]+)", run).group(1))
    with h5py.File(sorted(glob.glob(f"runs/{run}/output/*.h5"))[0]) as h:
        x = h["x"][()]; v = h["v"][()]; m = h["m"][()]; L0 = np.sum(m*(x[:, 0]*v[:, 1]-x[:, 1]*v[:, 0]))
    acc = {}
    for f in fs:
        with h5py.File(f) as h:
            x = h["x"][()]; r = np.linalg.norm(x - x.mean(0), axis=1); redge = np.percentile(r, 99.5)
            surf = r > redge - 3*d
            for k in ("tqDev", "tqWtT", "tqHLLC"):
                if k not in h: continue
                t = -h[k][()]/L0
                acc.setdefault(k, []).append((t[surf].sum(), t[~surf].sum(), np.abs(t).sum()))
    print(f"{run:24s} d={d}  time-mean dL/dt / L0:  [surface band 3d | bulk | sum|.|]  N={len(x)}")
    for k, a in acc.items():
        a = np.array(a); print(f"   {k:7s} surf {a[:,0].mean():+.2e}  bulk {a[:,1].mean():+.2e}  abs {a[:,2].mean():.2e}")
