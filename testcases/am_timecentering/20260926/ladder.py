#!/usr/bin/env python3
"""ladder.py [variants...]: AM error vs resolution for runs/L_<variant>_d<d> (rotating disc, omega 0.2, one rotation).
Metrics: max_t |dL/L| over the rotation, RMS_t, |dL/L| at t_end; plus rho_E/rho_K drift (bulk), r_gyr change, u_min, negE.
Order p = slope of log(metric) vs log(d), least squares over all finished rungs, and pairwise."""
import sys, os, re, glob, numpy as np, h5py
os.chdir(os.path.dirname(os.path.abspath(__file__)))
from analyze import dts
def lz(run):
    fs = sorted(glob.glob(f"runs/{run}/output/*.h5")); T, L, rr, rg, umin = [], [], [], [], []
    for f in fs:
        with h5py.File(f) as h:
            x = h["x"][()]; v = h["v"][()]; m = h["m"][()]
            T.append(h["time"][()][0]); L.append(np.sum(m*(x[:, 0]*v[:, 1]-x[:, 1]*v[:, 0])))
            r = np.linalg.norm(x, axis=1); b = r < 0.7
            if "rhoKernel" in h and h["rhoKernel"][()][b].min() > 0: rr.append(np.mean(h["rho"][()][b]/h["rhoKernel"][()][b]))
            else: rr.append(np.nan)
            rg.append(np.sqrt(np.sum(m*r**2)/m.sum())); umin.append(h["u"][()].min())
    return np.array(T), np.array(L), np.array(rr), np.array(rg), np.array(umin)
variants = sys.argv[1:] or sorted({re.sub(r"_d[0-9.]+$", "", os.path.basename(p))[2:] for p in glob.glob("runs/L_*")})
for v in variants:
    rows = []
    for p in glob.glob(f"runs/L_{v}_d*"):
        if re.sub(r"_d[0-9.]+$", "", os.path.basename(p))[2:] != v: continue
        if not os.path.exists(p + "/exit_code") or not glob.glob(p + "/output/*.h5"): continue
        d = float(re.search(r"_d([0-9.]+)$", p).group(1)); run = os.path.basename(p)
        T, L, rr, rg, umin = lz(run); rel = (L-L[0])/L[0]; dd, neg = dts(run)
        ok = T[-1] >= 30.99
        rows.append((d, np.abs(rel).max(), np.sqrt(np.mean(rel**2)), abs(rel[-1]), rr[-1]-1 if len(rr) else np.nan, rg[-1]/rg[0]-1, umin.min(), neg, ok, len(dd), open(p+"/exit_code").read().strip()))
    rows.sort(reverse=True)
    print(f"## {v}")
    print("%8s %10s %10s %10s %10s %10s %8s %6s %6s %5s" % ("d", "max|dL/L|", "rms", "|end|", "rhoE/rhoK-1", "rgyr-1", "umin", "negE", "steps", "exit"))
    for r in rows:
        print("%8.4f %10.3e %10.3e %10.3e %+10.2e %+10.2e %8.4f %6d %6d %5s%s" % (r[0], r[1], r[2], r[3], r[4], r[5], r[6], r[7], r[9], r[10], "" if r[8] else "  INCOMPLETE"))
    good = [r for r in rows if r[8] and r[10] == "0"]
    if len(good) >= 3:
        ld = np.log([r[0] for r in good])
        for k, name in ((1, "max"), (2, "rms"), (3, "end")):
            y = np.log([r[k] for r in good]); A = np.vstack([ld, np.ones_like(ld)]).T; (p, c), res, *_ = np.linalg.lstsq(A, y, rcond=None)
            r2 = 1 - np.sum((y - A@[p, c])**2)/np.sum((y-y.mean())**2)
            pw = " ".join("%.2f" % ((y[i]-y[i+1])/(ld[i]-ld[i+1])) for i in range(len(y)-1))
            print(f"   order[{name}] p = {p:.2f} (R2 {r2:.3f})   pairwise: {pw}")
    print()
