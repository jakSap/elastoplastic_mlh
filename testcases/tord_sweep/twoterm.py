#!/usr/bin/env python3
"""Is the dt error first order or second? Fit E = a*dt + c*dt^2 + b (linear LSQ)
and report whether the LINEAR coefficient survives. If the adiabatic face
correction restored second order, a must collapse to ~0."""
import subprocess, re, numpy as np, sys
txt = subprocess.run([sys.executable, "analyse.py"], capture_output=True, text=True).stdout
blocks, cur = {}, None
for ln in txt.splitlines():
    m = re.match(r"=== (\S+)\s+\((\d+) rungs\)", ln)
    if m: cur = m.group(1); blocks[cur] = []
    elif cur and re.match(r"\s*\d+\s+[\d.]+e", ln):
        f = ln.split(); blocks[cur].append((float(f[1]), float(f[2]), float(f[3])))
print("Model:  E = a*dt + c*dt^2 + b        (14 or 7 rungs, 3 params, overdetermined)")
print("%-14s %-9s %12s %12s %12s   %9s   %s" % ("case", "metric", "a (dt)", "c (dt^2)", "b (floor)", "rms res", "a*dt at CFL0.2-ish"))
for case in ("ideal", "idealAFC", "idealLimN", "idealLimNAFC", "idealN128", "idealAFCN128", "tillotson"):
    if case not in blocks or len(blocks[case]) < 4: continue
    arr = np.array(blocks[case]); dt = arr[:, 0]
    for j, lab in ((1, "1-a_rho"), (2, "L1(rho)")):
        E = arr[:, j]
        X = np.vstack([dt, dt ** 2, np.ones_like(dt)]).T
        coef, *_ = np.linalg.lstsq(X, E, rcond=None)
        r = E - X @ coef
        print("%-14s %-9s %12.4e %12.4e %12.6e   %9.2e   %.4e (%.0f%% of E_max)"
              % (case, lab, coef[0], coef[1], coef[2],
                 np.sqrt((r**2).mean())/abs(E).max(), coef[0]*dt.max(),
                 100*coef[0]*dt.max()/E.max()))
