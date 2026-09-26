#!/usr/bin/env python3
"""analyze.py run_tag [...]: L_z(t) drift, dt stats, and the kick/drift prediction w^2/2 * sum dt^2 (DRIFT_KDK only).
Prints per run: steps, mean dt, dL/L at the end and at t_check, rho_mean drift, u_min, negE, r_max/R."""
import sys, re, glob, json, numpy as np, h5py, os
os.chdir(os.path.dirname(os.path.abspath(__file__)))
def load(run):
    fs = sorted(glob.glob(f"runs/{run}/output/*.h5"))
    T, L, rho, umin, rmax, KE = [], [], [], [], [], []
    for f in fs:
        with h5py.File(f) as h:
            x = h["x"][()]; v = h["v"][()]; m = h["m"][()]
            T.append(h["time"][()][0]); L.append(np.sum(m*(x[:, 0]*v[:, 1]-x[:, 1]*v[:, 0])))
            rho.append(np.sum(m*h["rho"][()])/m.sum()); umin.append(h["u"][()].min())
            c = (m[:, None]*x).sum(0)/m.sum(); rmax.append(np.max(np.linalg.norm(x-c, axis=1)))
            KE.append(0.5*np.sum(m*(v**2).sum(1)))
    return np.array(T), np.array(L), np.array(rho), np.array(umin), np.array(rmax), np.array(KE)
def dts(run):
    s = open(f"runs/{run}/run.log", errors="ignore").read()
    return np.array([float(x) for x in re.findall(r"dt = ([0-9.eE+-]+) selected", s)]), s.count("Negative internal")
if __name__ == "__main__":
    for run in sys.argv[1:]:
        T, L, rho, umin, rmax, KE = load(run); d, neg = dts(run)
        w = float(h5py.File("ics/"+re.search(r"initFile ics/(\S+)", open(f"runs/{run}/config.info").read()).group(1)).attrs["omega"])
        ec = open(f"runs/{run}/exit_code").read().strip() if os.path.exists(f"runs/{run}/exit_code") else "run"
        rel = (L-L[0])/L[0]
        print(f"{run:24s} ex={ec:3s} t={T[-1]:6.2f} steps={len(d):6d} <dt>={d.mean() if len(d) else 0:.5f} "
              f"dL/L={rel[-1]:+.3e} pred_kdk={w*w/2*np.sum(d**2):+.3e} max|dL/L|={np.abs(rel).max():.3e} "
              f"drho={rho[-1]-rho[0]:+.2e} umin={umin.min():.3f} negE={neg} rmax={rmax.max():.4f}")
