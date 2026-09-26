#!/usr/bin/env python3
"""tqparse.py run [...]: integrated L_z per torque channel (TQDIAG intLz lines) at a few times + at the end, relative to L_z0."""
import sys, re, os, numpy as np, h5py, glob
os.chdir(os.path.dirname(os.path.abspath(__file__)))
def parse(run):
    s = open(f"runs/{run}/run.log", errors="ignore").read()
    t = [float(x) for x in re.findall(r"\[TQDIAG\] t=([0-9.eE+-]+)", s)]
    I = re.findall(r"intLz: hllc=(\S+) dev=(\S+) wtT=(\S+) wtR=(\S+) \| total=(\S+)", s)
    return np.array(t), np.array(I, dtype=float)
if __name__ == "__main__":
    for run in sys.argv[1:]:
        t, I = parse(run)
        L0 = None
        f0 = sorted(glob.glob(f"runs/{run}/output/*.h5"))[0]
        with h5py.File(f0) as h:
            x = h["x"][()]; v = h["v"][()]; m = h["m"][()]; L0 = np.sum(m*(x[:, 0]*v[:, 1]-x[:, 1]*v[:, 0]))
        print(run, "  (intLz / L0)   hllc      dev       wtT       total")
        for tt in (5, 10, 15, 20, 25, 30, t[-1]):
            k = np.searchsorted(t, tt) if tt < t[-1] else len(t)-1
            print(f"   t={t[k]:6.2f}   {I[k,0]/L0:+.2e} {I[k,1]/L0:+.2e} {I[k,2]/L0:+.2e} {I[k,4]/L0:+.2e}")
