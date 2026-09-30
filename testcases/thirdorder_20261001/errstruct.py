#!/usr/bin/env python3
"""Error structure of a soundwave run: L1 over time, and the error split into the fitted wave mode,
the rest of the y-averaged profile (other kx) and the y-dependent part. Usage: ./errstruct.py run_dir n"""
import h5py, numpy as np, glob, sys, os
k = 2*np.pi
d, n = sys.argv[1], int(sys.argv[2])
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from analysis import amplitude
A = amplitude(d + "/output")
for fn in sorted(glob.glob(d + "/output/*.h5")):
    with h5py.File(fn) as h:
        t = float(np.ravel(h['time'])[0]); x = h['x'][:]; v = h['v'][:]; rho = h['rho'][:]; u = h['u'][:]
    e = (rho - rho.mean()*(1 + A*np.sin(k*(x[:, 0] - t))))/A
    col = np.floor(x[:, 0]*n).astype(int) % n
    prof = np.bincount(col, e, n)/np.bincount(col, minlength=n)
    ydep = e - prof[col]
    F = np.fft.rfft(prof)/n
    top = np.argsort(-np.abs(F))[:4]
    print("t=%.3f L1 %.3e | profile L1 %.3e  (|F| kx: %s) | y-dependent L1 %.3e | max|vy|/A %.2e" % (
        t, np.mean(np.abs(e)), np.mean(np.abs(prof)), " ".join("%d:%.1e" % (kk, 2*abs(F[kk])) for kk in top),
        np.mean(np.abs(ydep)), np.abs(v[:, 1]).max()/A))
