#!/usr/bin/env python3
"""Error structure of the transverse (shear) wave: v_y error split into the y-averaged profile (Fourier kx) and the
y-dependent part, plus |v_x| and the density deviation. Usage: ./errstruct_shear.py run_dir n"""
import h5py, numpy as np, glob, sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from analysis import ic_attrs
k = 2*np.pi
d, n = sys.argv[1], int(sys.argv[2])
W = ic_attrs(d + "/output"); A, c = W["A"], W["c"]
for fn in sorted(glob.glob(d + "/output/*.h5")):
    with h5py.File(fn) as h:
        t = float(np.ravel(h['time'])[0]); x = h['x'][:]; v = h['v'][:]; rho = h['rho'][:]
    e = (v[:, 1] - c*A*np.sin(k*(x[:, 0] - c*t)))/(c*A)
    col = np.floor(x[:, 0]*n).astype(int) % n
    prof = np.bincount(col, e, n)/np.bincount(col, minlength=n)
    ydep = e - prof[col]
    F = np.fft.rfft(prof)/n
    top = np.argsort(-np.abs(F))[:5]
    print("t=%.3f L1 %.3e | profile %.3e (kx: %s) | y-dep %.3e | max|vx|/A %.1e | max|rho-1|/A %.1e" % (
        t, np.mean(np.abs(e)), np.mean(np.abs(prof)), " ".join("%d:%.1e" % (kk, 2*abs(F[kk])) for kk in top),
        np.mean(np.abs(ydep)), np.abs(v[:, 0]).max()/A, np.abs(rho - 1).max()/A))
