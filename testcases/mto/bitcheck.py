#!/usr/bin/env python3
"""Bit-identity check: final snapshot of runs/<A>/<T>/<R> vs runs/<B>/<T>/<R> (all datasets).
    ./bitcheck.py A B P1 w10"""
import sys, glob, h5py, numpy as np
A, B, T, R = sys.argv[1:5]
fa, fb = (sorted(glob.glob("runs/{}/{}/{}/output/*.h5".format(c, T, R)))[-1] for c in (A, B))
ha, hb = h5py.File(fa, "r"), h5py.File(fb, "r")
worst = 0.
for k in ha:
    if k in hb and hasattr(ha[k], "shape"):
        a, b = np.asarray(ha[k]), np.asarray(hb[k])
        if a.shape == b.shape and a.size:
            worst = max(worst, float(np.max(np.abs(a.astype(float) - b.astype(float)))))
print("{} vs {} {} {} ({} / {}): max|diff| over all datasets = {:.3e}".format(A, B, T, R, fa[-9:], fb[-9:], worst))
