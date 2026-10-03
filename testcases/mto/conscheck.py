#!/usr/bin/env python3
"""Conservation check of a two-body run: momentum of each body (index halves), total energy, and
L_z, over time.   ./conscheck.py runs/<cand>/<test>/<rung>"""
import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lib.analysis import snapshots, load
D = sys.argv[1]
sn = snapshots(os.path.join(D, "output"))
E0 = None
for t, p in sn[:: max(1, len(sn)//8)] + [sn[-1]]:
    d = load(p, ("x", "v", "m"))
    x, v, m = d["x"], d["v"], d["m"]; h = len(m)//2
    pA = (m[:h, None]*v[:h]).sum(0); pB = (m[h:, None]*v[h:]).sum(0)
    Lz = float(np.sum(m*(x[:, 0]*v[:, 1] - x[:, 1]*v[:, 0])))
    E0 = d["E"] if E0 is None else E0
    print("t={:6.2f}  pA=({:+.3e},{:+.3e})  pA+pB=({:+.1e},{:+.1e})  Lz={:+.3e}  dE/E0={:+.2e}".format(
        t, pA[0], pA[1], pA[0]+pB[0], pA[1]+pB[1], Lz, (d["E"]-E0)/abs(E0)))
