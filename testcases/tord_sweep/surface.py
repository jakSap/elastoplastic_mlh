#!/usr/bin/env python3
"""Fit the numerical sound speed as a surface c_num(dt, h) instead of two ladders.

Sweep 2 left one statement standing: the dt-dependent part of c_num - 1 is first
order in every variant, while its dt -> 0 floor is +3.83e-3 identically across all
of them.  The two-term ladder fit cannot separate those cleanly from a single
resolution -- it attributes everything that does not scale with dt to `b`.  With
three resolutions (n = 64, 128, 256) the model

    c_num(dt, h) - 1 = A dt^p + B h^q

is overdetermined: p and q are scanned, A and B are linear, and the residual says
whether the two-term ansatz is even the right shape.

Usage:  ROOT=~/Documents/tord_sweep python3 surface.py [tag ...]
"""
import json, os, sys
import numpy as np

ROOT = os.path.expanduser(os.environ.get("ROOT", "~/Documents/tord_sweep"))
SW = os.path.join(ROOT, "testcases", "soundwave_2d")
d = json.load(open(os.path.join(SW, "summary.json")))
TAGS = sys.argv[1:] or ["base", "kdk"]


def gather(tag):
    dt, h, c, damp = [], [], [], []
    for n in (64, 128, 256):
        key = "%s_N%d" % (tag, n)
        if key not in d:
            continue
        for r in d[key]["rows"]:
            dt.append(r["dt"]); h.append(1.0 / n)
            c.append(r["c_num"] - 1.0); damp.append(r["amp_def"])
    return map(np.asarray, (dt, h, c, damp))


def fit_surface(dt, h, y):
    """y = A dt^p + B h^q; p, q scanned, A and B linear."""
    best = None
    for p in np.arange(0.4, 3.01, 0.02):
        for q in np.arange(0.4, 3.01, 0.02):
            X = np.vstack([dt ** p, h ** q]).T
            coef, *_ = np.linalg.lstsq(X, y, rcond=None)
            r = y - X @ coef
            ss = float(r @ r)
            if best is None or ss < best[0]:
                best = (ss, p, q, coef[0], coef[1], r)
    ss, p, q, A, B, r = best
    return dict(p=float(p), q=float(q), A=float(A), B=float(B),
                rms_rel=float(np.sqrt(ss / len(y)) / max(abs(y).max(), 1e-300)))


for tag in TAGS:
    dt, h, c, damp = gather(tag)
    if len(dt) == 0:
        print("no data for", tag); continue
    ns = sorted({int(round(1 / x)) for x in h})
    print("\n=== %s :: %d points, n = %s ===" % (tag, len(dt), ns))
    for name, y in (("c_num - 1", c), ("damping 1 - a_rho", damp)):
        f = fit_surface(dt, h, y)
        print("  %-18s  A dt^%.2f + B h^%.2f   A = %+.4g  B = %+.4g   rms = %.1e"
              % (name, f["p"], f["q"], f["A"], f["B"], f["rms_rel"]))
