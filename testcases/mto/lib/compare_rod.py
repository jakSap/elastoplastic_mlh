#!/usr/bin/env python3
"""Side-by-side frames: demonstrator P1 run vs GIZMO control, same times.

    python3 lib/compare_rod.py CAND RUNG [-n 4] [--field rho]
"""
import argparse, os, sys
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from analysis import snapshots, load

ap = argparse.ArgumentParser(); ap.add_argument("cand"); ap.add_argument("rung")
ap.add_argument("-n", type=int, default=4); ap.add_argument("--field", default="rho")
a = ap.parse_args()
root = os.path.join(HERE, "..")
runs = [(a.cand, snapshots(os.path.join(root, "runs", a.cand, "P1", a.rung, "output"))),
        ("GIZMO", snapshots(os.path.join(root, "gizmo", "runs", a.rung, "conv")))]
tend = min(s[-1][0] for _, s in runs)
times = [tend * i / (a.n - 1) for i in range(a.n)]
fig, axs = plt.subplots(2, a.n, figsize=(3.0 * a.n, 4.6), squeeze=False)
for row, (name, snaps) in enumerate(runs):
    for col, t in enumerate(times):
        tt, p = min(snaps, key=lambda s: abs(s[0] - t))
        s = load(p, ("x", a.field)); x = s["x"]; ax = axs[row][col]
        sc = ax.scatter(x[:, 0], x[:, 1], c=s[a.field], s=max(0.3, 2500 / len(x)), cmap="turbo",
                        vmin=0.8, vmax=1.2, linewidths=0)
        ax.set_aspect("equal"); ax.set_xlim(-0.8, 0.8); ax.set_ylim(-0.5, 0.5)
        ax.set_title("{} t={:.2f}".format(name, tt), fontsize=8); ax.tick_params(labelsize=6)
fig.colorbar(sc, ax=axs.ravel().tolist(), shrink=0.7, label=a.field)
png = os.path.join(root, "runs", a.cand, "P1", a.rung, "vs_gizmo_{}.png".format(a.field))
fig.savefig(png, dpi=130, bbox_inches="tight"); print(png)
