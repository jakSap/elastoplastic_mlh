#!/usr/bin/env python3
"""
Side-by-side comparison of the 2D iron-impact benchmark:
  left columns : miluphcuda  (../../../miluphcuda_bare/plasticity_impact_2d)
  right columns: newMFM demonstrator (this directory's output/)

Both are the SAME physical setup (iron bullet -> iron block, dx=2.5mm, 1000 m/s,
Tillotson EOS + von Mises 2.5 GPa yield), so frames line up 1:1 in time
(frame N = N microseconds). Particles colored by speed |v| (shared field).
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import h5py

HERE = os.path.dirname(os.path.abspath(__file__))
MILU = os.path.join(HERE, "..", "..", "..", "miluphcuda_bare", "plasticity_impact_2d")
DEMO = os.path.join(HERE, "output")

TIMES_US = [0, 4, 8, 10]   # microseconds (rows); faithful window before the
                           # demonstrator's residual free-surface tensile
                           # instability ejects a few particles (t>~12 us)


def load_milu(fr):
    f = h5py.File(os.path.join(MILU, f"impact.{fr:04d}.h5"), "r")
    x = f["x"][:]
    v = f["v"][:]
    return x[:, 0], x[:, 1], np.sqrt((v**2).sum(1))


def load_demo(fr):
    f = h5py.File(os.path.join(DEMO, f"{fr:06d}.h5"), "r")
    x = f["x"][:]
    v = f["v"][:]
    return x[:, 0], x[:, 1], np.sqrt((v**2).sum(1))


def main():
    # keep only frames that exist on both sides
    times = []
    for t in TIMES_US:
        if (os.path.exists(os.path.join(MILU, f"impact.{t:04d}.h5"))
                and os.path.exists(os.path.join(DEMO, f"{t:06d}.h5"))):
            times.append(t)
    n = len(times)
    fig, axes = plt.subplots(n, 2, figsize=(8, 3.4 * n), dpi=120, squeeze=False)
    vmax = 1000.0
    for row, t in enumerate(times):
        for col, (loader, title) in enumerate(
                [(load_milu, "miluphcuda"), (load_demo, "newMFM demonstrator")]):
            ax = axes[row][col]
            xx, yy, sp = loader(t)
            sc = ax.scatter(xx, yy, c=sp, s=4, cmap="viridis", vmin=0, vmax=vmax)
            ax.set_aspect("equal")
            ax.set_xlim(-0.06, 0.06); ax.set_ylim(-0.09, 0.03)
            ax.set_title(f"{title}  t={t} us", fontsize=9)
            if col == 1:
                fig.colorbar(sc, ax=ax, fraction=0.046, label="|v| [m/s]")
    fig.suptitle("2D iron impact: miluphcuda vs newMFM demonstrator (MFM)", y=1.0)
    fig.tight_layout()
    out = os.path.join(HERE, "compare_milupch.png")
    fig.savefig(out, bbox_inches="tight")
    print("Saved", out)


if __name__ == "__main__":
    main()
