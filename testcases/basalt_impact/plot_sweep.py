#!/usr/bin/env python3
"""Compare the basalt-impact failure-model sweep. Produces two figures in the
sweep directory:
  sweep_stress.png : final-time svm = sqrt(2 J2) for all 7 combos
  sweep_damage.png : final-time damageTotal for the 3 fragmentation combos
"""
import glob, os
import numpy as np
import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = os.path.dirname(os.path.abspath(__file__))
NAMES = ["01_elastic", "02_vonmises", "03_mohrcoulomb", "04_druckerprager",
         "05_fragmentation", "06_fragmentation_damageS", "07_collins_damage"]


def last(d):
    fs = sorted(glob.glob(os.path.join(BASE, d, "output", "*.h5")))
    return fs[-1] if fs else None


def svm(f):
    Sxx, Syy, Sxy = f["Sxx"][:], f["Syy"][:], f["Sxy"][:]
    return np.sqrt(2 * (0.5 * (Sxx**2 + Syy**2 + 2 * Sxy**2)))


# --- stress figure (all 7) ---
fig, axes = plt.subplots(2, 4, figsize=(18, 9))
axes = axes.ravel()
for ax, name in zip(axes, NAMES):
    fn = last(name)
    if fn is None:
        ax.set_visible(False); continue
    f = h5py.File(fn, "r"); x = f["x"][:]; s = svm(f)
    t = float(np.ravel(f["time"])[0])
    sc = ax.scatter(x[:, 0], x[:, 1], c=s, s=3, cmap="viridis", vmin=0, vmax=0.05)
    ax.set_aspect("equal"); ax.set_title(f"{name}\nsvm  t={t:.1f}", fontsize=10)
    plt.colorbar(sc, ax=ax, fraction=0.046)
axes[-1].set_visible(False)
fig.suptitle("Basalt impact sweep — deviatoric stress sqrt(2 J2), final time", fontsize=13)
fig.tight_layout()
fig.savefig(os.path.join(BASE, "sweep_stress.png"), dpi=110)
print("wrote sweep_stress.png")

# --- damage figure (3 fragmentation) ---
FRAG = ["05_fragmentation", "06_fragmentation_damageS", "07_collins_damage"]
fig, axes = plt.subplots(1, 3, figsize=(16, 5.5))
for ax, name in zip(axes, FRAG):
    fn = last(name)
    f = h5py.File(fn, "r"); x = f["x"][:]; dt = f["damageTotal"][:]
    t = float(np.ravel(f["time"])[0])
    sc = ax.scatter(x[:, 0], x[:, 1], c=dt, s=3, cmap="inferno", vmin=0, vmax=1)
    ax.set_aspect("equal")
    ax.set_title(f"{name}\ndamageTotal t={t:.1f}  mean={dt.mean():.2f}", fontsize=10)
    plt.colorbar(sc, ax=ax, fraction=0.046)
fig.suptitle("Basalt impact sweep — total damage, final time", fontsize=13)
fig.tight_layout()
fig.savefig(os.path.join(BASE, "sweep_damage.png"), dpi=110)
print("wrote sweep_damage.png")
