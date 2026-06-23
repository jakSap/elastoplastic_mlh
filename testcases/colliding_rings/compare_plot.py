#!/usr/bin/env python3
"""Density plots comparing this codebase (lowres + hires) vs GIZMO reference.

Generates two figures:
  comparison.png       - 3-way side-by-side at t = 0, 10, 20, 30
  gizmo_evolution.png  - GIZMO-only full collision through t = 300
GIZMO uses GADGET-style HDF5 (PartType0/...), demonstrator uses flat keys (x, rho)."""
import h5py, numpy as np, matplotlib, os
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def load_gizmo(path):
    with h5py.File(path, 'r') as f:
        x = f['PartType0/Coordinates'][:]
        rho = f['PartType0/Density'][:]
        t = f['Header'].attrs['Time']
    return x[:, 0], x[:, 1], rho, t


def load_demo(path):
    with h5py.File(path, 'r') as f:
        x = f['x'][:]
        rho = f['rho'][:]
        t = f['time'][0]
    return x[:, 0], x[:, 1], rho, t


def panel(ax, X, Y, rho, t, title, vmin=0.7, vmax=1.1, s=4, xlim=(-5.5, 5.5), ylim=(-5.5, 5.5)):
    sc = ax.scatter(X, Y, c=rho, s=s, vmin=vmin, vmax=vmax, cmap='viridis', edgecolors='none')
    ax.set_aspect('equal')
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_title(f'{title}\nt={t:.2f}, N={len(rho)}', fontsize=9)
    return sc


gizmo_dir = '/home/sappler/Documents/GIZMO/ring_collision/output'
lowres_dir = '/home/sappler/Documents/newMFM/testcases/colliding_rings/no_envelope/output_tillotson'
hires_dir = '/home/sappler/Documents/newMFM/testcases/colliding_rings/no_envelope_hires/output_tillotson'

# ---- Figure 1: 3-way side-by-side at times where all runs have valid data ----
target_times = [0, 10, 20, 30]
fig, axes = plt.subplots(len(target_times), 3, figsize=(13, 3.4 * len(target_times)), dpi=120)
sc_last = None
for i, tt in enumerate(target_times):
    # Lowres: snap index == t (integer)
    p_lo = f'{lowres_dir}/{min(tt, 59):06d}.h5'
    if os.path.exists(p_lo):
        X, Y, rho, t = load_demo(p_lo)
        sc_last = panel(axes[i, 0], X, Y, rho, t, 'demonstrator lowres (deltap=0.05)', s=3)
    # Hires: snap index == t (integer), valid through ~t=20 before instability
    p_hi = f'{hires_dir}/{min(tt, 38):06d}.h5'
    if os.path.exists(p_hi):
        X, Y, rho, t = load_demo(p_hi)
        # for late hires snaps, particles disperse far - keep zoom but warn
        sc_last = panel(axes[i, 1], X, Y, rho, t, 'demonstrator hires (deltap=0.02)', s=1.5)
    # GIZMO: snap N at t = 5*N
    g_idx = min(tt // 5, 60)
    p_g = f'{gizmo_dir}/snapshot_{g_idx:03d}.hdf5'
    if os.path.exists(p_g):
        X, Y, rho, t = load_gizmo(p_g)
        sc_last = panel(axes[i, 2], X, Y, rho, t, 'GIZMO reference (Hopkins)', s=6)
fig.colorbar(sc_last, ax=axes.ravel().tolist(), shrink=0.7, label=r'$\rho$')
plt.suptitle('Colliding Rubber Rings: this codebase (lowres / hires) vs GIZMO reference', y=0.995)
out = '/home/sappler/Documents/newMFM/testcases/colliding_rings/comparison.png'
plt.savefig(out, bbox_inches='tight', dpi=120)
print(f'wrote {out}')
plt.close(fig)

# ---- Figure 2: GIZMO full evolution (only run that captured the actual bounce) ----
g_times = [0, 50, 100, 150, 200, 250, 300]
fig, axes = plt.subplots(2, 4, figsize=(16, 8), dpi=120)
axes = axes.ravel()
sc_last = None
for i, tt in enumerate(g_times):
    g_idx = tt // 5
    p_g = f'{gizmo_dir}/snapshot_{g_idx:03d}.hdf5'
    if os.path.exists(p_g):
        X, Y, rho, t = load_gizmo(p_g)
        sc_last = panel(axes[i], X, Y, rho, t, 'GIZMO', s=8, xlim=(-7, 7), ylim=(-5, 5))
# hide unused last panel
axes[-1].axis('off')
fig.colorbar(sc_last, ax=axes.tolist(), shrink=0.7, label=r'$\rho$')
plt.suptitle('GIZMO reference: full collision evolution (Hopkins canonical params)', y=0.995)
out = '/home/sappler/Documents/newMFM/testcases/colliding_rings/gizmo_evolution.png'
plt.savefig(out, bbox_inches='tight', dpi=120)
print(f'wrote {out}')
plt.close(fig)
