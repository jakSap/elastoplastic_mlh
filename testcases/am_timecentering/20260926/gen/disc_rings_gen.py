#!/usr/bin/env python3
"""Rigidly rotating elastic disc on CONCENTRIC RINGS -- the same test as disc_gen.py (A.1), but
with a particle layout that forms a smooth refinement family.

Why: disc_gen.py cuts a circle out of a square lattice. The cut differs from rung to rung in a
non-smooth way (L_z/L_z,rigid of the IC: d0.1 0.985, d0.07 1.00006, d0.05 0.994, d0.0177 1.00012),
i.e. the free surface carries an O(d) roughness that is not refined consistently. An order fit
over such rungs measures that roughness as much as the scheme. On concentric rings the surface is
the circle r = R at every resolution.

Layout: ring k = 0..K-1 at r_k = (k + 1/2) d, K = round(R/d), n_k = max(3, round(2 pi r_k / d))
particles, equally spaced, ring k rotated by k * golden angle so no radial spokes form. Every
particle of ring k gets the mass of its share of the annulus [k d, (k+1) d]: m = rho0 pi (2k+1) d^2 / n_k,
so the mass density is rho0 everywhere and the body edge sits at r = K d = R.

Material and state as disc_gen.py (Hopkins-canonical Tillotson, mu = 0.22, u = 1, S = 0).

Usage: disc_rings_gen.py -d 0.025 -w 0.2 -o disc_rings_d0.025_w0.2.h5
"""
import argparse
import h5py
import numpy as np

RHO0, MU, U_BAR = 1.0, 0.22, 1.0
GOLDEN = np.pi*(3. - np.sqrt(5.))


def make(radius, spacing, omega, out):
    K = int(round(radius/spacing)); d = radius/K
    pos, m = [], []
    for k in range(K):
        r = (k + .5)*d; n = max(3, int(round(2*np.pi*r/d)))
        th = k*GOLDEN + 2*np.pi*np.arange(n)/n
        pos.append(np.column_stack([r*np.cos(th), r*np.sin(th)]))
        m.append(np.full(n, RHO0*np.pi*(2*k + 1)*d*d/n))
    pos = np.vstack(pos); m = np.concatenate(m); N = len(m)
    vel = np.column_stack([-omega*pos[:, 1], omega*pos[:, 0]])
    Lz = float(np.sum(m*(pos[:, 0]*vel[:, 1] - pos[:, 1]*vel[:, 0])))
    M = float(m.sum()); Lrig = 0.5*M*radius**2*omega
    with h5py.File(out, "w") as f:
        f.create_dataset("x", data=pos); f.create_dataset("v", data=vel); f.create_dataset("m", data=m)
        f.create_dataset("u", data=np.full(N, U_BAR)); f.create_dataset("materialId", data=np.zeros(N, dtype=np.int8))
        for name in ("Sxx", "Sxy", "Syy"): f.create_dataset(name, data=np.zeros(N))
        for key, val in (("radius", radius), ("spacing", d), ("omega", omega), ("rho0", RHO0), ("mu", MU), ("N", N),
                         ("mass", M), ("Lz", Lz), ("Lz_rigid", Lrig), ("period", 2*np.pi/omega), ("layout", "rings")):
            f.attrs[key] = val
    print("Wrote {}: N = {}, rings = {}, d = {:g}, M = {:.6g} (pi = {:.6g}), L_z/L_rigid = {:.6f}".format(
        out, N, K, d, M, np.pi*RHO0*radius**2, Lz/Lrig))


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-R", "--radius", type=float, default=1.0)
    ap.add_argument("-d", "--spacing", type=float, default=0.025)
    ap.add_argument("-w", "--omega", type=float, default=0.2)
    ap.add_argument("-o", "--out", required=True)
    a = ap.parse_args(); make(a.radius, a.spacing, a.omega, a.out)
