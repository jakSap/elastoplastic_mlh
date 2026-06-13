"""Weibull flaw assignment for the Grady-Kipp damage model (Benz & Asphaug 1995).

The demonstrator reads two datasets from the IC HDF5 when FRAGMENTATION is on:

    /numFlaws : int   [N]              number of activation strains per particle
    /flaws    : float [N, MAX_NUM_FLAWS] the activation strains (padded with 0)

`max_flaws` here MUST match MAX_NUM_FLAWS in demonstrator/include/parameter.h.
With the default MAX_NUM_FLAWS = 1 only the single weakest flaw per particle is
stored (a coarse, near-binary damage model); raise it for a gradual Grady-Kipp
crack-growth response.

Assignment (Benz & Asphaug 1995): flaws are handed out one at a time. The j-th
flaw gets activation strain eps_j = (j / (k * V_tot))**(1/m) and is dropped on a
uniformly random particle. Flaws keep being assigned until every particle owns at
least one. k [length^-DIM] and m are the Weibull material parameters.
"""

import numpy as np


def assign_weibull_flaws(volumes, k, m, max_flaws=1, seed=0):
    """Return (numFlaws[N] int32, flaws[N, max_flaws] float64).

    volumes : per-particle volume (mass/density), length N.
    k, m    : Weibull parameters (k in length^-DIM, m dimensionless).
    """
    volumes = np.asarray(volumes, dtype=np.float64)
    N = volumes.size
    V_tot = float(volumes.sum())
    if N == 0:
        return np.zeros(0, dtype=np.int32), np.zeros((0, max_flaws))
    if k <= 0.0 or m <= 0.0:
        raise ValueError("Weibull k and m must be positive")

    rng = np.random.default_rng(seed)
    flaw_lists = [[] for _ in range(N)]
    covered = np.zeros(N, dtype=bool)
    n_covered = 0
    j = 0
    inv_m = 1.0 / m
    while n_covered < N:
        j += 1
        eps = (j / (k * V_tot)) ** inv_m
        p = int(rng.integers(0, N))
        flaw_lists[p].append(eps)
        if not covered[p]:
            covered[p] = True
            n_covered += 1

    numFlaws = np.zeros(N, dtype=np.int32)
    flaws = np.zeros((N, max_flaws), dtype=np.float64)
    for i in range(N):
        fl = sorted(flaw_lists[i])          # ascending: weakest flaw first
        nf = min(len(fl), max_flaws)
        numFlaws[i] = nf
        flaws[i, :nf] = fl[:nf]
    print("Weibull flaws: total assigned =", j,
          "| per-particle min/mean/max =",
          int(numFlaws.min()),
          round(float(numFlaws.mean()), 2),
          int(numFlaws.max()))
    return numFlaws, flaws


def write_flaws(h5f, volumes, k, m, max_flaws=1, seed=0):
    """Convenience: assign flaws and write /numFlaws and /flaws into an open
    h5py.File `h5f`."""
    numFlaws, flaws = assign_weibull_flaws(volumes, k, m, max_flaws, seed)
    h5f.create_dataset("numFlaws", data=numFlaws)
    h5f.create_dataset("flaws", data=flaws)
    return numFlaws, flaws
