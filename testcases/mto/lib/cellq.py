#!/usr/bin/env python3
"""Particle-arrangement quality under large deformation (P1 Taylor rods, 10-05).

A Lagrangian particle set that moves exactly with the material inherits its deformation: a square
lattice stretched by (l1, l2) has rectangular Voronoi cells of aspect l1/l2. GIZMO does not do that
beyond moderate strain; its particles leave the affine map and re-pack into regular cells (P1 w16:
interior cells with material distortion > 1 keep G/G_sq - 1 ~ 0.2 where the affine value is ~1.4).
This module measures how much of the material deformation the particle cells inherit.

Per particle i (body = rod, particle order = IC order in both codes):

  delta_i = ln(l1/l2)        material distortion: principal stretches of the least-squares deformation
                             gradient F_i = argmin sum_j |x_ij - F X_ij|^2 over the IC lattice neighbours
                             (|X_ij| <= 1.5 d0, same body)
  g_i = G(V_i)/G_sq - 1      cell irregularity: G(V) = int_V |x - x_i|^2 dA / (2 |V|^2), the dimensionless
                             second moment of the particle's Voronoi cell among its own body
                             (Conway & Sloane; G_sq = 1/12 square, 5/(36 sqrt 3) hexagon -> g >= -0.038).
                             Shape and centring in one number: an off-centre particle (a pair) raises G too.
                             Cells are clipped to |x - x_i| <= 2 Delta_i, Delta_i = sqrt(m_i/rho_i); a cell
                             that touches the clip disc is a free-surface cell and is left out.

  affine motion of the IC lattice  =>  g_i = cosh(delta_i) - 1  exactly.

Run score at t_end, over S = {interior i : delta_i > 1} (cells the affine map would stretch beyond aspect e):

  kappa = sum_S g_i / sum_S (cosh(delta_i) - 1)      inheritance: 1 = particles carry the deformation
                                                     (deformed lattice), 0 = cells fully re-regularised
  g99   = 99th percentile of g over all interior particles (outcome, no deformation reference)
"""
import numpy as np
from scipy.spatial import cKDTree

G_SQ = 1.0 / 12.0
DELTA_STAR = 1.0
_NPOLY = 48
_ang = np.linspace(0, 2 * np.pi, _NPOLY, endpoint=False)
_CIRC = np.column_stack([np.cos(_ang), np.sin(_ang)])


def _clip(P, n, c):
    """Part of the convex polygon P with n.x <= c (Sutherland-Hodgman, one half-plane)."""
    s = P @ n - c
    inside = s <= 0
    if inside.all():
        return P
    out = []
    k = len(P)
    for a in range(k):
        b = (a + 1) % k
        if inside[a]:
            out.append(P[a])
        if inside[a] != inside[b]:
            t = s[a] / (s[a] - s[b])
            out.append(P[a] + t * (P[b] - P[a]))
    return np.array(out) if len(out) >= 3 else None


def cell_g(x, m, rho, body):
    """g_i = G(V_i)/G_sq - 1 per particle; NaN for free-surface cells."""
    g = np.full(len(x), np.nan)
    D = np.sqrt(m / rho)
    for b in np.unique(body):
        idx = np.where(body == b)[0]
        xb = x[idx]
        tr = cKDTree(xb)
        dist, nb = tr.query(xb, k=min(40, len(idx)))
        for k, i in enumerate(idx):
            R = 2.0 * D[i]
            P = R * _CIRC
            rmax = R
            for dj, j in zip(dist[k, 1:], nb[k, 1:]):
                if dj > 2.0 * rmax:           # no farther particle can cut the cell
                    break
                if dj < 1e-15:
                    continue
                r = xb[j] - xb[k]
                P = _clip(P, r, 0.5 * dj * dj)
                if P is None:
                    break
                rmax = np.sqrt((P * P).sum(1).max())
            if P is None or rmax > 0.999 * R:
                continue
            q = np.roll(P, -1, axis=0)
            cr = P[:, 0] * q[:, 1] - P[:, 1] * q[:, 0]
            A = 0.5 * cr.sum()
            I = (cr * ((P * P).sum(1) + (q * q).sum(1) + (P * q).sum(1))).sum() / 12.0
            g[i] = I / (2 * A * A) / G_SQ - 1.0
    return g


def ic_pairs(X0, body, d0):
    pairs = cKDTree(X0).query_pairs(1.5 * d0, output_type="ndarray")
    return pairs[body[pairs[:, 0]] == body[pairs[:, 1]]]


def distortion(x, X0, pairs):
    """delta_i = ln(l1/l2) of the least-squares deformation gradient over the IC neighbours."""
    n = len(x)
    i, j = pairs[:, 0], pairs[:, 1]
    dX, dx = X0[j] - X0[i], x[j] - x[i]
    A = np.zeros((n, 2, 2)); B = np.zeros((n, 2, 2))
    o, gg = np.einsum("ka,kb->kab", dx, dX), np.einsum("ka,kb->kab", dX, dX)
    for idx in (i, j):
        np.add.at(A, idx, o); np.add.at(B, idx, gg)
    F = A @ np.linalg.inv(B)
    ev = np.linalg.eigvalsh(F @ np.transpose(F, (0, 2, 1)))
    return 0.5 * np.log(ev[:, 1] / ev[:, 0])


def score(x, m, rho, X0, body, d0, pairs=None):
    """kappa, g99 and the sizes of the sets for one snapshot."""
    if pairs is None:
        pairs = ic_pairs(X0, body, d0)
    g = cell_g(x, m, rho, body)
    d = distortion(x, X0, pairs)
    ok = ~np.isnan(g)
    S = ok & (d > DELTA_STAR)
    aff = np.cosh(d[S]) - 1.0
    return dict(kappa=float(g[S].sum() / aff.sum()) if S.any() else float("nan"),
                g99=float(np.percentile(g[ok], 99)), g_mean=float(g[ok].mean()),
                n_S=int(S.sum()), n_interior=int(ok.sum()), g=g, delta=d)


def rank(rungs, root, jobs=4):
    """kappa/g99 at t_end for every candidate's P1 run on the given rungs, plus the GIZMO controls."""
    import glob, json, os
    from concurrent.futures import ProcessPoolExecutor
    tasks = []
    for rung in rungs:
        tasks.append(("GIZMO", rung, os.path.join(root, "gizmo", "runs", rung, "conv"), None))
        for D in sorted(glob.glob(os.path.join(root, "runs", "*", "P1", rung))):
            tasks.append((D.split(os.sep)[-3], rung, os.path.join(D, "output"), D))
    with ProcessPoolExecutor(jobs) as ex:
        return [r for r in ex.map(_one, tasks) if r]


def _one(task):
    import json, os
    from analysis import snapshots, load
    cand, rung, out, D = task
    s = snapshots(out)
    if len(s) < 2:
        return None
    X0 = load(s[0][1], ("x",))["x"]
    body = (X0[:, 0] > 0).astype(int)
    meta = json.load(open(os.path.join(D, "meta.json"))) if D else {}
    res = json.load(open(os.path.join(D, "result.json"))) if D and os.path.exists(os.path.join(D, "result.json")) else {}
    import tests as T
    r = T.rung("P1", rung)
    d0, t_end = r["d"], r["t_end"]
    f = load(s[-1][1], ("x", "m", "rho"))
    if not np.all(np.isfinite(f["x"])):
        return dict(cand=cand, rung=rung, t=f["t"], t_end=t_end, kappa=float("nan"), g99=float("nan"),
                    g_mean=float("nan"), n_S=0, verdict=res.get("verdict", "?"), final_clumps=res.get("final_clumps"))
    sc = score(f["x"], f["m"], f["rho"], X0, body, d0)
    return dict(cand=cand, rung=rung, t=f["t"], t_end=t_end, kappa=sc["kappa"], g99=sc["g99"], g_mean=sc["g_mean"],
                n_S=sc["n_S"], verdict=res.get("verdict", "?"), final_clumps=res.get("final_clumps"))


if __name__ == "__main__":
    import argparse, json, os, sys
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("rungs", nargs="*", help="rank all P1 runs on these rungs (none: self-test)")
    ap.add_argument("-j", type=int, default=4)
    ap.add_argument("-o", help="write the table as JSON")
    a = ap.parse_args()
    if not a.rungs:
        for a_ in (1.0, 2.0, 5.0):   # stretched square lattice: interior g = cosh(ln a) - 1
            L = np.mgrid[0:21, 0:21].reshape(2, -1).T.astype(float)
            x = L * np.array([np.sqrt(a_), 1 / np.sqrt(a_)])
            g = cell_g(x, np.ones(len(x)), np.ones(len(x)), np.zeros(len(x), int))
            c = (L[:, 0] == 10) & (L[:, 1] == 10)
            print("a = %g: g %.4f, cosh(ln a) - 1 = %.4f, surface cells %d" % (a_, g[c][0], np.cosh(np.log(a_)) - 1, np.isnan(g).sum()))
        sys.exit()
    rows = rank(a.rungs, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."), a.j)
    if a.o:
        json.dump(rows, open(a.o, "w"), indent=1)
    for rung in a.rungs:
        print("\n%s  (kappa: share of the affine cell distortion kept where delta > %g; g99: 99th pct cell irregularity)" % (rung, DELTA_STAR))
        rr = sorted([r for r in rows if r["rung"] == rung], key=lambda r: (np.nan_to_num(r["kappa"], nan=9), r["g99"]))
        for r in rr:
            print("  %-34s kappa %6.3f  g99 %6.3f  n_S %4d  t %4.2f%s  %s clumps %s" % (
                r["cand"], r["kappa"], r["g99"], r["n_S"], r["t"], "" if r["t"] >= r["t_end"] - 1e-6 else " (unfinished)",
                r["verdict"], r["final_clumps"]))
