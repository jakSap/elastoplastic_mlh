#!/usr/bin/env python3
"""Ladder analysis for the time-order sweep (2026-09-24).

Metrics identical to mars_fdt_sweep/analyse.py so the numbers are directly
comparable to the 2026-09-23 sweep: dt is exact (= 1/M, one wave period), the
snapshot is picked by stored /time and never by file index, and the damping
1 - a_rho comes from the same column-profile sine fit.

Two fits per run set:
  * E = a dt^p + b   -- p is the headline number; p = 2 is the goal.
  * E = a dt + c dt^2 + b -- a is the linear coefficient; if the order really
    moved to 2, a must go to zero.  This is the fit that settled the question on
    2026-09-24 (a fell 5.9 -> 2.1 and then stood still).

Usage:  ROOT=~/Documents/tord_sweep python3 analyse.py [--json out.json]
"""
import glob, json, os, sys
import numpy as np, h5py

ROOT = os.path.expanduser(os.environ.get("ROOT", "~/Documents/tord_sweep"))
SW = os.path.join(ROOT, "testcases", "soundwave_2d")
sys.path.insert(0, SW)
os.chdir(SW)
from makeMovie import snapshots, column_profile, sine_fit  # noqa: E402

GAMMA = 5.0 / 3.0
TARGET = 1.0  # one wave period


def diag(outdir, ic, n):
    with h5py.File(ic, "r") as f:
        A, lam = float(f.attrs["amplitude"]), float(f.attrs["wavelength"])
    k = 2 * np.pi / lam
    files, times = snapshots(outdir)
    if not files:
        return None
    i = int(np.argmin([abs(t - TARGET) for t in times]))
    path, t = files[i], times[i]
    if abs(t - TARGET) > 1e-6:
        return ("OFFTARGET", t)
    with h5py.File(path, "r") as d:
        x, rho, vx, u, P = (d["x"][:, 0], d["rho"][:], d["v"][:, 0],
                            d["u"][:], d["P"][:])
    mean = rho.mean()
    wave = np.sin(k * (x - t))
    xc, prof, m = column_profile(path, n)
    a_rho, phi = sine_fit(xc, (prof - m) / (m * A), k, t)
    # spurious entropy: amplitude of the P/rho^gamma perturbation
    s = P / rho ** GAMMA
    return dict(t=float(t), a_rho=float(a_rho), amp_def=float(1 - a_rho),
                phi=float(phi),
                L1_rho=float(np.mean(np.abs(rho - mean * (1 + A * wave)))),
                L1_vx=float(np.mean(np.abs(vx - A * wave))),
                nonisen_amp=float((s.max() - s.min()) / (2 * s.mean())),
                c_num=float(1 + phi / (k * t)),
                c_exp=float(np.sqrt(GAMMA * (GAMMA - 1) * u.mean())))


def fit_pab(E, dt):
    """E = a dt^p + b by least squares over ALL rungs (p scanned, a,b linear)."""
    E, dt = np.asarray(E, float), np.asarray(dt, float)
    if len(E) < 3:
        return None
    best = None
    for p in np.arange(0.20, 4.0005, 0.001):
        Xd = np.vstack([dt ** p, np.ones_like(dt)]).T
        coef, *_ = np.linalg.lstsq(Xd, E, rcond=None)
        r = E - Xd @ coef
        ss = float(np.sum(r ** 2))
        if best is None or ss < best[0]:
            best = (ss, p, coef[0], coef[1], r)
    ss, p, a, b, r = best
    scale = max(abs(E).max(), 1e-300)
    return dict(p=float(p), a=float(a), b=float(b),
                rms_rel=float(np.sqrt(ss / len(E)) / scale),
                maxres_rel=float(np.abs(r).max() / scale))


def fit_two_term(E, dt):
    """E = a dt + c dt^2 + b.  a -> 0 is the signature of second order."""
    E, dt = np.asarray(E, float), np.asarray(dt, float)
    if len(E) < 4:
        return None
    X = np.vstack([dt, dt ** 2, np.ones_like(dt)]).T
    coef, *_ = np.linalg.lstsq(X, E, rcond=None)
    a, c, b = map(float, coef)
    r = E - X @ coef
    scale = max(abs(E).max(), 1e-300)
    return dict(a=a, c=c, b=b,
                lin_frac=float(abs(a * dt.max()) / max(abs(E).max(), 1e-300)),
                rms_rel=float(np.sqrt(np.sum(r ** 2) / len(E)) / scale))


rows_by_set = {}
for d in sorted(glob.glob(os.path.join(SW, "output_tord_*"))):
    name = os.path.basename(d)[len("output_tord_"):]
    tag, nstr, mstr = name.rsplit("_", 2)
    n, M = int(nstr[1:]), int(mstr[1:])
    ic = os.path.join(SW, "soundwave_ideal_N%d.h5" % n)
    r = diag(d, ic, n)
    key = "%s_N%d" % (tag, n)
    if r is None:
        print("  (no data)     %s" % name); continue
    if isinstance(r, tuple):
        print("  (off-target %.9g) %s" % (r[1], name)); continue
    r["M"], r["dt"] = M, TARGET / M
    rows_by_set.setdefault(key, []).append(r)

out = {}
hdr = ("%-16s %5s   %-6s %-11s %-9s   %-9s %-9s   %-10s" %
       ("set", "rungs", "p", "b (floor)", "rms", "a (2-term)", "c", "nonisen@fin"))
print("\n" + hdr)
print("-" * len(hdr))
for key in sorted(rows_by_set):
    rows = sorted(rows_by_set[key], key=lambda r: r["M"])
    dt = [r["dt"] for r in rows]
    E = [r["amp_def"] for r in rows]
    f1, f2 = fit_pab(E, dt), fit_two_term(E, dt)
    fv = fit_pab([r["L1_vx"] for r in rows], dt)
    out[key] = dict(rows=rows, fit_damp=f1, fit_damp_2term=f2, fit_L1vx=fv)
    if f1:
        print("%-16s %5d   %-6.3f %-11.4e %-9.1e   %-9.3f %-9.2f   %-10.3e"
              % (key, len(rows), f1["p"], f1["b"], f1["rms_rel"],
                 (f2 or {}).get("a", float("nan")), (f2 or {}).get("c", float("nan")),
                 rows[-1]["nonisen_amp"]))

if "--json" in sys.argv:
    path = sys.argv[sys.argv.index("--json") + 1]
    json.dump(out, open(path, "w"), indent=1)
    print("\nwrote", path)
