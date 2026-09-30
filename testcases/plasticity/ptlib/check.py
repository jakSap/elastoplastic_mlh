"""Checkers of the plasticity test plan. Each returns (verdict, lines, data); `pt check` writes
checks/<test>.txt (human) and checks/<test>.json (verdict + numbers).

Verdicts: PASS, FAIL, INCOMPLETE (runs missing), LOOK (automatic criteria pass, a rendered movie
must still be inspected -> `pt accept`), with FINDING lines for reported-but-not-gating results
(TC footprint, P-R4)."""
import datetime
import glob
import json
import math
import os

import h5py
import numpy as np

from .core import run_state
from .gen import drive_matrix


# ------------------------------------------------------------------------ helpers
def snap_files(S, rel):
    return sorted(glob.glob(S.p("runs", rel, "output", "*.h5")))


def read(f, keys=None):
    with h5py.File(f, "r") as h:
        ks = keys or list(h.keys())
        out = {k: h[k][()] for k in ks if k in h}
    if "time" in out:
        out["t"] = float(np.ravel(out["time"])[0])
    return out


def pressure(d, mat):
    return d["P"] if "P" in d else mat["A"] * (d["rho"] / mat["rho0"] - 1.0)


def mat_of(S, test):
    return {**S.M.MATERIAL, **S.M.TESTS[test].get("material", {})}


def model_of(S, variant):
    f = S.M.VARIANTS[variant].get("flags", {})
    if S.code(variant) == "gizmo":
        return 1, 1                                     # GIZMO: true 2D, no-sqrt return
    return int(f.get("SOLID_2D_MODEL", 0)), int(f.get("YIELD_RETURN", 0))


def runs_by(S, test):
    """{(variant, rung): rel} of finished (exit 0) runs, plus the list of missing ones."""
    ok, missing = {}, []
    for v, r, rel in S.runs(test):
        st = run_state(S, rel)
        if st == "done" and snap_files(S, rel):
            ok[(v, r)] = rel
        else:
            missing.append("%s (%s)" % (rel, st))
    return ok, missing


def finest(rungs, have):
    rs = [r for r in rungs if r in have]
    return rs[-1] if rs else None


def fmt(x):
    return "%.3e" % x if isinstance(x, float) else str(x)


# ---------------------------------------------------------------- P1: bit identity
def check_identical(S, test, a, b):
    runs, missing = runs_by(S, test)
    lines, data, ok = [], {}, True
    T = S.M.TESTS[test]
    for r in T["rungs"]:
        if (a, r) not in runs or (b, r) not in runs:
            return "INCOMPLETE", lines + ["missing: " + ", ".join(missing)], data
        fa, fb = snap_files(S, runs[(a, r)]), snap_files(S, runs[(b, r)])
        if len(fa) != len(fb):
            ok = False
            lines.append("rung %s: %d vs %d snapshots" % (r, len(fa), len(fb)))
            continue
        first, nfield, gp_max = None, 0, 0.0
        for i, (x, y) in enumerate(zip(fa, fb)):
            da, db = read(x), read(y)
            for k in sorted(set(da) & set(db)):
                if k in ("t",):
                    continue
                nfield += 1
                if da[k].shape != db[k].shape or da[k].tobytes() != db[k].tobytes():
                    if first is None:
                        diff = np.max(np.abs(np.asarray(da[k], float) - np.asarray(db[k], float))) \
                            if da[k].shape == db[k].shape else float("nan")
                        first = "snapshot %d field %s max|diff| %.3e" % (i, k, diff)
            for d in (da, db):
                if "plasticStrain" in d:
                    gp_max = max(gp_max, float(np.max(d["plasticStrain"])))
        same = first is None and gp_max == 0.0
        ok &= same
        lines.append("rung %s: %d snapshots, %d field comparisons: %s%s" % (
            r, len(fa), nfield, "bitwise identical" if first is None else "DIFFERENT, first at " + first,
            "" if gp_max == 0 else "; max plasticStrain %.3e (must be 0)" % gp_max))
        data[r] = dict(identical=first is None, first_diff=first, gp_max=gp_max, snapshots=len(fa))
        # margin to yield (P1c): max sqrt(2 J2)/Y of the elastic run b
        if "Y0" in S.M.MATERIAL and test.startswith("P1c"):
            Y = mat_of(S, test)["Y0"]
            m = max(float(np.max(np.sqrt(2 * (0.5 * (d["Sxx"]**2 + d["Syy"]**2) + d["Sxy"]**2)))) / Y
                    for d in (read(f, ["Sxx", "Syy", "Sxy"]) for f in fb))
            lines.append("rung %s: elastic run max sqrt(2 J2)/Y = %.3e (yield at 1)" % (r, m))
            data[r]["max_J2_ratio"] = m
    return ("PASS" if ok else "FAIL"), lines, data


# ------------------------------------------------------------------ P1d: smoke look
def check_smoke(S, test):
    runs, missing = runs_by(S, test)
    if missing:
        return "INCOMPLETE", ["missing: " + ", ".join(missing)], {}
    lines, ok = [], True
    for (v, r), rel in runs.items():
        log = open(S.p("runs", rel, "run.log"), errors="replace").read()
        nerr = sum("ERROR" in l for l in log.split("\n"))
        f = snap_files(S, rel)
        d0, d1 = read(f[0], ["x"]), read(f[-1], ["x"])
        ext0 = np.ptp(d0["x"], axis=0).max()
        ext1 = np.ptp(d1["x"], axis=0).max()
        bad = nerr > 0 or ext1 > 3 * ext0
        ok &= not bad
        lines.append("%s: ERROR lines %d, extent %.3f -> %.3f%s" % (rel, nerr, ext0, ext1,
                                                                    "  <- FAIL" if bad else ""))
    lines.append("Render and inspect: ./pt movie %s ; then ./pt accept %s \"<what you saw>\"" % (test, test))
    return ("LOOK" if ok else "FAIL"), lines, {}


# ------------------------------------------------------------------ P1d: wave speed
def wave_phase(S, rel, mode):
    ts, ph, amp = [], [], []
    for f in snap_files(S, rel):
        d = read(f, ["time", "x", "v"])
        x = np.mod(d["x"][:, 0], 1.0)
        comp = d["v"][:, 0] if mode == "p" else d["v"][:, 1]
        a = np.mean(comp * np.exp(-2j * np.pi * x))
        ts.append(d["t"]); ph.append(np.angle(a)); amp.append(abs(a))
    ph = np.unwrap(ph)
    return np.array(ts), ph, np.array(amp)


def check_wavespeed(S, test, mode, c, tol, at):
    runs, missing = runs_by(S, test)
    T = S.M.TESTS[test]
    lines, data = [], {}
    for (v, r), rel in sorted(runs.items(), key=lambda kv: (kv[0][0], int(kv[0][1]))):
        t, ph, amp = wave_phase(S, rel, mode)
        slope = np.polyfit(t, ph, 1)[0]
        cn = -slope / (2 * np.pi)
        data.setdefault(v, {})[r] = dict(c_num=cn, err=cn / c - 1, damping=amp[-1] / amp[0])
        lines.append("%-5s n=%-4s c_num = %.7f  c_num/c - 1 = %+.3e  amplitude ratio %.5f" % (
            v, r, cn, cn / c - 1, amp[-1] / amp[0]))
    if missing:
        return "INCOMPLETE", lines + ["missing: " + ", ".join(missing)], data
    ok = True
    for v, byr in data.items():
        rs = [r for r in T["rungs"] if r in byr]
        errs = [abs(byr[r]["err"]) for r in rs]
        orders = [math.log(errs[i] / errs[i + 1]) / math.log(float(rs[i + 1]) / float(rs[i]))
                  for i in range(len(rs) - 1) if errs[i + 1] > 0]
        lines.append("%s: |c_num/c - 1| order(s) %s" % (v, ", ".join("%.2f" % o for o in orders)))
        if v != "gz":
            ok &= at in byr and abs(byr[at]["err"]) < tol and all(o >= 1.7 for o in orders)
    if "gz" in data:
        rs = [r for r in T["rungs"] if r in data["gz"] and any(r in data[v] for v in data if v != "gz")]
        if rs:
            r = rs[-1]
            v = [k for k in data if k != "gz" and r in data[k]][0]
            dc = abs(data[v][r]["c_num"] - data["gz"][r]["c_num"]) / c
            lines.append("newMFM vs GIZMO at n=%s: |dc|/c = %.3e (plan: < 1e-3)" % (r, dc))
            ok &= dc < 1e-3
    return ("PASS" if ok else "FAIL"), lines, data


# ------------------------------------------------------------ P0b: material point ODE
def dev_rate(mu, L, model):
    a, b, c, d = L[0, 0], L[0, 1], L[1, 0], L[1, 1]       # dvx_dx, dvx_dy, dvy_dx, dvy_dy
    if model == 1:
        rxx = mu * (a - d); ryy = -rxx
    else:
        rxx = 4 / 3 * mu * a - 2 / 3 * mu * d; ryy = -2 / 3 * mu * a + 4 / 3 * mu * d
    return rxx, ryy, mu * (c + b)


def material_point(mat, L0, t_end, model, nsteps=40000, times=None):
    """Continuum reference of a homogeneous ballistic deformation x = (I + L0 t) X:
    L(t) = L0 (I + L0 t)^-1, rho = rho0/det, S from the model with RK2 substeps and an exact
    radial return (dt -> 0), u from du/dt = sigma:D/rho, gamma_p from the returns."""
    mu, Y, A, r0 = mat["mu"], mat["Y0"], mat["A"], mat["rho0"]
    dt = t_end / nsteps
    s = np.zeros(3)                 # sxx, syy, sxy
    u, gp, t = mat["u0"], 0.0, 0.0
    out = []
    times = list(times) if times is not None else []

    def Lof(t):
        F = np.eye(2) + L0 * t
        return L0 @ np.linalg.inv(F), r0 / np.linalg.det(F)

    def rate(s, L):
        rxx, ryy, rxy = dev_rate(mu, L, model)
        R = L[1, 0] - L[0, 1]
        return np.array([rxx - s[2] * R, ryy + s[2] * R, rxy + 0.5 * (s[0] - s[1]) * R])

    k = 0
    for n in range(nsteps + 1):
        while k < len(times) and times[k] <= t + 0.5 * dt:
            L, rho = Lof(t)
            out.append((times[k], s.copy(), u, gp, rho))
            k += 1
        if n == nsteps:
            break
        L, rho = Lof(t + 0.5 * dt)
        s_mid = s + 0.5 * dt * rate(s, Lof(t)[0])
        s = s + dt * rate(s_mid, L)
        J2 = 0.5 * (s[0]**2 + s[1]**2) + s[2]**2
        if J2 > Y * Y / 2:
            f = math.sqrt(Y * Y / 2 / J2)
            gp += (1 - f) * math.sqrt(J2) / mu
            s *= f
        D = 0.5 * (L + L.T)
        P = A * (rho / r0 - 1)
        power = -P * np.trace(D) + s[0] * D[0, 0] + s[1] * D[1, 1] + 2 * s[2] * D[0, 1]
        u += dt * power / rho
        t += dt
    return out


def interior_series(S, rel, inner):
    fs = snap_files(S, rel)
    x0 = read(fs[0], ["x"])["x"]
    sel = (np.abs(x0[:, 0]) <= inner) & (np.abs(x0[:, 1]) <= inner)
    rows = []
    for f in fs:
        d = read(f, ["time", "time_S", "Sxx", "Syy", "Sxy", "u", "plasticStrain"])
        t = float(d["time_S"][0]) if "time_S" in d else d["t"]     # GIZMO: stress is half a step behind
        rows.append(dict(t=t, S=np.array([d[k][sel].mean() for k in ("Sxx", "Syy", "Sxy")]),
                         Sstd=np.array([d[k][sel].std() for k in ("Sxx", "Syy", "Sxy")]),
                         u=float(d["u"][sel].mean()),
                         gp=float(d["plasticStrain"][sel].mean()) if "plasticStrain" in d else None))
    return rows, int(sel.sum())


def check_p0b(S, test, drive, rate, inner, tol_elastic, tol_plateau, tol_plateau_prr, tol_scatter,
              tol_gz, tol_tc):
    runs, missing = runs_by(S, test)
    T = S.M.TESTS[test]
    mat = mat_of(S, test)
    Y = mat["Y0"]
    L0 = drive_matrix(drive, rate)
    lines, data, ok, findings = [], {}, True, []
    series = {}
    # GIZMO's rotation term has the opposite sign (2026-09-30: simple-shear normal stresses = -Jaumann
    # solution), so with rotation (ssh) GIZMO is judged on S_xy only; the sign is reported as a FINDING.
    comps = lambda v: [2] if (v == "gz" and drive == "ssh") else [0, 1, 2]
    for (v, r), rel in runs.items():
        rows, nin = interior_series(S, rel, inner)
        model, _ = model_of(S, v)
        ref = material_point(mat, L0, T["t_end"] * 1.02, model, times=[w["t"] for w in rows])
        c = comps(v)
        refJ = [math.sqrt(2 * (0.5 * (q[1][0]**2 + q[1][1]**2) + q[1][2]**2)) / Y for q in ref]
        tY = next((q[0] for q, j in zip(ref, refJ) if j >= 0.999), None)
        el = [np.max(np.abs(w["S"] - q[1])[c]) / Y for w, q in zip(rows, ref) if tY and w["t"] < 0.9 * tY]
        pl = [np.max(np.abs(w["S"] - q[1])[c]) / Y for w, q in zip(rows, ref) if tY and w["t"] >= 1.5 * tY]
        sc = max(np.max(w["Sstd"][c]) / Y for w in rows)
        if drive == "ssh":
            num = [(w["S"][0], q[1][0]) for w, q in zip(rows, ref) if tY and 0.3 * tY < w["t"] < 0.9 * tY]
            if num:
                findings.append("FINDING: %s d=%s simple-shear S_xx / Jaumann solution = %+.3f (+1 = correct "
                                "rotation sign)" % (v, r, float(np.mean([a / b for a, b in num]))))
        m = dict(N_interior=nin, t_yield_ref=tY, err_elastic=max(el) if el else None,
                 err_plateau=max(pl) if pl else None, scatter=sc)
        if rows[-1]["gp"] is not None:
            m["gp_rel_err"] = abs(rows[-1]["gp"] - ref[-1][3]) / max(ref[-1][3], 1e-30)
            m["heat_rel_err"] = abs((rows[-1]["u"] - mat["u0"]) - (ref[-1][2] - mat["u0"])) / \
                abs(ref[-1][2] - mat["u0"])
        data.setdefault(v, {})[r] = m
        series[(v, r)] = rows
    # table
    for v in T["variants"]:
        for r in T["rungs"]:
            m = data.get(v, {}).get(r)
            if m:
                lines.append("%-5s d=%-7s Nin=%-5d elastic %s  plateau %s  scatter %s  gp %s  heat %s   (/Y)" % (
                    v, r, m["N_interior"], fmt(m["err_elastic"]), fmt(m["err_plateau"]), fmt(m["scatter"]),
                    fmt(m.get("gp_rel_err", "-")), fmt(m.get("heat_rel_err", "-"))))
    if missing:
        return "INCOMPLETE", lines + findings + ["missing: " + ", ".join(missing)], data
    rf = T["rungs"][-1]

    def gate(v, key, tol):
        val = data[v][rf][key]
        good = val is not None and val <= tol
        lines.append("  %-5s %-12s %s <= %g : %s" % (v, key, fmt(val), tol, "ok" if good else "FAIL"))
        return good

    lines.append("gates at d = %s:" % rf)
    ok &= gate("pref", "err_elastic", tol_elastic) & gate("pref", "err_plateau", tol_plateau) \
        & gate("pref", "scatter", tol_scatter)
    for key in ("gp_rel_err", "heat_rel_err"):
        ok &= gate("pref", key, 0.05)
    if "prr" in data:
        ok &= gate("prr", "err_plateau", tol_plateau_prr)
    if "gz" in data:
        ok &= gate("gz", "err_elastic", tol_elastic) & gate("gz", "err_plateau", tol_plateau)

    def pair_diff(a, b):
        ra, rb = series[(a, rf)], series[(b, rf)]
        tb = np.array([w["t"] for w in rb])
        tY = data["pref"][rf]["t_yield_ref"] or 0
        diffs = []
        for w in ra:
            if w["t"] < 1.5 * tY or w["t"] > tb[-1]:
                continue
            sb = np.array([np.interp(w["t"], tb, [q["S"][c] for q in rb]) for c in range(3)])
            diffs.append(np.max(np.abs(w["S"] - sb)[comps(b)]) / Y)
        return max(diffs) if diffs else None

    if "gz" in data:
        dg = pair_diff("pref", "gz")
        good = dg is not None and dg <= tol_gz
        ok &= good
        lines.append("  pref vs gz plateau |dS|/Y %s <= %g : %s" % (fmt(dg), tol_gz, "ok" if good else "FAIL"))
    if "ptc0" in data:
        dt_ = pair_diff("pref", "ptc0")
        lines.append("  pref vs ptc0 plateau |dS|/Y %s (TC footprint, plan tol %g)" % (fmt(dt_), tol_tc))
        if dt_ is not None and dt_ > tol_tc:
            findings.append("FINDING: the tensile correction changes homogeneous plastic flow by %.2f%% Y" % (100 * dt_))
    lines += findings
    return ("PASS" if ok else "FAIL"), lines, data


# ------------------------------------------------------------------ P2: exact plate
def exact_plate(mat, up, model):
    """Hypo-elastic precursor (S_xx = -mu_eff ln eta until |S_xx| = Y_ax) + Rankine-Hugoniot for both
    waves, frame of the undisturbed target. model 1: mu_eff = mu, Y_ax = Y/sqrt2; model 0 (old
    hybrid): S_xx = -4/3 mu ln eta, S_yy = 2/3 mu ln eta, J2 over 2x2 = Y^2/2 -> |S_xx| = 2Y/sqrt5."""
    from scipy.optimize import brentq
    A, r0, mu, Y = mat["A"], mat["rho0"], mat["mu"], mat["Y0"]
    mu_eff, Y_ax = (mu, Y / math.sqrt(2)) if model == 1 else (4 / 3 * mu, 2 * Y / math.sqrt(5))
    Sel = lambda e: A * (e - 1) + mu_eff * math.log(e)
    Spl = lambda e: A * (e - 1) + Y_ax
    eH = math.exp(Y_ax / mu_eff)
    SH = Sel(eH)

    def jump(ra, rb, Sa, Sb):
        return math.sqrt((Sb - Sa) * (1 / ra - 1 / rb)), math.sqrt((Sb - Sa) / (1 / ra - 1 / rb)) / ra

    uH, Ue = jump(r0, r0 * eH, 0.0, SH)
    eb = brentq(lambda e: uH + jump(r0 * eH, r0 * e, SH, Spl(e))[0] - up, eH * (1 + 1e-9), 3.0)
    du, Up_rel = jump(r0 * eH, r0 * eb, SH, Spl(eb))
    return dict(U_e=Ue, sigma_hel=SH, u_hel=uH, eta_hel=eH, U_p=uH + Up_rel, sigma_b=Spl(eb), eta_b=eb,
                W_p=eH * Up_rel)            # Lagrangian speed of the plastic shock


def profile(X, q, x0, x1, w):
    edges = np.arange(x0, x1 + 0.5 * w, w)
    idx = np.digitize(X, edges) - 1
    ok = (idx >= 0) & (idx < len(edges) - 1)
    cnt = np.bincount(idx[ok], minlength=len(edges) - 1)
    s = np.bincount(idx[ok], weights=q[ok], minlength=len(edges) - 1)
    c = 0.5 * (edges[1:] + edges[:-1])
    m = cnt > 0
    return c[m], s[m] / cnt[m]


def last_cross(Xc, q, level):
    """Largest X where the profile is above level (front of a wave running to +X), interpolated."""
    above = np.where(q > level)[0]
    if len(above) == 0 or above[-1] == len(q) - 1:
        return None
    i = above[-1]
    return Xc[i] + (Xc[i + 1] - Xc[i]) * (q[i] - level) / (q[i] - q[i + 1])


def rise_width(Xc, q, lo, hi):
    a, b = last_cross(Xc, q, lo + 0.9 * (hi - lo)), last_cross(Xc, q, lo + 0.1 * (hi - lo))
    return None if a is None or b is None else b - a


def energy(d):
    return float(np.sum(d["m"] * (d["u"] + 0.5 * np.sum(d["v"]**2, axis=1))))


def plate_metrics(S, rel, mat, ex, up, t_fit, d, L):
    fs = snap_files(S, rel)
    d0 = read(fs[0])
    X = d0["x"][:, 0]
    tgt = X > 0
    E0 = energy(d0)
    te, Xe, Xp = [], [], []
    for f in fs:
        q = read(f)
        if q["t"] < t_fit[0] - 1e-9 or q["t"] > t_fit[1] + 1e-9:
            continue
        Xc, vp = profile(X[tgt], q["v"][tgt, 0] + up, 0, L, d)
        e1 = last_cross(Xc, vp, 0.5 * ex["u_hel"])
        p1 = last_cross(Xc, vp, 0.5 * (ex["u_hel"] + up))
        if e1 is not None and p1 is not None:
            te.append(q["t"]); Xe.append(e1); Xp.append(p1)
    q = read(fs[-1])
    Xc, vp = profile(X[tgt], q["v"][tgt, 0] + up, 0, L, d)
    _, sig = profile(X[tgt], (pressure(q, mat) - q["Sxx"])[tgt], 0, L, d)
    _, rho = profile(X[tgt], q["rho"][tgt], 0, L, d)
    xe, xp = Xe[-1], Xp[-1]
    ahead = Xc > xe + 5 * d
    sig0 = float(np.median(sig[ahead])) if ahead.any() else 0.0
    hel = (Xc > xp + 3 * d) & (Xc < xe - 3 * d)
    beh = (Xc > 0.25 * xp) & (Xc < 0.75 * xp)
    SH = float(np.median(sig[hel])) - sig0
    Sb = float(np.median(sig[beh])) - sig0
    uH = float(np.median(vp[hel]))
    etaH = float(np.median(rho[hel])) / mat["rho0"]
    We = np.polyfit(te, Xe, 1)[0]
    Wp = np.polyfit(te, Xp, 1)[0]
    m = dict(U_e=We, U_p=uH + Wp / etaH, sigma_hel=SH, sigma_b=Sb, u_hel=uH,
             width_e=(rise_width(Xc, vp, 0, ex["u_hel"]) or float("nan")) / d,
             width_p=(rise_width(Xc, vp, ex["u_hel"], up) or float("nan")) / d,
             dE=abs(energy(q) - E0) / abs(E0), negE=int(np.sum(q["u"] < 0)), sigma0=sig0)
    for k in ("U_e", "U_p", "sigma_hel", "sigma_b"):
        m["err_" + k] = abs(m[k] / ex[k] - 1)
    return m


def check_p2(S, test, up, t_fit, tol_U, tol_hel, tol_sigma, margin_gz):
    runs, missing = runs_by(S, test)
    T = S.M.TESTS[test]
    mat = mat_of(S, test)
    L = T["ic"][1]["L"]
    data, lines = {}, []
    exs = {m: exact_plate(mat, up, m) for m in (0, 1)}
    lines.append("exact (true 2D): U_e %.5f  sigma_HEL %.6f  u_HEL %.6f  U_p %.5f  sigma_b %.6f" % tuple(
        exs[1][k] for k in ("U_e", "sigma_hel", "u_hel", "U_p", "sigma_b")))
    for (v, r), rel in runs.items():
        data.setdefault(v, {})[r] = plate_metrics(S, rel, mat, exs[model_of(S, v)[0]], up, t_fit, float(r), L)
    keys = ("U_e", "U_p", "sigma_hel", "sigma_b")
    for v in T["variants"]:
        for r in T["rungs"]:
            m = data.get(v, {}).get(r)
            if m:
                lines.append("%-5s d=%-6s " % (v, r) + "  ".join("%s %.3e" % (k, m["err_" + k]) for k in keys)
                             + "  widths %.1f/%.1f d  dE %.1e  negE %d" % (m["width_e"], m["width_p"], m["dE"],
                                                                           m["negE"]))
    if missing:
        return "INCOMPLETE", lines + ["missing: " + ", ".join(missing)], data
    ok, rf = True, T["rungs"][-1]
    tol = dict(U_e=tol_U, U_p=tol_U, sigma_hel=tol_hel, sigma_b=tol_sigma)
    lines.append("gates (pref, d = %s):" % rf)
    for k in keys:
        e = [data["pref"][r]["err_" + k] for r in T["rungs"]]
        good = e[-1] <= tol[k] and all(e[i + 1] <= e[i] for i in range(len(e) - 1))
        ok &= good
        lines.append("  %-10s %s <= %g, falling %s : %s" % (k, fmt(e[-1]), tol[k],
                                                           " > ".join("%.2e" % x for x in e), "ok" if good else "FAIL"))
    good = data["pref"][rf]["dE"] <= 1e-10 and all(data["pref"][r]["negE"] == 0 for r in T["rungs"])
    ok &= good
    lines.append("  energy %.1e <= 1e-10, negE 0 : %s" % (data["pref"][rf]["dE"], "ok" if good else "FAIL"))
    if "gz" in data:
        for r in T["rungs"]:
            for k in keys:
                a, b = data["pref"][r]["err_" + k], data["gz"][r]["err_" + k]
                if a > b + margin_gz:
                    ok = False
                    lines.append("  vs GIZMO d=%s %s: newMFM %.3e > GIZMO %.3e + %g : FAIL" % (r, k, a, b, margin_gz))
        lines.append("  vs GIZMO: newMFM error <= GIZMO error + %g at every rung for every metric: %s" % (
            margin_gz, "checked"))
    return ("PASS" if ok else "FAIL"), lines, data


# ------------------------------------------------------- P3: 1D Lagrangian reference
def ref1d(mat, Lf, Lt, vf, t_out, cells, cfl=0.3, c0=1.0, c1=0.1):
    """Uniaxial strain, true-2D deviator (dS_xx/dt = mu dv/dx, J2 = S_xx^2), exact radial return,
    P = A(rho/rho0 - 1), von Neumann-Richtmyer + linear viscosity. Returns per output time the zone
    centres X (Lagrangian) and Sigma = P - S_xx (compression positive), v."""
    A, r0, mu, Y = mat["A"], mat["rho0"], mat["mu"], mat["Y0"]
    Ltot = Lf + Lt
    nz = int(round(cells * Ltot / Lt))
    Xn = np.linspace(-Lf, Lt, nz + 1)
    x = Xn.copy()
    dX = Xn[1:] - Xn[:-1]
    mz = r0 * dX
    mn = np.zeros(nz + 1); mn[:-1] += 0.5 * mz; mn[1:] += 0.5 * mz
    v = np.where(Xn < 0, vf, 0.0); v[np.isclose(Xn, 0.0)] = 0.5 * vf
    Sx = np.zeros(nz)
    rho = np.full(nz, r0)
    t, k, out = 0.0, 0, []
    Xc = 0.5 * (Xn[1:] + Xn[:-1])
    Ya = Y / math.sqrt(2)
    while k < len(t_out):
        dx = x[1:] - x[:-1]
        rho = mz / dx
        dv = v[1:] - v[:-1]
        c = np.sqrt((A * rho / r0 + mu) / rho)
        q = np.where(dv < 0, rho * (c0**2 * dv**2 + c1 * c * np.abs(dv)), 0.0)
        P = A * (rho / r0 - 1)
        sig = -P + Sx - q                                   # tension positive
        if t + 1e-14 >= t_out[k]:
            out.append((t, Xc.copy(), P - Sx, 0.5 * (v[1:] + v[:-1])))
            k += 1
            continue
        dt = cfl * np.min(dx / (c + np.abs(dv) + 1e-30))
        dt = min(dt, t_out[k] - t) if t_out[k] - t > 1e-14 else dt
        f = np.zeros(nz + 1)
        f[:-1] -= sig; f[1:] += sig                           # free ends: outside stress 0
        v = v + dt * f / mn
        xn = x + dt * v
        dxn = xn[1:] - xn[:-1]
        eps = (v[1:] - v[:-1]) / (0.5 * (dx + dxn))
        Sx = Sx + dt * mu * eps
        over = np.abs(Sx) > Ya
        Sx[over] = np.sign(Sx[over]) * Ya
        x = xn
        t += dt
    return out


def check_p3(S, test, vf, Lf, Lt, gauges, tol_peak, tol_arrival, margin_gz, ref_cells):
    runs, missing = runs_by(S, test)
    T = S.M.TESTS[test]
    mat = mat_of(S, test)
    lines, data = [], {}
    if not runs:
        return "INCOMPLETE", ["missing: " + ", ".join(missing)], data
    any_rel = next(iter(runs.values()))
    times = [read(f, ["time"])["t"] for f in snap_files(S, any_rel)]
    cache = S.p("checks", "P3_ref1d_%d.npz" % ref_cells)
    if os.path.exists(cache):
        z = np.load(cache, allow_pickle=True)
        ref, ref_half = list(z["ref"]), list(z["ref_half"])
    else:
        ref = ref1d(mat, Lf, Lt, vf, times, ref_cells)
        ref_half = ref1d(mat, Lf, Lt, vf, times, ref_cells // 2)
        np.savez(cache, ref=np.array(ref, dtype=object), ref_half=np.array(ref_half, dtype=object))
    sHEL = exact_plate(mat, 0.5 * vf, 1)["sigma_hel"]

    def curve(series, d):
        """peak Sigma vs its Lagrangian depth, and front arrival times at the gauges."""
        pk, xpk, te, xe = [], [], [], []
        for t, Xc, sig, v in series:
            m = Xc > 0
            Xb, sb = profile(Xc[m], sig[m], 0, Lt, d)
            i = int(np.argmax(sb))
            pk.append(sb[i]); xpk.append(Xb[i])
            e = last_cross(Xb, sb, 0.5 * sHEL)
            te.append(t); xe.append(e if e is not None else np.nan)
        pk, xpk, te, xe = map(np.array, (pk, xpk, te, xe))
        good = np.isfinite(xe)
        peak_at = {g: float(np.interp(g, xpk[np.argsort(xpk)], pk[np.argsort(xpk)])) for g in gauges}
        arr = {g: float(np.interp(g, xe[good], te[good])) for g in gauges if xe[good].max() > g}
        return peak_at, arr

    def run_series(rel):
        fs = snap_files(S, rel)
        X = read(fs[0], ["x"])["x"][:, 0]
        out = []
        for f in fs:
            q = read(f)
            out.append((q["t"], X, pressure(q, mat) - q["Sxx"], q["v"][:, 0]))
        return out

    for r in T["rungs"]:
        d = float(r)
        rp, ra = curve(ref, d)
        hp, _ = curve(ref_half, d)
        conv = max(abs(hp[g] / rp[g] - 1) for g in gauges)
        for (v, rr), rel in runs.items():
            if rr != r:
                continue
            p, a = curve(run_series(rel), d)
            m = dict(err_peak=max(abs(p[g] / rp[g] - 1) for g in gauges),
                     err_arrival=max(abs(a[g] / ra[g] - 1) for g in gauges if g in a and g in ra),
                     peaks={str(g): p[g] for g in gauges}, ref_peaks={str(g): rp[g] for g in gauges},
                     ref_selfconv=conv)
            data.setdefault(v, {})[r] = m
            lines.append("%-5s d=%-6s max peak err %.3e  max arrival err %.3e   (ref N vs N/2: %.1e)" % (
                v, r, m["err_peak"], m["err_arrival"], conv))
    if missing:
        return "INCOMPLETE", lines + ["missing: " + ", ".join(missing)], data
    ok, rf = True, T["rungs"][-1]
    for k, tol in (("err_peak", tol_peak), ("err_arrival", tol_arrival)):
        e = [data["pref"][r][k] for r in T["rungs"]]
        good = e[-1] <= tol and all(e[i + 1] <= e[i] for i in range(len(e) - 1))
        ok &= good
        lines.append("  pref %-11s %s <= %g, falling %s : %s" % (k, fmt(e[-1]), tol,
                                                                " > ".join("%.2e" % x for x in e), "ok" if good else "FAIL"))
    if data["pref"][rf]["ref_selfconv"] > 0.002:
        lines.append("  WARNING: reference not converged (N vs N/2 %.2e > 0.2%%)" % data["pref"][rf]["ref_selfconv"])
    if "gz" in data:
        for r in T["rungs"]:
            for k in ("err_peak", "err_arrival"):
                a, b = data["pref"][r][k], data["gz"][r][k]
                if a > b + margin_gz:
                    ok = False
                    lines.append("  vs GIZMO d=%s %s: %.3e > %.3e + %g : FAIL" % (r, k, a, b, margin_gz))
    return ("PASS" if ok else "FAIL"), lines, data


# ------------------------------------------------------------------ P4: Taylor bars
def components(x, link):
    from scipy.sparse.csgraph import connected_components
    from scipy.spatial import cKDTree
    pairs = cKDTree(x).query_pairs(link, output_type="ndarray")
    import scipy.sparse as sp
    n = len(x)
    g = sp.coo_matrix((np.ones(len(pairs)), (pairs[:, 0], pairs[:, 1])), shape=(n, n))
    return connected_components(g, directed=False)[0]


def taylor_metrics(S, rel, d, L, W):
    fs = snap_files(S, rel)
    d0 = read(fs[0], ["x"])
    X0 = d0["x"]
    right = X0[:, 0] > 0
    gtot, ts = [], []
    for f in fs:
        q = read(f, ["time", "plasticStrain"])
        if "plasticStrain" in q:
            gtot.append(float(q["plasticStrain"].sum())); ts.append(q["t"])
    q = read(fs[-1])
    m = {}
    Ls, Ws, frag = [], [], []
    for side in (right, ~right):
        x = q["x"][side]
        Lf = np.ptp(x[:, 0]) + d
        near = np.abs(x[:, 0]) <= np.abs(x[:, 0]).min() + 0.1 * Lf
        Ls.append(Lf); Ws.append(np.ptp(x[near, 1]) + d)
        frag.append(components(x, 1.6 * d))
    m.update(L_f=float(np.mean(Ls)), W_f=float(np.mean(Ws)), fragments=frag,
             mirror=float(abs(q["x"][right, 0].mean() + q["x"][~right, 0].mean())),
             rebound=float(abs(q["v"][right, 0].mean())), negE=int(np.sum(q["u"] < 0)))
    if "plasticStrain" in q:
        gp = q["plasticStrain"]
        m["gp_max"] = float(gp.max())
        pz = right & (gp > 1e-3)
        m["plastic_zone"] = float(X0[pz, 0].max() - X0[pz, 0].min() + d) if pz.any() else 0.0
        if len(ts) > 2:
            rate = np.diff(gtot) / np.diff(ts)
            i = int(np.argmax(rate))
            after = np.where(rate[i:] < 0.01 * rate[i])[0]
            m["t_arrest"] = float(ts[i + after[0] + 1]) if len(after) else None
    return m


def check_p4(S, test, L, W, v, tol_Lconv, tol_Wconv, tol_Lgz, tol_Wgz):
    runs, missing = runs_by(S, test)
    T = S.M.TESTS[test]
    lines, data = [], {}
    lines.append("sanity (1D Wilkins-Guinan, no criterion): L_f/L0 ~ exp(-rho v^2/2Y) = %.3f" % math.exp(
        -mat_of(S, test)["rho0"] * v * v / (2 * mat_of(S, test)["Y0"])))
    for (va, r), rel in runs.items():
        data.setdefault(va, {})[r] = taylor_metrics(S, rel, float(r), L, W)
    for va in T["variants"]:
        for r in T["rungs"]:
            m = data.get(va, {}).get(r)
            if m:
                lines.append("%-5s d=%-7s L_f/L0 %.4f  W_f/W0 %.4f  frags %s  mirror %.1e  negE %d  gp_max %s  "
                             "zone %s  t_arrest %s" % (va, r, m["L_f"] / L, m["W_f"] / W, m["fragments"], m["mirror"],
                                                       m["negE"], fmt(m.get("gp_max", "-")),
                                                       fmt(m.get("plastic_zone", "-")), fmt(m.get("t_arrest", "-"))))
    if missing:
        return "INCOMPLETE", lines + ["missing: " + ", ".join(missing)], data
    ok = True
    rs = T["rungs"]
    p = data["pref"]
    dL = [abs(p[rs[i + 1]]["L_f"] - p[rs[i]]["L_f"]) / L for i in range(len(rs) - 1)]
    dW = [abs(p[rs[i + 1]]["W_f"] - p[rs[i]]["W_f"]) / W for i in range(len(rs) - 1)]
    for name, dd, tol in (("L_f", dL, tol_Lconv), ("W_f", dW, tol_Wconv)):
        good = dd[-1] < tol and all(dd[i + 1] < dd[i] for i in range(len(dd) - 1))
        ok &= good
        lines.append("  pref %s rung differences %s (last < %g, shrinking): %s" % (
            name, " > ".join("%.2e" % x for x in dd), tol, "ok" if good else "FAIL"))
    good = all(m["fragments"] == [1, 1] and m["negE"] == 0 and m["mirror"] <= 1e-10 for m in p.values())
    ok &= good
    lines.append("  pref one body per bar, negE 0, mirror <= 1e-10 at every rung: %s" % ("ok" if good else "FAIL"))
    if "gz" in data:
        rf = rs[-1]
        a, b = abs(p[rf]["L_f"] - data["gz"][rf]["L_f"]) / L, abs(p[rf]["W_f"] - data["gz"][rf]["W_f"]) / W
        good = a <= tol_Lgz and b <= tol_Wgz
        ok &= good
        lines.append("  vs GIZMO at d=%s: |dL|/L0 %.3e <= %g, |dW|/W0 %.3e <= %g : %s" % (
            rf, a, tol_Lgz, b, tol_Wgz, "ok" if good else "FAIL"))
    return ("PASS" if ok else "FAIL"), lines, data


# -------------------------------------------------------------------------- driver
CHECKERS = dict(identical=check_identical, smoke=check_smoke, wavespeed=check_wavespeed, p0b=check_p0b,
                p2=check_p2, p3=check_p3, p4=check_p4)


def run_check(S, test):
    name, kw = S.M.TESTS[test]["check"]
    os.makedirs(S.p("checks"), exist_ok=True)
    verdict, lines, data = CHECKERS[name](S, test, **kw)
    when = datetime.datetime.now().strftime("%m-%d %H:%M")
    txt = "%s  %s  (%s, commit %s)\n" % (test, verdict, when, S.meta["commit"][:7]) + "\n".join(lines) + "\n"
    open(S.p("checks", test + ".txt"), "w").write(txt)
    json.dump(dict(verdict=verdict, when=when, commit=S.meta["commit"][:7], data=data),
              open(S.p("checks", test + ".json"), "w"), indent=1, default=float)
    return verdict, txt


def cmd_check(args):
    from .core import find_study
    S = find_study(args.study)
    for test in S.tests(args.tests):
        verdict, txt = run_check(S, test)
        print(txt)


def cmd_accept(args):
    """Record a human/model inspection for a LOOK verdict (elastic plan S.1)."""
    from .core import find_study, die
    S = find_study(args.study)
    j = S.p("checks", args.test + ".json")
    if not os.path.exists(j):
        die("run ./pt check %s first" % args.test)
    v = json.load(open(j))
    if v["verdict"] != "LOOK":
        die("%s is %s, not LOOK -- nothing to accept" % (args.test, v["verdict"]))
    v.update(verdict="PASS", accepted=args.note, when=datetime.datetime.now().strftime("%m-%d %H:%M"))
    json.dump(v, open(j, "w"), indent=1)
    open(S.p("checks", args.test + ".txt"), "a").write("ACCEPTED after inspection: %s\n" % args.note)
    print("%s: PASS (accepted: %s)" % (args.test, args.note))
