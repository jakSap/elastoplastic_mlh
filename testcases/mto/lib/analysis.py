"""Analysis of mto runs: per-run metrics + verdict (result.json), scoreboard, frames.

Snapshots are always selected by their stored /time, never by file index
(MeshlessScheme dump-index off-by-one, see wave_validation/state.md).

Per-run verdicts
  T1 rings   FAIL on: crash, rip (> 2 fragments >= 2 % N, or the two largest
             fragments < 99.5 % N), sticking (closest cross-body distance at the end
             < 2 d0: the bodies did not separate again).
             Clumps (nearest neighbour < 0.5 d0) are reported, not failed.
  T2 blocks  FAIL on: crash, rip, clumps, interior lines in rho or u.
  P1 rods    FAIL on: crash, rip, rods stuck together, contact-plane clumps above
             the GIZMO control's own count (+10 %), buried surface particles (fold),
             outline differing from GIZMO's by more than P1_OUTLINE_TOL * W0.
  T1/T2/P1 also FAIL on symmetry breaking: the ICs are mirror- (T1, P1) or point-
  (T2) symmetric, so a 99th-percentile partner mismatch > ASYM_TOL d0 (max over time)
  is a numerical defect (the max is reported as asym_max).
  Clumps are always counted WITHIN a body: two surface rows meeting at a contact
  sit closer than 0.5 d0 in GIZMO too.
  A1, W*     INFO per rung (FAIL only on crash/rip); the verdict is the ORDER over
             the ladder, computed in scoreboard().

The interior-line detector (T2): residual of a particle's value against the
least-squares plane through its neighbours (r < 3 d0), on interior particles
only (full lattice neighbourhood) AND inside the central circle of each block
(radius = CORE_FRAC x edge length, around the block's current centre of mass;
user 2026-10-02: surface-row bands are a pass). Particles with |residual| > LINE_TOL are
clustered (link 1.5 d0); a cluster is a *line* when it has >= 4 particles, an
extent >= 4 d0 and principal-axis elongation >= 3. Thresholds are a first guess,
to be calibrated against the frames of ref (see the testplan).
"""

import glob
import json
import math
import os
import re

import h5py
import numpy as np
from scipy.spatial import cKDTree
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

LINE_TOL = {"rho": 0.01, "u": 0.01}     # absolute (rho0 = 1, u0 = 1)
CLUMP = 0.5                              # nearest neighbour < CLUMP * d0
FRAG_MIN = 0.02                          # fragments >= 2 % of N count
LINK = 1.5                               # connectivity link length / d0
ORDER_PASS = 1.8
ORDER_FALLBACK = 0.8                      # A1 fallback: everything else passes, AM converges at 1st order (user 10-03)
CORE_FRAC = 0.25                         # T2 line region: circle r = CORE_FRAC * edge (see interior_lines)
# Symmetry verdict on the 99th percentile of the mirror-partner mismatch (max over time), not on
# the maximum: one particle stranded in a void dominates the max. Calibrated on the user's eye
# (P1 2x2 at w24, 2026-10-02): p99 = 1.42 / 2.34 / 0.83 judged asymmetric, 0.14 judged symmetric.
ASYM_TOL = 0.25                          # p99 mirror-partner mismatch / d0 (GIZMO: <= 1e-5)
P1_OUTLINE_TOL = 0.05                     # mean |outline - GIZMO outline| / W0


# ------------------------------------------------------------------ io

def snapshots(outdir):
    out = []
    for f in glob.glob(os.path.join(outdir, "*.h5")):
        try:
            with h5py.File(f, "r") as d:
                out.append((float(np.ravel(d["time"])[0]), f))
        except (OSError, KeyError):
            continue
    return sorted(out)


def load(path, keys=("x", "v", "m", "rho", "u")):
    with h5py.File(path, "r") as d:
        out = {k: d[k][:] for k in keys if k in d}
        out["t"] = float(np.ravel(d["time"])[0])
        out["E"] = float(np.ravel(d["energy"])[0]) if "energy" in d else float("nan")
    return out


def log_stats(D):
    neg, inj, warn_sml = 0, 0.0, 0
    log = os.path.join(D, "run.log")
    if not os.path.exists(log):
        return dict(negE=0, limiter_E=0.0)
    with open(log, errors="replace") as f:
        for line in f:
            if "Negative internal energy" in line:
                neg += 1
            elif "ENERGY_HALF_LIMIT" in line:
                m = re.search(r"injected energy = ([0-9.eE+-]+)", line)
                if m:
                    inj = float(m.group(1))
            elif "bad state" in line:
                warn_sml += 1
    return dict(negE=neg, limiter_E=inj, sml_bad_lines=warn_sml)


def finished(D, meta):
    t_last = meta.get("t_last") or 0.0
    ok = meta.get("exit") == 0 and t_last >= meta["t_end"] * (1 - 1e-6) - 1e-9
    return ok, t_last


# ------------------------------------------------------------------ body metrics

def intra_nn(x, nbody):
    """Nearest-neighbour distance within each body (bodies = equal index blocks, as gen.py
    writes them). Cross-body pairs at a contact are not clumping: two surface rows meeting
    at impact sit closer than 0.5 d0 in GIZMO too (P1 control)."""
    nn = np.empty(len(x))
    for b in np.array_split(np.arange(len(x)), nbody):
        dist, _ = cKDTree(x[b]).query(x[b], k=2)
        nn[b] = dist[:, 1]
    return nn


# Symmetries of the testcase IC (exact on the lattice): a physically symmetric problem must
# stay symmetric, so the partner mismatch is a pure numerical defect (user, 2026-10-02:
# ref breaks the rod symmetry, GIZMO does not).
SYMMETRIES = {"T1": [(-1, 1), (1, -1)], "P1": [(-1, 1), (1, -1)], "T2": [(-1, -1)]}


def sym_partners(x0, test):
    """Index maps i -> partner(i) for each symmetry (sx, sy): x -> (sx x, sy y)."""
    tree = cKDTree(x0)
    maps = []
    for sx, sy in SYMMETRIES.get(test, []):
        d, j = tree.query(x0 * [sx, sy])
        if d.max() < 1e-9 * max(1.0, np.abs(x0).max()):
            maps.append(((sx, sy), j))
    return maps


def sym_error(x, maps, d0, q=99):
    """Mirror-partner mismatch |x_i - S x_partner(i)| / d0, worst symmetry: (q-th percentile, max)."""
    if not maps:
        return float("nan"), float("nan")
    errs = [np.linalg.norm(x - x[j] * [sx, sy], axis=1) / d0 for (sx, sy), j in maps]
    return (max(float(np.percentile(e, q)) for e in errs), max(float(e.max()) for e in errs))


def body_gap(x, d0):
    """Closest distance between the two bodies (index halves), in d0. A welded pair sits at
    ~1 d0; two bodies that rebounded are many d0 apart."""
    h = len(x) // 2
    d, _ = cKDTree(x[:h]).query(x[h:], k=1)
    return float(d.min() / d0)


def geometry_metrics(x, d0, N, nbody):
    tree = cKDTree(x)
    nn = intra_nn(x, nbody)
    pairs = tree.query_pairs(LINK * d0, output_type="ndarray")
    g = coo_matrix((np.ones(len(pairs)), (pairs[:, 0], pairs[:, 1])), shape=(N, N))
    ncomp, lab = connected_components(g, directed=False)
    sizes = np.sort(np.bincount(lab))[::-1]
    nfrag = int(np.sum(sizes >= FRAG_MIN * N))
    return dict(dmin=float(nn.min() / d0), nclump=int(np.sum(nn < CLUMP * d0)),
                nfrag=nfrag, main_frac=float(sizes[:nbody].sum() / N), tree=tree, nn=nn)


def interior_lines(x, fields, d0, tree, core=None):
    """Line-shaped clusters of outliers in each field, on interior particles.

    core: optional mask; only flagged particles inside it count (T2: the central circle of
    each block -- user, 2026-10-02: bands in the surface rows are a pass, only strips in the
    middle count). The user's circle r = edge/2 touches the face midpoints, where the
    third-row edge of the surface bands sits, and would fail ref, which the user judged
    a strong pass; CORE_FRAC = 1/4 is used instead."""
    N = len(x)
    pairs = tree.query_pairs(3.0 * d0, output_type="ndarray")
    i = np.concatenate([pairs[:, 0], pairs[:, 1]])
    j = np.concatenate([pairs[:, 1], pairs[:, 0]])
    cnt = np.bincount(i, minlength=N)
    full = np.percentile(cnt, 90)
    interior = cnt >= 0.95 * full
    dx = x[j] - x[i]
    M = np.zeros((N, 3, 3))
    basis = np.column_stack([np.ones(len(i)), dx[:, 0], dx[:, 1]])
    for a in range(3):
        for b in range(3):
            M[:, a, b] = np.bincount(i, weights=basis[:, a] * basis[:, b], minlength=N)
    M[~interior] = np.eye(3)
    out = {}
    for name, f in fields.items():
        df = f[j] - f[i]
        rhs = np.column_stack([np.bincount(i, weights=basis[:, a] * df, minlength=N)
                               for a in range(3)])
        rhs[~interior] = 0.0
        coef = np.linalg.solve(M, rhs[:, :, None])[:, :, 0]
        resid = -coef[:, 0]                       # f_i minus the neighbours' plane at x_i
        flag = interior & (np.abs(resid) > LINE_TOL[name])
        if core is not None:
            flag &= core
        idx = np.flatnonzero(flag)
        longest = 0.0
        if len(idx) >= 4:
            sub = cKDTree(x[idx]).query_pairs(LINK * d0, output_type="ndarray")
            g = coo_matrix((np.ones(len(sub)), (sub[:, 0], sub[:, 1])),
                           shape=(len(idx), len(idx)))
            nc, lab = connected_components(g, directed=False)
            for c in range(nc):
                pts = x[idx[lab == c]]
                if len(pts) < 4:
                    continue
                ev, evec = np.linalg.eigh(np.cov(pts.T))
                proj = (pts - pts.mean(0)) @ evec[:, -1]
                extent = (proj.max() - proj.min()) / d0
                elong = math.sqrt(ev[-1] / max(ev[0], 1e-30))
                if extent >= 4.0 and elong >= 3.0:
                    longest = max(longest, extent)
        out[name] = dict(nflag=int(flag.sum()), line=longest,
                         max_resid=float(np.abs(resid[interior]).max()) if interior.any() else 0.0)
    return out


def analyse_body(D, meta, test):
    snaps = snapshots(os.path.join(D, "output"))
    with h5py.File(meta.get("ic") or os.path.join(os.path.dirname(D), "..", "..", "..", "ics",
                                "{}_{}.h5".format(test, meta["tag"])), "r") as f:
        attrs = dict(f.attrs)
        x0 = f["x"][:]
    d0 = float(attrs["spacing"])
    N = len(x0)
    half = N // 2
    maps = sym_partners(x0, test)
    series = []
    lines = {"rho": 0.0, "u": 0.0}
    stride = 1 if test == "T2" else max(1, len(snaps) // 60)
    for k, (t, path) in enumerate(snaps):
        s = load(path)
        x = s["x"]
        g = geometry_metrics(x, d0, N, 2)
        row = dict(t=t, dmin=g["dmin"], nclump=g["nclump"], nfrag=g["nfrag"],
                   main_frac=g["main_frac"], rho_min=float(s["rho"].min()),
                   rho_max=float(s["rho"].max()), u_min=float(s["u"].min()), E=s["E"],
                   com_dist=float(np.linalg.norm(x[half:].mean(0) - x[:half].mean(0))),
                   asym=sym_error(x, maps, d0)[0], asym_max=sym_error(x, maps, d0)[1])
        m, v = s["m"], s["v"]
        row["Lz"] = float(np.sum(m * (x[:, 0] * v[:, 1] - x[:, 1] * v[:, 0])))
        if test == "T2" and k % stride == 0:
            core = np.zeros(N, bool)
            for blk in (slice(0, half), slice(half, N)):
                core[blk] = (np.linalg.norm(x[blk] - x[blk].mean(0), axis=1)
                             < CORE_FRAC * float(attrs["size"]))
            li = interior_lines(x, {"rho": s["rho"], "u": s["u"]}, d0, g["tree"], core)
            row["line_rho"], row["line_u"] = li["rho"]["line"], li["u"]["line"]
            for f in lines:
                lines[f] = max(lines[f], li[f]["line"])
        series.append(row)
    if not series:
        return dict(verdict="FAIL", summary="no snapshots", series=[])
    S = {k: np.array([r[k] for r in series]) for k in series[0] if k != "t"}
    t = np.array([r["t"] for r in series])
    res = dict(d0=d0, N=N, t=t.tolist(),
               dmin=float(S["dmin"].min()), nclump_max=int(S["nclump"].max()),
               nfrag_max=int(S["nfrag"].max()), main_frac_min=float(S["main_frac"].min()),
               rho_min=float(S["rho_min"].min()), rho_max=float(S["rho_max"].max()),
               u_min=float(S["u_min"].min()),
               dE_rel=float(abs(S["E"][-1] / S["E"][0] - 1)) if S["E"][0] else float("nan"),
               Lz_rel_end=float(S["Lz"][-1] / S["Lz"][0] - 1) if abs(S["Lz"][0]) > 1e-14 else None,
               com_start=float(S["com_dist"][0]), com_min=float(S["com_dist"].min()),
               com_end=float(S["com_dist"][-1]),
               line_rho=lines["rho"], line_u=lines["u"], asym=float(np.nanmax(S["asym"])),
               asym_max=float(np.nanmax(S["asym_max"])))
    res["rip"] = res["nfrag_max"] > 2 or res["main_frac_min"] < 0.995
    res["gap_end"] = body_gap(load(snaps[-1][1], ("x",))["x"], d0)
    res["stuck"] = res["gap_end"] < 2.0
    res["series"] = {k: v.tolist() for k, v in S.items()}
    return res


# ------------------------------------------------------------------ disc

def analyse_disc(D, meta):
    snaps = snapshots(os.path.join(D, "output"))
    L = []
    for t, path in snaps:
        s = load(path, ("x", "v", "m"))
        x, v, m = s["x"], s["v"], s["m"]
        L.append((t, float(np.sum(m * (x[:, 0] * v[:, 1] - x[:, 1] * v[:, 0])))))
    t = np.array([a for a, _ in L])
    Lz = np.array([b for _, b in L])
    rel = Lz / Lz[0] - 1
    s_end = load(snaps[-1][1], ("x",))
    d0 = meta["d"]
    g = geometry_metrics(s_end["x"], d0, len(s_end["x"]), 1)
    return dict(t=t.tolist(), dLz=rel.tolist(), dLz_T=float(rel[-1]),
                dLz_max=float(np.abs(rel).max()), rip=g["nfrag"] > 1 or g["main_frac"] < 0.995,
                dmin=g["dmin"])


# ------------------------------------------------------------------ waves

def analyse_wave(D, meta, test):
    ic = meta.get("ic") or os.path.join(os.path.dirname(D), "..", "..", "..", "ics",
                                        "{}_{}.h5".format(test, meta["tag"]))
    with h5py.File(ic, "r") as f:
        A = float(f.attrs["amplitude"])
        k = float(np.linalg.norm(f.attrs["k"]))
        c = float(f.attrs["c_wave"])
        mode = str(f.attrs["mode"])
        rho_bar = float(f.attrs["rho_bar"])
    Av = c * A
    comp = 0 if mode == "p" else 1
    # Material viscosity (MATERIAL_VISCOSITY): the reference is the damped linear solution, not the undamped wave.
    # The pair force is the Morris form MINUS the mean spin W r; to second order in r the spin term removes
    # half of the Laplacian and adds half a grad div, so the operator is div(eta D) = eta/2 (lap v + grad div v)
    # (D = sym grad v), not eta lap v: P mode eta k^2, S mode eta k^2 / 2 (measured 10-03: S rates 1.148,
    # 1.052, 1.007 -> eta k^2/4rho = 0.987, P rates -> eta k^2/2rho = 1.974).
    #   d/dt [v, s] = [[-G, ik/rho], [ik rho c^2, 0]] [v, s],  G = f eta k^2 / rho (f = 1 P, 1/2 S),  s0 = -rho c v0
    # with the same initial condition; amplitude exp(-G t / 2), phase speed sqrt(c^2 - (G/2k)^2).
    eta = float((meta.get("materials") or [{}])[0].get("eta", 0.0))
    visc = eta > 0.
    G = (1.0 if mode == "p" else 0.5) * eta * k * k / rho_bar

    def exact_v(t):
        """(amplitude, phase) of the exact v(x,t) = amp sin(kx + phase), amplitude as a fraction of Av."""
        if not visc:
            return 1.0, -k * c * t
        Mat = np.array([[-G, 1j * k / rho_bar], [1j * k * rho_bar * c * c, 0.0]])
        w, V = np.linalg.eig(Mat)
        y = V @ (np.exp(w * t) * np.linalg.solve(V, np.array([1.0, -rho_bar * c], dtype=complex)))
        return abs(y[0]), float(np.angle(y[0]))
    rows = []
    for t, path in snapshots(os.path.join(D, "output")):
        s = load(path, ("x", "v", "rho"))
        x = s["x"][:, 0]
        X = np.column_stack([np.ones_like(x), np.sin(k * x), np.cos(k * x),
                             np.sin(2 * k * x), np.cos(2 * k * x)])
        out = dict(t=t)
        for name, f, amp in (("v", s["v"][:, comp], Av), ("rho", s["rho"], rho_bar * A)):
            if name == "rho" and (mode == "s" or visc):
                continue
            coef, *_ = np.linalg.lstsq(X, f, rcond=None)
            a = math.hypot(coef[1], coef[2])
            psi = math.atan2(coef[2], coef[1])          # f ~ a sin(kx + psi)
            a_ex, psi_ex = exact_v(t) if name == "v" else (1.0, -k * c * t)
            dpsi = (psi - psi_ex + math.pi) % (2 * math.pi) - math.pi
            mean = coef[0] if name == "v" else rho_bar
            exact = mean + amp * a_ex * np.sin(k * x + psi_ex)
            out[name] = dict(a=a / amp, a_ex=a_ex, dphi=dpsi, harm2=math.hypot(coef[3], coef[4]) / amp,
                             L1=float(np.mean(np.abs(f - exact)) / amp))
        rows.append(out)
    t = np.array([r["t"] for r in rows])
    dphi = np.unwrap([r["v"]["dphi"] for r in rows])
    slope = np.polyfit(t, dphi, 1)[0] if len(t) > 2 else float("nan")
    fin = rows[-1]
    # c_num/c - 1, positive = wave too fast: from the phase after one traversal (c_err,
    # the user's criterion) and from the phase drift over all snapshots (freq_err).
    res = dict(t_final=float(t[-1]),
               c_err=float(-fin["v"]["dphi"] / (k * c * t[-1])),
               freq_err=float(-slope / (k * c)),
               damp_v=float(fin["v"]["a_ex"] - fin["v"]["a"]), dphi_v=float(fin["v"]["dphi"]),
               eta=eta, a_final=float(fin["v"]["a"]), a_exact=float(fin["v"]["a_ex"]),
               gamma_meas=float(-2 * math.log(max(fin["v"]["a"], 1e-300)) / t[-1] / 2) if eta > 0 else 0.0,
               gamma_th=float(G / 2),
               L1_v=float(fin["v"]["L1"]), harm2_v=float(fin["v"]["harm2"]),
               series=rows)
    if mode == "p" and not visc:
        res.update(damp_rho=float(1 - fin["rho"]["a"]), L1_rho=float(fin["rho"]["L1"]))
    return res


# ------------------------------------------------------------------ rods

def rod_shape(x, d0, half, W0):
    """Metrics of the right rod (second half of the particles, x > 0)."""
    xb = x[half:]
    xmin = xb[:, 0].min()
    length = xb[:, 0].max() - xmin
    near = xb[:, 0] - xmin < 3 * d0
    width = 2 * np.abs(xb[near, 1]).max() if near.any() else float("nan")
    bins = np.arange(0.0, length + d0, d0)
    idx = np.clip(((xb[:, 0] - xmin) / d0).astype(int), 0, len(bins) - 1)
    outline = np.zeros(len(bins))
    np.maximum.at(outline, idx, np.abs(xb[:, 1]))
    return dict(length=length, width=width, outline=outline, idx=idx, xmin=xmin)


def plane_metrics(snaps, d0, half):
    """Clumping in the contact plane (|x| < 3 d0) over the run + body COM distances."""
    plane_clumps, plane_dmin, com = 0, 9.9, []
    for t, path in snaps[:: max(1, len(snaps) // 40)] + [snaps[-1]]:
        x = load(path, ("x",))["x"]
        com.append(float(x[half:, 0].mean() - x[:half, 0].mean()))
        nn = intra_nn(x, 2) / d0
        plane = np.abs(x[:, 0]) < 3 * d0
        if plane.any():
            plane_clumps = max(plane_clumps, int(np.sum(nn[plane] < CLUMP)))
            plane_dmin = min(plane_dmin, float(nn[plane].min()))
    return plane_clumps, plane_dmin, com


def analyse_rod(D, meta):
    snaps = snapshots(os.path.join(D, "output"))
    ic = os.path.join(os.path.dirname(D), "..", "..", "..", "ics", "P1_{}.h5".format(meta["tag"]))
    with h5py.File(ic, "r") as f:
        attrs = dict(f.attrs)
        x0 = f["x"][:]
    d0, W0, L0 = float(attrs["spacing"]), float(attrs["width"]), float(attrs["length"])
    N = len(x0)
    half = N // 2
    surface0 = np.abs(x0[half:, 1]) >= W0 / 2 - 0.75 * d0
    plane_clumps, plane_dmin, com = plane_metrics(snaps, d0, half)
    maps = sym_partners(x0, "P1")
    se = [sym_error(load(p, ("x",))["x"], maps, d0) for _, p in snaps]
    asym, asym_max = max(a for a, _ in se), max(b for _, b in se)
    s = load(snaps[-1][1], ("x", "v", "m"))
    x = s["x"]
    shp = rod_shape(x, d0, half, W0)
    buried = int(np.sum(surface0 & (np.abs(x[half:, 1]) < shp["outline"][shp["idx"]] - 2.5 * d0)))
    g = geometry_metrics(x, d0, N, 2)
    res = dict(d0=d0, L_ratio=float(shp["length"] / L0), W_ratio=float(shp["width"] / W0),
               plane_clumps=plane_clumps, plane_dmin=plane_dmin, buried=buried,
               outline=(shp["outline"] / W0).tolist(), rip=g["nfrag"] > 2 or g["main_frac"] < 0.995,
               KE_end=float(0.5 * np.sum(s["m"] * np.sum(s["v"] ** 2, axis=1))),
               gap_end=body_gap(x, d0), stuck=bool(body_gap(x, d0) < 2.0), asym=asym,
               asym_max=asym_max)
    res["final_clumps"] = int(np.sum(intra_nn(x, 2) / d0 < CLUMP))
    eta = float((meta.get("materials") or [{}])[0].get("eta", 0.0))
    ctrl_name = meta["tag"] + ("_eta{:g}".format(eta) if eta > 0 else "")
    res["gizmo_control"] = ctrl_name
    # original GIZMO control, and the one with the objective Jaumann term (10-03: GIZMO's
    # S W - W S is not objective; the objective control widens the mushroom, W/W0 1.60 -> 1.74 at w16)
    for key, name in (("gizmo", ctrl_name), ("gizmo_obj", ctrl_name + "_objjau")):
        ctrl = os.path.join(os.path.dirname(__file__), "..", "gizmo", "runs", name, "conv")
        if not (os.path.isdir(ctrl) and snapshots(ctrl)):
            continue
        gs = snapshots(ctrl)
        tg, pg = min(gs, key=lambda a: abs(a[0] - snaps[-1][0]))
        xg = load(pg, ("x",))["x"]
        gsh = rod_shape(xg, d0, len(xg) // 2, W0)
        n = min(len(gsh["outline"]), len(shp["outline"]))
        gpc, gpd, gcom = plane_metrics(gs, d0, len(xg) // 2)
        gasym = max(sym_error(load(p, ("x",))["x"], maps, d0)[0] for _, p in gs)
        res[key] = dict(t=tg, L_ratio=float(gsh["length"] / L0), plane_clumps=gpc,
                        plane_dmin=gpd, asym=gasym, gap_end=body_gap(xg, d0),
                        stuck=bool(body_gap(xg, d0) < 2.0),
                        W_ratio=float(gsh["width"] / W0),
                        final_clumps=int(np.sum(intra_nn(xg, 2) / d0 < CLUMP)),
                        outline=(gsh["outline"] / W0).tolist())
        sfx = "" if key == "gizmo" else "_obj"
        res["dL_vs_gizmo" + sfx] = float(abs(shp["length"] - gsh["length"]) / L0)
        res["dW_vs_gizmo" + sfx] = float(abs(shp["width"] - gsh["width"]) / W0)
        res["outline_L1_vs_gizmo" + sfx] = float(np.mean(np.abs(shp["outline"][:n]
                                                                - gsh["outline"][:n])) / W0)
    return res


def p1_control_fails(res, key, sfx):
    """P1 checks relative to a GIZMO control: plane clumps at impact (GIZMO clumps there too),
    clumps left at the end (GIZMO: 0 at every rung, 10-03), outline."""
    g = res.get(key)
    out = []
    cmax = int(1.1 * g["plane_clumps"]) + 2 if g else 0
    if res.get("plane_clumps", 0) > cmax:
        out.append("plane_clumps={}>{}".format(res["plane_clumps"], cmax))
    if g and "final_clumps" in g and "final_clumps" in res:
        fmax = int(1.1 * g["final_clumps"]) + 4
        if res["final_clumps"] > fmax:
            out.append("final_clumps={}>{}".format(res["final_clumps"], fmax))
    if g and res.get("outline_L1_vs_gizmo" + sfx, 0) > P1_OUTLINE_TOL:
        out.append("outline_vs_gizmo{}={:.3f}".format(sfx, res["outline_L1_vs_gizmo" + sfx]))
    return out


# ------------------------------------------------------------------ dispatch

def analyse_run(D):
    meta = json.load(open(os.path.join(D, "meta.json")))
    test = meta["test"]
    ok, t_last = finished(D, meta)
    res = dict(test=test, tag=meta["tag"], cand=meta["cand"], exit=meta.get("exit"),
               t_last=t_last, wall_s=meta.get("wall_s"), N=meta.get("N"), **log_stats(D))
    try:
        if test in ("T1", "T2"):
            res.update(analyse_body(D, meta, test))
        elif test == "A1":
            res.update(analyse_disc(D, meta))
        elif test == "P1":
            res.update(analyse_rod(D, meta))
        else:
            res.update(analyse_wave(D, meta, test))
    except (IndexError, ValueError, OSError, KeyError) as e:
        res["analysis_error"] = repr(e)
    fails = []
    if not ok:
        fails.append("crash(exit={}, t={:.4g})".format(meta.get("exit"), t_last))
    if res.get("rip"):
        fails.append("rip")
    if test in ("T1", "T2", "P1") and res.get("asym", 0) > ASYM_TOL:
        fails.append("asym={:.2g}d0".format(res["asym"]))
    if test == "T1" and res.get("stuck"):
        fails.append("stuck")
    if test == "T2":
        if res.get("nclump_max", 0) > 0:
            fails.append("clumps={}".format(res["nclump_max"]))
        for f in ("rho", "u"):
            if res.get("line_" + f, 0) > 0:
                fails.append("line_{}={:.0f}d0".format(f, res["line_" + f]))
    if test == "P1":
        if res.get("stuck"):
            fails.append("stuck")
        if res.get("buried", 0) > 0:
            fails.append("fold(buried={})".format(res["buried"]))
        base = list(fails)
        fails += p1_control_fails(res, "gizmo", "")
        # verdict against the objective-Jaumann GIZMO control (where it exists)
        if res.get("gizmo_obj"):
            fo = base + p1_control_fails(res, "gizmo_obj", "_obj")
            if not ok:
                fo = ["crash"] + fo
            res["fails_obj"], res["verdict_obj"] = fo, ("FAIL" if fo else "PASS")
    if "analysis_error" in res:
        fails.append("analysis_error")
    if fails:
        res["verdict"] = "FAIL"
    elif test in ("T1", "T2", "P1"):
        res["verdict"] = "PASS"
    else:
        res["verdict"] = "INFO"
    res["fails"] = fails
    res["summary"] = summary_line(res)
    json.dump(res, open(os.path.join(D, "result.json"), "w"), indent=1, default=float)
    return res


def summary_line(r):
    t = r["test"]
    f = " ".join(r.get("fails", []))
    if t in ("T1", "T2"):
        s = "dmin={:.2f} clumps={} frag={} rho=[{:.2f},{:.2f}] negE={} limE={:.2g}".format(
            r.get("dmin", float("nan")), r.get("nclump_max"), r.get("nfrag_max"),
            r.get("rho_min", float("nan")), r.get("rho_max", float("nan")),
            r.get("negE"), r.get("limiter_E", 0))
        s += " asym={:.2g}d0".format(r.get("asym", float("nan")))
        if t == "T2":
            s += " lines rho/u={:.0f}/{:.0f}".format(r.get("line_rho", 0), r.get("line_u", 0))
    elif t == "A1":
        s = "dLz(T)={:+.3e} max={:.3e}".format(r.get("dLz_T", float("nan")),
                                                r.get("dLz_max", float("nan")))
    elif t == "P1":
        s = "gap={:.1f}d0 L/L0={:.3f} W/W0={:.3f} plane_clumps={} buried={} asym={:.2g}d0".format(
            r.get("gap_end", float("nan")), r.get("L_ratio", float("nan")),
            r.get("W_ratio", float("nan")), r.get("plane_clumps"), r.get("buried"),
            r.get("asym", float("nan")))
        if "gizmo" in r:
            s += " | GIZMO gap={:.1f}d0 L/L0={:.3f} W/W0={:.3f} clumps={} outline dL1={:.3f}".format(
                r["gizmo"]["gap_end"], r["gizmo"]["L_ratio"], r["gizmo"]["W_ratio"],
                r["gizmo"]["plane_clumps"], r["outline_L1_vs_gizmo"])
    else:
        s = "L1_v={:.3e} damp={:.3e} c_err={:+.3e} freq_err={:+.3e}".format(
            r.get("L1_v", float("nan")), r.get("damp_v", float("nan")),
            r.get("c_err", float("nan")), r.get("freq_err", float("nan")))
        if "L1_rho" in r:
            s += " L1_rho={:.3e}".format(r["L1_rho"])
    return (s + ("  [" + f + "]" if f else "")).strip()


# ------------------------------------------------------------------ scoreboard

def _load_results(RUNS, cand, test):
    out = {}
    for p in glob.glob(os.path.join(RUNS, cand, test, "*", "result.json")):
        r = json.load(open(p))
        out[r["tag"]] = r
    return out


def order_fit(hs, errs):
    hs, errs = np.asarray(hs, float), np.abs(np.asarray(errs, float))
    if len(hs) < 2 or np.any(errs <= 0) or np.any(~np.isfinite(errs)):
        return float("nan"), []
    pair = [math.log(errs[i] / errs[i + 1]) / math.log(hs[i] / hs[i + 1])
            for i in range(len(hs) - 1)]
    return float(np.polyfit(np.log(hs), np.log(errs), 1)[0]), pair


def verdict_test(RUNS, cand, test):
    if test == "T1r":
        # relaxed ring criterion (user priorities 10-03): one run at d0.05 that does not fall apart
        r = _load_results(RUNS, cand, "T1").get("d0.05")
        if not r:
            return "-", ""
        bad = (r.get("exit") != 0 or r.get("rip") or r.get("nfrag_max", 2) > 2
               or r.get("main_frac_min", 1.0) < 0.995)
        return ("FAIL" if bad else "PASS"), "d0.05 frag={} main_frac={:.4f}".format(
            r.get("nfrag_max"), r.get("main_frac_min", float("nan")))
    if test == "P1obj":
        # P1 against the objective-Jaumann GIZMO control
        R = {g: dict(r, verdict=r["verdict_obj"], fails=r.get("fails_obj", []))
             for g, r in _load_results(RUNS, cand, "P1").items() if "verdict_obj" in r}
        test = "P1"
    else:
        R = _load_results(RUNS, cand, test)
    if not R:
        return "-", ""
    if test in ("T1", "T2", "P1"):
        gate = {"T1": ["d0.1", "d0.08", "d0.06"],
                "T2": ["d0.05", "d0.025", "d0.0177", "d0.0125"],
                "P1": ["w10", "w16", "w24"]}[test]
        have = [g for g in gate if g in R]
        bad = [g for g in have if R[g]["verdict"] == "FAIL"]
        if bad:
            return "FAIL", "fails at " + ",".join(
                "{}({})".format(g, " ".join(R[g]["fails"])) for g in bad)
        if len(have) < len(gate):
            return "partial", "{}/{} rungs pass".format(len(have), len(gate))
        return "PASS", "all {} rungs".format(len(gate))
    if any(r["verdict"] == "FAIL" for r in R.values()):
        return "FAIL", "crash/rip at " + ",".join(g for g, r in R.items() if r["verdict"] == "FAIL")
    if test == "A1":
        tags = sorted(R, key=lambda g: -float(g[1:]))
        hs = [R[g].get("d0") or float(g[1:]) for g in tags]
        e = [R[g]["dLz_T"] for g in tags]
        p, pair = order_fit(hs, e)
        ok = len(tags) >= 3 and p >= ORDER_PASS
        return ("PASS" if ok else "FAIL" if len(tags) >= 3 else "partial",
                "p={:.2f} pairs={} |dLz|={}".format(
                    p, ",".join("{:.2f}".format(q) for q in pair),
                    ",".join("{:.1e}".format(abs(x)) for x in e)))
    tags = sorted([g for g in R if int(g[1:]) >= 64], key=lambda g: int(g[1:]))
    if len(tags) < 2:
        return "partial", "need n64/n128/n256"
    hs = [1.0 / int(g[1:]) for g in tags]
    metrics = ["L1_v", "damp_v", "c_err"] + (["L1_rho"] if "L1_rho" in R[tags[0]] else [])
    parts, ok = [], len(tags) >= 3
    for mname in metrics:
        p, pair = order_fit(hs, [R[g][mname] for g in tags])
        last = pair[-1] if pair else float("nan")
        ok = ok and last >= ORDER_PASS
        parts.append("{} p={:.2f}".format(mname, last))
    return ("PASS" if ok else "FAIL" if len(tags) >= 3 else "partial"), " ".join(parts)


def a1_order(RUNS, cand):
    """Fitted A1 order over the full ladder (>= 3 rungs, no crash), else nan."""
    R = _load_results(RUNS, cand, "A1")
    if len(R) < 3 or any(r["verdict"] == "FAIL" for r in R.values()):
        return float("nan")
    tags = sorted(R, key=lambda g: -float(g[1:]))
    return order_fit([R[g].get("d0") or float(g[1:]) for g in tags], [R[g]["dLz_T"] for g in tags])[0]


def scoreboard(RUNS, cands):
    tests = ["T1", "T1r", "T2", "P1", "P1obj", "A1", "W1i", "W1m", "W1t", "W2", "W3"]
    lines = ["| test | " + " | ".join(cands) + " |",
             "|---|" + "---|" * len(cands)]
    V = {c: {} for c in cands}
    for t in tests:
        cells = []
        for c in cands:
            v, why = verdict_test(RUNS, c, t)
            V[c][t] = v
            cells.append("**{}** {}".format(v, why) if v != "-" else "-")
        lines.append("| {} | {} |".format(t, " | ".join(cells)))
    # fallback (user 10-03): every test passes, only A1 converges at first order
    fb = [c for c in cands if all(V[c][t] == "PASS" for t in tests if t not in ("A1", "P1obj", "T1", "T2"))
          and (V[c]["A1"] == "PASS" or ORDER_FALLBACK <= a1_order(RUNS, c) < ORDER_PASS)]
    lines.append("")
    lines.append("Accepted under the 10-03 priorities (P1 strict, waves order 2, T1r, A1 p >= {}): {}".format(
        ORDER_FALLBACK, ", ".join(fb) if fb else "none"))
    return "\n".join(lines)


# ------------------------------------------------------------------ frames

def frames(D, n=6, field="rho"):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    meta = json.load(open(os.path.join(D, "meta.json")))
    snaps = snapshots(os.path.join(D, "output"))
    if not snaps:
        return None
    if meta["test"].startswith("W"):
        r = json.load(open(os.path.join(D, "result.json")))
        fig, ax = plt.subplots(figsize=(6, 3))
        s = load(snaps[-1][1], ("x", "v", "rho"))
        comp = 1 if meta["test"] == "W3" else 0
        ax.plot(s["x"][:, 0], s["v"][:, comp], ",", color="#2a78d6")
        ax.set_title("{} {} t={:.3f}  {}".format(meta["cand"], meta["tag"], snaps[-1][0],
                                                r["summary"][:60]), fontsize=7)
        png = os.path.join(D, "frames.png")
        fig.savefig(png, dpi=130, bbox_inches="tight")
        plt.close(fig)
        return png
    pick = [snaps[int(round(i * (len(snaps) - 1) / (n - 1)))] for i in range(n)]
    fig, axs = plt.subplots(1, n, figsize=(2.6 * n, 2.8), squeeze=False)
    lim = {"rho": (0.8, 1.2), "u": (0.9, 1.1)}.get(field, (None, None))
    for ax, (t, path) in zip(axs[0], pick):
        s = load(path, ("x", field))
        x = s["x"]
        sz = max(0.2, 3000.0 / len(x))
        sc = ax.scatter(x[:, 0], x[:, 1], c=s[field], s=sz, cmap="turbo", vmin=lim[0],
                        vmax=lim[1], linewidths=0)
        ax.set_aspect("equal")
        ax.set_title("t={:.3g}".format(t), fontsize=8)
        ax.tick_params(labelsize=6)
    fig.colorbar(sc, ax=axs[0].tolist(), shrink=0.8, label=field)
    fig.suptitle("{} {} {} ({})".format(meta["cand"], meta["test"], meta["tag"], field),
                 fontsize=9)
    png = os.path.join(D, "frames_{}.png".format(field))
    fig.savefig(png, dpi=130, bbox_inches="tight")
    plt.close(fig)
    return png
