#!/usr/bin/env python3
"""Metric-by-metric comparison of wave runs (thirdorder study, elastic part, 2026-10-01): every run gets the same
metrics, whatever study it comes from, so the 2026-09-25 second-order runs and today's ladders sit in one table.

Per run, after one period (wave parameters given explicitly, see RUNS):
  L1     L1(rho)/A against the linear reference          (longitudinal; shear: L1(v_y)/(c A))
  L1nl   L1(rho)/A against the weakly nonlinear O(A^2) reference (analysis.py)
  L1S    L1(S_xx + 4/3 mu ln rho_nl)/(4/3 mu A)          (shear: L1(S_xy + mu A sin)/(mu A))
  damp   1 - fitted amplitude of the y-averaged profile; phase = |fitted phase|
  harm   L1 of the y-averaged error with the fitted fundamental removed (= harmonics / non-smooth part)
  ydep   L1 of the y-dependent part of the error (lattice noise)
  vy     max |transverse velocity| / (c A)  (longitudinal: v_y, shear: v_x)
  dE     |E(T) - E(0)|/E(0);  dp = |sum m v(T) - sum m v(0)| / (M c A)
  wall   wall-clock minutes (run.log)
Usage: ./compare.py [long|shear|big]   -> prints a table and writes compare_<set>.txt"""
import os, sys, glob, re, numpy as np, h5py
HERE = os.path.dirname(os.path.abspath(__file__))
OLD = os.path.expanduser("~/Documents/newMFM/testcases/soundwave_2d/20260925_xsplit_eos")
k = 2*np.pi
MU = 6/13
CL = np.sqrt(1 + 4/3*MU); CS = np.sqrt(MU)
C1 = {"m": 5/42, "t": -4/21}


def snaps(d):
    out = []
    for fn in sorted(glob.glob(os.path.join(d, "*.h5"))):
        with h5py.File(fn) as h: out.append((float(np.ravel(h["time"])[0]), fn))
    return out


def load(fn):
    with h5py.File(fn) as h:
        return {q: h[q][:] for q in ("x", "v", "rho", "m", "u", "Sxx", "Sxy") if q in h}


def wall(run):
    try:
        m = re.findall(r"Total wall-clock time: ([0-9.]+) min", open(os.path.join(run, "run.log")).read())
        return float(m[-1])
    except Exception:
        return np.nan


def metrics(run, n, A, c, c1, mu, T, mode):
    S = snaps(os.path.join(run, "output"))
    f1 = [fn for t, fn in S if abs(t - T) < 1e-6]
    if not f1: return None
    s0, s1 = load(S[0][1]), load(f1[0]); t = T
    x = s1["x"][:, 0]; m = s1["m"]
    col = np.floor(x*n).astype(int) % n
    xi = k*(x - c*t)
    if mode == "shear":
        e = (s1["v"][:, 1] - c*A*np.sin(xi))/(c*A)
        eS = (s1["Sxy"] + mu*A*np.sin(xi))/(mu*A)
        L1 = L1nl = np.mean(np.abs(e)); L1S = np.mean(np.abs(eS)); fld = s1["v"][:, 1]/(c*A); vt = s1["v"][:, 0]
        ref = np.sin(xi)
    else:
        rho = s1["rho"]
        a = 2*c*A*np.sin(xi) + c*(c1 - 1)/2*A**2*np.sin(xi)**2 - (1 + c1)/2*t*(2*c*A*np.sin(xi))*(2*c*A*k*np.cos(xi))
        b = c*(1 - c1)/2*A**2*np.sin(k*(x + c*t))**2
        rbar = np.sum(m)/np.sum(m/rho)
        rho_nl = rbar*(1 + (a - b)/(2*c) + (1 - c1)*a**2/(8*c*c))
        L1 = np.mean(np.abs(rho - rho.mean()*(1 + A*np.sin(xi))))/A
        e = (rho - rho_nl)/A
        L1nl = np.mean(np.abs(e))
        L1S = np.mean(np.abs(s1["Sxx"] + 4/3*mu*np.log(rho_nl)))/(4/3*mu*A) if (mu > 0 and "Sxx" in s1) else np.nan
        fld = (rho - rho.mean())/(rho.mean()*A); vt = s1["v"][:, 1]
        ref = np.sin(xi)
    # fitted mode of the y-averaged field, harmonics = profile error minus the fitted-mode error
    xc = (np.arange(n) + .5)/n
    prof = np.bincount(col, fld, n)/np.bincount(col, minlength=n)
    B = np.column_stack([np.sin(k*(xc - c*t)), np.cos(k*(xc - c*t))])
    (p, q), *_ = np.linalg.lstsq(B, prof, rcond=None)
    eprof = np.bincount(col, e, n)/np.bincount(col, minlength=n)
    (pe, qe), *_ = np.linalg.lstsq(B, eprof, rcond=None)
    harm = np.mean(np.abs(eprof - B @ np.array([pe, qe])))
    ydep = np.mean(np.abs(e - eprof[col]))
    E = lambda s: np.sum(s["m"]*(s["u"] + .5*np.sum(s["v"]**2, 1)))
    dE = abs(E(s1) - E(s0))/abs(E(s0))
    P = lambda s: np.sum(s["m"][:, None]*s["v"], 0)
    dp = np.max(np.abs(P(s1) - P(s0)))/(np.sum(m)*c*A)
    return dict(L1=L1, L1nl=L1nl, L1S=L1S, damp=1 - np.hypot(p, q), phase=abs(np.arctan2(-q, p)), harm=harm,
                ydep=ydep, vy=np.abs(vt).max()/(c*A), dE=dE, dp=dp, wall=wall(run))


def R(tag, n, eos, mode="long", A=1e-6, base=HERE, name=None):
    c = CS if mode == "shear" else CL
    run = os.path.join(base, name or "run_%s_N%d" % (tag, n))
    return (tag, n, run, A, c, C1[eos], MU, 1/c, mode)


SETS = {
    "long": [R(t, n, e) for t, e in (("e1mc", "m"), ("t3pm", "m"), ("t3pmx1", "m"), ("t3pmx2", "m"), ("q23r4m", "m"))
             for n in (64, 128, 256, 512)]
            # the 2026-09-25 second-order binaries (build_xs_murn/till: DRIFT_KDK, GO 1, limiters), re-run here on the
            # mass-corrected ICs (their own snapshots were deleted on 09-28)
            + [R("0925m", n, "m") for n in (64, 128, 256)] + [R("0925t", n, "t") for n in (64, 128, 256)]
            + [R("e1mcx2", n, "m") for n in (64, 128, 256)],
    "shear": [R(t, n, "m", "shear") for t in ("se1mc", "se1mcx2", "st3pm", "st3pmx1", "st3pmx2") for n in (32, 64, 128, 256)],
    "big": [R("%s_%s" % (t, a), n, "m", A=amp) for a, amp in (("A2", 1e-2), ("A1", 1e-1))
            for t in ("e1mc", "t3pmx2", "q23r4m") for n in (64, 128)],
}

if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "long"
    keys = ("L1", "L1nl", "L1S", "damp", "phase", "harm", "ydep", "vy", "dE", "dp", "wall")
    out = ["%-12s %4s " % ("run", "n") + " ".join("%9s" % kk for kk in keys)]
    for tag, n, run, A, c, c1, mu, T, mode in SETS[which]:
        if not os.path.exists(os.path.join(run, "output")): continue
        r = metrics(run, n, A, c, c1, mu, T, mode)
        if r is None: continue
        out.append("%-12s %4d " % (tag, n) + " ".join("%9.2e" % r[kk] for kk in keys))
    print("\n".join(out))
    open(os.path.join(HERE, "compare_%s.txt" % which), "w").write("\n".join(out) + "\n")
