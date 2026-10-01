#!/usr/bin/env python3
"""Convergence of the 2D soundwave (ideal gas) or elastic wave (generateIC_elastic.py, mu > 0) after one period.
Usage: ./analysis.py [--t T] tag [tag ...]   (run dirs run_<tag>_N<n>/output; T defaults to the IC's period)
Metrics (all / A): L1(rho) = <|rho_i - <rho>(1 + A sin k(x_i - ct))|>, L1(v_x), damping 1 - a and |phase| from the
sinus fit of the y-averaged rho profile, plus dE/E and a finite check. Order p between consecutive rungs.
c, A, the period and the nonlinearity c1 = (rho/c) dc/drho come from the IC attributes (ideal gas: c = 1,
c1 = (gamma-1)/2). Elastic runs add L1S = <|S_xx + 4/3 mu ln(rho_nl/rho0)|>/(4/3 mu A) against the O(A^2) reference."""
import numpy as np, h5py, os, sys, glob
HERE = os.path.dirname(os.path.abspath(__file__))
A = 1e-6; k = 2*np.pi; c = 1.0
def snaps(d):
    out = []
    for fn in sorted(glob.glob(os.path.join(d, "*.h5"))):
        with h5py.File(fn) as h: out.append((float(np.ravel(h['time'])[0]), fn))
    return out
def load(fn):
    with h5py.File(fn) as h:
        return {q: h[q][:] for q in ("x", "v", "rho", "m", "u", "Sxx", "Sxy") if q in h}
def ic_attrs(d):
    """Wave parameters from the run's IC: A, c, c1, period, mu (ideal-gas defaults if absent)."""
    out = dict(A=1e-6, c=1.0, c1=(5/3 - 1)/2, period=1.0, mu=0.0, mode="long")
    try:
        cfg = open(os.path.join(os.path.dirname(d.rstrip("/")), "config.info")).read()
        ic = [l.split()[1] for l in cfg.splitlines() if l.startswith("initFile")][0]
        with h5py.File(os.path.join(HERE, ic)) as h:
            at = h.attrs
            out["A"] = float(at.get("amplitude", 1e-6))
            if str(at.get("mode", "long")) == "shear":
                out.update(c=float(at["c_s"]), c1=0., period=float(at["period"]), mu=float(at["mu"]), mode="shear")
            elif "c_l" in at:
                out.update(c=float(at["c_l"]), c1=float(at["c1"]), period=float(at["period"]), mu=float(at["mu"]))
    except Exception:
        pass
    return out
def metrics(d, n, T):
    global A, c
    W = ic_attrs(d); A, c, c1, mu = W["A"], W["c"], W["c1"], W["mu"]
    if T is None: T = W["period"]
    S = snaps(d); f1 = [fn for t, fn in S if abs(t - T) < 1e-7]
    if not f1: return None
    s0, s1 = load(S[0][1]), load(f1[0]); t = T
    if W["mode"] == "shear": return metrics_shear(s0, s1, n, t, mu)
    x = s1["x"][:, 0]; rho = s1["rho"]; vx = s1["v"][:, 0]/c; wave = np.sin(k*(x - c*t))
    L1r = np.mean(np.abs(rho - rho.mean()*(1 + A*wave)))/A
    L1v = np.mean(np.abs(vx - A*wave))/A
    col = np.floor(x*n).astype(int) % n
    prof = np.bincount(col, rho, n)/np.bincount(col, minlength=n)
    m = rho.mean(); y = (prof - m)/(m*A); xc = (np.arange(n) + .5)/n
    (p, q), *_ = np.linalg.lstsq(np.column_stack([np.sin(k*(xc-c*t)), np.cos(k*(xc-c*t))]), y, rcond=None)
    E = lambda s: np.sum(s["m"]*(s["u"] + .5*np.sum(s["v"]**2, 1)))
    dE = (E(s1) - E(s0))/abs(E(s0))
    ok = all(np.isfinite(s1[q_]).all() for q_ in ("rho", "v", "u"))
    # Weakly nonlinear reference to O(A^2) via the Riemann invariants J+- = v +- h(rho), h = int c/rho drho, of any
    # barotropic 1D flow (ideal gas on the isentrope; elastic: Pi = P + 4/3 mu ln rho). With c = c0 (1 + c1 delta),
    # delta = rho/rhobar - 1: h = c0 (delta + (c1-1) delta^2/2), wave speeds v +- c = +-c0 + a (1 +- c1)/2 + ...
    # The forward part a = J+ - J+bar steepens, a = a0(xi) - (1+c1)/2 t a0 a0'(xi), xi = x - c t; the IC (linear
    # eigenmode, not an exact simple wave) launches a backward O(A^2) wave b0 = c (1-c1)/2 A^2 sin^2(kx), moving
    # with -c. rho/rhobar = 1 + (a-b)/(2c) + (1-c1) a^2/(8 c^2), v = (a+b)/2; rhobar = total mass / total volume.
    # Ideal gas: c1 = (gamma-1)/2, i.e. (gamma+1)/4, (3-gamma)/4 and (3-gamma)/16 as before.
    s_ = np.sin(k*(x - c*t)); c_ = np.cos(k*(x - c*t))
    a = 2*c*A*s_ + c*(c1 - 1)/2*A**2*s_**2 - (1 + c1)/2*t*(2*c*A*s_)*(2*c*A*k*c_)
    b = c*(1 - c1)/2*A**2*np.sin(k*(x + c*t))**2
    rbar = np.sum(s1["m"])/np.sum(s1["m"]/rho)
    rho_nl = rbar*(1 + (a - b)/(2*c) + (1 - c1)*a**2/(8*c*c))
    L1nl = np.mean(np.abs(rho - rho_nl))/A
    L1vnl = np.mean(np.abs(s1["v"][:, 0] - (a + b)/2))/(c*A)
    L1S = np.nan
    if mu > 0:   # slaved deviator S_xx = -4/3 mu ln(rho/rho0), rho0 = 1
        L1S = np.mean(np.abs(s1["Sxx"] + 4/3*mu*np.log(rho_nl)))/(4/3*mu*A)
    return L1r, L1v, 1 - np.hypot(p, q), abs(np.arctan2(-q, p)), dE, ok, L1nl, L1vnl, L1S
def metrics_shear(s0, s1, n, t, mu):
    """Transverse wave: columns L1(rho) -> L1(v_y)/(c_s A), L1(v_x) -> L1(S_xy)/(mu A), damping/phase of the
    y-averaged v_y profile; the L1nl column repeats L1(v_y) (linear to O(A^2)), L1S = L1(S_xy)."""
    x = s1["x"][:, 0]; wave = np.sin(k*(x - c*t))
    L1vy = np.mean(np.abs(s1["v"][:, 1] - c*A*wave))/(c*A)
    L1Sxy = np.mean(np.abs(s1["Sxy"] + mu*A*wave))/(mu*A)
    col = np.floor(x*n).astype(int) % n
    y = np.bincount(col, s1["v"][:, 1], n)/np.bincount(col, minlength=n)/(c*A); xc = (np.arange(n) + .5)/n
    (p, q), *_ = np.linalg.lstsq(np.column_stack([np.sin(k*(xc-c*t)), np.cos(k*(xc-c*t))]), y, rcond=None)
    E = lambda s: np.sum(s["m"]*(s["u"] + .5*np.sum(s["v"]**2, 1)))
    dE = (E(s1) - E(s0))/abs(E(s0))
    ok = all(np.isfinite(s1[q_]).all() for q_ in ("rho", "v", "u", "Sxy"))
    return L1vy, L1Sxy, 1 - np.hypot(p, q), abs(np.arctan2(-q, p)), dE, ok, L1vy, L1vy, L1Sxy
def main(args):
    T = None
    if args and args[0] == "--t": T = float(args[1]); args = args[2:]
    out = ["t = %s; amplitudes / A; p = log2(err(n/2)/err(n))" % ("one period" if T is None else "%g" % T),
           "%-12s %4s  %10s %10s %11s %10s %9s %10s %10s   p: L1rho L1vx damp phase L1nl L1S" % ("arm", "n", "L1(rho)", "L1(vx)", "1-a_rho", "|phase|", "dE/E", "L1nl(rho)", "L1S")]
    for tag in args:
        prev = None
        for n in (16, 32, 64, 128, 256, 512, 1024):
            d = os.path.join(HERE, "run_%s_N%d/output" % (tag, n))
            if not os.path.isdir(d): continue
            r = metrics(d, n, T)
            if r is None: continue
            s = "%-12s %4d  %.4e %.4e %+.4e %.4e %+.1e %.4e %.4e" % (tag, n, *r[:5], r[6], r[8]) + ("" if r[5] else " NaN!")
            if prev: s += "   %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f" % tuple(np.log2(abs(prev[j]/r[j])) for j in (0, 1, 2, 3, 6, 8))
            out.append(s); prev = r
        out.append("")
    print("\n".join(out))
    open(os.path.join(HERE, "analysis_%s.txt" % "_".join(args)[:80]), "w").write("\n".join(out) + "\n")
if __name__ == "__main__":
    main(sys.argv[1:])
