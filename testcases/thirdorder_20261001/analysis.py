#!/usr/bin/env python3
"""Convergence of the 2D ideal-gas soundwave at t = T (default 1).
Usage: ./analysis.py [--t T] tag [tag ...]   (run dirs run_<tag>_N<n>/output)
Metrics (all / A): L1(rho) = <|rho_i - <rho>(1 + A sin k(x_i - ct))|>, L1(v_x), damping 1 - a and |phase| from the
sinus fit of the y-averaged rho profile, plus dE/E and a finite check. Order p between consecutive rungs."""
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
        return {q: h[q][:] for q in ("x", "v", "rho", "m", "u")}
def metrics(d, n, T):
    S = snaps(d); f1 = [fn for t, fn in S if abs(t - T) < 1e-7]
    if not f1: return None
    s0, s1 = load(S[0][1]), load(f1[0]); t = T
    x = s1["x"][:, 0]; rho = s1["rho"]; vx = s1["v"][:, 0]/c; wave = np.sin(k*(x - c*t))
    L1r = np.mean(np.abs(rho - rho.mean()*(1 + A*wave)))/A
    L1v = np.mean(np.abs(vx - A*wave))/A
    col = np.floor(x*n).astype(int) % n
    prof = np.bincount(col, rho, n)/np.bincount(col, minlength=n)
    m = rho.mean(); y = (prof - m)/(m*A); xc = (np.arange(n) + .5)/n
    (p, q), *_ = np.linalg.lstsq(np.column_stack([np.sin(k*(xc-c*t)), np.cos(k*(xc-c*t))]), y, rcond=None)
    E = lambda s: np.sum(s["m"]*(s["u"] + .5*np.sum(s["v"]**2, 1)))
    dE = (E(s1) - E(s0))/abs(E(s0))
    ok = all(np.isfinite(s1[q]).all() for q in ("rho", "v", "u"))
    return L1r, L1v, 1 - np.hypot(p, q), abs(np.arctan2(-q, p)), dE, ok
def main(args):
    T = 1.0
    if args and args[0] == "--t": T = float(args[1]); args = args[2:]
    out = ["t = %g; amplitudes / A; p = log2(err(n/2)/err(n))" % T,
           "%-10s %4s  %10s %10s %11s %10s %9s   p: L1rho L1vx damp phase" % ("arm", "n", "L1(rho)", "L1(vx)", "1-a_rho", "|phase|", "dE/E")]
    for tag in args:
        prev = None
        for n in (16, 32, 64, 128, 256, 512, 1024):
            d = os.path.join(HERE, "run_%s_N%d/output" % (tag, n))
            if not os.path.isdir(d): continue
            r = metrics(d, n, T)
            if r is None: out.append("%-10s %4d  (no t=%g snapshot)" % (tag, n, T)); prev = None; continue
            s = "%-10s %4d  %.4e %.4e %+.4e %.4e %+.1e" % (tag, n, *r[:5]) + ("" if r[5] else " NaN!")
            if prev: s += "   %5.2f %5.2f %5.2f %5.2f" % tuple(np.log2(abs(prev[j]/r[j])) for j in range(4))
            out.append(s); prev = r
        out.append("")
    print("\n".join(out))
    open(os.path.join(HERE, "analysis_%s.txt" % "_".join(args)[:80]), "w").write("\n".join(out) + "\n")
if __name__ == "__main__":
    main(sys.argv[1:])
