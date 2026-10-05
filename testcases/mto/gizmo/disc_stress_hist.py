#!/usr/bin/env python3
"""Stress-magnitude histograms of the GIZMO rotating-disc series (runs/disc_<rung>_f0.5).

    ./disc_stress_hist.py OUTDIR [--prefix 2026-10-05-gizmo-scheibe-stress]

Per run: a movie (one frame per snapshot) of the histogram of |S| = Frobenius norm of the 3x3 StressTensor
per particle, with a least-squares Gauss fit A exp(-(x-mu)^2 / 2 sigma^2) to the density histogram.
Next to it, the radial profile of |S| (mean |S| per radial shell, NSHELL shells on [0, RMAX]) with a Gauss fit
B exp(-(r-r0)^2 / 2 s^2); r0 is the radial position of the maximum of the fit.
Over all runs: histogram-fit height A(t), and radial position r0(t).
Bins are common to all runs and frames (80 bins on [0, 1.02 max|S|]), so the heights are comparable.
Frames whose spread is below one bin width are not resolved by the histogram (the fit is meaningless there):
they are drawn as open markers. For r0, open markers mean r0 lies outside the disc (profile rises to the edge).
"""
import argparse, glob, os, sys
import h5py
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import animation
from scipy.optimize import curve_fit

HERE = os.path.dirname(os.path.abspath(__file__))
OFFSET = 5.0
NBINS = 80
NSHELL = 25
RMAX = 1.0125


def load(run):
    out = []
    for s in sorted(glob.glob(os.path.join(HERE, "runs", run, "output", "snapshot_*.hdf5"))):
        with h5py.File(s, "r") as f:
            g = f["PartType0"]
            S = g["StressTensor"][:].reshape(-1, 9)
            xy = g["Coordinates"][:][:, :2] - OFFSET
            out.append((float(f["Header"].attrs["Time"]), np.linalg.norm(S, axis=1), np.hypot(xy[:, 0], xy[:, 1])))
    return out


def gauss(x, A, mu, sig):
    return A * np.exp(-0.5 * ((x - mu) / sig) ** 2)


def fit(dens, centres, w):
    """Least-squares Gauss fit to the density histogram. Returns (A, mu, sigma, resolved)."""
    mu0 = np.sum(dens * centres) / np.sum(dens)
    sig0 = max(np.sqrt(np.sum(dens * (centres - mu0) ** 2) / np.sum(dens)), w / 2)
    try:
        p, _ = curve_fit(gauss, centres, dens, p0=[dens.max(), mu0, sig0], maxfev=20000)
        p[2] = abs(p[2])
    except RuntimeError:
        p = np.array([np.nan, np.nan, np.nan])
    return p[0], p[1], p[2], sig0 > w


def fit_radial(prof, rc):
    """Gauss fit to the shell-mean profile. Returns (B, r0, sigma, inside); NaN if the profile is ~0 (t = 0)."""
    if not np.nanmax(prof) > 1e-8:
        return np.nan, np.nan, np.nan, False
    m = np.isfinite(prof)
    i = np.nanargmax(prof)
    try:
        p, _ = curve_fit(gauss, rc[m], prof[m], p0=[prof[i], rc[i], 0.3], maxfev=20000)
        p[2] = abs(p[2])
    except RuntimeError:
        return np.nan, np.nan, np.nan, False
    return p[0], p[1], p[2], 0 <= p[1] <= RMAX


def profile(s, rad, redges):
    idx = np.digitize(rad, redges) - 1
    return np.array([s[idx == i].mean() if np.any(idx == i) else np.nan for i in range(len(redges) - 1)])


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("outdir")
    ap.add_argument("--prefix", default="2026-10-05-gizmo-scheibe-stress")
    a = ap.parse_args()

    runs = sorted((os.path.basename(d) for d in glob.glob(os.path.join(HERE, "runs", "disc_d0.0*_f0.5"))),
                  key=lambda r: -float(r.split("_")[1][1:]))
    data = {r: load(r) for r in runs}
    smax = 1.02 * max(s.max() for r in runs for _, s, _ in data[r])
    edges = np.linspace(0, smax, NBINS + 1)
    centres = 0.5 * (edges[1:] + edges[:-1])
    w = edges[1] - edges[0]
    redges = np.linspace(0, RMAX, NSHELL + 1)
    rc = 0.5 * (redges[1:] + redges[:-1])

    res = {}
    for r in runs:
        rows = []
        for t, s, rad in data[r]:
            dens, _ = np.histogram(s, bins=edges, density=True)
            A, mu, sig, ok = fit(dens, centres, w)
            B, r0, sr, inside = fit_radial(profile(s, rad, redges), rc)
            rows.append((t, A, mu, sig, ok, r0, len(s), B, sr, inside))
        res[r] = np.array(rows, dtype=float)

    # ---- one movie per run
    for r in runs:
        d0 = r.split("_")[1][1:]
        fig, (ax, ar) = plt.subplots(1, 2, figsize=(13, 4.8), dpi=120)

        def frame(k, r=r, d0=d0):
            ax.clear(); ar.clear()
            t, s, rad = data[r][k]
            dens, _ = np.histogram(s, bins=edges, density=True)
            A, mu, sig, ok = res[r][k, 1:5]
            ax.bar(centres, dens, width=w, color="#4c78a8", alpha=0.8, label="Histogramm $|S|$")
            if np.isfinite(A):
                xx = np.linspace(0, smax, 600)
                ax.plot(xx, gauss(xx, A, mu, sig), color="#e45756", lw=2,
                        label="Gauß-Fit: $A$=%.0f, $\\mu$=%.2e, $\\sigma$=%.2e%s" % (A, mu, sig, "" if ok else " (nicht aufgelöst)"))
            ax.set_xlim(0, smax)
            ax.set_xlabel("$|S|$ = Frobenius-Norm des Spannungstensors")
            ax.set_ylabel("Dichte")
            ax.set_title("GIZMO Scheibe d0=%s (N=%d)   t = %.2f  (T/2 = 15.71)" % (d0, len(s), t))
            ax.legend(loc="upper right", fontsize=8)
            ax.grid(alpha=0.25)
            prof = profile(s, rad, redges)
            B, sr, inside = res[r][k, 7:10]
            r0 = res[r][k, 5]
            ar.plot(rc, prof, "o", color="#4c78a8", label="mittleres $|S|$ je Radiusschale")
            if np.isfinite(B):
                rr = np.linspace(0, 1.3, 400)
                ar.plot(rr, gauss(rr, B, r0, sr), color="#e45756", lw=2,
                        label="Gauß-Fit: $r_0$=%.3f, $\\sigma$=%.3f%s" % (r0, sr, "" if inside else " ($r_0$ außerhalb)"))
                ar.axvline(r0, color="#e45756", ls=":", lw=1)
            ar.set_xlim(0, 1.3)
            ar.set_xlabel("Radius $r$")
            ar.set_ylabel("mittleres $|S|$")
            ar.set_title("Radialprofil von $|S|$")
            ar.legend(loc="upper left", fontsize=8)
            ar.grid(alpha=0.25)
            fig.tight_layout()

        out = os.path.join(a.outdir, "%s-d%s.mp4" % (a.prefix, d0.replace(".", "")))
        ani = animation.FuncAnimation(fig, frame, frames=len(data[r]))
        ani.save(out, writer=animation.FFMpegWriter(fps=2, codec="libx264", extra_args=["-pix_fmt", "yuv420p"]))
        plt.close(fig)
        print("wrote", out)

    # ---- combined plots
    cmap = plt.get_cmap("viridis")
    col = {r: cmap(i / (len(runs) - 1)) for i, r in enumerate(runs)}

    def lab(r):
        return "d0=%s (N=%d)" % (r.split("_")[1][1:], res[r][0, 6])

    def draw(ax, key, ylab):
        for r in runs:
            R = res[r]
            ok = R[:, 4] > 0
            ax.plot(R[:, 0], R[:, key], "-", color=col[r], lw=1.2, alpha=0.8, label=lab(r))
            ax.plot(R[ok, 0], R[ok, key], "o", color=col[r], ms=5)
            ax.plot(R[~ok, 0], R[~ok, key], "o", mfc="none", color=col[r], ms=5)
        ax.set_ylabel(ylab)
        ax.grid(alpha=0.25)

    fig, ax = plt.subplots(figsize=(8, 5), dpi=120)
    draw(ax, 1, "Höhe $A$ des Gauß-Fits (Dichte)")
    ax.set_yscale("log")
    ax.set_xlabel("$t$"); ax.legend(fontsize=8, title="offen: Breite < 1 Bin")
    ax.set_title("Höhe der Gauß-Glocke von $|S|$ über der Zeit")
    fig.tight_layout(); p = os.path.join(a.outdir, a.prefix + "-hoehe.png"); fig.savefig(p); print("wrote", p)

    fig, ax = plt.subplots(figsize=(8, 5), dpi=120)
    for r in runs:
        R = res[r]
        good = np.isfinite(R[:, 5])
        ax.plot(R[good, 0], R[good, 5], "-", color=col[r], lw=1.2, alpha=0.8, label=lab(r))
        ins = good & (R[:, 9] > 0)
        ax.plot(R[ins, 0], R[ins, 5], "o", color=col[r], ms=5)
        ax.plot(R[good & ~ins, 0], R[good & ~ins, 5], "o", mfc="none", color=col[r], ms=5)
    ax.axhline(RMAX, color="k", lw=0.8, ls="--"); ax.text(0.2, RMAX, " Scheibenrand", va="bottom", fontsize=8)
    ax.set_ylim(0.5, 2.0)
    ax.set_xlabel("$t$"); ax.set_ylabel("radiale Lage $r_0$ des Gauß-Maximums von $|S|(r)$")
    ax.grid(alpha=0.25); ax.legend(fontsize=8, title="offen: $r_0$ außerhalb der Scheibe", loc="center left", bbox_to_anchor=(1.01, 0.5))
    ax.set_title("Radiale Lage des Maximums der Gauß-Glocke über der Zeit")
    fig.set_size_inches(10.5, 5); fig.tight_layout(); p = os.path.join(a.outdir, a.prefix + "-radius.png"); fig.savefig(p); print("wrote", p)

    print("\nrung        t_end   A      mu       sigma    r0     sigma_r")
    for r in runs:
        t, A, mu, sig, ok, r0, n, B, sr, ins = res[r][-1]
        print("%-8s %7.2f %7.1f %.3e %.3e %.3f %.3f%s" % (r.split("_")[1], t, A, mu, sig, r0, sr, "" if ins else "  (außerhalb)"))


if __name__ == "__main__":
    main()
