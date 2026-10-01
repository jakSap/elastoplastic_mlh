#!/usr/bin/env python3
"""Elastic-branch figure of the thirdorder study (mu = 6/13 > 0, Murnaghan and Tillotson, one period):
(a) longitudinal wave, L1(rho)/A against the O(A^2) reference vs n: today's scheme, the third-order recipe as is
    (tangential traction picked by sign(S*)), with the shear-wave star state (XSPLIT_SHEAR_STAR 2), and the
    fourth-order variant with the stress reconstructed quadratically or linearly;
(b) Fourier spectrum of the y-averaged density error at n = 128 (the side switch's even-harmonic comb);
(c) transverse (shear) wave, L1(v_y)/(c_s A)."""
import os, sys, glob, numpy as np, h5py
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from analysis import metrics, ic_attrs
HERE = os.path.dirname(os.path.abspath(__file__))
SURF, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e4e3df"
C = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948"]


def series(tag):
    ns, e = [], []
    for n in (32, 64, 128, 256, 512):
        d = os.path.join(HERE, "run_%s_N%d/output" % (tag, n))
        if not os.path.isdir(d) or not os.path.exists(os.path.join(HERE, "run_%s_N%d/exit_code" % (tag, n))):
            continue
        r = metrics(d, n, None)
        if r is None: continue
        ns.append(n); e.append(r[6])
    return np.array(ns), np.array(e)


def spectrum(tag, n):
    d = os.path.join(HERE, "run_%s_N%d/output" % (tag, n))
    W = ic_attrs(d); A, c, T = W["A"], W["c"], W["period"]
    for fn in sorted(glob.glob(d + "/*.h5")):
        with h5py.File(fn) as h:
            t = float(np.ravel(h["time"])[0])
            if abs(t - T) > 1e-7: continue
            x = h["x"][:, 0]; rho = h["rho"][:]
    e = (rho - rho.mean()*(1 + A*np.sin(2*np.pi*(x - c*T))))/A
    col = np.floor(x*n).astype(int) % n
    prof = np.bincount(col, e, n)/np.bincount(col, minlength=n)
    return 2*np.abs(np.fft.rfft(prof)/n)


plt.rcParams.update({"font.size": 10, "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2,
                     "ytick.color": INK2, "text.color": INK})
fig, axs = plt.subplots(1, 3, figsize=(15.5, 4.8), facecolor=SURF, gridspec_kw=dict(width_ratios=[1.15, 0.9, 1]))
for ax in axs:
    ax.set_facecolor(SURF); ax.grid(True, which="major", color=GRID, lw=0.8); ax.set_axisbelow(True)
    for s in ("top", "right"): ax.spines[s].set_visible(False)


def ladder(ax, S, ylim, guides):
    n0 = np.array([32, 512])
    for p, y0, lab in guides:
        ax.plot(n0, y0*(n0/32.)**(-p), ls=(0, (4, 3)), color=INK2, lw=1)
        ax.text(n0[1]*1.05, y0*(n0[1]/32.)**(-p), lab, color=INK2, va="center", fontsize=9)
    for tagm, tagt, lab, col in S:
        for tag, filled in ((tagm, True), (tagt, False)):
            if tag is None: continue
            ns, e = series(tag)
            if len(ns) == 0: continue
            ax.plot(ns, e, "-" if filled else "none", color=col, lw=2, marker="o" if filled else "s",
                    ms=7 if filled else 9, mfc=col if filled else "none", mec=SURF if filled else col,
                    mew=1.5, label=lab if filled else None, zorder=3 if filled else 4)
    ax.set_xscale("log", base=2); ax.set_yscale("log")
    ax.set_xticks([32, 64, 128, 256, 512]); ax.set_xticklabels(["32", "64", "128", "256", "512"])
    ax.set_xlim(26, 800); ax.set_ylim(*ylim)
    ax.set_xlabel("particles per wavelength n")


ladder(axs[0], [("e1mc", "e1tc", "today's scheme (GEF 0, KDK, limiters)", C[1]),
                ("t3pm", "t3pt", "third-order recipe as is (t3p)", C[0]),
                ("t3pmx2", "t3ptx2", "t3p + shear-wave star state (XSPLIT_SHEAR_STAR 2)", C[2]),
                ("q23r4mls", None, "RK4 + 2/3-quadratic recon, S linear", C[4]),
                ("q23r4m", "q23r4t", "RK4 + 2/3-quadratic recon incl. S", C[3])],
       (5e-9, 0.1), ((2, 2e-2, "p = 2"), (3, 2.5e-3, "p = 3"), (4, 1.2e-4, "p = 4")))
axs[0].plot([], [], "s", mfc="none", mec=INK2, ms=9, label="hollow squares: Tillotson (filled: Murnaghan)")
axs[0].legend(frameon=False, fontsize=8.2, loc="lower left")
axs[0].set_ylabel("L1(rho) / A after one period (vs O(A^2) reference)")
axs[0].set_title("(a) longitudinal wave, mu = 6/13", loc="left", color=INK, fontsize=11)

ax = axs[1]
kx = np.arange(1, 25)
for tag, lab, col, off in (("t3pm", "t3p as is", C[0], -0.2), ("t3pmx2", "+ XSPLIT_SHEAR_STAR 2", C[2], 0.2)):
    if not os.path.exists(os.path.join(HERE, "run_%s_N128/exit_code" % tag)): continue
    F = spectrum(tag, 128)[1:25]
    ax.vlines(kx + off, 1e-12, F, color=col, lw=2.5, label=lab)
    ax.plot(kx + off, F, "o", color=col, ms=4)
ax.axhline(3.5e-6, color=INK2, ls=(0, (4, 3)), lw=1)
ax.text(24.5, 3.5e-6*1.4, "physical 2nd harmonic", color=INK2, fontsize=8.5, ha="right")
ax.set_yscale("log"); ax.set_ylim(1e-11, 2e-4); ax.set_xlim(0, 25)
ax.set_xlabel("harmonic kx of the y-averaged density error")
ax.set_ylabel("Fourier amplitude / A")
ax.set_title("(b) error spectrum, Murnaghan, n = 128", loc="left", color=INK, fontsize=11)
ax.legend(frameon=False, fontsize=8.5, loc="upper right")

ladder(axs[2], [("se1mc", "se1tc", "today's scheme", C[1]),
                ("st3pm", "st3pt", "t3p as is (side by sign S*)", C[0]),
                ("st3pmx1", None, "t3p + central S_nt (XSPLIT_SHEAR_STAR 1)", C[6]),
                ("st3pmx2", "st3ptx2", "t3p + shear-wave star state (2)", C[2])],
       (5e-7, 0.1), ((2, 1.4e-2, "p = 2"), (3, 3.5e-3, "p = 3")))
axs[2].legend(frameon=False, fontsize=8.2, loc="lower left")
axs[2].set_ylabel("L1(v_y) / (c_s A) after one period")
axs[2].set_title("(c) transverse wave, c_s = sqrt(mu/rho)", loc="left", color=INK, fontsize=11)

fig.suptitle("Elastic branch at third order: Murnaghan / Tillotson, mu = 6/13, xSplitElasticHLLC (GEF 0), fixed CFL",
             x=0.01, ha="left", color=INK, fontsize=11)
fig.tight_layout()
out = os.path.join(HERE, "fig_elastic.png"); fig.savefig(out, dpi=150, facecolor=SURF); print(out)
