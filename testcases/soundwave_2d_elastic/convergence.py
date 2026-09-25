#!/usr/bin/env python3
"""Convergence analysis of the 2D elastic travelling-wave runs.

For each resolution, at every wave passing t = T, 2T, 3T (T = lambda/c_l,
snapshots selected by their stored /time, never by file index -- see the
dump-index caveat in makeConfigs.py) the exact solution is the IC again:

  L1_rho    1/N sum |rho_i - <rho>(1 + A sin(k(x_i - c_l t)))|
  L1_vx     1/N sum |vx_i - c_l A sin(k(x_i - c_l t))|
  L1_Sxx    1/N sum |Sxx_i + (4/3) mu A sin(k(x_i - c_l t))|

and from a sine fit of the y-averaged profiles, y ~ a sin(k x - k c_l t - phi):

  a_rho, a_Sxx   amplitude in units of A resp. (4/3) mu A
  phi            phase shift (+ = ahead of c_l)
  c_num/c_l - 1  measured wave speed relative to sqrt((K + 4/3 mu)/rho0)

Plus two checks specific to mu > 0:

  slaving    max |S_xx + (4/3) mu ln(rho/rho0)| / ((4/3) mu A) -- the hypo-elastic
             stress must stay slaved to the density for this mode
  L1_Sxy     1/N sum |S_xy|, which must stay at round-off (no shear is excited)

Writes plots/L1_vs_N.png, plots/amp_phase_vs_t.png, plots/residual.png,
plots/summary.txt and plots/results.json.
"""

import json
import os
import numpy as np
import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from makeMovie import snapshots, ic_name

RESOLUTIONS = (64, 128, 256)
PASSINGS = (1, 2, 3)

SURFACE, INK, INK_2, MUTED, GRID, AXIS = "#fcfcfb", "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7"
# resolution is ORDINAL -> single-hue blue ramp (250/450/650), as in soundwave_2d
RES_COLOR = {64: "#86b6ef", 128: "#2a78d6", 256: "#104281"}
SERIES_2 = "#eb6834"

plt.rcParams.update({
    "font.family": "sans-serif", "font.size": 9, "text.color": INK,
    "axes.labelcolor": INK_2, "axes.edgecolor": AXIS, "axes.linewidth": 0.8,
    "xtick.color": MUTED, "ytick.color": MUTED, "xtick.labelcolor": INK_2, "ytick.labelcolor": INK_2,
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
})


def stem(n):
    return "elastic_N{}_corr".format(n)


def meta(n):
    with h5py.File(ic_name(stem(n)), "r") as f:
        return {k: float(f.attrs[k]) for k in
                ("amplitude", "wavelength", "mu", "c_l", "period", "kernel_offset")}


def per_snapshot(path, md, n):
    """All diagnostics of one snapshot."""
    A, k, mu, c_l = md["amplitude"], 2 * np.pi / md["wavelength"], md["mu"], md["c_l"]
    S_SCALE = (4.0 / 3.0) * mu * A
    with h5py.File(path, "r") as d:
        t = float(np.ravel(d["time"])[0])
        x, vx = d["x"][:, 0], d["v"][:, 0]
        rho, Sxx, Sxy, Syy = d["rho"][:], d["Sxx"][:], d["Sxy"][:], d["Syy"][:]
        cons = {key: float(np.ravel(d[key])[0])
                for key in ("totalMass", "energy", "xMomentum", "yMomentum")}
    phase = k * c_l * t
    rho_m = rho.mean()
    wave = np.sin(k * x - phase)

    out = dict(t=t, rho_mean=rho_m, **cons)
    out["L1_rho"] = np.mean(np.abs(rho - rho_m * (1.0 + A * wave)))
    out["L1_vx"] = np.mean(np.abs(vx - c_l * A * wave))
    out["L1_Sxx"] = np.mean(np.abs(Sxx + S_SCALE * wave))
    out["L1_Sxy"] = np.mean(np.abs(Sxy))
    # the hypo-elastic stress must stay slaved to the density
    out["slaving"] = np.max(np.abs(Sxx + (4.0 / 3.0) * mu * np.log(rho / 1.0))) / S_SCALE
    # deviatoric trace check: S_yy = -S_xx/2 for this mode
    out["trace_dev"] = np.max(np.abs(Syy + 0.5 * Sxx)) / S_SCALE

    # y-averaged column profiles + sine fits
    col = np.floor(x * n).astype(int) % n
    cnt = np.maximum(np.bincount(col, minlength=n), 1)
    xc = (np.arange(n) + 0.5) / n
    M = np.column_stack([np.sin(k * xc - phase), np.cos(k * xc - phase)])
    for name, vals, scale in (("rho", (np.bincount(col, weights=rho, minlength=n) / cnt - rho_m)
                                      / (rho_m * A), 1.0),
                              ("Sxx", -np.bincount(col, weights=Sxx, minlength=n) / cnt / S_SCALE, 1.0)):
        (p, q), *_ = np.linalg.lstsq(M, vals * scale, rcond=None)
        out["a_" + name] = float(np.hypot(p, q))
        out["phi_" + name] = float(np.arctan2(-q, p))
    out["prof_x"] = xc
    out["prof_res"] = ((np.bincount(col, weights=rho, minlength=n) / cnt - rho_m)
                       / (rho_m * A) - np.sin(k * xc - phase))
    return out


def analyse(n):
    md = meta(n)
    files, times = snapshots("output_" + stem(n))
    if not files:
        return None
    series = [per_snapshot(f, md, n) for f in files]
    t0 = series[0]
    # pick the snapshot closest to each period boundary
    at = {}
    for p in PASSINGS:
        target = p * md["period"]
        i = int(np.argmin([abs(s["t"] - target) for s in series]))
        if abs(series[i]["t"] - target) > 1e-9 * max(1.0, target):
            print("  warning: n = {}, {} T: nearest snapshot is t = {:.9f}, target {:.9f}"
                  .format(n, p, series[i]["t"], target))
        at[p] = series[i]
    return dict(meta=md, series=series, at=at, t0=t0)


def order(vals, ns):
    return [np.log(vals[i] / vals[i + 1]) / np.log(ns[i + 1] / ns[i]) for i in range(len(ns) - 1)]


if __name__ == "__main__":
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    os.makedirs("plots", exist_ok=True)

    res = {n: analyse(n) for n in RESOLUTIONS}
    have = [n for n in RESOLUTIONS if res[n]]
    if not have:
        raise SystemExit("no runs found")
    md = res[have[0]]["meta"]
    T, c_l, mu, A = md["period"], md["c_l"], md["mu"], md["amplitude"]
    k = 2 * np.pi / md["wavelength"]

    lines = []
    def out(s=""):
        print(s)
        lines.append(s)

    out("2D elastic travelling wave, Tillotson (K = 1), mu = {:.9g} (Poisson 0.3)".format(mu))
    out("c_l = sqrt((K + 4/3 mu)/rho0) = {:.9f}, period T = {:.9f}, A = {:g}".format(c_l, T, A))
    out("masses corrected for the kernel-sum offset {:+.6e}".format(md["kernel_offset"]))
    out()

    for p in PASSINGS:
        out("--- t = {} T = {:.6f} ({} wave passing{}) ---".format(p, p * T, p, "" if p == 1 else "s"))
        out("  {:>5s} {:>11s} {:>11s} {:>11s} {:>9s} {:>9s} {:>11s} {:>11s}".format(
            "n", "L1_rho", "L1_vx", "L1_Sxx", "a_rho", "a_Sxx", "c_num/c_l-1", "slaving"))
        for n in have:
            s = res[n]["at"][p]
            cnum = s["phi_rho"] / (k * s["t"]) / c_l
            out("  {:>5d} {:11.4e} {:11.4e} {:11.4e} {:9.6f} {:9.6f} {:+11.3e} {:11.3e}".format(
                n, s["L1_rho"], s["L1_vx"], s["L1_Sxx"], s["a_rho"], s["a_Sxx"], cnum, s["slaving"]))
        if len(have) > 1:
            for key in ("L1_rho", "L1_vx", "L1_Sxx"):
                o = order([res[n]["at"][p][key] for n in have], have)
                out("  order {:<8s} {}".format(key, "  ".join("{:.2f}".format(v) for v in o)))
        out()

    out("--- conservation and shear cleanliness (over the whole run) ---")
    out("  {:>5s} {:>12s} {:>12s} {:>12s} {:>11s} {:>11s}".format(
        "n", "dM/M", "dE/E", "dPx", "max L1_Sxy", "max trace"))
    for n in have:
        s0, ser = res[n]["t0"], res[n]["series"]
        dM = max(abs(s["totalMass"] / s0["totalMass"] - 1) for s in ser)
        dE = max(abs(s["energy"] / s0["energy"] - 1) for s in ser)
        dP = max(abs(s["xMomentum"] - s0["xMomentum"]) for s in ser)
        out("  {:>5d} {:12.3e} {:12.3e} {:12.3e} {:11.3e} {:11.3e}".format(
            n, dM, dE, dP, max(s["L1_Sxy"] for s in ser), max(s["trace_dev"] for s in ser)))
    out()

    # ---------------- plots ----------------
    # 1. L1 vs n, one line per passing
    fig, axes = plt.subplots(1, 3, figsize=(10.5, 3.6))
    for ax, key, lab in zip(axes, ("L1_rho", "L1_vx", "L1_Sxx"),
                            (r"$L_1(\rho)$", r"$L_1(v_x)$", r"$L_1(S_{xx})$")):
        # passings are ORDINAL too -> same single-hue blue ramp, light = 1 T
        for p, col, mk in zip(PASSINGS, ("#86b6ef", "#2a78d6", "#104281"), ("o", "s", "^")):
            y = [res[n]["at"][p][key] for n in have]
            ax.plot(have, y, "-", marker=mk, ms=5, lw=2, color=col,
                    label="{} T".format(p))
        # order-2 guide, offset below the data: at exactly order 2 it would
        # otherwise sit on top of the 3 T line and read as part of it
        ref = 0.45 * res[have[0]]["at"][PASSINGS[-1]][key]
        guide = [ref * (have[0] / n) ** 2 for n in have]
        ax.plot(have, guide, color=AXIS, lw=1.2, zorder=1)
        ax.text(have[len(have) // 2], guide[len(have) // 2], "slope −2\n", color=MUTED,
                fontsize=8, va="center", ha="center", rotation=-38,
                rotation_mode="anchor")
        ax.set_xscale("log"); ax.set_yscale("log")
        # the log locator would add "6x10^1" style minor labels on top of ours
        ax.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
        ax.set_xticks(have); ax.set_xticklabels([str(n) for n in have])
        ax.set_xlabel("particles per wavelength $n$")
        ax.set_title(lab, color=INK_2, fontsize=9, loc="left")
        ax.grid(True, which="both", color=GRID, lw=0.6); ax.set_axisbelow(True)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
    axes[0].legend(frameon=False, fontsize=8, labelcolor=INK_2, title="after", title_fontsize=8)
    fig.suptitle("Elastic wave: convergence after 1, 2 and 3 wave passings",
                 x=0.02, ha="left", fontsize=11, color=INK)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig("plots/L1_vs_N.png", dpi=170)
    plt.close(fig)

    # 2. amplitude and phase vs t
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 5.4), sharex=True,
                                   gridspec_kw={"hspace": 0.12})
    for n in have:
        ser = res[n]["series"]
        tt = np.array([s["t"] for s in ser]) / T
        ax1.plot(tt, [s["a_rho"] for s in ser], lw=1.8, color=RES_COLOR[n], label="n = {}".format(n))
        ax1.plot(tt, [s["a_Sxx"] for s in ser], lw=1.0, ls="--", color=RES_COLOR[n])
        ax2.plot(tt, [s["phi_rho"] for s in ser], lw=1.8, color=RES_COLOR[n])
    ax1.set_ylabel("amplitude  [$A$]")
    ax1.text(0.99, 0.06, "solid: $\\rho$    dashed: $S_{xx}$", transform=ax1.transAxes,
             fontsize=8, color=MUTED, ha="right")
    ax2.set_ylabel("phase shift  [rad]\n(+ = ahead of $c_l$)")
    ax2.set_xlabel("wave passings  $t/T$")
    for ax in (ax1, ax2):
        ax.grid(True, color=GRID, lw=0.6); ax.set_axisbelow(True)
        ax.axhline(1.0 if ax is ax1 else 0.0, color=AXIS, lw=0.8)
        for p in PASSINGS:
            ax.axvline(p, color=AXIS, lw=0.8, ls=":")
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
    ax1.legend(frameon=False, fontsize=8, labelcolor=INK_2, loc="lower left")
    fig.suptitle("Elastic wave: amplitude and phase over 3 passings",
                 x=0.02, ha="left", fontsize=11, color=INK)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig("plots/amp_phase_vs_t.png", dpi=170)
    plt.close(fig)

    # 3. density residual at 3 T
    fig, ax = plt.subplots(figsize=(8, 3.4))
    for n in have:
        s = res[n]["at"][PASSINGS[-1]]
        ax.plot(s["prof_x"], s["prof_res"], lw=1.4, color=RES_COLOR[n], label="n = {}".format(n))
    ax.set_xlim(0, 1)
    ax.set_xlabel("x"); ax.set_ylabel("numerical − analytic  [$A$]")
    ax.axhline(0, color=AXIS, lw=0.8)
    ax.grid(True, color=GRID, lw=0.6); ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.legend(frameon=False, fontsize=8, labelcolor=INK_2)
    fig.suptitle("Density residual after 3 wave passings", x=0.02, ha="left",
                 fontsize=11, color=INK)
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    fig.savefig("plots/residual.png", dpi=170)
    plt.close(fig)

    with open("plots/summary.txt", "w") as f:
        f.write("\n".join(lines) + "\n")
    with open("plots/results.json", "w") as f:
        json.dump({str(n): {"meta": res[n]["meta"],
                            "at": {str(p): {kk: vv for kk, vv in res[n]["at"][p].items()
                                            if not kk.startswith("prof_")}
                                   for p in PASSINGS}} for n in have}, f, indent=1)
    out("wrote plots/L1_vs_N.png, plots/amp_phase_vs_t.png, plots/residual.png, "
        "plots/summary.txt, plots/results.json")
