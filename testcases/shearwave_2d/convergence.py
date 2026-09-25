#!/usr/bin/env python3
"""Convergence analysis of the 2D shear-wave runs.

Structured like testcases/soundwave_2d_elastic/convergence.py so the two tables can be read
side by side, but for the transverse mode: at every passing t = T, 2T, 3T (T = lambda/c_s,
snapshots picked by their stored /time, never by file index -- dump-index caveat in
makeConfigs.py) the exact solution is the IC again,

  L1_vy     1/N sum |vy_i - c_s A sin(k(x_i - c_s t))|
  L1_Sxy    1/N sum |Sxy_i + rho_0 c_s^2 A sin(k(x_i - c_s t))|

and from a sine fit of the y-averaged profiles:

  a_vy, a_Sxy     amplitude in units of c_s A resp. rho_0 c_s^2 A
  phi             phase shift (+ = ahead of c_s)
  c_num/c_s - 1   measured shear-wave speed relative to sqrt(mu/rho_0)

Checks specific to this mode -- it is the one whose restoring force is entirely deviatoric,
and which to leading order carries no density perturbation:

  coupling    max |S_xy + rho c_s v_y| / (rho_0 c_s^2 A), the travelling-wave relation; the
              analogue of the "slaving" defect of the P-wave, and the quantity that exposed
              the O(dt) operator-splitting offset on 2026-09-18
  drho        max |rho/<rho> - 1| / A, the longitudinal contamination. A shear wave perturbs
              the density only at O(A^2) = 1e-12 relative, so anything at O(A) is spurious.
  L1_vx       the same contamination in the velocity
  S_norm      max |S_xx|, |S_yy| / (rho_0 c_s^2 A): zero to leading order, fed at O(A^2) by
              the Jaumann terms (R_xy S_xy with R_xy = d_x v_y), so ~1e-6 in these units is
              the prediction, not an error

Writes plots/L1_vs_N.png, plots/amp_phase_vs_t.png, plots/residual.png, plots/summary.txt
and plots/results.json.
"""

import glob
import json
import os

import h5py
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RESOLUTIONS = (64, 128, 256)
PASSINGS = (1, 2, 3)

SURFACE, INK, INK_2, MUTED, GRID, AXIS = "#fcfcfb", "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7"
RES_COLOR = {64: "#86b6ef", 128: "#2a78d6", 256: "#104281"}

plt.rcParams.update({
    "font.family": "sans-serif", "font.size": 9, "text.color": INK,
    "axes.labelcolor": INK_2, "axes.edgecolor": AXIS, "axes.linewidth": 0.8,
    "xtick.color": MUTED, "ytick.color": MUTED, "xtick.labelcolor": INK_2,
    "ytick.labelcolor": INK_2, "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
    "savefig.facecolor": SURFACE,
})


def stem(n):
    return "shear_N{}_corr".format(n)


def ic_name(n):
    return "shearwave_N{}_corr.h5".format(n)


def snapshots(outdir):
    files = glob.glob(os.path.join(outdir, "*.h5"))
    times = []
    for f in files:
        with h5py.File(f, "r") as d:
            times.append(float(np.ravel(d["time"])[0]))
    order = np.argsort(times)
    return [files[i] for i in order], [times[i] for i in order]


def meta(n):
    with h5py.File(ic_name(n), "r") as f:
        return {k: float(f.attrs[k]) for k in
                ("amplitude", "wavelength", "mu", "c_wave", "c_l", "period", "kernel_offset",
                 "rho_bar")}


def per_snapshot(path, md, n):
    A, k = md["amplitude"], 2 * np.pi / md["wavelength"]
    c_s, rho0 = md["c_wave"], md["rho_bar"]
    V_SCALE = c_s * A                      # velocity amplitude
    S_SCALE = rho0 * c_s * c_s * A         # stress amplitude, = mu A
    with h5py.File(path, "r") as d:
        t = float(np.ravel(d["time"])[0])
        x, vx, vy = d["x"][:, 0], d["v"][:, 0], d["v"][:, 1]
        rho, Sxx, Sxy, Syy = d["rho"][:], d["Sxx"][:], d["Sxy"][:], d["Syy"][:]
        cons = {key: float(np.ravel(d[key])[0])
                for key in ("totalMass", "energy", "xMomentum", "yMomentum")}
    phase = k * c_s * t
    wave = np.sin(k * x - phase)
    rho_m = rho.mean()

    out = dict(t=t, rho_mean=rho_m, **cons)
    out["L1_vy"] = np.mean(np.abs(vy - V_SCALE * wave))
    out["L1_Sxy"] = np.mean(np.abs(Sxy + S_SCALE * wave))
    # travelling-wave coupling S_xy = -rho c_s v_y
    out["coupling"] = np.max(np.abs(Sxy + rho * c_s * vy)) / S_SCALE
    # longitudinal contamination: must stay at O(A^2)
    out["drho"] = np.max(np.abs(rho / rho_m - 1.0)) / A
    out["L1_vx"] = np.mean(np.abs(vx)) / V_SCALE
    out["S_norm"] = max(np.max(np.abs(Sxx)), np.max(np.abs(Syy))) / S_SCALE

    # y-averaged column profiles + sine fits
    col = np.floor(x * n).astype(int) % n
    cnt = np.maximum(np.bincount(col, minlength=n), 1)
    xc = (np.arange(n) + 0.5) / n
    M = np.column_stack([np.sin(k * xc - phase), np.cos(k * xc - phase)])
    prof_vy = np.bincount(col, weights=vy, minlength=n) / cnt / V_SCALE
    prof_Sxy = -np.bincount(col, weights=Sxy, minlength=n) / cnt / S_SCALE
    for name, vals in (("vy", prof_vy), ("Sxy", prof_Sxy)):
        (p, q), *_ = np.linalg.lstsq(M, vals, rcond=None)
        out["a_" + name] = float(np.hypot(p, q))
        out["phi_" + name] = float(np.arctan2(-q, p))
    out["prof_x"] = xc
    out["prof_res"] = prof_vy - np.sin(k * xc - phase)
    return out


def analyse(n):
    if not os.path.isdir("output_" + stem(n)):
        return None
    md = meta(n)
    files, _ = snapshots("output_" + stem(n))
    if not files:
        return None
    series = [per_snapshot(f, md, n) for f in files]
    at = {}
    for p in PASSINGS:
        target = p * md["period"]
        i = int(np.argmin([abs(s["t"] - target) for s in series]))
        # a run still in progress has no snapshot at the later boundaries; including its
        # nearest one would put a number from the wrong time into the convergence table
        if abs(series[i]["t"] - target) > 1e-6 * max(1.0, target):
            print("  skipping n = {}: no snapshot at {} T (nearest t = {:.6f}, target {:.6f})"
                  " -- run still in progress?".format(n, p, series[i]["t"], target))
            return None
        at[p] = series[i]
    return dict(meta=md, series=series, at=at, t0=series[0])


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
    T, c_s, c_l, mu, A = md["period"], md["c_wave"], md["c_l"], md["mu"], md["amplitude"]
    k = 2 * np.pi / md["wavelength"]

    lines = []

    def out(s=""):
        print(s)
        lines.append(s)

    out("2D travelling SHEAR wave, Tillotson (K = 1), mu = {:.9g} (Poisson 0.3)".format(mu))
    out("c_s = sqrt(mu/rho0) = {:.9f} (c_l = {:.9f}), period T = {:.9f}, A = {:g}".format(
        c_s, c_l, T, A))
    out("masses corrected for the kernel-sum offset {:+.6e}".format(md["kernel_offset"]))
    out("restoring force 100 % deviatoric; rho uniform at t = 0")
    out()

    for p in PASSINGS:
        out("--- t = {} T = {:.6f} ({} passing{}) ---".format(p, p * T, p, "" if p == 1 else "s"))
        out("  {:>5s} {:>11s} {:>11s} {:>9s} {:>9s} {:>12s} {:>11s} {:>10s}".format(
            "n", "L1_vy", "L1_Sxy", "a_vy", "a_Sxy", "c_num/c_s-1", "coupling", "drho [A]"))
        for n in have:
            s = res[n]["at"][p]
            cnum = s["phi_vy"] / (k * s["t"]) / c_s
            out("  {:>5d} {:11.4e} {:11.4e} {:9.6f} {:9.6f} {:+12.3e} {:11.3e} {:10.3e}".format(
                n, s["L1_vy"], s["L1_Sxy"], s["a_vy"], s["a_Sxy"], cnum,
                s["coupling"], s["drho"]))
        if len(have) > 1:
            for key in ("L1_vy", "L1_Sxy", "coupling"):
                o = order([res[n]["at"][p][key] for n in have], have)
                out("  order {:<9s} {}".format(key, "  ".join("{:.2f}".format(v) for v in o)))
        out()

    out("--- conservation and mode purity (over the whole run) ---")
    out("  {:>5s} {:>12s} {:>12s} {:>12s} {:>12s} {:>11s} {:>11s}".format(
        "n", "dM/M", "dE/E", "dPy", "max drho [A]", "max L1_vx", "max S_norm"))
    for n in have:
        s0, ser = res[n]["t0"], res[n]["series"]
        dM = max(abs(s["totalMass"] / s0["totalMass"] - 1) for s in ser)
        dE = max(abs(s["energy"] / s0["energy"] - 1) for s in ser)
        dP = max(abs(s["yMomentum"] - s0["yMomentum"]) for s in ser)
        out("  {:>5d} {:12.3e} {:12.3e} {:12.3e} {:12.3e} {:11.3e} {:11.3e}".format(
            n, dM, dE, dP, max(s["drho"] for s in ser), max(s["L1_vx"] for s in ser),
            max(s["S_norm"] for s in ser)))
    out()

    # ---------------- plots ----------------
    fig, axes = plt.subplots(1, 3, figsize=(10.5, 3.6))
    for ax, key, lab in zip(axes, ("L1_vy", "L1_Sxy", "coupling"),
                            (r"$L_1(v_y)$", r"$L_1(S_{xy})$",
                             r"coupling defect $|S_{xy}+\rho c_s v_y|$")):
        for p, col, mk in zip(PASSINGS, ("#86b6ef", "#2a78d6", "#104281"), ("o", "s", "^")):
            y = [res[n]["at"][p][key] for n in have]
            ax.plot(have, y, "-", marker=mk, ms=5, lw=2, color=col, label="{} T".format(p))
        ref = 0.45 * res[have[0]]["at"][PASSINGS[-1]][key]
        guide = [ref * (have[0] / n) ** 2 for n in have]
        ax.plot(have, guide, color=AXIS, lw=1.2, zorder=1)
        ax.text(have[len(have) // 2], guide[len(have) // 2], "slope −2\n", color=MUTED,
                fontsize=8, va="center", ha="center", rotation=-38, rotation_mode="anchor")
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
        ax.set_xticks(have); ax.set_xticklabels([str(n) for n in have])
        ax.set_xlabel("particles per wavelength $n$")
        ax.set_title(lab, color=INK_2, fontsize=9, loc="left")
        ax.grid(True, which="both", color=GRID, lw=0.6); ax.set_axisbelow(True)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
    axes[0].legend(frameon=False, fontsize=8, labelcolor=INK_2, title="after", title_fontsize=8)
    fig.suptitle("Shear wave: convergence after 1, 2 and 3 passings",
                 x=0.02, ha="left", fontsize=11, color=INK)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig("plots/L1_vs_N.png", dpi=170)
    plt.close(fig)

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 5.4), sharex=True,
                                   gridspec_kw={"hspace": 0.12})
    for n in have:
        ser = res[n]["series"]
        tt = np.array([s["t"] for s in ser]) / T
        ax1.plot(tt, [s["a_vy"] for s in ser], lw=1.8, color=RES_COLOR[n], label="n = {}".format(n))
        ax1.plot(tt, [s["a_Sxy"] for s in ser], lw=1.0, ls="--", color=RES_COLOR[n])
        ax2.plot(tt, [s["phi_vy"] for s in ser], lw=1.8, color=RES_COLOR[n])
    ax1.set_ylabel("amplitude  [$c_s A$, $\\mu A$]")
    ax1.text(0.99, 0.06, "solid: $v_y$    dashed: $S_{xy}$", transform=ax1.transAxes,
             fontsize=8, color=MUTED, ha="right")
    ax2.set_ylabel("phase shift  [rad]\n(+ = ahead of $c_s$)")
    ax2.set_xlabel("passings  $t/T$")
    for ax in (ax1, ax2):
        ax.grid(True, color=GRID, lw=0.6); ax.set_axisbelow(True)
        ax.axhline(1.0 if ax is ax1 else 0.0, color=AXIS, lw=0.8)
        for p in PASSINGS:
            ax.axvline(p, color=AXIS, lw=0.8, ls=":")
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
    ax1.legend(frameon=False, fontsize=8, labelcolor=INK_2, loc="lower left")
    fig.suptitle("Shear wave: amplitude and phase over 3 passings",
                 x=0.02, ha="left", fontsize=11, color=INK)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig("plots/amp_phase_vs_t.png", dpi=170)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8, 3.4))
    for n in have:
        s = res[n]["at"][PASSINGS[-1]]
        ax.plot(s["prof_x"], s["prof_res"], lw=1.4, color=RES_COLOR[n], label="n = {}".format(n))
    ax.set_xlim(0, 1)
    ax.set_xlabel("x"); ax.set_ylabel("numerical − analytic  [$c_s A$]")
    ax.axhline(0, color=AXIS, lw=0.8)
    ax.grid(True, color=GRID, lw=0.6); ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.legend(frameon=False, fontsize=8, labelcolor=INK_2)
    fig.suptitle("$v_y$ residual after 3 passings", x=0.02, ha="left", fontsize=11, color=INK)
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
