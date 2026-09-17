#!/usr/bin/env python3
"""Convergence analysis of the 2D soundwave runs (Hopkins 2015, sec. 4.2.1).

Per run, at t = 1 (one period; the solution equals the IC):

  L1_rho_nom   1/N sum |rho_i - (1 + A sin(k(x_i - t)))|         Hopkins eq. 31, nominal rho_bar = 1
  L1_rho       1/N sum |rho_i - <rho>(1 + A sin(k(x_i - t)))|    same, kernel-sum offset removed
  L1_rho_self  1/N sum |rho_i(t=1) - rho_i(t=0)|                 against the run's own t = 0 state
  L1_vx        1/N sum |vx_i - A sin(k(x_i - t))|

and from a sine fit of the y-averaged profiles, y ~ a sin(k(x - t) - phi):
  a_rho(0), a_rho(1), a_vx(1)   amplitude in units of A
  phi(1)                        phase shift (+ = ahead of c = 1)
  c_num - 1 = phi / (k t)       measured wave speed
  c_exp - 1                     speed expected from the realised mean state
                                (ideal: sqrt(gamma(gamma-1)<u>), Murnaghan: sqrt(<rho>), Tillotson: 1)
  dc_disp = c_num - c_exp       numerical dispersion

Also: conservation (mass, energy, x-momentum) and the kernel-sum density offset
from an ideal-lattice Wendland C2 sum at the solver's h/dx, compared with <rho> - 1.

GIZMO runs on the same ICs (output_<eos>_N<n>_gizmo, ideal gas and Tillotson, from
~/Documents/GIZMO/soundwave_2d) are added to the table, the orders and to
plots/convergence.png and plots/L1_vs_N.png (dashed, diamonds) where present.

Writes plots/convergence.png, plots/L1_vs_N.png, plots/amp_phase_vs_t.png, plots/residual_t1.png,
plots/convergence_summary.txt, plots/convergence_results.json.
"""

import json
import os
import numpy as np
import h5py
import matplotlib
import matplotlib.ticker
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from makeMovie import snapshots, column_profile, sine_fit

EOS_NAMES = ("ideal", "murnaghan", "tillotson")
EOS_LABEL = {"ideal": "ideal gas", "murnaghan": "Murnaghan", "tillotson": "Tillotson"}
RUNS = [(e, n, 32, False) for e in EOS_NAMES for n in (64, 128, 256)] + \
       [(e, n, 64, False) for e in EOS_NAMES for n in (64, 128)] + \
       [("murnaghan", n, 32, True) for n in (64, 128, 256)]
# figure columns: (eos, masses corrected)
COLUMNS = [("ideal", False), ("murnaghan", False), ("murnaghan", True), ("tillotson", False)]
COL_LABEL = {("ideal", False): "ideal gas", ("murnaghan", False): "Murnaghan",
             ("murnaghan", True): "Murnaghan, masses corrected", ("tillotson", False): "Tillotson"}
GAMMA = 5.0 / 3.0

SURFACE, INK, INK_2, MUTED, GRID, AXIS = "#fcfcfb", "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7"
EOS_COLOR = {"ideal": "#2a78d6", "murnaghan": "#eb6834", "tillotson": "#1baf7a"}   # categorical slots 1-3
RES_COLOR = {64: "#86b6ef", 128: "#2a78d6", 256: "#104281"}                     # ordinal blue 250/450/650

plt.rcParams.update({
    "font.family": "sans-serif", "font.size": 9, "text.color": INK,
    "axes.labelcolor": INK_2, "axes.edgecolor": AXIS, "axes.linewidth": 0.8,
    "xtick.color": MUTED, "ytick.color": MUTED, "xtick.labelcolor": INK_2, "ytick.labelcolor": INK_2,
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
})


def stem(eos, n, nnn, corr=False, suffix=""):
    return "{}_N{}".format(eos, n) + ("" if nnn == 32 else "_nnn{}".format(nnn)) + ("_corr" if corr else "") + suffix


def wendland_c2_lattice_offset(h_over_dx, R=10):
    i = np.arange(-R, R + 1)
    X, Y = np.meshgrid(i, i)
    u = np.hypot(X, Y) / h_over_dx
    t = np.clip(1.0 - u, 0.0, None)
    return (7.0 / np.pi / h_over_dx**2 * t**4 * (1.0 + 4.0 * u)).sum() - 1.0


def analyse(eos, n, nnn, corr=False, suffix=""):
    s = stem(eos, n, nnn, corr, suffix)
    with h5py.File("soundwave_{}_N{}.h5".format(eos, n), "r") as f:
        A, lam = float(f.attrs["amplitude"]), float(f.attrs["wavelength"])
    k = 2.0 * np.pi / lam
    files, times = snapshots("output_" + s)
    if abs(times[-1] - 1.0) > 1e-9:
        raise SystemExit("{}: last snapshot at t = {}".format(s, times[-1]))

    def load(path):
        with h5py.File(path, "r") as d:
            return {key: d[key][()] for key in ("x", "v", "rho", "u", "P", "sml", "totalMass",
                                                "energy", "xMomentum", "yMomentum", "time")}
    d0, d1 = load(files[0]), load(files[-1])
    x, rho, vx, t = d1["x"][:, 0], d1["rho"], d1["v"][:, 0], float(d1["time"][0])
    wave = np.sin(k * (x - t))
    mean = rho.mean()

    r = {"run": s, "eos": eos, "n": n, "nnn": nnn, "corr": corr, "N": int(rho.size), "t": t,
         "rho_mean_m1_t0": float(d0["rho"].mean() - 1.0), "rho_mean_m1_t1": float(mean - 1.0),
         "h_over_dx": float(d0["sml"].mean() * n),
         "L1_rho_nom": float(np.mean(np.abs(rho - (1.0 + A * wave)))),
         "L1_rho": float(np.mean(np.abs(rho - mean * (1.0 + A * wave)))),
         "L1_rho_self": float(np.mean(np.abs(rho - d0["rho"]))),
         "L1_vx": float(np.mean(np.abs(vx - A * wave)))}
    r["lattice_offset"] = float(wendland_c2_lattice_offset(r["h_over_dx"]))

    # amplitude / phase time series of the y-averaged density and velocity
    ts, a_rho, phi_rho, a_vx, phi_vx = [], [], [], [], []
    for f, tt in zip(files, times):
        xc, prof, m = column_profile(f, n)
        a, p = sine_fit(xc, (prof - m) / (m * A), k, tt)
        with h5py.File(f, "r") as d:
            xv, vv = d["x"][:, 0], d["v"][:, 0]
        col = np.floor(xv * n).astype(int) % n
        vprof = np.bincount(col, weights=vv, minlength=n) / np.bincount(col, minlength=n)
        av, pv = sine_fit(xc, vprof / A, k, tt)
        ts.append(tt); a_rho.append(a); phi_rho.append(p); a_vx.append(av); phi_vx.append(pv)
    r["series"] = {"t": ts, "a_rho": a_rho, "phi_rho": phi_rho, "a_vx": a_vx, "phi_vx": phi_vx}
    r["a_rho_t0"], r["a_rho_t1"], r["a_vx_t1"] = a_rho[0], a_rho[-1], a_vx[-1]
    r["phi_rho_t1"], r["phi_vx_t1"] = phi_rho[-1], phi_vx[-1]
    r["c_num_m1"] = phi_rho[-1] / (k * t)

    if eos == "ideal":
        c_exp = np.sqrt(GAMMA * (GAMMA - 1.0) * d1["u"].mean())
    elif eos == "murnaghan":
        c_exp = np.sqrt(mean)          # K0 = rho0 = 1, n = 2
    else:
        c_exp = 1.0                    # A = rho0 = 1, a = b = 0
    r["c_exp_m1"] = float(c_exp - 1.0)
    r["dc_disp"] = r["c_num_m1"] - r["c_exp_m1"]

    for key, name in (("totalMass", "dM"), ("energy", "dE")):
        r[name] = float((d1[key][0] - d0[key][0]) / d0[key][0])
    r["dpx"] = float(d1["xMomentum"][0] - d0["xMomentum"][0])
    r["dpy"] = float(d1["yMomentum"][0] - d0["yMomentum"][0])
    # non-isentropic content: entropy proxy K = P / rho^gamma (ideal gas), u (solids;
    # P does not depend on u there). Amplitude of the y-averaged, mean-removed
    # perturbation at t = 0 and t = 1, in units of A (sin and cos parts, so any phase).
    for tag, d in (("t0", d0), ("t1", d1)):
        xd = d["x"][:, 0]
        if eos == "ideal":
            q = d["P"] / d["rho"] ** GAMMA
            q = q / q.mean() - 1.0
        else:
            q = (d["u"] - d["u"].mean()) / d["u"].mean()
        col = np.floor(xd * n).astype(int) % n
        qc = np.bincount(col, weights=q, minlength=n) / np.bincount(col, minlength=n)
        xc_ = (np.arange(n) + 0.5) / n
        M = np.column_stack([np.sin(k * xc_), np.cos(k * xc_)])
        (p1, p2), *_ = np.linalg.lstsq(M, qc / A, rcond=None)
        r["nonisen_amp_" + tag] = float(np.hypot(p1, p2))
    xc, prof, m = column_profile(files[-1], n)
    r["residual"] = {"xc": xc.tolist(), "res": ((prof - m) / (m * A) - np.sin(k * (xc - t))).tolist()}
    return r


def order(v1, v2):
    return np.log2(v1 / v2) if v1 > 0 and v2 > 0 else float("nan")


def style_axes(ax, logx=True, logy=True):
    if logx:
        ax.set_xscale("log", base=2)
    if logy:
        ax.set_yscale("log")
    ax.grid(True, which="major", color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.tick_params(length=2.5, width=0.6)


if __name__ == "__main__":
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    os.makedirs("plots", exist_ok=True)
    res = {stem(*run): analyse(*run) for run in RUNS}
    # GIZMO on the same uncorrected ICs, N_NN = 32, where converted (~/Documents/GIZMO/soundwave_2d)
    GZ = {eos: [n for n in (64, 128, 256) if os.path.isdir("output_{}_N{}_gizmo".format(eos, n))]
          for eos in ("ideal", "tillotson")}
    res_gz = {stem(eos, n, 32, False, "_gizmo"): analyse(eos, n, 32, False, "_gizmo")
              for eos, gns in GZ.items() for n in gns}

    # ---- table ---------------------------------------------------------------
    lines = []
    hdr = ("{:22s} {:>6s} {:>10s} {:>10s} {:>10s} {:>10s} {:>10s} {:>8s} {:>8s} {:>8s} {:>10s} {:>10s} {:>10s} {:>9s} {:>9s}"
           .format("run", "N", "<rho>-1", "L1_nom", "L1_rho", "L1_self", "L1_vx", "a_rho0", "a_rho1",
                   "a_vx1", "phi1", "c_exp-1", "dc_disp", "dE/E", "dM/M"))
    lines.append(hdr)
    for s, r in {**res, **res_gz}.items():
        lines.append("{:22s} {:6d} {:+10.3e} {:10.3e} {:10.3e} {:10.3e} {:10.3e} {:8.5f} {:8.5f} {:8.5f} {:+10.3e} {:+10.3e} {:+10.3e} {:+9.1e} {:+9.1e}"
                     .format(s, r["N"], r["rho_mean_m1_t1"], r["L1_rho_nom"], r["L1_rho"], r["L1_rho_self"],
                             r["L1_vx"], r["a_rho_t0"], r["a_rho_t1"], r["a_vx_t1"], r["phi_rho_t1"],
                             r["c_exp_m1"], r["dc_disp"], r["dE"], r["dM"]))
    lines.append("")
    lines.append("Non-isentropic perturbation amplitude [A]: ideal dK/K (K = P/rho^gamma), solids du/<u>; t=0 -> t=1")
    for s_, r in {**res, **res_gz}.items():
        lines.append("  {:22s} {:.3e} -> {:.3e}".format(s_, r["nonisen_amp_t0"], r["nonisen_amp_t1"]))
    lines.append("")
    lines.append("Convergence orders p = log2(err(n/2)/err(n)), N_NN = 32 [64->128, 128->256]; N_NN = 64 [64->128]")
    metrics = [("L1_rho", "L1_rho"), ("L1_rho_self", "L1_self"), ("L1_vx", "L1_vx"),
               ("amp_def", "1-a_rho1"), ("amp_def_vx", "1-a_vx1"), ("dc_disp", "|dc_disp|")]

    def metric(r, key):
        if key == "amp_def":
            return 1.0 - r["a_rho_t1"]
        if key == "amp_def_vx":
            return 1.0 - r["a_vx_t1"]
        if key == "dc_disp":
            return abs(r["dc_disp"])
        return r[key]
    for eos in EOS_NAMES:
        parts = []
        for key, label in metrics:
            e = [metric(res[stem(eos, n, 32)], key) for n in (64, 128, 256)]
            f = [metric(res[stem(eos, n, 64)], key) for n in (64, 128)]
            parts.append("{} {:.2f},{:.2f} | {:.2f}".format(label, order(e[0], e[1]), order(e[1], e[2]),
                                                             order(f[0], f[1])))
        lines.append("  {:10s} ".format(eos) + ";  ".join(parts))
    parts = []
    for key, label in metrics:
        e = [metric(res[stem("murnaghan", n, 32, True)], key) for n in (64, 128, 256)]
        parts.append("{} {:.2f},{:.2f}".format(label, order(e[0], e[1]), order(e[1], e[2])))
    lines.append("  {:10s} ".format("murn_corr") + ";  ".join(parts))
    for eos, gns in GZ.items():
        if len(gns) < 2:
            continue
        parts = []
        for key, label in metrics:
            e = [metric(res_gz[stem(eos, n, 32, False, "_gizmo")], key) for n in gns]
            parts.append("{} ".format(label) + ",".join("{:.2f}".format(order(abs(a), abs(b))) for a, b in zip(e, e[1:])))
        lines.append("  {:10s} ".format(eos + "_gz") + ";  ".join(parts) + "   (GIZMO, n = {}, abs values)".format(gns))
    lines.append("")
    lines.append("Kernel-sum density offset: ideal-lattice Wendland C2 sum at the solver's h/dx vs measured <rho>-1 (t=0)")
    for nnn in (32, 64):
        r = res[stem("ideal", 64, nnn)]
        lines.append("  N_NN = {}: h/dx = {:.5f}, lattice sum {:+.7e}, measured {:+.7e}".format(
            nnn, r["h_over_dx"], r["lattice_offset"], r["rho_mean_m1_t0"]))
    lines.append("  c_exp - 1 for Murnaghan (sqrt(<rho>) - 1): N_NN = 32 {:+.4e}, N_NN = 64 {:+.4e}".format(
        res["murnaghan_N64"]["c_exp_m1"], res["murnaghan_N64_nnn64"]["c_exp_m1"]))
    text = "\n".join(lines)
    with open("plots/convergence_summary.txt", "w") as f:
        f.write(text + "\n")
    print(text)
    with open("plots/convergence_results.json", "w") as f:
        json.dump(res, f, indent=1)

    # ---- figure 1: convergence -------------------------------------------------
    panels = [("L1_rho", r"$L_1(\rho)$, offset removed"), ("L1_vx", r"$L_1(v_x)$"),
              ("amp_def", r"amplitude deficit $1 - a_\rho(t{=}1)$"),
              ("dc_disp", r"dispersion $|c_{num} - c_{exp}|$")]
    fig, axes = plt.subplots(1, 4, figsize=(13, 3.6), gridspec_kw={"wspace": 0.36})
    ns = np.array([64, 128, 256])
    for ax, (key, title) in zip(axes, panels):
        style_axes(ax)
        allv = []
        for eos in EOS_NAMES:
            c = EOS_COLOR[eos]
            v32 = [metric(res[stem(eos, n, 32)], key) for n in ns]
            v64 = [metric(res[stem(eos, n, 64)], key) for n in ns[:2]]
            ax.plot(ns, v32, "-o", color=c, lw=2, ms=7, mec=SURFACE, mew=1.5, zorder=3)
            ax.plot(ns[:2], v64, "-o", color=c, lw=1.2, ms=7, mfc=SURFACE, mec=c, mew=1.5, zorder=3)
            allv += v32 + v64
        vc = [metric(res[stem("murnaghan", n, 32, True)], key) for n in ns]
        ax.plot(ns, vc, "-s", color=EOS_COLOR["murnaghan"], lw=2, ms=7, mec=SURFACE, mew=1.5, zorder=4)
        allv += vc
        for eos, gns in GZ.items():
            if gns:
                vg = np.array([metric(res_gz[stem(eos, n, 32, False, "_gizmo")], key) for n in gns])
                # amplitude growth (a > 1) gives a negative deficit: plot |1 - a|, hollow markers
                ax.plot(gns, np.abs(vg), "--", color=EOS_COLOR[eos], lw=1.4, zorder=5)
                for n_, v_ in zip(gns, vg):
                    ax.plot(n_, abs(v_), "D", ms=6, color=EOS_COLOR[eos], mew=1.2, zorder=6,
                            mfc=EOS_COLOR[eos] if v_ >= 0 else SURFACE, mec=SURFACE if v_ >= 0 else EOS_COLOR[eos])
                if key == "amp_def" and (vg < 0).any():
                    ax.text(0.97, 0.97, "hollow ◇: growth, |1 − a| plotted", transform=ax.transAxes, ha="right", va="top",
                            fontsize=7.5, color=INK_2)
        ref0 = np.exp(np.mean(np.log([v for v in allv if v > 0])))
        ax.plot(ns, ref0 * 1.8 * (ns / 128.0) ** -2, color=MUTED, lw=1.0, zorder=2)
        ax.text(256, ref0 * 1.8 * 0.25, r"  $\propto n^{-2}$", color=MUTED, fontsize=8, va="center")
        ax.set_xticks(ns)
        ax.set_xticklabels([str(n) for n in ns])
        ax.set_xlim(52, 380)
        ax.set_xlabel("particles per wavelength n")
        ax.set_title(title, fontsize=9.5, color=INK, loc="left")
    handles = [Line2D([], [], color=EOS_COLOR[e], lw=2, marker="o", ms=7, mec=SURFACE, label=EOS_LABEL[e])
               for e in EOS_NAMES]
    handles += [Line2D([], [], color=EOS_COLOR["murnaghan"], lw=2, marker="s", ms=7, mec=SURFACE,
                       label="Murnaghan, masses corrected"),
                Line2D([], [], color=INK_2, lw=2, marker="o", ms=7, mec=SURFACE, label=r"$N_{NN} = 32$"),
                Line2D([], [], color=INK_2, lw=1.2, marker="o", ms=7, mfc=SURFACE, mec=INK_2, label=r"$N_{NN} = 64$")]
    if any(GZ.values()):
        handles.append(Line2D([], [], color=INK_2, lw=1.4, ls="--", marker="D", ms=6, mec=SURFACE,
                              label=r"GIZMO, $N_{NN} = 32$"))
    fig.legend(handles=handles, loc="upper center", ncol=len(handles), frameon=False, fontsize=8.5,
               bbox_to_anchor=(0.5, 1.07), labelcolor=INK_2)
    fig.savefig("plots/convergence.png", dpi=200, bbox_inches="tight")
    plt.close(fig)

    # ---- figure 0: Hopkins (2015) Fig. 2 bottom analogue, L1 vs N ---------------
    # x: particles per wavelength n, the 1D-equivalent of Hopkins' N (second order
    # stays L1 ~ n^-2); top axis: total 2D particle count N = n * n/4.
    fig, ax = plt.subplots(figsize=(5.6, 4.4))
    style_axes(ax)
    for eos in EOS_NAMES:
        c = EOS_COLOR[eos]
        ax.plot(ns, [res[stem(eos, n, 32)]["L1_rho"] for n in ns], "-o", color=c, lw=2, ms=8,
                mec=SURFACE, mew=1.5, zorder=3)
        ax.plot(ns[:2], [res[stem(eos, n, 64)]["L1_rho"] for n in ns[:2]], "-o", color=c, lw=1.2,
                ms=8, mfc=SURFACE, mec=c, mew=1.5, zorder=3)
    # corrected Murnaghan coincides with Tillotson: draw it underneath, wider, so both stay visible
    vc = [res[stem("murnaghan", n, 32, True)]["L1_rho"] for n in ns]
    ax.plot(ns, vc, "-s", color=EOS_COLOR["murnaghan"], lw=5, ms=13, mec=SURFACE, mew=1.5, zorder=2.5)
    ax.annotate("Murnaghan corrected\n≈ Tillotson", xy=(100, np.exp(np.interp(np.log(100), np.log(ns), np.log(vc)))),
                xytext=(66, vc[1] * 0.72), fontsize=8, color=INK_2, ha="left", va="top",
                arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.8))
    # GIZMO on the same ICs (if converted: ~/Documents/GIZMO/soundwave_2d)
    for eos, gns in GZ.items():
        if gns:
            ax.plot(gns, [res_gz[stem(eos, n, 32, False, "_gizmo")]["L1_rho"] for n in gns], "--D",
                    color=EOS_COLOR[eos], lw=1.4, ms=7, mec=SURFACE, mew=1.2, zorder=4)
    nref = np.array([48, 340])
    l1ref = res[stem("tillotson", 128, 32)]["L1_rho"] * 1.6
    ax.plot(nref, l1ref * (nref / 128.0) ** -2, ":", color=MUTED, lw=1.4, zorder=2)
    ax.text(200, l1ref * (200 / 128.0) ** -2 * 1.25, r"$\propto N^{-2}$", color=MUTED, fontsize=9,
            ha="left", va="bottom")
    ax.set_xticks(ns)
    ax.set_xticklabels([str(n) for n in ns])
    ax.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
    ax.set_xlim(48, 340)
    ax.set_xlabel("N (particles per wavelength)")
    ax.set_ylabel(r"$L_1$ error norm, $\rho$")
    top = ax.secondary_xaxis("top")
    top.set_xticks(ns)
    top.set_xticklabels([str(n * n // 4) for n in ns])
    top.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
    top.set_xlabel("total particles in 2D", color=INK_2)
    top.tick_params(length=2.5, width=0.6, colors=MUTED, labelcolor=INK_2)
    handles = [Line2D([], [], color=EOS_COLOR[e], lw=2, marker="o", ms=8, mec=SURFACE, label=EOS_LABEL[e])
               for e in EOS_NAMES]
    handles += [Line2D([], [], color=EOS_COLOR["murnaghan"], lw=4, marker="s", ms=10, mec=SURFACE,
                       label="Murnaghan, corrected"),
                Line2D([], [], color=INK_2, lw=2, marker="o", ms=8, mec=SURFACE, label=r"$N_{NN} = 32$"),
                Line2D([], [], color=INK_2, lw=1.2, marker="o", ms=8, mfc=SURFACE, mec=INK_2, label=r"$N_{NN} = 64$"),
                Line2D([], [], color=MUTED, lw=1.4, ls=":", label="second order")]
    if any(GZ.values()):
        handles.insert(4, Line2D([], [], color=INK_2, lw=1.4, ls="--", marker="D", ms=7, mec=SURFACE,
                                 label="GIZMO (same ICs)"))
    ax.legend(handles=handles, loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False, fontsize=8.5,
              labelcolor=INK_2, ncol=1)
    ax.set_title("2D soundwave, one period (MFM; demonstrator xSplit HLLC, GIZMO)", fontsize=10, color=INK, loc="left", pad=26)
    fig.savefig("plots/L1_vs_N.png", dpi=200, bbox_inches="tight")
    plt.close(fig)

    # ---- figure 2: amplitude and phase vs time ----------------------------------
    fig, axes = plt.subplots(2, 4, figsize=(15, 5.6), sharex=True, gridspec_kw={"hspace": 0.18, "wspace": 0.25})
    for col, (eos, corr) in enumerate(COLUMNS):
        for nnn in ((32,) if corr else (32, 64)):
            for n in ((64, 128, 256) if nnn == 32 else (64, 128)):
                r = res[stem(eos, n, nnn, corr)]
                sr = r["series"]
                kw = dict(color=RES_COLOR[n], lw=1.8 if nnn == 32 else 1.1, zorder=3)
                mk = {} if nnn == 32 else dict(marker="o", markevery=10, ms=5, mfc=SURFACE, mec=RES_COLOR[n])
                axes[0, col].plot(sr["t"], 1.0 - np.array(sr["a_rho"]), **kw, **mk)
                axes[1, col].plot(sr["t"], sr["phi_rho"], **kw, **mk)
        axes[0, col].set_title(COL_LABEL[(eos, corr)], fontsize=10, color=INK, loc="left")
        for row in (0, 1):
            style_axes(axes[row, col], logx=False, logy=False)
            axes[row, col].axhline(0, color=AXIS, lw=0.8, zorder=1)
        axes[1, col].set_xlabel("t")
        axes[1, col].set_xlim(0, 1)
    axes[0, 0].set_ylabel(r"$1 - a_\rho(t)$  [A]")
    axes[1, 0].set_ylabel(r"phase shift $\varphi(t)$ [rad], + = ahead")
    handles = [Line2D([], [], color=RES_COLOR[n], lw=1.8, label="n = {}".format(n)) for n in (64, 128, 256)]
    handles += [Line2D([], [], color=INK_2, lw=1.8, label=r"$N_{NN} = 32$"),
                Line2D([], [], color=INK_2, lw=1.1, marker="o", ms=5, mfc=SURFACE, mec=INK_2, label=r"$N_{NN} = 64$")]
    fig.legend(handles=handles, loc="upper center", ncol=5, frameon=False, fontsize=8.5,
               bbox_to_anchor=(0.5, 1.0), labelcolor=INK_2)
    fig.savefig("plots/amp_phase_vs_t.png", dpi=200, bbox_inches="tight")
    plt.close(fig)

    # ---- figure 3: density residual at t = 1 ------------------------------------
    fig, axes = plt.subplots(2, 4, figsize=(15, 5.2), sharex=True, gridspec_kw={"hspace": 0.2, "wspace": 0.25})
    for col, (eos, corr) in enumerate(COLUMNS):
        for row, nnn in enumerate((32, 64)):
            ax = axes[row, col]
            if corr and nnn == 64:
                ax.axis("off")
                continue
            style_axes(ax, logx=False, logy=False)
            ax.axhline(0, color=AXIS, lw=0.8, zorder=1)
            for n in ((64, 128, 256) if nnn == 32 else (64, 128)):
                rr = res[stem(eos, n, nnn, corr)]["residual"]
                ax.plot(rr["xc"], rr["res"], color=RES_COLOR[n], lw=1.4, zorder=3, label="n = {}".format(n))
            ax.set_xlim(0, 1)
            if row == 0:
                ax.set_title(COL_LABEL[(eos, corr)], fontsize=10, color=INK, loc="left")
            if col == 0:
                ax.set_ylabel(r"$N_{{NN}} = {}$" "\n" r"$\delta\rho/A$ − analytic".format(nnn))
            if row == 1 or corr:
                ax.set_xlabel("x")
    axes[0, 0].legend(loc="lower left", frameon=False, fontsize=8, labelcolor=INK_2)
    fig.savefig("plots/residual_t1.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print("Wrote plots/convergence.png, plots/amp_phase_vs_t.png, plots/residual_t1.png")
