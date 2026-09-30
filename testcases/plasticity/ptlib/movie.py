"""Movies of plasticity runs: particles coloured by a field, one mp4 per run (elastic plan rule 7:
a movie counts only after a frame has been looked at -- `pt movie` also writes a PNG of the last
frame next to the mp4 for exactly that)."""
import glob
import os
import shutil
import subprocess
from concurrent.futures import ProcessPoolExecutor

import numpy as np

FIELDS = {
    "gp":    ("plastic shear strain $\\gamma_p$", lambda d, m: d.get("plasticStrain", np.zeros(len(d["m"])))),
    "j2":    ("$\\sqrt{2 J_2}/Y$", lambda d, m: np.sqrt(d["Sxx"]**2 + d["Syy"]**2 + 2 * d["Sxy"]**2) / m["Y0"]),
    "sigma": ("$\\Sigma = P - S_{xx}$", lambda d, m: (d["P"] if "P" in d else m["A"] * (d["rho"] / m["rho0"] - 1))
              - d["Sxx"]),
    "rho":   ("$\\rho/\\rho_0$", lambda d, m: d["rho"] / m["rho0"]),
    "fab":   ("tensile-correction factor f", lambda d, m: d.get("fab", d.get("fabParticle", np.zeros(len(d["m"]))))),
    "u":     ("specific internal energy u", lambda d, m: d["u"]),
}


def _frame(job):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from .check import read
    f, out, field, mat, lims, title, ps = job
    d = read(f)
    label, fn = FIELDS[field]
    c = fn(d, mat)
    fig, ax = plt.subplots(figsize=(10, 5.6), dpi=110)
    sc = ax.scatter(d["x"][:, 0], d["x"][:, 1], c=c, s=ps, cmap="viridis", vmin=lims[0], vmax=lims[1],
                    linewidths=0)
    ax.set_xlim(lims[2], lims[3]); ax.set_ylim(lims[4], lims[5]); ax.set_aspect("equal")
    fig.colorbar(sc, ax=ax, label=label, shrink=0.8)
    ax.set_title("%s   t = %.4f" % (title, d["t"]))
    fig.savefig(out); plt.close(fig)


def render(S, rel, field, jobs=6):
    from .check import read, snap_files, mat_of
    test = rel.split("/")[0]
    fs = snap_files(S, rel)
    if not fs:
        return None
    mat = mat_of(S, test)
    label, fn = FIELDS[field]
    vals = np.concatenate([fn(read(f), mat) for f in fs[::max(1, len(fs) // 8)] + [fs[-1]]])
    lo, hi = np.percentile(vals, 1), np.percentile(vals, 99.5)
    if hi <= lo:
        hi = lo + 1e-12
    xs = np.concatenate([read(f, ["x"])["x"] for f in (fs[0], fs[-1])])
    pad = 0.05 * np.ptp(xs, axis=0).max()
    lims = (lo, hi, xs[:, 0].min() - pad, xs[:, 0].max() + pad, xs[:, 1].min() - pad, xs[:, 1].max() + pad)
    N = len(read(fs[0], ["m"])["m"])
    ps = max(0.2, min(12.0, 4e4 / N))
    tmp = S.p("movies", "_frames_" + rel.replace("/", "_") + "_" + field)
    shutil.rmtree(tmp, ignore_errors=True)
    os.makedirs(tmp)
    title = "%s  (%s)" % (rel, field)
    with ProcessPoolExecutor(jobs) as ex:
        list(ex.map(_frame, [(f, os.path.join(tmp, "%05d.png" % i), field, mat, lims, title, ps)
                             for i, f in enumerate(fs)]))
    name = rel.replace("/", "_") + "_" + field
    mp4 = S.p("movies", name + ".mp4")
    subprocess.run(["ffmpeg", "-y", "-loglevel", "error", "-framerate", "12", "-pattern_type", "glob", "-i",
                    tmp + "/*.png", "-vf", "pad=ceil(iw/2)*2:ceil(ih/2)*2", "-c:v", "libx264", "-pix_fmt",
                    "yuv420p", mp4], check=True)
    last = sorted(glob.glob(tmp + "/*.png"))[-1]
    shutil.copy(last, S.p("movies", name + "_last.png"))
    shutil.rmtree(tmp)
    return mp4


def _nearest(fs_times, t):
    return int(np.argmin(np.abs(np.asarray(fs_times) - t)))


def _encode(S, tmp, name, keep_mid=True):
    mp4 = S.p("movies", name + ".mp4")
    subprocess.run(["ffmpeg", "-y", "-loglevel", "error", "-framerate", "12", "-pattern_type", "glob", "-i",
                    tmp + "/*.png", "-vf", "pad=ceil(iw/2)*2:ceil(ih/2)*2", "-c:v", "libx264", "-pix_fmt",
                    "yuv420p", mp4], check=True)
    frames = sorted(glob.glob(tmp + "/*.png"))
    shutil.copy(frames[-1], S.p("movies", name + "_last.png"))
    if keep_mid:
        shutil.copy(frames[len(frames) // 2], S.p("movies", name + "_mid.png"))
    shutil.rmtree(tmp)
    return mp4


def _cmp_frame(job):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from .check import read
    out, panels, t0, mat, suptitle = job
    n = len(panels)
    ncol = 2 if n == 4 else n
    nrow = (n + ncol - 1) // ncol
    fig, axs = plt.subplots(nrow, ncol, figsize=(6.2 * ncol, 4.6 * nrow), dpi=100, squeeze=False)
    for ax, (f, field, label, lims, ps) in zip(axs.ravel(), panels):
        d = read(f)
        flabel, fn = FIELDS[field]
        sc = ax.scatter(d["x"][:, 0], d["x"][:, 1], c=fn(d, mat), s=ps, cmap="viridis", vmin=lims[0], vmax=lims[1],
                        linewidths=0)
        ax.set_xlim(lims[2], lims[3]); ax.set_ylim(lims[4], lims[5]); ax.set_aspect("equal")
        fig.colorbar(sc, ax=ax, label=flabel, shrink=0.85)
        ax.set_title("%s   t = %.3f" % (label, d["t"]), fontsize=10)
    for ax in axs.ravel()[n:]:
        ax.axis("off")
    fig.suptitle(suptitle, fontsize=11)
    fig.tight_layout()
    fig.savefig(out); plt.close(fig)


def render_compare(S, test, panels, name, jobs=6):
    """panels: [(rel, field, label)] -> one mp4, frames matched to the first panel's snapshot times; each panel
    has its own colour scale (1st-99.5th percentile over its run)."""
    from .check import read, snap_files, mat_of
    mat = mat_of(S, test)
    runs = []
    for rel, field, label in panels:
        fs = snap_files(S, rel)
        ts = [read(f, ["time"])["t"] for f in fs]
        fn = FIELDS[field][1]
        vals = np.concatenate([fn(read(f), mat) for f in fs[::max(1, len(fs) // 8)] + [fs[-1]]])
        lo, hi = np.percentile(vals, 1), np.percentile(vals, 99.5)
        hi = hi if hi > lo else lo + 1e-12
        runs.append((fs, ts, field, label, lo, hi))
    xs = np.concatenate([read(f, ["x"])["x"] for r in runs for f in (r[0][0], r[0][-1])])
    pad = 0.05 * np.ptp(xs, axis=0).max()
    box = (xs[:, 0].min() - pad, xs[:, 0].max() + pad, xs[:, 1].min() - pad, xs[:, 1].max() + pad)
    N = len(read(runs[0][0][0], ["m"])["m"])
    ps = max(0.3, min(14.0, 2.5e4 / N))
    tmp = S.p("movies", "_frames_" + name)
    shutil.rmtree(tmp, ignore_errors=True)
    os.makedirs(tmp)
    title = "%s  (commit %s)" % (test, S.meta["commit"][:7])
    jobs_ = []
    for i, t in enumerate(runs[0][1]):
        pan = [(r[0][_nearest(r[1], t)], r[2], r[3], (r[4], r[5]) + box, ps) for r in runs]
        jobs_.append((os.path.join(tmp, "%05d.png" % i), pan, t, mat, title))
    with ProcessPoolExecutor(jobs) as ex:
        list(ex.map(_cmp_frame, jobs_))
    return _encode(S, tmp, name)


def _prof_frame(job):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    out, t, curves, ref, xlim, ylims, title, vlabel = job
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(11, 7), dpi=100, sharex=True)
    for (X, sig, v, label, st) in curves:
        a1.plot(X, sig, st, lw=1.2, label=label)
        a2.plot(X, v, st, lw=1.2, label=label)
    if ref is not None:
        X, sig, v, label = ref
        a1.plot(X, sig, "k-", lw=0.9, alpha=0.7, label=label)
        a2.plot(X, v, "k-", lw=0.9, alpha=0.7, label=label)
    a1.set_ylabel("$\\Sigma = P - S_{xx}$ (rel. to undisturbed)"); a2.set_ylabel(vlabel)
    a2.set_xlabel("Lagrangian position X (initial x)")
    a1.set_xlim(*xlim); a1.set_ylim(*ylims[0]); a2.set_ylim(*ylims[1])
    a1.legend(loc="upper right", fontsize=9); a1.grid(alpha=0.3); a2.grid(alpha=0.3)
    a1.set_title("%s   t = %.3f" % (title, t))
    fig.tight_layout(); fig.savefig(out); plt.close(fig)


def render_profile(S, test, runs_, name, jobs=6):
    """P2/P3: y-averaged Sigma and v_x against the Lagrangian coordinate, every run as a line, plus the exact
    two-wave solution (P2, target frame) or the 1D reference (P3)."""
    from .check import read, snap_files, mat_of, profile, exact_plate, ref1d, pressure
    T = S.M.TESTS[test]
    mat = mat_of(S, test)
    kind, a = T["ic"]
    styles = ["C0-", "C3-", "C2-", "C1-"]
    series = []
    for k, (rel, label) in enumerate(runs_):
        fs = snap_files(S, rel)
        d0 = read(fs[0])
        X = d0["x"][:, 0]
        sel = X > 0 if kind == "slabs" else np.ones(len(X), bool)
        s0 = float(np.median((pressure(d0, mat) - d0["Sxx"])[(X > 1.0) & (X < 3.0)]))
        dd = float(rel.split("_d")[-1])
        fr = []
        for f in fs:
            q = read(f)
            v = q["v"][:, 0] + (a["up"] if kind == "slabs" else 0.0)
            Xc, sg = profile(X[sel], (pressure(q, mat) - q["Sxx"])[sel] - s0, X[sel].min(), X[sel].max(), dd)
            _, vv = profile(X[sel], v[sel], X[sel].min(), X[sel].max(), dd)
            fr.append((q["t"], Xc, sg, vv))
        series.append((fr, label, styles[k % 4]))
    times = [f[0] for f in series[0][0]]
    if kind == "slabs":
        ex = exact_plate(mat, a["up"], 1)
        xlim = (0, a["L"])
        def ref_at(t):
            X = np.linspace(0, a["L"], 2000)
            sig = np.where(X < ex["W_p"] * t, ex["sigma_b"], np.where(X < ex["U_e"] * t, ex["sigma_hel"], 0.0))
            v = np.where(X < ex["W_p"] * t, a["up"], np.where(X < ex["U_e"] * t, ex["u_hel"], 0.0))
            return X, sig, v, "exact (true 2D)"
        vlabel = "$v_x$ in the frame of the undisturbed target"
    else:
        rr = ref1d(mat, a["Lf"], a["Lt"], a["vf"], times, 4000)
        xlim = (-a["Lf"], a["Lt"])
        def ref_at(t):
            k = _nearest([r[0] for r in rr], t)
            return rr[k][1], rr[k][2], rr[k][3], "1D reference (4000 cells)"
        vlabel = "$v_x$ (lab)"
    allsig = np.concatenate([f[2] for s in series for f in s[0]])
    allv = np.concatenate([f[3] for s in series for f in s[0]])
    ylims = ((min(-0.01, allsig.min() * 1.1), allsig.max() * 1.15), (min(-0.01, allv.min() * 1.1), allv.max() * 1.15))
    tmp = S.p("movies", "_frames_" + name)
    shutil.rmtree(tmp, ignore_errors=True)
    os.makedirs(tmp)
    title = "%s  (commit %s)" % (test, S.meta["commit"][:7])
    jobs_ = []
    for i, t in enumerate(times):
        curves = []
        for fr, label, st in series:
            k = _nearest([f[0] for f in fr], t)
            curves.append((fr[k][1], fr[k][2], fr[k][3], label, st))
        jobs_.append((os.path.join(tmp, "%05d.png" % i), t, curves, ref_at(t), xlim, ylims, title, vlabel))
    with ProcessPoolExecutor(jobs) as ex:
        list(ex.map(_prof_frame, jobs_))
    return _encode(S, tmp, name)


def cmd_movie(args):
    from .core import find_study
    S = find_study(args.study)
    os.makedirs(S.p("movies"), exist_ok=True)
    variants = args.variants.split(",") if args.variants else None
    rungs = args.rungs.split(",") if args.rungs else None
    if args.compare or args.profile:
        for test in S.tests(args.tests):
            T = S.M.TESTS[test]
            for r in T["rungs"]:
                sel = [(v, rel) for v, rr, rel in S.runs(test, variants, [r])
                       if rr == r and os.path.exists(S.p("runs", rel, "exit_code"))]
                if len(sel) < 1:
                    continue
                tag = "%s_d%s" % (test, r)
                if args.profile:
                    out = render_profile(S, test, [(rel, v) for v, rel in sel], tag + "_profiles")
                else:
                    fields = args.field.split(",")
                    panels = [(rel, fields[min(i, len(fields) - 1)], v) for i, (v, rel) in enumerate(sel)]
                    out = render_compare(S, test, panels, tag + "_" + "-".join(dict.fromkeys(fields)))
                print("  ", out, "(inspect the _mid/_last.png frames before using it)")
        return
    for test in S.tests(args.tests):
        for v, r, rel in S.runs(test, variants, rungs):
            if os.path.exists(S.p("runs", rel, "exit_code")):
                out = render(S, rel, args.field)
                if out:
                    print("  ", out, "(inspect the _last.png frame before using it)")
