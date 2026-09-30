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


def cmd_movie(args):
    from .core import find_study
    S = find_study(args.study)
    os.makedirs(S.p("movies"), exist_ok=True)
    variants = args.variants.split(",") if args.variants else None
    rungs = args.rungs.split(",") if args.rungs else None
    for test in S.tests(args.tests):
        for v, r, rel in S.runs(test, variants, rungs):
            if os.path.exists(S.p("runs", rel, "exit_code")):
                out = render(S, rel, args.field)
                if out:
                    print("  ", out, "(inspect the _last.png frame before using it)")
