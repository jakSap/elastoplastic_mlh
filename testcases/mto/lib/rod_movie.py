#!/usr/bin/env python3
"""P1 movie: a candidate's Taylor rods above the GIZMO control, same times, coloured by density.

    python3 lib/rod_movie.py CAND RUNG [-o OUT.mp4] [--fps 8]

Particles whose nearest neighbour in the same rod is closer than 0.5 d0 (the clump criterion of
lib/analysis.py) are drawn black. The caption carries the mto verdict and summary of the run
(result.json). The GIZMO control is the one the analysis compared with (result.json 'gizmo_control').
"""
import argparse, json, os, shutil, subprocess, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from analysis import snapshots, load, intra_nn

ROOT = os.path.join(HERE, "..")
INK, INK2, SURF = "#0b0b0b", "#52514e", "#fcfcfb"
CMAP = LinearSegmentedColormap.from_list("div", ["#104281", "#3987e5", "#f0efec", "#e34948", "#8f1d1d"])


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("cand"); ap.add_argument("rung")
    ap.add_argument("-o", "--out"); ap.add_argument("--fps", type=float, default=8)
    a = ap.parse_args()
    D = os.path.join(ROOT, "runs", a.cand, "P1", a.rung)
    res = json.load(open(os.path.join(D, "result.json")))
    d0 = res.get("d0") or json.load(open(os.path.join(D, "meta.json")))["d"]
    ctrl = res.get("gizmo_control") or a.rung
    panels = [("{}  ({})".format(a.cand, res.get("verdict", "?")), snapshots(os.path.join(D, "output"))),
              ("GIZMO control ({})".format(ctrl), snapshots(os.path.join(ROOT, "gizmo", "runs", ctrl, "conv")))]
    times = [t for t, _ in panels[0][1] if t <= panels[1][1][-1][0] + 1e-9]
    out = a.out or os.path.join(D, "p1_{}_{}.mp4".format(a.cand, a.rung))

    plt.rcParams.update({"font.size": 10, "text.color": INK, "figure.facecolor": SURF, "axes.facecolor": SURF})
    X, Y = 1.25, 0.30
    fig, axs = plt.subplots(2, 1, figsize=(12, 7.2), dpi=100)
    fig.subplots_adjust(left=0.02, right=0.90, bottom=0.10, top=0.90, hspace=0.18)
    cax = fig.add_axes([0.92, 0.15, 0.015, 0.7])
    wpt = axs[0].get_window_extent().width * 72 / fig.dpi
    s = (0.85 * d0 / (2 * X) * wpt) ** 2
    summ = res.get("summary", "").split(" | ")
    caption = "mto: {}   ours:  {}".format(res.get("verdict", "?"), summ[0].split("  [")[0])
    if len(summ) > 1:
        caption += "\n           GIZMO: " + summ[1].split("  [")[0]
    if res.get("fails"):
        caption += "\n           fails: " + "; ".join(res["fails"])
    tmp = out + "_frames"; shutil.rmtree(tmp, ignore_errors=True); os.makedirs(tmp)
    for k, t in enumerate(times):
        for ax, (lab, sn) in zip(axs, panels):
            ax.cla()
            tt, path = min(sn, key=lambda q: abs(q[0] - t))
            f = load(path, ("x", "rho")); x = f["x"]
            cl = intra_nn(x, 2) / d0 < 0.5
            sc = ax.scatter(x[:, 0], x[:, 1], c=f["rho"], s=s, cmap=CMAP, vmin=0.8, vmax=1.2, linewidths=0)
            ax.scatter(x[cl, 0], x[cl, 1], s=2.2 * s, c=INK, linewidths=0)
            ax.set_xlim(-X, X); ax.set_ylim(-Y, Y); ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
            for sp in ax.spines.values():
                sp.set_visible(False)
            ax.set_title("{}   t = {:.2f}   clumped (black): {}".format(lab, tt, int(cl.sum())), fontsize=11)
        cb = fig.colorbar(sc, cax=cax); cb.set_label(r"density $\rho$", color=INK2); cb.outline.set_visible(False)
        fig.suptitle("P1 Taylor rods {}, {} vs GIZMO".format(a.rung, a.cand), fontsize=12)
        fig.texts = [tx for tx in fig.texts if tx.get_gid() != "cap"]
        fig.text(0.02, 0.02, caption, fontsize=9, family="monospace", color=INK2, gid="cap")
        fig.savefig(os.path.join(tmp, "f{:04d}.png".format(k)))
    plt.close(fig)
    subprocess.run(["ffmpeg", "-y", "-loglevel", "error", "-framerate", str(a.fps), "-i", os.path.join(tmp, "f%04d.png"),
                    "-vf", "tpad=stop_mode=clone:stop_duration=2,scale=trunc(iw/2)*2:trunc(ih/2)*2",
                    "-c:v", "libx264", "-pix_fmt", "yuv420p", "-crf", "20", out], check=True)
    shutil.rmtree(tmp, ignore_errors=True)
    print(out)


if __name__ == "__main__":
    main()
