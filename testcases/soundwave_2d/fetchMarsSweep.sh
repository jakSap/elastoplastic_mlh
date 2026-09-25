#!/usr/bin/env bash
# Pull the 2026-09-23 fixed-dt sweep back from Mars.
#
# Only the snapshot AT the target time plus run.log per run -- each run holds ~70
# dumps (13 MB), so a blind rsync would drag gigabytes over the line. The target
# snapshot is selected on Mars by the stored /time, never by the file index:
# in the fixed-dt path `step` gets an extra ++ inside the dump block, so with
# dump interval I the j-th dump sits at t = j*(I-1)*dt and indices are not
# equidistant in time. See MA-Obsidian/2026-09-23-mars-fixed-dt-sweep.md.
set -euo pipefail
DEST="${1:-$HOME/Documents/elastoplastic_mlh/testcases/mars_fdt_sweep}"
mkdir -p "$DEST"

echo "=== selecting target-time snapshots on Mars ==="
ssh Mars 'cd ~/Documents/cfl_sweep && python3 - <<PY
import glob, os, re, shutil, h5py
T = {"ideal":1.0, "idealLimN":1.0, "tillotson":1.0, "murnaghan":1.0,
     "idealN128":1.0, "elastic":3*0.7867957924694431, "elastN128":3*0.7867957924694431,
     "shear":3*1.4719601443879746}
stage = os.path.expanduser("~/Documents/cfl_sweep/_stage")
shutil.rmtree(stage, ignore_errors=True); os.makedirs(stage)
rows = []
for d in sorted(glob.glob("testcases/*/output_fdt_*")):
    tag = os.path.basename(d)[len("output_fdt_"):]
    case = tag.rsplit("_M", 1)[0]
    if case not in T: continue
    best = None
    for f in glob.glob(d + "/*.h5"):
        try: t = float(h5py.File(f)["time"][()].ravel()[0])
        except Exception: continue
        if best is None or abs(t - T[case]) < abs(best[1] - T[case]): best = (f, t)
    if best is None:
        rows.append((tag, "NO SNAPSHOT", "")); continue
    o = os.path.join(stage, tag); os.makedirs(o, exist_ok=True)
    shutil.copy2(best[0], os.path.join(o, "target.h5"))
    for extra in ("run.log",):
        if os.path.exists(os.path.join(d, extra)): shutil.copy2(os.path.join(d, extra), o)
    cfg = "testcases/%s/config_fdt_%s.info" % (d.split("/")[1], tag)
    if os.path.exists(cfg): shutil.copy2(cfg, o)
    rows.append((tag, "%.17g" % best[1], "%.3e" % (best[1] - T[case])))
with open(os.path.join(stage, "MANIFEST.tsv"), "w") as fh:
    fh.write("tag\tt_snapshot\tdev_from_target\n")
    for r in rows: fh.write("\t".join(map(str, r)) + "\n")
print("staged %d runs" % len(rows))
bad = [r for r in rows if r[2] and abs(float(r[2])) > 1e-9]
print("OFF-TARGET:", bad if bad else "none")
PY'
echo "=== copying ==="
rsync -az --info=progress2 Mars:~/Documents/cfl_sweep/_stage/ "$DEST/"
rsync -az Mars:~/Documents/cfl_sweep/progress.log Mars:~/Documents/cfl_sweep/jobs.txt "$DEST/"
echo "=== done -> $DEST ==="
du -sh "$DEST"; column -t -s$'\t' "$DEST/MANIFEST.tsv" | head -20
