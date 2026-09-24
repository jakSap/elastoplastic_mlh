#!/usr/bin/env bash
# Pull back only what the analysis needs: the t = 1.0 snapshot of every run plus
# its run.log and config -- not the ~25 dumps per run.  Run from the laptop.
set -euo pipefail
DEST="${DEST:-$HOME/Documents/elastoplastic_mlh/testcases/tord_sweep/data}"
mkdir -p "$DEST"
rsync -av Mars:'~/Documents/tord_sweep/{progress.log,go.log,jobs.txt,variants.txt,summary.json,summary.txt}' "$DEST/" || true
# target snapshots: identified on Mars by stored /time, listed into targets.txt
ssh Mars 'cd ~/Documents/tord_sweep && ROOT=$PWD python3 - <<PY
import glob,os,h5py,numpy as np
sw=os.path.expanduser("~/Documents/tord_sweep/testcases/soundwave_2d")
out=[]
for d in sorted(glob.glob(sw+"/output_tord_*")):
    best=None
    for f in glob.glob(d+"/*.h5"):
        t=float(np.ravel(h5py.File(f,"r")["time"])[0])
        if abs(t-1.0)<1e-6 and (best is None or abs(t-1.0)<best[0]): best=(abs(t-1.0),f)
    if best: out.append(os.path.relpath(best[1], os.path.expanduser("~/Documents/tord_sweep")))
open(os.path.expanduser("~/Documents/tord_sweep/targets.txt"),"w").write("\n".join(out)+"\n")
print(len(out),"target snapshots")
PY'
rsync -av --files-from=<(ssh Mars 'cat ~/Documents/tord_sweep/targets.txt') \
      Mars:'~/Documents/tord_sweep/' "$DEST/tree/"
echo "fetched into $DEST"
