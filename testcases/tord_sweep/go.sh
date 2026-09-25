#!/usr/bin/env bash
# One entry point, meant for:  setsid nohup ./go.sh > go.log 2>&1 < /dev/null &
# Survives the ssh session going away.
set -euo pipefail
export ROOT="${ROOT:-$HOME/Documents/tord_sweep}"
cd "$ROOT"
[ -s progress.log ] && mv progress.log "progress_$(date +%Y%m%d_%H%M%S).log"
: > progress.log

echo "### setup $(date)"
./setup.sh

echo "### jobs $(date)"
python3 makeJobs.py

echo "### smoke test $(date)"
# One cheap run of `base` must land a snapshot on t = 1.0.  If the dump interval
# rule is wrong, every one of the ~170 jobs is wasted, so check it first.
SW="$ROOT/testcases/soundwave_2d"
( cd "$SW" && ./build_t_w3exact/bin/mlh -c config_tord_w3exact_N64_M192.info > output_tord_w3exact_N64_M192/run.log 2>&1 )
python3 - <<'PY'
import glob, os, sys, h5py
import numpy as np
sw = os.path.expanduser(os.environ["ROOT"]) + "/testcases/soundwave_2d"
d = sw + "/output_tord_w3exact_N64_M192"
ts = []
for f in sorted(glob.glob(d + "/*.h5")):
    with h5py.File(f, "r") as h:
        ts.append((float(np.ravel(h["time"])[0]), os.path.basename(f)))
hit = [x for x in ts if abs(x[0] - 1.0) < 1e-9]
print("snapshots:", ts)
if not hit:
    print("SMOKE FAILED: no snapshot at t = 1.0", file=sys.stderr); sys.exit(1)
print("SMOKE OK ->", hit)
PY

echo "### queue $(date)"
./runQueue.sh
echo "### all done $(date)"
