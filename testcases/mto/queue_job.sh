#!/bin/sh
# One queue entry: "CANDIDATE TEST:RUNG [note]" or "@sh <command ...>" (see runqueue.py).
# A job whose build fails (e.g. the source tree was mid-edit) is handed back to the queue: its line is
# removed from the state file after a pause, so the runner starts it again later.
cd "$(dirname "$0")"
if [ "$1" = "@sh" ]; then shift; sh -c "$*" >> runs/queue_jobs.log 2>&1; exit; fi
entry="$*"
echo "$(date +%H:%M:%S) start $1 $2" >> runs/queue_jobs.log
out=$(./mto.py run "$1" --only "$2" -j 1 2>&1)
rc=$?
echo "$out" >> runs/queue_jobs.log
if echo "$out" | grep -q "^build .* failed"; then
    echo "$(date +%H:%M:%S) build failed: $entry handed back to the queue (retry in 2 min)" >> runs/queue_jobs.log
    sleep 120
    python3 - "$entry" <<'PY'
import sys
S = "runs/queue_state.txt"; e = " ".join(sys.argv[1].split())
lines = [l for l in open(S).read().splitlines() if l.strip() != e]
open(S, "w").write("\n".join(lines) + "\n")
PY
    exit
fi
echo "$(date +%H:%M:%S) done  $1 $2 (exit $rc)" >> runs/queue_jobs.log
