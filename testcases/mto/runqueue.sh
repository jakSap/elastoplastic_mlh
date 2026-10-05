#!/bin/sh
# Run a ranked queue file (queue/*.txt) with at most 4 concurrent solvers, in file order.
# Runs on battery and in Low Power Mode (user 10-05; timings are ~10x the estimates then).
cd "$(dirname "$0")"
Q=${1:-queue/2026-10-03.txt}
echo "$(date +%H:%M:%S) queue $Q started" >> runs/queue_jobs.log
grep -v '^#' "$Q" | grep -v '^ *$' | xargs -P 4 -L 1 ./queue_job.sh
echo "$(date +%H:%M:%S) queue $Q finished" >> runs/queue_jobs.log
