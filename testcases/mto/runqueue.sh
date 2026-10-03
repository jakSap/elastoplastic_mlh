#!/bin/sh
# Run a ranked queue file (queue/*.txt) with at most 4 concurrent solvers, in file order.
# Waits until the laptop is on mains power with Low Power Mode off (timings are ~10x off otherwise).
cd "$(dirname "$0")"
Q=${1:-queue/2026-10-03.txt}
until pmset -g batt | grep -q "AC Power" && pmset -g | grep -q "lowpowermode *0"; do
    echo "$(date +%H:%M:%S) waiting for AC + Low Power Mode off" >> runs/queue_jobs.log; sleep 120
done
echo "$(date +%H:%M:%S) queue $Q started" >> runs/queue_jobs.log
grep -v '^#' "$Q" | grep -v '^ *$' | xargs -P 4 -L 1 ./queue_job.sh
echo "$(date +%H:%M:%S) queue $Q finished" >> runs/queue_jobs.log
