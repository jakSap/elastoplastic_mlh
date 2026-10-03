#!/bin/sh
# One queue entry: "CANDIDATE TEST:RUNG" or "@sh <command ...>" (see runqueue.sh).
cd "$(dirname "$0")"
if [ "$1" = "@sh" ]; then shift; sh -c "$*" >> runs/queue_jobs.log 2>&1; exit; fi
echo "$(date +%H:%M:%S) start $1 $2" >> runs/queue_jobs.log
./mto.py run "$1" --only "$2" -j 1 >> runs/queue_jobs.log 2>&1
echo "$(date +%H:%M:%S) done  $1 $2 (exit $?)" >> runs/queue_jobs.log
