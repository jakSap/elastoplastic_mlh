#!/usr/bin/env bash
# queue_mixed.sh jobs MAXJ : lines "L build ic tag tend ddt" (disc, run.sh) or "W build test N tag" (wave, runw.sh)
cd "$(dirname "$0")"; J=$1; MAXJ=${2:-8}
while read -r kind a b c d e; do
  [ -z "$kind" ] && continue
  if [ "$kind" = L ]; then [ -f runs/$c/exit_code ] && continue; else [ -f runsw/$b/$d/exit_code ] && continue; fi
  while [ "$(pgrep -c -x mlh)" -ge "$MAXJ" ]; do sleep 5; done
  echo "$(date +%T) launch $kind $c $d" >> queue.log
  if [ "$kind" = L ]; then ./run.sh $a $b $c $d $e & else ./runw.sh $a $b $c $d & fi
  sleep 2
done < "$J"
wait; echo "$(date +%T) QUEUE DONE $J" >> queue.log
