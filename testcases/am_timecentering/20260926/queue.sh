#!/usr/bin/env bash
# queue.sh jobs.txt MAXJ : each line "build ic tag tend dumpdt"; runs up to MAXJ at once, in file order; skips tags with exit_code
cd "$(dirname "$0")"; J=$1; MAXJ=${2:-8}
while read -r b ic tag tend ddt; do
  [ -z "$b" ] && continue; [[ "$b" == \#* ]] && continue
  [ -f runs/$tag/exit_code ] && continue
  while [ "$(pgrep -c -x mlh)" -ge "$MAXJ" ]; do sleep 5; done
  avail=$(df --output=avail -k . | tail -1); [ "$avail" -lt 4000000 ] && { echo "$(date +%T) low disk, pausing" >> queue.log; while [ "$(df --output=avail -k . | tail -1)" -lt 4000000 ]; do sleep 60; done; }
  echo "$(date +%T) launch $tag" >> queue.log
  ./run.sh $b $ic $tag $tend $ddt &
  sleep 2
done < "$J"
wait; echo "$(date +%T) QUEUE DONE $J" >> queue.log
