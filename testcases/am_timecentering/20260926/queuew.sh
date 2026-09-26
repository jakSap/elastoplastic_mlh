#!/usr/bin/env bash
# queuew.sh jobs MAXJ : lines "build test N tag" -> runw.sh; skips existing exit_code
cd "$(dirname "$0")"; J=$1; MAXJ=${2:-8}
while read -r b test n tag; do
  [ -z "$b" ] && continue; [[ "$b" == \#* ]] && continue
  [ -f runsw/$test/$tag/exit_code ] && continue
  while [ "$(pgrep -c -x mlh)" -ge "$MAXJ" ]; do sleep 5; done
  echo "$(date +%T) launch W $test/$tag" >> queue.log
  ./runw.sh $b $test $n $tag &
  sleep 2
done < "$J"
wait; echo "$(date +%T) QUEUE DONE $J" >> queue.log
