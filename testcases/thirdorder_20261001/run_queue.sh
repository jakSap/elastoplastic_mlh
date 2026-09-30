#!/usr/bin/env bash
# Queue: args are run dirs holding config.info + BUILD. A job starts only while fewer than MAXBUSY
# mlh processes of THIS study run (other sessions' processes are not ours to count or touch).
cd "$(dirname "$0")"; MAXBUSY=${MAXBUSY:-4}; STUDY="$(pwd)"
busy() { for p in $(pgrep -x mlh); do readlink /proc/$p/cwd; done 2>/dev/null | grep -c "^$STUDY\$"; }
for r in "$@"; do
    r=${r%/}
    while [ "$(busy)" -ge "$MAXBUSY" ]; do sleep 5; done
    ( t0=$(date +%s); rc=0; rm -f $r/output/*.h5
      ./$(cat $r/BUILD)/bin/mlh -c $r/config.info > $r/run.log 2>&1 || rc=$?
      echo $rc > $r/exit_code
      echo "$r exit=$rc snaps=$(ls $r/output/*.h5 2>/dev/null | wc -l) wall=$(( $(date +%s)-t0 ))s $(date +%T)" >> summary.txt ) &
    echo "START $r $(date +%T)" >> queue.log; sleep 2
done
wait; echo "ALLDONE $(date +%T)" >> queue.log
