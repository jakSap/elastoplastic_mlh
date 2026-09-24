#!/usr/bin/env bash
# 6 workers (Mars has 6 physical cores) over jobs.txt, tier order preserved.
# Guard: skip a job if free disk drops below 5 G -- Mars sits at 95 % full.
ROOT="${ROOT:-$HOME/Documents/tord_sweep}"
cd "$ROOT"
export LC_ALL=C

runone() {
    IFS=$'\t' read -r dd build cfg out tag <<< "$1"
    avail=$(df --output=avail -BG "$HOME" | tail -1 | tr -dc '0-9')
    if [ "${avail:-0}" -lt 5 ]; then
        echo "$(date +%H:%M:%S) SKIP $tag (only ${avail}G free)" >> "$ROOT/progress.log"; return 0
    fi
    s=$(date +%s)
    if ( cd "$dd" && ./"$build"/bin/mlh -c "$cfg" > "$out/run.log" 2>&1 ); then st=OK; else st=FAIL; fi
    n=$(ls "$dd/$out"/*.h5 2>/dev/null | wc -l)
    echo "$(date +%H:%M:%S) $st $tag  $(( $(date +%s) - s ))s  ${n}files" >> "$ROOT/progress.log"
}
export -f runone
export ROOT

echo "START $(date)" >> progress.log
tr '\n' '\0' < jobs.txt | xargs -0 -P "${NPROC:-6}" -I{} bash -c 'runone "$@"' _ {}
echo "DONE $(date)" >> progress.log
df -h "$HOME" | tail -1 >> progress.log
