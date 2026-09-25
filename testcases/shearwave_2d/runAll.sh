#!/usr/bin/env bash
# Run the 2D shear-wave matrix from this directory as a job queue.
#
#   default: n in {64, 128, 256}, N_NN = 32, masses corrected, smallest first
#
# Each run logs to output_<run>/run.log; afterwards its y-averaged density
# movie is rendered to movies/<run>.mp4. One status line per finished run is
# appended to runs_status.log ("DONE ..." / "FAILED ..."), then "ALLDONE".
#
# Usage: ./runAll.sh [-j JOBS] [RUN ...]
#   RUN is a config stem, e.g. shear_N128_corr.
#   -j JOBS parallel runs (default 3 = one per resolution).

set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"

JOBS=3
if [ "${1:-}" = "-j" ]; then JOBS="$2"; shift 2; fi

RUNS=("$@")
if [ ${#RUNS[@]} -eq 0 ]; then
    for n in 64 128 256; do RUNS+=("shear_N${n}_corr"); done
fi

STATUS=runs_status.log
echo "QUEUE $(date '+%F %T') jobs=${JOBS} runs=${RUNS[*]}" >> "$STATUS"

run_one() {
    local run="$1" t0 rc movie
    t0=$(date +%s)
    mkdir -p "output_${run}"
    rc=0
    build_shear/bin/mlh -c "config_${run}.info" > "output_${run}/run.log" 2>&1 || rc=$?
    local el=$(( $(date +%s) - t0 ))
    local nerr nwarn
    nerr=$(grep -c "ERROR" "output_${run}/run.log" || true)
    nwarn=$(grep -c "WARN" "output_${run}/run.log" || true)
    if [ "$rc" -eq 0 ] && grep -q "FINISHED" "output_${run}/run.log"; then
        movie="none"
        if [ -f makeMovie.py ]; then
            movie="movies/${run}.mp4"
            python3 makeMovie.py "$run" > "output_${run}/movie.log" 2>&1 || movie="movie_failed"
        fi
        echo "DONE ${run} runtime=${el}s errors=${nerr} warns=${nwarn} movie=${movie} $(date '+%T')" >> runs_status.log
    else
        echo "FAILED ${run} exit=${rc} runtime=${el}s errors=${nerr} warns=${nwarn} $(date '+%T')" >> runs_status.log
    fi
}
export -f run_one

printf '%s\n' "${RUNS[@]}" | xargs -P "$JOBS" -n 1 bash -c 'run_one "$1"' _
echo "ALLDONE $(date '+%F %T')" >> "$STATUS"
