#!/usr/bin/env bash
# Run the 2D soundwave matrix from this directory as a job queue.
#
#   default: 3 EOS x N{64,128,256,512} x N_NN{32,64}, smallest first
#
# Memory, not CPU, limits the parallelism for the large runs (peak RSS, 16 GB
# machine): N64 ~80 MB, N128 ~250 MB, N256 ~0.9 GB, N512 ~3.6 GB (extrapolated,
# ~56 kB per particle). Run N512 with -j 2.
#
# Each run logs to output_<run>/run.log; afterwards its y-averaged density
# movie is rendered to movies/<run>.mp4. One status line per finished run is
# appended to runs_status.log ("DONE ..." / "FAILED ..."), then "ALLDONE".
#
# Usage: ./runAll.sh [-j JOBS] [RUN ...]
#   RUN is a config stem, e.g. ideal_N64 or tillotson_N128_nnn64.
#   -j JOBS parallel runs (default 5 = ~2/3 of the 8 logical CPUs).

set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"

JOBS=5
if [ "${1:-}" = "-j" ]; then JOBS="$2"; shift 2; fi

RUNS=("$@")
if [ ${#RUNS[@]} -eq 0 ]; then
    for n in 64 128 256 512; do
        for sfx in "" "_nnn64"; do
            for eos in ideal murnaghan tillotson; do RUNS+=("${eos}_N${n}${sfx}"); done
        done
    done
fi

STATUS=runs_status.log
echo "QUEUE $(date '+%F %T') jobs=${JOBS} runs=${RUNS[*]}" >> "$STATUS"

run_one() {
    local run="$1" eos="${1%%_*}" t0 rc movie
    t0=$(date +%s)
    mkdir -p "output_${run}"
    rc=0
    "build_${eos}/bin/mlh" -c "config_${run}.info" > "output_${run}/run.log" 2>&1 || rc=$?
    local el=$(( $(date +%s) - t0 ))
    local nerr nwarn
    nerr=$(grep -c "ERROR" "output_${run}/run.log" || true)
    nwarn=$(grep -c "WARN" "output_${run}/run.log" || true)
    if [ "$rc" -eq 0 ] && grep -q "FINISHED" "output_${run}/run.log"; then
        movie="movies/${run}.mp4"
        python3 makeMovie.py "$run" > "output_${run}/movie.log" 2>&1 || movie="movie_failed"
        echo "DONE ${run} runtime=${el}s errors=${nerr} warns=${nwarn} movie=${movie} $(date '+%T')" >> runs_status.log
    else
        echo "FAILED ${run} exit=${rc} runtime=${el}s errors=${nerr} warns=${nwarn} $(date '+%T')" >> runs_status.log
    fi
}
export -f run_one

printf '%s\n' "${RUNS[@]}" | xargs -P "$JOBS" -n 1 bash -c 'run_one "$1"' _
echo "ALLDONE $(date '+%F %T')" >> "$STATUS"
