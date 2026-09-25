#!/usr/bin/env bash
# Is the first-order-in-dt damping of the ideal gas caused by the slope limiter?
#
# Hypothesis (2) from 2026-09-22-cfl-dump-und-signalgeschwindigkeit §8: the MUSCL
# half-step (Particles.cpp:4311-4334) consumes the LIMITED gradients. Clipping
# them breaks the time-centering and leaves an O(dt) jump at the face for the
# upwind flux to dissipate -- which would show up exactly as the measured
# p = 0.97 damping, 80% of the error at CFL 0.2.
#
# Test: the same three-point CFL ladder as cflTest.sh, but with both limiters
# OFF (SLOPE_LIMITING 0, PAIRWISE_LIMITER 0, as in limiterTest.sh's limNone).
#   p collapses to ~2  -> it is the limiter.
#   p stays ~1         -> it is not; next suspect is the adiabatic face
#                         correction / S_M frame (hypothesis 1).
#
# Builds: limiter_builds/build_<eos>_limNone_cfl<tag>; runs: output_<stem>_limNone[_cfl<tag>].
#
# Usage: ./limiterCflTest.sh [STEM ...]      (default: ideal_N64)

set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
HERE="$(pwd)"

CFLS=(.2 .05 .0125)

set_flag() {
    local file="$1" name="$2" value="$3"
    grep -qE "^#define ${name} " "$file" || { echo "flag ${name} not found in ${file}" >&2; exit 1; }
    sed -i.bak -E "s|^#define ${name} .*$|#define ${name} ${value}|" "$file"
    rm -f "${file}.bak"
}

make_build() {  # make_build <eos> <cfl>
    local eos="$1" cfl="$2" tag B P
    tag="cfl${cfl#.}"
    B="$HERE/limiter_builds/build_${eos}_limNone_${tag}"
    [ -x "$B/bin/mlh" ] && return 0
    local REF="$HERE/build_${eos}"
    echo "=== build_${eos}_limNone_${tag} ==="
    rm -rf "$B"; mkdir -p "$B"
    cp -R "$REF/Makefile" "$REF/include" "$REF/src" "$REF/tests" "$B/"
    rm -f "$B/include/parameter.h.cmp_backup"
    ln -s "$(readlink "$REF/holib")" "$B/holib"
    P="$B/include/parameter.h"
    set_flag "$P" SLOPE_LIMITING 0
    set_flag "$P" PAIRWISE_LIMITER 0
    set_flag "$P" CFL "$cfl"
    (cd "$B" && make -j4 > make.log 2>&1) || { grep -E "error:" "$B/make.log" >&2; exit 1; }
    [ -x "$B/bin/mlh" ] || { echo "$B: bin/mlh missing" >&2; exit 1; }
}

STEMS=("$@"); [ ${#STEMS[@]} -eq 0 ] && STEMS=(ideal_N64)
mkdir -p limiter_builds
for stem in "${STEMS[@]}"; do for c in "${CFLS[@]}"; do make_build "${stem%%_*}" "$c"; done; done

echo "QUEUE $(date '+%F %T') limNone ladder stems=${STEMS[*]} cfls=${CFLS[*]}" >> cfl_status.log

run_one() {  # run_one <stem> <cfl>
    local stem="$1" cfl="$2" eos="${1%%_*}" tag run t0 el rc=0
    tag="cfl${cfl#.}"
    # CFL 0.2 is the reference: no _cfl suffix, matching cflTest.sh's naming
    if [ "$cfl" = ".2" ]; then run="${stem}_limNone"; else run="${stem}_limNone_${tag}"; fi
    mkdir -p "output_${run}"
    sed "s|^outDir.*|outDir   output_${run}|" "config_${stem}.info" > "config_${run}.info"
    t0=$(date +%s)
    "limiter_builds/build_${eos}_limNone_${tag}/bin/mlh" -c "config_${run}.info" \
        > "output_${run}/run.log" 2>&1 || rc=$?
    el=$(( $(date +%s) - t0 ))
    if [ "$rc" -eq 0 ] && grep -q "FINISHED" "output_${run}/run.log"; then
        echo "DONE ${run} runtime=${el}s steps=$(grep -c 'selected' output_${run}/run.log)" >> cfl_status.log
    else
        echo "FAILED ${run} exit=${rc} runtime=${el}s" >> cfl_status.log
    fi
}
export -f run_one

for stem in "${STEMS[@]}"; do
    printf '%s\n' "${CFLS[@]}" | xargs -P 3 -I{} bash -c 'run_one "$1" "$2"' _ "$stem" {}
done
echo "ALLDONE $(date '+%F %T')" >> cfl_status.log

python3 cflAnalysis.py "${STEMS[@]/%/_limNone}" | tee plots/cfl_limNone_summary.txt
