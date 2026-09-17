#!/usr/bin/env bash
# Limiter experiment, demonstrator side: are the error spikes at the wave extrema
# (and their factor ~2 over GIZMO) caused by the slope/pairwise limiters?
#
# Tillotson, N_NN = 32, n = 64 and 128, the same uncorrected ICs as the
# reference runs (output_tillotson_N<n>). Each variant is a copy of build_tillotson/
# (the snapshot the reference runs were made with) differing ONLY in:
#
#   limNone   SLOPE_LIMITING 0, PAIRWISE_LIMITER 0
#   limSL     SLOPE_LIMITING 1, PAIRWISE_LIMITER 0
#   limPW     SLOPE_LIMITING 0, PAIRWISE_LIMITER 1
#   limSLfix  SLOPE_LIMITING 1, PAIRWISE_LIMITER 1, slope limiter as in Hopkins (2015)
#             eq. B2: grad *= min(1, BETA*alpha). The demonstrator applies grad *= alpha
#             when BETA*alpha < 1 (Particles::slopeLimiter), i.e. BETA times stronger
#             and discontinuous at the threshold.
#
# Builds: limiter_builds/build_tillotson_<variant>; runs: output_tillotson_N<n>_<variant>.
# The GIZMO counterpart is ~/Documents/GIZMO/soundwave_2d/limiterTest.sh.
#
# Usage: ./limiterTest.sh [-j JOBS]

set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
HERE="$(pwd)"
JOBS=5
if [ "${1:-}" = "-j" ]; then JOBS="$2"; shift 2; fi

REF="$HERE/build_tillotson"
VARIANTS=(limNone limSL limPW limSLfix)

set_flag() {  # set_flag FILE NAME VALUE
    local file="$1" name="$2" value="$3"
    grep -qE "^#define ${name} " "$file" || { echo "flag ${name} not found in ${file}" >&2; exit 1; }
    sed -i.bak -E "s|^#define ${name} .*$|#define ${name} ${value}|" "$file"
    rm -f "${file}.bak"
}

mkdir -p limiter_builds
for v in "${VARIANTS[@]}"; do
    B="$HERE/limiter_builds/build_tillotson_${v}"
    echo "=== $B ==="
    rm -rf "$B"; mkdir -p "$B"
    cp -R "$REF/Makefile" "$REF/include" "$REF/src" "$REF/tests" "$B/"
    ln -s "$(readlink "$REF/holib")" "$B/holib"
    P="$B/include/parameter.h"
    case "$v" in
        limNone)  set_flag "$P" SLOPE_LIMITING 0; set_flag "$P" PAIRWISE_LIMITER 0 ;;
        limSL)    set_flag "$P" SLOPE_LIMITING 1; set_flag "$P" PAIRWISE_LIMITER 0 ;;
        limPW)    set_flag "$P" SLOPE_LIMITING 0; set_flag "$P" PAIRWISE_LIMITER 1 ;;
        limSLfix) set_flag "$P" SLOPE_LIMITING 1; set_flag "$P" PAIRWISE_LIMITER 1
                  sed -i.bak -E 's/grad\[i\]\[([012])\] \*= alpha(Min|Max);/grad[i][\1] *= BETA*alpha\2;/' "$B/src/Particles.cpp"
                  rm -f "$B/src/Particles.cpp.bak"
                  [ "$(grep -c 'grad\[i\]\[[012]\] \*= BETA\*alpha' "$B/src/Particles.cpp")" -eq 6 ] \
                      || { echo "slope limiter patch did not apply" >&2; exit 1; } ;;
    esac
    { echo "limiter variant ${v}, copy of ${REF} taken $(date '+%Y-%m-%d %H:%M')"
      echo "parameter.h diff vs build_tillotson:"; diff "$REF/include/parameter.h" "$P" || true
      echo "Particles.cpp diff vs build_tillotson:"; diff "$REF/src/Particles.cpp" "$B/src/Particles.cpp" || true
    } > "$B/SNAPSHOT.txt"
    (cd "$B" && make -j8 > make.log 2>&1) || { grep -E "error:|ERROR" "$B/make.log" >&2; exit 1; }
    [ -x "$B/bin/mlh" ] || { echo "${B}: bin/mlh missing" >&2; exit 1; }
done

RUNS=()
for v in "${VARIANTS[@]}"; do
    for n in 64 128; do
        run="tillotson_N${n}_${v}"
        { echo "; LIMITER EXPERIMENT (limiterTest.sh): variant ${v}, build limiter_builds/build_tillotson_${v}"
          sed -E "s|^outDir .*$|outDir   output_${run}|" "config_tillotson_N${n}.info"
        } > "config_${run}.info"
        RUNS+=("$run")
    done
done

STATUS=limiter_status.log
echo "QUEUE $(date '+%F %T') jobs=${JOBS} runs=${RUNS[*]}" >> "$STATUS"
run_one() {
    local run="$1" v="${1##*_}" t0 rc
    t0=$(date +%s); mkdir -p "output_${run}"; rc=0
    "limiter_builds/build_tillotson_${v}/bin/mlh" -c "config_${run}.info" > "output_${run}/run.log" 2>&1 || rc=$?
    local el=$(( $(date +%s) - t0 ))
    local nerr; nerr=$(grep -c "ERROR" "output_${run}/run.log" || true)
    if [ "$rc" -eq 0 ] && grep -q "FINISHED" "output_${run}/run.log"; then
        echo "DONE ${run} runtime=${el}s errors=${nerr} $(date '+%T')" >> limiter_status.log
    else
        echo "FAILED ${run} exit=${rc} runtime=${el}s errors=${nerr} $(date '+%T')" >> limiter_status.log
    fi
}
export -f run_one
printf '%s\n' "${RUNS[@]}" | xargs -P "$JOBS" -n 1 bash -c 'run_one "$1"' _
echo "ALLDONE $(date '+%F %T')" >> "$STATUS"
