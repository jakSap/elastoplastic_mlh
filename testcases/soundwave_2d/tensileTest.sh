#!/usr/bin/env bash
# Tensile-correction experiment, demonstrator side: why does GIZMO's corrected-mass
# Murnaghan soundwave run 1.8 % slow while the demonstrator's does not?
#
# Murnaghan (n = 2, K0 = rho0 = 1), corrected masses (P crosses zero, so the
# dummy-pressure/tensile-correction branch acts), N_NN = 32, n = 64/128/256.
# Each variant is a copy of build_murnaghan/ (the snapshot of the reference
# runs output_murnaghan_N<n>_corr) differing ONLY in:
#
#   tcOff       TENSILE_CORRECTION 0
#   tcMonaghan  TENSILE_CORRECTION_LOCAL 0, TENSILE_CORRECTION_PREFAC 0.2,
#               TENSILE_CORRECTION_POWER 4  ->  f = 0.2 (W(r/h)/W(dp))^4, the form
#               GIZMO uses (solids/elastic_physics.c)
#
# The reference build has TENSILE_CORRECTION_LOCAL 1: f = 4 max(0, 1 - s_ab),
# s_ab = r_ab / mean((m/rho)^(1/DIM)).
#
# Builds: tensile_builds/build_murnaghan_<variant>; runs: output_murnaghan_N<n>_corr_<variant>.
# Usage: ./tensileTest.sh [-j JOBS]

set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
HERE="$(pwd)"
JOBS=6
if [ "${1:-}" = "-j" ]; then JOBS="$2"; shift 2; fi

REF="$HERE/build_murnaghan"
VARIANTS=(tcOff tcMonaghan)

set_flag() {  # set_flag FILE NAME VALUE
    local file="$1" name="$2" value="$3"
    grep -qE "^#define ${name} " "$file" || { echo "flag ${name} not found in ${file}" >&2; exit 1; }
    sed -i.bak -E "s|^#define ${name} .*$|#define ${name} ${value}|" "$file"
    rm -f "${file}.bak"
}

mkdir -p tensile_builds
for v in "${VARIANTS[@]}"; do
    B="$HERE/tensile_builds/build_murnaghan_${v}"
    echo "=== $B ==="
    rm -rf "$B"; mkdir -p "$B"
    cp -R "$REF/Makefile" "$REF/include" "$REF/src" "$REF/tests" "$B/"
    ln -s "$(readlink "$REF/holib")" "$B/holib"
    P="$B/include/parameter.h"
    case "$v" in
        tcOff)      set_flag "$P" TENSILE_CORRECTION 0 ;;
        tcMonaghan) set_flag "$P" TENSILE_CORRECTION_LOCAL 0
                    set_flag "$P" TENSILE_CORRECTION_PREFAC 0.2
                    set_flag "$P" TENSILE_CORRECTION_POWER 4. ;;
    esac
    { echo "tensile variant ${v}, copy of ${REF} taken $(date '+%Y-%m-%d %H:%M')"
      echo "parameter.h diff vs build_murnaghan:"; diff "$REF/include/parameter.h" "$P" || true
    } > "$B/SNAPSHOT.txt"
    (cd "$B" && make -j8 > make.log 2>&1) || { grep -E "error:|ERROR" "$B/make.log" >&2; exit 1; }
    [ -x "$B/bin/mlh" ] || { echo "${B}: bin/mlh missing" >&2; exit 1; }
done

RUNS=()
for v in "${VARIANTS[@]}"; do
    for n in 64 128 256; do
        run="murnaghan_N${n}_corr_${v}"
        { echo "; TENSILE EXPERIMENT (tensileTest.sh): variant ${v}, build tensile_builds/build_murnaghan_${v}"
          sed -E "s|^outDir .*$|outDir   output_${run}|" "config_murnaghan_N${n}_corr.info"
        } > "config_${run}.info"
        RUNS+=("$run")
    done
done

STATUS=tensile_status.log
echo "QUEUE $(date '+%F %T') jobs=${JOBS} runs=${RUNS[*]}" >> "$STATUS"
run_one() {
    local run="$1" v="${1##*_}" t0 rc
    t0=$(date +%s); mkdir -p "output_${run}"; rc=0
    "tensile_builds/build_murnaghan_${v}/bin/mlh" -c "config_${run}.info" > "output_${run}/run.log" 2>&1 || rc=$?
    local el=$(( $(date +%s) - t0 ))
    local nerr; nerr=$(grep -c "ERROR" "output_${run}/run.log" || true)
    if [ "$rc" -eq 0 ] && grep -q "FINISHED" "output_${run}/run.log"; then
        echo "DONE ${run} runtime=${el}s errors=${nerr} $(date '+%T')" >> tensile_status.log
    else
        echo "FAILED ${run} exit=${rc} runtime=${el}s errors=${nerr} $(date '+%T')" >> tensile_status.log
    fi
}
export -f run_one
printf '%s\n' "${RUNS[@]}" | xargs -P "$JOBS" -n 1 bash -c 'run_one "$1"' _
echo "ALLDONE $(date '+%F %T')" >> "$STATUS"
