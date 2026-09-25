#!/usr/bin/env bash
# Separate the TIME-integration error from the SPATIAL error in the 2D soundwave
# (mu = 0), for all three EOS.
#
# Why: the 17.09. numbers (damping deficit 1 - a_rho, dispersion c_num - c_exp,
# L1) were taken from a pure n-scan at fixed CFL 0.2. Under the CFL condition
# dt = CFL * h / v_sig ~ h, so a first-order-in-dt error is indistinguishable
# from a first-order-in-h spatial error in such a scan -- this is what bit the
# elastic wave on 18.09. (soundwave_2d_elastic/cflTest.sh: 93% of the stress
# slaving defect was O(dt)). The ideal gas damping at order 1 while GIZMO gives
# order 2 is the obvious candidate here.
#
# This holds n fixed and cuts CFL (parameter.h, compile-time) from 0.2 to 0.05.
# A quantity that drops by the dt ratio (~4) is time integration; one that stays
# put is spatial.
#
# Two CFL points cannot tell a dt term from a dt^2 term: both drop when dt does.
# Use at least three (0.2 / 0.05 / 0.0125) and let cflAnalysis.py fit the exponent.
#
# Usage: ./cflTest.sh [-c CFL] [STEM ...]
#   -c CFL       low CFL to compare against the 0.2 reference (default .05);
#                outputs go to output_<STEM>_cfl<digits>, e.g. _cfl05, _cfl0125
#   default stems: ideal_N64 murnaghan_N64_corr tillotson_N64
#   Each STEM needs an existing config_<STEM>.info (the CFL 0.2 reference run).

set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"

CFL_LOW=.05
if [ "${1:-}" = "-c" ]; then CFL_LOW="$2"; shift 2; fi
TAG="cfl${CFL_LOW#.}"
DEMO="$(cd ../../demonstrator && pwd)"

eos_id() {
    case "$1" in
        ideal)     echo 0 ;;
        murnaghan) echo 1 ;;
        tillotson) echo 2 ;;
        *) echo "unknown EOS '$1'" >&2; exit 1 ;;
    esac
}

make_build() {  # make_build <eos>
    local eos="$1" id B P
    id="$(eos_id "$eos")"
    B="build_${eos}_${TAG}"
    [ -x "$B/bin/mlh" ] && return 0
    echo "=== ${B} (as build_${eos}, but CFL ${CFL_LOW}) ==="
    # NOTE: build from a FRESH copy of demonstrator/, never by copying an
    # existing build dir -- the stale obj/ and dep/ there make `make` a no-op
    # and you silently re-run the old binary (cost us a round on 18.09.).
    rm -rf "$B"; mkdir -p "$B"
    cp -R "$DEMO/Makefile" "$DEMO/include" "$DEMO/src" "$DEMO/tests" "$B/"
    rm -f "$B/include/parameter.h.cmp_backup"
    ln -s "$DEMO/holib" "$B/holib"
    P="$B/include/parameter.h"
    # same flags as makeBuilds.sh, plus CFL
    for kv in "EOS ${id}" "ELASTIC 1" "GIZMO_ELASTIC_FLUX 0" "RECON_FACE_DV 0" \
              "ELASTIC_HLLC_EP 0" "SURFACE_VOLCORR 0" "TENSILE_CORRECTION 1" \
              "EXPLICIT_VOL_INTEGRATION 1" "PERIODIC_BOUNDARIES 1" "DIM 2" \
              "CFL ${CFL_LOW}"; do
        grep -qE "^#define ${kv%% *} " "$P" || { echo "flag ${kv%% *} not found" >&2; exit 1; }
        sed -i.bak -E "s|^#define ${kv%% *} .*$|#define ${kv%% *} ${kv#* }|" "$P"
    done
    rm -f "$P.bak"
    {
        echo "snapshot of demonstrator/ taken $(date '+%Y-%m-%d %H:%M')"
        git -C "$DEMO" rev-parse HEAD
        git -C "$DEMO" status --short -- . | sed 's/^/  /'
        echo
        echo "parameter.h diff vs demonstrator/include/parameter.h:"
        diff "$DEMO/include/parameter.h" "$P" || true
    } > "$B/SNAPSHOT.txt"
    (cd "$B" && make -j4 > make.log 2>&1) || { grep -E "error:" "$B/make.log" >&2; exit 1; }
    [ -x "$B/bin/mlh" ] || { echo "${B}: bin/mlh missing, see make.log" >&2; exit 1; }
}

STEMS=("$@")
[ ${#STEMS[@]} -eq 0 ] && STEMS=(ideal_N64 murnaghan_N64_corr tillotson_N64)

STATUS=cfl_status.log
echo "QUEUE $(date '+%F %T') cfl=${CFL_LOW} tag=${TAG} stems=${STEMS[*]}" >> "$STATUS"

for stem in "${STEMS[@]}"; do make_build "${stem%%_*}"; done

run_one() {
    local stem="$1" eos="${1%%_*}" run="${1}_${TAG}" t0 el rc=0
    [ -f "config_${stem}.info" ] || { echo "MISSING config_${stem}.info" >> cfl_status.log; return; }
    mkdir -p "output_${run}"
    sed "s|^outDir.*|outDir   output_${run}|" "config_${stem}.info" > "config_${run}.info"
    t0=$(date +%s)
    "build_${eos}_${TAG}/bin/mlh" -c "config_${run}.info" > "output_${run}/run.log" 2>&1 || rc=$?
    el=$(( $(date +%s) - t0 ))
    local nerr steps
    nerr=$(grep -c "ERROR" "output_${run}/run.log" || true)
    steps=$(grep -c "dt = .* selected" "output_${run}/run.log" || true)
    if [ "$rc" -eq 0 ] && grep -q "FINISHED" "output_${run}/run.log"; then
        echo "DONE ${run} runtime=${el}s steps=${steps} errors=${nerr} $(date '+%T')" >> cfl_status.log
    else
        echo "FAILED ${run} exit=${rc} runtime=${el}s steps=${steps} errors=${nerr} $(date '+%T')" >> cfl_status.log
    fi
}
export -f run_one
export TAG

printf '%s\n' "${STEMS[@]}" | xargs -P 3 -n 1 bash -c 'run_one "$1"' _
echo "ALLDONE $(date '+%F %T')" >> "$STATUS"

python3 cflAnalysis.py "${STEMS[@]}" | tee plots/cfl_summary.txt
