#!/usr/bin/env bash
# Separate the TIME-splitting error from the SPATIAL error in the elastic wave.
#
# Why: in the n-refinement the deviatoric stress's deviation from its slaved
# value, max |S_xx + 4/3 mu ln(rho/rho0)| / (4/3 mu A), converges at order ~1
# while L1(rho) and L1(v_x) converge at order 2. Under the CFL condition
# dt ~ h, so a first-order-in-dt splitting error is indistinguishable from a
# first-order-in-h spatial error as long as only n is varied.
#
# This holds n fixed and cuts CFL (parameter.h, compile-time) from 0.2 to 0.05.
# A defect that drops by ~4 is O(dt); one that stays put is spatial.
#
# Measured (n = 64, after 3 passings):
#   CFL 0.20   max|S-defect| = 2.028e-02   L1(rho) = 2.358e-08   a_rho = 0.98320
#   CFL 0.05   max|S-defect| = 6.089e-03   L1(rho) = 1.763e-08   a_rho = 0.98863
# i.e. ratio 3.33 of the possible 4 -> ~93% of the defect is the time splitting,
# ~1.4e-3 is spatial. Note the density amplitude and phase improve too, so part
# of the wave damping is also first order in dt.
#
# Two CFL points cannot tell a dt term from a dt^2 term: both shrink when dt
# does. Use three (0.2 / 0.05 / 0.0125) and let cflAnalysis.py fit the exponent.
#
# Usage: ./cflTest.sh [-c CFL] [n ...]     (default: CFL .05, n = 64)

set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"

CFL_LOW=.05
if [ "${1:-}" = "-c" ]; then CFL_LOW="$2"; shift 2; fi
TAG="cfl${CFL_LOW#.}"
DEMO="$(cd ../../demonstrator && pwd)"
B="build_${TAG}"

if [ ! -x "$B/bin/mlh" ]; then
    echo "=== ${B} (as build_elastic, but CFL ${CFL_LOW}) ==="
    # NOTE: build from a FRESH copy of demonstrator/, never by copying an
    # existing build dir -- the stale obj/ and dep/ there make `make` a no-op
    # and you silently re-run the old binary.
    rm -rf "$B"; mkdir -p "$B"
    cp -R "$DEMO/Makefile" "$DEMO/include" "$DEMO/src" "$DEMO/tests" "$B/"
    rm -f "$B/include/parameter.h.cmp_backup"
    ln -s "$DEMO/holib" "$B/holib"
    for kv in "EOS 2" "ELASTIC 1" "GIZMO_ELASTIC_FLUX 0" "RECON_FACE_DV 0" \
              "ELASTIC_HLLC_EP 0" "SURFACE_VOLCORR 0" "TENSILE_CORRECTION 1" \
              "TENSILE_CORRECTION_LOCAL 1" "EXPLICIT_VOL_INTEGRATION 1" \
              "PERIODIC_BOUNDARIES 1" "DIM 2" "CFL ${CFL_LOW}"; do
        sed -i.bak -E "s|^#define ${kv%% *} .*$|#define ${kv%% *} ${kv#* }|" \
            "$B/include/parameter.h"
    done
    rm -f "$B/include/parameter.h.bak"
    (cd "$B" && make -j8 > make.log 2>&1) || { grep -E "error:" "$B/make.log" >&2; exit 1; }
fi

NS=("$@"); [ ${#NS[@]} -eq 0 ] && NS=(64)
for n in "${NS[@]}"; do
    run="elastic_N${n}_corr_${TAG}"
    mkdir -p "output_${run}"
    sed "s|^outDir.*|outDir   output_${run}|" "config_elastic_N${n}_corr.info" \
        > "config_${run}.info"
    echo "=== ${run} ==="
    "$B/bin/mlh" -c "config_${run}.info" > "output_${run}/run.log" 2>&1
    grep -c "STEP:" "output_${run}/run.log" | sed 's/^/  steps: /'
done

python3 cflAnalysis.py "${NS[@]}"
