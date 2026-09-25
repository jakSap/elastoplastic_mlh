#!/usr/bin/env bash
# Separate the TIME-splitting error from the SPATIAL error in the SHEAR wave.
#
# Why here: on 2026-09-18 the elastic P-wave showed that the deviatoric stress's
# deviation from its slaved value falls with the timestep, not with the grid --
# 93% O(dt), 7% O(h) -- and that dt ~ h under the CFL condition hides this in a
# pure n-scan. The P-wave's restoring traction is only 38% deviatoric, so that
# finding was measured on a mode where the stress is a minority partner.
#
# The shear wave is the sharp version of the same test: its restoring force is
# 100% deviatoric and it carries no density perturbation at leading order. If
# the splitting offset is real, the coupling defect max|S_xy + rho c_s v_y|
# must be at least as dt-dominated here as it was there.
#
# This holds n fixed and cuts CFL (parameter.h, compile-time) from 0.2 to 0.05.
# A defect that drops by ~4 is O(dt); one that stays put is spatial.
#
# Usage: ./cflTest.sh [n ...]     (default: 64)

set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"

CFL_LOW=.05
DEMO="$(cd ../../demonstrator && pwd)"
B=build_cfl

if [ ! -x "$B/bin/mlh" ]; then
    echo "=== ${B} (as build_shear, but CFL ${CFL_LOW}) ==="
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
    run="shear_N${n}_corr_cfl05"
    mkdir -p "output_${run}"
    sed "s|^outDir.*|outDir   output_${run}|" "config_shear_N${n}_corr.info" \
        > "config_${run}.info"
    echo "=== ${run} ==="
    "$B/bin/mlh" -c "config_${run}.info" > "output_${run}/run.log" 2>&1
    grep -c "STEP:" "output_${run}/run.log" | sed 's/^/  steps: /'
done

python3 cflAnalysis.py "${NS[@]}"
