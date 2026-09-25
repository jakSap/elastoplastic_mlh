#!/usr/bin/env bash
# Does the adiabatic (entropic) face correction change convergence?
#
# Builds two snapshots of demonstrator/ that differ in EXACTLY ONE flag,
# ADIABATIC_FACE_CORRECTION (see demonstrator/include/parameter.h; Hopkins 2015,
# GIZMO hydro/hydro_core_meshless.h:335-375), and runs the ideal-gas 2D
# travelling soundwave at N = 64, 128, 256 through both.
#
# Everything else matches build_ideal/ from makeBuilds.sh, so afc0 reproduces
# the existing output_ideal_N<n> runs bit for bit:
#   EOS 0, ELASTIC 1 (mu = 0, so pure hydro), GIZMO_ELASTIC_FLUX 0
#   -> Riemann::HLLCFlux uses xSplitElasticHLLC, which is the path the
#   correction hooks into. RECON_FACE_DV 0, SURFACE_VOLCORR 0,
#   TENSILE_CORRECTION 1, EXPLICIT_VOL_INTEGRATION 1, PERIODIC_BOUNDARIES 1.
#
# Output dirs are output_ideal_N<n>_afc0 / _afc1, which convergence.py's
# analyse(..., suffix=) already understands; adiabaticFaceAnalysis.py reads them.
#
# Usage: ./adiabaticFaceTest.sh [-j JOBS] [N ...]        (default: 64 128 256)

set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"

HERE="$(pwd)"
DEMO="$(cd ../../demonstrator && pwd)"

JOBS=3
if [ "${1:-}" = "-j" ]; then JOBS="$2"; shift 2; fi
NS=("$@")
[ ${#NS[@]} -eq 0 ] && NS=(64 128 256)

set_flag() {  # set_flag FILE NAME VALUE
    local file="$1" name="$2" value="$3"
    if ! grep -qE "^#define ${name} " "$file"; then
        echo "flag ${name} not found in ${file}" >&2; exit 1
    fi
    sed -i.bak -E "s|^#define ${name} .*$|#define ${name} ${value}|" "$file"
    rm -f "${file}.bak"
}

build_one() {  # build_one AFC
    local afc="$1" B="$HERE/build_ideal_afc$1"
    echo "=== build_ideal_afc${afc} (ADIABATIC_FACE_CORRECTION ${afc}) ==="
    rm -rf "$B"; mkdir -p "$B"
    cp -R "$DEMO/Makefile" "$DEMO/include" "$DEMO/src" "$DEMO/tests" "$B/"
    rm -f "$B/include/parameter.h.cmp_backup"
    ln -s "$DEMO/holib" "$B/holib"

    local P="$B/include/parameter.h"
    set_flag "$P" EOS 0
    set_flag "$P" ELASTIC 1
    set_flag "$P" GIZMO_ELASTIC_FLUX 0
    set_flag "$P" RECON_FACE_DV 0
    set_flag "$P" ELASTIC_HLLC_EP 0
    set_flag "$P" SURFACE_VOLCORR 0
    set_flag "$P" TENSILE_CORRECTION 1
    set_flag "$P" EXPLICIT_VOL_INTEGRATION 1
    set_flag "$P" PERIODIC_BOUNDARIES 1
    set_flag "$P" DIM 2
    set_flag "$P" ADIABATIC_FACE_CORRECTION "$afc"
    set_flag "$P" ADIABATIC_FACE_DIAG 0

    {
        echo "snapshot of demonstrator/ taken $(date '+%Y-%m-%d %H:%M')"
        git -C "$DEMO" rev-parse HEAD
        git -C "$DEMO" status --short -- . | sed 's/^/  /'
        echo
        echo "parameter.h diff vs demonstrator/include/parameter.h:"
        diff "$DEMO/include/parameter.h" "$P" || true
    } > "$B/SNAPSHOT.txt"

    (cd "$B" && make -j8 > make.log 2>&1) || { grep -E "error:" "$B/make.log" >&2; exit 1; }
    [ -x "$B/bin/mlh" ] || { echo "build_ideal_afc${afc}: bin/mlh missing" >&2; exit 1; }
}

for afc in 0 1; do build_one "$afc"; done

# The two binaries must differ only through the flag: a diff of the two
# parameter.h files that shows anything but ADIABATIC_FACE_CORRECTION means the
# comparison is not controlled.
echo "=== parameter.h delta between the two builds ==="
diff "$HERE/build_ideal_afc0/include/parameter.h" \
     "$HERE/build_ideal_afc1/include/parameter.h" || true

run_one() {  # run_one N AFC
    local n="$1" afc="$2" run="ideal_N${n}_afc${afc}" cfg="config_ideal_N${n}_afc${afc}.info"
    sed -e "s|^outDir .*|outDir output_${run}|" "config_ideal_N${n}.info" > "$cfg"
    mkdir -p "output_${run}"
    local t0 rc; t0=$(date +%s); rc=0
    "build_ideal_afc${afc}/bin/mlh" -c "$cfg" > "output_${run}/run.log" 2>&1 || rc=$?
    if [ $rc -eq 0 ]; then
        echo "DONE ${run} $(( $(date +%s) - t0 ))s"
    else
        echo "FAILED ${run} rc=${rc} -- see output_${run}/run.log" >&2
    fi
}
export -f run_one

for n in "${NS[@]}"; do
    for afc in 0 1; do echo "$n $afc"; done
done | xargs -P "$JOBS" -n 2 bash -c 'run_one "$0" "$1"'

echo "ALLDONE -- now run ./adiabaticFaceAnalysis.py"
