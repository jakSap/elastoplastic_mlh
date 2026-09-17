#!/usr/bin/env bash
# Create and compile the three build directories for the 2D soundwave test.
#
# Each build is a snapshot copy of demonstrator/ (Makefile, include/, src/,
# tests/; holib/ is symlinked) with include/parameter.h changed in exactly these
# flags relative to the demonstrator tree:
#
#   EOS                 0 | 1 | 2   (ideal gas | Murnaghan | Tillotson)
#   GIZMO_ELASTIC_FLUX  0           -> Riemann::HLLCFlux uses xSplitElasticHLLC
#   RECON_FACE_DV       0           (#error guard: needs GIZMO_ELASTIC_FLUX)
#   SURFACE_VOLCORR     0
#
# ELASTIC stays 1: xSplitElasticHLLC only exists on the ELASTIC branch, and
# EXPLICIT_VOL_INTEGRATION requires it. Solid dynamics are switched off through
# the shear modulus instead: mu = 0 in the material tables (EOS 1/2), and
# the EOS 0 path hard-codes mu = 0. With S = 0 at t = 0 and mu = 0 the
# deviatoric stress stays identically zero.
#
# Enforced to the requested values (already set that way in the tree):
#   ELASTIC 1, PERIODIC_BOUNDARIES 1, TENSILE_CORRECTION 1,
#   EXPLICIT_VOL_INTEGRATION 1, ELASTIC_HLLC_EP 0, DIM 2
#
# Usage: ./makeBuilds.sh [ideal|murnaghan|tillotson ...]   (default: all three)

set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DEMO="$(cd "$HERE/../../demonstrator" && pwd)"

eos_id() {
    case "$1" in
        ideal)     echo 0 ;;
        murnaghan) echo 1 ;;
        tillotson) echo 2 ;;
        *) echo "unknown EOS '$1'" >&2; exit 1 ;;
    esac
}

set_flag() {  # set_flag FILE NAME VALUE
    local file="$1" name="$2" value="$3"
    if ! grep -qE "^#define ${name} " "$file"; then
        echo "flag ${name} not found in ${file}" >&2; exit 1
    fi
    sed -i.bak -E "s|^#define ${name} .*$|#define ${name} ${value}|" "$file"
    rm -f "${file}.bak"
}

EOS_LIST=("$@")
[ ${#EOS_LIST[@]} -eq 0 ] && EOS_LIST=(ideal murnaghan tillotson)

for eos in "${EOS_LIST[@]}"; do
    id="$(eos_id "$eos")"
    B="$HERE/build_${eos}"
    echo "=== build_${eos} (EOS ${id}) ==="
    rm -rf "$B"
    mkdir -p "$B"
    cp -R "$DEMO/Makefile" "$DEMO/include" "$DEMO/src" "$DEMO/tests" "$B/"
    rm -f "$B/include/parameter.h.cmp_backup"
    ln -s "$DEMO/holib" "$B/holib"

    P="$B/include/parameter.h"
    set_flag "$P" EOS "$id"
    set_flag "$P" ELASTIC 1
    set_flag "$P" GIZMO_ELASTIC_FLUX 0
    set_flag "$P" RECON_FACE_DV 0
    set_flag "$P" ELASTIC_HLLC_EP 0
    set_flag "$P" SURFACE_VOLCORR 0
    set_flag "$P" TENSILE_CORRECTION 1
    set_flag "$P" EXPLICIT_VOL_INTEGRATION 1
    set_flag "$P" PERIODIC_BOUNDARIES 1
    set_flag "$P" DIM 2

    # record provenance next to the snapshot
    {
        echo "snapshot of demonstrator/ taken $(date '+%Y-%m-%d %H:%M')"
        git -C "$DEMO" rev-parse HEAD
        git -C "$DEMO" status --short -- . | sed 's/^/  /'
        echo
        echo "parameter.h diff vs demonstrator/include/parameter.h:"
        diff "$DEMO/include/parameter.h" "$P" || true
    } > "$B/SNAPSHOT.txt"

    (cd "$B" && make -j8 > make.log 2>&1) || { grep -E "error:|ERROR|aborted" "$B/make.log" >&2; exit 1; }
    [ -x "$B/bin/mlh" ] || { echo "build_${eos}: bin/mlh missing, see make.log" >&2; exit 1; }
    grep -E "flags are referenced|Linking" "$B/make.log"
done
