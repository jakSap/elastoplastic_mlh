#!/usr/bin/env bash
# Create and compile the build directory for the 2D elastic travelling wave.
#
# One build (Tillotson, EOS 2). It is a snapshot copy of demonstrator/ (Makefile,
# include/, src/, tests/; holib/ is symlinked) with include/parameter.h changed
# in exactly these flags relative to the demonstrator tree:
#
#   EOS                 2           (Tillotson)
#   GIZMO_ELASTIC_FLUX  0           -> Riemann::HLLCFlux uses xSplitElasticHLLC
#   RECON_FACE_DV       0           (#error guard: needs GIZMO_ELASTIC_FLUX)
#   SURFACE_VOLCORR     0
#
# Same flags as testcases/soundwave_2d, with one difference that is the whole
# point of this case: the shear modulus in the material table is mu > 0, so the
# deviatoric stress is live. The IC seeds S (see generateIC.py).
#
# Enforced to the requested values (already set that way in the tree):
#   ELASTIC 1, PERIODIC_BOUNDARIES 1, TENSILE_CORRECTION 1,
#   TENSILE_CORRECTION_LOCAL 1, EXPLICIT_VOL_INTEGRATION 1, ELASTIC_HLLC_EP 0, DIM 2
#
# Usage: ./makeBuilds.sh

set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DEMO="$(cd "$HERE/../../demonstrator" && pwd)"

set_flag() {  # set_flag FILE NAME VALUE
    local file="$1" name="$2" value="$3"
    if ! grep -qE "^#define ${name} " "$file"; then
        echo "flag ${name} not found in ${file}" >&2; exit 1
    fi
    sed -i.bak -E "s|^#define ${name} .*$|#define ${name} ${value}|" "$file"
    rm -f "${file}.bak"
}

B="$HERE/build_elastic"
echo "=== build_elastic (EOS 2, Tillotson, mu > 0) ==="
rm -rf "$B"
mkdir -p "$B"
cp -R "$DEMO/Makefile" "$DEMO/include" "$DEMO/src" "$DEMO/tests" "$B/"
rm -f "$B/include/parameter.h.cmp_backup"
ln -s "$DEMO/holib" "$B/holib"

P="$B/include/parameter.h"
set_flag "$P" EOS 2
set_flag "$P" ELASTIC 1
set_flag "$P" GIZMO_ELASTIC_FLUX 0
set_flag "$P" RECON_FACE_DV 0
set_flag "$P" ELASTIC_HLLC_EP 0
set_flag "$P" SURFACE_VOLCORR 0
set_flag "$P" TENSILE_CORRECTION 1
set_flag "$P" TENSILE_CORRECTION_LOCAL 1
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
[ -x "$B/bin/mlh" ] || { echo "build_elastic: bin/mlh missing, see make.log" >&2; exit 1; }
grep -E "flags are referenced|Linking" "$B/make.log"
