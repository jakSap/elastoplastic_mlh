#!/usr/bin/env bash
# Time-order sweep, 2026-09-24.  Runs ON MARS in $ROOT (an rsync'd source tree).
#
# Question: the fixed-dt ladder of 2026-09-24 leaves p = 1.0 over six
# configurations with ADIABATIC_FACE_CORRECTION on -- only the prefactor moves
# (5.9 -> 2.1).  This sweep runs one build per candidate cause, each differing
# from `base` in exactly the flags listed in variants.txt, over the same dt
# ladder, ideal gas only.
#
# dt path: ADAPTIVE_TIMESTEP 0, dt = t_target/M exactly (identical to the
# 2026-09-23 sweep, so `base` must reproduce its numbers: p = 0.971,
# b = 7.2943e-3 at n = 64).
set -euo pipefail
ROOT="${ROOT:-$HOME/Documents/tord_sweep}"
cd "$ROOT"
DEMO="$ROOT/demonstrator"
SW="$ROOT/testcases/soundwave_2d"

set_flag() {  # set_flag FILE NAME VALUE
    local file="$1" name="$2" value="$3"
    grep -qE "^#define ${name} " "$file" || { echo "flag ${name} missing in ${file}" >&2; exit 1; }
    sed -i -E "s|^#define ${name} .*$|#define ${name} ${value}|" "$file"
}

echo "=== 1. fixed dt in the shared tree (every build inherits it) ==="
set_flag "$DEMO/include/parameter.h" ADAPTIVE_TIMESTEP 0
grep -nE '^#define (ADAPTIVE_TIMESTEP|MFM_ENERGY_STARFIX|DRIFT_KDK|RECON_TIME_PRED|ADIABATIC_FACE_CORRECTION) ' "$DEMO/include/parameter.h"

echo "=== 2. builds ==="
# Column order in variants.txt after the tag.
FLAGS=(MFM_ENERGY_STARFIX DRIFT_KDK ADIABATIC_FACE_CORRECTION RECON_TIME_PRED SLOPE_LIMITING FRAME_HALFSTEP TIME_DKD)
while IFS=$'\t' read -r tag v0 v1 v2 v3 v4 v5 v6; do
    case "$tag" in ''|\#*) continue ;; esac
    B="$SW/build_t_${tag}"
    echo "--- build_t_${tag} ---"
    rm -rf "$B"; mkdir -p "$B"
    cp -R "$DEMO/Makefile" "$DEMO/include" "$DEMO/src" "$DEMO/tests" "$B/"
    rm -f "$B/include/parameter.h.cmp_backup"
    ln -s "$DEMO/holib" "$B/holib"
    P="$B/include/parameter.h"
    # build_ideal baseline (makeBuilds.sh), then the per-variant flags
    for kv in "EOS 0" "ELASTIC 1" "GIZMO_ELASTIC_FLUX 0" "RECON_FACE_DV 0" \
              "ELASTIC_HLLC_EP 0" "SURFACE_VOLCORR 0" "TENSILE_CORRECTION 1" \
              "EXPLICIT_VOL_INTEGRATION 1" "PERIODIC_BOUNDARIES 1" "DIM 2" \
              "ADAPTIVE_TIMESTEP 0"; do
        set_flag "$P" "${kv% *}" "${kv#* }"
    done
    vals=("$v0" "$v1" "$v2" "$v3" "$v4" "${v5:-0}" "${v6:-0}")
    for k in 0 1 2 3 4 5 6; do set_flag "$P" "${FLAGS[$k]}" "${vals[$k]}"; done
    # PAIRWISE_LIMITER follows SLOPE_LIMITING (the limNone twin of the old sweep
    # turned both off together).
    [ "$v4" = "0" ] && set_flag "$P" PAIRWISE_LIMITER 0
    { echo "variant ${tag}  $(date '+%Y-%m-%d %H:%M')"
      for k in 0 1 2 3 4 5 6; do echo "  ${FLAGS[$k]} = ${vals[$k]}"; done
      echo; echo "parameter.h diff vs build_t_base:"; } > "$B/SNAPSHOT.txt"
    [ -f "$SW/build_t_base/include/parameter.h" ] && \
        diff "$SW/build_t_base/include/parameter.h" "$P" >> "$B/SNAPSHOT.txt" || true
    ( cd "$B" && make -j4 > make.log 2>&1 ) || { grep -E "error:" "$B/make.log" | head -20 >&2; echo "BUILD FAILED: $tag" >&2; exit 1; }
    [ -x "$B/bin/mlh" ] || { echo "build_t_${tag}: bin/mlh missing" >&2; exit 1; }
done < variants.txt

echo "=== 3. the single-flag diffs (this is the controlled part of the comparison) ==="
for B in "$SW"/build_t_*; do
    t="$(basename "$B" | sed 's/^build_t_//')"
    [ "$t" = base ] && continue
    printf '%-10s : ' "$t"
    { diff "$SW/build_t_base/include/parameter.h" "$B/include/parameter.h" || true; } \
      | grep -E '^[<>] #define' | tr '\n' ' ' || true ; echo
done

echo "=== 4. initial conditions (ideal gas, n = 64 and 128) ==="
( cd "$SW" && python3 generateIC.py --all > ic.log 2>&1 ) || { tail -20 "$SW/ic.log"; exit 1; }
ls -l "$SW"/soundwave_ideal_N64.h5 "$SW"/soundwave_ideal_N128.h5

echo "SETUP OK"
