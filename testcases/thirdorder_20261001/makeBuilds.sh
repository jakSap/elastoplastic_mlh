#!/usr/bin/env bash
# Builds for the thirdorder study (2026-10-01). Snapshot of ../../demonstrator (branch thirdorder).
# Common to all: PERIODIC_BOUNDARIES 1, ADAPTIVE_TIMESTEP 0, DIM 2, EOS 0 (ideal gas), ELASTIC 1 (default),
# GIZMO_ELASTIC_FLUX 0 (xSplitElasticHLLC carries everything; mu = 0 in the config) and the flags that are
# only valid with GEF 1 switched off (RECON_FACE_DV*, ELASTIC_FACE_MODE, TC gate/particle/damp-limit),
# POLY_RECONSTRUCTION 0 unless a variant sets it (periodic guard in parameter.h).
# Usage: ./makeBuilds.sh tag [tag ...]; variants in V[] below.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DEMO="$(cd "$HERE/../../demonstrator" && pwd)"
set_flag() { grep -qE "^#define $2 " "$1" || { echo "missing $2" >&2; exit 1; }
             sed -i -E "s|^#define $2 .*$|#define $2 $3|" "$1"; }
COMMON="PERIODIC_BOUNDARIES:1 ADAPTIVE_TIMESTEP:0 DIM:2 EOS:0 GIZMO_ELASTIC_FLUX:0 RECON_FACE_DV:0 RECON_FACE_DV_SPATIAL:0
 ELASTIC_FACE_MODE:0 TENSILE_CORRECTION_TRACTION_GATE:0 TENSILE_CORRECTION_PARTICLE:0 TENSILE_CORRECTION_DAMP_LIMIT:0
 POLY_RECONSTRUCTION:0"
declare -A V
V[b1]="GRADIENT_ORDER:1"
V[b1chk]="GRADIENT_ORDER:1"                   # = b1 rebuilt from the thirdorder tree (bit-identity check)                       # baseline = today's defaults (TIME_KDK) on the GEF 0 path
V[b2]="GRADIENT_ORDER:2"
V[b3]="GRADIENT_ORDER:3"
V[b1nl]="GRADIENT_ORDER:1 SLOPE_LIMITING:0"
V[b3nl]="GRADIENT_ORDER:3 SLOPE_LIMITING:0"
# --- thirdorder ingredients (SURFACE_GRAD_BLEND off: periodic box has no surface, and the blend is broken on
#     the periodic path -- gradBlend is never set there and the stress pass drops the ghost neighbours) ---
RK="TIME_RK3:1 TIME_KDK:0 TIME_KDK_FIELD_DRIFT:0"
V[b3fix]="GRADIENT_ORDER:3 SLOPE_LIMITING:0 SURFACE_GRAD_BLEND:0"                         # GO3 alone, blend bug removed
V[r1]="GRADIENT_ORDER:1 SLOPE_LIMITING:0 $RK"                                              # RK3 alone
V[k3f]="GRADIENT_ORDER:3 SLOPE_LIMITING:0 SURFACE_GRAD_BLEND:0 FACE_GRAD_WEIGHTS:1 EXPLICIT_VOL_SEED_IC:1"   # spatial, KDK
V[t3]="GRADIENT_ORDER:3 SLOPE_LIMITING:0 SURFACE_GRAD_BLEND:0 FACE_GRAD_WEIGHTS:1 EXPLICIT_VOL_SEED_IC:1 $RK" # all
V[t3ns]="GRADIENT_ORDER:3 SLOPE_LIMITING:0 SURFACE_GRAD_BLEND:0 FACE_GRAD_WEIGHTS:1 $RK"   # all but the IC seed
V[t3nf]="GRADIENT_ORDER:3 SLOPE_LIMITING:0 SURFACE_GRAD_BLEND:0 EXPLICIT_VOL_SEED_IC:1 $RK" # all but the faces
SP="GRADIENT_ORDER:3 SLOPE_LIMITING:0 PAIRWISE_LIMITER:0 SURFACE_GRAD_BLEND:0 EXPLICIT_VOL_SEED_IC:1 MOMENTUM_PREF:1"
V[t3p]="$SP FACE_GRAD_WEIGHTS:1 $RK"          # all + P_ref
V[k3p]="$SP FACE_GRAD_WEIGHTS:1"              # all + P_ref, KDK
V[t3pf1]="$SP $RK"                            # P_ref with first-order faces, linear recon
SP="GRADIENT_ORDER:3 SLOPE_LIMITING:0 PAIRWISE_LIMITER:0 SURFACE_GRAD_BLEND:0 EXPLICIT_VOL_SEED_IC:1 MOMENTUM_PREF:1"
V[t3p]="$SP FACE_GRAD_WEIGHTS:1 $RK"          # all + P_ref
V[k3p]="$SP FACE_GRAD_WEIGHTS:1"              # all + P_ref, KDK
V[t3pf1]="$SP $RK"                            # P_ref with first-order faces, linear recon
# ablation 2026-10-01: the full scheme t3p minus one ingredient each
FULL="$SP FACE_GRAD_WEIGHTS:1 $RK"
V[a_kdk]="${FULL/TIME_RK3:1 TIME_KDK:0 TIME_KDK_FIELD_DRIFT:0/}"   # time integrator back to TIME_KDK (default)
V[a_f1]="${FULL/FACE_GRAD_WEIGHTS:1/}"                             # first-order faces
V[a_seed]="${FULL/EXPLICIT_VOL_SEED_IC:1/}"                         # rhoExplicit seeded from the kernel sum
V[a_pref]="${FULL/MOMENTUM_PREF:1/}"                                # background pressure kept in the momentum flux
V[a_go2]="${FULL/GRADIENT_ORDER:3/GRADIENT_ORDER:2}"                # quadratic-exact weights
V[a_plim]="${FULL/PAIRWISE_LIMITER:0/PAIRWISE_LIMITER:1}"           # pairwise limiter on
V[a_slim]="${FULL/SLOPE_LIMITING:0/SLOPE_LIMITING:1}"               # slope limiter on
V[a_f1pref]="${V[a_f1]/MOMENTUM_PREF:1/}"                           # GO3 + seed + RK3 only (default faces, P0 kept)
V[hop]="GRADIENT_ORDER:3 SLOPE_LIMITING:0 PAIRWISE_LIMITER:0 SURFACE_GRAD_BLEND:0"   # Hopkins' recipe alone
# quadratic face reconstruction (RECON_HESSIAN): point-value corrected (2/3) vs Taylor (1)
V[q23]="$FULL RECON_HESSIAN:1"                                     # full + corrected quadratic recon
V[qT]="$FULL RECON_HESSIAN:1 RECON_HESSIAN_FAC:1.0"                # full + Taylor quadratic recon
V[f1q23]="${V[a_f1]} RECON_HESSIAN:1"                              # first-order faces + P_ref + corrected recon
[ -f "$HERE/variants.sh" ] && source "$HERE/variants.sh"
for t in "$@"; do
    [ -n "${V[$t]+x}" ] || { echo "unknown variant $t" >&2; exit 1; }
    B="$HERE/build_$t"; rm -rf "$B"; mkdir -p "$B"
    cp -R "$DEMO/Makefile" "$DEMO/include" "$DEMO/src" "$B/"
    [ -d "$DEMO/tests" ] && cp -R "$DEMO/tests" "$B/"
    rm -f "$B/include/parameter.h.cmp_backup"; [ -e "$B/holib" ] || ln -s "$DEMO/holib" "$B/holib"
    P="$B/include/parameter.h"
    for kv in $COMMON ${V[$t]}; do set_flag "$P" "${kv%:*}" "${kv#*:}"; done
    { echo "snapshot of demonstrator/ $(date '+%F %T')"; git -C "$DEMO" rev-parse HEAD
      git -C "$DEMO" status --short -- . ; echo; echo "variant $t: ${V[$t]}"; echo "parameter.h diff vs demonstrator:"
      diff "$DEMO/include/parameter.h" "$P" || true; git -C "$DEMO" diff -- . ; } > "$B/SNAPSHOT.txt"
    (cd "$B" && make -j4 > make.log 2>&1) || { grep -E "error" "$B/make.log" | head; exit 1; }
    [ -x "$B/bin/mlh" ] && echo "built $t"
done
