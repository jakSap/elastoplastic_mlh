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
# fourth order in time (TIME_RK4) on top of the corrected quadratic reconstruction
RK4="TIME_RK4:1 TIME_KDK:0 TIME_KDK_FIELD_DRIFT:0"
V[q23r4]="${V[q23]/TIME_RK3:1 TIME_KDK:0 TIME_KDK_FIELD_DRIFT:0/$RK4}"
V[f1q23r4]="${V[f1q23]/TIME_RK3:1 TIME_KDK:0 TIME_KDK_FIELD_DRIFT:0/$RK4}"
# elastic branch, mu > 0 (2026-10-01 afternoon): Murnaghan (EOS 1) / Tillotson (EOS 2), ICs generateIC_elastic.py.
# Suffix m = Murnaghan, t = Tillotson; nc = TENSILE_CORRECTION 0.
for e in m:1 t:2; do
    s=${e%:*}; E="EOS:${e#*:}"
    V[e1$s]="GRADIENT_ORDER:1 $E"                    # today's GEF 0 scheme (2nd-order baseline of 2026-09-25)
    V[t3p$s]="$FULL $E"                              # third-order recipe as for the ideal gas
    V[t3p${s}nc]="$FULL $E TENSILE_CORRECTION:0"
    V[e1${s}x2]="GRADIENT_ORDER:1 $E XSPLIT_SHEAR_STAR:2"   # today's GEF 0 scheme + the tangential-traction fix
    V[q23r4$s]="${V[q23r4]} $E"                      # fourth order: RK4 + 2/3-quadratic recon incl. S
    V[q23r4${s}ls]="${V[q23r4]} $E RECON_HESSIAN_STRESS:0"   # ablation: S reconstructed linearly
    V[t3p${s}x1]="$FULL $E XSPLIT_SHEAR_STAR:1"      # tangential traction central (no side switch)
    V[t3p${s}x2]="$FULL $E XSPLIT_SHEAR_STAR:2"      # tangential traction from the shear-wave star state
    V[q23r4${s}x2]="${V[q23r4]} $E XSPLIT_SHEAR_STAR:2"
done
V[t3pchk]="$FULL"                                    # t3p rebuilt after the elastic changes (bit-identity check)
V[t3pmchk]="$FULL EOS:1"                             # t3pm rebuilt from the new tree (RK Euler stress stage, flag 0)
# rotating disc (test plan A.1: Tillotson, mu = 0.22, free surface, fixed dt), 2026-10-01 evening. Non-periodic, so no
# MOMENTUM_PREF (guarded: at a free surface P_ref would be an ambient pressure; the free-surface value is 0).
# NOCOMMON = demonstrator defaults without the periodic COMMON block.
DSK="PERIODIC_BOUNDARIES:0 EOS:2"
V[dstd]="NOCOMMON PERIODIC_BOUNDARIES:0 ADAPTIVE_TIMESTEP:0"          # standard scheme (GEF 1, defaults), fixed dt
V[dxs]="GRADIENT_ORDER:1 $DSK"                                         # today's GEF 0 scheme (2nd order)
V[d3]="${FULL/MOMENTUM_PREF:1/} $DSK XSPLIT_SHEAR_STAR:2 SURFACE_GRAD_BLEND:1"   # third order + surface fade
V[d3n]="${FULL/MOMENTUM_PREF:1/} $DSK XSPLIT_SHEAR_STAR:2"                      # third order, hard fallback only
V[d3x0]="${FULL/MOMENTUM_PREF:1/} $DSK SURFACE_GRAD_BLEND:1"                    # third order without the fix
V[d4]="${V[q23r4]/MOMENTUM_PREF:1/} $DSK XSPLIT_SHEAR_STAR:2 SURFACE_GRAD_BLEND:1" # fourth order (RK4 + quad recon)
# disc ablation of d3 (the free surface goes unstable after ~0.15 rotations at d0.05): one ingredient changed each
V[d3f1]="${V[d3]/FACE_GRAD_WEIGHTS:1/}"                                # first-order faces
V[d3lim]="${V[d3]/SLOPE_LIMITING:0 PAIRWISE_LIMITER:0/SLOPE_LIMITING:1 PAIRWISE_LIMITER:1}"   # both limiters on
V[d3ns]="${V[d3]/EXPLICIT_VOL_SEED_IC:1/}"                             # kernel-sum density seed
V[d3g1]="${V[d3f1]/GRADIENT_ORDER:3/GRADIENT_ORDER:1}"                 # linear-exact gradients (faces first order)
V[d3kdk]="${V[d3]/TIME_RK3:1 TIME_KDK:0 TIME_KDK_FIELD_DRIFT:0/}"     # TIME_KDK instead of RK3
# fix: cubic faces faded to first order at the surface with gradBlend (FACE_GRAD_WEIGHTS 2)
V[d3fb]="${V[d3]/FACE_GRAD_WEIGHTS:1/FACE_GRAD_WEIGHTS:2}"             # third order, faded faces
V[d4fb]="${V[d4]/FACE_GRAD_WEIGHTS:1/FACE_GRAD_WEIGHTS:2}"             # fourth order, faded faces
V[d4f1]="${V[d4]/FACE_GRAD_WEIGHTS:1/}"                                # fourth order, first-order faces (f1q23r4 route)
[ -f "$HERE/variants.sh" ] && source "$HERE/variants.sh"
for t in "$@"; do
    [ -n "${V[$t]+x}" ] || { echo "unknown variant $t" >&2; exit 1; }
    B="$HERE/build_$t"; rm -rf "$B"; mkdir -p "$B"
    cp -R "$DEMO/Makefile" "$DEMO/include" "$DEMO/src" "$B/"
    [ -d "$DEMO/tests" ] && cp -R "$DEMO/tests" "$B/"
    rm -f "$B/include/parameter.h.cmp_backup"; [ -e "$B/holib" ] || ln -s "$DEMO/holib" "$B/holib"
    P="$B/include/parameter.h"
    C="$COMMON"; VV="${V[$t]}"; [[ $VV == NOCOMMON* ]] && { C=""; VV="${VV#NOCOMMON}"; }
    for kv in $C $VV; do set_flag "$P" "${kv%:*}" "${kv#*:}"; done
    { echo "snapshot of demonstrator/ $(date '+%F %T')"; git -C "$DEMO" rev-parse HEAD
      git -C "$DEMO" status --short -- . ; echo; echo "variant $t: ${V[$t]}"; echo "parameter.h diff vs demonstrator:"
      diff "$DEMO/include/parameter.h" "$P" || true; git -C "$DEMO" diff -- . ; } > "$B/SNAPSHOT.txt"
    (cd "$B" && make -j4 > make.log 2>&1) || { grep -E "error" "$B/make.log" | head; exit 1; }
    [ -x "$B/bin/mlh" ] && echo "built $t"
done
