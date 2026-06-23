#!/bin/bash
# Basalt-impact failure-model sweep, rerun WITH the GIZMO-faithful elastic flux +
# explicit-volume seed fix on top of each failure model:
#   GIZMO_ELASTIC_FLUX=1  EXPLICIT_VOL_INTEGRATION=1  EXPLICIT_VOL_SEED_SKIP=1
# The stress flux is damage-aware (1-D under DAMAGE_ACTS_ON_S). Reruns in-place
# into the 0X_* dirs (baseline summary PNGs were copied to *_baseline.png first).
# Usage: ./run_sweep_gizmoflux.sh [timeEnd]   (default 8.0)
set -e
ROOT=/home/sappler/Documents/newMFM
DEM=$ROOT/demonstrator
BASE=$ROOT/testcases/basalt_impact
PARAM=$DEM/include/parameter.h
TIMEEND=${1:-8.0}

NAMES=(01_elastic 02_vonmises 03_mohrcoulomb 04_druckerprager \
       05_fragmentation 06_fragmentation_damageS 07_collins_damage)
FLAGS=(""                                                   \
       "VON_MISES_PLASTICITY"                               \
       "MOHR_COULOMB_PLASTICITY"                            \
       "DRUCKER_PRAGER_PLASTICITY"                          \
       "FRAGMENTATION"                                      \
       "FRAGMENTATION DAMAGE_ACTS_ON_S"                     \
       "FRAGMENTATION DAMAGE_ACTS_ON_S COLLINS_PLASTICITY")

# new GIZMO-faithful elastic flags applied to every variant
enable_gizmo () {
    sed -i \
        -e 's/^#define GIZMO_ELASTIC_FLUX [01]$/#define GIZMO_ELASTIC_FLUX 1/' \
        -e 's/^#define EXPLICIT_VOL_INTEGRATION [01]$/#define EXPLICIT_VOL_INTEGRATION 1/' \
        -e 's/^#define EXPLICIT_VOL_SEED_SKIP [0-9]*$/#define EXPLICIT_VOL_SEED_SKIP 1/' \
        "$PARAM"
}

cp "$PARAM" /tmp/param.master.bak
cd "$DEM"
> "$BASE/sweep_gizmoflux.status"
for i in "${!NAMES[@]}"; do
    name=${NAMES[$i]}; flags=${FLAGS[$i]}
    cp /tmp/param.master.bak "$PARAM"
    for f in $flags; do sed -i "s/#define $f 0/#define $f 1/" "$PARAM"; done
    if [[ "$flags" == *FRAGMENTATION* ]]; then
        # High-velocity impact: contact is near t=0, so the damage gate only needs to
        # clear the startup transient (the seed fix heals t=0 by step 1), not wait for
        # the (now immediate) contact. Small value keeps damage active through the impact.
        sed -i 's/#define DAMAGE_START_TIME 0.0/#define DAMAGE_START_TIME 0.05/' "$PARAM"
    fi
    enable_gizmo
    echo "=== building $name [${flags:-elastic baseline}] + gizmoflux ==="
    rm -f "$DEM"/obj/*.o "$DEM"/dep/*.d   # Makefile .d targets are unreliable on header-only changes
    if ! make -j"$(nproc)" >/dev/null 2>"$BASE/$name.build.err"; then
        echo "BUILD FAILED for $name"; tail -20 "$BASE/$name.build.err"; cp /tmp/param.master.bak "$PARAM"; exit 1
    fi
    mkdir -p "$BASE/$name/output"; rm -f "$BASE/$name/output"/*.h5
    cp bin/mlh "$BASE/$name/mlh"
    cp "$PARAM" "$BASE/$name/parameter.h"
    sed -e "s#__OUTDIR__#$BASE/$name/output#" -e "s/__TIMEEND__/$TIMEEND/" \
        "$BASE/config_base.info" > "$BASE/$name/config.info"
    rm -f "$BASE/$name.build.err"
done
cp /tmp/param.master.bak "$PARAM"   # restore repo defaults

echo "=== all builds done; launching ${#NAMES[@]} runs in parallel (timeEnd=$TIMEEND) ==="
for name in "${NAMES[@]}"; do
    ( cd "$BASE/$name" && ./mlh --config config.info >run.log 2>&1; \
      echo "$name finished rc=$?" >> "$BASE/sweep_gizmoflux.status" ) &
done
wait
echo "=== ALL RUNS DONE ==="
cat "$BASE/sweep_gizmoflux.status"
