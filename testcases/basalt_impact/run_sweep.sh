#!/bin/bash
# Basalt-impact failure-model sweep.
# Builds one binary per compile-time flag combo (simple plasticity -> fragmentation
# -> Collins), each into its own result folder, then runs them all in parallel on
# the SAME IC. Usage: ./run_sweep.sh [timeEnd]   (default 8.0)
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

cp "$PARAM" /tmp/param.master.bak
cd "$DEM"
for i in "${!NAMES[@]}"; do
    name=${NAMES[$i]}; flags=${FLAGS[$i]}
    cp /tmp/param.master.bak "$PARAM"
    for f in $flags; do sed -i "s/#define $f 0/#define $f 1/" "$PARAM"; done
    # damage runs: skip the startup transient so it doesn't latch every flaw
    if [[ "$flags" == *FRAGMENTATION* ]]; then
        sed -i 's/#define DAMAGE_START_TIME 0.0/#define DAMAGE_START_TIME 0.6/' "$PARAM"
    fi
    echo "=== building $name [${flags:-elastic baseline}] ==="
    if ! make >/dev/null 2>"$BASE/$name.build.err"; then
        echo "BUILD FAILED for $name"; tail -20 "$BASE/$name.build.err"; cp /tmp/param.master.bak "$PARAM"; exit 1
    fi
    mkdir -p "$BASE/$name/output"
    cp bin/mlh "$BASE/$name/mlh"
    sed -e "s#__OUTDIR__#$BASE/$name/output#" -e "s/__TIMEEND__/$TIMEEND/" \
        "$BASE/config_base.info" > "$BASE/$name/config.info"
    rm -f "$BASE/$name.build.err"
done
cp /tmp/param.master.bak "$PARAM"
make >/dev/null 2>&1   # restore default (all-off) binary

echo "=== all builds done; launching $((${#NAMES[@]})) runs in parallel (timeEnd=$TIMEEND) ==="
for name in "${NAMES[@]}"; do
    ( cd "$BASE/$name" && ./mlh --config config.info >run.log 2>&1; \
      echo "$name finished rc=$?" >> "$BASE/sweep.status" ) &
done
wait
echo "=== ALL RUNS DONE ==="
cat "$BASE/sweep.status"
