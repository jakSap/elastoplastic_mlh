#!/bin/bash
# Tensile-correction sweep for the colliding rubber rings (IC: 20260611 Hopkins glass,
# config: 20260611/run_1). Phase 1 builds one binary per flag combination (compile-time
# flags in parameter.h), phase 2 runs all variants in parallel from demonstrator/.
set -e

BASE="$(cd "$(dirname "$0")" && pwd)"
DEMO="$(cd "$BASE/../../../demonstrator" && pwd)"

# name  KERNEL_FUNCTION  TENSILE_CORRECTION  TC_1  TC_2  TC_3
VARIANTS="\
cubic_tc_none 3 0 0 0 0
cubic_tc1     3 1 1 0 0
cubic_tc2     3 1 0 1 0
cubic_tc3     3 1 0 0 1
cubic_tc1_tc2 3 1 1 1 0
cubic_tc1_tc3 3 1 1 0 1
wc2_tc_none   6 0 0 0 0
wc2_tc1       6 1 1 0 0
wc2_tc2       6 1 0 1 0
wc2_tc3       6 1 0 0 1
wc2_tc1_tc2   6 1 1 1 0
wc2_tc1_tc3   6 1 1 0 1
wc4_tc1_tc3   7 1 1 0 1"

set_flags () { # k tc tc1 tc2 tc3
    sed -i \
        -e "s/^#define KERNEL_FUNCTION .*/#define KERNEL_FUNCTION $1/" \
        -e "s/^#define TENSILE_CORRECTION [01]\$/#define TENSILE_CORRECTION $2/" \
        -e "s/^#define TENSILE_CORRECTION_1 [01]\$/#define TENSILE_CORRECTION_1 $3/" \
        -e "s/^#define TENSILE_CORRECTION_2 [01]\$/#define TENSILE_CORRECTION_2 $4/" \
        -e "s/^#define TENSILE_CORRECTION_3 [01]\$/#define TENSILE_CORRECTION_3 $5/" \
        "$DEMO/include/parameter.h"
}

echo "=== Phase 1: building variant binaries ==="
while read -r name k tc tc1 tc2 tc3; do
    [ -z "$name" ] && continue
    echo "  > $name (KERNEL_FUNCTION=$k TC=$tc TC1=$tc1 TC2=$tc2 TC3=$tc3)"
    mkdir -p "$BASE/$name/output"
    set_flags "$k" "$tc" "$tc1" "$tc2" "$tc3"
    make -C "$DEMO" -j"$(nproc)" > "$BASE/$name/build.log" 2>&1
    cp "$DEMO/bin/mlh" "$DEMO/bin/mlh_$name"
    cp "$DEMO/include/parameter.h" "$BASE/$name/parameter.h"
done <<< "$VARIANTS"

# restore repository defaults
set_flags 3 1 1 0 0

echo "=== Phase 2: running all variants in parallel ==="
cd "$DEMO"
while read -r name _; do
    [ -z "$name" ] && continue
    (
        "bin/mlh_$name" -v -c "$BASE/$name/config.info" \
            > "$BASE/$name/log.txt" 2>&1
        echo "$?" > "$BASE/$name/exit_code"
    ) &
done <<< "$VARIANTS"
wait

echo "=== Summary ==="
while read -r name _; do
    [ -z "$name" ] && continue
    rc="$(cat "$BASE/$name/exit_code" 2>/dev/null || echo '?')"
    nsnap="$(ls "$BASE/$name/output" 2>/dev/null | wc -l)"
    echo "  $name: exit=$rc snapshots=$nsnap"
done <<< "$VARIANTS" | tee "$BASE/summary.txt"
