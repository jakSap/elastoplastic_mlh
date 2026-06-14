#!/bin/bash
# Build + run the one missing combo wc4_tc1_tc2 (KERNEL_FUNCTION=7, TC1+TC2).
set -e
BASE="$(cd "$(dirname "$0")" && pwd)"
DEMO="$(cd "$BASE/../../../demonstrator" && pwd)"
NAME=wc4_tc1_tc2

sed -i \
    -e "s/^#define KERNEL_FUNCTION .*/#define KERNEL_FUNCTION 7/" \
    -e "s/^#define TENSILE_CORRECTION [01]\$/#define TENSILE_CORRECTION 1/" \
    -e "s/^#define TENSILE_CORRECTION_1 [01]\$/#define TENSILE_CORRECTION_1 1/" \
    -e "s/^#define TENSILE_CORRECTION_2 [01]\$/#define TENSILE_CORRECTION_2 1/" \
    -e "s/^#define TENSILE_CORRECTION_3 [01]\$/#define TENSILE_CORRECTION_3 0/" \
    "$DEMO/include/parameter.h"

make -C "$DEMO" -j"$(nproc)" > "$BASE/$NAME/build.log" 2>&1
cp "$DEMO/bin/mlh" "$DEMO/bin/mlh_$NAME"
cp "$DEMO/include/parameter.h" "$BASE/$NAME/parameter.h"

# restore repository defaults
sed -i \
    -e "s/^#define KERNEL_FUNCTION .*/#define KERNEL_FUNCTION 3/" \
    -e "s/^#define TENSILE_CORRECTION [01]\$/#define TENSILE_CORRECTION 1/" \
    -e "s/^#define TENSILE_CORRECTION_1 [01]\$/#define TENSILE_CORRECTION_1 1/" \
    -e "s/^#define TENSILE_CORRECTION_2 [01]\$/#define TENSILE_CORRECTION_2 0/" \
    -e "s/^#define TENSILE_CORRECTION_3 [01]\$/#define TENSILE_CORRECTION_3 0/" \
    "$DEMO/include/parameter.h"

cd "$DEMO"
"bin/mlh_$NAME" -v -c "$BASE/$NAME/config.info" > "$BASE/$NAME/log.txt" 2>&1
echo "$?" > "$BASE/$NAME/exit_code"
echo "wc4_tc1_tc2 done: exit=$(cat "$BASE/$NAME/exit_code") snapshots=$(ls "$BASE/$NAME/output" | wc -l)"
