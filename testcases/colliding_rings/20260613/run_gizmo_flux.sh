#!/bin/bash
# Build + run the GIZMO-faithful elastic-flux variant(s).
# Usage: run_gizmo_flux.sh <name> <KERNEL_FUNCTION> <GIZMO_ELASTIC_FLUX>
set -e
BASE="$(cd "$(dirname "$0")" && pwd)"
DEMO="$(cd "$BASE/../../../demonstrator" && pwd)"
NAME="$1"; KF="${2:-6}"; GEF="${3:-1}"

set_flags () { # KERNEL_FUNCTION TC TC1 TC2 TC3 GIZMO_ELASTIC_FLUX
    sed -i \
        -e "s/^#define KERNEL_FUNCTION .*/#define KERNEL_FUNCTION $1/" \
        -e "s/^#define TENSILE_CORRECTION [01]\$/#define TENSILE_CORRECTION $2/" \
        -e "s/^#define TENSILE_CORRECTION_1 [01]\$/#define TENSILE_CORRECTION_1 $3/" \
        -e "s/^#define TENSILE_CORRECTION_2 [01]\$/#define TENSILE_CORRECTION_2 $4/" \
        -e "s/^#define TENSILE_CORRECTION_3 [01]\$/#define TENSILE_CORRECTION_3 $5/" \
        -e "s/^#define GIZMO_ELASTIC_FLUX [01]\$/#define GIZMO_ELASTIC_FLUX $6/" \
        "$DEMO/include/parameter.h"
}

mkdir -p "$BASE/$NAME/output"
set_flags "$KF" 1 1 0 0 "$GEF"
make -C "$DEMO" -j"$(nproc)" > "$BASE/$NAME/build.log" 2>&1
cp "$DEMO/bin/mlh" "$DEMO/bin/mlh_$NAME"
cp "$DEMO/include/parameter.h" "$BASE/$NAME/parameter.h"

# restore repository defaults
set_flags 3 1 1 0 0 0

cd "$DEMO"
"bin/mlh_$NAME" -v -c "$BASE/$NAME/config.info" > "$BASE/$NAME/log.txt" 2>&1
echo "$?" > "$BASE/$NAME/exit_code"
echo "$NAME done: exit=$(cat "$BASE/$NAME/exit_code") snapshots=$(ls "$BASE/$NAME/output" | wc -l)"
