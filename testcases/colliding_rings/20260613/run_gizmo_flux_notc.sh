#!/bin/bash
# GIZMO-faithful elastic flux with tensile correction OFF (ablation): isolates the
# transverse shear-wave HLL dissipation as the ring stabilizer.
set -e
BASE="$(cd "$(dirname "$0")" && pwd)"
DEMO="$(cd "$BASE/../../../demonstrator" && pwd)"
NAME=gizmo_flux_wc2_notc

set_flags () { # KF TC TC1 TC2 TC3 GEF
    sed -i \
        -e "s/^#define KERNEL_FUNCTION .*/#define KERNEL_FUNCTION $1/" \
        -e "s/^#define TENSILE_CORRECTION [01]\$/#define TENSILE_CORRECTION $2/" \
        -e "s/^#define TENSILE_CORRECTION_1 [01]\$/#define TENSILE_CORRECTION_1 $3/" \
        -e "s/^#define TENSILE_CORRECTION_2 [01]\$/#define TENSILE_CORRECTION_2 $4/" \
        -e "s/^#define TENSILE_CORRECTION_3 [01]\$/#define TENSILE_CORRECTION_3 $5/" \
        -e "s/^#define GIZMO_ELASTIC_FLUX [01]\$/#define GIZMO_ELASTIC_FLUX $6/" \
        "$DEMO/include/parameter.h"
}

set_flags 6 0 0 0 0 1
# Force a clean recompile: the Makefile's .d targets use bare 'Particles.o' (not
# 'obj/Particles.o'), so header-only (parameter.h) changes don't reliably trigger
# rebuilds. Removing objects/deps guarantees the binary matches these flags.
rm -f "$DEMO"/obj/*.o "$DEMO"/dep/*.d
make -C "$DEMO" -j"$(nproc)" > "$BASE/$NAME/build.log" 2>&1
cp "$DEMO/bin/mlh" "$DEMO/bin/mlh_$NAME"
cp "$DEMO/include/parameter.h" "$BASE/$NAME/parameter.h"
set_flags 3 1 1 0 0 0  # restore defaults

cd "$DEMO"
"bin/mlh_$NAME" -v -c "$BASE/$NAME/config.info" > "$BASE/$NAME/log.txt" 2>&1
echo "$?" > "$BASE/$NAME/exit_code"
echo "$NAME done: exit=$(cat "$BASE/$NAME/exit_code") snapshots=$(ls "$BASE/$NAME/output" | wc -l)"
