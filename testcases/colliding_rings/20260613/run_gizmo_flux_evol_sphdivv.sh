#!/bin/bash
# GIZMO-faithful flux + explicit volume integration driven by GIZMO's SPH
# kernel-weighted Particle_DivVel (EXPLICIT_VOL_SPH_DIVV=1). Tests whether the
# robust surface divergence removes the bloating/drift seen with the MFM-grad trace.
set -e
BASE="$(cd "$(dirname "$0")" && pwd)"
DEMO="$(cd "$BASE/../../../demonstrator" && pwd)"
NAME=gizmo_flux_wc2_evol_sphdivv

set_flags () { # KF GEF EVOL SPHDIVV
    sed -i \
        -e "s/^#define KERNEL_FUNCTION .*/#define KERNEL_FUNCTION $1/" \
        -e "s/^#define TENSILE_CORRECTION [01]\$/#define TENSILE_CORRECTION 1/" \
        -e "s/^#define TENSILE_CORRECTION_1 [01]\$/#define TENSILE_CORRECTION_1 1/" \
        -e "s/^#define TENSILE_CORRECTION_2 [01]\$/#define TENSILE_CORRECTION_2 0/" \
        -e "s/^#define TENSILE_CORRECTION_3 [01]\$/#define TENSILE_CORRECTION_3 0/" \
        -e "s/^#define GIZMO_ELASTIC_FLUX [01]\$/#define GIZMO_ELASTIC_FLUX $2/" \
        -e "s/^#define EXPLICIT_VOL_INTEGRATION [01]\$/#define EXPLICIT_VOL_INTEGRATION $3/" \
        -e "s/^#define EXPLICIT_VOL_SPH_DIVV [01]\$/#define EXPLICIT_VOL_SPH_DIVV $4/" \
        "$DEMO/include/parameter.h"
}

set_flags 6 1 1 1
rm -f "$DEMO"/obj/*.o "$DEMO"/dep/*.d
make -C "$DEMO" -j"$(nproc)" > "$BASE/$NAME/build.log" 2>&1
cp "$DEMO/bin/mlh" "$DEMO/bin/mlh_$NAME"
cp "$DEMO/include/parameter.h" "$BASE/$NAME/parameter.h"
set_flags 3 0 0 1  # restore defaults (KERNEL=3, GEF=0, EVOL=0, SPHDIVV default 1)

cd "$DEMO"
"bin/mlh_$NAME" -v -c "$BASE/$NAME/config.info" > "$BASE/$NAME/log.txt" 2>&1
echo "$?" > "$BASE/$NAME/exit_code"
echo "$NAME done: exit=$(cat "$BASE/$NAME/exit_code") snapshots=$(ls "$BASE/$NAME/output" | wc -l)"
