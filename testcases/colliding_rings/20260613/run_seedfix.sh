#!/bin/bash
set -e
BASE="$(cd "$(dirname "$0")" && pwd)"; DEMO="$(cd "$BASE/../../../demonstrator" && pwd)"; NAME=gizmo_flux_wc2_evol_seedfix
sed -i -e 's/^#define KERNEL_FUNCTION .*/#define KERNEL_FUNCTION 6/' \
  -e 's/^#define TENSILE_CORRECTION [01]$/#define TENSILE_CORRECTION 1/' \
  -e 's/^#define TENSILE_CORRECTION_1 [01]$/#define TENSILE_CORRECTION_1 1/' \
  -e 's/^#define TENSILE_CORRECTION_2 [01]$/#define TENSILE_CORRECTION_2 0/' \
  -e 's/^#define TENSILE_CORRECTION_3 [01]$/#define TENSILE_CORRECTION_3 0/' \
  -e 's/^#define GIZMO_ELASTIC_FLUX [01]$/#define GIZMO_ELASTIC_FLUX 1/' \
  -e 's/^#define EXPLICIT_VOL_INTEGRATION [01]$/#define EXPLICIT_VOL_INTEGRATION 1/' \
  -e 's/^#define EXPLICIT_VOL_SPH_DIVV [01]$/#define EXPLICIT_VOL_SPH_DIVV 0/' \
  -e 's/^#define EXPLICIT_VOL_SEED_SKIP [0-9]*$/#define EXPLICIT_VOL_SEED_SKIP 1/' \
  "$DEMO/include/parameter.h"
rm -f "$DEMO"/obj/*.o "$DEMO"/dep/*.d
make -C "$DEMO" -j"$(nproc)" > "$BASE/$NAME/build.log" 2>&1
cp "$DEMO/bin/mlh" "$DEMO/bin/mlh_$NAME"; cp "$DEMO/include/parameter.h" "$BASE/$NAME/parameter.h"
sed -i -e 's/^#define KERNEL_FUNCTION .*/#define KERNEL_FUNCTION 3/' \
  -e 's/^#define GIZMO_ELASTIC_FLUX [01]$/#define GIZMO_ELASTIC_FLUX 0/' \
  -e 's/^#define EXPLICIT_VOL_INTEGRATION [01]$/#define EXPLICIT_VOL_INTEGRATION 0/' \
  -e 's/^#define EXPLICIT_VOL_SPH_DIVV [01]$/#define EXPLICIT_VOL_SPH_DIVV 1/' "$DEMO/include/parameter.h"
cd "$DEMO"; "bin/mlh_$NAME" -v -c "$BASE/$NAME/config.info" > "$BASE/$NAME/log.txt" 2>&1
echo "$?" > "$BASE/$NAME/exit_code"
echo "$NAME done: exit=$(cat "$BASE/$NAME/exit_code") snapshots=$(ls "$BASE/$NAME/output"|wc -l)"
