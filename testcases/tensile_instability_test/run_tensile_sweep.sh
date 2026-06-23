#!/bin/bash
# Parameter sweep for the tensile-instability stretch testcase.
#
# Grid = strain-rate (IC axis) x compile-flag combo (kernel x BC mode).
#   - strain-rate changes only the IC h5 -> generated once per rate (Phase 0).
#   - kernel + BC mode are compile-time flags in parameter.h -> one binary per
#     combo (Phase 1). Phase 2 runs combo x strain-rate; Phase 3 summarizes.
#
# Fixed for all runs: MFM (MESHLESS_FINITE_MASS 1), Tillotson (EOS 2), ELASTIC 1,
# adaptive SML (VARIABLE_SML 1) with NNNTarget 32, no tensile correction,
# delta_p 0.05, timeStep 0.01, h5DumpInterval 10, timeEnd 10.
#
# BC modes:
#   periodic : PERIODIC_BOUNDARIES 1, SURFACE_VOLCORR 0, EXPLICIT_VOL_INTEGRATION 0
#   free     : PERIODIC_BOUNDARIES 0, SURFACE_VOLCORR 1, EXPLICIT_VOL_INTEGRATION 1
#     (EXPLICIT_VOL_INTEGRATION/SURFACE_VOLCORR error under PERIODIC and require
#      ELASTIC 1 -- see guards in demonstrator/include/parameter.h. PERIODIC is a
#      global flag, so the free mode is non-periodic in y too: free top/bottom.)
#
# Extending this script:
#   - new strain-rate           -> add to STRAINRATES
#   - new compile-flag axis      -> add a column to COMBOS + one line in set_flags
#   - new test problem           -> change TESTNAME, the Phase-0 IC command, and
#                                   config_base.info placeholders; phases 1-3 stay.
#
# Usage: ./run_tensile_sweep.sh [timeEnd]   (default 10)
set -e

# --- configurable head -----------------------------------------------------
TESTNAME=tensile_instability_test
BASE="$(cd "$(dirname "$0")" && pwd)"
DEMO="$(cd "$BASE/../../demonstrator" && pwd)"
PARAM="$DEMO/include/parameter.h"
DATE="$(date +%Y%m%d)"
OUTROOT="$BASE/$DATE"

STRAINRATES=(0.01 0.1 0.5 1)        # IC axis: vx = a*x  (--strain-rate)
DELTAP=0.05                          # IC resolution
TIMEEND=${1:-10}
BUFFER_SAFETY=1.3                    # x-box half-width = ceil(10*(1+a*tEnd)*safety)

# compile-flag combos: "name KERNEL_FUNCTION PERIODIC SURFACE_VOLCORR EXPLICIT_VOL"
COMBOS=(
  "periodic_cubic 3 1 0 0"
  "periodic_wc2   6 1 0 0"
  "free_cubic     3 0 1 1"
  "free_wc2       6 0 1 1"
)

# --- helpers ---------------------------------------------------------------
# Set every swept/fixed compile-time flag for one combo. [01]$ anchors keep
# TENSILE_CORRECTION/SURFACE_VOLCORR from matching their _1/_A... siblings.
set_flags () { # KERNEL PERIODIC SURFVOLCORR EXPLICITVOL
    sed -i \
        -e 's/^#define EOS .*/#define EOS 2/' \
        -e 's/^#define ELASTIC .*/#define ELASTIC 1/' \
        -e 's/^#define MESHLESS_FINITE_MASS .*/#define MESHLESS_FINITE_MASS 1/' \
        -e 's/^#define VARIABLE_SML .*/#define VARIABLE_SML 1/' \
        -e 's/^#define TENSILE_CORRECTION [01]$/#define TENSILE_CORRECTION 0/' \
        -e "s/^#define KERNEL_FUNCTION .*/#define KERNEL_FUNCTION $1/" \
        -e "s/^#define PERIODIC_BOUNDARIES .*/#define PERIODIC_BOUNDARIES $2/" \
        -e "s/^#define SURFACE_VOLCORR [01]$/#define SURFACE_VOLCORR $3/" \
        -e "s/^#define EXPLICIT_VOL_INTEGRATION [01]$/#define EXPLICIT_VOL_INTEGRATION $4/" \
        "$PARAM"
}

xbox () { # strain-rate -> ceil(10*(1+a*tEnd)*safety)
    awk -v a="$1" -v t="$TIMEEND" -v s="$BUFFER_SAFETY" \
        'BEGIN{x=10*(1+a*t)*s; printf "%d", (x==int(x))?x:int(x)+1}'
}

mkdir -p "$OUTROOT"

# --- Phase 0: generate ICs (one per strain-rate, idempotent) ---------------
echo "=== Phase 0: generating ICs ==="
for a in "${STRAINRATES[@]}"; do
    ic="$BASE/tensile_stretch_ic_sr${a}_dp${DELTAP}.h5"
    if [ -f "$ic" ]; then
        echo "  > sr=$a: $(basename "$ic") exists, skipping"
    else
        echo "  > sr=$a: generating $(basename "$ic")"
        python3 "$BASE/generateIC_tensile_stretch.py" --strain-rate "$a" --delta-p "$DELTAP"
    fi
done

# --- Phase 1: build one binary per compile-flag combo ----------------------
cp "$PARAM" /tmp/param.master.bak
trap 'cp /tmp/param.master.bak "$PARAM"' EXIT   # always restore repo defaults
cd "$DEMO"
echo "=== Phase 1: building ${#COMBOS[@]} combo binaries ==="
for combo in "${COMBOS[@]}"; do
    read -r name k p s e <<< "$combo"
    echo "  > $name (KERNEL=$k PERIODIC=$p SURFVOLCORR=$s EXPLICITVOL=$e)"
    cp /tmp/param.master.bak "$PARAM"
    set_flags "$k" "$p" "$s" "$e"
    rm -f "$DEMO"/obj/*.o "$DEMO"/dep/*.d   # parameter.h is a header; .d targets unreliable
    if ! make -j"$(nproc)" > "$OUTROOT/$name.build.log" 2>&1; then
        echo "BUILD FAILED for $name"; tail -20 "$OUTROOT/$name.build.log"; exit 1
    fi
    cp "$DEMO/bin/mlh" "$DEMO/bin/mlh_$name"
    cp "$PARAM" "$OUTROOT/$name.parameter.h"
done
cp /tmp/param.master.bak "$PARAM"   # restore repo defaults (also via trap)

# --- Phase 2: run combo x strain-rate in parallel --------------------------
echo "=== Phase 2: launching $(( ${#COMBOS[@]} * ${#STRAINRATES[@]} )) runs (timeEnd=$TIMEEND) ==="
for combo in "${COMBOS[@]}"; do
    read -r name _ <<< "$combo"
    for a in "${STRAINRATES[@]}"; do
        run="${name}_sr${a}"
        rundir="$OUTROOT/$run"
        mkdir -p "$rundir/output"; rm -f "$rundir/output"/*.h5
        ux="$(xbox "$a")"
        ic="$BASE/tensile_stretch_ic_sr${a}_dp${DELTAP}.h5"
        sed -e "s#__INITFILE__#$ic#" \
            -e "s#__OUTDIR__#$rundir/output#" \
            -e "s/__TIMEEND__/$TIMEEND/" \
            -e "s/__LOWERX__/-$ux/" \
            -e "s/__UPPERX__/$ux/" \
            "$BASE/config_base.info" > "$rundir/config.info"
        (
            set +e   # capture mlh's real exit code instead of aborting the subshell
            cd "$DEMO"
            "bin/mlh_$name" -c "$rundir/config.info" > "$rundir/log.txt" 2>&1
            echo "$?" > "$rundir/exit_code"
        ) &
    done
done
wait

# --- Phase 3: summary ------------------------------------------------------
echo "=== Summary ==="
{
    for combo in "${COMBOS[@]}"; do
        read -r name _ <<< "$combo"
        for a in "${STRAINRATES[@]}"; do
            run="${name}_sr${a}"
            rc="$(cat "$OUTROOT/$run/exit_code" 2>/dev/null || echo '?')"
            nsnap="$(ls "$OUTROOT/$run/output" 2>/dev/null | wc -l)"
            printf "  %-22s exit=%s snapshots=%s\n" "$run" "$rc" "$nsnap"
        done
    done
} | tee "$OUTROOT/summary.txt"
echo "=== ALL DONE ($OUTROOT) ==="
