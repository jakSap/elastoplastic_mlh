#!/usr/bin/env bash
# runw.sh <build_tag> <W1_sound|W2_pwave|W3_shear> <N> <run_tag>
# Wave run with the config of testplan_20260926 runs/<test>/ref[_tillotson]_N<N> (IC, box, dt, t_end, dumps), own build.
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; cd "$HERE"
b=$1; test=$2; n=$3; tag=$4
src=../../testplan_20260926/runs/$test/ref_N$n; [ $test = W1_sound ] && src=../../testplan_20260926/runs/$test/ref_tillotson_N$n
[ $test = W1_ideal ] && src=../../testplan_20260926/runs/W1_sound/ref_ideal_N$n
R=runsw/$test/$tag; rm -rf "$R"; mkdir -p "$R/output"
cp build_$b/include/parameter.h "$R/"; echo "$b" > "$R/BUILD"
sed -e "s|^initFile ics/|initFile ics_tp/|" -e "s|^outDir .*|outDir $R/output|" -e "1s|.*|; am_timecentering wave $test $tag (build_$b, config from $src)|" $src/config.info > "$R/config.info"
start=$(date +%s)
./build_$b/bin/mlh -c "$R/config.info" > "$R/run.log" 2>&1; ec=$?
echo $ec > "$R/exit_code"
echo "W $test/$tag build=$b exit=$ec wall=$(( $(date +%s)-start ))s nsnap=$(ls $R/output | wc -l) negE=$(grep -c 'Negative internal' $R/run.log)" >> summary.txt
