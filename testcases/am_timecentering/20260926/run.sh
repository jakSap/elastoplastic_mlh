#!/usr/bin/env bash
# run.sh <build_tag> <ic.h5 (in ics/)> <run_tag> <t_end> <dump_dt> [spacing]
# -> runs/<run_tag>/{config.info,parameter.h,BUILD,output/,run.log,exit_code}
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; cd "$HERE"
b=$1; ic=$2; tag=$3; tend=$4; ddt=$5
d=$(python3 -c "import h5py;print(h5py.File('ics/$ic').attrs['spacing'])")
R=runs/$tag; rm -rf "$R"; mkdir -p "$R/output"
cp build_$b/include/parameter.h "$R/"; echo "$b" > "$R/BUILD"
ks=$(python3 -c "print(0.07978845608028655/0.025*$d)")
sed -e "s|^initFile .*|initFile ics/$ic|" -e "s|^outDir .*|outDir $R/output|" -e "s|^timeStep .*|timeStep $ddt|" \
    -e "s|^timeEnd .*|timeEnd $tend|" -e "s|^kernelSize .*|kernelSize $ks|" -e "1s|.*|; am_timecentering $tag (build_$b, ic $ic)|" \
    config_template.info > "$R/config.info"
start=$(date +%s)
./build_$b/bin/mlh -c "$R/config.info" > "$R/run.log" 2>&1; ec=$?
echo $ec > "$R/exit_code"
echo "$tag build=$b ic=$ic exit=$ec wall=$(( $(date +%s)-start ))s nsnap=$(ls $R/output | wc -l) negE=$(grep -c 'Negative internal' $R/run.log)" >> summary.txt
