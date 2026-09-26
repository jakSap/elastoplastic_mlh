#!/usr/bin/env bash
# mkbuild.sh <tag> [FLAG:VAL ...]  -> build_<tag>/ = snapshot of demonstrator/ + only the listed flags
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DEMO="$(cd "$HERE/../../../demonstrator" && pwd)"
t=$1; shift; B="$HERE/build_$t"
set_flag() { grep -qE "^#define $2( |$)" "$1" || { echo "missing $2" >&2; exit 1; }
             sed -i -E "s|^#define $2( .*)?$|#define $2 $3|" "$1"; }
rm -rf "$B"; mkdir -p "$B"; cp -R "$DEMO/Makefile" "$DEMO/include" "$DEMO/src" "$DEMO/tests" "$B/"
ln -s "$DEMO/holib" "$B/holib"
for kv in "$@"; do set_flag "$B/include/parameter.h" "${kv%%:*}" "${kv#*:}"; done
{ echo "build_$t  $(date '+%F %T')  flags: $*"; git -C "$DEMO" rev-parse HEAD; git -C "$DEMO" status --short -- .
  echo "md5(src+include): $(cat "$B"/src/* "$B"/include/* | md5sum | cut -c1-12)"
  echo; echo "diff vs demonstrator/:"; diff -r "$DEMO/include" "$B/include" || true; diff -r "$DEMO/src" "$B/src" || true; } > "$B/SNAPSHOT.txt"
(cd "$B" && make -j6 > make.log 2>&1) && [ -x "$B/bin/mlh" ] && echo "built $t" || { echo "FAILED $t"; grep -m5 -i error "$B/make.log"; exit 1; }
