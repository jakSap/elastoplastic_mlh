#!/usr/bin/env bash
# Render density frames + encode a movie for the named variant(s) only.
# Usage: make_movies_new.sh <variant_dir> [<variant_dir> ...]
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PLOTTER="$ROOT/../densityPlotter.py"
MOVIEDIR="$ROOT/movies"
mkdir -p "$MOVIEDIR"

for series in "$@"; do
    out="$ROOT/$series/output"
    [ -d "$out" ] || { echo "skip $series (no output dir)"; continue; }
    ls "$out"/*.h5 >/dev/null 2>&1 || { echo "skip $series (no h5)"; continue; }
    echo "=== [$series] rendering frames ==="
    python3 "$PLOTTER" -d "$out" -o "$out" --openBorders
    echo "=== [$series] encoding movie ==="
    ffmpeg -y -framerate 30 -pattern_type glob -i "$out/*.png" \
        -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" \
        -c:v libx264 -pix_fmt yuv420p \
        "$MOVIEDIR/$series.mp4" </dev/null
    echo "=== [$series] -> $MOVIEDIR/$series.mp4 ==="
done
ls -la "$MOVIEDIR"
