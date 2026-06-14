#!/usr/bin/env bash
# Render density plots (--openBorders) for every series under 20260613/ and
# combine each series into one movie, collected in 20260613/movies/.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PLOTTER="$ROOT/../densityPlotter.py"
MOVIEDIR="$ROOT/movies"
mkdir -p "$MOVIEDIR"

for d in "$ROOT"/*/; do
    series="$(basename "$d")"
    out="$d/output"
    [ -d "$out" ] || continue
    ls "$out"/*.h5 >/dev/null 2>&1 || continue

    echo "=== [$series] rendering frames ==="
    python3 "$PLOTTER" -d "$out" -o "$out" --openBorders

    echo "=== [$series] encoding movie ==="
    ffmpeg -y -framerate 30 -pattern_type glob -i "$out/*.png" \
        -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" \
        -c:v libx264 -pix_fmt yuv420p \
        "$MOVIEDIR/$series.mp4" </dev/null
done

echo "=== all done; movies in $MOVIEDIR ==="
ls -la "$MOVIEDIR"
