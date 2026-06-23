#!/usr/bin/env bash
# Render combined (rho,P,u,Sxx,Sxy,Syy) frames for each failure-model variant and
# encode one movie per variant into basalt_impact/movies/.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PLOTTER="$ROOT/../colliding_rings/densityPlotter.py"
MOVIEDIR="$ROOT/movies"; mkdir -p "$MOVIEDIR"

for d in "$ROOT"/0*/; do
    series="$(basename "$d")"; out="$d/output"
    ls "$out"/*.h5 >/dev/null 2>&1 || { echo "skip $series"; continue; }
    echo "=== [$series] rendering combined frames ==="
    python3 "$PLOTTER" -d "$out" -o "$out" --combined --openBorders -w 1
    echo "=== [$series] encoding movie ==="
    ffmpeg -y -framerate 12 -pattern_type glob -i "$out/*.png" \
        -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" -c:v libx264 -pix_fmt yuv420p \
        "$MOVIEDIR/$series.mp4" </dev/null
done
echo "=== done ==="; ls -la "$MOVIEDIR"
