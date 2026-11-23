#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
HOTSPOT="/home/czh/HotSpot-master/hotspot"
SCRIPTS="/home/czh/HotSpot-master/scripts"
OUTPUT_DIR="$SCRIPT_DIR/output"
POWER_TRACE="$SCRIPT_DIR/test25tier_original.ptrace"
PREFIX="test25tier_original"

mkdir -p "$OUTPUT_DIR"
cd "$SCRIPT_DIR"

rm -f \
  "$OUTPUT_DIR"/"$PREFIX".steady \
  "$OUTPUT_DIR"/"$PREFIX".ttrace \
  "$OUTPUT_DIR"/"${PREFIX}_grid.steady" \
  "$OUTPUT_DIR"/"${PREFIX}_grid.ttrace" \
  "$OUTPUT_DIR"/"${PREFIX}_grid_layer"*.grid.steady \
  "$OUTPUT_DIR"/"${PREFIX}_heatmap.png" \
  "$OUTPUT_DIR"/"${PREFIX}_heatmap.svg" \
  "$OUTPUT_DIR"/test25tier_floorplan.png

# Block model: one temperature value per chiplet.
"$HOTSPOT" \
  -f test25tier.flp \
  -p "$POWER_TRACE" \
  -c test25tier.config \
  -model_type block \
  -o "$OUTPUT_DIR/${PREFIX}.ttrace" \
  -steady_file "$OUTPUT_DIR/${PREFIX}.steady"

# Grid model: used to produce a 64x64 spatial thermal map.
"$HOTSPOT" \
  -f test25tier.flp \
  -p "$POWER_TRACE" \
  -c test25tier.config \
  -model_type grid \
  -grid_rows 64 \
  -grid_cols 64 \
  -o "$OUTPUT_DIR/${PREFIX}_grid.ttrace" \
  -steady_file "$OUTPUT_DIR/${PREFIX}_grid.steady" \
  -grid_steady_file "$OUTPUT_DIR/${PREFIX}_grid.steady"

python3 "$SCRIPTS/split_grid_steady.py" \
  "$OUTPUT_DIR/${PREFIX}_grid.steady" 4 64 64

python3 "$SCRIPTS/visualize_floorplan.py" \
  test25tier.flp \
  "$OUTPUT_DIR/test25tier_floorplan.png"
python3 "$SCRIPTS/grid_thermal_map.py" \
  test25tier.flp \
  "$OUTPUT_DIR/${PREFIX}_grid_layer0.grid.steady" \
  64 64 \
  "$OUTPUT_DIR/${PREFIX}_heatmap.png"
perl "$SCRIPTS/grid_thermal_map.pl" \
  test25tier.flp \
  "$OUTPUT_DIR/${PREFIX}_grid_layer0.grid.steady" \
  64 64 \
  > "$OUTPUT_DIR/${PREFIX}_heatmap.svg"

printf '\nOriginal-power results saved in: %s\n' "$OUTPUT_DIR"
