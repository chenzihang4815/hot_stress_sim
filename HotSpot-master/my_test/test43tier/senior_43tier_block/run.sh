#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
HOTSPOT="/home/czh/HotSpot-master/hotspot"
SCRIPTS="/home/czh/HotSpot-master/scripts"
OUTPUT_DIR="/home/czh/HotSpot-master/my_test/output"

mkdir -p "$OUTPUT_DIR"
cd "$SCRIPT_DIR"

rm -f \
  "$OUTPUT_DIR"/senior_43tier.steady \
  "$OUTPUT_DIR"/senior_43tier.ttrace \
  "$OUTPUT_DIR"/senior_43tier_grid.steady \
  "$OUTPUT_DIR"/senior_43tier_grid.ttrace \
  "$OUTPUT_DIR"/senior_43tier_layer*.grid.steady \
  "$OUTPUT_DIR"/senior_43tier_heatmap.png \
  "$OUTPUT_DIR"/senior_43tier_heatmap.svg \
  "$OUTPUT_DIR"/senior_43tier_floorplan.png

# Block model: one temperature value per tier.
"$HOTSPOT" \
  -f senior_43tier.flp \
  -p senior_43tier.ptrace \
  -c senior_43tier.config \
  -model_type block \
  -o "$OUTPUT_DIR/senior_43tier.ttrace" \
  -steady_file "$OUTPUT_DIR/senior_43tier.steady"

# Grid model: used to produce a 64x64 spatial thermal map.
"$HOTSPOT" \
  -f senior_43tier.flp \
  -p senior_43tier.ptrace \
  -c senior_43tier.config \
  -model_type grid \
  -grid_rows 64 \
  -grid_cols 64 \
  -o "$OUTPUT_DIR/senior_43tier_grid.ttrace" \
  -steady_file "$OUTPUT_DIR/senior_43tier_grid.steady" \
  -grid_steady_file "$OUTPUT_DIR/senior_43tier_grid.steady"

# Split the combined grid output into chip, TIM, spreader, and heatsink layers.
python3 "$SCRIPTS/split_grid_steady.py" \
  "$OUTPUT_DIR/senior_43tier_grid.steady" 4 64 64

# Generate a floorplan image and a chip-layer thermal map.
python3 "$SCRIPTS/visualize_floorplan.py" \
  senior_43tier.flp \
  "$OUTPUT_DIR/senior_43tier_floorplan.png"
python3 "$SCRIPTS/grid_thermal_map.py" \
  senior_43tier.flp \
  "$OUTPUT_DIR/senior_43tier_grid_layer0.grid.steady" \
  64 64 \
  "$OUTPUT_DIR/senior_43tier_heatmap.png"
perl "$SCRIPTS/grid_thermal_map.pl" \
  senior_43tier.flp \
  "$OUTPUT_DIR/senior_43tier_grid_layer0.grid.steady" \
  64 64 \
  > "$OUTPUT_DIR/senior_43tier_heatmap.svg"

printf '\nResults saved in: %s\n' "$OUTPUT_DIR"
