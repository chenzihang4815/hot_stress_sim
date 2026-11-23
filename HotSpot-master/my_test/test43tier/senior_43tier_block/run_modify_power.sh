#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
HOTSPOT="/home/czh/HotSpot-master/hotspot"
SCRIPTS="/home/czh/HotSpot-master/scripts"
OUTPUT_DIR="/home/czh/HotSpot-master/my_test/output"
POWER_TRACE="$SCRIPT_DIR/senior_43tier_modified_10W.ptrace"

mkdir -p "$OUTPUT_DIR"
cd "$SCRIPT_DIR"

rm -f \
  "$OUTPUT_DIR"/senior_43tier_modified_10W.steady \
  "$OUTPUT_DIR"/senior_43tier_modified_10W.ttrace \
  "$OUTPUT_DIR"/senior_43tier_modified_10W_grid.steady \
  "$OUTPUT_DIR"/senior_43tier_modified_10W_grid.ttrace \
  "$OUTPUT_DIR"/senior_43tier_modified_10W_grid_layer*.grid.steady \
  "$OUTPUT_DIR"/senior_43tier_modified_10W_heatmap.png \
  "$OUTPUT_DIR"/senior_43tier_modified_10W_heatmap.svg

# Block model: one temperature value per tier.
"$HOTSPOT" \
  -f senior_43tier.flp \
  -p "$POWER_TRACE" \
  -c senior_43tier.config \
  -model_type block \
  -o "$OUTPUT_DIR/senior_43tier_modified_10W.ttrace" \
  -steady_file "$OUTPUT_DIR/senior_43tier_modified_10W.steady"

# Grid model: used to produce a 64x64 spatial thermal map.
"$HOTSPOT" \
  -f senior_43tier.flp \
  -p "$POWER_TRACE" \
  -c senior_43tier.config \
  -model_type grid \
  -grid_rows 64 \
  -grid_cols 64 \
  -o "$OUTPUT_DIR/senior_43tier_modified_10W_grid.ttrace" \
  -steady_file "$OUTPUT_DIR/senior_43tier_modified_10W_grid.steady" \
  -grid_steady_file "$OUTPUT_DIR/senior_43tier_modified_10W_grid.steady"

python3 "$SCRIPTS/split_grid_steady.py" \
  "$OUTPUT_DIR/senior_43tier_modified_10W_grid.steady" 4 64 64

python3 "$SCRIPTS/grid_thermal_map.py" \
  senior_43tier.flp \
  "$OUTPUT_DIR/senior_43tier_modified_10W_grid_layer0.grid.steady" \
  64 64 \
  "$OUTPUT_DIR/senior_43tier_modified_10W_heatmap.png"

perl "$SCRIPTS/grid_thermal_map.pl" \
  senior_43tier.flp \
  "$OUTPUT_DIR/senior_43tier_modified_10W_grid_layer0.grid.steady" \
  64 64 \
  > "$OUTPUT_DIR/senior_43tier_modified_10W_heatmap.svg"

printf '\nModified-power results saved in: %s\n' "$OUTPUT_DIR"
