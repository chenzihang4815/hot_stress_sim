#!/usr/bin/env bash
set -e
ROOT="$(cd "$(dirname "$0")/../HotSpot-master/HotSpot-master" && pwd)"
mkdir -p outputs
HS="$ROOT/hotspot"
if [ ! -x "$HS" ]; then echo "compile HotSpot: cd $ROOT && make"; exit 1; fi
"$HS" -c case3_3chip.config -f case3_3chip.flp -p case3_3chip.ptrace -materials_file case3_3chip.materials -model_type grid -grid_rows 64 -grid_cols 64 -steady_file outputs/case3_3chip.steady -grid_steady_file outputs/case3_3chip.grid.steady
