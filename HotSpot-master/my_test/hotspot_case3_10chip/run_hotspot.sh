#!/usr/bin/env bash
set -e
ROOT="$(cd "$(dirname "$0")/../HotSpot-master/HotSpot-master" && pwd)"
mkdir -p outputs
HS="$ROOT/hotspot"
if [ ! -x "$HS" ]; then echo "compile HotSpot: cd $ROOT && make"; exit 1; fi
"$HS" -c case3_10chip.config -f case3_10chip.flp -p case3_10chip.ptrace -materials_file case3_10chip.materials -model_type grid -grid_rows 64 -grid_cols 64 -steady_file outputs/case3_10chip.steady -grid_steady_file outputs/case3_10chip.grid.steady
