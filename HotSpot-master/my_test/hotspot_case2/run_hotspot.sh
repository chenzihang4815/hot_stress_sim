#!/usr/bin/env bash
set -e
ROOT="$(cd "$(dirname "$0")/../HotSpot-master/HotSpot-master" && pwd)"
mkdir -p outputs
HS="$ROOT/hotspot"
if [ ! -x "$HS" ]; then echo "compile HotSpot: cd $ROOT && make"; exit 1; fi
"$HS" -c case2_package.config -f case2_package.flp -p case2_package.ptrace -materials_file case2_package.materials -model_type grid -grid_rows 64 -grid_cols 64 -steady_file outputs/case2_package.steady -grid_steady_file outputs/case2_package.grid.steady
