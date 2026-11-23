#!/usr/bin/env bash
set -e
ROOT="$(cd "$(dirname "$0")/../HotSpot-master/HotSpot-master" && pwd)"
mkdir -p outputs
HS="$ROOT/hotspot"
if [ ! -x "$HS" ]; then echo "compile HotSpot: cd $ROOT && make"; exit 1; fi
"$HS" -c lbfgs_L2__jia.config -f lbfgs_L2__jia.flp -p lbfgs_L2__jia.ptrace -materials_file lbfgs_L2__jia.materials -model_type grid -grid_rows 64 -grid_cols 64 -steady_file outputs/lbfgs_L2__jia.steady -grid_steady_file outputs/lbfgs_L2__jia.grid.steady
