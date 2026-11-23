#!/usr/bin/env bash
set -e
ROOT="$(cd "$(dirname "$0")/../HotSpot-master/HotSpot-master" && pwd)"
mkdir -p outputs
HS="$ROOT/hotspot"
if [ ! -x "$HS" ]; then echo "compile HotSpot: cd $ROOT && make"; exit 1; fi
"$HS" -c lbfgs_starter__jia.config -f lbfgs_starter__jia.flp -p lbfgs_starter__jia.ptrace -materials_file lbfgs_starter__jia.materials -model_type grid -grid_rows 64 -grid_cols 64 -steady_file outputs/lbfgs_starter__jia.steady -grid_steady_file outputs/lbfgs_starter__jia.grid.steady
