#!/usr/bin/env bash
set -e
ROOT="$(cd "$(dirname "$0")/../HotSpot-master/HotSpot-master" && pwd)"
mkdir -p outputs
HS="$ROOT/hotspot"
if [ ! -x "$HS" ]; then echo "compile HotSpot: cd $ROOT && make"; exit 1; fi
"$HS" -c lbfgs_atomic__wang.config -f lbfgs_atomic__wang.flp -p lbfgs_atomic__wang.ptrace -materials_file lbfgs_atomic__wang.materials -model_type grid -grid_rows 64 -grid_cols 64 -steady_file outputs/lbfgs_atomic__wang.steady -grid_steady_file outputs/lbfgs_atomic__wang.grid.steady
