#!/usr/bin/env bash
set -e
ROOT="$(cd "$(dirname "$0")/../HotSpot-master/HotSpot-master" && pwd)"
mkdir -p outputs
HS="$ROOT/hotspot"
if [ ! -x "$HS" ]; then echo "compile HotSpot: cd $ROOT && make"; exit 1; fi
"$HS" -c lbfgs_starter__colonnade.config -f lbfgs_starter__colonnade.flp -p lbfgs_starter__colonnade.ptrace -materials_file lbfgs_starter__colonnade.materials -model_type grid -grid_rows 64 -grid_cols 64 -steady_file outputs/lbfgs_starter__colonnade.steady -grid_steady_file outputs/lbfgs_starter__colonnade.grid.steady
