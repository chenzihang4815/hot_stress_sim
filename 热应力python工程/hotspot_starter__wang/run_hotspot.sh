#!/usr/bin/env bash
set -e
ROOT="$(cd "$(dirname "$0")/../HotSpot-master/HotSpot-master" && pwd)"
mkdir -p outputs
HS="$ROOT/hotspot"
if [ ! -x "$HS" ]; then echo "compile HotSpot: cd $ROOT && make"; exit 1; fi
"$HS" -c cpu_gpu_cim_starter__wang.config -f cpu_gpu_cim_starter__wang.flp -p cpu_gpu_cim_starter__wang.ptrace -materials_file cpu_gpu_cim_starter__wang.materials -model_type grid -grid_rows 64 -grid_cols 64 -steady_file outputs/cpu_gpu_cim_starter__wang.steady -grid_steady_file outputs/cpu_gpu_cim_starter__wang.grid.steady
