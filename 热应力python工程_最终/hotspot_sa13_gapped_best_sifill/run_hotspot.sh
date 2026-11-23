#!/usr/bin/env bash
set -e
ROOT="$(cd "$(dirname "$0")/../HotSpot-master/HotSpot-master" && pwd)"
mkdir -p outputs
HS="$ROOT/hotspot"
if [ ! -x "$HS" ]; then echo "compile HotSpot: cd $ROOT && make"; exit 1; fi
"$HS" -c sa13_gapped_best_sifill.config -f sa13_gapped_best_sifill.flp -p sa13_gapped_best_sifill.ptrace -materials_file sa13_gapped_best_sifill.materials -model_type grid -grid_rows 64 -grid_cols 64 -steady_file outputs/sa13_gapped_best_sifill.steady -grid_steady_file outputs/sa13_gapped_best_sifill.grid.steady
