#!/usr/bin/env bash
# Re-run the official gem5 ARM hello (SE + TimingSimpleCPU + L1/L2).
set -euo pipefail
export PATH="${HOME}/miniconda3/bin:${PATH}"
export LD_LIBRARY_PATH="${HOME}/miniconda3/lib:${LD_LIBRARY_PATH:-}"
GEM5_ROOT="/home/czh/任务仿真产生功耗文件工具/gem5-stable/gem5-stable"
HERE="$(cd "$(dirname "$0")/.." && pwd)"
cd "$GEM5_ROOT"
"$GEM5_ROOT/build/ARM/gem5.opt" -d "$HERE/output" \
  configs/deprecated/example/se.py \
  --cpu-type TimingSimpleCPU --caches --l2cache \
  --cmd tests/test-progs/hello/bin/arm/linux/hello
