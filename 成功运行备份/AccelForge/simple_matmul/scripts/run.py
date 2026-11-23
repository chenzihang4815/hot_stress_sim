#!/usr/bin/env python3
"""Run official toy: simple.yaml + 64x32 matmul. Requires conda env accelforge."""
import json
from pathlib import Path
import accelforge as af
from accelforge.model.main import evaluate_mapping

ROOT = Path("/home/czh/任务仿真产生功耗文件工具/accelforge-main/accelforge-main")
EX = ROOT / "examples"
HERE = Path(__file__).resolve().parents[1]
OUT = HERE / "output"
OUT.mkdir(exist_ok=True)

M, KN = 64, 32
spec = af.Spec.from_yaml(
    EX / "arches" / "simple.yaml",
    EX / "workloads" / "basic" / "matmuls.yaml",
    EX / "mappings" / "unfused_matmuls_to_simple.yaml",
    jinja_parse_data={"N_EINSUMS": 1, "M": M, "KN": KN},
)
result = evaluate_mapping(spec)
summary = {
    "energy": float(result.energy()),
    "latency": float(result.latency()),
    "n_computes": float(result.n_computes()),
    "energy_per_component": {
        str(k): float(v) for k, v in result.energy(per_component=True).items()
    },
}
(OUT / "result.json").write_text(json.dumps(summary, indent=2))
print(json.dumps(summary, indent=2))
