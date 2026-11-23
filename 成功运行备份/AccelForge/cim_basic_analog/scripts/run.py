#!/usr/bin/env python3
"""Run official CIM basic_analog. Requires conda env accelforge."""
import json
import sys
from pathlib import Path
import accelforge as af

CIM = Path(
    "/home/czh/任务仿真产生功耗文件工具/accelforge-main/accelforge-main"
    "/examples/arches/compute_in_memory"
)
sys.path.insert(0, str(CIM))
from _load_spec import get_spec  # noqa: E402

HERE = Path(__file__).resolve().parents[1]
OUT = HERE / "output"
OUT.mkdir(exist_ok=True)

basic_path = CIM / "basic_analog.yaml"
workload = af.Workload.from_yaml(basic_path, top_key="workload")
renames = af.Renames.from_yaml(basic_path, top_key="renames")
spec = get_spec("basic_analog", add_dummy_main_memory=True, n_macros=1)
spec.workload = workload
spec.renames = renames
spec.mapper.metrics = af.Metrics.ENERGY
spec.mapper.explore_loop_orders = False
result = spec.map_workload_to_arch(print_progress=False)
result = result.drop_components_with_zero_energy_and_latency()
evaluated = spec.calculate_component_costs()
n = float(result.n_computes())
e = float(result.energy())
summary = {
    "energy_J": e,
    "latency_s": float(result.latency()),
    "n_computes": n,
    "energy_per_MAC_J": e / n,
    "energy_per_component_J": {
        str(k): float(v)
        for k, v in result.energy(per_component=True).items()
        if float(v) != 0
    },
    "area_m2": {
        k: float(v) for k, v in evaluated.arch.per_component_total_area.items()
    },
}
(OUT / "result.json").write_text(json.dumps(summary, indent=2))
print(json.dumps(summary, indent=2))
