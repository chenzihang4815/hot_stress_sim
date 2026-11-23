#!/usr/bin/env python3
"""给线长内环的接受步补 ATSim 应力。

被拒绝的步布局没有变，不重复算。应力是正常芯片内部 von Mises 最大。
布局按记录重放，重放结果先和已存的最优布局对一下。
"""

from __future__ import annotations

import json
import os
import random
import sys
from pathlib import Path

ROOT = Path("/home/czh/任务仿真产生功耗文件工具")
SRC = ROOT / "热仿真流水线" / "src"
sys.path.insert(0, str(SRC))

import run_closed_loop_multidie as cl
import run_closed_loop_suite as suite

SUITE = ROOT / "热仿真流水线" / "runs" / "closed_loop_wire"
CASES = tuple(s for s in ("n5_silicon", "n5_epoxy", "n8_silicon", "n8_epoxy", "n11_epoxy", "n11_silicon") if s in os.environ.get("STRESS_ONLY", "n5_silicon,n5_epoxy,n8_silicon,n8_epoxy,n11_epoxy,n11_silicon").split(","))


def replay(name: str, trace: list[dict]) -> dict[int, list[dict]]:
    cat = cl.catalog()
    current = suite.shelf(cat, name.split("_")[0])
    saved = {0: cl.pack(current)}
    rng_m = random.Random(suite.SEED)
    for step in range(1, suite.N_INNER + 1):
        roll = rng_m.random()
        if roll < 0.60:
            trial, op = cl.nudge(current, rng_m)
        elif roll < 0.80:
            trial, op = cl.rotate(current, rng_m)
        else:
            trial, op = cl.swap_hw(current, cat, rng_m)
        rec = next(r for r in trace if r.get("step") == step)
        if rec.get("op") != op:
            raise RuntimeError(f"{name} 第 {step} 步对不上: {op} / {rec.get('op')}")
        if rec.get("accept"):
            current = trial
            saved[step] = cl.pack(current)
    return saved


def main() -> None:
    for name in CASES:
        out = SUITE / name
        trace = json.loads((out / "trace.json").read_text(encoding="utf-8"))
        saved = replay(name, trace)
        best = json.loads((out / "layouts" / "best.json").read_text(encoding="utf-8"))
        last = saved[max(saved)]
        if [(d["name"], d["j0"], d["ib0"], d["w"], d["h"]) for d in last] != [
            (d["name"], d["j0"], d["ib0"], d["w"], d["h"]) for d in best
        ]:
            raise RuntimeError(f"{name} 重放和已存最优布局对不上")
        print(f"{name} replay ok, {len(saved)} accepted", flush=True)
        fill = name.split("_")[1]
        (out / "anneal_atsim").mkdir(exist_ok=True)
        stress_path = out / "anneal_stress.json"
        done = {}
        if stress_path.exists():
            done = {r["step"]: r for r in json.loads(stress_path.read_text(encoding="utf-8"))}
        rows = []
        for step in sorted(saved):
            if step in done:
                rows.append(done[step])
                continue
            result_path = out / "anneal_atsim" / f"step_{step:02d}_result.json"
            if not result_path.exists():
                suite.run_atsim(
                    saved[step], fill,
                    out / "anneal_atsim" / f"step_{step:02d}.json",
                    result_path,
                )
            result = json.loads(result_path.read_text(encoding="utf-8"))
            rows.append({
                "step": step,
                "die_cell_max_C": result["die_cell_max_C"],
                "vm_MPa": result["vm_normal_chip_max_MPa"],
                "vm_name": result["vm_normal_chip_name"],
            })
            stress_path.write_text(json.dumps(rows, indent=2) + "\n", encoding="utf-8")
            print(f"  step {step}  vm {rows[-1]['vm_MPa']:.1f}", flush=True)


if __name__ == "__main__":
    main()
