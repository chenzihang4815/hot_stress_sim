#!/usr/bin/env python3
"""老化第一步：垂直热阻变为 2 倍和 3 倍，坐标不动，只换填充。

垂直路径按模型里的串联热阻：整层芯片厚度加上整层 TIM。
硅的导热系数不改，只把 TIM 的导热系数降到让这段串联热阻正好是 2 倍、3 倍。
均热片仍是 400。温度是 python 场，不是 HotSpot。
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path("/home/czh/任务仿真产生功耗文件工具")
SRC = ROOT / "热仿真流水线" / "src"
sys.path.insert(0, str(SRC))

import run_factor_controls as fc  # noqa: E402
import run_hw_area_field as hw  # noqa: E402

OUT = ROOT / "热仿真流水线" / "runs" / "aging_fill"
T_CHIP = 150e-6
T_TIM = 20e-6
K_SI = 130.0
K_TIM0 = 4.0
R0 = T_CHIP / K_SI + T_TIM / K_TIM0

# 现成材料的导热系数，不新编。
FILLS = (
    ("silicon", 130.0),
    ("underfill", 1.0),
    ("epoxy", 0.8),
    ("organic", 0.5),
    ("air", 0.028),
)


def k_tim_for(factor: float) -> float:
    r_tim = factor * R0 - T_CHIP / K_SI
    return T_TIM / r_tim


def die_max(T: np.ndarray, die: dict) -> float:
    vals = [T[i, j] for i, j in hw.cells(die)]
    return float(max(vals) - 273.15)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    dies, _ = fc.base_dies()
    rows = []
    fields = {}
    for factor in (1.0, 2.0, 3.0):
        k_tim = k_tim_for(factor)
        fields[factor] = {}
        for name, k_fill in FILLS:
            st = fc.truth(dies, k_fill, k_tim, 400.0)
            T = st[0]
            rec = {
                "vertical_factor": factor,
                "k_tim": k_tim,
                "fill": name,
                "k_fill": k_fill,
                "python_Tmax_C": fc.tmax(T),
                "dies_C": {d["name"]: die_max(T, d) for d in dies},
            }
            rows.append(rec)
            fields[factor][name] = T
            print(f"  x{factor:.0f} {name:10} k {k_fill:7.3f}  Tmax {rec['python_Tmax_C']:.3f}", flush=True)
    # 温度差：同一垂直热阻下，相对环氧。
    deltas = []
    for factor in (1.0, 2.0, 3.0):
        base = next(r for r in rows if r["vertical_factor"] == factor and r["fill"] == "epoxy")
        Tb = fields[factor]["epoxy"]
        for name, _k in FILLS:
            r = next(x for x in rows if x["vertical_factor"] == factor and x["fill"] == name)
            dT = fields[factor][name] - Tb
            deltas.append({
                "vertical_factor": factor,
                "fill": name,
                "dTmax_vs_epoxy_K": r["python_Tmax_C"] - base["python_Tmax_C"],
                "dT_field_max_K": float(np.max(np.abs(dT))),
                "dies_dTmax_K": {
                    n: r["dies_C"][n] - base["dies_C"][n] for n in r["dies_C"]
                },
            })
    report = {
        "layout": [fc.sm.dump_die(d) for d in dies],
        "note": (
            "坐标是仓库起点，不移动。垂直热阻是芯片厚度/130 加上 TIM 厚度/k_TIM，"
            "乘 2 和 3 时只改 k_TIM，均热片保持 400。"
            "python 峰值不是 HotSpot。功耗是各模板自己的 McPAT 380 K 总功耗和同一份 backprop。"
        ),
        "R0_m2K_W": R0,
        "k_tim": {str(f): k_tim_for(f) for f in (1.0, 2.0, 3.0)},
        "rows": rows,
        "delta_vs_epoxy": deltas,
    }
    (OUT / "report.json").write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"写出 {OUT / 'report.json'}", flush=True)


if __name__ == "__main__":
    main()
