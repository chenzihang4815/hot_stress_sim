#!/usr/bin/env python3
"""芯片是硅，芯片以外的填充换成别的材料。

不改流水线默认。HotSpot / ATSim 仍然把芯片层当成连续硅。
这里只回答：填充不用硅之后，场叠加还差多少、峰值温度动多少。

填充取三种，覆盖“随便”：
- 环氧模塑料 0.8 W/m·K（模封料常见值）
- 底部填充取 1.0，和模塑料同一量级
- 空气 0.028，HotSpot package.h 里的下界

只跑最小的一套（Penryn + 3070）和最大的一套（Xeon + TITAN X）。
布局和功耗与 hw_area_field 相同。
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from scipy.sparse.linalg import splu

ROOT = Path("/home/czh/任务仿真产生功耗文件工具")
SRC = ROOT / "热仿真流水线" / "src"
sys.path.insert(0, str(SRC))

import run_fill_superposition as hs  # noqa: E402
import run_hw_area_field as hw  # noqa: E402

OUT = ROOT / "热仿真流水线" / "runs" / "hw_fill_material"
K_SI = hs.K_SI
T_AMB = hs.T_AMB

# 芯片永远是硅。下面只换芯片之外的填充。
FILLS = (
    {"name": "emc", "note": "环氧模塑料", "k": 0.8},
    {"name": "underfill", "note": "底部填充，和模塑料同一量级", "k": 1.0},
    {"name": "air", "note": "空气，HotSpot package.h 的 0.028", "k": 0.028},
)


def one(dies: list[dict], k_fill: float, fac_si, g_si) -> dict:
    k = hw.paint(dies, k_fill)
    p = hw.power_map(dies)
    truth = hw.solve(k, p)[0]
    parts = []
    for d in dies:
        alone = hw.paint([d], k_fill)
        parts.append(hw.solve(alone, hw.power_map([d]))[0] - T_AMB)
    naive = np.sum(parts, axis=0) + T_AMB
    mat, g = hs.build_operator(k)
    fac = splu(mat.tocsc())
    frozen = np.zeros_like(truth)
    for d in dies:
        frozen += hw.solve(k, hw.power_map([d]), fac, g)[0] - T_AMB
    frozen += T_AMB
    silicon = hw.solve(np.full_like(k, K_SI), p, fac_si, g_si)[0]
    err = naive - truth
    return {
        "k_fill": k_fill,
        "Tmax_C": float(truth.max() - 273.15),
        "Tmax_silicon_fill_C": float(silicon.max() - 273.15),
        "fill_minus_silicon_K": float(truth.max() - silicon.max()),
        "naive_max_abs_K": float(np.abs(err).max()),
        "naive_dTmax_K": float(naive.max() - truth.max()),
        "frozen_max_abs_K": float(np.abs(frozen - truth).max()),
    }


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    cpus, gpus, hbm, cim = hw.place_all()
    picks = [(cpus[0], gpus[0]), (cpus[-1], gpus[-1])]
    k_si = np.full((hs.N, hs.N), K_SI)
    mat_si, g_si = hs.build_operator(k_si)
    fac_si = splu(mat_si.tocsc())
    rows = []
    for cpu, gpu in picks:
        dies = [cpu, gpu, hbm, cim]
        for fill in FILLS:
            rec = one(dies, fill["k"], fac_si, g_si)
            rec.update({
                "cpu": cpu["name"], "gpu": gpu["name"],
                "fill": fill["name"], "fill_note": fill["note"],
                "cpu_power_W": cpu["power_W"], "gpu_power_W": gpu["power_W"],
            })
            rows.append(rec)
            print(
                f"{cpu['name']:7s} {gpu['name']:10s} {fill['name']:10s} k={fill['k']:<6}  "
                f"Tmax {rec['Tmax_C']:.3f} °C  比全硅 {rec['fill_minus_silicon_K']:+.4f} K  "
                f"naive {rec['naive_max_abs_K']:.4f} K  frozen {rec['frozen_max_abs_K']:.2e}",
                flush=True,
            )
            (OUT / "cases.json").write_text(json.dumps(rows, indent=2) + "\n", encoding="utf-8")
    report = {
        "chip_k": K_SI,
        "chip_material": "silicon",
        "fills": list(FILLS),
        "note": (
            "芯片足迹是硅 k=130。芯片以外的芯片层换成表里的填充。"
            "绝对温度不是 HotSpot。默认流水线的填充没有改。"
        ),
        "cases": rows,
    }
    (OUT / "report.json").write_text(
        json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"写出 {OUT / 'report.json'}", flush=True)


if __name__ == "__main__":
    main()
