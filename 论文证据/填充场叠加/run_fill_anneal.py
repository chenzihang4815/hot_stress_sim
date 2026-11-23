#!/usr/bin/env python3
"""自己的三芯粒方案，不同填充，材料场更新对不更新。

芯片是 sa_heat_vs_coupled_pkg_1cpu_1gpu_1cim 的 CPU_matmul + GPU_3070 + CIM_jia，
72 mm 中介层，起点是那一批的 layout_start.json。

求解用 run_interposer_map 的五层冷边算子：中介层 100 µm 硅、芯片层 150 µm、
TIM 20 µm、铜均热片 1 mm、热沉 6.9 mm，板边接到 45 °C，底面绝热。芯片层缝里的
导热系数按填充材料改，芯片占的格子固定是硅 130。这是 python 场，不是 HotSpot，
也不是 ATSim。

两种做法走同一批提议。exact 每步按最新布局重画材料图、重分解、重解。frozen
只在起点分解一次，后面每步把三颗芯片的功耗分开，在起点的材料图上回代再相加，
材料图不再更新。接受只按 exact 的峰值温度决定，两条记录的是同一条布局序列上
frozen 对 exact 的误差，所以误差是材料图过期造成的，不是两条退火分叉造成的。

移动用 sa_thermal_feedback.propose：平移、跳动、对调、旋转，比例 40/25/20/15，
种子 22，温度从 2.0 起每步乘 0.95，120 步。
"""

from __future__ import annotations

import json
import os
import random
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path("/home/czh/任务仿真产生功耗文件工具")
SRC = ROOT / "热仿真流水线" / "src"
sys.path.insert(0, str(SRC))

import run_hw_area_field as hw  # noqa: E402
import run_interposer_map as ip  # noqa: E402
import sa_thermal_feedback as sa  # noqa: E402

N = ip.N
T_AMB = ip.T_AMB
CW = hw.hs.CW
K_SI = 130.0
SIDE = 0.072

FILLS = {
    "air": 0.028,
    "epoxy": 0.8,
    "underfill": 1.0,
    "silicon": 130.0,
}

STEPS = int(os.environ.get("FILL_STEPS", "120"))
SEED = 22
T0 = 2.0
ALPHA = 0.95


def to_dies(chips) -> list[dict]:
    """米制坐标换成 64 网格的格子。j0 是列，ib0 是从下往上的行。"""
    dies = []
    for c in chips:
        j0 = int(round(c.x / CW))
        ib0 = int(round(c.y / CW))
        w = max(1, int(round(c.w / CW)))
        h = max(1, int(round(c.h / CW)))
        dies.append({
            "name": c.name, "kind": c.kind, "power_W": c.power_W,
            "j0": j0, "ib0": ib0, "w": w, "h": h,
        })
    return dies


def paint(dies, k_fill: float) -> np.ndarray:
    k = np.full((N, N), k_fill)
    for d in dies:
        for i, j in hw.cells(d):
            k[i, j] = K_SI
    return k


def one_power(dies, idx: int) -> np.ndarray:
    p = np.zeros((N, N))
    d = dies[idx]
    cc = hw.cells(d)
    share = d["power_W"] / len(cc)
    for i, j in cc:
        p[i, j] = share
    return p


def die_max_c(T: np.ndarray, dies) -> dict:
    return {d["name"]: float(max(T[i, j] for i, j in hw.cells(d)) - 273.15) for d in dies}


def layout_record(chips) -> list[dict]:
    return [
        {"name": c.name, "x_mm": c.x * 1e3, "y_mm": c.y * 1e3,
         "w_mm": c.w * 1e3, "h_mm": c.h * 1e3, "power_W": c.power_W}
        for c in chips
    ]


def run(fill: str, out: Path) -> None:
    k_fill = FILLS[fill]
    out.mkdir(parents=True, exist_ok=True)
    src = ROOT / "热仿真流水线" / "runs" / "sa_heat_vs_coupled_pkg_1cpu_1gpu_1cim" / "layout_start.json"
    chips = sa.load_chips(src)

    rng = random.Random(SEED)
    cur = chips
    dies0 = to_dies(cur)
    k0 = paint(dies0, k_fill)

    t0 = time.perf_counter()
    fac, g, _ = ip.factor_edge(k0)
    t_factor = time.perf_counter() - t0

    def frozen_of(dies) -> np.ndarray:
        parts = [ip.apply_edge(fac, g, one_power(dies, i)) for i in range(len(dies))]
        return sum(p - T_AMB for p in parts) + T_AMB

    exact_secs = 0.0
    frozen_secs = 0.0
    exact_solves = 0
    best_peak = float("inf")
    best_layout = layout_record(cur)
    best_step = 0
    temp = T0
    rows = []
    err_max = np.zeros((N, N))

    T_exact0 = ip.solve_edge(k0, hw.power_map(dies0))
    exact_solves += 1
    T_frozen0 = frozen_of(dies0)
    d0 = T_frozen0 - T_exact0
    cur_peak = float(T_exact0.max() - 273.15)
    best_peak = cur_peak
    rows.append({
        "step": 0, "op": "start", "accept": True,
        "exact_Tmax_C": cur_peak,
        "frozen_Tmax_C": float(T_frozen0.max() - 273.15),
        "peak_err_K": float(T_frozen0.max() - T_exact0.max()),
        "grid_max_K": float(np.max(np.abs(d0))),
        "grid_mean_K": float(np.mean(np.abs(d0))),
    })

    grids = {"step0_exact": T_exact0 - 273.15}

    for step in range(1, STEPS + 1):
        trial, op = sa.propose(cur, SIDE, SIDE, rng)
        accepted = False
        if sa.legal(trial, SIDE, SIDE):
            dies = to_dies(trial)
            t1 = time.perf_counter()
            T_exact = ip.solve_edge(paint(dies, k_fill), hw.power_map(dies))
            exact_secs += time.perf_counter() - t1
            exact_solves += 1
            t1 = time.perf_counter()
            T_frozen = frozen_of(dies)
            frozen_secs += time.perf_counter() - t1
            peak = float(T_exact.max() - 273.15)
            diff = peak - cur_peak
            accepted = diff <= 0 or rng.random() < np.exp(-diff / temp)
            d = T_frozen - T_exact
            err_max = np.maximum(err_max, np.abs(d))
            rows.append({
                "step": step, "op": op, "accept": bool(accepted),
                "exact_Tmax_C": peak,
                "frozen_Tmax_C": float(T_frozen.max() - 273.15),
                "peak_err_K": float(T_frozen.max() - T_exact.max()),
                "grid_max_K": float(np.max(np.abs(d))),
                "grid_mean_K": float(np.mean(np.abs(d))),
            })
            if accepted:
                cur = trial
                cur_peak = peak
                if peak < best_peak:
                    best_peak = peak
                    best_layout = layout_record(cur)
                    best_step = step
                    grids["best_exact"] = T_exact - 273.15
                    grids["best_frozen"] = T_frozen - 273.15
                    grids["best_err"] = d
        else:
            rows.append({"step": step, "op": op, "accept": False, "illegal": True})
        temp *= ALPHA
        if step % 20 == 0:
            print(f"  {fill} 第 {step} 步 峰值 {cur_peak:.2f} 最优 {best_peak:.2f}", flush=True)

    summary = {
        "fill": fill,
        "k_fill": k_fill,
        "steps": STEPS,
        "seed": SEED,
        "start_factor_s": t_factor,
        "exact_solve_s": exact_secs,
        "frozen_solve_s": frozen_secs,
        "exact_solves": exact_solves,
        "start_Tmax_C": rows[0]["exact_Tmax_C"],
        "best_Tmax_C": best_peak,
        "best_step": best_step,
        "err_max_over_steps_K": float(err_max.max()),
        "note": (
            "120 步，种子 22，T0 2.0，每步乘 0.95。提议是生产退火的平移、跳动、对调、旋转。"
            "接受由 exact 决定，两条走同一批布局。frozen 用起点材料图的分解回代，不重分解。"
            "exact 每步按最新材料图重分解重解。err_max_over_steps 是每步误差场逐格取绝对值的最大。"
            "冷边 45 C，底面绝热。python，不是 HotSpot，也不是 ATSim。"
        ),
    }
    (out / "summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    (out / "trace.json").write_text(
        json.dumps(rows, indent=1, ensure_ascii=False) + "\n", encoding="utf-8")
    (out / "best_layout.json").write_text(
        json.dumps(best_layout, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    np.savez_compressed(out / "grids.npz", err_max_over_steps=err_max, **grids)
    print(json.dumps({k: summary[k] for k in (
        "fill", "start_factor_s", "exact_solve_s", "frozen_solve_s",
        "start_Tmax_C", "best_Tmax_C", "err_max_over_steps_K")}), flush=True)


if __name__ == "__main__":
    fill = os.environ["FILL_CASE"]
    out = Path(os.environ["FILL_OUT"])
    run(fill, out)
