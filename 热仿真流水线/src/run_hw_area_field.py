#!/usr/bin/env python3
"""同一件任务换硬件：面积变了，填充被硅替掉，场还差多少。

封装里 CIM 和 HBM 不动。CPU 三套、GPU 三套，任务各一件：

- CPU 任务：x86 192×192 matmul。功耗不是把这趟 matmul 在三块芯片上重跑，
  用的是该模板自己的 McPAT 380 K 总功耗（动态 + 亚阈 + 栅）。
- GPU 任务：AccelWattch 的 backprop，三张卡各一份 kernel 平均功耗。
  面积是裸片面积，不是仿真器算出来的。

每套都是一颗 CPU、一颗 GPU、一颗 HBM3、一颗 CIM_jia。另外四颗 CIM
小于一格（1.27 mm²），不进材料图。左下角对齐，往上、往右长。
CIM 放在最高那颗 CPU 上方两格，HBM 放在最高那颗 GPU 上方两格。
所以芯片越矮，和上面那颗之间的填充越厚；芯片越宽，和旁边那颗之间的填充越窄。

比较四件事，都在同一个离散里，不把绝对温度当成 HotSpot：

- naive：每颗芯片单独求解（别的位置仍是环氧），温升相加。
- frozen：材料图冻成“现在这四颗都是硅”，只轮流开功耗。应当只剩舍入误差。
- ring：功耗图不变，把“比最小那套多出来的格子”从硅改回环氧。
  这是填充被硅替掉本身的温差。
- 横向：从芯片顶边流向上面那颗、从芯片侧边流向旁边那颗的热，单位 W。
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

OUT = ROOT / "热仿真流水线" / "runs" / "hw_area_field"
N, CW = hs.N, hs.CW
T_AMB, K_SI, K_FILL = hs.T_AMB, hs.K_SI, hs.K_FILL
CELL_MM2 = (CW * 1e3) ** 2
GAP = 2  # 格。最大的那套和上面芯片之间留 2.25 mm 环氧。

# McPAT 官方模板，380 K。漏电不另做温度迭代：这一问看的是面积和材料图。
CPUS = (
    {"name": "Penryn", "area_mm2": 91.1348, "power_W": 47.9315 + 10.7217 + 0.588125},
    {"name": "T2", "area_mm2": 273.523, "power_W": 50.421 + 15.0682 + 3.50762},
    {"name": "Xeon", "area_mm2": 421.676, "power_W": 66.8485 + 35.4132 + 1.68429},
)
# 同一份 backprop。面积是流水线里用的裸片面积。
GPUS = (
    {"name": "GPU_3070", "area_mm2": 392.0, "power_W": 66.89},
    {"name": "GPU_2060S", "area_mm2": 445.0, "power_W": 71.12},
    {"name": "GPU_TITANX", "area_mm2": 471.0, "power_W": 67.41},
)
# 一颗 HBM3，一条 CPU 轨迹加一条 GPU 轨迹。公式同 run_full_chain.hbm_power(1, 1)。
HBM_POWER = (0.088 + (1.035 - 0.088) + (0.844 - 0.088)) * 16
HBM = {"name": "HBM3", "area_mm2": 121.0, "power_W": HBM_POWER}
CIM = {"name": "CIM_jia", "area_mm2": 4.281, "power_W": 0.002773}
TINY_CIM = (
    ("CIM_colonnade", 0.2314, 0.0412),
    ("CIM_wang", 0.122, 0.026),
    ("CIM_basic", 0.04672, 0.002625),
    ("CIM_sinangil", 0.002577, 0.0007948),
)


def snap(area_mm2: float) -> tuple[int, int]:
    """取最接近的矩形，长宽比不超过 1.6，并且让高不小于宽，方便往上长。"""
    target = area_mm2 / CELL_MM2
    best = None
    for w in range(1, 40):
        for h in range(w, int(w * 1.6) + 1):
            err = abs(w * h - target)
            if best is None or err < best[0]:
                best = (err, w, h)
    assert best is not None
    return best[1], best[2]


def prepare(spec: dict) -> dict:
    w, h = snap(spec["area_mm2"])
    out = dict(spec)
    out["w"] = w
    out["h"] = h
    out["grid_area_mm2"] = w * h * CELL_MM2
    return out


def place_all() -> tuple[list[dict], list[dict], dict, dict]:
    cpus = [prepare(s) for s in CPUS]
    gpus = [prepare(s) for s in GPUS]
    hbm = prepare(HBM)
    cim = prepare(CIM)
    max_cpu_w = max(c["w"] for c in cpus)
    max_cpu_h = max(c["h"] for c in cpus)
    max_gpu_h = max(g["h"] for g in gpus)
    max_gpu_w = max(g["w"] for g in gpus)
    cpu_j0, cpu_ib0 = 4, 4
    cim["j0"] = cpu_j0
    cim["ib0"] = cpu_ib0 + max_cpu_h + GAP
    gpu_j0 = cpu_j0 + max_cpu_w + GAP
    gpu_ib0 = 4
    hbm["j0"] = gpu_j0
    hbm["ib0"] = gpu_ib0 + max_gpu_h + GAP
    for c in cpus:
        c["j0"], c["ib0"] = cpu_j0, cpu_ib0
    for g in gpus:
        g["j0"], g["ib0"] = gpu_j0, gpu_ib0
    for d in cpus + gpus + [hbm, cim]:
        if d["j0"] + d["w"] > N or d["ib0"] + d["h"] > N:
            raise SystemExit(f"{d['name']} 超出 72 mm 网格")
    return cpus, gpus, hbm, cim


def cells(d: dict) -> list[tuple[int, int]]:
    """HotSpot 行号，第 0 行在 y 最大。"""
    i0 = N - (d["ib0"] + d["h"])
    i1 = N - d["ib0"]
    out = []
    for i in range(i0, i1):
        for j in range(d["j0"], d["j0"] + d["w"]):
            out.append((i, j))
    return out


def paint(dies: list[dict], k_fill: float, skip: set[tuple[int, int]] | None = None) -> np.ndarray:
    k = np.full((N, N), k_fill, dtype=float)
    for d in dies:
        for i, j in cells(d):
            if skip and (i, j) in skip:
                continue
            k[i, j] = K_SI
    return k


def power_map(dies: list[dict]) -> np.ndarray:
    p = np.zeros((N, N), dtype=float)
    for d in dies:
        cc = cells(d)
        share = d["power_W"] / len(cc)
        for i, j in cc:
            p[i, j] += share
    return p


def solve(k: np.ndarray, power: np.ndarray, factor=None, g_amb: float | None = None) -> np.ndarray:
    if factor is None:
        mat, g_amb = hs.build_operator(k)
        factor = splu(mat.tocsc())
    return hs.solve_state(k, power, factor=factor, g_amb=g_amb)


def stats(truth: np.ndarray, summed: np.ndarray) -> dict:
    err = summed - truth
    return {
        "Tmax_truth_C": float(truth.max() - 273.15),
        "Tmax_sum_C": float(summed.max() - 273.15),
        "dTmax_K": float(summed.max() - truth.max()),
        "max_abs_K": float(np.abs(err).max()),
        "rmse_K": float(np.sqrt(np.mean(err ** 2))),
    }


def rx(k: float) -> float:
    return CW / (k * CW * hs.LAYERS[0][1])


def edge_flow(T: np.ndarray, k: np.ndarray, die_cells: set[tuple[int, int]], di: int, dj: int) -> float:
    """从这颗芯片跨出指定方向的净热，W。正号表示热离开芯片。"""
    total = 0.0
    for i, j in die_cells:
        i2, j2 = i + di, j + dj
        if not (0 <= i2 < N and 0 <= j2 < N):
            continue
        if (i2, j2) in die_cells:
            continue
        r = 0.5 * rx(k[i, j]) + 0.5 * rx(k[i2, j2])
        total += (T[i, j] - T[i2, j2]) / r
    return float(total)


def up_flow(state: np.ndarray, k: np.ndarray, die_cells: set[tuple[int, int]]) -> float:
    """芯片层进入 TIM 的热，W。"""
    t = hs.LAYERS[0][1]
    total = 0.0
    for i, j in die_cells:
        rz = t / (k[i, j] * CW * CW)
        total += (state[0, i, j] - state[1, i, j]) / rz
    return float(total)


def gap_mm(a: dict, b: dict) -> dict:
    """a 在下、b 在上时的竖直间隙，以及左右相邻时的水平间隙。单位 mm。"""
    a_top = a["ib0"] + a["h"]
    b_bot = b["ib0"]
    a_right = a["j0"] + a["w"]
    vertical = (b_bot - a_top) * CW * 1e3
    horizontal = (b["j0"] - a_right) * CW * 1e3
    overlap_w = min(a["j0"] + a["w"], b["j0"] + b["w"]) - max(a["j0"], b["j0"])
    overlap_h = min(a["ib0"] + a["h"], b["ib0"] + b["h"]) - max(a["ib0"], b["ib0"])
    return {
        "vertical_gap_mm": vertical,
        "horizontal_gap_mm": horizontal,
        "facing_width_mm": max(0, overlap_w) * CW * 1e3,
        "facing_height_mm": max(0, overlap_h) * CW * 1e3,
    }


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    cpus, gpus, hbm, cim = place_all()
    ref_cpu, ref_gpu = cpus[0], gpus[0]
    ref_cells = set(cells(ref_cpu)) | set(cells(ref_gpu)) | set(cells(hbm)) | set(cells(cim))

    k_si = np.full((N, N), K_SI)
    mat_si, g_si = hs.build_operator(k_si)
    fac_si = splu(mat_si.tocsc())

    rows = []
    for cpu in cpus:
        for gpu in gpus:
            dies = [cpu, gpu, hbm, cim]
            cc = {d["name"]: set(cells(d)) for d in dies}
            now = set()
            for s in cc.values():
                now |= s
            extra = now - ref_cells
            k_ep = paint(dies, K_FILL)
            k_ring = k_ep.copy()
            for i, j in extra:
                k_ring[i, j] = K_FILL
            # k_ep 已经是硅在足迹上。环上改回环氧才是“多出来的仍是填充”。
            k_ep = paint(dies, K_FILL)
            k_ring = paint(dies, K_FILL, skip=extra)
            p = power_map(dies)

            st_ep = solve(k_ep, p)
            st_si = solve(k_si, p, fac_si, g_si)
            st_ring = solve(k_ring, p)

            parts = []
            for d in dies:
                alone = paint([d], K_FILL)
                field = solve(alone, power_map([d]))[0]
                parts.append(field - T_AMB)
            naive = np.sum(parts, axis=0) + T_AMB

            mat_fr, g_fr = hs.build_operator(k_ep)
            fac_fr = splu(mat_fr.tocsc())
            frozen = np.zeros((N, N))
            for d in dies:
                frozen += solve(k_ep, power_map([d]), fac_fr, g_fr)[0] - T_AMB
            frozen += T_AMB

            chip = st_ep[0]
            cpu_c, gpu_c = cc[cpu["name"]], cc[gpu["name"]]
            rec = {
                "cpu": cpu["name"],
                "gpu": gpu["name"],
                "cpu_area_mm2": cpu["area_mm2"],
                "gpu_area_mm2": gpu["area_mm2"],
                "cpu_grid_mm2": cpu["grid_area_mm2"],
                "gpu_grid_mm2": gpu["grid_area_mm2"],
                "cpu_wh_cells": [cpu["w"], cpu["h"]],
                "gpu_wh_cells": [gpu["w"], gpu["h"]],
                "cpu_power_W": cpu["power_W"],
                "gpu_power_W": gpu["power_W"],
                "spreader_contact_mm2": cpu["grid_area_mm2"] + gpu["grid_area_mm2"],
                "fill_replaced_mm2": len(extra) * CELL_MM2,
                "cpu_to_cim": gap_mm(cpu, cim),
                "gpu_to_hbm": gap_mm(gpu, hbm),
                "cpu_to_gpu_gap_mm": (gpu["j0"] - (cpu["j0"] + cpu["w"])) * CW * 1e3,
                "Tmax_epoxy_C": float(chip.max() - 273.15),
                "Tmax_silicon_C": float(st_si[0].max() - 273.15),
                "Tmax_cpu_C": float(max(chip[i, j] for i, j in cpu_c) - 273.15),
                "Tmax_gpu_C": float(max(chip[i, j] for i, j in gpu_c) - 273.15),
                "naive": stats(chip, naive),
                "frozen": stats(chip, frozen),
                "ring_epoxy_minus_silicon_dTmax_K": float(st_ring[0].max() - chip.max()),
                "ring_max_abs_K": float(np.abs(st_ring[0] - chip).max()),
                "Q_cpu_up_W": up_flow(st_ep, k_ep, cpu_c),
                "Q_gpu_up_W": up_flow(st_ep, k_ep, gpu_c),
                "Q_cpu_toward_cim_W": edge_flow(chip, k_ep, cpu_c, di=-1, dj=0),
                "Q_gpu_toward_hbm_W": edge_flow(chip, k_ep, gpu_c, di=-1, dj=0),
                "Q_cpu_toward_gpu_W": edge_flow(chip, k_ep, cpu_c, di=0, dj=1),
                "Q_gpu_toward_cpu_W": edge_flow(chip, k_ep, gpu_c, di=0, dj=-1),
            }
            for key in ("Q_cpu_toward_cim_W", "Q_gpu_toward_hbm_W", "Q_cpu_toward_gpu_W"):
                pass
            rec["cpu_lateral_frac"] = (cpu["power_W"] - rec["Q_cpu_up_W"]) / cpu["power_W"]
            rec["gpu_lateral_frac"] = (gpu["power_W"] - rec["Q_gpu_up_W"]) / gpu["power_W"]
            rows.append(rec)
            tag = f"{cpu['name']}_{gpu['name']}"
            np.save(OUT / f"{tag}_epoxy.npy", chip)
            np.save(OUT / f"{tag}_naive_err.npy", naive - chip)
            print(
                f"{cpu['name']:7s} {gpu['name']:10s}  "
                f"替掉 {rec['fill_replaced_mm2']:.0f} mm²  "
                f"CPU-CIM 间隙 {rec['cpu_to_cim']['vertical_gap_mm']:.2f} mm  "
                f"CPU-GPU 间隙 {rec['cpu_to_gpu_gap_mm']:.2f} mm  "
                f"naive {rec['naive']['max_abs_K']:.4f} K  "
                f"ring {rec['ring_max_abs_K']:.4f} K  "
                f"向上CIM {rec['Q_cpu_toward_cim_W']*1e3:.2f} mW  "
                f"侧向 {rec['Q_cpu_toward_gpu_W']*1e3:.1f} mW  "
                f"Tmax {rec['Tmax_epoxy_C']:.2f} °C",
                flush=True,
            )
            (OUT / "cases.json").write_text(json.dumps(rows, indent=2) + "\n", encoding="utf-8")

    report = {
        "task_cpu": "x86 192x192 matmul",
        "task_gpu": "AccelWattch backprop",
        "k_fill": K_FILL,
        "k_silicon": K_SI,
        "cell_mm": CW * 1e3,
        "gap_cells_at_largest": GAP,
        "hbm_power_W": HBM_POWER,
        "cim": {"name": CIM["name"], "area_mm2": CIM["area_mm2"], "power_W": CIM["power_W"],
                "grid_mm2": cim["grid_area_mm2"], "wh": [cim["w"], cim["h"]]},
        "tiny_cim_not_on_grid": [
            {"name": n, "area_mm2": a, "power_W": p} for n, a, p in TINY_CIM
        ],
        "note": (
            "绝对温度是没有 spreader 外伸环的 Python 模板，不是 HotSpot。"
            "naive 是每颗芯片单独求解再相加。ring 是多出来的格子改回环氧、功耗图不变。"
            "CPU 功耗是各模板 380 K 的 McPAT 总功耗，不是 matmul 在三块芯片上重跑。"
        ),
        "cpus": [{k: c[k] for k in ("name", "area_mm2", "grid_area_mm2", "power_W", "w", "h")} for c in cpus],
        "gpus": [{k: g[k] for k in ("name", "area_mm2", "grid_area_mm2", "power_W", "w", "h")} for g in gpus],
        "cases": rows,
    }
    (OUT / "report.json").write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"写出 {OUT / 'report.json'}", flush=True)


if __name__ == "__main__":
    main()
