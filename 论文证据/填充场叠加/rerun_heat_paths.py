#!/usr/bin/env python3
"""复现 2026-11-13 的散热路径实验。

当时没有单独的设计记录，脚本是在命令行里用 heredoc 写的，跑完没有存。
这里按会话记录里留下的那一版重写，口径和数字都对得上
热仿真流水线/runs/interposer_map/bottom_exit.json 和 fill_with_both_exits.json。

模型是 run_interposer_map.solve_with_interposer 的五层：中介层 100 µm 硅、
芯片层 150 µm、TIM 20 µm（k=4）、铜均热片 1 mm、热沉 6.9 mm。环境 45 °C。
芯片是 run_factor_controls.base_dies 的四颗：Penryn、GPU_3070、HBM3、CIM_jia，
总功耗 154.790098 W。变的是芯片层缝里的填充，芯片占的格子固定是硅 130。
python 场，不是 HotSpot，也不是 ATSim。

热流从求解后的温度场积分，不另解。顶上是热沉每一格的对流。底面是中介层每一格
的底面电导乘温差。侧边是中介层最外一圈的边界电导 1e4 W/K 乘温差。三条加起来
等于总功耗。

电导：顶上对流系数 h = 1 / (R_conv · 0.096²) = 2777.78 W/m²K。底面挂和顶上
相同的散热器时，每格电导是 h 乘格子面积。底面先过 1 mm、k=0.5 的有机层时，
每格电导是格子面积除以 (1e-3/0.5 + 1/h)。
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path("/home/czh/任务仿真产生功耗文件工具")
sys.path.insert(0, str(ROOT / "热仿真流水线" / "src"))

import run_factor_controls as fc  # noqa: E402
import run_hw_area_field as hw  # noqa: E402
import run_interposer_map as im  # noqa: E402

FILLS = [("air", 0.028), ("organic", 0.5), ("epoxy", 0.8),
         ("underfill", 1.0), ("silicon", 130.0)]
G_EDGE = 1.0e4
T_CHIP = 150e-6

AREA = im.hs.CW * im.hs.CW
H_TOP = 1.0 / (im.hs.R_CONVEC * im.hs.S_SINK * im.hs.S_SINK)
G_SAME = H_TOP * AREA
G_ORG = AREA / (1.0e-3 / 0.5 + 1.0 / H_TOP)
R_SINK = im.hs._rz(400.0, 6.9e-3) + im.hs.R_CONVEC * (im.hs.S_SINK ** 2) / AREA
G_AMB = 1.0 / R_SINK
G_DOWN = 1.0 / im.hs._rz(im.K_IP, im.T_IP)


def flows(st, k_chip, g_bottom, edge, silicon):
    ip, chip, sink = st[0], st[1], st[4]
    q_top = float((G_AMB * (sink - im.T_AMB)).sum())
    q_bot = float((g_bottom * (ip - im.T_AMB)).sum())
    border = np.zeros_like(ip, dtype=bool)
    border[0, :] = border[-1, :] = border[:, 0] = border[:, -1] = True
    q_edge = float((G_EDGE * (ip - im.T_AMB))[border].sum()) if edge else 0.0
    q_down = G_DOWN * (chip - ip)
    q_lat = 0.0
    for i in range(im.N):
        for j in range(im.N):
            if not silicon[i, j]:
                continue
            for di, dj in ((0, 1), (1, 0), (0, -1), (-1, 0)):
                ni, nj = i + di, j + dj
                if not (0 <= ni < im.N and 0 <= nj < im.N) or silicon[ni, nj]:
                    continue
                r = (0.5 * im.hs._rx(float(k_chip[i, j]), T_CHIP)
                     + 0.5 * im.hs._rx(float(k_chip[ni, nj]), T_CHIP))
                q_lat += (chip[i, j] - chip[ni, nj]) / r
    return {
        "Q_top_W": q_top,
        "Q_bottom_W": q_bot,
        "Q_edge_W": q_edge,
        "Q_sum_W": q_top + q_bot + q_edge,
        "Q_down_under_dies_W": float(q_down[silicon].sum()),
        "Q_down_under_gaps_W": float(q_down[~silicon].sum()),
        "Q_lateral_die_to_gap_W": q_lat,
    }


def one(dies, power, silicon, name, k_fill, g_bottom, edge):
    k_chip = hw.paint(dies, k_fill)
    st = im.solve_with_interposer(k_chip, power, 4.0, g_bottom=g_bottom,
                                  edge_to_ambient=edge)
    chip = st[1]
    rec = {
        "name": name, "k_fill": k_fill,
        "g_bottom_W_per_K": g_bottom, "edge_to_ambient": edge,
        "python_Tmax_C": float(chip.max() - 273.15),
        "dies_C": {d["name"]: im.die_max(chip, d) for d in dies},
    }
    rec.update(flows(st, k_chip, g_bottom, edge, silicon))
    return rec


def check(tag, rec, archived):
    keys = ("python_Tmax_C", "Q_top_W", "Q_bottom_W", "Q_edge_W")
    worst = max(abs(rec[k] - archived[k]) for k in keys)
    print(f"  {tag:22s} 和存档最大差 {worst:.2e}", flush=True)
    return worst


def main():
    dies, _ = fc.base_dies()
    power = hw.power_map(dies)
    silicon = np.zeros((im.N, im.N), dtype=bool)
    for d in dies:
        for i, j in hw.cells(d):
            silicon[i, j] = True
    print(f"总功耗 {power.sum():.6f} W  顶上对流 {H_TOP:.2f} W/m²K  "
          f"底面相同散热器 {G_SAME:.6e}  底面有机层 {G_ORG:.6e}", flush=True)

    out = Path(__file__).resolve().parent / "复现"
    out.mkdir(exist_ok=True)
    archived = ROOT / "热仿真流水线" / "runs" / "interposer_map"
    old_exit = {r["name"]: r for r in json.loads(
        (archived / "bottom_exit.json").read_text())["rows"]}
    old_fill = {r["fill"]: r for r in json.loads(
        (archived / "fill_with_both_exits.json").read_text())["rows"]}

    exits = [
        ("adiabatic", 0.8, 0.0, False),
        ("edge_to_ambient", 0.8, 0.0, True),
        ("bottom_same_h", 0.8, G_SAME, False),
        ("bottom_organic_1mm", 0.8, G_ORG, False),
        ("bottom_same_h_air", 0.028, G_SAME, False),
        ("bottom_organic_air", 0.028, G_ORG, False),
        ("edge_plus_organic", 0.8, G_ORG, True),
    ]
    exit_rows = []
    for name, k, g, edge in exits:
        rec = one(dies, power, silicon, name, k, g, edge)
        exit_rows.append(rec)
        check(name, rec, old_exit[name])
    (out / "exits.json").write_text(
        json.dumps(exit_rows, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    fill_rows = []
    for name, k in FILLS:
        rec = one(dies, power, silicon, name, k, G_ORG, True)
        fill_rows.append(rec)
        check(name, rec, old_fill[name])
    (out / "fills.json").write_text(
        json.dumps(fill_rows, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"写到 {out}", flush=True)


if __name__ == "__main__":
    main()
