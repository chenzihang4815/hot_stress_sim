#!/usr/bin/env python3
"""温度对线长的帕累托。

两组布局合在一起：400 个随机布局，加上两颗热源从贴合到拉开的 15 个。
都是 5 芯粒、环氧、芯片固定、频率固定，只动位置。温度和线长都取 Python 场，
不跑 ATSim。
"""

from __future__ import annotations

import json
import sys
from copy import deepcopy
from pathlib import Path

SRC = Path("/home/czh/任务仿真产生功耗文件工具/热仿真流水线/src")
sys.path.insert(0, str(SRC))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import font_manager
from scipy.sparse.linalg import splu

import run_closed_loop_multidie as cl
import run_fill_superposition as hs
import run_hw_area_field as hw
import run_wirelength as wire

HERE = Path("/home/czh/任务仿真产生功耗文件工具/论文证据/帕累托_布局")
OUT = HERE / "runs"
FONT = "/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc"
BASE = json.loads((OUT / "layouts" / "init.json").read_text(encoding="utf-8"))


def separated(step: int) -> list[dict] | None:
    """五颗从左下角一簇起步，每步把每颗芯粒沿对角线往外推一格。

    线长随步数单调增加，热源同时在分开，温度下降。
    """
    dies = deepcopy(BASE)
    by = {d["slot"]: d for d in dies}
    home = {"Ca": (0, 0), "Ga": (17, 0), "HBM3": (0, 22),
            "CIM_colonnade_m512": (17, 22), "CIM_jia_m64": (26, 22)}
    push = {"Ca": (0, 0), "Ga": (1, 0), "HBM3": (0, 1),
            "CIM_colonnade_m512": (1, 1), "CIM_jia_m64": (1, 1)}
    for slot, die in by.items():
        j, i = home[slot]
        dj, di = push[slot]
        die["j0"] = j + dj * step
        die["ib0"] = i + di * step
        if die["j0"] + die["w"] > hs.N or die["ib0"] + die["h"] > hs.N:
            return None
    return dies


def score(dies: list[dict], fac, g_amb, matrix: list[list[int]]) -> dict:
    cl.apply_power(dies)
    temp = hw.solve(hw.paint(dies, hs.K_SI), hw.power_map(dies), fac, g_amb)[0]
    ev = cl.metrics(temp)
    return {"T_C": ev["Tmax_C"], "wl_mm": wire.wirelength_cells(dies, matrix, cl.CW_MM)}


def main() -> None:
    rows = json.loads((OUT / "scan.json").read_text(encoding="utf-8"))
    for r in rows:
        r["T_C"] = r.pop("T_py_C")
        r["group"] = "随机布局"

    k_si = hw.paint(deepcopy(BASE), hs.K_SI)
    mat, g_amb = hs.build_operator(k_si)
    fac = splu(mat.tocsc())
    matrix = wire.assumed_matrix(BASE)
    for step in range(0, 40, 2):
        dies = separated(step)
        if dies is None:
            break
        row = score(dies, fac, g_amb, matrix)
        row.update({"tag": f"spread{step:02d}", "group": "逐步打散"})
        rows.append(row)
        print(f"{row['tag']:10} T {row['T_C']:.2f}  wl {row['wl_mm']:.2f}", flush=True)

    def dominates(a, b) -> bool:
        return (a["T_C"] <= b["T_C"] and a["wl_mm"] <= b["wl_mm"]
                and (a["T_C"] < b["T_C"] or a["wl_mm"] < b["wl_mm"]))

    for r in rows:
        r["on_front"] = not any(dominates(o, r) for o in rows if o is not r)
    front = sorted((r for r in rows if r["on_front"]), key=lambda r: r["T_C"])
    print(f"front {len(front)} / {len(rows)}", flush=True)
    (OUT / "pareto_tw.json").write_text(
        json.dumps({"n": len(rows), "n_front": len(front), "rows": rows},
                   indent=1, ensure_ascii=False) + "\n", encoding="utf-8")

    font_manager.fontManager.addfont(FONT)
    plt.rcParams["font.family"] = "Noto Sans CJK JP"
    plt.rcParams["axes.unicode_minus"] = False
    fig, ax = plt.subplots(figsize=(7.2, 5.0))
    rest = [r for r in rows if not r["on_front"] and r["group"] == "随机布局"]
    ax.scatter([r["wl_mm"] for r in rest], [r["T_C"] for r in rest],
               s=12, color="#c5c1b9", label="被支配", zorder=2)
    ax.plot([r["wl_mm"] for r in front], [r["T_C"] for r in front],
            "-o", color="#1f4e79", ms=5, label="非支配", zorder=3)
    ax.set_xlabel("线长，mm/线")
    ax.set_ylabel("峰值温度，°C")
    ax.set_title("5 芯粒环氧，只动位置", loc="left")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(HERE / "帕累托_温度线长.png", dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    main()
