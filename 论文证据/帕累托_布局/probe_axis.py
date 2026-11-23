#!/usr/bin/env python3
"""沿一条轴探温度和应力的方向：两颗高功耗芯粒从贴合到拉开。

只算 Python 场，不跑 ATSim。Xeon 固定在左下，GPU 从 Xeon 右边一格开始
往右挪，其余三颗放在右上角，离热源最远。
"""

from __future__ import annotations

import json
import sys
from copy import deepcopy
from pathlib import Path

SRC = Path("/home/czh/任务仿真产生功耗文件工具/热仿真流水线/src")
sys.path.insert(0, str(SRC))

import numpy as np
from scipy.sparse.linalg import splu

import run_closed_loop_multidie as cl
import run_fill_superposition as hs
import run_hw_area_field as hw

HERE = Path("/home/czh/任务仿真产生功耗文件工具/论文证据/帕累托_布局")
BASE = json.loads((HERE / "runs" / "layouts" / "init.json").read_text(encoding="utf-8"))


def place(gap: int) -> list[dict] | None:
    dies = deepcopy(BASE)
    by = {d["slot"]: d for d in dies}
    cpu, gpu = by["Ca"], by["Ga"]
    cpu["j0"], cpu["ib0"] = 0, 0
    gpu["j0"] = cpu["w"] + gap
    gpu["ib0"] = 0
    if gpu["j0"] + gpu["w"] > hs.N:
        return None
    # 冷的三颗堆到右上角，远离两颗热源。
    x = hs.N
    for slot in ("HBM3", "CIM_jia_m64", "CIM_colonnade_m512"):
        die = by[slot]
        x -= die["w"]
        die["j0"] = x
        die["ib0"] = hs.N - die["h"]
        x -= 1
    return dies


def main() -> None:
    dies0 = deepcopy(BASE)
    k_si = hw.paint(dies0, hs.K_SI)
    mat, g_amb = hs.build_operator(k_si)
    fac = splu(mat.tocsc())
    print(f"{'gap':5} {'Tpy':8} {'g_metric':8} {'g_max':8} {'g_p99':8}")
    for gap in range(0, 30, 2):
        dies = place(gap)
        if dies is None:
            print(f"{gap:5} 超出板边")
            break
        cl.apply_power(dies)
        temp = hw.solve(k_si, hw.power_map(dies), fac, g_amb)[0]
        ev = cl.metrics(temp)
        g = np.gradient(temp)
        gmag = np.hypot(g[0], g[1]) / cl.CW_MM
        print(f"{gap:5} {ev['Tmax_C']:8.2f} {ev['grad_max']:8.3f} "
              f"{float(gmag.max()):8.3f} {float(np.percentile(gmag, 99)):8.3f}", flush=True)


if __name__ == "__main__":
    main()
