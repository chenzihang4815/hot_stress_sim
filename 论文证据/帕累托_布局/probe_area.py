#!/usr/bin/env python3
"""沿芯粒面积探温度方向。功耗固定，Xeon 和 GPU 的面积一起缩放。

只算 Python 场。两颗热源固定贴在左下角，其余三颗固定在右上角，
所以变的只有功率密度。面积从 0.4 倍到 1.6 倍，步长 0.1。
"""

from __future__ import annotations

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
import json
BASE = json.loads((HERE / "runs" / "layouts" / "init.json").read_text(encoding="utf-8"))


def resized(scale: float) -> list[dict] | None:
    dies = deepcopy(BASE)
    by = {d["slot"]: d for d in dies}
    for slot in ("Ca", "Ga"):
        die = by[slot]
        w = max(2, round(die["w"] * scale ** 0.5))
        h = max(2, round(die["h"] * scale ** 0.5))
        die["w"], die["h"] = w, h
    cpu, gpu = by["Ca"], by["Ga"]
    cpu["j0"], cpu["ib0"] = 0, 0
    gpu["j0"], gpu["ib0"] = cpu["w"] + 1, 0
    if gpu["j0"] + gpu["w"] > hs.N:
        return None
    x = hs.N
    for slot in ("HBM3", "CIM_jia_m64", "CIM_colonnade_m512"):
        die = by[slot]
        x -= die["w"]
        if x < gpu["j0"] + gpu["w"]:
            return None
        die["j0"] = x
        die["ib0"] = hs.N - die["h"]
        x -= 1
    return dies


def main() -> None:
    k_si = hw.paint(deepcopy(BASE), hs.K_SI)
    mat, g_amb = hs.build_operator(k_si)
    fac = splu(mat.tocsc())
    print(f"{'scale':7} {'Xeon mm2':9} {'GPU mm2':9} {'Tpy':8}")
    for k in range(4, 17):
        scale = k / 10
        dies = resized(scale)
        if dies is None:
            print(f"{scale:7.1f} 超出板边")
            continue
        cl.apply_power(dies)
        temp = hw.solve(k_si, hw.power_map(dies), fac, g_amb)[0]
        ev = cl.metrics(temp)
        by = {d["slot"]: d for d in dies}
        print(f"{scale:7.1f} {by['Ca']['w']*by['Ca']['h']*1.2656:9.1f} "
              f"{by['Ga']['w']*by['Ga']['h']*1.2656:9.1f} {ev['Tmax_C']:8.2f}", flush=True)


if __name__ == "__main__":
    main()
