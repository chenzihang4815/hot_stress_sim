#!/usr/bin/env python3
"""用退火找温度-线长权衡的两个端点和中间的折中点。

目标是 a*温度归一化 + (1-a)*线长归一化，a 从 0 到 1 取 9 个值。
a=0 只压线长，a=1 只压温度。每个 a 跑 300 步，只动位置，不换芯片。
归一化的尺度取自先前 400 个随机布局的最小最大值。
"""

from __future__ import annotations

import json
import random
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
import run_sa_hw_smooth as sm
import run_wirelength as wire

HERE = Path("/home/czh/任务仿真产生功耗文件工具/论文证据/帕累托_n9")
OUT = HERE / "runs"
BASE = json.loads((OUT / "layouts" / "init.json").read_text(encoding="utf-8"))
STEPS = 250
ALPHAS = [0.35, 0.45, 0.55]
# 400 个随机布局的范围，用来把两项拉到同一尺度。
T_LO, T_HI = 80.8, 91.95
W_LO, W_HI = 35.28, 55.84


def move(dies: list[dict], rng: random.Random) -> list[dict]:
    """连试几次，直到走出一个合法的布局。11 颗里大颗多，单次扰动多半重叠。"""
    for _ in range(25):
        trial = deepcopy(dies)
        i = rng.randrange(len(trial))
        roll = rng.random()
        if roll < 0.55:
            step = rng.choice([1, 2, 3, 4, 6, 8])
            dj, di = rng.choice([(1, 0), (-1, 0), (0, 1), (0, -1)])
            trial[i]["j0"] = min(max(trial[i]["j0"] + dj * step, 0), hs.N - trial[i]["w"])
            trial[i]["ib0"] = min(max(trial[i]["ib0"] + di * step, 0), hs.N - trial[i]["h"])
        elif roll < 0.85:
            trial[i]["j0"] = rng.randrange(hs.N - trial[i]["w"] + 1)
            trial[i]["ib0"] = rng.randrange(hs.N - trial[i]["h"] + 1)
        else:
            trial[i]["w"], trial[i]["h"] = trial[i]["h"], trial[i]["w"]
        if sm.legal(trial):
            return trial
    return dies


def norm_cost(t: float, w: float, alpha: float) -> float:
    return alpha * (t - T_LO) / (T_HI - T_LO) + (1 - alpha) * (w - W_LO) / (W_HI - W_LO)


def main() -> None:
    k_si = hw.paint(deepcopy(BASE), hs.K_SI)
    mat, g_amb = hs.build_operator(k_si)
    fac = splu(mat.tocsc())
    matrix = wire.assumed_matrix(BASE)

    def evaluate(dies: list[dict]) -> tuple[float, float]:
        cl.apply_power(dies)
        temp = hw.solve(k_si, hw.power_map(dies), fac, g_amb)[0]
        return cl.metrics(temp)["Tmax_C"], wire.wirelength_cells(dies, matrix, cl.CW_MM)

    rows = []
    for alpha in ALPHAS:
        rng = random.Random(20261117 + int(alpha * 1000))
        current = deepcopy(BASE)
        t, w = evaluate(current)
        j = norm_cost(t, w, alpha)
        best = (j, t, w, deepcopy(current))
        temp_sa = 0.05
        for step in range(STEPS):
            trial = move(current, rng)
            if not sm.legal(trial):
                continue
            tt, ww = evaluate(trial)
            jj = norm_cost(tt, ww, alpha)
            if jj <= j or rng.random() < np.exp((j - jj) / max(temp_sa, 1e-9)):
                current, t, w, j = trial, tt, ww, jj
                if jj < best[0]:
                    best = (jj, tt, ww, deepcopy(trial))
            temp_sa *= 0.985
        _, bt, bw, _ = best
        rows.append({"alpha": alpha, "T_C": bt, "wl_mm": bw})
        print(f"a {alpha:.3f}  T {bt:6.2f}  wl {bw:6.2f}", flush=True)
        (OUT / "anneal_ends.json").write_text(
            json.dumps(rows, indent=1) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
