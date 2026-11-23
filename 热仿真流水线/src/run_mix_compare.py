#!/usr/bin/env python3
"""四种动作占比，在现在的芯片案例上比。

同一颗种子先生成一次初始布局，四种占比从这一份出发。
代价沿用闭环的三项各三分之一，越小越好。画出来的是相对起点
下降的比例，所以曲线越高越好。

前 60% 步和后 40% 步分开记，用来看哪段该用哪种占比。
不跑 ATSim，不改已有的案例。
"""

from __future__ import annotations

import json
import random
import sys
from copy import deepcopy
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import font_manager
from scipy.sparse.linalg import splu

ROOT = Path("/home/czh/任务仿真产生功耗文件工具")
SRC = ROOT / "热仿真流水线" / "src"
sys.path.insert(0, str(SRC))

import run_closed_loop_multidie as cl
import run_closed_loop_suite as suite
import run_fill_superposition as hs
import run_hw_area_field as hw
import run_sa_hw_smooth as sm
import run_wirelength as wire

OUT = ROOT / "热仿真流水线" / "runs" / "mix_compare"
STEPS = 60
SEEDS = (11, 22, 33)
KINDS = ("n5", "n8", "n11")
# (微调, 大跳, 交换位置, 旋转)
MIXES = {
    "多跳": (0.10, 0.60, 0.20, 0.10),
    "多交换": (0.20, 0.20, 0.50, 0.10),
    "传统": (0.40, 0.25, 0.20, 0.15),
    "多微调": (0.70, 0.10, 0.10, 0.10),
}
COLORS = {"多跳": "#b85c38", "多交换": "#1f4e79", "传统": "#888780", "多微调": "#2f6f4e"}


def price(dies: list[dict], fac, g_amb, matrix: list[list[int]], scales: dict) -> tuple[float, dict]:
    cl.apply_power(dies)
    field = hw.solve(np.full((hs.N, hs.N), hs.K_SI), hw.power_map(dies), fac, g_amb)[0]
    ev = cl.metrics(field)
    wl = wire.wirelength_cells(dies, matrix, cl.CW_MM)
    return suite.cost_of(ev, wl, scales), ev | {"wl": wl}


def search(start: list[dict], cat: dict, matrix: list[list[int]], scales: dict,
           fac, g_amb, mix: tuple[float, float, float, float], seed: int) -> dict:
    rng_m = random.Random(seed)
    rng_a = random.Random(seed + 7919)
    current = deepcopy(start)
    j_cur, ev = price(current, fac, g_amb, matrix, scales)
    best_j = j_cur
    curve = [j_cur]
    anneal_t = 0.20
    for _ in range(STEPS):
        trial, _, _ = suite.move_of(current, rng_m, mix)
        if not sm.legal(trial):
            curve.append(best_j)
            anneal_t *= suite.ALPHA
            continue
        j_new, ev_new = price(trial, fac, g_amb, matrix, scales)
        take, _ = sm.accept(j_cur, j_new, anneal_t, rng_a)
        if take:
            current = trial
            j_cur = j_new
            best_j = min(best_j, j_new)
        curve.append(best_j)
        anneal_t *= suite.ALPHA
    return {"curve": curve, "final": best_j}


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    font_manager.fontManager.addfont("/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc")
    plt.rcParams["font.family"] = "Noto Sans CJK JP"
    plt.rcParams["axes.unicode_minus"] = False

    k = np.full((hs.N, hs.N), hs.K_SI)
    mat, g_amb = hs.build_operator(k)
    fac = splu(mat.tocsc())
    cat = cl.catalog()
    records = []
    for kind in KINDS:
        members = suite.members(cat, kind)
        for seed in SEEDS:
            suite.SEED = seed
            start, name, _ = suite.initial_layout(members, hs.K_SI)
            matrix = wire.assumed_matrix(start)
            _, ev0 = price(start, fac, g_amb, matrix, {"dT": 1.0, "g0": 1.0, "wl0": 1.0})
            scales = {"dT": max(ev0["Tmax_C"] - cl.T_AMB_C, 1e-6),
                      "g0": max(ev0["grad_max"], 1e-9), "wl0": max(ev0["wl"], 1e-9)}
            print(f"{kind} seed {seed} init {name} J0 1.000", flush=True)
            for mix_name, mix in MIXES.items():
                got = search(start, cat, matrix, scales, fac, g_amb, mix, seed)
                records.append({"kind": kind, "seed": seed, "init": name, "mix": mix_name,
                                "j0": 1.0, **got})
                print(f"  {mix_name} final {got['final']:.4f}", flush=True)
    (OUT / "records.json").write_text(json.dumps(records) + "\n", encoding="utf-8")
    draw(records)


def draw(records: list[dict]) -> None:
    """曲线是 (J0 - J) / J0，越高越好。"""
    fig, axes = plt.subplots(len(KINDS), len(SEEDS), figsize=(12, 8), sharey=True)
    for row, kind in enumerate(KINDS):
        for col, seed in enumerate(SEEDS):
            ax = axes[row][col]
            for mix_name in MIXES:
                rec = next(r for r in records if r["kind"] == kind and r["seed"] == seed and r["mix"] == mix_name)
                gain = [(1.0 - v) for v in rec["curve"]]
                ax.plot(gain, color=COLORS[mix_name], label=mix_name, lw=1.2)
            ax.axvline(int(0.6 * STEPS), color="#ccc", lw=0.6)
            ax.set_title(f"{kind} 种子 {seed}", loc="left", fontsize=9)
            if row == len(KINDS) - 1:
                ax.set_xlabel("步")
            if col == 0:
                ax.set_ylabel("下降比例")
    axes[0][0].legend(frameon=False, fontsize=7, ncol=2)
    fig.suptitle("同一起点、四种占比。曲线是相对起点下降的比例，越高越好。竖线是 60% 步。", fontsize=11)
    fig.tight_layout()
    fig.savefig(OUT / "mix_curves.png", dpi=130)
    plt.close(fig)
    print(f"wrote {OUT / 'mix_curves.png'}", flush=True)


if __name__ == "__main__":
    main()
