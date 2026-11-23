#!/usr/bin/env python3
"""Record the best temperature at every step for the chosen search, and save
the 32-grid and 64-grid thermal maps of one layout."""

from __future__ import annotations

import importlib.util
import json
import math
import random
import time
from copy import deepcopy
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
import numpy as np

import sa_thermal_feedback as sa

ROOT = Path("/home/czh/任务仿真产生功耗文件工具")
SRC = ROOT / "热仿真流水线" / "src"
BASE = ROOT / "热仿真流水线"
OUT = BASE / "runs" / "superposition" / "chosen"
spec = importlib.util.spec_from_file_location("f", SRC / "run_sa_field_120.py")
f = importlib.util.module_from_spec(spec)
spec.loader.exec_module(f)
f.OUT = OUT

INTP = 0.072
STEPS = 120
SEED = 22
JUMP = (0.10, 0.60, 0.20, 0.10)
UNIFORM = (0.40, 0.25, 0.20, 0.15)

font_manager.fontManager.addfont("/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc")
plt.rcParams["font.family"] = "Noto Sans CJK JP"
plt.rcParams["axes.unicode_minus"] = False


def trace(chips, mix, coarse_frac):
    rng = random.Random(SEED)
    cache32, cache64 = {}, {}
    n_coarse = round(STEPS * coarse_frac)
    cur = best = deepcopy(chips)
    f.ROWS = 64
    tmax, _, _ = f.field_temp(cur, cache64, "start")
    tmax_cur = tmax_best = tmax
    T = sa.T0
    best_at = [tmax]
    for step in range(1, STEPS + 1):
        coarse = step <= n_coarse
        trial, op = sa.propose(cur, INTP, INTP, rng, mix)
        tries = 1
        while not sa.legal(trial, INTP, INTP) and tries < 25:
            trial, op = sa.propose(cur, INTP, INTP, rng, mix)
            tries += 1
        if not sa.legal(trial, INTP, INTP):
            T *= sa.ALPHA
            best_at.append(tmax_best)
            continue
        f.ROWS = 32 if coarse else 64
        tmax_new, _, _ = f.field_temp(trial, cache32 if coarse else cache64, f"s{step}")
        j_cur, j_new = (tmax_cur - 45) / 25, (tmax_new - 45) / 25
        if rng.random() < (1.0 if j_new <= j_cur else math.exp((j_cur - j_new) / max(T, 1e-9))):
            cur, tmax_cur = trial, tmax_new
            if tmax_new < tmax_best:
                best, tmax_best = deepcopy(trial), tmax_new
        T *= sa.ALPHA
        best_at.append(tmax_best)
    return best, best_at


def grid_map(chips, rows):
    f.ROWS = rows
    work_tag = f"map{rows}"
    g = f.solve(chips, work_tag)
    return g


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    chips = [sa.Chip(c["name"], c["kind"], c["x_m"], c["y_m"], c["w_m"], c["h_m"],
                     c["power_W"], c["area_mm2"])
             for c in json.loads((BASE / "runs/with_hbm_2/Lstar.json").read_text())]
    best, best_at = trace(chips, JUMP, 0.984)
    (OUT / "trace.json").write_text(json.dumps({"best_C": best_at}))
    print(f"起点 {best_at[0]:.2f} 终点 {best_at[-1]:.2f}", flush=True)

    g64 = grid_map(best, 64)
    g32 = grid_map(best, 32)
    lo = min(g64.min(), g32.min()) - 273.15
    hi = max(g64.max(), g32.max()) - 273.15
    fig, axes = plt.subplots(1, 2, figsize=(11.6, 5.6), constrained_layout=True)
    im = None
    for ax, g, rows, title in ((axes[0], g64, 64, "64 网格"), (axes[1], g32, 32, "32 网格")):
        im = ax.imshow(g.reshape(rows, rows)[::-1] - 273.15, origin="lower",
                       extent=(0, 72, 0, 72), vmin=lo, vmax=hi, cmap="inferno")
        ax.set_title(f"{title}\n峰值 {float(g.max()-273.15):.2f} °C", fontsize=13)
        ax.set_xlabel("x / mm")
        ax.set_ylabel("y / mm")
    fig.colorbar(im, ax=axes, fraction=0.046, pad=0.02).set_label("温度 / °C")
    fig.suptitle("同一个最优布局，两种网格", fontsize=14)
    fig.savefig(OUT / "grid_compare.png", dpi=140)
    print("WROTE", OUT / "grid_compare.png", flush=True)
    print(f"64 网格 {float(g64.max()-273.15):.2f} C  32 网格 {float(g32.max()-273.15):.2f} C", flush=True)


if __name__ == "__main__":
    main()
