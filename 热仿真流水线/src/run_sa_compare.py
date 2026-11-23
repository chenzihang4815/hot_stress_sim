#!/usr/bin/env python3
"""Same SA twice: once priced by a full solve, once by cached fields.

Both start from the same layout and the same seed, and both use the fixed
72 mm frame with no snap back to the origin, so the only thing that differs
is where the temperature of a move comes from.
"""

from __future__ import annotations

import importlib.util
import json
import math
import random
import shutil
import time
from copy import deepcopy
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import Rectangle
import numpy as np

import sa_thermal_feedback as sa
from hotspot_io import run_hotspot, write_hotspot_inputs

ROOT = Path("/home/czh/任务仿真产生功耗文件工具")
SRC = ROOT / "热仿真流水线" / "src"
HS = ROOT / "HotSpot-master" / "hotspot"
BASE = ROOT / "热仿真流水线"
OUT = BASE / "runs" / "superposition" / "sa_compare"
AMB = 318.15
ROWS = 64
INTP = 0.072
STEPS = 120
SEED = 22

spec = importlib.util.spec_from_file_location("f", SRC / "run_sa_field_120.py")
f = importlib.util.module_from_spec(spec)
spec.loader.exec_module(f)
f.OUT = OUT

font_manager.fontManager.addfont("/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc")
plt.rcParams["font.family"] = "Noto Sans CJK JP"
plt.rcParams["axes.unicode_minus"] = False


def load():
    chips = [sa.Chip(c["name"], c["kind"], c["x_m"], c["y_m"], c["w_m"], c["h_m"],
                     c["power_W"], c["area_mm2"])
             for c in json.loads((BASE / "runs/with_hbm_2/Lstar.json").read_text())]
    return chips if sa.legal(chips, INTP, INTP) else sa.spread_to_gap(chips, INTP, INTP)


def search(chips, use_field: bool):
    rng = random.Random(SEED)
    cache = {}
    cur = best = deepcopy(chips)
    if use_field:
        tmax, dt, fresh = f.field_temp(cur, cache, "init")
        t_used = dt
        solved = fresh
    else:
        t0 = time.time()
        g = f.solve(cur, "full_init")
        tmax = float(g.max() - 273.15)
        t_used = time.time() - t0
        solved = 1
    tmax_cur = tmax_best = tmax
    T = sa.T0
    accepted = 0
    for step in range(1, STEPS + 1):
        trial, op = sa.propose(cur, INTP, INTP, rng)
        tries = 1
        while not sa.legal(trial, INTP, INTP) and tries < 25:
            trial, op = sa.propose(cur, INTP, INTP, rng)
            tries += 1
        if not sa.legal(trial, INTP, INTP):
            T *= sa.ALPHA
            continue
        if use_field:
            tmax_new, dt, fresh = f.field_temp(trial, cache, f"s{step}")
            t_used += dt
            solved += fresh
        else:
            t0 = time.time()
            g = f.solve(trial, "full_step")
            tmax_new = float(g.max() - 273.15)
            t_used += time.time() - t0
            solved += 1
        j_cur, j_new = (tmax_cur - 45) / 25, (tmax_new - 45) / 25
        if rng.random() < (1.0 if j_new <= j_cur else math.exp((j_cur - j_new) / max(T, 1e-9))):
            cur, tmax_cur = trial, tmax_new
            accepted += 1
            if tmax_new < tmax_best:
                best, tmax_best = deepcopy(trial), tmax_new
        T *= sa.ALPHA
        if step % 30 == 0:
            print(f"  步 {step} 接受 {accepted} 最优 {tmax_best:.2f}", flush=True)
    return best, round(tmax_best, 2), round(t_used, 1), solved, accepted


def render(path, panels):
    lo = min(g.min() for _, g, _ in panels) - 273.15
    hi = max(g.max() for _, g, _ in panels) - 273.15
    fig, axes = plt.subplots(1, 3, figsize=(16.2, 5.8), constrained_layout=True)
    im = None
    for ax, (title, g, lay) in zip(axes, panels):
        im = ax.imshow(g.reshape(64, 64)[::-1] - 273.15, origin="lower",
                       extent=(0, 72, 0, 72), vmin=lo, vmax=hi, cmap="inferno")
        for c in lay:
            ax.add_patch(Rectangle((c.x * 1e3, c.y * 1e3), c.w * 1e3, c.h * 1e3,
                                   fill=False, edgecolor="white", lw=0.8))
            ax.text((c.x + c.w / 2) * 1e3, (c.y + c.h / 2) * 1e3, c.kind,
                    ha="center", va="center", fontsize=7, color="white")
        ax.set_title(f"{title}\n峰值 {float(g.max() - 273.15):.2f} °C", fontsize=12)
        ax.set_xlabel("x / mm")
        ax.set_ylabel("y / mm")
    fig.colorbar(im, ax=axes, fraction=0.046, pad=0.02).set_label("温度 / °C")
    fig.suptitle("同一起点、同一随机种子、120 步，只是每步温度的来源不同", fontsize=14)
    fig.savefig(path, dpi=140)
    print("WROTE", path, flush=True)


def main() -> None:
    if OUT.exists():
        shutil.rmtree(OUT)
    OUT.mkdir(parents=True)
    chips = load()
    print("完整重解驱动", flush=True)
    best_full, t_full, time_full, solves_full, acc_full = search(chips, False)
    print("场叠加驱动", flush=True)
    best_field, t_field, time_field, solves_field, acc_field = search(chips, True)
    g_start = f.solve(chips, "map_start")
    g_full = f.solve(best_full, "map_full")
    g_field = f.solve(best_field, "map_field")
    render(OUT / "sa_full_vs_field.png", [
        ("搜索前", g_start, chips),
        (f"完整重解驱动\n{time_full:.0f}s，接受 {acc_full} 步", g_full, best_full),
        (f"场叠加驱动\n{time_field:.0f}s，接受 {acc_field} 步", g_field, best_field),
    ])
    # Both final layouts re-scored by a full solve, so the peaks are comparable.
    print(f"完整重解驱动 峰值 {float(g_full.max()-273.15):.2f} C  用时 {time_full}s", flush=True)
    print(f"场叠加驱动   峰值 {float(g_field.max()-273.15):.2f} C  用时 {time_field}s", flush=True)
    (OUT / "summary.json").write_text(json.dumps({
        "start_C": round(float(g_start.max()-273.15), 2),
        "full": {"best_C": round(float(g_full.max()-273.15), 2), "secs": time_full,
                 "solves": solves_full, "accepted": acc_full},
        "field": {"best_C": round(float(g_field.max()-273.15), 2), "secs": time_field,
                  "solves": solves_field, "accepted": acc_field},
    }, indent=2))


if __name__ == "__main__":
    main()
