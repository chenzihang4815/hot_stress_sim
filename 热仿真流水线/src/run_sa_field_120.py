#!/usr/bin/env python3
"""120-step SA timed both ways, plus a before/after thermal map.

The search itself is driven by the cached fields. Every accepted-candidate
step is also solved in full, so the time and the temperature can be compared
on the same sequence of moves.
"""

from __future__ import annotations

import json
import math
import random
import shutil
import time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
import numpy as np

import sa_thermal_feedback as sa
from hotspot_io import run_hotspot, write_hotspot_inputs

ROOT = Path("/home/czh/任务仿真产生功耗文件工具")
HS = ROOT / "HotSpot-master" / "hotspot"
BASE = ROOT / "热仿真流水线"
OUT = BASE / "runs" / "superposition" / "sa120"
AMB = 318.15
ROWS = 64
INTP = 0.072
STEPS = 120
SEED = 22

font_manager.fontManager.addfont("/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc")
plt.rcParams["font.family"] = "Noto Sans CJK JP"
plt.rcParams["axes.unicode_minus"] = False


def grid_of(work: Path) -> np.ndarray:
    vals = []
    for line in (work / "grid.steady").read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("Layer"):
            continue
        vals.append(float(line.split()[-1]))
        if len(vals) == 4096:
            break
    return np.array(vals)


def solve(chips, tag) -> np.ndarray:
    work = OUT / tag
    if work.exists():
        shutil.rmtree(work)
    write_hotspot_inputs(sa.to_rects(chips), work, rows=ROWS, cols=ROWS)
    run_hotspot(HS, work)
    return grid_of(work)


def solo(chips, i):
    return [sa.Chip(c.name, c.kind, c.x, c.y, c.w, c.h,
                    chips[i].power_W if k == i else 0.0, c.area_mm2)
            for k, c in enumerate(chips)]


def field_temp(chips, cache, tag):
    """Peak temperature from cached per-chip fields. Returns (Tmax_C, seconds, new solves)."""
    t0 = time.time()
    acc = None
    fresh = 0
    for i in range(len(chips)):
        key = (i, round(chips[i].x, 6), round(chips[i].y, 6),
               round(chips[i].w, 6), round(chips[i].h, 6))
        if key not in cache:
            cache[key] = solve(solo(chips, i), f"{tag}_{i}")
            fresh += 1
        acc = cache[key] if acc is None else acc + cache[key]
    acc = acc - (len(chips) - 1) * AMB
    return float(acc.max() - 273.15), time.time() - t0, fresh


def main() -> None:
    chips = [sa.Chip(c["name"], c["kind"], c["x_m"], c["y_m"], c["w_m"], c["h_m"],
                     c["power_W"], c["area_mm2"])
             for c in json.loads((BASE / "runs/with_hbm_2/Lstar.json").read_text())]
    if not sa.legal(chips, INTP, INTP):
        chips = sa.spread_to_gap(chips, INTP, INTP)
    if OUT.exists():
        shutil.rmtree(OUT)
    OUT.mkdir(parents=True)

    cache = {}
    t_field = t_full = 0.0
    solved = 0
    tmax, dt, fresh = field_temp(chips, cache, "init")
    t_field += dt
    solved += fresh
    g0 = solve(chips, "map_before")
    t_full += time.time() - time.time()  # placeholder, replaced below

    # Time the initial full solve separately so the map is not double counted.
    t0 = time.time()
    g0 = solve(chips, "map_before")
    t_init_full = time.time() - t0

    rng = random.Random(SEED)
    cur = chips
    tmax_cur = tmax
    best = chips
    tmax_best = tmax
    T = sa.T0
    accepted = 0
    errs = []
    for step in range(1, STEPS + 1):
        trial, op = sa.propose(cur, INTP, INTP, rng)
        tries = 1
        while not sa.legal(trial, INTP, INTP) and tries < 25:
            trial, op = sa.propose(cur, INTP, INTP, rng)
            tries += 1
        if not sa.legal(trial, INTP, INTP):
            T *= sa.ALPHA
            continue
        tmax_new, dt, fresh = field_temp(trial, cache, f"s{step}")
        t_field += dt
        solved += fresh
        t0 = time.time()
        full = solve(trial, "check")
        t_full += time.time() - t0
        errs.append(abs(tmax_new - float(full.max() - 273.15)))
        j_cur, j_new = (tmax_cur - 45) / 25, (tmax_new - 45) / 25
        take = rng.random() < (1.0 if j_new <= j_cur else math.exp((j_cur - j_new) / max(T, 1e-9)))
        if take:
            cur, tmax_cur = trial, tmax_new
            accepted += 1
            if tmax_new < tmax_best:
                best, tmax_best = trial, tmax_new
        T *= sa.ALPHA
        if step % 20 == 0:
            print(f"步 {step} 接受 {accepted} 当前 {tmax_cur:.2f} 最优 {tmax_best:.2f} "
                  f"误差 {errs[-1]:.3f} K", flush=True)

    g1 = solve(best, "map_after")
    print(f"\n场叠加 {t_field:.1f}s  完整重解 {t_init_full + t_full:.1f}s  "
          f"新求解 {solved} 次  接受 {accepted} 步", flush=True)
    print(f"峰值误差 最大 {max(errs):.3f} K  中位 {sorted(errs)[len(errs)//2]:.3f} K", flush=True)
    print(f"起点 {float(g0.max()-273.15):.2f} C  最优 {float(g1.max()-273.15):.2f} C", flush=True)

    lo = min(g0.min(), g1.min()) - 273.15
    hi = max(g0.max(), g1.max()) - 273.15
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 5.6), constrained_layout=True)
    im = None
    for ax, g, title in ((axes[0], g0, "搜索前"), (axes[1], g1, "搜索后（120 步最优）")):
        im = ax.imshow(g.reshape(64, 64)[::-1] - 273.15, origin="lower",
                       extent=(0, 72, 0, 72), vmin=lo, vmax=hi, cmap="inferno")
        ax.set_title(f"{title}\n峰值 {float(g.max()-273.15):.2f} °C", fontsize=13)
        ax.set_xlabel("x / mm")
        ax.set_ylabel("y / mm")
    fig.colorbar(im, ax=axes, fraction=0.046, pad=0.02).set_label("温度 / °C")
    fig.suptitle("7 芯布局，固定 72 mm 芯片层，场叠加驱动的 120 步 SA", fontsize=14)
    png = OUT / "sa120_before_after.png"
    fig.savefig(png, dpi=140)
    print("WROTE", png, flush=True)
    (OUT / "summary.json").write_text(json.dumps({
        "steps": STEPS, "accepted": accepted, "field_s": round(t_field, 1),
        "full_s": round(t_init_full + t_full, 1), "solved": solved,
        "err_max_K": round(max(errs), 3), "start_C": round(float(g0.max()-273.15), 2),
        "best_C": round(float(g1.max()-273.15), 2)}, indent=2))


if __name__ == "__main__":
    main()
