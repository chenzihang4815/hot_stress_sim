#!/usr/bin/env python3
"""Can stored per-chip fields replace a full HotSpot solve inside SA?

The package is fixed at 72 mm, so every legal layout shares one mesh and a
chip's field depends only on its own position. Cache one field per chip per
position. A move then costs one solve for each chip that moved; the rest come
from the cache. Compare against solving the whole layout every step.
"""

from __future__ import annotations

import json
import random
import shutil
import time
from pathlib import Path

import numpy as np

import sa_thermal_feedback as sa
from hotspot_io import run_hotspot, write_hotspot_inputs

ROOT = Path("/home/czh/任务仿真产生功耗文件工具")
HS = ROOT / "HotSpot-master" / "hotspot"
BASE = ROOT / "热仿真流水线"
OUT = Path("/tmp/sa_field")
AMB = 318.15
ROWS = 64
INTP = 0.072
STEPS = 40


def grid(work: Path) -> np.ndarray:
    lines = (work / "grid.steady").read_text().splitlines()
    vals = []
    for ln in lines[1:4097]:
        parts = ln.split()
        if len(parts) == 2:
            vals.append(float(parts[1]))
    return np.array(vals)


def solve(chips, tag) -> np.ndarray:
    work = OUT / tag
    if work.exists():
        shutil.rmtree(work)
    write_hotspot_inputs(sa.to_rects(chips), work, rows=ROWS, cols=ROWS)
    run_hotspot(HS, work)
    return grid(work)


def solo(chips, i) -> list[sa.Chip]:
    return [sa.Chip(c.name, c.kind, c.x, c.y, c.w, c.h,
                    chips[i].power_W if k == i else 0.0, c.area_mm2)
            for k, c in enumerate(chips)]


def run_case(label, path) -> dict:
    chips = [sa.Chip(c["name"], c["kind"], c["x_m"], c["y_m"], c["w_m"], c["h_m"],
                     c["power_W"], c["area_mm2"])
             for c in json.loads((BASE / path).read_text())]
    n = len(chips)
    if not sa.legal(chips, INTP, INTP):
        chips = sa.spread_to_gap(chips, INTP, INTP)
        if not sa.legal(chips, INTP, INTP):
            print(f"{label} 无法在 72 mm 内排开, 跳过", flush=True)
            return None
        print(f"{label} 初始布局间隙不足, 已推开", flush=True)
    global OUT
    OUT = Path("/tmp/sa_field") / label
    OUT.mkdir(parents=True, exist_ok=True)

    # Precompute: one field per chip at its starting position.
    t0 = time.time()
    cache = {}
    for i in range(n):
        key = (i, round(chips[i].x, 6), round(chips[i].y, 6), round(chips[i].w, 6), round(chips[i].h, 6))
        cache[key] = solve(solo(chips, i), f"pre{i}")
    t_pre = time.time() - t0
    print(f"{label} 预计算 {n} 个场 {t_pre:.1f}s", flush=True)

    rng = random.Random(22)
    cur = chips
    t_full = t_field = 0.0
    solved = cached = 0
    errs = []
    for step in range(1, STEPS + 1):
        trial, op = sa.propose(cur, INTP, INTP, rng)
        tries = 1
        while not sa.legal(trial, INTP, INTP) and tries < 25:
            trial, op = sa.propose(cur, INTP, INTP, rng)
            tries += 1
        if not sa.legal(trial, INTP, INTP):
            continue

        t0 = time.time()
        full = solve(trial, f"full{step}")
        t_full += time.time() - t0

        moved = [i for i in range(n)
                 if abs(trial[i].x - cur[i].x) > 1e-9 or abs(trial[i].y - cur[i].y) > 1e-9
                 or abs(trial[i].w - cur[i].w) > 1e-9 or abs(trial[i].h - cur[i].h) > 1e-9]
        # The key is the whole layout, not the one chip. propose() snaps every
        # chip back to the origin, so a chip that did not move still shifts,
        # and a field computed on the old grid cannot be added to the new one.
        t0 = time.time()
        acc = np.zeros_like(full)
        for i in range(n):
            key = (i, round(trial[i].x, 6), round(trial[i].y, 6),
                   round(trial[i].w, 6), round(trial[i].h, 6))
            if key not in cache:
                cache[key] = solve(solo(trial, i), f"s{step}_{i}")
                solved += 1
            else:
                cached += 1
            acc += cache[key]
        acc -= (n - 1) * AMB
        t_field += time.time() - t0
        errs.append(float(np.max(np.abs(acc - full))))
        cur = trial
        if step % 10 == 0:
            print(f"  {label} 步 {step} 误差 {errs[-1]:.4f} K", flush=True)

    rec = {"label": label, "chips": n, "steps": len(errs), "pre_s": round(t_pre, 1),
           "full_s": round(t_full, 1), "field_s": round(t_field, 1),
           "solved": solved, "cached": cached,
           "err_max_K": round(max(errs), 4), "err_median_K": round(sorted(errs)[len(errs)//2], 4)}
    print(json.dumps(rec, ensure_ascii=False), flush=True)
    return rec


CASES = [
    ("C1_4芯", "runs/sa_paper_cases/C1_1T2_1GPU/field_old/layout.json"),
    ("big_6芯", "runs/big_balanced_2/lbfgs_init.json"),
    ("hbm_7芯", "runs/with_hbm_2/Lstar.json"),
    ("hbm2_7芯", "runs/with_hbm_2_cpu_plus_gpu/lbfgs_init.json"),
    ("t2_9芯", "runs/full_chain/pkg_2t2_2gpu/Lstar.json"),
    ("mix_9芯", "runs/full_chain/pkg_1xeon_1t2_2gpu/Lstar.json"),
    ("xeon2_9芯", "runs/full_chain/pkg_2xeon_2gpu/Lstar.json"),
    ("xeon3g_10芯", "runs/full_chain/pkg_2xeon_3gpu/Lstar.json"),
    ("xeon3_11芯", "runs/full_chain/pkg_3xeon_3gpu/Lstar.json"),
]


def main() -> None:
    rows = [r for r in (run_case(label, path) for label, path in CASES) if r]
    out = BASE / "runs/superposition/sa_field_cases.json"
    out.write_text(json.dumps(rows, indent=2, ensure_ascii=False))
    print("saved", out, flush=True)


if __name__ == "__main__":
    main()
