#!/usr/bin/env python3
"""Traditional SA against the current fast chain, same layouts, same 120 steps.

Traditional: every candidate is a full HotSpot solve on the 64 grid, with the
default move mix (nudge 0.40, jump 0.25, swap 0.20, rotate 0.15).

Fast: every candidate is priced from cached per-chip fields on the 16 grid,
with the big-jump mix (0.10, 0.60, 0.20, 0.10). Its winner is re-scored once
on the 64 grid, which is the number that gets reported.

Both arms share the seed, so the only things that differ are the grid, the
pricing and the move mix. One HotSpot at a time, work dir removed after every
solve, resumable from speedup.jsonl.
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

import numpy as np

ROOT = Path("/home/czh/任务仿真产生功耗文件工具")
SRC = ROOT / "热仿真流水线" / "src"
BASE = ROOT / "热仿真流水线"
HS = ROOT / "HotSpot-master" / "hotspot"

import sa_thermal_feedback as sa  # noqa: E402
from hotspot_io import run_hotspot, write_hotspot_inputs  # noqa: E402

spec = importlib.util.spec_from_file_location("m", SRC / "run_sa_multicase.py")
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)

OUT = BASE / "runs" / "superposition" / "speedup"
WORK = OUT / "_work"
AMB = 318.15
INTP = 0.072
STEPS = 120
SEED = 22
JUMP = (0.10, 0.60, 0.20, 0.10)


def full_peak(chips, rows: int) -> tuple[float, float]:
    """One full solve. Returns (peak_C, seconds) and deletes the work dir."""
    if WORK.exists():
        shutil.rmtree(WORK)
    t0 = time.time()
    write_hotspot_inputs(sa.to_rects(chips), WORK, rows=rows, cols=rows)
    run_hotspot(HS, WORK)
    vals = []
    for line in (WORK / "grid.steady").read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("Layer"):
            continue
        vals.append(float(line.split()[-1]))
        if len(vals) == rows * rows:
            break
    dt = time.time() - t0
    shutil.rmtree(WORK, ignore_errors=True)
    return max(vals) - 273.15, dt


def field_peak(chips, cache, rows: int) -> tuple[float, float, int]:
    """Peak from cached single-chip fields on a fixed grid."""
    n = rows * rows
    acc = None
    fresh = 0
    t0 = time.time()
    for i, c in enumerate(chips):
        key = (i, round(c.x, 6), round(c.y, 6), round(c.w, 6), round(c.h, 6))
        if key not in cache:
            if WORK.exists():
                shutil.rmtree(WORK)
            solo = [sa.Chip(d.name, d.kind, d.x, d.y, d.w, d.h,
                            c.power_W if k == i else 0.0, d.area_mm2)
                    for k, d in enumerate(chips)]
            write_hotspot_inputs(sa.to_rects(solo), WORK, rows=rows, cols=rows)
            run_hotspot(HS, WORK)
            vals = []
            for line in (WORK / "grid.steady").read_text().splitlines():
                line = line.strip()
                if not line or line.startswith("Layer"):
                    continue
                vals.append(float(line.split()[-1]))
                if len(vals) == n:
                    break
            cache[key] = np.array(vals)
            fresh += 1
            shutil.rmtree(WORK, ignore_errors=True)
        acc = cache[key] if acc is None else acc + cache[key]
    acc = acc - (len(chips) - 1) * AMB
    return float(acc.max() - 273.15), time.time() - t0, fresh


def search(chips, rows, mix, use_field):
    """120 steps. Returns (best chips, best coarse peak, solve seconds, solves)."""
    rng = random.Random(SEED)
    cache = {}
    cur = best = deepcopy(chips)
    price = (lambda cs: field_peak(cs, cache, rows)) if use_field else (lambda cs: (*full_peak(cs, rows), 1))
    tmax, dt, fresh = price(cur)
    tmax_cur = tmax_best = tmax
    solved, t_used = fresh, dt
    T = sa.T0
    for step in range(1, STEPS + 1):
        trial, _ = sa.propose(cur, INTP, INTP, rng, mix)
        tries = 1
        while not sa.legal(trial, INTP, INTP) and tries < 25:
            trial, _ = sa.propose(cur, INTP, INTP, rng, mix)
            tries += 1
        if not sa.legal(trial, INTP, INTP):
            T *= sa.ALPHA
            continue
        tmax_new, dt, fresh = price(trial)
        t_used += dt
        solved += fresh
        j_cur, j_new = (tmax_cur - 45) / 25, (tmax_new - 45) / 25
        take = rng.random() < (1.0 if j_new <= j_cur else math.exp((j_cur - j_new) / max(T, 1e-9)))
        if take:
            cur, tmax_cur = trial, tmax_new
            if tmax_new < tmax_best:
                best, tmax_best = deepcopy(trial), tmax_new
        T *= sa.ALPHA
    return best, tmax_best, t_used, solved


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    log = OUT / "speedup.jsonl"
    done = set()
    if log.exists():
        for line in log.read_text().splitlines():
            if line.strip():
                done.add(json.loads(line)["label"])
    for label, path in m.CASES:
        if label in done:
            continue
        chips = m.load(path)
        t0 = time.time()
        best_old, coarse_old, secs_old, solves_old = search(chips, 64, sa.propose.__defaults__[0], False)
        wall_old = time.time() - t0
        final_old, _ = full_peak(best_old, 64)

        t0 = time.time()
        best_new, coarse_new, secs_new, solves_new = search(chips, 16, JUMP, True)
        final_new, resc = full_peak(best_new, 64)
        wall_new = time.time() - t0

        rec = {
            "label": label, "chips": len(chips),
            "trad_final_64_C": round(final_old, 3), "trad_secs": round(secs_old, 1),
            "trad_solves": solves_old, "trad_wall_s": round(wall_old, 1),
            "fast_coarse_C": round(coarse_new, 3), "fast_final_64_C": round(final_new, 3),
            "fast_secs": round(secs_new + resc, 1), "fast_solves": solves_new,
            "fast_wall_s": round(wall_new, 1),
            "speedup": round(secs_old / (secs_new + resc), 2),
            "delta_C": round(final_new - final_old, 3),
        }
        with log.open("a") as fh:
            fh.write(json.dumps(rec, ensure_ascii=False) + "\n")
        print(json.dumps(rec, ensure_ascii=False), flush=True)


if __name__ == "__main__":
    main()
