#!/usr/bin/env python3
"""Move-mix sweep on one layout, coarse for 75% then fine, field driven."""

from __future__ import annotations

import importlib.util
import json
import math
import random
import time
from copy import deepcopy
from pathlib import Path

import sa_thermal_feedback as sa

ROOT = Path("/home/czh/任务仿真产生功耗文件工具")
SRC = ROOT / "热仿真流水线" / "src"
BASE = ROOT / "热仿真流水线"
spec = importlib.util.spec_from_file_location("f", SRC / "run_sa_field_120.py")
f = importlib.util.module_from_spec(spec)
spec.loader.exec_module(f)
spec = importlib.util.spec_from_file_location("c", SRC / "run_sa_coarse.py")
c = importlib.util.module_from_spec(spec)
spec.loader.exec_module(c)

INTP = 0.072
STEPS = 120
SEED = 22
UNIFORM = (0.40, 0.25, 0.20, 0.15)
MIXES = {
    "传统 40/25/20/15": UNIFORM,
    "粗段多跳 10/60/20/10": (0.10, 0.60, 0.20, 0.10),
    "粗段多微调 70/10/10/10": (0.70, 0.10, 0.10, 0.10),
    "粗段多交换 20/20/50/10": (0.20, 0.20, 0.50, 0.10),
}


def run(chips, mix):
    rng = random.Random(SEED)
    cache32, cache64 = {}, {}
    n_coarse = round(STEPS * 0.75)
    cur = best = deepcopy(chips)
    f.ROWS = 64
    tmax, _, _ = f.field_temp(cur, cache64, "start")
    tmax_cur = tmax_best = tmax
    T = sa.T0
    t_used = 0.0
    drawn = {}
    for step in range(1, STEPS + 1):
        trial, op = sa.propose(cur, INTP, INTP, rng, mix)
        tries = 1
        while not sa.legal(trial, INTP, INTP) and tries < 25:
            trial, op = sa.propose(cur, INTP, INTP, rng, mix)
            tries += 1
        if not sa.legal(trial, INTP, INTP):
            T *= sa.ALPHA
            continue
        kind = op.split()[0]
        drawn[kind] = drawn.get(kind, 0) + 1
        coarse = step <= n_coarse
        f.ROWS = 32 if coarse else 64
        t0 = time.time()
        tmax_new, _, _ = f.field_temp(trial, cache32 if coarse else cache64, f"s{step}")
        t_used += time.time() - t0
        j_cur, j_new = (tmax_cur - 45) / 25, (tmax_new - 45) / 25
        if rng.random() < (1.0 if j_new <= j_cur else math.exp((j_cur - j_new) / max(T, 1e-9))):
            cur, tmax_cur = trial, tmax_new
            if tmax_new < tmax_best:
                best, tmax_best = deepcopy(trial), tmax_new
        T *= sa.ALPHA
    f.ROWS = 64
    final, _, _ = f.field_temp(best, {}, "final")
    return round(final, 3), round(t_used, 1), drawn


def main() -> None:
    chips = c.load("runs/with_hbm_2/Lstar.json")
    for name, mix in MIXES.items():
        final, secs, drawn = run(chips, mix)
        print(json.dumps({"mix": name, "final_64_C": final, "secs": secs, "moves": drawn},
                         ensure_ascii=False), flush=True)


if __name__ == "__main__":
    main()
