#!/usr/bin/env python3
"""Coarse-then-fine SA against an all-fine SA, both driven by cached fields.

The fine run prices every step on the 64 grid. The coarse run prices the first
part on the 32 grid and the rest on the 64 grid. Both end by re-scoring the
best layout on the 64 grid, so the peaks are comparable. One HotSpot at a time.
"""

from __future__ import annotations

import importlib.util
import json
import math
import random
import time
from copy import deepcopy
from pathlib import Path

import numpy as np

import sa_thermal_feedback as sa

ROOT = Path("/home/czh/任务仿真产生功耗文件工具")
SRC = ROOT / "热仿真流水线" / "src"
BASE = ROOT / "热仿真流水线"
OUT = BASE / "runs" / "superposition" / "coarse_sa"
INTP = 0.072
STEPS = 120
SEED = 22

spec = importlib.util.spec_from_file_location("f", SRC / "run_sa_field_120.py")
f = importlib.util.module_from_spec(spec)
spec.loader.exec_module(f)

CASES = [
    ("C1_4芯", "runs/sa_paper_cases/C1_1T2_1GPU/field_old/layout.json"),
    ("big_6芯", "runs/big_balanced_2/lbfgs_init.json"),
    ("hbm_7芯", "runs/with_hbm_2/Lstar.json"),
    ("t2_9芯", "runs/full_chain/pkg_2t2_2gpu/Lstar.json"),
    ("xeon2_9芯", "runs/full_chain/pkg_2xeon_2gpu/Lstar.json"),
    ("xeon3_11芯", "runs/full_chain/pkg_3xeon_3gpu/Lstar.json"),
]
FRACTIONS = [0.5, 0.75]


def load(path):
    chips = [sa.Chip(c["name"], c["kind"], c["x_m"], c["y_m"], c["w_m"], c["h_m"],
                     c["power_W"], c["area_mm2"]) for c in json.loads((BASE / path).read_text())]
    return chips if sa.legal(chips, INTP, INTP) else sa.spread_to_gap(chips, INTP, INTP)


def run(chips, frac: float):
    """frac of the steps priced on the 32 grid, the rest on 64. Returns the best layout."""
    rng = random.Random(SEED)
    cache32, cache64 = {}, {}
    n_coarse = round(STEPS * frac)
    cur = best = deepcopy(chips)
    f.ROWS = 64
    tmax, _, _ = f.field_temp(cur, cache64, "start")
    tmax_cur = tmax_best = tmax
    T = sa.T0
    t_used = 0.0
    for step in range(1, STEPS + 1):
        trial, op = sa.propose(cur, INTP, INTP, rng)
        tries = 1
        while not sa.legal(trial, INTP, INTP) and tries < 25:
            trial, op = sa.propose(cur, INTP, INTP, rng)
            tries += 1
        if not sa.legal(trial, INTP, INTP):
            T *= sa.ALPHA
            continue
        coarse = step <= n_coarse
        f.ROWS = 32 if coarse else 64
        cache = cache32 if coarse else cache64
        t0 = time.time()
        tmax_new, dt, _ = f.field_temp(trial, cache, f"s{step}")
        t_used += time.time() - t0
        j_cur, j_new = (tmax_cur - 45) / 25, (tmax_new - 45) / 25
        if rng.random() < (1.0 if j_new <= j_cur else math.exp((j_cur - j_new) / max(T, 1e-9))):
            cur, tmax_cur = trial, tmax_new
            if coarse:
                # The coarse temperature only ranks moves. Track the best by it,
                # then re-score that layout on the fine grid at the end.
                if tmax_new < tmax_best:
                    best, tmax_best = deepcopy(trial), tmax_new
            elif tmax_new < tmax_best:
                best, tmax_best = deepcopy(trial), tmax_new
        T *= sa.ALPHA
    f.ROWS = 64
    final, _, _ = f.field_temp(best, {}, "final")
    return best, round(final, 3), round(t_used, 1)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    log = OUT / "results.jsonl"
    for label, path in CASES:
        chips = load(path)
        for frac in [0.0] + FRACTIONS:
            best, final, secs = run(chips, frac)
            rec = {"label": label, "chips": len(chips), "coarse_frac": frac,
                   "final_64_C": final, "secs": secs}
            with log.open("a") as fh:
                fh.write(json.dumps(rec) + "\n")
            print(json.dumps(rec, ensure_ascii=False), flush=True)


if __name__ == "__main__":
    main()
