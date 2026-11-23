#!/usr/bin/env python3
"""Fine-phase move mix, with the coarse mix held fixed. 7-chip layout."""

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
JUMP = (0.10, 0.60, 0.20, 0.10)
NUDGE = (0.70, 0.10, 0.10, 0.10)
SWAP = (0.20, 0.20, 0.50, 0.10)
FINE = {"传统": UNIFORM, "多跳": JUMP, "多微调": NUDGE, "多交换": SWAP}


def run(chips, coarse_mix, fine_mix):
    rng = random.Random(SEED)
    cache32, cache64 = {}, {}
    n_coarse = round(STEPS * 0.75)
    cur = best = deepcopy(chips)
    f.ROWS = 64
    tmax, _, _ = f.field_temp(cur, cache64, "start")
    tmax_cur = tmax_best = tmax
    T = sa.T0
    t_used = 0.0
    for step in range(1, STEPS + 1):
        coarse = step <= n_coarse
        trial, op = sa.propose(cur, INTP, INTP, rng, coarse_mix if coarse else fine_mix)
        tries = 1
        while not sa.legal(trial, INTP, INTP) and tries < 25:
            trial, op = sa.propose(cur, INTP, INTP, rng, coarse_mix if coarse else fine_mix)
            tries += 1
        if not sa.legal(trial, INTP, INTP):
            T *= sa.ALPHA
            continue
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
    return round(final, 3), round(t_used, 1)


def main() -> None:
    chips = c.load("runs/with_hbm_2/Lstar.json")
    for cname, cmix in (("传统", UNIFORM), ("多跳", JUMP)):
        for fname, fmix in FINE.items():
            final, secs = run(chips, cmix, fmix)
            print(json.dumps({"粗阶段": cname, "细阶段": fname, "final_64_C": final, "secs": secs},
                             ensure_ascii=False), flush=True)


if __name__ == "__main__":
    main()
