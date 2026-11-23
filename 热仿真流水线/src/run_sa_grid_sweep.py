#!/usr/bin/env python3
"""Coarse-then-fine SA sweep.

Fixed budget of N_SA steps. The first `frac` of them run on a coarse grid, the
rest on the 64 grid, with the annealing temperature continuing unbroken across
the switch. Every final layout is re-scored on the 64 grid, so all configs are
compared at the same accuracy. One row per run is appended to results.jsonl.
"""

from __future__ import annotations

import json
import math
import random
import shutil
import sys
import time
from copy import deepcopy
from pathlib import Path

ROOT = Path("/home/czh/任务仿真产生功耗文件工具")
SRC = ROOT / "热仿真流水线" / "src"
sys.path.insert(0, str(SRC))

import sa_thermal_feedback as sa  # noqa: E402
from geometry import packing_metrics  # noqa: E402
from hotspot_io import run_hotspot, write_hotspot_inputs  # noqa: E402

HS = Path("/home/czh/任务仿真产生功耗文件工具/HotSpot-master/hotspot")
INIT = ROOT / "热仿真流水线" / "runs" / "with_hbm_2" / "lbfgs_init.json"
OUT = ROOT / "热仿真流水线" / "runs" / "sa_grid_sweep"
FINE = 64
STEPS = 120
SEEDS = [22, 23, 24]
COARSE_GRIDS = [16, 32]
FRACS = [0.25, 0.5, 0.75, 1.0]


def starting_layout():
    chips = [sa.Chip(c["name"], c["kind"], c["x_m"], c["y_m"], c["w_m"], c["h_m"],
                     c["power_W"], c["area_mm2"]) for c in json.loads(INIT.read_text())]
    sa.shift_to_origin(chips)
    m = packing_metrics(sa.to_rects(chips))
    n_gpu = sum(1 for c in chips if c.kind == "GPU")
    w = m["width_m"] + sa.GAP_GPU_M * max(n_gpu, 1) + 0.008
    h = m["height_m"] + sa.GAP_GPU_M * max(n_gpu, 1) + 0.008
    return sa.spread_to_gap(chips, w, h), w, h


def evaluate(chips, rows, work: Path) -> float:
    if work.exists():
        shutil.rmtree(work)
    work.mkdir(parents=True)
    write_hotspot_inputs(sa.to_rects(chips), work, rows=rows, cols=rows)
    return float(run_hotspot(HS, work)["Tmax_C"])


def anneal(start, steps, rng, rows, T, w, h, work: Path):
    cur = best = deepcopy(start)
    t_cur = t_best = evaluate(cur, rows, work)
    evals = accepted = 1
    for _ in range(steps):
        trial, _ = sa.propose(cur, w, h, rng)
        if not sa.legal(trial, w, h):
            T *= sa.ALPHA
            continue
        try:
            t_new = evaluate(trial, rows, work)
        except Exception:
            T *= sa.ALPHA
            continue
        evals += 1
        j_cur, j_new = (t_cur - 45) / 25, (t_new - 45) / 25
        if rng.random() < (1.0 if j_new <= j_cur else math.exp((j_cur - j_new) / max(T, 1e-9))):
            cur, t_cur = trial, t_new
            accepted += 1
            if t_new < t_best - 1e-12:
                best, t_best = deepcopy(trial), t_new
        T *= sa.ALPHA
    return best, t_best, evals, accepted, T


def run_one(coarse, frac, seed, base, w, h, work: Path) -> dict:
    rng = random.Random(seed)
    n_coarse = round(STEPS * frac)
    t0 = time.time()
    lay, t_coarse, e1, a1, T = anneal(base, n_coarse, rng, coarse, sa.T0, w, h, work)
    lay, t_fine, e2, a2, T = anneal(lay, STEPS - n_coarse, rng, FINE, T, w, h, work)
    t_final = evaluate(lay, FINE, work / "final")
    return {
        "coarse_grid": coarse, "frac_coarse": frac, "seed": seed,
        "secs": round(time.time() - t0, 1),
        "evals_coarse": e1, "evals_fine": e2,
        "best_native_C": round(t_fine if n_coarse < STEPS else t_coarse, 3),
        "final_64_C": round(t_final, 3),
    }


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    log = OUT / "results.jsonl"
    base, w, h = starting_layout()
    jobs = [(FINE, 0.0, s) for s in SEEDS]
    jobs += [(g, f, s) for g in COARSE_GRIDS for f in FRACS for s in SEEDS]
    for i, (g, f, s) in enumerate(jobs, 1):
        row = run_one(g, f, s, base, w, h, OUT / "work")
        with log.open("a") as fh:
            fh.write(json.dumps(row) + "\n")
        print(f"[{i}/{len(jobs)}] {row}", flush=True)


if __name__ == "__main__":
    main()
