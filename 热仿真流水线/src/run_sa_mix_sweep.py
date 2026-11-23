#!/usr/bin/env python3
"""Stage-dependent move mix vs the uniform mix.

Coarse fraction is fixed at 0.75 (chosen from the grid sweep). The coarse stage
runs on a 16 or 32 grid, the fine stage on 64, and every final layout is
re-scored on 64. Two mixes:

- uniform: 40% nudge, 25% jump, 20% swap, 15% rotate, both stages. This is what
  sa.propose does, so these rows reproduce the grid sweep and are the baseline.
- staged: coarse stage jumps more (20/45/25/10), fine stage nudges more
  (80/5/5/10).

One row per run is appended to results.jsonl.
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
OUT = ROOT / "热仿真流水线" / "runs" / "sa_mix_sweep"
FINE = 64
STEPS = 120
FRAC = 0.75
SEEDS = [22, 23, 24]
COARSE_GRIDS = [16, 32]

# (nudge, jump, swap, rotate). Cumulative edges are computed below.
UNIFORM = (0.40, 0.25, 0.20, 0.15)
STAGED_COARSE = (0.20, 0.45, 0.25, 0.10)
STAGED_FINE = (0.80, 0.05, 0.05, 0.10)


def starting_layout():
    chips = [sa.Chip(c["name"], c["kind"], c["x_m"], c["y_m"], c["w_m"], c["h_m"],
                     c["power_W"], c["area_mm2"]) for c in json.loads(INIT.read_text())]
    sa.shift_to_origin(chips)
    m = packing_metrics(sa.to_rects(chips))
    n_gpu = sum(1 for c in chips if c.kind == "GPU")
    w = m["width_m"] + sa.GAP_GPU_M * max(n_gpu, 1) + 0.008
    h = m["height_m"] + sa.GAP_GPU_M * max(n_gpu, 1) + 0.008
    return sa.spread_to_gap(chips, w, h), w, h


def propose(chips, intp_w, intp_h, rng, mix):
    """Same four moves as sa.propose, only the probabilities change."""
    nudge, jump, swap, _rotate = mix
    e_nudge, e_jump, e_swap = nudge, nudge + jump, nudge + jump + swap
    n = len(chips)
    trial = deepcopy(chips)
    gpus = [k for k, c in enumerate(chips) if c.kind == "GPU"]
    r = rng.random()
    if r < e_nudge:
        i = rng.choice(gpus) if gpus and rng.random() < 0.55 else rng.randrange(n)
        steps = rng.choice([1, 1, 2, 2, 3, 4])
        dx, dy = rng.choice([(sa.GRAN_M, 0), (-sa.GRAN_M, 0), (0, sa.GRAN_M), (0, -sa.GRAN_M)])
        trial[i].x += dx * steps
        trial[i].y += dy * steps
        op = "nudge"
    elif r < e_jump:
        i = rng.choice(gpus) if gpus and rng.random() < 0.6 else rng.randrange(n)
        max_x = max(0, int((intp_w - trial[i].w) / sa.GRAN_M))
        max_y = max(0, int((intp_h - trial[i].h) / sa.GRAN_M))
        trial[i].x = rng.randint(0, max_x) * sa.GRAN_M
        trial[i].y = rng.randint(0, max_y) * sa.GRAN_M
        if rng.random() < 0.4:
            trial[i].rotate90()
            op = "jump+rot"
        else:
            op = "jump"
    elif r < e_swap:
        i, j = rng.sample(range(n), 2)
        trial[i].x, trial[j].x = trial[j].x, trial[i].x
        trial[i].y, trial[j].y = trial[j].y, trial[i].y
        op = "swap"
    else:
        i = rng.randrange(n)
        trial[i].rotate90()
        op = "rotate"
    sa.shift_to_origin(trial)
    return trial, op


def evaluate(chips, rows, work: Path) -> float:
    if work.exists():
        shutil.rmtree(work)
    work.mkdir(parents=True)
    write_hotspot_inputs(sa.to_rects(chips), work, rows=rows, cols=rows)
    return float(run_hotspot(HS, work)["Tmax_C"])


def anneal(start, steps, rng, rows, mix, T, w, h, work: Path):
    cur = best = deepcopy(start)
    t_cur = t_best = evaluate(cur, rows, work)
    evals = accepted = 1
    ops = {"nudge": 0, "jump": 0, "jump+rot": 0, "swap": 0, "rotate": 0}
    for _ in range(steps):
        trial, op = propose(cur, w, h, rng, mix)
        ops[op] = ops.get(op, 0) + 1
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
    return best, t_best, evals, accepted, T, ops


def run_one(coarse, mix_name, seed, base, w, h, work: Path) -> dict:
    rng = random.Random(seed)
    n_coarse = round(STEPS * FRAC)
    if mix_name == "uniform":
        mix_c = mix_f = UNIFORM
    else:
        mix_c, mix_f = STAGED_COARSE, STAGED_FINE
    t0 = time.time()
    lay, _, e1, a1, T, ops_c = anneal(base, n_coarse, rng, coarse, mix_c, sa.T0, w, h, work)
    lay, _, e2, a2, T, ops_f = anneal(lay, STEPS - n_coarse, rng, FINE, mix_f, T, w, h, work)
    t_final = evaluate(lay, FINE, work / "final")
    return {
        "coarse_grid": coarse, "frac_coarse": FRAC, "mix": mix_name, "seed": seed,
        "secs": round(time.time() - t0, 1),
        "evals_coarse": e1, "evals_fine": e2,
        "accepted_coarse": a1, "accepted_fine": a2,
        "ops_coarse": ops_c, "ops_fine": ops_f,
        "final_64_C": round(t_final, 3),
    }


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    log = OUT / "results.jsonl"
    base, w, h = starting_layout()
    jobs = [(g, m, s) for g in COARSE_GRIDS for m in ("uniform", "staged") for s in SEEDS]
    for i, (g, m, s) in enumerate(jobs, 1):
        row = run_one(g, m, s, base, w, h, OUT / "work")
        with log.open("a") as fh:
            fh.write(json.dumps(row) + "\n")
        print(f"[{i}/{len(jobs)}] grid {g} {m} seed {s}: {row['final_64_C']:.3f} C "
              f"{row['secs']}s evals {row['evals_coarse']}+{row['evals_fine']}", flush=True)


if __name__ == "__main__":
    main()
