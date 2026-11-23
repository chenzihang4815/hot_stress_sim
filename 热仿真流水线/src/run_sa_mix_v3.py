#!/usr/bin/env python3
"""Staged move mixes on top of redraw.

Every coarse proposal is re-drawn until legal (cap 25), so a mix that jumps
more no longer buys that with wasted steps. Coarse fraction 0.75, coarse grid
32, fine grid 64, final layout re-scored on 64. Three variants, three seeds:

- redraw_uniform: both stages use 40/25/20/15. The v2 winner, the baseline.
- redraw_staged: coarse stage jumps more (20/45/25/10), fine stage nudges more
  (80/5/5/10).
- redraw_jump: coarse stage jumps harder still (10/60/20/10), same fine stage.
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
OUT = ROOT / "热仿真流水线" / "runs" / "sa_mix_v3"
FINE = 64
COARSE = 32
STEPS = 120
FRAC = 0.75
SEEDS = [22, 23, 24]
REDRAW_CAP = 25

UNIFORM = (0.40, 0.25, 0.20, 0.15)
STAGED_COARSE = (0.20, 0.45, 0.25, 0.10)
STAGED_FINE = (0.80, 0.05, 0.05, 0.10)
JUMP_COARSE = (0.10, 0.60, 0.20, 0.10)


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


def anneal(start, steps, rng, rows, mix, T, w, h, work: Path, redraw: bool):
    cur = best = deepcopy(start)
    t_cur = t_best = evaluate(cur, rows, work)
    evals = accepted = 1
    drawn = {"nudge": 0, "jump": 0, "swap": 0, "rotate": 0}
    for _ in range(steps):
        trial, op = propose(cur, w, h, rng, mix)
        tries = 1
        while redraw and not sa.legal(trial, w, h) and tries < REDRAW_CAP:
            trial, op = propose(cur, w, h, rng, mix)
            tries += 1
        if not sa.legal(trial, w, h):
            T *= sa.ALPHA
            continue
        drawn[op] += 1
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
    return best, evals, accepted, drawn, T


def run_one(variant, seed, base, w, h, work: Path) -> dict:
    rng = random.Random(seed)
    n_coarse = round(STEPS * FRAC)
    if variant == "redraw_uniform":
        coarse, fine = UNIFORM, UNIFORM
    elif variant == "redraw_staged":
        coarse, fine = STAGED_COARSE, STAGED_FINE
    elif variant == "redraw_jump":
        coarse, fine = JUMP_COARSE, STAGED_FINE
    else:
        raise ValueError(variant)
    t0 = time.time()
    lay, e1, a1, d1, T = anneal(base, n_coarse, rng, COARSE, coarse, sa.T0, w, h, work, True)
    lay, e2, a2, d2, T = anneal(lay, STEPS - n_coarse, rng, FINE, fine, T, w, h, work, False)
    t_final = evaluate(lay, FINE, work / "final")
    return {
        "variant": variant, "seed": seed, "secs": round(time.time() - t0, 1),
        "evals_coarse": e1, "accepted_coarse": a1, "drawn_coarse": d1,
        "evals_fine": e2, "accepted_fine": a2, "drawn_fine": d2,
        "final_64_C": round(t_final, 3),
    }


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    log = OUT / "results.jsonl"
    base, w, h = starting_layout()
    variants = ["redraw_uniform", "redraw_staged", "redraw_jump"]
    jobs = [(v, s) for v in variants for s in SEEDS]
    for i, (v, s) in enumerate(jobs, 1):
        row = run_one(v, s, base, w, h, OUT / "work")
        with log.open("a") as fh:
            fh.write(json.dumps(row) + "\n")
        print(f"[{i}/{len(jobs)}] {v} seed {s}: {row['final_64_C']:.3f} C {row['secs']}s "
              f"粗{row['drawn_coarse']} 细{row['drawn_fine']}", flush=True)


if __name__ == "__main__":
    main()
