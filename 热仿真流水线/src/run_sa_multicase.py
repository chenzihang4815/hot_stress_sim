#!/usr/bin/env python3
"""Multi-case check of the chosen search, one HotSpot solve at a time.

Two outputs, both resumable from results.jsonl:

* traces: the best temperature at every step, for the chosen setting
  (60% big jumps, first 98.4% of steps priced on the 32 grid), over every
  layout. This is what trace.png was, extended past the single 7-chip case.

* mixes: the final 64-grid temperature of the move mixes, for every layout.
  A mix that runs through all 120 steps, and the same mix used only for the
  last 25% after a fixed coarse phase. This is what mix_compare.png was.

Work dirs are deleted after each solve. A field is reusable only while that
chip has not moved, so the cache stays small and is dropped between runs.
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

import sa_thermal_feedback as sa

ROOT = Path("/home/czh/任务仿真产生功耗文件工具")
SRC = ROOT / "热仿真流水线" / "src"
BASE = ROOT / "热仿真流水线"
OUT = BASE / "runs" / "superposition" / "multicase"

spec = importlib.util.spec_from_file_location("f", SRC / "run_sa_field_120.py")
f = importlib.util.module_from_spec(spec)
spec.loader.exec_module(f)
f.OUT = OUT / "_work"

INTP = 0.072
STEPS = 120
SEED = 22
AMB = 318.15
COARSE_FRAC = 0.984

JUMP = (0.10, 0.60, 0.20, 0.10)
UNIFORM = (0.40, 0.25, 0.20, 0.15)
NUDGE = (0.70, 0.10, 0.10, 0.10)
SWAP = (0.20, 0.20, 0.50, 0.10)
MIXES = [("传统", UNIFORM), ("多跳", JUMP), ("多微调", NUDGE), ("多交换", SWAP)]

CASES = [
    ("C1_4芯", "runs/sa_paper_cases/C1_1T2_1GPU/field_old/layout.json"),
    ("big_6芯", "runs/big_balanced_2/lbfgs_init.json"),
    ("hbm_7芯", "runs/with_hbm_2/Lstar.json"),
    ("t2_9芯", "runs/full_chain/pkg_2t2_2gpu/Lstar.json"),
    ("xeon2_9芯", "runs/full_chain/pkg_2xeon_2gpu/Lstar.json"),
    ("xeon3_11芯", "runs/full_chain/pkg_3xeon_3gpu/Lstar.json"),
]


def load(path):
    chips = [sa.Chip(c["name"], c["kind"], c["x_m"], c["y_m"], c["w_m"], c["h_m"],
                     c["power_W"], c["area_mm2"]) for c in json.loads((BASE / path).read_text())]
    return chips if sa.legal(chips, INTP, INTP) else sa.spread_to_gap(chips, INTP, INTP)


def field_temp(chips, cache, rows):
    """Peak from cached per-chip fields. Deletes the work dir after every solve."""
    f.ROWS = rows
    n = rows * rows
    acc = None
    fresh = 0
    t0 = time.time()
    for i, c in enumerate(chips):
        key = (i, round(c.x, 6), round(c.y, 6), round(c.w, 6), round(c.h, 6))
        if key not in cache:
            g = f.solve(f.solo(chips, i), "w")
            cache[key] = g[:n].copy()
            fresh += 1
            shutil.rmtree(f.OUT / "w", ignore_errors=True)
        acc = cache[key] if acc is None else acc + cache[key]
    acc = acc - (len(chips) - 1) * AMB
    return float(acc.max() - 273.15), time.time() - t0, fresh


def search(chips, coarse_mix, fine_mix, coarse_frac, want_trace):
    """One 120-step search. The two mixes split the coarse and the fine phase."""
    rng = random.Random(SEED)
    cache32, cache64 = {}, {}
    n_coarse = round(STEPS * coarse_frac)
    cur = best = deepcopy(chips)
    tmax, _, _ = field_temp(cur, cache64, 64)
    tmax_cur = tmax_best = tmax
    best_at = [tmax] if want_trace else None
    T = sa.T0
    t_used = 0.0
    for step in range(1, STEPS + 1):
        coarse = step <= n_coarse
        mix = coarse_mix if coarse else fine_mix
        trial, _ = sa.propose(cur, INTP, INTP, rng, mix)
        tries = 1
        while not sa.legal(trial, INTP, INTP) and tries < 25:
            trial, _ = sa.propose(cur, INTP, INTP, rng, mix)
            tries += 1
        if not sa.legal(trial, INTP, INTP):
            T *= sa.ALPHA
            if want_trace:
                best_at.append(tmax_best)
            continue
        t0 = time.time()
        tmax_new, _, _ = field_temp(trial, cache32 if coarse else cache64, 32 if coarse else 64)
        t_used += time.time() - t0
        j_cur, j_new = (tmax_cur - 45) / 25, (tmax_new - 45) / 25
        if rng.random() < (1.0 if j_new <= j_cur else math.exp((j_cur - j_new) / max(T, 1e-9))):
            cur, tmax_cur = trial, tmax_new
            if tmax_new < tmax_best:
                best, tmax_best = deepcopy(trial), tmax_new
        T *= sa.ALPHA
        if want_trace:
            best_at.append(tmax_best)
    final, _, _ = field_temp(best, {}, 64)
    return best_at, round(final, 3), round(t_used, 1)


# Coarse fractions to compare against the 98.4% trace, same big-jump mix.
# 0.0, 0.5, 0.75 and 0.984 are already stored, so only the missing one remains.
TRACE_FRACS = [1.0]


def traces_main() -> None:
    """Best-temperature curve per step, for several coarse fractions."""
    OUT.mkdir(parents=True, exist_ok=True)
    log = OUT / "results.jsonl"
    done = set()
    if log.exists():
        for line in log.read_text().splitlines():
            r = json.loads(line)
            if r["kind"] == "trace":
                done.add((r["label"], r["coarse_frac"]))
    for label, path in CASES:
        chips = load(path)
        for frac in TRACE_FRACS:
            if (label, frac) in done:
                continue
            t0 = time.time()
            best_at, final, secs = search(chips, JUMP, JUMP, frac, True)
            rec = {"kind": "trace", "label": label, "chips": len(chips),
                   "coarse": "多跳", "fine": "多跳", "coarse_frac": frac,
                   "final_64_C": final, "secs": secs, "wall_s": round(time.time() - t0, 1),
                   "best_C": [round(v, 3) for v in best_at]}
            with log.open("a") as fh:
                fh.write(json.dumps(rec, ensure_ascii=False) + "\n")
            print(json.dumps({k: v for k, v in rec.items() if k != "best_C"},
                             ensure_ascii=False), flush=True)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    log = OUT / "results.jsonl"
    done = set()
    if log.exists():
        for line in log.read_text().splitlines():
            r = json.loads(line)
            done.add((r["kind"], r["label"], r["coarse"], r["fine"]))

    def emit(rec):
        with log.open("a") as fh:
            fh.write(json.dumps(rec, ensure_ascii=False) + "\n")
        print(json.dumps({k: v for k, v in rec.items() if k != "best_C"},
                         ensure_ascii=False), flush=True)

    for label, path in CASES:
        chips = load(path)
        key = ("trace", label, "多跳", "多跳")
        if key not in done:
            t0 = time.time()
            best_at, final, secs = search(chips, JUMP, JUMP, COARSE_FRAC, True)
            emit({"kind": "trace", "label": label, "chips": len(chips),
                  "coarse": "多跳", "fine": "多跳", "coarse_frac": COARSE_FRAC,
                  "final_64_C": final, "secs": secs, "wall_s": round(time.time() - t0, 1),
                  "best_C": [round(v, 3) for v in best_at]})
        for name, mix in MIXES:
            key = ("whole", label, name, name)
            if key not in done:
                t0 = time.time()
                _, final, secs = search(chips, mix, mix, COARSE_FRAC, False)
                emit({"kind": "whole", "label": label, "chips": len(chips),
                      "coarse": name, "fine": name, "coarse_frac": COARSE_FRAC,
                      "final_64_C": final, "secs": secs, "wall_s": round(time.time() - t0, 1)})
        for cname, cmix in (("传统", UNIFORM), ("多跳", JUMP)):
            for fname, fmix in MIXES:
                key = ("split", label, cname, fname)
                if key not in done:
                    t0 = time.time()
                    _, final, secs = search(chips, cmix, fmix, 0.75, False)
                    emit({"kind": "split", "label": label, "chips": len(chips),
                          "coarse": cname, "fine": fname, "coarse_frac": 0.75,
                          "final_64_C": final, "secs": secs,
                          "wall_s": round(time.time() - t0, 1)})


if __name__ == "__main__":
    import sys
    traces_main() if "--traces" in sys.argv else main()
