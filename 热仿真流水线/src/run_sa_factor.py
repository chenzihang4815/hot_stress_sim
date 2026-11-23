#!/usr/bin/env python3
"""Two questions about the coarse/fine search, kept separate.

1. Decision risk. The coarse grid mis-estimates a fixed layout's temperature by
   some error (measured in runs/sa_paper_cases/grid_error). A move only risks
   being accepted or rejected wrongly when its temperature change is smaller
   than that error. Run the coarse phase and record every move's temperature
   change, then count how many fall inside the error.

2. Result quality. Sweep the coarse fraction and the move mix of each phase.
   Report the final layout's peak temperature, re-scored on the 64 grid, and
   the time. The mix and the fraction affect which layout is found; the grid
   error above says whether the coarse grid can be trusted to choose it.
"""

from __future__ import annotations

import importlib.util
import json
import random
import sys
import time
from pathlib import Path

ROOT = Path("/home/czh/任务仿真产生功耗文件工具")
SRC = ROOT / "热仿真流水线" / "src"
sys.path.insert(0, str(SRC))

import sa_thermal_feedback as sa  # noqa: E402

spec = importlib.util.spec_from_file_location("v3", SRC / "run_sa_mix_v3.py")
v3 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(v3)

spec = importlib.util.spec_from_file_location("pc", SRC / "run_sa_paper_cases.py")
pc = importlib.util.module_from_spec(spec)
spec.loader.exec_module(pc)

from hotspot_io import run_hotspot, write_hotspot_inputs  # noqa: E402

HS = Path("/home/czh/任务仿真产生功耗文件工具/HotSpot-master/hotspot")
OUT = ROOT / "热仿真流水线" / "runs" / "sa_factor"
STEPS = 120
SEEDS = [22, 23, 24]

# C1 small, C4 mid (where the new search was slightly worse), C6 large.
CASE_INDEX = [0, 1, 2, 3, 4, 5]
FRACTIONS = [0.25, 0.5, 0.75, 1.0]
COARSE_MIX = {"uniform": v3.UNIFORM, "jump": v3.JUMP_COARSE}
FINE_MIX = {"uniform": v3.UNIFORM, "nudge": v3.STAGED_FINE}


def score(chips, rows, work: Path) -> float:
    return v3.evaluate(chips, rows, work)


def case_list():
    chosen = []
    for k in CASE_INDEX:
        name, n_cpu, n_gpu, path = pc.CASES[k]
        chosen.append((name, pc.load(path, n_cpu, n_gpu)))
    return chosen


def decision_risk(name, base, w, h, thresh: float) -> dict:
    """Coarse phase only. Record each evaluated move's temperature change."""
    rng = random.Random(22)
    n_coarse = round(STEPS * 0.75)
    cur = base
    t_cur = score(cur, 32, OUT / "risk" / name)
    deltas = []
    T = sa.T0
    for step in range(1, n_coarse + 1):
        trial, _ = v3.propose(cur, w, h, rng, v3.JUMP_COARSE)
        tries = 1
        while not sa.legal(trial, w, h) and tries < v3.REDRAW_CAP:
            trial, _ = v3.propose(cur, w, h, rng, v3.JUMP_COARSE)
            tries += 1
        if not sa.legal(trial, w, h):
            T *= sa.ALPHA
            continue
        t_new = score(trial, 32, OUT / "risk" / name)
        deltas.append(t_new - t_cur)
        j_cur, j_new = (t_cur - 45) / 25, (t_new - 45) / 25
        if rng.random() < (1.0 if j_new <= j_cur else __import__("math").exp((j_cur - j_new) / max(T, 1e-9))):
            cur, t_cur = trial, t_new
        T *= sa.ALPHA
    n = len(deltas)
    half = n // 2

    def frac(vals):
        return round(sum(1 for d in vals if abs(d) < thresh) / len(vals), 3) if vals else None

    return {
        "case": name, "moves": n, "thresh_K": round(thresh, 3),
        "risk_all": frac(deltas), "risk_first_half": frac(deltas[:half]),
        "risk_second_half": frac(deltas[half:]),
        "median_abs_dT": round(sorted(abs(d) for d in deltas)[n // 2], 3) if n else None,
    }


def run_search(base, w, h, seed, frac, cname, fname, work: Path):
    rng = random.Random(seed)
    t0 = time.time()
    if frac <= 0:
        lay, *_ = v3.anneal(base, STEPS, rng, 64, v3.UNIFORM, sa.T0, w, h, work, False)
    else:
        n_coarse = round(STEPS * frac)
        lay, _, _, _, T = v3.anneal(base, n_coarse, rng, 32, COARSE_MIX[cname], sa.T0, w, h, work, True)
        lay, *_ = v3.anneal(lay, STEPS - n_coarse, rng, 64, FINE_MIX[fname], T, w, h, work, False)
    final = score(lay, 64, work / "final")
    return round(final, 3), round(time.time() - t0, 1)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    err = {r["case"]: abs(r["err32_K"]) for r in json.loads(
        (ROOT / "热仿真流水线" / "runs" / "sa_paper_cases" / "grid_error" / "grid_error.json").read_text())}
    cases = case_list()

    print("== 决策风险：温度变化小于网格误差的步占比", flush=True)
    risks = []
    for name, (base, w, h) in cases:
        r = decision_risk(name, base, w, h, err[name])
        risks.append(r)
        print(f"  {name}: 误差阈值 {r['thresh_K']} K, {r['moves']} 步, "
              f"可能被翻转 {r['risk_all']*100:.0f}% (前半 {r['risk_first_half']*100:.0f}% "
              f"后半 {r['risk_second_half']*100:.0f}%), 温度变化中位 {r['median_abs_dT']} K", flush=True)
    (OUT / "decision_risk.json").write_text(json.dumps(risks, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    print("== 粗占比二分", flush=True)
    # For each mix pair, bisect the coarse fraction. The two endpoints and the
    # midpoint come first; each round halves the side whose peak is lower,
    # until the two fractions still in play are less than 2% apart.
    log = OUT / "results.jsonl"
    done_frac = {}

    def cached(name, seed, frac, cname, fname):
        key = (name, seed, round(frac, 4), cname, fname)
        if key in done_frac:
            return done_frac[key]
        work = OUT / "work"
        final, secs = run_search(base_of[name], w_of[name], h_of[name], seed, frac, cname, fname, work)
        done_frac[key] = (final, secs)
        row = {"case": name, "frac": round(frac, 4), "coarse_mix": cname, "fine_mix": fname,
               "seed": seed, "final_64_C": final, "secs": secs}
        with log.open("a") as fh:
            fh.write(json.dumps(row) + "\n")
        return final, secs

    base_of, w_of, h_of = {}, {}, {}
    for name, (base, w, h) in cases:
        base_of[name], w_of[name], h_of[name] = base, w, h

    mixes = [(c, fn) for c in COARSE_MIX for fn in FINE_MIX]
    for mi, (cname, fname) in enumerate(mixes, 1):
        print(f"-- 比例组合 {mi}/{len(mixes)}: 粗 {cname} 细 {fname}", flush=True)
        for name, _ in cases:
            lo, hi = 0.0, 1.0
            best_frac, best_temp = None, None
            while hi - lo > 0.02:
                probes = sorted({round(lo, 4), round((lo + hi) / 2, 4), round(hi, 4)})
                temps = {}
                for frac in probes:
                    vals = [cached(name, seed, frac, cname, fname)[0] for seed in SEEDS]
                    temps[frac] = sum(vals) / len(vals)
                    print(f"   {name} 占比 {frac:.3f} 平均 {temps[frac]:.3f} C", flush=True)
                winner = min(probes, key=lambda f: temps[f])
                if best_temp is None or temps[winner] < best_temp:
                    best_frac, best_temp = winner, temps[winner]
                # Keep the winner inside the next bracket.
                if winner == probes[0]:
                    hi = probes[1]
                elif winner == probes[-1]:
                    lo = probes[-2]
                else:
                    lo, hi = probes[0], probes[-1]
                    # Narrow toward the better neighbour so the bracket still halves.
                    left, right = (probes[0] + winner) / 2, (winner + probes[-1]) / 2
                    lo, hi = (winner, probes[-1]) if temps[probes[-1]] < temps[probes[0]] else (probes[0], winner)
            print(f"   {name} 收敛: 占比 {best_frac:.3f} 峰值 {best_temp:.3f} C", flush=True)


if __name__ == "__main__":
    main()
