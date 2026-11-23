#!/usr/bin/env python3
"""Bisect the coarse fraction upward until the final temperature rises.

For each layout, keep the largest fraction of steps priced on the 32 grid that
still finishes within 0.05 C of an all-fine search, measured by re-scoring the
best layout on the 64 grid. The bracket stops once it is under 2% wide.
"""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import sa_thermal_feedback as sa

ROOT = Path("/home/czh/任务仿真产生功耗文件工具")
SRC = ROOT / "热仿真流水线" / "src"
BASE = ROOT / "热仿真流水线"
OUT = BASE / "runs" / "superposition" / "bisect"
spec = importlib.util.spec_from_file_location("c", SRC / "run_sa_coarse.py")
c = importlib.util.module_from_spec(spec)
spec.loader.exec_module(c)

TOL = 0.05  # C. A coarse run is acceptable while it stays within this of all-fine.


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    log = OUT / "results.jsonl"
    done = {}
    if log.exists():
        for line in log.read_text().splitlines():
            r = json.loads(line)
            done[(r["label"], r["coarse_frac"])] = r

    def one(label, chips, frac):
        key = (label, round(frac, 4))
        if key in done:
            return done[key]
        best, final, secs = c.run(chips, frac)
        rec = {"label": label, "chips": len(chips), "coarse_frac": round(frac, 4),
               "final_64_C": final, "secs": secs}
        done[key] = rec
        with log.open("a") as fh:
            fh.write(json.dumps(rec) + "\n")
        print(json.dumps(rec, ensure_ascii=False), flush=True)
        return rec

    for label, path in c.CASES:
        chips = c.load(path)
        base = one(label, chips, 0.0)["final_64_C"]
        lo, hi = 0.75, 1.0
        best_frac = 0.75
        while hi - lo > 0.02:
            mid = round((lo + hi) / 2, 4)
            rec = one(label, chips, mid)
            if rec["final_64_C"] <= base + TOL:
                best_frac, lo = mid, mid
            else:
                hi = mid
        print(f"{label} 收敛: 粗占比 {best_frac:.3f} 仍不升温", flush=True)


if __name__ == "__main__":
    main()
