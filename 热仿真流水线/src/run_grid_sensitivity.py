#!/usr/bin/env python3
"""When does the 32 grid stop matching the 64 grid.

One hot die of a chosen size and power, centred in the fixed 72 mm frame, with
three cold neighbours so the mesh matches a real package. The peak is read on
the 64 grid and the 32 grid; the difference is the error a coarse search would
rank its moves by. One HotSpot at a time, work dir removed after each solve.
"""

from __future__ import annotations

import json
import shutil
import sys
from pathlib import Path

import numpy as np

ROOT = Path("/home/czh/任务仿真产生功耗文件工具")
SRC = ROOT / "热仿真流水线" / "src"
sys.path.insert(0, str(SRC))

import sa_thermal_feedback as sa  # noqa: E402
from hotspot_io import run_hotspot, write_hotspot_inputs  # noqa: E402

HS = ROOT / "HotSpot-master" / "hotspot"
OUT = ROOT / "热仿真流水线" / "runs" / "superposition" / "grid_sensitivity"
WORK = OUT / "_work"

# Side length in mm. 32-grid cells are 2.25 mm, so this runs from a die smaller
# than one cell up to the CPUs already shown to be insensitive.
SIDES = [1, 2, 3, 4, 6, 8, 11, 16, 22]
# Power density in W/mm^2. The current chips sit at 0.14 to 0.25.
DENS = [0.25, 1, 2, 4]


def chips_for(side_mm: float, dens: float):
    s = side_mm / 1e3
    hot = sa.Chip("HOT", "GPU", 0.036 - s / 2, 0.036 - s / 2, s, s,
                  dens * side_mm * side_mm, side_mm * side_mm)
    cold = []
    for i, (x, y) in enumerate(((0.004, 0.004), (0.050, 0.004), (0.004, 0.050))):
        cold.append(sa.Chip(f"COLD{i}", "CPU", x, y, 0.012, 0.012, 0.0, 144.0))
    return [hot] + cold


def peak(chips, rows: int) -> float:
    if WORK.exists():
        shutil.rmtree(WORK)
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
    shutil.rmtree(WORK, ignore_errors=True)
    return max(vals) - 273.15


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    log = OUT / "results.jsonl"
    done = {}
    if log.exists():
        for line in log.read_text().splitlines():
            r = json.loads(line)
            done[(r["side_mm"], r["dens"])] = r
    for side in SIDES:
        for dens in DENS:
            rec = done.get((side, dens))
            if rec is None or "t16" not in rec:
                cs = chips_for(side, dens)
                t64 = rec["t64"] if rec else round(peak(cs, 64), 3)
                t32 = rec["t32"] if rec else round(peak(cs, 32), 3)
                t16 = round(peak(cs, 16), 3)
                rec = {"side_mm": side, "dens": dens, "power_W": round(dens * side * side, 2),
                       "cells32": round(side / 2.25, 2), "cells16": round(side / 4.5, 2),
                       "t64": t64, "t32": t32, "t16": t16,
                       "err": round(t32 - t64, 3), "err16": round(t16 - t64, 3)}
                done[(side, dens)] = rec
                with log.open("a") as fh:
                    fh.write(json.dumps(rec) + "\n")
            print(json.dumps(rec), flush=True)


if __name__ == "__main__":
    main()
