#!/usr/bin/env python3
"""Grid error on the fixed 72 mm frame. One layout, three grids, one at a time."""

from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np

import sa_thermal_feedback as sa
from hotspot_io import run_hotspot, write_hotspot_inputs

ROOT = Path("/home/czh/任务仿真产生功耗文件工具")
HS = ROOT / "HotSpot-master" / "hotspot"
BASE = ROOT / "热仿真流水线"
OUT = BASE / "runs" / "superposition" / "grid_error"

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
    return chips if sa.legal(chips, 0.072, 0.072) else sa.spread_to_gap(chips, 0.072, 0.072)


def peak(chips, rows) -> float:
    work = OUT / "work"
    if work.exists():
        import shutil
        shutil.rmtree(work)
    write_hotspot_inputs(sa.to_rects(chips), work, rows=rows, cols=rows)
    return float(run_hotspot(HS, work)["Tmax_C"])


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    rows = []
    for label, path in CASES:
        chips = load(path)
        rec = {"label": label, "chips": len(chips)}
        base = None
        for g in (64, 32, 16):
            t0 = time.time()
            t = peak(chips, g)
            rec[f"g{g}"] = round(t, 3)
            rec[f"s{g}"] = round(time.time() - t0, 2)
            if g == 64:
                base = t
            else:
                rec[f"err{g}"] = round(t - base, 3)
        rows.append(rec)
        print(json.dumps(rec, ensure_ascii=False), flush=True)
    (OUT / "grid_error.json").write_text(json.dumps(rows, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
