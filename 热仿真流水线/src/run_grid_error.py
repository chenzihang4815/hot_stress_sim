#!/usr/bin/env python3
"""Grid error on a fixed layout.

Take the layout the traditional SA found for each of the six packages and
re-score that exact layout on a 32 grid and on a 64 grid. The difference is the
temperature error of the coarse grid, with the search taken out of it.
"""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

ROOT = Path("/home/czh/任务仿真产生功耗文件工具")
SRC = ROOT / "热仿真流水线" / "src"
sys.path.insert(0, str(SRC))

import sa_thermal_feedback as sa  # noqa: E402
from hotspot_io import run_hotspot, write_hotspot_inputs  # noqa: E402

HS = Path("/home/czh/任务仿真产生功耗文件工具/HotSpot-master/hotspot")
CASES = ROOT / "热仿真流水线" / "runs" / "sa_paper_cases"
OUT = CASES / "grid_error"


def score(chips, rows, work: Path) -> float:
    if work.exists():
        import shutil
        shutil.rmtree(work)
    work.mkdir(parents=True)
    write_hotspot_inputs(sa.to_rects(chips), work, rows=rows, cols=rows)
    return float(run_hotspot(HS, work)["Tmax_C"])


def load_layout(path: Path):
    return [sa.Chip(c["name"], c["kind"], c["x_m"], c["y_m"], c["w_m"], c["h_m"],
                    c["power_W"], c["area_mm2"]) for c in json.loads(path.read_text())]


def main() -> None:
    spec = importlib.util.spec_from_file_location("p", SRC / "run_sa_paper_cases.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    OUT.mkdir(parents=True, exist_ok=True)
    rows = []
    for name, n_cpu, n_gpu, path in m.CASES:
        case = CASES / name
        lay = case / "field_old" / "layout.json"
        if not lay.exists():
            print(f"{name}: 没有传统 SA 的布局，跳过", flush=True)
            continue
        chips = load_layout(lay)
        t64 = score(chips, 64, OUT / name / "g64")
        t32 = score(chips, 32, OUT / name / "g32")
        t16 = score(chips, 16, OUT / name / "g16")
        row = {"case": name, "chips": len(chips),
               "grid64_C": round(t64, 3), "grid32_C": round(t32, 3), "grid16_C": round(t16, 3),
               "err32_K": round(t32 - t64, 3), "err16_K": round(t16 - t64, 3)}
        rows.append(row)
        print(f"{name}: 64格 {t64:.3f}  32格 {t32:.3f} (误差 {t32-t64:+.3f})  "
              f"16格 {t16:.3f} (误差 {t16-t64:+.3f})", flush=True)
    (OUT / "grid_error.json").write_text(
        json.dumps(rows, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
