#!/usr/bin/env python3
"""默认路径：三颗 CPU/GPU/CIM → L-BFGS 初始布局 → HotSpot。"""

from __future__ import annotations

import json
import shutil
import sys
from pathlib import Path

ROOT = Path("/home/czh/任务仿真产生功耗文件工具")
PIPE = ROOT / "热仿真流水线"
SRC = PIPE / "src"
sys.path.insert(0, str(SRC))

from combos import all_mixes  # noqa: E402
from defaults import INITIAL_LAYOUT, INITIAL_LAYOUT_SEED  # noqa: E402
from hotspot_io import run_hotspot, write_hotspot_inputs  # noqa: E402
from lbfgs import random_lbfgs_layout  # noqa: E402
from visualize import draw_floorplan  # noqa: E402

HOTSPOT_BIN = ROOT / "HotSpot-master" / "hotspot"
RUNS = PIPE / "runs" / "lbfgs_default"
BACKUP = ROOT / "成功运行备份" / "HotSpot_初始布局" / "lbfgs_default"

LAYOUT_FNS = {"lbfgs": random_lbfgs_layout}


def _copy_run(src: Path, dst: Path) -> None:
    if dst.exists():
        shutil.rmtree(dst)
    shutil.copytree(src, dst)


def run_one(name: str, blocks, seed: int) -> dict:
    if INITIAL_LAYOUT not in LAYOUT_FNS:
        raise SystemExit(f"默认布局 {INITIAL_LAYOUT} 未接入 run_pipeline")
    out = RUNS / name
    if out.exists():
        shutil.rmtree(out)
    out.mkdir(parents=True)
    rects, meta = LAYOUT_FNS[INITIAL_LAYOUT](blocks, seed=seed)
    write_hotspot_inputs(rects, out)
    thermal = run_hotspot(HOTSPOT_BIN, out)
    draw_floorplan(rects, out / "floorplan.png", f"{name}  {INITIAL_LAYOUT}")
    summary = {
        "name": name,
        "algorithm": INITIAL_LAYOUT,
        "n_blocks": len(blocks),
        "layout": meta,
        "thermal": {
            "Tmax_C": thermal.get("Tmax_C"),
            "Tmin_C": thermal.get("Tmin_C"),
            "temperatures_K": thermal.get("temperatures_K"),
            "returncode": thermal.get("returncode"),
        },
        "power_total_W": sum(b.power_W for b in blocks),
        "area_total_mm2": sum(b.area_mm2 for b in blocks),
    }
    (out / "summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    _copy_run(out, BACKUP / name)
    return summary


def main() -> None:
    if not HOTSPOT_BIN.exists():
        raise SystemExit(f"找不到 HotSpot: {HOTSPOT_BIN}")
    if INITIAL_LAYOUT != "lbfgs":
        raise SystemExit(f"已拍板初始布局为 lbfgs，当前 defaults.py 是 {INITIAL_LAYOUT}")
    RUNS.mkdir(parents=True, exist_ok=True)
    BACKUP.mkdir(parents=True, exist_ok=True)

    mixes = all_mixes()
    reports = []
    for mix in mixes:
        kinds = [b.kind for b in mix.blocks]
        if kinds != ["CPU", "GPU", "CIM"]:
            raise SystemExit(f"{mix.id} 不是三颗 CPU/GPU/CIM: {kinds}")
        print(f"==> {mix.id}  {INITIAL_LAYOUT}", flush=True)
        reports.append(run_one(mix.id, mix.blocks, INITIAL_LAYOUT_SEED))

    table = [
        {
            "case": r["name"],
            "algorithm": r["algorithm"],
            "outline_mm2": r["layout"]["outline_mm2"],
            "deadspace_ratio": r["layout"]["deadspace_ratio"],
            "overlap_m2": r["layout"]["overlap_m2"],
            "Tmax_C": r["thermal"]["Tmax_C"],
            "Tmin_C": r["thermal"]["Tmin_C"],
            "power_total_W": r["power_total_W"],
        }
        for r in reports
    ]
    (PIPE / "lbfgs_default_summary.json").write_text(
        json.dumps({"algorithm": INITIAL_LAYOUT, "cases": table, "details": reports}, indent=2, ensure_ascii=False)
        + "\n",
        encoding="utf-8",
    )
    shutil.copy2(PIPE / "lbfgs_default_summary.json", BACKUP / "lbfgs_default_summary.json")
    print(json.dumps(table, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
