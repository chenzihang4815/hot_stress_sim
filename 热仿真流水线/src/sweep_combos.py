#!/usr/bin/env python3
"""对照实验：三颗 CPU/GPU/CIM × bstree / mcts / lbfgs。

日常初始布局已拍板为 lbfgs，见 defaults.py 和 run_pipeline.py。
本脚本只用来复现 2026-10-19 的算法对比，不是默认入口。
"""

from __future__ import annotations

import csv
import json
import shutil
import sys
from pathlib import Path

ROOT = Path("/home/czh/任务仿真产生功耗文件工具")
PIPE = ROOT / "热仿真流水线"
SRC = PIPE / "src"
sys.path.insert(0, str(SRC))

from bstree import random_bstree_layout  # noqa: E402
from combos import all_mixes  # noqa: E402
from hotspot_io import run_hotspot, write_hotspot_inputs  # noqa: E402
from lbfgs import random_lbfgs_layout  # noqa: E402
from mcts import random_mcts_layout  # noqa: E402
from visualize import draw_compare_bars, draw_floorplan, draw_tmax_grouped  # noqa: E402

HOTSPOT_BIN = ROOT / "HotSpot-master" / "hotspot"
RUNS = PIPE / "runs" / "combo_sweep"
BACKUP = ROOT / "成功运行备份" / "HotSpot_初始布局" / "combo_sweep"

ALGS = [
    ("bstree", random_bstree_layout, 7),
    ("mcts", random_mcts_layout, 11),
    ("lbfgs", random_lbfgs_layout, 13),
]


def _copy(src: Path, dst: Path) -> None:
    if dst.exists():
        shutil.rmtree(dst)
    shutil.copytree(src, dst)


def run_one(mix_id: str, alg: str, blocks, fn, seed: int) -> dict:
    name = f"{mix_id}__{alg}"
    out = RUNS / name
    if out.exists():
        shutil.rmtree(out)
    out.mkdir(parents=True)
    rects, meta = fn(blocks, seed=seed)
    write_hotspot_inputs(rects, out)
    thermal = run_hotspot(HOTSPOT_BIN, out)
    draw_floorplan(rects, out / "floorplan.png", f"{mix_id} / {alg}")
    temps = thermal.get("temperatures_K") or {}
    block_c = {k: v - 273.15 for k, v in temps.items() if not str(k).startswith("hsp_")}
    summary = {
        "mix": mix_id,
        "algorithm": alg,
        "n_blocks": len(blocks),
        "layout": meta,
        "thermal": {
            "Tmax_C": thermal.get("Tmax_C"),
            "Tmin_C": thermal.get("Tmin_C"),
            "block_T_C": block_c,
            "returncode": thermal.get("returncode"),
        },
        "power_total_W": sum(b.power_W for b in blocks),
        "area_total_mm2": sum(b.area_mm2 for b in blocks),
        "blocks": [
            {
                "name": b.name,
                "kind": b.kind,
                "area_mm2": b.area_mm2,
                "power_W": b.power_W,
                "notes": b.notes,
            }
            for b in blocks
        ],
    }
    (out / "summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    _copy(out, BACKUP / name)
    return summary


def main() -> None:
    if not HOTSPOT_BIN.exists():
        raise SystemExit(f"找不到 HotSpot: {HOTSPOT_BIN}")
    if RUNS.exists():
        shutil.rmtree(RUNS)
    RUNS.mkdir(parents=True)
    if BACKUP.exists():
        shutil.rmtree(BACKUP)
    BACKUP.mkdir(parents=True)

    mixes = all_mixes()
    reports = []
    for mix in mixes:
        kinds = [b.kind for b in mix.blocks]
        names = [b.name for b in mix.blocks]
        if len(mix.blocks) != 3 or kinds != ["CPU", "GPU", "CIM"]:
            raise SystemExit(f"{mix.id} 不是三颗 CPU/GPU/CIM: {names} {kinds}")
        print(f"== mix {mix.id}: {mix.title}", flush=True)
        for alg, fn, seed in ALGS:
            print(f"   {alg}", flush=True)
            reports.append(run_one(mix.id, alg, mix.blocks, fn, seed))

    rows = []
    for r in reports:
        rows.append(
            {
                "mix": r["mix"],
                "algorithm": r["algorithm"],
                "n_blocks": r["n_blocks"],
                "area_total_mm2": r["area_total_mm2"],
                "power_total_W": r["power_total_W"],
                "outline_mm2": r["layout"]["outline_mm2"],
                "deadspace_ratio": r["layout"]["deadspace_ratio"],
                "overlap_m2": r["layout"]["overlap_m2"],
                "Tmax_C": r["thermal"]["Tmax_C"],
                "Tmin_C": r["thermal"]["Tmin_C"],
                "dT_C": (r["thermal"]["Tmax_C"] or 0) - (r["thermal"]["Tmin_C"] or 0),
            }
        )

    (PIPE / "combo_sweep_summary.json").write_text(
        json.dumps(
            {
                "mixes": [
                    {
                        "id": m.id,
                        "title": m.title,
                        "why": m.why,
                        "blocks": [
                            {
                                "name": b.name,
                                "kind": b.kind,
                                "area_mm2": b.area_mm2,
                                "power_W": b.power_W,
                            }
                            for b in m.blocks
                        ],
                    }
                    for m in mixes
                ],
                "table": rows,
                "details": reports,
            },
            indent=2,
            ensure_ascii=False,
        )
        + "\n",
        encoding="utf-8",
    )
    with (PIPE / "combo_sweep_table.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    titles = {m.id: m.id for m in mixes}
    fig_dir = PIPE / "compare"
    if fig_dir.exists():
        shutil.rmtree(fig_dir)
    fig_dir.mkdir()
    draw_tmax_grouped(rows, titles, fig_dir / "tmax_by_mix.png")
    draw_compare_bars(rows, titles, "outline_mm2", "outline area (mm2)", fig_dir / "outline_by_mix.png")
    draw_compare_bars(rows, titles, "deadspace_ratio", "deadspace ratio", fig_dir / "deadspace_by_mix.png")
    draw_compare_bars(rows, titles, "dT_C", "delta T (C)", fig_dir / "dT_by_mix.png")

    shutil.copy2(PIPE / "combo_sweep_summary.json", BACKUP / "combo_sweep_summary.json")
    shutil.copy2(PIPE / "combo_sweep_table.csv", BACKUP / "combo_sweep_table.csv")
    shutil.copytree(fig_dir, BACKUP / "compare")
    print(json.dumps(rows, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
