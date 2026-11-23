#!/usr/bin/env python3
"""Layouts for the multi-scheme heat-vs-coupled SA study.

Schemes:
  pkg_5cpu_3gpu_5cim  — existing 13-chip large package (reuse layout.json)
  pkg_5cpu_4gpu_5cim  — existing 14-chip / 4-GPU package (reuse layout.json)
  pkg_4cpu_2gpu_4cim  — new mid-size mix (L-BFGS, seed 13)
"""

from __future__ import annotations

import json
import shutil
import sys
from dataclasses import replace
from pathlib import Path

ROOT = Path("/home/czh/任务仿真产生功耗文件工具")
SRC = ROOT / "热仿真流水线" / "src"
sys.path.insert(0, str(SRC))

from blocks import BACKUP, accelforge_from_summary, cim_basic_analog, cpu_from_gem5_mcpat  # noqa: E402
from defaults import INITIAL_LAYOUT_SEED  # noqa: E402
from lbfgs import random_lbfgs_layout  # noqa: E402
from run_large_package import gpu  # noqa: E402

OUT = ROOT / "热仿真流水线" / "runs" / "sa_compare_schemes"


def cpu(txt: Path, name: str, notes: str):
    b = cpu_from_gem5_mcpat(txt, name=name, notes=notes)
    return replace(b, name=name, kind="CPU")


def cim(combo: str, name: str):
    b = accelforge_from_summary(combo, kind="CIM", name=name)
    return replace(b, name=name, kind="CIM")


def rects_to_layout(rects) -> list[dict]:
    return [
        {
            "name": r.name,
            "kind": r.kind,
            "x_m": r.x,
            "y_m": r.y,
            "w_m": r.w,
            "h_m": r.h,
            "area_mm2": r.area_mm2,
            "power_W": r.power_W,
        }
        for r in rects
    ]


def main() -> None:
    if OUT.exists():
        shutil.rmtree(OUT)
    OUT.mkdir(parents=True)

    gdir = BACKUP / "Accel-Sim/小面积GPU_backprop"
    large = ROOT / "热仿真流水线" / "runs" / "large_package"

    copied = []
    for sid, title, why in (
        (
            "pkg_5cpu_3gpu_5cim",
            "5 CPU + 3 GPU + 5 CIM（13 颗）",
            "已有大封装；对照实验的基准规模",
        ),
        (
            "pkg_5cpu_4gpu_5cim",
            "5 CPU + 4 GPU + 5 CIM（14 颗）",
            "多一颗同规格 3070，看同类 GPU 间距是否被应力拉开",
        ),
    ):
        src = large / sid / "layout.json"
        if not src.exists():
            raise SystemExit(f"missing {src}")
        dest = OUT / sid
        dest.mkdir()
        shutil.copy2(src, dest / "layout.json")
        meta = {
            "id": sid,
            "title": title,
            "why": why,
            "layout": str(dest / "layout.json"),
            "source": str(src),
            "n_chiplets": len(json.loads(src.read_text())),
            "atsim_pitch_die_mm": 2.5 if sid.endswith("4gpu_5cim") else 2.0,
            "atsim_pitch_coarse_mm": 10.0 if sid.endswith("4gpu_5cim") else 8.0,
            "atsim_nz": {"nz_die": 2, "nz_spreader": 2, "nz_sink": 3, "nz_tim": 1, "nz_substrate": 2, "nz_interposer": 2}
            if sid.endswith("4gpu_5cim")
            else None,
        }
        (dest / "scheme.json").write_text(json.dumps(meta, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
        copied.append(meta)
        print("copied", sid, "n=", meta["n_chiplets"])

    cpu_matmul = cpu(BACKUP / "gem5/x86_tasks/matmul_timing_l2/mcpat.txt", "CPU_matmul", "x86 Timing+L1/L2 矩阵乘")
    cpu_threads = cpu(BACKUP / "gem5/x86_tasks/threads_timing_l2/mcpat.txt", "CPU_threads", "x86 Timing+L1/L2 threads")
    cpu_arm = cpu(BACKUP / "gem5/output/mcpat_from_gem5.txt", "CPU_armL2", "ARM Timing+L1/L2 hello")
    cpu_starter = cpu(BACKUP / "gem5/多案例/03_starter_se_atomic/mcpat.txt", "CPU_starter", "ARM starter atomic hello")
    gpu_3070 = gpu(gdir / "SM86_RTX3070/accelwattch_power_report.log", "GPU_3070", 392.0, "RTX3070 backprop")
    gpu_titanx = gpu(gdir / "SM6_TITANX/accelwattch_power_report.log", "GPU_TITANX", 471.0, "TitanX backprop")
    blocks = [
        cpu_matmul,
        cpu_threads,
        cpu_arm,
        cpu_starter,
        gpu_3070,
        gpu_titanx,
        cim("jia_jssc_2020+builtin", "CIM_jia"),
        cim("colonnade_jssc_2021+builtin", "CIM_colonnade"),
        cim("wang_vlsi_2022+builtin", "CIM_wang"),
        replace(cim_basic_analog(), name="CIM_basic", kind="CIM"),
    ]
    sid = "pkg_4cpu_2gpu_4cim"
    print("==", sid, "n=", len(blocks), "lbfgs seed", INITIAL_LAYOUT_SEED, flush=True)
    rects, pack = random_lbfgs_layout(blocks, seed=INITIAL_LAYOUT_SEED)
    dest = OUT / sid
    dest.mkdir()
    layout = rects_to_layout(rects)
    (dest / "layout.json").write_text(json.dumps(layout, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    meta = {
        "id": sid,
        "title": "4 CPU + 2 GPU + 4 CIM（10 颗，新方案）",
        "why": "中等规模：两颗不同 GPU、四颗 CPU、四颗 CIM，比 13/14 芯省内存，仍有异构旋钮",
        "layout": str(dest / "layout.json"),
        "source": "lbfgs",
        "n_chiplets": len(blocks),
        "kinds": {"CPU": 4, "GPU": 2, "CIM": 4},
        "power_total_W": sum(b.power_W for b in blocks),
        "area_total_mm2": sum(b.area_mm2 for b in blocks),
        "outline_mm2": pack["outline_mm2"],
        "atsim_pitch_die_mm": 2.0,
        "atsim_pitch_coarse_mm": 8.0,
        "atsim_nz": None,
        "chips": [
            {"name": b.name, "kind": b.kind, "area_mm2": b.area_mm2, "power_W": b.power_W, "notes": b.notes}
            for b in blocks
        ],
    }
    (dest / "scheme.json").write_text(json.dumps(meta, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    copied.append(meta)
    (OUT / "index.json").write_text(json.dumps(copied, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print("wrote", OUT)


if __name__ == "__main__":
    main()
