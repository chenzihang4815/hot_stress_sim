#!/usr/bin/env python3
"""Add 3/4/11/12/15/16-chip L-BFGS layouts into sa_compare_schemes (keep existing)."""

from __future__ import annotations

import json
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


def coarse_nz():
    return {
        "nz_die": 2,
        "nz_spreader": 2,
        "nz_sink": 3,
        "nz_tim": 1,
        "nz_substrate": 2,
        "nz_interposer": 2,
    }


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    gdir = BACKUP / "Accel-Sim/小面积GPU_backprop"
    cpu_matmul = cpu(BACKUP / "gem5/x86_tasks/matmul_timing_l2/mcpat.txt", "CPU_matmul", "x86 Timing+L1/L2 矩阵乘")
    cpu_threads = cpu(BACKUP / "gem5/x86_tasks/threads_timing_l2/mcpat.txt", "CPU_threads", "x86 Timing+L1/L2 threads")
    cpu_arm = cpu(BACKUP / "gem5/output/mcpat_from_gem5.txt", "CPU_armL2", "ARM Timing+L1/L2 hello")
    cpu_atomic = cpu(BACKUP / "gem5/多案例/02_se_atomic_无cache/mcpat.txt", "CPU_atomic", "ARM Atomic hello")
    cpu_starter = cpu(BACKUP / "gem5/多案例/03_starter_se_atomic/mcpat.txt", "CPU_starter", "ARM starter atomic hello")
    gpu_3070 = gpu(gdir / "SM86_RTX3070/accelwattch_power_report.log", "GPU_3070", 392.0, "RTX3070 backprop")
    gpu_2060 = gpu(gdir / "SM75_RTX2060_S/accelwattch_power_report.log", "GPU_2060S", 445.0, "RTX2060S backprop")
    gpu_titanx = gpu(gdir / "SM6_TITANX/accelwattch_power_report.log", "GPU_TITANX", 471.0, "TitanX backprop")
    gpu_3070b = replace(gpu_3070, name="GPU_3070b", notes=gpu_3070.notes + "；第二颗 3070")
    gpu_3070c = replace(gpu_3070, name="GPU_3070c", notes=gpu_3070.notes + "；第三颗 3070")
    cim_jia = cim("jia_jssc_2020+builtin", "CIM_jia")
    cim_col = cim("colonnade_jssc_2021+builtin", "CIM_colonnade")
    cim_wang = cim("wang_vlsi_2022+builtin", "CIM_wang")
    cim_basic = replace(cim_basic_analog(), name="CIM_basic", kind="CIM")
    cim_sin = cim("sinangil_jssc_2021+builtin", "CIM_sinangil")
    cim_jia_b = replace(cim_jia, name="CIM_jia_b", notes=cim_jia.notes + "；第二颗 jia")

    schemes = [
        {
            "id": "pkg_1cpu_1gpu_1cim",
            "title": "1 CPU + 1 GPU + 1 CIM（3 颗）",
            "why": "最小异构：matmul + 3070 + jia",
            "blocks": [cpu_matmul, gpu_3070, cim_jia],
            "pitch": (2.0, 8.0),
            "nz": None,
        },
        {
            "id": "pkg_2cpu_1gpu_1cim",
            "title": "2 CPU + 1 GPU + 1 CIM（4 颗）",
            "why": "两颗 x86 CPU + 一颗 GPU",
            "blocks": [cpu_matmul, cpu_threads, gpu_3070, cim_jia],
            "pitch": (2.0, 8.0),
            "nz": None,
        },
        {
            "id": "pkg_4cpu_3gpu_4cim",
            "title": "4 CPU + 3 GPU + 4 CIM（11 颗）",
            "why": "夹在 10 与 13 之间",
            "blocks": [
                cpu_matmul, cpu_threads, cpu_arm, cpu_starter,
                gpu_3070, gpu_2060, gpu_titanx,
                cim_jia, cim_col, cim_wang, cim_basic,
            ],
            "pitch": (2.0, 8.0),
            "nz": None,
        },
        {
            "id": "pkg_4cpu_3gpu_5cim",
            "title": "4 CPU + 3 GPU + 5 CIM（12 颗）",
            "why": "11 芯再加 sinangil",
            "blocks": [
                cpu_matmul, cpu_threads, cpu_arm, cpu_starter,
                gpu_3070, gpu_2060, gpu_titanx,
                cim_jia, cim_col, cim_wang, cim_basic, cim_sin,
            ],
            "pitch": (2.0, 8.0),
            "nz": None,
        },
        {
            "id": "pkg_5cpu_5gpu_5cim",
            "title": "5 CPU + 5 GPU + 5 CIM（15 颗）",
            "why": "14 芯再加第三颗 3070",
            "blocks": [
                cpu_matmul, cpu_threads, cpu_arm, cpu_atomic, cpu_starter,
                gpu_3070, gpu_2060, gpu_titanx, gpu_3070b, gpu_3070c,
                cim_jia, cim_col, cim_wang, cim_basic, cim_sin,
            ],
            "pitch": (2.5, 10.0),
            "nz": coarse_nz(),
        },
        {
            "id": "pkg_5cpu_5gpu_6cim",
            "title": "5 CPU + 5 GPU + 6 CIM（16 颗）",
            "why": "15 芯再加第二颗 jia",
            "blocks": [
                cpu_matmul, cpu_threads, cpu_arm, cpu_atomic, cpu_starter,
                gpu_3070, gpu_2060, gpu_titanx, gpu_3070b, gpu_3070c,
                cim_jia, cim_col, cim_wang, cim_basic, cim_sin, cim_jia_b,
            ],
            "pitch": (2.5, 10.0),
            "nz": coarse_nz(),
        },
    ]

    index_path = OUT / "index.json"
    index = json.loads(index_path.read_text(encoding="utf-8")) if index_path.exists() else []
    by_id = {r["id"]: r for r in index}

    for sch in schemes:
        sid = sch["id"]
        blocks = sch["blocks"]
        names = [b.name for b in blocks]
        if len(names) != len(set(names)):
            raise SystemExit(f"dup names {names}")
        print(f"== {sid} n={len(blocks)} lbfgs seed {INITIAL_LAYOUT_SEED}", flush=True)
        rects, pack = random_lbfgs_layout(blocks, seed=INITIAL_LAYOUT_SEED)
        dest = OUT / sid
        dest.mkdir(exist_ok=True)
        layout = rects_to_layout(rects)
        (dest / "layout.json").write_text(json.dumps(layout, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
        kinds = {k: sum(1 for b in blocks if b.kind == k) for k in ("CPU", "GPU", "CIM")}
        meta = {
            "id": sid,
            "title": sch["title"],
            "why": sch["why"],
            "layout": str(dest / "layout.json"),
            "source": "lbfgs",
            "n_chiplets": len(blocks),
            "kinds": kinds,
            "power_total_W": sum(b.power_W for b in blocks),
            "area_total_mm2": sum(b.area_mm2 for b in blocks),
            "outline_mm2": pack["outline_mm2"],
            "atsim_pitch_die_mm": sch["pitch"][0],
            "atsim_pitch_coarse_mm": sch["pitch"][1],
            "atsim_nz": sch["nz"],
            "chips": [
                {"name": b.name, "kind": b.kind, "area_mm2": b.area_mm2, "power_W": b.power_W, "notes": b.notes}
                for b in blocks
            ],
        }
        (dest / "scheme.json").write_text(json.dumps(meta, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
        by_id[sid] = meta
        print(f"   P={meta['power_total_W']:.1f}W area={meta['area_total_mm2']:.0f}mm2", flush=True)

    index_path.write_text(json.dumps(list(by_id.values()), indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print("index n=", len(by_id))


if __name__ == "__main__":
    main()
