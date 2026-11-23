#!/usr/bin/env python3
"""5 CPU + ≥3 GPU + 5 CIM 大封装，L-BFGS，HotSpot 官方热图。"""

from __future__ import annotations

import json
import shutil
import subprocess
import sys
from dataclasses import replace
from pathlib import Path

ROOT = Path("/home/czh/任务仿真产生功耗文件工具")
SRC = ROOT / "热仿真流水线" / "src"
sys.path.insert(0, str(SRC))

from blocks import BACKUP, Block, accelforge_from_summary, cpu_from_gem5_mcpat  # noqa: E402
from defaults import INITIAL_LAYOUT_SEED  # noqa: E402
from hotspot_io import run_hotspot, write_hotspot_inputs  # noqa: E402
from lbfgs import random_lbfgs_layout  # noqa: E402

HS_SCRIPTS = ROOT / "HotSpot-master" / "scripts"
HS_BIN = ROOT / "HotSpot-master" / "hotspot"
OUT = ROOT / "热仿真流水线" / "runs" / "large_package"
BACKUP_OUT = ROOT / "成功运行备份" / "大五颗CPU_多GPU_五CIM_2026-10-21"
ATSIM_CASES = ROOT / "热应力python工程" / "cases" / "large_package"


def tot_power(report: Path) -> float:
    tot = last = None
    for line in report.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith("kernel_avg_power"):
            last = float(line.split("=")[1])
        if line.startswith("gpu_tot_avg_power"):
            tot = float(line.split("=")[1])
    if tot is not None:
        return tot
    if last is not None:
        return last
    raise ValueError(report)


def cpu(txt: Path, name: str, notes: str) -> Block:
    b = cpu_from_gem5_mcpat(txt, name=name, notes=notes)
    return replace(b, name=name, kind="CPU")


def gpu(report: Path, name: str, area: float, notes: str) -> Block:
    p = tot_power(report)
    return Block(name, "GPU", str(report), area, p, notes + f"；P={p:.4g}W 文献面积 {area} mm²")


def cim(combo: str, name: str) -> Block:
    b = accelforge_from_summary(combo, kind="CIM", name=name)
    return replace(b, name=name, kind="CIM")


def official_map(flp: Path, grid: Path, png: Path) -> None:
    work = png.parent
    all_grid = work / "all.grid.steady"
    if grid.resolve() != all_grid.resolve():
        shutil.copy2(grid, all_grid)
    split = subprocess.run(
        [sys.executable, str(HS_SCRIPTS / "split_grid_steady.py"), str(all_grid), "4", "64", "64"],
        cwd=str(work),
        text=True,
        capture_output=True,
        check=False,
    )
    if split.returncode != 0:
        raise RuntimeError(split.stderr + split.stdout)
    mapped = subprocess.run(
        [
            sys.executable,
            str(HS_SCRIPTS / "grid_thermal_map.py"),
            str(flp),
            str(work / "all_layer0.grid.steady"),
            "64",
            "64",
            str(png),
        ],
        cwd=str(work),
        text=True,
        capture_output=True,
        check=False,
    )
    if mapped.returncode != 0:
        raise RuntimeError(mapped.stderr + mapped.stdout)


def write_atsim(scheme_id: str, rects, title: str) -> Path:
    ATSIM_CASES.mkdir(parents=True, exist_ok=True)
    payload = {
        "name": scheme_id,
        "description": title,
        "chips": [
            {
                "name": r.name,
                "x_mm": r.x * 1e3,
                "y_mm": r.y * 1e3,
                "w_mm": r.w * 1e3,
                "h_mm": r.h * 1e3,
                "t_um": 150.0,
                "power_W": r.power_W,
            }
            for r in rects
        ],
        "package": {
            "tim_um": 20,
            "t_spreader_mm": 1.0,
            "t_sink_mm": 6.9,
            "t_substrate_mm": 1.0,
            "t_interposer_um": 100.0,
            "r_convec_K_per_W": 0.1,
            "ambient_C": 45,
            "pitch_die_mm": 1.0,
            "pitch_coarse_mm": 4.0,
            "fill_die": False,
        },
        "layout_algorithm": "lbfgs",
        "layout_seed": INITIAL_LAYOUT_SEED,
    }
    path = ATSIM_CASES / f"{scheme_id}.json"
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return path


def main() -> None:
    gdir = BACKUP / "Accel-Sim/小面积GPU_backprop"
    cpus = [
        cpu(BACKUP / "gem5/x86_tasks/matmul_timing_l2/mcpat.txt", "CPU_matmul", "x86 Timing+L1/L2 矩阵乘"),
        cpu(BACKUP / "gem5/x86_tasks/threads_timing_l2/mcpat.txt", "CPU_threads", "x86 Timing+L1/L2 threads"),
        cpu(BACKUP / "gem5/output/mcpat_from_gem5.txt", "CPU_armL2", "ARM Timing+L1/L2 hello"),
        cpu(BACKUP / "gem5/多案例/02_se_atomic_无cache/mcpat.txt", "CPU_atomic", "ARM Atomic hello"),
        cpu(BACKUP / "gem5/多案例/03_starter_se_atomic/mcpat.txt", "CPU_starter", "ARM starter atomic hello"),
    ]
    gpus3 = [
        gpu(gdir / "SM86_RTX3070/accelwattch_power_report.log", "GPU_3070", 392.0, "RTX3070 backprop 重跑"),
        gpu(gdir / "SM75_RTX2060_S/accelwattch_power_report.log", "GPU_2060S", 445.0, "RTX2060S backprop 重跑"),
        gpu(gdir / "SM6_TITANX/accelwattch_power_report.log", "GPU_TITANX", 471.0, "TitanX backprop 重跑"),
    ]
    gpu3070b = replace(gpus3[0], name="GPU_3070b", notes=gpus3[0].notes + "；第二颗同规格 3070")
    from blocks import cim_basic_analog

    cim_basic = replace(cim_basic_analog(), name="CIM_basic", kind="CIM")
    cims = [
        cim("jia_jssc_2020+builtin", "CIM_jia"),
        cim("colonnade_jssc_2021+builtin", "CIM_colonnade"),
        cim("wang_vlsi_2022+builtin", "CIM_wang"),
        cim_basic,
        cim("sinangil_jssc_2021+builtin", "CIM_sinangil"),
    ]

    schemes = [
        {
            "id": "pkg_5cpu_3gpu_5cim",
            "title": "5 CPU + 3 GPU + 5 CIM（13 颗）",
            "why": "CPU 五种任务/机器，三颗不同面积 GPU，五颗 CIM 密度差很大，给热应力回馈足够旋钮",
            "blocks": cpus + gpus3 + cims,
        },
        {
            "id": "pkg_5cpu_4gpu_5cim",
            "title": "5 CPU + 4 GPU + 5 CIM（14 颗，多一颗同规格 3070）",
            "why": "GPU 不少于 3 且再加一颗相同 3070，可回馈两颗同类 GPU 的间距",
            "blocks": cpus + gpus3 + [gpu3070b] + cims,
        },
    ]

    if OUT.exists():
        shutil.rmtree(OUT)
    OUT.mkdir(parents=True)
    rows = []
    for sch in schemes:
        sid = sch["id"]
        blocks = sch["blocks"]
        names = [b.name for b in blocks]
        if len(names) != len(set(names)):
            raise SystemExit(f"dup names {names}")
        print(f"== {sid} n={len(blocks)}", flush=True)
        for b in blocks:
            print(f"   {b.name:16s} {b.kind:3s} {b.area_mm2:10.4f} mm2 {b.power_W:10.4f} W", flush=True)
        rects, meta = random_lbfgs_layout(blocks, seed=INITIAL_LAYOUT_SEED)
        case_dir = OUT / sid
        case_dir.mkdir()
        write_hotspot_inputs(rects, case_dir)
        thermal = run_hotspot(HS_BIN, case_dir)
        png = case_dir / f"{sid}_hotspot.png"
        official_map(case_dir / "floorplan.flp", case_dir / "grid.steady", png)
        shutil.copy2(png, OUT / f"{sid}.png")
        json_path = write_atsim(sid, rects, sch["title"])
        chips = []
        for b in blocks:
            dens = b.power_W / b.area_mm2 if b.area_mm2 else 0.0
            chips.append(
                {
                    "name": b.name,
                    "kind": b.kind,
                    "area_mm2": b.area_mm2,
                    "power_W": b.power_W,
                    "W_per_cm2": dens * 100.0,
                    "notes": b.notes,
                }
            )
        row = {
            "id": sid,
            "title": sch["title"],
            "why": sch["why"],
            "n_chiplets": len(blocks),
            "kinds": {k: sum(1 for b in blocks if b.kind == k) for k in ("CPU", "GPU", "CIM")},
            "power_total_W": sum(b.power_W for b in blocks),
            "area_total_mm2": sum(b.area_mm2 for b in blocks),
            "outline_mm2": meta["outline_mm2"],
            "deadspace_ratio": meta["deadspace_ratio"],
            "overlap_m2": meta["overlap_m2"],
            "HotSpot_Tmax_C": thermal.get("Tmax_C"),
            "HotSpot_Tmin_C": thermal.get("Tmin_C"),
            "chips": chips,
            "atsim_json": str(json_path),
            "heatmap": str(png),
        }
        rows.append(row)
        (case_dir / "scheme.json").write_text(json.dumps(row, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
        print(
            f"   P={row['power_total_W']:.1f}W area={row['area_total_mm2']:.0f}mm2 "
            f"Tmax={row['HotSpot_Tmax_C']:.2f}C overlap={row['overlap_m2']}",
            flush=True,
        )

    (OUT / "summary.json").write_text(json.dumps(rows, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    if BACKUP_OUT.exists():
        shutil.rmtree(BACKUP_OUT)
    shutil.copytree(OUT, BACKUP_OUT)
    print("wrote", OUT)


if __name__ == "__main__":
    main()
