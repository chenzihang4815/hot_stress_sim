#!/usr/bin/env python3
"""多颗 CPU/GPU/CIM 组成一套封装，给后续热应力回馈用。热图用 HotSpot 官方脚本。"""

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
OUT = ROOT / "热仿真流水线" / "runs" / "multi_chiplet_schemes"
BACKUP_OUT = ROOT / "成功运行备份" / "多CPU_GPU_CIM方案热图_2026-10-21"
ATSIM_CASES = ROOT / "热应力python工程" / "cases" / "multi_chiplet_schemes"


def parse_tot_power(report: Path) -> float:
    tot = last_k = None
    for line in report.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith("kernel_avg_power"):
            last_k = float(line.split("=")[1])
        if line.startswith("gpu_tot_avg_power"):
            tot = float(line.split("=")[1])
    if tot is not None:
        return tot
    if last_k is not None:
        return last_k
    raise ValueError(f"no power in {report}")


def named_cpu(txt: Path, name: str, notes: str) -> Block:
    b = cpu_from_gem5_mcpat(txt, name=name, notes=notes)
    return replace(b, name=name, kind="CPU")


def named_gpu(report: Path, name: str, area_mm2: float, notes: str) -> Block:
    p = parse_tot_power(report)
    return Block(name, "GPU", str(report), area_mm2, p, notes + f"；P={p:.4g} W；文献面积 {area_mm2} mm²")


def named_cim(combo: str, name: str) -> Block:
    b = accelforge_from_summary(combo, kind="CIM", name=name)
    return replace(b, name=name, kind="CIM")


def split_and_map(flp: Path, grid: Path, png: Path) -> None:
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
    layer0 = work / "all_layer0.grid.steady"
    mapped = subprocess.run(
        [
            sys.executable,
            str(HS_SCRIPTS / "grid_thermal_map.py"),
            str(flp),
            str(layer0),
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


def write_atsim_json(scheme_id: str, rects, title: str) -> Path:
    ATSIM_CASES.mkdir(parents=True, exist_ok=True)
    chips = [
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
    ]
    payload = {
        "name": scheme_id,
        "description": title,
        "chips": chips,
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
    gpu_dir = BACKUP / "Accel-Sim/小面积GPU_backprop"
    cpu_matmul = named_cpu(
        BACKUP / "gem5/x86_tasks/matmul_timing_l2/mcpat.txt",
        "CPU_matmul",
        "x86 Timing+L1/L2，64x64 matmul",
    )
    cpu_threads = named_cpu(
        BACKUP / "gem5/x86_tasks/threads_timing_l2/mcpat.txt",
        "CPU_threads",
        "x86 Timing+L1/L2，gem5 threads",
    )
    cpu_starter = named_cpu(
        BACKUP / "gem5/多案例/03_starter_se_atomic/mcpat.txt",
        "CPU_starter",
        "ARM starter_se atomic，hello（高功率密度对照）",
    )
    gpu_3070 = named_gpu(
        gpu_dir / "SM86_RTX3070/accelwattch_power_report.log",
        "GPU_3070",
        392.0,
        "RTX3070，backprop 按该配置重跑 AccelWattch",
    )
    gpu_titanx = named_gpu(
        gpu_dir / "SM6_TITANX/accelwattch_power_report.log",
        "GPU_TITANX",
        471.0,
        "Pascal Titan X，backprop 按该配置重跑",
    )
    gpu_2060 = named_gpu(
        gpu_dir / "SM75_RTX2060_S/accelwattch_power_report.log",
        "GPU_2060S",
        445.0,
        "RTX2060 Super，backprop 按该配置重跑",
    )
    cim_jia = named_cim("jia_jssc_2020+builtin", "CIM_jia")
    cim_col = named_cim("colonnade_jssc_2021+builtin", "CIM_colonnade")
    cim_wang = named_cim("wang_vlsi_2022+builtin", "CIM_wang")

    schemes = [
        {
            "id": "pkg_2cpu_1gpu_2cim",
            "title": "2 CPU (x86 matmul + threads) + RTX3070 + 2 CIM (jia 大而冷, colonnade 小而热)",
            "why": "后续回馈的基线：两颗 CPU 任务不同、两颗 CIM 密度差一个数量级、GPU 用较小的 392 mm²",
            "blocks": [cpu_matmul, cpu_threads, gpu_3070, cim_jia, cim_col],
        },
        {
            "id": "pkg_3cpu_1gpu_3cim",
            "title": "3 CPU (x86 matmul, threads, ARM starter) + TitanX + 3 CIM (jia/colonnade/wang)",
            "why": "更多频率/架构旋钮：三颗 CPU 密度不同；三颗 CIM 面积和功率密度都不同",
            "blocks": [cpu_matmul, cpu_threads, cpu_starter, gpu_titanx, cim_jia, cim_col, cim_wang],
        },
        {
            "id": "pkg_2cpu_2gpu_2cim",
            "title": "2 CPU (x86 matmul + ARM starter) + 2 GPU (3070 + 2060S) + 2 CIM (jia + wang)",
            "why": "两颗不同面积/功耗的 GPU，布局回馈可以挪 GPU 相对位置；CPU 一颗带 L2、一颗小核高密度",
            "blocks": [cpu_matmul, cpu_starter, gpu_3070, gpu_2060, cim_jia, cim_wang],
        },
    ]

    if OUT.exists():
        shutil.rmtree(OUT)
    OUT.mkdir(parents=True)

    summaries = []
    for sch in schemes:
        sid = sch["id"]
        blocks = sch["blocks"]
        names = [b.name for b in blocks]
        if len(names) != len(set(names)):
            raise SystemExit(f"duplicate names in {sid}: {names}")
        print(f"== {sid}  n={len(blocks)}  {names}", flush=True)
        rects, meta = random_lbfgs_layout(blocks, seed=INITIAL_LAYOUT_SEED)
        case_dir = OUT / sid
        case_dir.mkdir()
        write_hotspot_inputs(rects, case_dir)
        thermal = run_hotspot(HS_BIN, case_dir)
        png = case_dir / f"{sid}_hotspot.png"
        split_and_map(case_dir / "floorplan.flp", case_dir / "grid.steady", png)
        shutil.copy2(png, OUT / f"{sid}.png")
        json_path = write_atsim_json(sid, rects, sch["title"])
        chips = []
        for b in blocks:
            dens = b.power_W / b.area_mm2 if b.area_mm2 else 0.0
            chips.append(
                {
                    "name": b.name,
                    "kind": b.kind,
                    "area_mm2": b.area_mm2,
                    "power_W": b.power_W,
                    "W_per_mm2": dens,
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
        summaries.append(row)
        (case_dir / "scheme.json").write_text(
            json.dumps(row, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
        )
        print(
            f"   P={row['power_total_W']:.2f}W  area={row['area_total_mm2']:.1f}mm2  "
            f"Tmax={row['HotSpot_Tmax_C']:.2f}C  overlap={row['overlap_m2']}",
            flush=True,
        )

    (OUT / "summary.json").write_text(
        json.dumps(summaries, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    if BACKUP_OUT.exists():
        shutil.rmtree(BACKUP_OUT)
    shutil.copytree(OUT, BACKUP_OUT)
    print("wrote", OUT)


if __name__ == "__main__":
    main()
