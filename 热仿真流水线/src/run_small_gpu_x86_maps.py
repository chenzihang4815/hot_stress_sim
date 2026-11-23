#!/usr/bin/env python3
"""x86 gem5 tasks + smaller GPUs (recomputed AccelWattch) + CIM jia → L-BFGS → HotSpot official maps."""

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
from geometry import Rect  # noqa: E402
from hotspot_io import run_hotspot, write_hotspot_inputs  # noqa: E402
from lbfgs import random_lbfgs_layout  # noqa: E402

HS_SCRIPTS = ROOT / "HotSpot-master" / "scripts"
HS_BIN = ROOT / "HotSpot-master" / "hotspot"
PIPE = ROOT / "热仿真流水线"
OUT = PIPE / "runs" / "x86_smallgpu_3chip"
BACKUP_OUT = ROOT / "成功运行备份" / "x86与小面积GPU三芯热图_2026-10-21"


def parse_tot_power(report: Path) -> float:
    tot = None
    last_k = None
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


def gpu_block(name: str, report: Path, area_mm2: float, notes: str) -> Block:
    p = parse_tot_power(report)
    return Block(
        name="GPU",
        kind="GPU",
        source=str(report),
        area_mm2=area_mm2,
        power_W=p,
        notes=notes + f"；AccelWattch gpu_tot_avg_power={p:.4g} W；面积为文献裸片 {area_mm2} mm²",
    )


def split_and_map(flp: Path, grid: Path, png: Path) -> None:
    work = png.parent
    all_grid = work / "all.grid.steady"
    if grid.resolve() != all_grid.resolve():
        shutil.copy2(grid, all_grid)
    flp_copy = work / "floorplan.flp"
    if flp.resolve() != flp_copy.resolve():
        shutil.copy2(flp, flp_copy)
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
            str(work / "floorplan.flp"),
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


def main() -> None:
    if OUT.exists():
        shutil.rmtree(OUT)
    OUT.mkdir(parents=True)

    cpu_matmul = replace(
        cpu_from_gem5_mcpat(
            BACKUP / "gem5/x86_tasks/matmul_timing_l2/mcpat.txt",
            name="CPU",
            notes="x86 TimingSimpleCPU + L1/L2，64x64 matmul（静态 gcc -static）",
        ),
        name="CPU",
        kind="CPU",
    )
    cpu_threads = replace(
        cpu_from_gem5_mcpat(
            BACKUP / "gem5/x86_tasks/threads_timing_l2/mcpat.txt",
            name="CPU",
            notes="x86 TimingSimpleCPU + L1/L2，gem5 自带 threads",
        ),
        name="CPU",
        kind="CPU",
    )
    cim = replace(
        accelforge_from_summary("jia_jssc_2020+builtin", kind="CIM", name="CIM"),
        name="CIM",
        kind="CIM",
    )

    gpu_dir = BACKUP / "Accel-Sim/小面积GPU_backprop"
    gpus = [
        gpu_block(
            "TITANX",
            gpu_dir / "SM6_TITANX/accelwattch_power_report.log",
            471.0,
            "Pascal Titan X / SM6_TITANX，同一份 backprop 轨迹重跑功耗",
        ),
        gpu_block(
            "RTX2060S",
            gpu_dir / "SM75_RTX2060_S/accelwattch_power_report.log",
            445.0,
            "Turing RTX 2060 Super / SM75_RTX2060_S，同一份 backprop 重跑功耗",
        ),
        gpu_block(
            "RTX3070",
            gpu_dir / "SM86_RTX3070/accelwattch_power_report.log",
            392.0,
            "Ampere RTX 3070 / SM86_RTX3070，同一份 backprop 重跑功耗（AccelWattch XML 仍是 Volta 标定）",
        ),
    ]

    mixes = []
    for gpu in gpus:
        mixes.append((f"x86matmul__{gpu.source.split('/')[-2]}", cpu_matmul, gpu, cim))
    mixes.append(("x86threads__SM86_RTX3070", cpu_threads, gpus[-1], cim))

    rows = []
    for mix_id, cpu, gpu, cim_b in mixes:
        blocks = [cpu, gpu, cim_b]
        print("==", mix_id, flush=True)
        print(
            f"   CPU {cpu.area_mm2:.3f}mm2 {cpu.power_W:.3f}W | "
            f"GPU {gpu.area_mm2:.0f}mm2 {gpu.power_W:.3f}W | "
            f"CIM {cim_b.area_mm2:.3f}mm2 {cim_b.power_W*1e3:.2f}mW",
            flush=True,
        )
        rects, meta = random_lbfgs_layout(blocks, seed=INITIAL_LAYOUT_SEED)
        case_dir = OUT / mix_id
        case_dir.mkdir()
        write_hotspot_inputs(rects, case_dir)
        thermal = run_hotspot(HS_BIN, case_dir)
        png = case_dir / f"{mix_id}_hotspot.png"
        split_and_map(case_dir / "floorplan.flp", case_dir / "grid.steady", png)
        shutil.copy2(png, OUT / f"{mix_id}.png")
        row = {
            "mix": mix_id,
            "CPU_task": cpu.notes,
            "CPU_mm2": cpu.area_mm2,
            "CPU_W": cpu.power_W,
            "GPU": gpu.notes,
            "GPU_mm2": gpu.area_mm2,
            "GPU_W": gpu.power_W,
            "CIM": "jia_jssc_2020",
            "CIM_mm2": cim_b.area_mm2,
            "CIM_W": cim_b.power_W,
            "outline_mm2": meta["outline_mm2"],
            "HotSpot_Tmax_C": thermal.get("Tmax_C"),
            "HotSpot_Tmin_C": thermal.get("Tmin_C"),
            "map": str(png),
        }
        rows.append(row)
        (case_dir / "combo.json").write_text(json.dumps(row, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    (OUT / "summary.json").write_text(json.dumps(rows, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    if BACKUP_OUT.exists():
        shutil.rmtree(BACKUP_OUT)
    shutil.copytree(OUT, BACKUP_OUT)
    print(json.dumps(rows, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
