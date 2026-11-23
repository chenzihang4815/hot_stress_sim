#!/usr/bin/env python3
"""13-chip package with all 5 CPUs replaced by McPAT official Xeon.

Does not overwrite pkg_5cpu_3gpu_5cim (small TimingSimple CPUs).
"""

from __future__ import annotations

import json
import shutil
import subprocess
import sys
from dataclasses import replace
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path("/home/czh/任务仿真产生功耗文件工具")
SRC = ROOT / "热仿真流水线" / "src"
sys.path.insert(0, str(SRC))

from blocks import BACKUP, Block, accelforge_from_summary, cim_basic_analog, cpu_from_gem5_mcpat  # noqa: E402
from defaults import INITIAL_LAYOUT_SEED  # noqa: E402
from hotspot_io import run_hotspot, write_hotspot_inputs  # noqa: E402
from lbfgs import random_lbfgs_layout  # noqa: E402
from run_large_package import gpu, official_map  # noqa: E402

HS_BIN = ROOT / "HotSpot-master" / "hotspot"
SID = "pkg_5xeon_3gpu_5cim"
OUT = ROOT / "热仿真流水线" / "runs" / SID
SCHEMES = ROOT / "热仿真流水线" / "runs" / "sa_compare_schemes"
BACKUP_OUT = ROOT / "成功运行备份" / "13芯_5Xeon_3GPU_5CIM_2026-10-26"
XEON_TXT = BACKUP / "McPAT/official_templates/Xeon.txt"

CPU_NAMES = [
    ("CPU_matmul", "原 matmul 槽位，现为 McPAT 官方 Xeon（2核 O3+L3，65 nm）"),
    ("CPU_threads", "原 threads 槽位，现为 McPAT 官方 Xeon"),
    ("CPU_armL2", "原 ARM L2 槽位，现为 McPAT 官方 Xeon"),
    ("CPU_atomic", "原 atomic 槽位，现为 McPAT 官方 Xeon"),
    ("CPU_starter", "原 starter 槽位，现为 McPAT 官方 Xeon"),
]


def main() -> None:
    gdir = BACKUP / "Accel-Sim/小面积GPU_backprop"
    cpus = [
        replace(
            cpu_from_gem5_mcpat(XEON_TXT, name=name, notes=notes),
            name=name,
            kind="CPU",
        )
        for name, notes in CPU_NAMES
    ]
    gpus = [
        gpu(gdir / "SM86_RTX3070/accelwattch_power_report.log", "GPU_3070", 392.0, "RTX3070 backprop"),
        gpu(gdir / "SM75_RTX2060_S/accelwattch_power_report.log", "GPU_2060S", 445.0, "RTX2060S backprop"),
        gpu(gdir / "SM6_TITANX/accelwattch_power_report.log", "GPU_TITANX", 471.0, "TitanX backprop"),
    ]
    cims = [
        replace(accelforge_from_summary("jia_jssc_2020+builtin", kind="CIM", name="CIM_jia"), name="CIM_jia", kind="CIM"),
        replace(accelforge_from_summary("colonnade_jssc_2021+builtin", kind="CIM", name="CIM_colonnade"), name="CIM_colonnade", kind="CIM"),
        replace(accelforge_from_summary("wang_vlsi_2022+builtin", kind="CIM", name="CIM_wang"), name="CIM_wang", kind="CIM"),
        replace(cim_basic_analog(), name="CIM_basic", kind="CIM"),
        replace(accelforge_from_summary("sinangil_jssc_2021+builtin", kind="CIM", name="CIM_sinangil"), name="CIM_sinangil", kind="CIM"),
    ]
    blocks = cpus + gpus + cims
    names = [b.name for b in blocks]
    if len(names) != len(set(names)):
        raise SystemExit(f"dup names {names}")

    if OUT.exists():
        shutil.rmtree(OUT)
    OUT.mkdir(parents=True)

    md = OUT / "EXPERIMENT.md"
    lines = [
        "# 13芯：5 颗 CPU 全部换成 McPAT 官方 Xeon",
        "",
        f"- started (UTC): {datetime.now(timezone.utc).isoformat()}",
        "- 旧方案 `pkg_5cpu_3gpu_5cim` **不覆盖**（仍是 TimingSimple 13.6 mm²）",
        "- 新方案 id：`pkg_5xeon_3gpu_5cim`",
        "- 每颗 CPU：Xeon.xml 官方结果 Area=421.676 mm²，P=RuntimeDynamic 66.8485 + Leakage 37.0975 = **103.946 W**",
        "- GPU / CIM 与 13 芯原方案相同",
        "- 初始布局：L-BFGS seed=13",
        "- 名字仍用 CPU_matmul 等槽位，便于和旧布局对照，硬件已不是 gem5 小核",
        "",
        "## 芯粒",
        "",
        "| 名称 | 种类 | 面积 mm² | 功耗 W |",
        "|---|---|---:|---:|",
    ]
    for b in blocks:
        lines.append(f"| {b.name} | {b.kind} | {b.area_mm2:.4g} | {b.power_W:.4g} |")
    md.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("blocks:", flush=True)
    for b in blocks:
        print(f"  {b.name:16s} {b.kind:3s} {b.area_mm2:10.4f} mm2 {b.power_W:10.4f} W", flush=True)

    rects, meta = random_lbfgs_layout(blocks, seed=INITIAL_LAYOUT_SEED)
    write_hotspot_inputs(rects, OUT)
    thermal = run_hotspot(HS_BIN, OUT)
    png = OUT / f"{SID}_hotspot.png"
    official_map(OUT / "floorplan.flp", OUT / "grid.steady", png)

    chips = [
        {
            "name": b.name,
            "kind": b.kind,
            "area_mm2": b.area_mm2,
            "power_W": b.power_W,
            "W_per_cm2": (b.power_W / b.area_mm2 * 100.0) if b.area_mm2 else 0.0,
            "notes": b.notes,
        }
        for b in blocks
    ]
    row = {
        "id": SID,
        "title": "5 Xeon + 3 GPU + 5 CIM（13 颗，CPU 全换大核）",
        "why": "用户要求 5 颗 CPU 全部换成面积能和 GPU 比的硬件；用 McPAT 官方 Xeon 模板，不再用 gem5 TimingSimple",
        "n_chiplets": 13,
        "kinds": {"CPU": 5, "GPU": 3, "CIM": 5},
        "cpu_model": "McPAT official Xeon.xml (2-core O3 + L3, 65 nm)",
        "cpu_area_mm2": cpus[0].area_mm2,
        "cpu_power_W": cpus[0].power_W,
        "power_total_W": sum(b.power_W for b in blocks),
        "area_total_mm2": sum(b.area_mm2 for b in blocks),
        "outline_mm2": meta["outline_mm2"],
        "deadspace_ratio": meta["deadspace_ratio"],
        "overlap_m2": meta["overlap_m2"],
        "HotSpot_Tmax_C": thermal.get("Tmax_C"),
        "HotSpot_Tmin_C": thermal.get("Tmin_C"),
        "layout_algorithm": "lbfgs",
        "layout_seed": INITIAL_LAYOUT_SEED,
        "chips": chips,
        "heatmap": str(png),
        "atsim_pitch_die_mm": 2.5,
        "atsim_pitch_coarse_mm": 10.0,
        "atsim_nz": {
            "nz_die": 2,
            "nz_spreader": 2,
            "nz_sink": 3,
            "nz_tim": 1,
            "nz_substrate": 2,
            "nz_interposer": 2,
        },
    }
    (OUT / "scheme.json").write_text(json.dumps(row, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    # register for later SA without touching old 13-chip
    slot = SCHEMES / SID
    slot.mkdir(parents=True, exist_ok=True)
    shutil.copy2(OUT / "layout.json", slot / "layout.json")
    (slot / "scheme.json").write_text(json.dumps(row, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    idx_path = SCHEMES / "index.json"
    idx = json.loads(idx_path.read_text(encoding="utf-8")) if idx_path.exists() else []
    idx = [e for e in idx if e.get("id") != SID]
    idx.append(
        {
            "id": SID,
            "title": row["title"],
            "why": row["why"],
            "layout": str(slot / "layout.json"),
            "source": str(OUT / "layout.json"),
            "n_chiplets": 13,
            "atsim_pitch_die_mm": 2.5,
            "atsim_pitch_coarse_mm": 10.0,
            "atsim_nz": row["atsim_nz"],
        }
    )
    idx_path.write_text(json.dumps(idx, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    with md.open("a", encoding="utf-8") as f:
        f.write("\n## L-BFGS + HotSpot\n\n")
        f.write(f"- 硅面积合计 **{row['area_total_mm2']:.0f} mm²**（旧 13 芯约 1340 mm² 量级，CPU 只有 ~47 mm²）\n")
        f.write(f"- 封装功耗合计 **{row['power_total_W']:.1f} W**（旧 13 芯约 211 W，CPU 只有 ~5 W）\n")
        f.write(f"- 轮廓 {row['outline_mm2']:.0f} mm²  deadspace={row['deadspace_ratio']:.3f}\n")
        f.write(f"- HotSpot Tmax = **{row['HotSpot_Tmax_C']:.3f} °C**  Tmin={row['HotSpot_Tmin_C']:.3f} °C\n")
        f.write(f"- 热图：`{png}`\n")
        f.write(f"- 布局：`{OUT / 'layout.json'}`\n")
        f.write(f"- 已登记到 `sa_compare_schemes/{SID}` 供后续 SA\n")

    if BACKUP_OUT.exists():
        shutil.rmtree(BACKUP_OUT)
    shutil.copytree(OUT, BACKUP_OUT, ignore=shutil.ignore_patterns("all_layer*", "grid.steady", "all.grid.steady"))
    print(json.dumps({k: row[k] for k in ("id","area_total_mm2","power_total_W","outline_mm2","HotSpot_Tmax_C")}, indent=2))
    print("wrote", OUT)


if __name__ == "__main__":
    main()
