#!/usr/bin/env python3
"""芯粒面积轴上的温度-应力权衡。

功耗固定，Xeon 和 GPU 的面积从 0.4 倍到 1.2 倍，步长 0.1，共 9 个点。
位置固定：两颗热源贴在左下角，三颗冷的放在右上角。每个点跑一次 ATSim。
"""

from __future__ import annotations

import json
import sys
from copy import deepcopy
from pathlib import Path

SRC = Path("/home/czh/任务仿真产生功耗文件工具/热仿真流水线/src")
sys.path.insert(0, str(SRC))

import run_closed_loop_multidie as cl
import run_closed_loop_suite as suite
import run_fill_superposition as hs

HERE = Path("/home/czh/任务仿真产生功耗文件工具/论文证据/帕累托_布局")
OUT = HERE / "runs" / "area"
BASE = json.loads((HERE / "runs" / "layouts" / "init.json").read_text(encoding="utf-8"))
SCALES = [k / 10 for k in range(4, 13)]


def resized(scale: float) -> list[dict]:
    dies = deepcopy(BASE)
    by = {d["slot"]: d for d in dies}
    for slot in ("Ca", "Ga"):
        die = by[slot]
        die["w"] = max(2, round(die["w"] * scale ** 0.5))
        die["h"] = max(2, round(die["h"] * scale ** 0.5))
        die["area_mm2"] = die.get("area_mm2", 0) * scale
    cpu, gpu = by["Ca"], by["Ga"]
    cpu["j0"], cpu["ib0"] = 0, 0
    gpu["j0"], gpu["ib0"] = cpu["w"] + 1, 0
    x = hs.N
    for slot in ("HBM3", "CIM_jia_m64", "CIM_colonnade_m512"):
        die = by[slot]
        x -= die["w"]
        die["j0"] = x
        die["ib0"] = hs.N - die["h"]
        x -= 1
    return dies


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "atsim").mkdir(exist_ok=True)
    rows = []
    for scale in SCALES:
        tag = f"a{scale:.1f}"
        dies = resized(scale)
        cl.apply_power(dies)
        (OUT / f"{tag}.json").write_text(json.dumps(cl.pack(dies), indent=1) + "\n", encoding="utf-8")
        out_path = OUT / "atsim" / f"{tag}_result.json"
        if out_path.exists():
            result = json.loads(out_path.read_text(encoding="utf-8"))
        else:
            result = suite.run_atsim(dies, "epoxy", OUT / "atsim" / f"{tag}.json", out_path)
        by = {d["slot"]: d for d in dies}
        row = {
            "scale": scale,
            "xeon_mm2": round(by["Ca"]["w"] * by["Ca"]["h"] * 1.265625, 1),
            "gpu_mm2": round(by["Ga"]["w"] * by["Ga"]["h"] * 1.265625, 1),
            "T_C": result["die_cell_max_C"],
            "vm_MPa": result["vm_normal_chip_max_MPa"],
            "vm_name": result["vm_normal_chip_name"],
        }
        rows.append(row)
        print(f"{tag}  Xeon {row['xeon_mm2']:7.1f}  GPU {row['gpu_mm2']:7.1f}  "
              f"T {row['T_C']:.2f}  vm {row['vm_MPa']:.1f} ({row['vm_name']})", flush=True)
        (OUT / "rows.json").write_text(json.dumps(rows, indent=1, ensure_ascii=False) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
