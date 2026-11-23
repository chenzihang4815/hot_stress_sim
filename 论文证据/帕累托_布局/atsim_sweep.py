#!/usr/bin/env python3
"""对 sweep.py 挑出的布局各跑一次 ATSim，补温度和应力。

已有结果的布局直接读，不重跑。一次一个，OMP 限 2 线程，给主会话的验证留核。
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

SRC = Path("/home/czh/任务仿真产生功耗文件工具/热仿真流水线/src")
sys.path.insert(0, str(SRC))

import run_closed_loop_suite as suite

HERE = Path("/home/czh/任务仿真产生功耗文件工具/论文证据/帕累托_布局")
OUT = HERE / "runs"
FILL = "epoxy"


def main() -> None:
    tags = json.loads((OUT / "chosen.json").read_text(encoding="utf-8"))
    scan = {r["tag"]: r for r in json.loads((OUT / "scan.json").read_text(encoding="utf-8"))}
    (OUT / "atsim").mkdir(exist_ok=True)
    rows = []
    for i, tag in enumerate(tags, 1):
        dies = json.loads((OUT / "layouts" / f"{tag}.json").read_text(encoding="utf-8"))
        case_path = OUT / "atsim" / f"{tag}.json"
        out_path = OUT / "atsim" / f"{tag}_result.json"
        if out_path.exists():
            result = json.loads(out_path.read_text(encoding="utf-8"))
        else:
            result = suite.run_atsim(dies, FILL, case_path, out_path)
        row = {
            "tag": tag,
            "T_py_C": scan[tag]["T_py_C"],
            "grad_py": scan[tag]["grad"],
            "wl_mm": scan[tag]["wl_mm"],
            "T_atsim_C": result["die_cell_max_C"],
            "vm_MPa": result["vm_normal_chip_max_MPa"],
            "vm_name": result["vm_normal_chip_name"],
        }
        rows.append(row)
        print(
            f"{i:02d}/{len(tags)} {tag:8} T {row['T_atsim_C']:.2f}  "
            f"vm {row['vm_MPa']:.1f}  wl {row['wl_mm']:.2f}",
            flush=True,
        )
        (OUT / "atsim_rows.json").write_text(
            json.dumps(rows, indent=1, ensure_ascii=False) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
