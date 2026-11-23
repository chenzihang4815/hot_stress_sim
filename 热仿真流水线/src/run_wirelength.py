#!/usr/bin/env python3
"""按 TAP-2.5D 的方式算线长，并在它自己的案例上对一下。

TAP-2.5D/fastSA.py 的 compute_wirelength：每条有向连线是两端中心的
曼哈顿距离乘 connection_matrix[i][j]，再除以 wire_count。
wire_count 是矩阵里所有非对角正数之和。坐标单位是 mm，结果是 mm/wire。

这里不导入 TAP-2.5D。它的退火要拉 bstree 和 CPLEX。公式照抄，单独对。
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path("/home/czh/任务仿真产生功耗文件工具")
SRC = ROOT / "热仿真流水线" / "src"
sys.path.insert(0, str(SRC))

import run_closed_loop_multidie as cl
import run_closed_loop_suite as suite

TAP = ROOT / "TAP-2.5D" / "configs"
SAVED = ROOT / "热仿真流水线" / "runs" / "superposition" / "tap_field" / "sfull"
OUT = ROOT / "热仿真流水线" / "runs" / "wirelength"


def wirelength_cells(dies: list[dict], matrix: list[list[int]], cw_mm: float) -> float:
    """芯粒坐标是格。中心换成 mm 后走 TAP 的公式，返回 mm/wire。"""
    placed = [{
        "x_mm": d["j0"] * cw_mm,
        "y_mm": d["ib0"] * cw_mm,
        "w_mm": d["w"] * cw_mm,
        "h_mm": d["h"] * cw_mm,
    } for d in dies]
    return tap_wirelength(placed, matrix)["mm_per_wire"]


def tap_wirelength(dies: list[dict], matrix: list[list[int]]) -> dict:
    """dies 的顺序必须和矩阵的行顺序一致。"""
    if len(dies) != len(matrix):
        raise ValueError(f"{len(dies)} 颗对不上 {len(matrix)} 行")
    total = 0.0
    wire_count = 0
    nets = []
    for i, row in enumerate(matrix):
        for j, weight in enumerate(row):
            if i == j or weight <= 0:
                continue
            a, b = dies[i], dies[j]
            dist = (
                abs(a["x_mm"] + a["w_mm"] / 2 - (b["x_mm"] + b["w_mm"] / 2))
                + abs(a["y_mm"] + a["h_mm"] / 2 - (b["y_mm"] + b["h_mm"] / 2))
            )
            total += dist * weight
            wire_count += weight
            nets.append({"from": i, "to": j, "weight": weight, "manhattan_mm": dist})
    return {
        "mm_per_wire": total / wire_count,
        "total_mm_weight": total,
        "wire_count": wire_count,
        "nets": len(nets),
    }


def parse_cfg(path: Path) -> tuple[list[dict], list[list[int]]]:
    text = path.read_text(encoding="utf-8")
    section = None
    fields: dict[str, str] = {}
    for raw in text.splitlines():
        line = raw.split("#", 1)[0].strip()
        if not line:
            continue
        if line.startswith("[") and line.endswith("]"):
            section = line[1:-1]
            continue
        if "=" not in line:
            if section == "chiplets" and "connections" in fields:
                fields["connections"] += line
            continue
        key, value = (part.strip() for part in line.split("=", 1))
        fields[key] = value
    widths = [float(v) for v in fields["widths"].split(",") if v.strip()]
    heights = [float(v) for v in fields["heights"].split(",") if v.strip()]
    matrix = [
        [int(v) for v in row.split(",") if v.strip()]
        for row in fields["connections"].split(";")
        if row.strip()
    ]
    dies = [
        {"i": i, "w_mm": widths[i], "h_mm": heights[i]}
        for i in range(len(widths))
    ]
    return dies, matrix


def reference() -> None:
    """手算一个 2 颗的例子，确认和 TAP 公式一致。"""
    dies = [
        {"x_mm": 0.0, "y_mm": 0.0, "w_mm": 10.0, "h_mm": 10.0},
        {"x_mm": 20.0, "y_mm": 0.0, "w_mm": 10.0, "h_mm": 6.0},
    ]
    # 中心 (5, 5) 和 (25, 3)，曼哈顿 22。权重 4 和 6，平均仍是 22。
    got = tap_wirelength(dies, [[0, 4], [6, 0]])
    if abs(got["mm_per_wire"] - 22.0) > 1e-9 or got["wire_count"] != 10:
        raise RuntimeError(f"手算对不上: {got}")


def assumed_matrix(dies: list[dict]) -> list[list[int]]:
    n = len(dies)
    matrix = [[0] * n for _ in range(n)]
    hbm = next(i for i, d in enumerate(dies) if d["role"] == "HBM")
    cpus = [i for i, d in enumerate(dies) if d["role"] == "CPU"]
    for i, die in enumerate(dies):
        if die["role"] in ("CPU", "GPU"):
            matrix[i][hbm] = 1024
            matrix[hbm][i] = 1024
        elif die["role"] == "CIM":
            matrix[i][cpus[0]] = 128
            matrix[cpus[0]][i] = 128
    return matrix


def our_cases() -> list[dict]:
    """假设的网表，只为能算。不是 gem5、不是 McPAT 给的连接。"""
    cat = cl.catalog()
    rows = []
    for kind, assumed in (
        ("n5", [[0, 256, 64, 1024, 0],
                [256, 0, 0, 1024, 0],
                [64, 0, 0, 0, 128],
                [1024, 1024, 0, 0, 0],
                [0, 0, 128, 0, 0]]),
        ("n8", None),
        ("n11", None),
    ):
        dies = suite.shelf(cat, kind)
        if assumed is None:
            assumed = assumed_matrix(dies)
        placed = [{
            "name": d["name"],
            "x_mm": d["j0"] * cl.CW_MM,
            "y_mm": d["ib0"] * cl.CW_MM,
            "w_mm": d["w"] * cl.CW_MM,
            "h_mm": d["h"] * cl.CW_MM,
        } for d in dies]
        wl = tap_wirelength(placed, assumed)
        rows.append({
            "case": kind,
            "order": [d["name"] for d in dies],
            "assumed_connection": assumed,
            "mm_per_wire": wl["mm_per_wire"],
            "wire_count": wl["wire_count"],
            "note": "网表是假设的。CPU、GPU 连 HBM，CIM 连最近的一颗 CPU。",
        })
    return rows


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    reference()
    checks = []
    for name in ("ascend910", "micro150", "multigpu"):
        spec, matrix = parse_cfg(TAP / f"sys_{name}.cfg")
        saved = json.loads((SAVED / name / "init.json").read_text(encoding="utf-8"))
        dies = []
        for row in spec:
            hit = next(d for d in saved if d["i"] == row["i"])
            dies.append({**row, "x_mm": hit["x_mm"], "y_mm": hit["y_mm"]})
        wl = tap_wirelength(dies, matrix)
        checks.append({"case": name, "source": "TAP config + 已存 init 坐标", **wl})
        print(f"{name}  {wl['mm_per_wire']:.4f} mm/wire  wires {wl['wire_count']}", flush=True)
    report = {
        "formula": (
            "每条 i->j 是两端中心的曼哈顿距离乘 connection_matrix[i][j]，"
            "全部加起来再除以矩阵正权重之和。和 TAP-2.5D/fastSA.py 的 compute_wirelength 一样。"
            "结果单位 mm/wire。没有跑 TAP 的退火，也没有用 CPLEX。"
        ),
        "hand_check": "两颗中心曼哈顿 22 mm，权重 4 和 6，结果 22 mm/wire。",
        "tap_cases": checks,
        "our_cases_assumed_netlist": our_cases(),
    }
    (OUT / "report.json").write_text(
        json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"wrote {OUT / 'report.json'}", flush=True)


if __name__ == "__main__":
    main()
