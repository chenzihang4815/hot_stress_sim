#!/usr/bin/env python3
"""13-chip closed loop: freeze SA-best layout, feed HotSpot Tmax into McPAT
temperature each round, update CPU leakage, re-run HotSpot until the McPAT
10 K bin stops changing.

gem5 ops / sim time do not change with McPAT T.
GPU / CIM powers stay frozen.
"""

from __future__ import annotations

import json
import re
import shutil
import subprocess
import sys
from copy import deepcopy
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path("/home/czh/任务仿真产生功耗文件工具")
SRC = ROOT / "热仿真流水线" / "src"
sys.path.insert(0, str(SRC))

from geometry import Rect  # noqa: E402
from hotspot_io import run_hotspot, write_hotspot_inputs  # noqa: E402

HS_BIN = ROOT / "HotSpot-master" / "hotspot"
MCPAT_BIN = ROOT / "成功运行备份" / "McPAT" / "bin" / "mcpat"
LAYOUT_BEST = (
    ROOT
    / "热仿真流水线"
    / "runs"
    / "sa_heat_vs_coupled_pkg_5cpu_3gpu_5cim"
    / "heat_only"
    / "layout_best"
    / "layout.json"
)
XML_DIR = ROOT / "热仿真流水线" / "runs" / "mcpat_13chip_tmax_feedback"
OUT = ROOT / "热仿真流水线" / "runs" / "mcpat_13chip_tmax_loop"
BACKUP = ROOT / "成功运行备份" / "McPAT温度闭环_13芯_2026-10-25"

MAX_ROUNDS = 8
T_MIN_K = 300
T_MAX_K = 400

# gem5 stats (do not depend on McPAT T). From the 13-chip one-shot experiment.
GEM5 = {
    "CPU_matmul": {"task": "x86 64x64 matmul Timing+L1/L2", "N_op": 1297322.0, "t_sim_s": 0.002706},
    "CPU_threads": {"task": "x86 threads Timing+L1/L2", "N_op": 1808180.0, "t_sim_s": 0.004124},
    "CPU_armL2": {"task": "ARM Timing+L1/L2 hello", "N_op": 4994.0, "t_sim_s": 3.5e-05},
    "CPU_atomic": {"task": "ARM Atomic hello", "N_op": 5020.0, "t_sim_s": 3e-06},
    "CPU_starter": {"task": "ARM starter atomic hello", "N_op": 5020.0, "t_sim_s": 1e-06},
}
CPU_XML = {
    "CPU_matmul": XML_DIR / "CPU_matmul_base_T350.xml",
    "CPU_threads": XML_DIR / "CPU_threads_base_T350.xml",
    "CPU_armL2": XML_DIR / "CPU_armL2_base_T350.xml",
    "CPU_atomic": XML_DIR / "CPU_atomic_base_T350.xml",
    "CPU_starter": XML_DIR / "CPU_starter_base_T350.xml",
}


def log_md(path: Path, line: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8") as f:
        f.write(line.rstrip() + "\n")
    print(line, flush=True)


def round10_k(T_C: float) -> int:
    k = T_C + 273.15
    r = int(round(k / 10.0) * 10)
    return max(T_MIN_K, min(T_MAX_K, r))


def parse_mcpat_txt(path: Path) -> dict:
    area = dyn = leak = sub = gate = None
    in_processor = False
    for raw in path.read_text(encoding="utf-8", errors="replace").splitlines():
        line = raw.strip()
        if line.startswith("Processor:"):
            in_processor = True
            continue
        if in_processor and (
            line.startswith("Total Cores")
            or line.startswith("Total L2")
            or line.startswith("Total First")
            or line.startswith("Total NoC")
            or line.startswith("****")
        ):
            break
        if not in_processor:
            continue
        if line.startswith("Area"):
            area = float(line.split("=")[1].split("mm")[0])
        elif line.startswith("Runtime Dynamic"):
            dyn = float(line.split("=")[1].split("W")[0])
        elif line.startswith("Total Leakage"):
            leak = float(line.split("=")[1].split("W")[0])
        elif line.startswith("Subthreshold Leakage"):
            sub = float(line.split("=")[1].split("W")[0])
        elif line.startswith("Gate Leakage"):
            gate = float(line.split("=")[1].split("W")[0])
    if None in (dyn, leak):
        raise ValueError(f"McPAT parse failed: {path}")
    return {
        "area_mm2": area,
        "dyn_W": dyn,
        "leak_W": leak,
        "sub_W": sub,
        "gate_W": gate,
        "tot_W": dyn + leak,
    }


def set_xml_temperature(src: Path, dst: Path, T_K: int) -> None:
    text = src.read_text(encoding="utf-8")
    text2, n = re.subn(
        r'(<param name="temperature" value=")(\d+)("/>)',
        rf"\g<1>{T_K}\3",
        text,
        count=1,
    )
    if n != 1:
        raise ValueError(f"temperature param not unique in {src} (n={n})")
    text2 = re.sub(r"T=\d+K", f"T={T_K}K", text2, count=1)
    dst.write_text(text2, encoding="utf-8")


def run_mcpat(xml: Path, txt: Path, err: Path) -> dict:
    cmd = [str(MCPAT_BIN), "-infile", str(xml), "-print_level", "5", "-opt_for_clk", "1"]
    completed = subprocess.run(cmd, text=True, capture_output=True, check=False)
    txt.write_text(completed.stdout, encoding="utf-8")
    err.write_text(completed.stderr, encoding="utf-8")
    if completed.returncode != 0:
        raise RuntimeError(f"McPAT rc={completed.returncode}\n{completed.stderr[-800:]}")
    return parse_mcpat_txt(txt)


def load_rects(path: Path) -> list[Rect]:
    rows = json.loads(path.read_text(encoding="utf-8"))
    return [
        Rect(
            name=r["name"],
            x=r["x_m"],
            y=r["y_m"],
            w=r["w_m"],
            h=r["h_m"],
            power_W=float(r["power_W"]),
            area_mm2=float(r["area_mm2"]),
            kind=r["kind"],
        )
        for r in rows
    ]


def apply_cpu_power(rects: list[Rect], cpu_power: dict[str, float]) -> list[Rect]:
    out = deepcopy(rects)
    for r in out:
        if r.name in cpu_power:
            r.power_W = cpu_power[r.name]
    return out


def hotspot_tmax(rects: list[Rect], work: Path) -> dict:
    if work.exists():
        shutil.rmtree(work)
    work.mkdir(parents=True)
    write_hotspot_inputs(rects, work)
    thermal = run_hotspot(HS_BIN, work)
    return thermal


def main() -> None:
    if OUT.exists():
        shutil.rmtree(OUT)
    OUT.mkdir(parents=True)
    md = OUT / "EXPERIMENT.md"
    log_md(md, "# 13芯 McPAT 温度闭环（冻结 SA 最优布局）")
    log_md(md, "")
    log_md(md, f"- started (UTC): {datetime.now(timezone.utc).isoformat()}")
    log_md(md, "- 封装：`pkg_5cpu_3gpu_5cim`")
    log_md(md, "- 布局：**冻结** heat-only SA 最优（Tmax 原 70.04 °C，功耗来自默认 350 K McPAT）")
    log_md(md, "- 每轮：用本轮 HotSpot **封装 Tmax** 圆整到 McPAT 10 K 档，改 5 颗 CPU 的 XML `temperature`，重跑 McPAT")
    log_md(md, "- 用新的 Runtime Dynamic + Total Leakage 替换 CPU 功耗，**GPU/CIM 不变**，再跑 HotSpot")
    log_md(md, "- 停止：McPAT 温度档与上一轮相同，或达到 8 轮")
    log_md(md, "- gem5 操作数 / 仿真时间 **不随 McPAT T 变**")
    log_md(md, "")

    base_rects = load_rects(LAYOUT_BEST)
    rounds = []
    T_mcpat = 350  # official default, not from layout
    prev_T_bin = None

    for i in range(MAX_ROUNDS):
        rnd_dir = OUT / f"round_{i:02d}_T{T_mcpat}"
        rnd_dir.mkdir()
        log_md(md, f"## Round {i}: McPAT T = {T_mcpat} K")
        cpu_stats = {}
        cpu_power = {}
        leak_sum = dyn_sum = 0.0
        for name, xml0 in CPU_XML.items():
            xml = rnd_dir / f"{name}_T{T_mcpat}.xml"
            txt = rnd_dir / f"{name}_T{T_mcpat}.txt"
            err = rnd_dir / f"{name}_T{T_mcpat}.err"
            set_xml_temperature(xml0, xml, T_mcpat)
            parsed = run_mcpat(xml, txt, err)
            g = GEM5[name]
            tot = parsed["tot_W"]
            j_per_op = tot * g["t_sim_s"] / g["N_op"] if g["N_op"] else None
            ops_per_j = (1.0 / j_per_op) if j_per_op else None
            leak_frac = parsed["leak_W"] / tot if tot else None
            rec = {
                "name": name,
                "task": g["task"],
                "N_op": g["N_op"],
                "t_sim_s": g["t_sim_s"],
                "ops_per_s": g["N_op"] / g["t_sim_s"],
                "McPAT_T_K": T_mcpat,
                **parsed,
                "leak_frac": leak_frac,
                "J_per_op": j_per_op,
                "ops_per_J": ops_per_j,
            }
            cpu_stats[name] = rec
            cpu_power[name] = tot
            leak_sum += parsed["leak_W"]
            dyn_sum += parsed["dyn_W"]
            log_md(
                md,
                f"- {name}: dyn={parsed['dyn_W']:.6g} W  leak={parsed['leak_W']:.6g} W  "
                f"tot={tot:.6g} W  leak_frac={leak_frac:.3f}  N_op={g['N_op']:.0f} (unchanged)",
            )

        rects = apply_cpu_power(base_rects, cpu_power)
        p_cpu = sum(r.power_W for r in rects if r.kind == "CPU")
        p_gpu = sum(r.power_W for r in rects if r.kind == "GPU")
        p_cim = sum(r.power_W for r in rects if r.kind == "CIM")
        p_tot = p_cpu + p_gpu + p_cim
        hs = hotspot_tmax(rects, rnd_dir / "hotspot")
        tmax_C = float(hs["Tmax_C"])
        next_T = round10_k(tmax_C)
        rec_round = {
            "round": i,
            "McPAT_T_in_K": T_mcpat,
            "cpu": cpu_stats,
            "P_cpu_dyn_W": dyn_sum,
            "P_cpu_leak_W": leak_sum,
            "P_cpu_W": p_cpu,
            "P_gpu_W": p_gpu,
            "P_cim_W": p_cim,
            "P_pkg_W": p_tot,
            "Tmax_C": tmax_C,
            "Tmax_K": float(hs["Tmax_K"]),
            "McPAT_T_out_K": next_T,
            "converged": next_T == T_mcpat,
        }
        rounds.append(rec_round)
        log_md(
            md,
            f"- 封装功耗 CPU={p_cpu:.4f} W (dyn {dyn_sum:.4f}+leak {leak_sum:.4f})  "
            f"GPU={p_gpu:.3f} W  CIM={p_cim:.4f} W  tot={p_tot:.3f} W",
        )
        log_md(md, f"- HotSpot Tmax = **{tmax_C:.3f} °C** ({tmax_C+273.15:.2f} K) → 下一档 McPAT **{next_T} K**")
        if next_T == T_mcpat:
            log_md(md, f"- **收敛**：下一档仍是 {T_mcpat} K，停止。")
            break
        if prev_T_bin is not None and next_T == prev_T_bin:
            # oscillation 350↔340 etc: still continue one more? user asked replace each round.
            pass
        prev_T_bin = T_mcpat
        T_mcpat = next_T
        log_md(md, "")

    # summary table
    log_md(md, "")
    log_md(md, "## 各轮对照")
    log_md(md, "")
    log_md(
        md,
        "| 轮 | McPAT T | CPU dyn | CPU leak | CPU 合计 | 封装合计 | HotSpot Tmax | 下一档 | 操作数 |",
    )
    log_md(md, "|---:|---:|---:|---:|---:|---:|---:|---:|---|")
    for r in rounds:
        log_md(
            md,
            f"| {r['round']} | {r['McPAT_T_in_K']} K | {r['P_cpu_dyn_W']:.4f} | {r['P_cpu_leak_W']:.4f} | "
            f"{r['P_cpu_W']:.4f} | {r['P_pkg_W']:.3f} | {r['Tmax_C']:.3f} °C | {r['McPAT_T_out_K']} K | 五颗均不变 |",
        )
    log_md(md, "")
    log_md(md, "## 规律（写给下一轮查阅）")
    if len(rounds) >= 2:
        d_leak = rounds[-1]["P_cpu_leak_W"] - rounds[0]["P_cpu_leak_W"]
        d_t = rounds[-1]["Tmax_C"] - rounds[0]["Tmax_C"]
        log_md(
            md,
            f"1. 相对官方 350 K：末轮 CPU 漏电 {d_leak:+.4f} W，封装 Tmax {d_t:+.4f} K。",
        )
    else:
        log_md(md, "1. 只跑了一轮（输入档与输出档相同）。")
    log_md(md, "2. 操作数来自 gem5，McPAT 温度改不了 `simInsts` / `simSeconds`。")
    log_md(md, "3. Runtime Dynamic 若各轮相同，说明 McPAT 动态功耗不吃 temperature。")
    log_md(md, "4. GPU/CIM 功耗未改；封装 Tmax 主要由 GPU 决定，CPU 漏电变化对 Tmax 可能很小。")
    log_md(md, "")
    result = {
        "scheme": "pkg_5cpu_3gpu_5cim",
        "layout": str(LAYOUT_BEST),
        "layout_frozen": True,
        "McPAT_T0_K": 350,
        "note": "Closed loop: HotSpot Tmax -> McPAT temperature (10K grid) -> CPU leak+dyn -> HotSpot. gem5 ops unchanged.",
        "rounds": rounds,
        "finished_utc": datetime.now(timezone.utc).isoformat(),
    }
    (OUT / "result.json").write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    log_md(md, f"日志目录：`{OUT}`")
    if BACKUP.exists():
        shutil.rmtree(BACKUP)
    shutil.copytree(OUT, BACKUP)
    log_md(md, f"备份：`{BACKUP}`")
    print("WROTE", OUT)


if __name__ == "__main__":
    main()
