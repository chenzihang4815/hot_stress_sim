#!/usr/bin/env python3
"""Pre-tapeout chiplet mix vs post-tapeout frequency, stage 1: screen.

Stage 1 enumerates the mix and screens it with HotSpot only, one solve per
(combo, s) on the L-BFGS layout. Real ATSim is 60x slower and is reserved for
the combos the screen puts on the window boundary (stage 2).

The mix is what cannot change after tape-out: CPU model and count, GPU model
and count, CIM count. Frequency s is what can. A combo is "in the window"
when some s in S_LIST keeps HotSpot Tmax under T_LIM_C. Because ATSim runs
5-8 C hotter than HotSpot, the screen uses a band, not the line:

  T_LO = 82 C   below this, ATSim has never crossed 90 in the 120-step chain
  T_HI = 90 C   above this, already out even before ATSim's penalty

so a combo landing between the two needs the real solve before it can be
called either way.
"""

from __future__ import annotations

import json
import shutil
import sys
import time
from datetime import datetime, timezone
from itertools import combinations_with_replacement
from pathlib import Path

ROOT = Path("/home/czh/任务仿真产生功耗文件工具")
SRC = ROOT / "热仿真流水线" / "src"
sys.path.insert(0, str(SRC))

import sa_thermal_feedback as sa  # noqa: E402
from hotspot_io import run_hotspot, write_hotspot_inputs  # noqa: E402
from run_full_chain import (  # noqa: E402
    GPU_AREA, OFFICIAL_380, S_LIST, T_LIM_C, T_START_K, apply_power,
    build_chips, leak_of, round10_k, stash, theta,
)

OUT = ROOT / "热仿真流水线" / "runs" / "combo_screen"
BACKUP = ROOT / "成功运行备份" / "芯粒组合筛选_2026-10-29"

CPU_MODELS = ["T2", "Xeon"]
GPU_MODELS = ["GPU_3070", "GPU_2060S", "GPU_TITANX"]
N_CPU = (1, 2, 3)
N_GPU = (1, 2, 3)
N_CIM = 5

# HotSpot-to-ATSim margin, measured on the 120-step chain: ATSim Tmax runs
# 4.4 to 7.0 C above HotSpot at the same layout and power.
T_LO_C = T_LIM_C - 8.0
MAX_INNER = 5


def combos() -> list[dict]:
    out = []
    for nc in N_CPU:
        for cpus in combinations_with_replacement(CPU_MODELS, nc):
            for ng in N_GPU:
                for gpus in combinations_with_replacement(GPU_MODELS, ng):
                    cpus, gpus = list(cpus), list(gpus)
                    cid = "c%d_%s__g%d_%s" % (
                        nc, "-".join(cpus), ng, "-".join(g.replace("GPU_", "") for g in gpus))
                    out.append({"id": cid, "cpus": cpus, "gpus": gpus, "n_cim": N_CIM})
    return out


def screen_one(spec: dict, work: Path) -> dict:
    """Close the 10 K bin at each s on the fixed L-BFGS layout. HotSpot only."""
    chips, cpu_xml = build_chips(spec)
    stash(chips)
    per_s = []
    for s in S_LIST:
        T = T_START_K
        for i in range(MAX_INNER):
            br = apply_power(chips, s, cpu_xml, T)
            tmax = sa.eval_hotspot(chips, work / f"s{s}_i{i}")
            nxt = round10_k(tmax)
            if nxt == T:
                break
            T = nxt
        per_s.append({
            "s": s, "McPAT_T_K": T, "Tmax_C": tmax,
            "cpu_leak_W": br["cpu_leak_W"], "cpu_dyn_W": br["cpu_dyn_W"],
            "cpu_leak_frac": br["cpu_leak_frac"],
            "cpu_dyn_over_leak": br["cpu_dyn_over_leak"],
            "pkg_tot_W": br["pkg_tot_W"], "N_op": br["N_op"], "Theta_pkg": theta(s),
        })
    shutil.rmtree(work, ignore_errors=True)

    def band(t: float) -> str:
        if t < T_LO_C:
            return "in"
        if t < T_LIM_C:
            return "border"
        return "out"

    for r in per_s:
        r["band"] = band(r["Tmax_C"])
    feasible = [r for r in per_s if r["band"] == "in"]
    best = min(feasible, key=lambda r: r["cpu_leak_frac"]) if feasible else None
    return {
        "id": spec["id"], "cpus": spec["cpus"], "gpus": spec["gpus"],
        "n_cpu": len(spec["cpus"]), "n_gpu": len(spec["gpus"]),
        "cpu_area_mm2": sum(OFFICIAL_380[m]["area_mm2"] for m in spec["cpus"]),
        "gpu_area_mm2": sum(GPU_AREA[g] for g in spec["gpus"]),
        "per_s": per_s,
        "best_s": None if best is None else best["s"],
        "best_leak_frac": None if best is None else best["cpu_leak_frac"],
        "best_theta": None if best is None else best["Theta_pkg"],
        "best_Tmax_C": None if best is None else best["Tmax_C"],
        "s1_band": next(r["band"] for r in per_s if r["s"] == 1.0),
        "s1_Tmax_C": next(r["Tmax_C"] for r in per_s if r["s"] == 1.0),
        "s1_leak_frac": next(r["cpu_leak_frac"] for r in per_s if r["s"] == 1.0),
        "needs_atsim": any(r["band"] == "border" for r in per_s),
    }


def main() -> None:
    specs = combos()
    OUT.mkdir(parents=True, exist_ok=True)
    md = OUT / "EXPERIMENT.md"
    done = {}
    res_path = OUT / "screen.json"
    if res_path.exists():
        done = {r["id"]: r for r in json.loads(res_path.read_text())}
    lines = [
        "# 芯粒组合筛选：流片前的选择 vs 流片后只能调频率",
        "",
        f"- started (UTC): {datetime.now(timezone.utc).isoformat()}",
        f"- 组合：CPU {{T2, Xeon}} × 1..3 颗，GPU {{3070, 2060S, TITANX}} × 1..3 颗，CIM 固定 {N_CIM} 颗。共 {len(specs)} 套。",
        "- 这一阶段只用 HotSpot，布局是 L-BFGS 初始布局，不跑 SA，不跑 ATSim。",
        f"- 分档：HotSpot Tmax < {T_LO_C:.0f} °C 记为窗内（ATSim 比 HotSpot 高 4–7 °C，留 8 °C 余量），"
        f"{T_LO_C:.0f}–{T_LIM_C:.0f} °C 记为边界，需要真 ATSim 才能判定，> {T_LIM_C:.0f} °C 记为窗外。",
        "- 每个 s 都从 350 K 闭合 10 K 温度档，漏电按 ITRS HP 表缩放。",
        "",
    ]
    md.write_text("\n".join(lines), encoding="utf-8")

    rows = []
    t0 = time.time()
    for k, spec in enumerate(specs, 1):
        if spec["id"] in done:
            rows.append(done[spec["id"]])
            continue
        t = time.time()
        rec = screen_one(spec, OUT / "_work" / spec["id"])
        rows.append(rec)
        res_path.write_text(json.dumps(rows, indent=1, ensure_ascii=False) + "\n", encoding="utf-8")
        flag = {("in", "in"): "不降频", ("out", "in"): "需降频", ("out", "out"): "无解"}.get(
            (rec["s1_band"], "in" if rec["best_s"] else "out"), "待定")
        print(f"[{k:2d}/{len(specs)}] {spec['id']:34s} s1 {rec['s1_Tmax_C']:6.2f} {rec['s1_band']:6s} "
              f"best_s {str(rec['best_s']):4s} {flag:4s}  {time.time()-t:.0f}s  "
              f"elapsed {(time.time()-t0)/60:.1f}min", flush=True)

    n_in = sum(r["s1_band"] == "in" for r in rows)
    n_rescue = sum(r["s1_band"] != "in" and r["best_s"] is not None for r in rows)
    n_none = sum(r["best_s"] is None and not r["needs_atsim"] for r in rows)
    n_border = sum(r["needs_atsim"] for r in rows)
    with md.open("a", encoding="utf-8") as f:
        f.write(f"\n## 结果\n\n- 共 {len(rows)} 套，用时 {(time.time()-t0)/60:.1f} min\n")
        f.write(f"- s=1 就在窗内：{n_in} 套。这些不需要流片后降频。\n")
        f.write(f"- s=1 出窗、降频后能进：{n_rescue} 套。\n")
        f.write(f"- 降到 s=0.5 仍出窗：{n_none} 套。\n")
        f.write(f"- 落在边界带、需要真 ATSim 才能判定：{n_border} 套。\n\n")
        f.write("| 组合 | CPU mm² | GPU mm² | s=1 Tmax | s=1 档 | s=1 漏电占比 | 窗内最高 s | 该 s 漏电占比 | 该 s Tmax |\n")
        f.write("|---|---:|---:|---:|---|---:|---:|---:|---:|\n")
        for r in rows:
            f.write(f"| {r['id']} | {r['cpu_area_mm2']:.0f} | {r['gpu_area_mm2']:.0f} | "
                    f"{r['s1_Tmax_C']:.2f} | {r['s1_band']} | {r['s1_leak_frac']:.3f} | "
                    f"{'—' if r['best_s'] is None else f'{r['best_s']:g}'} | "
                    f"{'—' if r['best_leak_frac'] is None else f'{r['best_leak_frac']:.3f}'} | "
                    f"{'—' if r['best_Tmax_C'] is None else f'{r['best_Tmax_C']:.2f}'} |\n")
    if BACKUP.exists():
        shutil.rmtree(BACKUP)
    shutil.copytree(OUT, BACKUP, ignore=shutil.ignore_patterns("_work"))
    print("WROTE", OUT)


if __name__ == "__main__":
    main()
