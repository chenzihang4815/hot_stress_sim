#!/usr/bin/env python3
"""Closed loop: s → power → HotSpot → McPAT T → leak → heat-only SA → repeat.

Four large-CPU packages. Layout starts from existing L-BFGS.
Heat-only SA 16 steps, no ATSim.

Leakage is NOT a fresh mcpat binary run. Official templates are device_type=0
(HP): the only temperature dependence in McPAT is I_off_n[g_ip->temp-300]
(cacti/technology.cc). Gate leakage is flat for HP, area and both dynamic
powers do not depend on temperature, and CACTI's partition search ignores
leakage (weight 0). One template run is >1 h and repeats that search to
multiply one number. So each 10 K bin scales the official 380 K subthreshold
by the ITRS HP table ratio and keeps gate fixed — the number the binary emits.
"""

from __future__ import annotations

import json
import math
import random
import re
import shutil
import subprocess
import sys
import time
from copy import deepcopy
from datetime import datetime, timezone
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager

ROOT = Path("/home/czh/任务仿真产生功耗文件工具")
SRC = ROOT / "热仿真流水线" / "src"
sys.path.insert(0, str(SRC))

import sa_thermal_feedback as sa  # noqa: E402
from geometry import packing_metrics  # noqa: E402
from hotspot_io import run_hotspot, write_hotspot_inputs  # noqa: E402

_font_path = "/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc"
font_manager.fontManager.addfont(_font_path)
_prop = font_manager.FontProperties(fname=_font_path)
plt.rcParams["font.family"] = _prop.get_name()
plt.rcParams["axes.unicode_minus"] = False

HS_BIN = ROOT / "HotSpot-master" / "hotspot"
MCPAT_BIN = ROOT / "成功运行备份" / "McPAT" / "bin" / "mcpat"
TPL = ROOT / "成功运行备份" / "McPAT" / "official_templates"
OUT = ROOT / "热仿真流水线" / "runs" / "closed_loop_bigcpu"
BACKUP = ROOT / "成功运行备份" / "大CPU闭环_热McPAT_SA_2026-10-26"

S_LIST = [0.5, 1.0, 1.5, 2.0]
N_SA = 16
SEED = 22
T_LIM_C = 90.0
MAX_INNER = 5
MAX_OUTER = 3
T_MIN_K, T_MAX_K = 300, 400

GPU_STATIC_W = {
    "GPU_3070": 22.5847,
    "GPU_2060S": 24.4467,
    "GPU_TITANX": 29.3748,
}
GPU_CONST_W = {
    "GPU_3070": 32.3252,
    "GPU_2060S": 32.3252,
    "GPU_TITANX": 32.3252,
}

XML = {
    "Penryn": TPL / "Penryn.xml",
    "T2": TPL / "Niagara2.xml",
    "Xeon": TPL / "Xeon.xml",
}
SEED380 = {
    "Penryn": TPL / "Penryn.txt",
    "T2": TPL / "T2.txt",
    "Xeon": TPL / "Xeon.txt",
}

# Official ExampleResults at 380 K. Subthreshold and gate split out so only
# subthreshold is scaled. Numbers match Processor: in the .txt files.
OFFICIAL_380 = {
    "Penryn": {"area_mm2": 91.1348, "dyn_W": 47.9315, "sub_W": 10.7217, "gate_W": 0.588125},
    "T2": {"area_mm2": 273.523, "dyn_W": 50.421, "sub_W": 15.0682, "gate_W": 3.50762},
    "Xeon": {"area_mm2": 421.676, "dyn_W": 66.8485, "sub_W": 35.4132, "gate_W": 1.68429},
}
# HP I_off_n, index = T_K - 300. 45 nm = Penryn, 65 nm = T2 and Xeon.
# cacti/technology.cc I_off_n[0]. Gate I_g_on_n[0] is constant, so not tabulated.
IOFF_HP = {
    "Penryn": {
        0: 2.8e-7, 10: 3.28e-7, 20: 3.81e-7, 30: 4.39e-7, 40: 5.02e-7,
        50: 5.69e-7, 60: 6.42e-7, 70: 7.2e-7, 80: 8.03e-7, 90: 8.91e-7, 100: 9.84e-7,
    },
    "T2": {
        0: 1.96e-7, 10: 2.29e-7, 20: 2.66e-7, 30: 3.05e-7, 40: 3.49e-7,
        50: 3.95e-7, 60: 4.45e-7, 70: 4.97e-7, 80: 5.48e-7, 90: 5.94e-7, 100: 6.3e-7,
    },
    "Xeon": {
        0: 1.96e-7, 10: 2.29e-7, 20: 2.66e-7, 30: 3.05e-7, 40: 3.49e-7,
        50: 3.95e-7, 60: 4.45e-7, 70: 4.97e-7, 80: 5.48e-7, 90: 5.94e-7, 100: 6.3e-7,
    },
}

SCHEMES = [
    {
        "id": "pkg_5penryn_3gpu_5cim",
        "label": "5×Penryn",
        "layout": ROOT / "热仿真流水线" / "runs" / "bigcpu_freq_schemes" / "pkg_5penryn_3gpu_5cim" / "layout.json",
        "cpu_xml": {n: "Penryn" for n in ["CPU_matmul", "CPU_threads", "CPU_armL2", "CPU_atomic", "CPU_starter"]},
    },
    {
        "id": "pkg_5t2_3gpu_5cim",
        "label": "5×T2",
        "layout": ROOT / "热仿真流水线" / "runs" / "bigcpu_freq_schemes" / "pkg_5t2_3gpu_5cim" / "layout.json",
        "cpu_xml": {n: "T2" for n in ["CPU_matmul", "CPU_threads", "CPU_armL2", "CPU_atomic", "CPU_starter"]},
    },
    {
        "id": "pkg_2penryn_2t2_1xeon_3gpu_5cim",
        "label": "2P+2T+1X",
        "layout": ROOT
        / "热仿真流水线"
        / "runs"
        / "bigcpu_freq_schemes"
        / "pkg_2penryn_2t2_1xeon_3gpu_5cim"
        / "layout.json",
        "cpu_xml": {
            "CPU_matmul": "Penryn",
            "CPU_threads": "Penryn",
            "CPU_armL2": "T2",
            "CPU_atomic": "T2",
            "CPU_starter": "Xeon",
        },
    },
    {
        "id": "pkg_5xeon_3gpu_5cim",
        "label": "5×Xeon",
        "layout": ROOT / "热仿真流水线" / "runs" / "pkg_5xeon_3gpu_5cim" / "layout.json",
        "cpu_xml": {n: "Xeon" for n in ["CPU_matmul", "CPU_threads", "CPU_armL2", "CPU_atomic", "CPU_starter"]},
    },
]

HEAT = "#2a78d6"
COUP = "#eb6834"
AQUA = "#1baf7a"
YEL = "#eda100"
COLORS = {
    "pkg_5penryn_3gpu_5cim": HEAT,
    "pkg_5t2_3gpu_5cim": AQUA,
    "pkg_2penryn_2t2_1xeon_3gpu_5cim": YEL,
    "pkg_5xeon_3gpu_5cim": COUP,
}
SURFACE = "#fcfcfb"
TEXT = "#0b0b0b"
MUTED = "#52514e"
GRID = "#e8e7e3"


def log_md(path: Path, line: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8") as f:
        f.write(line.rstrip() + "\n")
    print(line, flush=True)


def round10_k(T_C: float) -> int:
    r = int(round((T_C + 273.15) / 10.0) * 10)
    return max(T_MIN_K, min(T_MAX_K, r))


def parse_processor(text: str) -> dict:
    dyn = leak = area = None
    in_p = False
    for raw in text.splitlines():
        line = raw.strip()
        if line.startswith("Processor:"):
            in_p = True
            continue
        if in_p and (
            line.startswith("Total Cores")
            or line.startswith("Total L2")
            or line.startswith("Total L3")
            or line.startswith("Total First")
            or line.startswith("Total NoC")
            or line.startswith("****")
        ):
            break
        if not in_p:
            continue
        if line.startswith("Area"):
            area = float(line.split("=")[1].split("mm")[0])
        elif line.startswith("Runtime Dynamic"):
            dyn = float(line.split("=")[1].split("W")[0])
        elif line.startswith("Total Leakage"):
            leak = float(line.split("=")[1].split("W")[0])
    if None in (dyn, leak, area):
        raise ValueError("McPAT processor parse failed")
    return {"area_mm2": area, "dyn_W": dyn, "leak_W": leak, "tot_W": dyn + leak}


class McPATCache:
    """Official 380 K result, subthreshold scaled by the HP I_off table.

    Does not launch mcpat. At 380 K the scaled leak equals the official
    Total Leakage, so the seed parse is a check, not an input.
    """

    def __init__(self, work: Path):
        self.work = work
        self.work.mkdir(parents=True, exist_ok=True)
        self.mem: dict[tuple[str, int], dict] = {}
        for name, txt in SEED380.items():
            parsed = parse_processor(txt.read_text(encoding="utf-8", errors="replace"))
            scaled = self._scale(name, 380)
            if abs(parsed["leak_W"] - scaled["leak_W"]) > 1e-3:
                raise ValueError(
                    f"{name} 380 K table leak {scaled['leak_W']:.4f} != official {parsed['leak_W']:.4f}"
                )
            if abs(parsed["dyn_W"] - scaled["dyn_W"]) > 1e-3:
                raise ValueError(f"{name} dyn mismatch")

    def _scale(self, model: str, T_K: int) -> dict:
        off = OFFICIAL_380[model]
        tab = IOFF_HP[model]
        if (T_K - 300) not in tab:
            raise ValueError(f"{model} T={T_K} not on the 10 K ITRS grid")
        ratio = tab[T_K - 300] / tab[80]
        sub = off["sub_W"] * ratio
        leak = sub + off["gate_W"]
        return {
            "area_mm2": off["area_mm2"],
            "dyn_W": off["dyn_W"],
            "sub_W": sub,
            "gate_W": off["gate_W"],
            "leak_W": leak,
            "tot_W": off["dyn_W"] + leak,
            "Ioff_ratio": ratio,
            "T_K": T_K,
            "from_cache": "itrs_hp_table",
        }

    def get(self, model: str, T_K: int) -> dict:
        key = (model, T_K)
        if key not in self.mem:
            parsed = self._scale(model, T_K)
            self.mem[key] = parsed
            slot = self.work / model
            slot.mkdir(exist_ok=True)
            (slot / f"T{T_K}.json").write_text(
                json.dumps(parsed, indent=2) + "\n", encoding="utf-8"
            )
            print(
                f"    McPAT表 {model} T={T_K}  dyn={parsed['dyn_W']:.4g} "
                f"leak={parsed['leak_W']:.4g}  Ioff×{parsed['Ioff_ratio']:.3f}",
                flush=True,
            )
        return self.mem[key]


def apply_power(chips: list[sa.Chip], s: float, cpu_xml: dict[str, str], cache: McPATCache, T_K: int) -> dict:
    cpu_dyn = cpu_leak = gpu_dyn = gpu_leak = gpu_const = cim = 0.0
    per = {}
    for c in chips:
        if c.kind == "CPU":
            m = cache.get(cpu_xml[c.name], T_K)
            dyn = s * m["dyn_W"]
            leak = m["leak_W"]
            c.power_W = dyn + leak
            cpu_dyn += dyn
            cpu_leak += leak
            per[c.name] = {"dyn": dyn, "leak": leak, "model": cpu_xml[c.name], "McPAT_T": T_K}
        elif c.kind == "GPU":
            st = GPU_STATIC_W[c.name]
            cst = GPU_CONST_W[c.name]
            # recover dyn0 from current? use stored s=1 from initial layout power
            # we stash dyn0 on first call via extra attr... use layout original later
            dyn0 = getattr(c, "_dyn0", None)
            if dyn0 is None:
                dyn0 = max(0.0, getattr(c, "_p0", c.power_W) - st - cst)
                c._dyn0 = dyn0
            dyn = s * dyn0
            c.power_W = dyn + st + cst
            gpu_dyn += dyn
            gpu_leak += st
            gpu_const += cst
        else:
            if not hasattr(c, "_p0"):
                c._p0 = c.power_W
            c.power_W = c._p0
            cim += c.power_W
    pkg_dyn = cpu_dyn + gpu_dyn
    pkg_leak = cpu_leak + gpu_leak
    pkg = pkg_dyn + pkg_leak + gpu_const + cim
    cpu_tot = cpu_dyn + cpu_leak
    return {
        "s": s,
        "McPAT_T_K": T_K,
        "cpu_dyn_W": cpu_dyn,
        "cpu_leak_W": cpu_leak,
        "cpu_tot_W": cpu_tot,
        "cpu_dyn_over_leak": cpu_dyn / cpu_leak if cpu_leak else None,
        "cpu_leak_frac": cpu_leak / cpu_tot if cpu_tot else None,
        "pkg_dyn_W": pkg_dyn,
        "pkg_leak_W": pkg_leak,
        "pkg_tot_W": pkg,
        "pkg_dyn_over_leak": pkg_dyn / pkg_leak if pkg_leak else None,
        "per_cpu": per,
    }


def stash_p0(chips: list[sa.Chip]) -> None:
    for c in chips:
        c._p0 = c.power_W
        if c.kind == "GPU":
            st = GPU_STATIC_W[c.name]
            cst = GPU_CONST_W[c.name]
            c._dyn0 = max(0.0, c.power_W - st - cst)


def heat_sa(chips0: list[sa.Chip], work: Path, log) -> tuple[list[sa.Chip], float, dict]:
    if work.exists():
        shutil.rmtree(work)
    work.mkdir(parents=True)
    rng = random.Random(SEED)
    sa.shift_to_origin(chips0)
    m_raw = packing_metrics(sa.to_rects(chips0))
    n_gpu = sum(1 for c in chips0 if c.kind == "GPU")
    intp_w = m_raw["width_m"] + sa.GAP_GPU_M * max(n_gpu, 1) + 0.008
    intp_h = m_raw["height_m"] + sa.GAP_GPU_M * max(n_gpu, 1) + 0.008
    chips = sa.spread_to_gap(deepcopy(chips0), intp_w, intp_h)
    tmax0 = sa.eval_hotspot(chips, work / "eval_init")
    current = deepcopy(chips)
    best = deepcopy(chips)
    tmax_cur = tmax_best = tmax0
    j_cur = j_best = (tmax0 - 45.0) / 25.0
    T = sa.T0
    accepted = 0
    for step in range(1, N_SA + 1):
        trial, op = sa.propose(current, intp_w, intp_h, rng)
        if not sa.legal(trial, intp_w, intp_h):
            T *= sa.ALPHA
            continue
        try:
            tmax_new = sa.eval_hotspot(trial, work / "eval_try")
        except Exception:
            T *= sa.ALPHA
            continue
        j_new = (tmax_new - 45.0) / 25.0
        delta = j_cur - j_new
        ap = 1.0 if delta >= 0 else math.exp(delta / max(T, 1e-9))
        if rng.random() < ap:
            current = trial
            tmax_cur, j_cur = tmax_new, j_new
            accepted += 1
            if j_new < j_best - 1e-12:
                best = deepcopy(trial)
                tmax_best, j_best = tmax_new, j_new
                log(f"    SA step {step:2d} BEST {op} Tmax={tmax_new:.3f}")
        T *= sa.ALPHA
    write_hotspot_inputs(sa.to_rects(best), work / "layout_best")
    run_hotspot(HS_BIN, work / "layout_best")
    return best, tmax_best, {"Tmax_init_C": tmax0, "Tmax_best_C": tmax_best, "accepted": accepted}


def inner_mcpat(chips, s, cpu_xml, cache, T_K, work, log) -> tuple[int, float, dict]:
    T = T_K
    last_br = {}
    tmax = None
    for i in range(MAX_INNER):
        br = apply_power(chips, s, cpu_xml, cache, T)
        tmax = sa.eval_hotspot(chips, work / f"hs_i{i}_T{T}")
        nxt = round10_k(tmax)
        log(f"    inner {i}: McPAT {T} K → Tmax={tmax:.3f} °C → {nxt} K  "
            f"CPU dyn/leak={br['cpu_dyn_over_leak']:.3f} leak_frac={br['cpu_leak_frac']:.3f}")
        last_br = br
        if nxt == T:
            return T, tmax, last_br
        T = nxt
    return T, tmax, last_br


def eval_scheme(spec: dict, cache: McPATCache, md: Path) -> dict:
    sid = spec["id"]
    chips0 = sa.load_chips(spec["layout"])
    stash_p0(chips0)
    cpu_xml = spec["cpu_xml"]
    rows = []
    log_md(md, f"## {spec['label']}  (`{sid}`)")
    for s in S_LIST:
        log_md(md, f"### s={s:g}")
        slot = OUT / sid / f"s_{s}"
        if slot.exists():
            shutil.rmtree(slot)
        slot.mkdir(parents=True)
        chips = deepcopy(chips0)
        T = 380
        tmax = None
        br = {}
        sa_info = None
        for outer in range(MAX_OUTER):
            log_md(md, f"- outer {outer}: McPAT 闭环（布局暂冻）")
            T, tmax, br = inner_mcpat(chips, s, cpu_xml, cache, T, slot / f"outer{outer}_inner", lambda m: log_md(md, m))
            log_md(md, f"- outer {outer}: 只热 SA 16 步")
            chips, tmax_sa, sa_info = heat_sa(chips, slot / f"outer{outer}_sa", lambda m: log_md(md, m))
            nxt = round10_k(tmax_sa)
            log_md(md, f"    SA Tmax {sa_info['Tmax_init_C']:.3f}→{tmax_sa:.3f} °C  next McPAT {nxt} K  accepted {sa_info['accepted']}/{N_SA}")
            tmax = tmax_sa
            if nxt == T:
                log_md(md, f"- **收敛** T 档 {T} K，布局+漏电一致")
                break
            T = nxt
        else:
            log_md(md, f"- 外环用尽，停在 McPAT {T} K")
        # final power/T on SA best
        br = apply_power(chips, s, cpu_xml, cache, T)
        tmax = sa.eval_hotspot(chips, slot / "final_hs")
        rec = {
            "s": s,
            "McPAT_T_K": T,
            "Tmax_C": tmax,
            "in_window": tmax < T_LIM_C,
            "sa": sa_info,
            **br,
        }
        flag = "OK" if rec["in_window"] else f"T={tmax:.1f}>90"
        log_md(
            md,
            f"- 结果 s={s:g}  McPAT={T} K  Tmax={tmax:.2f} °C  "
            f"CPU dyn/leak={br['cpu_dyn_over_leak']:.3f}  leak_frac_cpu={br['cpu_leak_frac']:.3f}  "
            f"P_pkg={br['pkg_tot_W']:.1f} W  {flag}",
        )
        (slot / "eval.json").write_text(json.dumps(rec, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
        (slot / "layout_best.json").write_text(
            json.dumps(sa.snapshot(chips), indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
        )
        rows.append(rec)
        log_md(md, "")
    return {"id": sid, "label": spec["label"], "freq": rows}


def main() -> None:
    if OUT.exists():
        shutil.rmtree(OUT)
    OUT.mkdir(parents=True)
    md = OUT / "EXPERIMENT.md"
    log_md(md, "# 大 CPU 闭环：频率 → 热 → McPAT 漏电 → SA 布局")
    log_md(md, "")
    log_md(md, f"- started (UTC): {datetime.now(timezone.utc).isoformat()}")
    log_md(md, "- 每套、每个 s：McPAT 温度档 ↔ HotSpot 直到 10 K 档稳定，然后 **只热 SA 16 步**，再回写 T；外环最多 3 次")
    log_md(md, "- CPU 动态 = s × McPAT Runtime Dynamic(T)；漏电 = McPAT Total Leakage(T)，不随 s")
    log_md(md, "- GPU 动态 ×s，STATIC 当漏电不乘；CIM 不乘")
    log_md(md, "- 本闭环 **不跑 ATSim**（12 GB）；布局只服务 Tmax/漏电")
    log_md(md, "- 漏电 **不是** 现场跑 mcpat。官方模板 device_type=0（HP），McPAT 里温度只改 I_off（cacti/technology.cc），HP 栅漏电不随温度，面积和动态功耗也不随温度。一档模板要跑 >1 h，重复的是与温度无关的 CACTI 分区搜索。")
    log_md(md, "- 所以每档用官方 380 K ExampleResults：亚阈值漏电 × ITRS HP I_off 表比值，栅漏电固定。380 K 处与官方 Total Leakage 逐位一致（脚本启动时核对）。这就是二进制会打印的数。")
    log_md(md, f"- 窗：Tmax < {T_LIM_C} °C。窗内赢家 = 最低 CPU leak_frac（任务占比最高）")
    log_md(md, "")

    cache = McPATCache(OUT / "mcpat_cache")
    results = []
    for spec in SCHEMES:
        results.append(eval_scheme(spec, cache, md))

    log_md(md, "## 总表")
    log_md(md, "")
    log_md(md, "| 方案 | s | McPAT T | Tmax | CPU dyn/leak | CPU 漏电占比 | 封装 W | 窗 |")
    log_md(md, "|---|---:|---:|---:|---:|---:|---:|---|")
    winners = {}
    for R in results:
        feas = [r for r in R["freq"] if r["in_window"]]
        if feas:
            w = min(feas, key=lambda r: r["cpu_leak_frac"])
            winners[R["id"]] = w
        for r in R["freq"]:
            log_md(
                md,
                f"| {R['label']} | {r['s']:g} | {r['McPAT_T_K']} K | {r['Tmax_C']:.2f} | "
                f"{r['cpu_dyn_over_leak']:.3f} | {r['cpu_leak_frac']:.3f} | {r['pkg_tot_W']:.1f} | "
                f"{'是' if r['in_window'] else '否'} |",
            )
    log_md(md, "")
    log_md(md, "## 窗内更好的方案")
    log_md(md, "")
    if winners:
        best_id = min(winners, key=lambda i: winners[i]["cpu_leak_frac"])
        for R in results:
            if R["id"] in winners:
                w = winners[R["id"]]
                mark = " ← 赢家" if R["id"] == best_id else ""
                log_md(
                    md,
                    f"- {R['label']}: s={w['s']:g}  Tmax={w['Tmax_C']:.2f}  "
                    f"leak_frac={w['cpu_leak_frac']:.3f}  dyn/leak={w['cpu_dyn_over_leak']:.3f}{mark}",
                )
    else:
        log_md(md, "四套在所扫 s 上都没有 Tmax<90。闭环仍然给出了漏电随 T 档变化的曲线，可再往下降 s。")

    # plots
    fig, axes = plt.subplots(1, 2, figsize=(12.4, 4.6))
    fig.set_facecolor(SURFACE)
    for R in results:
        xs = [r["s"] for r in R["freq"]]
        axes[0].plot(xs, [r["Tmax_C"] for r in R["freq"]], "-o", color=COLORS[R["id"]], label=R["label"], lw=2, ms=6)
        axes[1].plot(xs, [r["cpu_dyn_over_leak"] for r in R["freq"]], "-o", color=COLORS[R["id"]], label=R["label"], lw=2, ms=6)
    axes[0].axhline(T_LIM_C, color="#e34948", ls="--", lw=1)
    axes[0].set_title("闭环后 Tmax", loc="left", color=TEXT, fontsize=12)
    axes[1].set_title("闭环后 CPU 任务/漏电（可弯）", loc="left", color=TEXT, fontsize=12)
    for ax, ylab in ((axes[0], "Tmax (°C)"), (axes[1], "CPU P_dyn / P_leak")):
        ax.set_xlabel("s = f / f0", color=MUTED)
        ax.set_ylabel(ylab, color=MUTED)
        ax.set_facecolor(SURFACE)
        ax.yaxis.grid(True, color=GRID, lw=1)
        ax.set_axisbelow(True)
        ax.tick_params(colors=MUTED, length=0)
        ax.legend(frameon=False, fontsize=8)
        for sp in ax.spines.values():
            sp.set_visible(False)
        ax.spines["bottom"].set_visible(True)
        ax.spines["bottom"].set_color(GRID)
    fig.tight_layout()
    fig.savefig(OUT / "closed_tmax_ratio.png", dpi=150, bbox_inches="tight", facecolor=SURFACE)
    plt.close(fig)

    payload = {
        "s_list": S_LIST,
        "n_sa": N_SA,
        "T_lim_C": T_LIM_C,
        "results": results,
        "winners": {k: v for k, v in winners.items()},
        "finished_utc": datetime.now(timezone.utc).isoformat(),
    }
    (OUT / "result.json").write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    if BACKUP.exists():
        shutil.rmtree(BACKUP)
    shutil.copytree(OUT, BACKUP, ignore=shutil.ignore_patterns("outer*", "hs_*", "eval_*", "final_hs", "mcpat_cache"))
    log_md(md, f"日志：`{OUT}`")
    print("WROTE", OUT)


if __name__ == "__main__":
    main()
