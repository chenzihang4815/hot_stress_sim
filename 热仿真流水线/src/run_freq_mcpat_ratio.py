#!/usr/bin/env python3
"""Sweep frequency scale s on frozen 13-chip L*, compare task/leak ratios.

Two leak treatments:
  frozen350 — CPU leak locked at official McPAT 350 K (no T feedback)
  mcpatT    — HotSpot Tmax → McPAT 10 K bin → CPU leak, iterate until bin stable

CPU/GPU dynamic scale with s; CPU leak from McPAT; GPU STATIC=leak (unscaled),
CONST unscaled; CIM analog unscaled. gem5 ops unchanged.
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
OUT = ROOT / "热仿真流水线" / "runs" / "freq_mcpat_ratio_13chip"
BACKUP = ROOT / "成功运行备份" / "频率McPAT占比_13芯_2026-10-25"

S_LIST = [0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 3.0, 4.0]
T_MIN_K, T_MAX_K = 300, 400
MAX_INNER = 6

CPU_XML = {
    "CPU_matmul": XML_DIR / "CPU_matmul_base_T350.xml",
    "CPU_threads": XML_DIR / "CPU_threads_base_T350.xml",
    "CPU_armL2": XML_DIR / "CPU_armL2_base_T350.xml",
    "CPU_atomic": XML_DIR / "CPU_atomic_base_T350.xml",
    "CPU_starter": XML_DIR / "CPU_starter_base_T350.xml",
}
GPU_STATIC_W = {
    "GPU_3070": 22.5847,
    "GPU_3070b": 22.5847,
    "GPU_2060S": 24.4467,
    "GPU_TITANX": 29.3748,
}
GPU_CONST_W = {
    "GPU_3070": 32.3252,
    "GPU_3070b": 32.3252,
    "GPU_2060S": 32.3252,
    "GPU_TITANX": 32.3252,
}


def log_md(path: Path, line: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8") as f:
        f.write(line.rstrip() + "\n")
    print(line, flush=True)


def round10_k(T_C: float) -> int:
    r = int(round((T_C + 273.15) / 10.0) * 10)
    return max(T_MIN_K, min(T_MAX_K, r))


def parse_mcpat_txt(path: Path) -> dict:
    dyn = leak = sub = gate = peak = area = None
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
        elif line.startswith("Peak Power"):
            peak = float(line.split("=")[1].split("W")[0])
    if None in (dyn, leak):
        raise ValueError(f"McPAT parse failed: {path}")
    return {
        "area_mm2": area,
        "dyn_W": dyn,
        "leak_W": leak,
        "sub_W": sub,
        "gate_W": gate,
        "peak_W": peak,
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
    dst.write_text(text2, encoding="utf-8")


class McPATCache:
    def __init__(self, work: Path):
        self.work = work
        self.work.mkdir(parents=True, exist_ok=True)
        self.by_T: dict[int, dict[str, dict]] = {}

    def at(self, T_K: int) -> dict[str, dict]:
        if T_K in self.by_T:
            return self.by_T[T_K]
        slot = self.work / f"T{T_K}"
        slot.mkdir(exist_ok=True)
        out = {}
        for name, xml0 in CPU_XML.items():
            xml = slot / f"{name}.xml"
            txt = slot / f"{name}.txt"
            err = slot / f"{name}.err"
            set_xml_temperature(xml0, xml, T_K)
            cmd = [str(MCPAT_BIN), "-infile", str(xml), "-print_level", "5", "-opt_for_clk", "1"]
            completed = subprocess.run(cmd, text=True, capture_output=True, check=False)
            txt.write_text(completed.stdout, encoding="utf-8")
            err.write_text(completed.stderr, encoding="utf-8")
            if completed.returncode != 0:
                raise RuntimeError(f"McPAT {name} T={T_K} rc={completed.returncode}")
            out[name] = parse_mcpat_txt(txt)
        self.by_T[T_K] = out
        return out


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


def gpu_parts(r: Rect) -> tuple[float, float, float]:
    st = GPU_STATIC_W[r.name]
    cst = GPU_CONST_W[r.name]
    dyn = max(0.0, r.power_W - st - cst)
    return dyn, st, cst


def apply(rects0: list[Rect], s: float, cpu_mcpat: dict[str, dict]) -> tuple[list[Rect], dict]:
    out = deepcopy(rects0)
    cpu_dyn = cpu_leak = gpu_dyn = gpu_leak = gpu_const = cim = 0.0
    for r in out:
        if r.kind == "CPU":
            m = cpu_mcpat[r.name]
            dyn = s * m["dyn_W"]
            leak = m["leak_W"]
            r.power_W = dyn + leak
            cpu_dyn += dyn
            cpu_leak += leak
        elif r.kind == "GPU":
            dyn0, st, cst = gpu_parts(r)
            dyn = s * dyn0
            r.power_W = dyn + st + cst
            gpu_dyn += dyn
            gpu_leak += st
            gpu_const += cst
        else:
            cim += r.power_W
    pkg_dyn = cpu_dyn + gpu_dyn
    pkg_leak = cpu_leak + gpu_leak
    pkg_tot = pkg_dyn + pkg_leak + gpu_const + cim
    cpu_tot = cpu_dyn + cpu_leak
    br = {
        "s": s,
        "cpu_dyn_W": cpu_dyn,
        "cpu_leak_W": cpu_leak,
        "cpu_tot_W": cpu_tot,
        "cpu_dyn_over_leak": (cpu_dyn / cpu_leak) if cpu_leak else None,
        "cpu_leak_frac": (cpu_leak / cpu_tot) if cpu_tot else None,
        "cpu_task_frac": (cpu_dyn / cpu_tot) if cpu_tot else None,
        "gpu_dyn_W": gpu_dyn,
        "gpu_leak_W": gpu_leak,
        "gpu_const_W": gpu_const,
        "cim_W": cim,
        "pkg_dyn_W": pkg_dyn,
        "pkg_leak_W": pkg_leak,
        "pkg_tot_W": pkg_tot,
        "pkg_dyn_over_leak": (pkg_dyn / pkg_leak) if pkg_leak else None,
        "pkg_leak_frac": (pkg_leak / pkg_tot) if pkg_tot else None,
        "pkg_task_frac": (pkg_dyn / pkg_tot) if pkg_tot else None,
    }
    return out, br


def hotspot(rects: list[Rect], work: Path) -> float:
    if work.exists():
        shutil.rmtree(work)
    work.mkdir(parents=True)
    write_hotspot_inputs(rects, work)
    return float(run_hotspot(HS_BIN, work)["Tmax_C"])


def eval_frozen350(rects0, s, cache: McPATCache, slot: Path) -> dict:
    m = cache.at(350)
    rects, br = apply(rects0, s, m)
    tmax = hotspot(rects, slot / "hotspot")
    return {
        "arm": "frozen350",
        "McPAT_T_K": 350,
        "inner_rounds": 1,
        "Tmax_C": tmax,
        "McPAT_T_from_Tmax_K": round10_k(tmax),
        **br,
    }


def eval_mcpatT(rects0, s, cache: McPATCache, slot: Path, t0: int) -> dict:
    T = t0
    last = None
    hist = []
    for i in range(MAX_INNER):
        m = cache.at(T)
        rects, br = apply(rects0, s, m)
        tmax = hotspot(rects, slot / f"hotspot_i{i}_T{T}")
        nxt = round10_k(tmax)
        rec = {"i": i, "T_in": T, "Tmax_C": tmax, "T_out": nxt, **br}
        hist.append(rec)
        last = rec
        if nxt == T:
            break
        T = nxt
    assert last is not None
    return {
        "arm": "mcpatT",
        "McPAT_T_K": last["T_in"],
        "inner_rounds": len(hist),
        "inner": hist,
        "Tmax_C": last["Tmax_C"],
        "McPAT_T_from_Tmax_K": last["T_out"],
        **{k: last[k] for k in last if k not in ("i", "T_in", "T_out")},
    }


def fmt(x, nd=3):
    if x is None:
        return "—"
    return f"{x:.{nd}f}"


def main() -> None:
    if OUT.exists():
        shutil.rmtree(OUT)
    OUT.mkdir(parents=True)
    md = OUT / "EXPERIMENT.md"
    log_md(md, "# 13芯：频率 s × McPAT 温度 → 任务/漏电占比")
    log_md(md, "")
    log_md(md, f"- started (UTC): {datetime.now(timezone.utc).isoformat()}")
    log_md(md, "- 布局冻结：`pkg_5cpu_3gpu_5cim` heat-only SA 最优")
    log_md(md, "- s = f/f0，CPU 与 GPU **动态**同乘 s；CIM 不乘")
    log_md(md, "- CPU 漏电：McPAT Total Leakage（改 XML `temperature`）")
    log_md(md, "- GPU 漏电：AccelWattch STATICP，不随 s、也不进 McPAT")
    log_md(md, "- GPU CONSTP、CIM 功耗固定")
    log_md(md, "- 两条臂：`frozen350` 漏电锁 350 K；`mcpatT` 用 HotSpot Tmax 回写 McPAT 直到 10 K 档稳定")
    log_md(md, "- 看的量：任务/漏电 = P_dyn/P_leak；漏电占比 = P_leak/(P_dyn+P_leak)（CPU）或 /P_pkg")
    log_md(md, "")

    rects0 = load_rects(LAYOUT_BEST)
    cache = McPATCache(OUT / "mcpat_cache")
    rows_f = []
    rows_m = []
    t_guess = 350
    for s in S_LIST:
        log_md(md, f"## s = {s}")
        sf = eval_frozen350(rects0, s, cache, OUT / f"frozen_s{s}")
        sm = eval_mcpatT(rects0, s, cache, OUT / f"mcpat_s{s}", t_guess)
        t_guess = sm["McPAT_T_K"]
        rows_f.append(sf)
        rows_m.append(sm)
        log_md(
            md,
            f"- frozen350: Tmax={sf['Tmax_C']:.3f}C  CPU dyn/leak={sf['cpu_dyn_over_leak']:.3f}  "
            f"leak_frac_cpu={sf['cpu_leak_frac']:.3f}  pkg dyn/leak={sf['pkg_dyn_over_leak']:.3f}",
        )
        log_md(
            md,
            f"- mcpatT T={sm['McPAT_T_K']}K ({sm['inner_rounds']} inner): Tmax={sm['Tmax_C']:.3f}C  "
            f"CPU dyn/leak={sm['cpu_dyn_over_leak']:.3f}  leak_frac_cpu={sm['cpu_leak_frac']:.3f}  "
            f"CPU leak={sm['cpu_leak_W']:.4f}W",
        )
        log_md(md, "")

    def table(title: str, rows: list[dict]) -> None:
        log_md(md, f"## {title}")
        log_md(md, "")
        log_md(
            md,
            "| s | McPAT T | Tmax | CPU dyn | CPU leak | CPU dyn/leak | CPU 漏电占比 | "
            "封装 dyn | 封装 leak | 封装 dyn/leak | 封装漏电占比 |",
        )
        log_md(md, "|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
        for r in rows:
            log_md(
                md,
                f"| {r['s']:g} | {r['McPAT_T_K']} K | {r['Tmax_C']:.2f} °C | "
                f"{r['cpu_dyn_W']:.3f} | {r['cpu_leak_W']:.3f} | {r['cpu_dyn_over_leak']:.3f} | "
                f"{r['cpu_leak_frac']:.3f} | {r['pkg_dyn_W']:.2f} | {r['pkg_leak_W']:.2f} | "
                f"{r['pkg_dyn_over_leak']:.3f} | {r['pkg_leak_frac']:.3f} |",
            )
        log_md(md, "")

    table("frozen350（漏电不随温度）", rows_f)
    table("mcpatT（Tmax 回写 350K）", rows_m)

    log_md(md, "## 升频会不会让任务/漏电比下降？")
    log_md(md, "")
    dips_cpu = []
    dips_pkg = []
    for a, b in zip(rows_m, rows_m[1:]):
        if b["cpu_dyn_over_leak"] < a["cpu_dyn_over_leak"] - 1e-9:
            dips_cpu.append((a["s"], b["s"], a["cpu_dyn_over_leak"], b["cpu_dyn_over_leak"], a["McPAT_T_K"], b["McPAT_T_K"]))
        if b["pkg_dyn_over_leak"] < a["pkg_dyn_over_leak"] - 1e-9:
            dips_pkg.append((a["s"], b["s"], a["pkg_dyn_over_leak"], b["pkg_dyn_over_leak"]))
    if dips_cpu:
        log_md(md, "CPU **出现过** 升频后 dyn/leak 下降（McPAT 温度档跳变、漏电上台阶）：")
        for s0, s1, r0, r1, t0, t1 in dips_cpu:
            log_md(md, f"- s={s0:g}→{s1:g}：CPU dyn/leak {r0:.3f}→{r1:.3f}（McPAT {t0}→{t1} K）")
    else:
        log_md(md, "在本扫点上，CPU dyn/leak **没有** 因升频下降：漏电台阶相对动态的 s 倍增不够大。")
    if dips_pkg:
        log_md(md, "封装 dyn/leak 升频下降：")
        for s0, s1, r0, r1 in dips_pkg:
            log_md(md, f"- s={s0:g}→{s1:g}：{r0:.3f}→{r1:.3f}")
    else:
        log_md(md, "封装 dyn/leak 本扫点上也是升频单调升（GPU 动态跟着 s 涨，STATIC 漏电固定）。")
    log_md(md, "")
    log_md(md, "对照：`frozen350` 里漏电常数，dyn∝s，**dyn/leak 必随 s 单调升**。")
    log_md(md, "`mcpatT` 里漏电按 10 K 档跳，理论上可以在档边界把比压下去；是否发生取决于 Tmax 有没有跨档。")
    log_md(md, "")

    result = {
        "scheme": "pkg_5cpu_3gpu_5cim",
        "layout_frozen": True,
        "s_list": S_LIST,
        "frozen350": rows_f,
        "mcpatT": rows_m,
        "cpu_dyn_over_leak_dips": [
            {"s_from": x[0], "s_to": x[1], "ratio_from": x[2], "ratio_to": x[3], "T_from": x[4], "T_to": x[5]}
            for x in dips_cpu
        ],
        "finished_utc": datetime.now(timezone.utc).isoformat(),
    }
    (OUT / "result.json").write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    if BACKUP.exists():
        shutil.rmtree(BACKUP)
    shutil.copytree(OUT, BACKUP, ignore=shutil.ignore_patterns("hotspot*", "frozen_s*", "mcpat_s*"))
    log_md(md, f"日志：`{OUT}`")
    log_md(md, f"备份：`{BACKUP}`")


if __name__ == "__main__":
    main()
