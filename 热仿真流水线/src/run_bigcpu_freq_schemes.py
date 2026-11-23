#!/usr/bin/env python3
"""Large-CPU packages (Penryn / Niagara2 / Xeon / mix) + frequency sweep.

Does not overwrite the small-TimingSimple 13-chip runs.
Xeon-only layout is reused if already built.
"""

from __future__ import annotations

import json
import shutil
import sys
from copy import deepcopy
from dataclasses import replace
from datetime import datetime, timezone
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager

ROOT = Path("/home/czh/任务仿真产生功耗文件工具")
SRC = ROOT / "热仿真流水线" / "src"
sys.path.insert(0, str(SRC))

from blocks import BACKUP, accelforge_from_summary, cim_basic_analog, cpu_from_gem5_mcpat  # noqa: E402
from defaults import INITIAL_LAYOUT_SEED  # noqa: E402
from geometry import Rect  # noqa: E402
from hotspot_io import run_hotspot, write_hotspot_inputs  # noqa: E402
from lbfgs import random_lbfgs_layout  # noqa: E402
from run_large_package import gpu, official_map  # noqa: E402

_font_path = "/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc"
font_manager.fontManager.addfont(_font_path)
_prop = font_manager.FontProperties(fname=_font_path)
plt.rcParams["font.family"] = _prop.get_name()
plt.rcParams["axes.unicode_minus"] = False

HS_BIN = ROOT / "HotSpot-master" / "hotspot"
OUT = ROOT / "热仿真流水线" / "runs" / "bigcpu_freq_schemes"
BACKUP_OUT = ROOT / "成功运行备份" / "大CPU多方案频率对照_2026-10-26"
SCHEMES_DIR = ROOT / "热仿真流水线" / "runs" / "sa_compare_schemes"
TPL = BACKUP / "McPAT/official_templates"

S_LIST = [0.5, 0.75, 1.0, 1.25, 1.5, 2.0]
T_LIM_C = 90.0

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

HEAT = "#2a78d6"
COUP = "#eb6834"
AQUA = "#1baf7a"
YEL = "#eda100"
SURFACE = "#fcfcfb"
TEXT = "#0b0b0b"
MUTED = "#52514e"
GRID = "#e8e7e3"


def log_md(path: Path, line: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8") as f:
        f.write(line.rstrip() + "\n")
    print(line, flush=True)


def parse_mcpat(txt: Path) -> tuple[float, float, float]:
    b = cpu_from_gem5_mcpat(txt, name="tmp", notes="x")
    # recover dyn/leak from notes is fragile; reparse
    dyn = leak = None
    in_p = False
    for raw in txt.read_text(encoding="utf-8", errors="replace").splitlines():
        line = raw.strip()
        if line.startswith("Processor:"):
            in_p = True
            continue
        if in_p and (
            line.startswith("Total Cores")
            or line.startswith("Total L2")
            or line.startswith("Total First")
            or line.startswith("Total NoC")
            or line.startswith("****")
        ):
            break
        if not in_p:
            continue
        if line.startswith("Runtime Dynamic"):
            dyn = float(line.split("=")[1].split("W")[0])
        elif line.startswith("Total Leakage"):
            leak = float(line.split("=")[1].split("W")[0])
    assert dyn is not None and leak is not None
    return b.area_mm2, dyn, leak


def cpu_block(txt: Path, name: str, notes: str):
    b = cpu_from_gem5_mcpat(txt, name=name, notes=notes)
    return replace(b, name=name, kind="CPU")


def make_gpus_cims():
    gdir = BACKUP / "Accel-Sim/小面积GPU_backprop"
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
    return gpus, cims


def build_layout(sid: str, title: str, why: str, cpus, gpus, cims, pitch=(2.5, 10.0)) -> dict:
    blocks = list(cpus) + list(gpus) + list(cims)
    names = [b.name for b in blocks]
    if len(names) != len(set(names)):
        raise SystemExit(f"dup {names}")
    slot = OUT / sid
    if slot.exists():
        shutil.rmtree(slot)
    slot.mkdir(parents=True)
    print(f"== L-BFGS {sid}", flush=True)
    for b in blocks:
        print(f"   {b.name:16s} {b.kind:3s} {b.area_mm2:10.3f} {b.power_W:10.3f}", flush=True)
    rects, meta = random_lbfgs_layout(blocks, seed=INITIAL_LAYOUT_SEED)
    write_hotspot_inputs(rects, slot)
    thermal = run_hotspot(HS_BIN, slot)
    png = slot / f"{sid}_hotspot.png"
    official_map(slot / "floorplan.flp", slot / "grid.steady", png)
    row = {
        "id": sid,
        "title": title,
        "why": why,
        "n_chiplets": len(blocks),
        "kinds": {k: sum(1 for b in blocks if b.kind == k) for k in ("CPU", "GPU", "CIM")},
        "power_total_W": sum(b.power_W for b in blocks),
        "area_total_mm2": sum(b.area_mm2 for b in blocks),
        "cpu_area_mm2": sum(b.area_mm2 for b in cpus),
        "cpu_power_W": sum(b.power_W for b in cpus),
        "outline_mm2": meta["outline_mm2"],
        "deadspace_ratio": meta["deadspace_ratio"],
        "HotSpot_Tmax_C": thermal.get("Tmax_C"),
        "HotSpot_Tmin_C": thermal.get("Tmin_C"),
        "layout_algorithm": "lbfgs",
        "layout_seed": INITIAL_LAYOUT_SEED,
        "layout": str(slot / "layout.json"),
        "heatmap": str(png),
        "atsim_pitch_die_mm": pitch[0],
        "atsim_pitch_coarse_mm": pitch[1],
        "chips": [
            {
                "name": b.name,
                "kind": b.kind,
                "area_mm2": b.area_mm2,
                "power_W": b.power_W,
                "notes": b.notes,
            }
            for b in blocks
        ],
    }
    (slot / "scheme.json").write_text(json.dumps(row, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    reg = SCHEMES_DIR / sid
    reg.mkdir(parents=True, exist_ok=True)
    shutil.copy2(slot / "layout.json", reg / "layout.json")
    (reg / "scheme.json").write_text(json.dumps(row, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return row


def load_rects(path: Path) -> list[Rect]:
    rows = json.loads(path.read_text(encoding="utf-8"))
    return [
        Rect(r["name"], r["x_m"], r["y_m"], r["w_m"], r["h_m"], float(r["power_W"]), float(r["area_mm2"]), r["kind"])
        for r in rows
    ]


def split_cpu(r: Rect, cpu_parts: dict) -> tuple[float, float]:
    dyn, leak = cpu_parts[r.name]
    return dyn, leak


def split_gpu(r: Rect) -> tuple[float, float, float]:
    st = GPU_STATIC_W[r.name]
    cst = GPU_CONST_W[r.name]
    dyn = max(0.0, r.power_W - st - cst)
    return dyn, st, cst


def apply_s(rects0: list[Rect], s: float, cpu_parts: dict) -> tuple[list[Rect], dict]:
    out = deepcopy(rects0)
    cpu_dyn = cpu_leak = gpu_dyn = gpu_leak = gpu_const = cim = 0.0
    for r in out:
        if r.kind == "CPU":
            dyn0, leak = split_cpu(r, cpu_parts)
            dyn = s * dyn0
            r.power_W = dyn + leak
            cpu_dyn += dyn
            cpu_leak += leak
        elif r.kind == "GPU":
            dyn0, st, cst = split_gpu(r)
            dyn = s * dyn0
            r.power_W = dyn + st + cst
            gpu_dyn += dyn
            gpu_leak += st
            gpu_const += cst
        else:
            cim += r.power_W
    pkg_dyn = cpu_dyn + gpu_dyn
    pkg_leak = cpu_leak + gpu_leak
    pkg = pkg_dyn + pkg_leak + gpu_const + cim
    cpu_tot = cpu_dyn + cpu_leak
    return out, {
        "s": s,
        "cpu_dyn_W": cpu_dyn,
        "cpu_leak_W": cpu_leak,
        "cpu_tot_W": cpu_tot,
        "cpu_dyn_over_leak": cpu_dyn / cpu_leak if cpu_leak else None,
        "cpu_leak_frac": cpu_leak / cpu_tot if cpu_tot else None,
        "gpu_dyn_W": gpu_dyn,
        "gpu_leak_W": gpu_leak,
        "pkg_dyn_W": pkg_dyn,
        "pkg_leak_W": pkg_leak,
        "pkg_tot_W": pkg,
        "pkg_dyn_over_leak": pkg_dyn / pkg_leak if pkg_leak else None,
        "pkg_leak_frac": pkg_leak / pkg if pkg else None,
    }


def hotspot(rects: list[Rect], work: Path) -> float:
    if work.exists():
        shutil.rmtree(work)
    work.mkdir(parents=True)
    write_hotspot_inputs(rects, work)
    return float(run_hotspot(HS_BIN, work)["Tmax_C"])


def cpu_parts_from_layout(rects: list[Rect], models: dict[str, tuple[float, float, float]]) -> dict:
    """name -> (dyn, leak) using model lookup by area match or provided map."""
    out = {}
    for r in rects:
        if r.kind != "CPU":
            continue
        if r.name in models:
            _a, dyn, leak = models[r.name]
            out[r.name] = (dyn, leak)
            continue
        # fallback: power_W at s=1 is dyn+leak; need stored
        raise KeyError(r.name)
    return out


def register_index(entries: list[dict]) -> None:
    idx_path = SCHEMES_DIR / "index.json"
    idx = json.loads(idx_path.read_text(encoding="utf-8")) if idx_path.exists() else []
    ids = {e["id"] for e in entries}
    idx = [e for e in idx if e.get("id") not in ids]
    for row in entries:
        idx.append(
            {
                "id": row["id"],
                "title": row["title"],
                "why": row["why"],
                "layout": str(SCHEMES_DIR / row["id"] / "layout.json"),
                "source": row["layout"],
                "n_chiplets": row["n_chiplets"],
                "atsim_pitch_die_mm": row["atsim_pitch_die_mm"],
                "atsim_pitch_coarse_mm": row["atsim_pitch_coarse_mm"],
            }
        )
    idx_path.write_text(json.dumps(idx, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def main() -> None:
    if OUT.exists():
        shutil.rmtree(OUT)
    OUT.mkdir(parents=True)
    md = OUT / "EXPERIMENT.md"
    log_md(md, "# 大 CPU 多方案 × 升/降频")
    log_md(md, "")
    log_md(md, f"- started (UTC): {datetime.now(timezone.utc).isoformat()}")
    log_md(md, "- GPU/CIM 与 13 芯相同；CPU 换成 McPAT 官方大核模板（不是 gem5 任务）")
    log_md(md, "- Penryn 91.1 mm²，P=47.93+11.31=59.24 W；Niagara2/T2 273.5 mm²，P=50.42+18.58=69.00 W；Xeon 421.7 mm²，P=66.85+37.10=103.95 W")
    log_md(md, "- 没有真正 500+ mm² 的官方模板；混搭里最大的一颗用 Xeon 422 mm²")
    log_md(md, "- s=f/f0，CPU/GPU **动态** ×s，漏电与 GPU CONST/STATIC、CIM 不乘")
    log_md(md, "- 布局冻结后只重跑 HotSpot（12 GB，本对照不做 ATSim）")
    log_md(md, f"- 窗口记录：Tmax < {T_LIM_C} °C（仅标注，不作为唯一目标）")
    log_md(md, "")

    area_p, dyn_p, leak_p = parse_mcpat(TPL / "Penryn.txt")
    area_t, dyn_t, leak_t = parse_mcpat(TPL / "T2.txt")
    area_x, dyn_x, leak_x = parse_mcpat(TPL / "Xeon.txt")
    gpus, cims = make_gpus_cims()

    names5 = ["CPU_matmul", "CPU_threads", "CPU_armL2", "CPU_atomic", "CPU_starter"]

    def five(txt, tag):
        return [cpu_block(txt, n, f"McPAT official {tag}") for n in names5]

    schemes_spec = []

    # 1 all Penryn
    schemes_spec.append(
        (
            "pkg_5penryn_3gpu_5cim",
            "5×Penryn + 3 GPU + 5 CIM",
            "五颗都换成 91 mm² Penryn",
            five(TPL / "Penryn.txt", "Penryn"),
            {n: (area_p, dyn_p, leak_p) for n in names5},
        )
    )
    # 2 all T2
    schemes_spec.append(
        (
            "pkg_5t2_3gpu_5cim",
            "5×Niagara2 + 3 GPU + 5 CIM",
            "五颗都换成 274 mm² Niagara2/T2",
            five(TPL / "T2.txt", "Niagara2/T2"),
            {n: (area_t, dyn_t, leak_t) for n in names5},
        )
    )
    # 3 mix 2 Penryn + 2 T2 + 1 Xeon
    mix_cpus = [
        cpu_block(TPL / "Penryn.txt", "CPU_matmul", "Penryn 91 mm²"),
        cpu_block(TPL / "Penryn.txt", "CPU_threads", "Penryn 91 mm²"),
        cpu_block(TPL / "T2.txt", "CPU_armL2", "Niagara2 274 mm²"),
        cpu_block(TPL / "T2.txt", "CPU_atomic", "Niagara2 274 mm²"),
        cpu_block(TPL / "Xeon.txt", "CPU_starter", "Xeon 422 mm²"),
    ]
    mix_models = {
        "CPU_matmul": (area_p, dyn_p, leak_p),
        "CPU_threads": (area_p, dyn_p, leak_p),
        "CPU_armL2": (area_t, dyn_t, leak_t),
        "CPU_atomic": (area_t, dyn_t, leak_t),
        "CPU_starter": (area_x, dyn_x, leak_x),
    }
    schemes_spec.append(
        (
            "pkg_2penryn_2t2_1xeon_3gpu_5cim",
            "2×Penryn + 2×T2 + 1×Xeon + 3 GPU + 5 CIM",
            "用户点的混搭：两颗91、两颗270、一颗420",
            mix_cpus,
            mix_models,
        )
    )

    built = []
    models_by_sid = {}
    for sid, title, why, cpus, models in schemes_spec:
        row = build_layout(sid, title, why, cpus, gpus, cims)
        built.append(row)
        models_by_sid[sid] = models
        log_md(
            md,
            f"- 建成 `{sid}`  CPU面积={row['cpu_area_mm2']:.0f} mm²  CPU功耗={row['cpu_power_W']:.1f} W  "
            f"封装 {row['power_total_W']:.1f} W / {row['area_total_mm2']:.0f} mm²  "
            f"Tmax={row['HotSpot_Tmax_C']:.2f} °C",
        )

    # reuse xeon layout
    xeon_dir = ROOT / "热仿真流水线" / "runs" / "pkg_5xeon_3gpu_5cim"
    xeon_scheme = json.loads((xeon_dir / "scheme.json").read_text(encoding="utf-8"))
    built.append(
        {
            **xeon_scheme,
            "id": "pkg_5xeon_3gpu_5cim",
            "layout": str(xeon_dir / "layout.json"),
        }
    )
    models_by_sid["pkg_5xeon_3gpu_5cim"] = {n: (area_x, dyn_x, leak_x) for n in names5}
    log_md(
        md,
        f"- 复用 `{xeon_scheme['id']}`  CPU面积={xeon_scheme['cpu_area_mm2']:.0f} mm²  "
        f"Tmax={xeon_scheme['HotSpot_Tmax_C']:.2f} °C",
    )
    register_index(built[:-1])

    # frequency sweep
    log_md(md, "")
    log_md(md, "## 频率扫描")
    log_md(md, "")
    all_freq = {}
    for row in built:
        sid = row["id"]
        rects0 = load_rects(Path(row["layout"]))
        parts = {n: (m[1], m[2]) for n, m in models_by_sid[sid].items()}
        recs = []
        log_md(md, f"### {sid}")
        for s in S_LIST:
            scaled, br = apply_s(rects0, s, parts)
            tmax = hotspot(scaled, OUT / "freq" / sid / f"s_{s}")
            rec = {**br, "Tmax_C": tmax, "in_window": tmax < T_LIM_C}
            recs.append(rec)
            flag = "OK" if rec["in_window"] else f"T={tmax:.1f}>90"
            log_md(
                md,
                f"- s={s:g}  Tmax={tmax:.2f} °C  P_pkg={br['pkg_tot_W']:.1f} W  "
                f"CPU dyn/leak={br['cpu_dyn_over_leak']:.3f}  leak_frac_cpu={br['cpu_leak_frac']:.3f}  {flag}",
            )
        feas = [r for r in recs if r["in_window"]]
        if feas:
            best = max(feas, key=lambda r: r["s"])
            log_md(md, f"- 窗内最高 s = **{best['s']:g}**（Tmax={best['Tmax_C']:.2f} °C）")
        else:
            log_md(md, "- 本扫点 **没有** s 落在 Tmax<90（s=0.5 仍超）")
        all_freq[sid] = recs
        log_md(md, "")

    # summary table
    log_md(md, "## 总表（s=1）")
    log_md(md, "")
    log_md(md, "| 方案 | CPU 模型 | CPU 面积 | CPU 功耗 s=1 | 封装功耗 | Tmax s=1 | 窗内最高 s |")
    log_md(md, "|---|---|---:|---:|---:|---:|---:|")
    labels = {
        "pkg_5penryn_3gpu_5cim": "5×Penryn",
        "pkg_5t2_3gpu_5cim": "5×T2",
        "pkg_2penryn_2t2_1xeon_3gpu_5cim": "2P+2T+1X",
        "pkg_5xeon_3gpu_5cim": "5×Xeon",
    }
    for row in built:
        sid = row["id"]
        recs = all_freq[sid]
        r1 = next(r for r in recs if abs(r["s"] - 1) < 1e-9)
        feas = [r for r in recs if r["in_window"]]
        smax = max(r["s"] for r in feas) if feas else None
        log_md(
            md,
            f"| {labels[sid]} | {row['title'].split('+')[0].strip()} | {row['cpu_area_mm2']:.0f} | "
            f"{r1['cpu_tot_W']:.1f} | {r1['pkg_tot_W']:.1f} | {r1['Tmax_C']:.1f} | "
            f"{smax if smax is not None else '无'} |",
        )

    # plots
    fig, axes = plt.subplots(1, 2, figsize=(12.4, 4.6))
    fig.set_facecolor(SURFACE)
    colors = {
        "pkg_5penryn_3gpu_5cim": HEAT,
        "pkg_5t2_3gpu_5cim": AQUA,
        "pkg_2penryn_2t2_1xeon_3gpu_5cim": YEL,
        "pkg_5xeon_3gpu_5cim": COUP,
    }
    for row in built:
        sid = row["id"]
        recs = all_freq[sid]
        xs = [r["s"] for r in recs]
        axes[0].plot(xs, [r["Tmax_C"] for r in recs], "-o", color=colors[sid], label=labels[sid], lw=2, ms=6)
        axes[1].plot(xs, [r["cpu_dyn_over_leak"] for r in recs], "-o", color=colors[sid], label=labels[sid], lw=2, ms=6)
    axes[0].axhline(T_LIM_C, color="#e34948", ls="--", lw=1, label="90 °C")
    axes[0].set_title("升/降频 · 封装 Tmax", loc="left", color=TEXT, fontsize=12)
    axes[1].set_title("升/降频 · CPU 任务/漏电", loc="left", color=TEXT, fontsize=12)
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
    fig.savefig(OUT / "freq_tmax_and_ratio.png", dpi=150, bbox_inches="tight", facecolor=SURFACE)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8.2, 4.4))
    fig.set_facecolor(SURFACE)
    x = list(range(len(built)))
    t1 = [next(r["Tmax_C"] for r in all_freq[row["id"]] if abs(r["s"] - 1) < 1e-9) for row in built]
    labs = [labels[row["id"]] for row in built]
    ax.bar(x, t1, 0.55, color=[colors[row["id"]] for row in built], zorder=3)
    ax.axhline(T_LIM_C, color="#e34948", ls="--", lw=1)
    ax.set_xticks(x)
    ax.set_xticklabels(labs, color=TEXT)
    ax.set_ylabel("Tmax at s=1 (°C)", color=MUTED)
    ax.set_title("标称频率下四套大 CPU 封装的 Tmax", loc="left", color=TEXT, fontsize=12)
    ax.set_facecolor(SURFACE)
    ax.yaxis.grid(True, color=GRID, lw=1)
    ax.set_axisbelow(True)
    ax.tick_params(colors=MUTED, length=0)
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.spines["bottom"].set_visible(True)
    ax.spines["bottom"].set_color(GRID)
    for i, v in enumerate(t1):
        ax.text(i, v, f"{v:.1f}", ha="center", va="bottom", color=TEXT, fontsize=9)
    fig.tight_layout()
    fig.savefig(OUT / "s1_tmax.png", dpi=150, bbox_inches="tight", facecolor=SURFACE)
    plt.close(fig)

    payload = {
        "schemes": [{k: row[k] for k in row if k != "chips"} | {"chips": row.get("chips")} for row in built],
        "freq": all_freq,
        "s_list": S_LIST,
        "T_lim_C": T_LIM_C,
        "finished_utc": datetime.now(timezone.utc).isoformat(),
        "note": "CPU dyn/leak from official McPAT templates; GPU STATIC=leak unscaled; layout frozen.",
    }
    (OUT / "result.json").write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    log_md(md, "")
    log_md(md, f"图：`{OUT / 'freq_tmax_and_ratio.png'}`  `{OUT / 's1_tmax.png'}`")
    if BACKUP_OUT.exists():
        shutil.rmtree(BACKUP_OUT)
    shutil.copytree(OUT, BACKUP_OUT, ignore=shutil.ignore_patterns("freq", "all_layer*", "grid.steady", "all.grid.steady"))
    print("wrote", OUT)


if __name__ == "__main__":
    main()
