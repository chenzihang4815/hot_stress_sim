#!/usr/bin/env python3
"""Aggregate heat-only vs coupled SA across all schemes and draw comparison charts."""

from __future__ import annotations

import json
import shutil
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
import numpy as np

_font_path = "/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc"
font_manager.fontManager.addfont(_font_path)
_prop = font_manager.FontProperties(fname=_font_path)
plt.rcParams["font.family"] = _prop.get_name()
plt.rcParams["axes.unicode_minus"] = False

ROOT = Path("/home/czh/任务仿真产生功耗文件工具")
OUT = ROOT / "热仿真流水线" / "runs" / "sa_heat_vs_coupled_multi"
BACKUP = ROOT / "成功运行备份" / "SA热vs热应力多方案对照_2026-10-23"

HEAT = "#2a78d6"
COUP = "#eb6834"
SURFACE = "#fcfcfb"
TEXT = "#0b0b0b"
MUTED = "#52514e"
GRID = "#e8e7e3"
ZERO = "#8a8983"

SCHEMES = [
    {"id": "pkg_1cpu_1gpu_1cim", "label": "3芯", "short": "3芯 1C+1G+1M"},
    {"id": "pkg_2cpu_1gpu_1cim", "label": "4芯", "short": "4芯 2C+1G+1M"},
    {"id": "pkg_2cpu_1gpu_2cim", "label": "5芯", "short": "5芯 2C+1G+2M"},
    {"id": "pkg_2cpu_2gpu_2cim", "label": "6芯", "short": "6芯 2C+2G+2M"},
    {"id": "pkg_3cpu_1gpu_3cim", "label": "7芯", "short": "7芯 3C+1G+3M"},
    {"id": "pkg_4cpu_2gpu_4cim", "label": "10芯", "short": "10芯 4C+2G+4M"},
    {"id": "pkg_4cpu_3gpu_4cim", "label": "11芯", "short": "11芯 4C+3G+4M"},
    {"id": "pkg_4cpu_3gpu_5cim", "label": "12芯", "short": "12芯 4C+3G+5M"},
    {"id": "pkg_5cpu_3gpu_5cim", "label": "13芯", "short": "13芯 5C+3G+5M"},
    {"id": "pkg_5cpu_4gpu_5cim", "label": "14芯", "short": "14芯 5C+4G+5M"},
    {"id": "pkg_5cpu_5gpu_5cim", "label": "15芯", "short": "15芯 5C+5G+5M"},
    {"id": "pkg_5cpu_5gpu_6cim", "label": "16芯", "short": "16芯 5C+5G+6M"},
]


def _vm_best(arm: dict) -> float:
    return float(arm.get("vm_chip_max_best_MPa") if arm.get("vm_chip_max_best_MPa") is not None else arm["vm_mean_best_MPa"])


def _vm_init(arm: dict) -> float:
    return float(arm.get("vm_chip_max_init_MPa") if arm.get("vm_chip_max_init_MPa") is not None else arm["vm_mean_init_MPa"])


def _dvm(d: dict) -> float:
    if d.get("d_vm_chip_max_MPa") is not None:
        return float(d["d_vm_chip_max_MPa"])
    return float(d["d_vm_mean_MPa"])


def load_rows() -> list[dict]:
    rows = []
    for s in SCHEMES:
        path = ROOT / "热仿真流水线" / "runs" / f"sa_heat_vs_coupled_{s['id']}" / "compare_table.json"
        t = json.loads(path.read_text(encoding="utf-8"))
        h, c = t["heat_only"], t["coupled"]
        d = t["delta_coupled_minus_heat"]
        rows.append(
            {
                **s,
                "n_chiplets": t.get("n_chiplets"),
                "n_steps": t.get("n_steps", 16),
                "heat": h,
                "coupled": c,
                "delta": d,
                "vm_h": _vm_best(h),
                "vm_c": _vm_best(c),
                "vm_h_init": _vm_init(h),
                "vm_c_init": _vm_init(c),
                "dvm": _dvm(d),
                "vm_name_h": h.get("vm_chip_name_best"),
                "vm_name_c": c.get("vm_chip_name_best"),
                "wall_s": t.get("wall_s"),
                "metric": t.get("vm_metric"),
            }
        )
    return rows


def _style_ax(ax):
    ax.set_facecolor(SURFACE)
    ax.yaxis.grid(True, color=GRID, linewidth=1, zorder=0)
    ax.set_axisbelow(True)
    ax.tick_params(colors=MUTED, length=0)
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.spines["bottom"].set_visible(True)
    ax.spines["bottom"].set_color(GRID)


def grouped_bars(ax, labels, heat_vals, coup_vals, ylabel, title):
    x = np.arange(len(labels))
    w = 0.36
    gap = 0.03
    b1 = ax.bar(x - (w + gap) / 2, heat_vals, w, color=HEAT, label="只热 SA", zorder=3)
    b2 = ax.bar(x + (w + gap) / 2, coup_vals, w, color=COUP, label="热+应力 SA", zorder=3)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, color=TEXT, fontsize=8)
    ax.set_ylabel(ylabel, color=MUTED, fontsize=9)
    ax.set_title(title, loc="left", color=TEXT, fontsize=11, pad=8)
    _style_ax(ax)
    ax.legend(frameon=False, fontsize=8, loc="upper left")
    return b1, b2


def delta_bars(ax, labels, deltas, ylabel, title, pos_note, neg_note):
    x = np.arange(len(labels))
    colors = [COUP if v < -1e-9 else HEAT if v > 1e-9 else "#b4b2a9" for v in deltas]
    ax.bar(x, deltas, 0.55, color=colors, zorder=3)
    ax.axhline(0.0, color=ZERO, linewidth=1, zorder=2)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, color=TEXT, fontsize=8)
    ax.set_ylabel(ylabel, color=MUTED, fontsize=9)
    ax.set_title(title, loc="left", color=TEXT, fontsize=11, pad=8)
    _style_ax(ax)
    ax.text(0.0, 1.02, f"{neg_note}  ·  {pos_note}", transform=ax.transAxes, fontsize=8, color=MUTED, ha="left")
    for i, v in enumerate(deltas):
        if abs(v) < 0.005:
            continue
        va = "bottom" if v >= 0 else "top"
        ax.text(i, v, f"{v:+.2f}", ha="center", va=va, fontsize=7, color=TEXT)


def main() -> None:
    rows = load_rows()
    OUT.mkdir(parents=True, exist_ok=True)
    labels = [r["label"] for r in rows]

    # 1) best layout grouped
    fig, axes = plt.subplots(1, 3, figsize=(16.4, 4.8))
    fig.set_facecolor(SURFACE)
    grouped_bars(
        axes[0], labels,
        [r["heat"]["Tmax_best_C"] for r in rows],
        [r["coupled"]["Tmax_best_C"] for r in rows],
        "Tmax (°C)",
        "A  选出的布局 · 峰值温度",
    )
    grouped_bars(
        axes[1], labels,
        [r["vm_h"] for r in rows],
        [r["vm_c"] for r in rows],
        "最差芯粒内部 von Mises (MPa)",
        "B  选出的布局 · 峰值应力",
    )
    grouped_bars(
        axes[2], labels,
        [r["heat"]["w_pv_best_um"] for r in rows],
        [r["coupled"]["w_pv_best_um"] for r in rows],
        "中介层翘曲峰–谷 (μm)",
        "C  选出的布局 · 翘曲",
    )
    fig.suptitle(
        "只热 SA vs 热+应力 SA  ·  12 套封装，16 步，seed=22",
        color=TEXT, fontsize=13, x=0.01, ha="left", y=1.03,
    )
    fig.tight_layout()
    fig.savefig(OUT / "compare_best_grouped.png", dpi=150, bbox_inches="tight", facecolor=SURFACE)
    plt.close(fig)

    # 2) delta
    fig, axes = plt.subplots(1, 3, figsize=(16.4, 4.6))
    fig.set_facecolor(SURFACE)
    delta_bars(
        axes[0], labels, [r["delta"]["dTmax_C"] for r in rows],
        "ΔTmax (°C)", "D  热+应力 − 只热  ·  Tmax",
        "正 = 更热", "负 = 更凉",
    )
    delta_bars(
        axes[1], labels, [r["dvm"] for r in rows],
        "Δvm (MPa)", "E  热+应力 − 只热  ·  峰值应力",
        "正 = 应力更高", "负 = 应力更低（应力进成本的收益）",
    )
    delta_bars(
        axes[2], labels, [r["delta"]["d_w_pv_um"] for r in rows],
        "Δw_pv (μm)", "F  热+应力 − 只热  ·  翘曲",
        "正 = 翘曲更大", "负 = 翘曲更小",
    )
    fig.tight_layout()
    fig.savefig(OUT / "compare_delta.png", dpi=150, bbox_inches="tight", facecolor=SURFACE)
    plt.close(fig)

    # 3) vs init — the intuitive story: heat-only often RAISES stress
    fig, axes = plt.subplots(1, 2, figsize=(13.2, 4.6))
    fig.set_facecolor(SURFACE)
    x = np.arange(len(labels))
    w = 0.36
    dt_h = [r["heat"]["Tmax_best_C"] - r["heat"]["Tmax_init_C"] for r in rows]
    dt_c = [r["coupled"]["Tmax_best_C"] - r["coupled"]["Tmax_init_C"] for r in rows]
    dv_h = [r["vm_h"] - r["vm_h_init"] for r in rows]
    dv_c = [r["vm_c"] - r["vm_c_init"] for r in rows]
    axes[0].bar(x - 0.2, dt_h, w, color=HEAT, label="只热", zorder=3)
    axes[0].bar(x + 0.2, dt_c, w, color=COUP, label="热+应力", zorder=3)
    axes[0].axhline(0, color=ZERO, linewidth=1)
    axes[0].set_title("G  相对初值 · ΔTmax（负=降温）", loc="left", color=TEXT, fontsize=11)
    axes[1].bar(x - 0.2, dv_h, w, color=HEAT, label="只热", zorder=3)
    axes[1].bar(x + 0.2, dv_c, w, color=COUP, label="热+应力", zorder=3)
    axes[1].axhline(0, color=ZERO, linewidth=1)
    axes[1].set_title("H  相对初值 · Δvm（负=应力降）", loc="left", color=TEXT, fontsize=11)
    for ax in axes:
        ax.set_xticks(x)
        ax.set_xticklabels(labels, color=TEXT, fontsize=8)
        _style_ax(ax)
        ax.legend(frameon=False, fontsize=8)
    axes[0].set_ylabel("°C", color=MUTED)
    axes[1].set_ylabel("MPa", color=MUTED)
    fig.tight_layout()
    fig.savefig(OUT / "compare_vs_init.png", dpi=150, bbox_inches="tight", facecolor=SURFACE)
    plt.close(fig)

    # 4) hero scatter: ΔTmax vs Δvm  — one point per scheme
    fig, ax = plt.subplots(figsize=(8.6, 6.2))
    fig.set_facecolor(SURFACE)
    ax.set_facecolor(SURFACE)
    ax.axhline(0, color=ZERO, linewidth=1, zorder=1)
    ax.axvline(0, color=ZERO, linewidth=1, zorder=1)
    ax.fill_between([-0.02, 0.22], -32, 0, color="#1baf7a", alpha=0.08, zorder=0)
    ax.text(0.14, -30.4, "用一点温度换应力下降", color="#008300", fontsize=9, ha="center")
    dts = [r["delta"]["dTmax_C"] for r in rows]
    dvms = [r["dvm"] for r in rows]
    ax.scatter(dts, dvms, s=78, color=COUP, zorder=4, edgecolors=SURFACE, linewidths=1.2)
    # data-space offsets with a short leader so near-origin / 3+6 overlap stay readable
    OFF = {
        "3芯": (0.018, 1.8),
        "4芯": (0.012, -2.2),
        "5芯": (-0.018, -2.4),
        "6芯": (0.018, -1.6),
        "7芯": (-0.018, -1.8),
        "10芯": (0.016, 1.8),
        "11芯": (-0.012, 2.2),
        "12芯": (0.012, -2.4),
        "13芯": (-0.018, 1.8),
        "14芯": (0.016, 1.8),
        "15芯": (0.012, 1.8),
        "16芯": (0.018, -2.0),
    }
    for r, dt, dvm in zip(rows, dts, dvms):
        ox, oy = OFF[r["label"]]
        ax.annotate(
            r["label"],
            xy=(dt, dvm),
            xytext=(dt + ox, dvm + oy),
            textcoords="data",
            color=TEXT,
            fontsize=9,
            arrowprops=dict(arrowstyle="-", color="#b4b2a9", lw=0.8),
            zorder=5,
        )
    ax.set_xlabel("ΔTmax  （热+应力 − 只热，°C）  正=更热", color=MUTED, fontsize=10)
    ax.set_ylabel("Δvm  （热+应力 − 只热，MPa）  负=应力更低", color=MUTED, fontsize=10)
    ax.set_title("应力进 SA 成本之后，布局往哪偏", loc="left", color=TEXT, fontsize=13, pad=10)
    ax.xaxis.grid(True, color=GRID, linewidth=1, zorder=0)
    ax.yaxis.grid(True, color=GRID, linewidth=1, zorder=0)
    ax.set_axisbelow(True)
    ax.tick_params(colors=MUTED, length=0)
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.spines["bottom"].set_visible(True)
    ax.spines["left"].set_visible(True)
    ax.spines["bottom"].set_color(GRID)
    ax.spines["left"].set_color(GRID)
    ax.set_xlim(-0.02, 0.16)
    fig.tight_layout()
    fig.savefig(OUT / "compare_tradeoff_scatter.png", dpi=150, bbox_inches="tight", facecolor=SURFACE)
    plt.close(fig)

    payload = []
    for r in rows:
        payload.append(
            {
                "id": r["id"],
                "short": r["short"],
                "n_chiplets": r["n_chiplets"],
                "n_steps": r["n_steps"],
                "vm_metric": r["metric"],
                "heat_only_init_Tmax_C": r["heat"]["Tmax_init_C"],
                "coupled_init_Tmax_C": r["coupled"]["Tmax_init_C"],
                "heat_only_best_Tmax_C": r["heat"]["Tmax_best_C"],
                "coupled_best_Tmax_C": r["coupled"]["Tmax_best_C"],
                "heat_only_init_vm_MPa": r["vm_h_init"],
                "coupled_init_vm_MPa": r["vm_c_init"],
                "heat_only_best_vm_MPa": r["vm_h"],
                "coupled_best_vm_MPa": r["vm_c"],
                "heat_only_vm_chip": r["vm_name_h"],
                "coupled_vm_chip": r["vm_name_c"],
                "heat_only_best_wpv_um": r["heat"]["w_pv_best_um"],
                "coupled_best_wpv_um": r["coupled"]["w_pv_best_um"],
                "dTmax_C": r["delta"]["dTmax_C"],
                "d_vm_MPa": r["dvm"],
                "d_wpv_um": r["delta"]["d_w_pv_um"],
                "heat_dT_vs_init": r["heat"]["Tmax_best_C"] - r["heat"]["Tmax_init_C"],
                "coup_dT_vs_init": r["coupled"]["Tmax_best_C"] - r["coupled"]["Tmax_init_C"],
                "heat_dvm_vs_init": r["vm_h"] - r["vm_h_init"],
                "coup_dvm_vs_init": r["vm_c"] - r["vm_c_init"],
                "heat_accept": r["heat"]["accepted"],
                "coupled_accept": r["coupled"]["accepted"],
            }
        )
    (OUT / "all_schemes.json").write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print("wrote", OUT)


if __name__ == "__main__":
    main()
