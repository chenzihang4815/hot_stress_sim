"""布局图和算法对比图。"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

from geometry import Rect

SURFACE = "#fcfcfb"
TEXT = "#0b0b0b"
TEXT2 = "#52514e"
GRID = "#e8e7e3"

ALG_COLOR = {
    "bstree": "#2a78d6",
    "mcts": "#eb6834",
    "lbfgs": "#1baf7a",
}
ALG_ORDER = ["bstree", "mcts", "lbfgs"]

KIND_COLOR = {
    "CPU": "#2a78d6",
    "GPU": "#eb6834",
    "CIM": "#1baf7a",
    "ACCEL": "#eda100",
}


def _style(ax) -> None:
    ax.set_facecolor(SURFACE)
    ax.figure.set_facecolor(SURFACE)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color("#d9d8d3")
    ax.spines["bottom"].set_color("#d9d8d3")
    ax.tick_params(colors=TEXT2, labelsize=8)
    ax.yaxis.grid(True, color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)


def draw_floorplan(rects: list[Rect], path: Path, title: str) -> None:
    fig, ax = plt.subplots(figsize=(7.2, 5.6))
    ax.set_facecolor(SURFACE)
    fig.set_facecolor(SURFACE)
    for r in rects:
        color = KIND_COLOR.get(r.kind, "#9D755D")
        ax.add_patch(
            Rectangle(
                (r.x * 1e3, r.y * 1e3),
                r.w * 1e3,
                r.h * 1e3,
                facecolor=color,
                edgecolor=SURFACE,
                linewidth=1.2,
                alpha=0.92,
            )
        )
        ax.add_patch(
            Rectangle(
                (r.x * 1e3, r.y * 1e3),
                r.w * 1e3,
                r.h * 1e3,
                facecolor="none",
                edgecolor=TEXT,
                linewidth=0.6,
            )
        )
        ax.text(
            (r.x + r.w / 2) * 1e3,
            (r.y + r.h / 2) * 1e3,
            f"{r.name}\n{r.power_W:.3g} W\n{r.area_mm2:.3g} mm²",
            ha="center",
            va="center",
            fontsize=7,
            color=TEXT,
        )
    ax.set_xlabel("x (mm)", color=TEXT2, fontsize=9)
    ax.set_ylabel("y (mm)", color=TEXT2, fontsize=9)
    ax.set_title(title, color=TEXT, fontsize=11, loc="left")
    ax.set_aspect("equal", adjustable="box")
    xmax = max(r.x + r.w for r in rects) * 1e3
    ymax = max(r.y + r.h for r in rects) * 1e3
    ax.set_xlim(-0.02 * xmax, 1.08 * xmax)
    ax.set_ylim(-0.02 * ymax, 1.08 * ymax)
    kinds = []
    for r in rects:
        if r.kind not in kinds:
            kinds.append(r.kind)
    handles = [
        Rectangle((0, 0), 1, 1, facecolor=KIND_COLOR.get(k, "#999"), edgecolor=TEXT, linewidth=0.4)
        for k in kinds
    ]
    ax.legend(handles, kinds, frameon=False, fontsize=8, loc="upper right")
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=140)
    plt.close(fig)


def _mix_order(rows: list[dict]) -> list[str]:
    seen = []
    for r in rows:
        if r["mix"] not in seen:
            seen.append(r["mix"])
    return seen


def draw_tmax_grouped(rows: list[dict], titles: dict[str, str], path: Path) -> None:
    draw_compare_bars(rows, titles, "Tmax_C", "Tmax (C)", path, baseline=45.0)


def draw_compare_bars(
    rows: list[dict],
    titles: dict[str, str],
    field: str,
    ylabel: str,
    path: Path,
    baseline: float | None = None,
) -> None:
    mixes = _mix_order(rows)
    by = {(r["mix"], r["algorithm"]): r[field] for r in rows}
    n = len(mixes)
    width = 0.24
    x = list(range(n))
    fig, ax = plt.subplots(figsize=(11.2, 5.4))
    _style(ax)
    for i, alg in enumerate(ALG_ORDER):
        xs = [xi + (i - 1) * (width + 0.02) for xi in x]
        ys = [by[(m, alg)] for m in mixes]
        ax.bar(
            xs,
            ys,
            width=width,
            color=ALG_COLOR[alg],
            label=alg,
            linewidth=0,
            zorder=2,
        )
    ax.set_xticks(x)
    labels = [titles.get(m, m) for m in mixes]
    ax.set_xticklabels(labels, rotation=22, ha="right", fontsize=8, color=TEXT)
    ax.set_ylabel(ylabel, color=TEXT2, fontsize=9)
    ax.set_title(ylabel + "  ·  bstree / mcts / lbfgs", color=TEXT, fontsize=12, loc="left")
    ax.legend(frameon=False, fontsize=8, title=None, loc="upper right")
    if field in {"Tmax_C", "Tmin_C", "dT_C"}:
        lo = min(r[field] for r in rows)
        hi = max(r[field] for r in rows)
        pad = max(0.4, 0.08 * (hi - lo + 1e-9))
        ax.set_ylim(lo - pad, hi + pad)
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=140)
    plt.close(fig)
