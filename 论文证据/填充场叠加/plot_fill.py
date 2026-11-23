#!/usr/bin/env python3
"""四种填充的误差场。每张图并排起点和最优布局，色标是 frozen 减 exact。"""

import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
import numpy as np

font_manager.fontManager.addfont(
    "/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc")
plt.rcParams["font.family"] = "Noto Sans CJK JP"
plt.rcParams["axes.unicode_minus"] = False

ROOT = Path("/home/czh/任务仿真产生功耗文件工具")
sys.path.insert(0, str(ROOT / "热仿真流水线" / "src"))
import run_hw_area_field as hw  # noqa: E402

HERE = Path(__file__).resolve().parent
RUNS = HERE / "runs"
FIG = HERE / "图"
FIG.mkdir(exist_ok=True)

NAMES = {"air": "空气 0.028", "epoxy": "环氧 0.8", "underfill": "底部填充 1.0",
         "silicon": "硅 130"}
SIDE = 72.0


def boxes(layout_path: Path):
    rows = json.loads(layout_path.read_text(encoding="utf-8"))
    return [(r["x_mm"], r["y_mm"], r["w_mm"], r["h_mm"]) for r in rows]


def start_boxes():
    src = (ROOT / "热仿真流水线" / "runs"
           / "sa_heat_vs_coupled_pkg_1cpu_1gpu_1cim" / "layout_start.json")
    rows = json.loads(src.read_text(encoding="utf-8"))
    return [(r["x_m"] * 1e3, r["y_m"] * 1e3, r["w_m"] * 1e3, r["h_m"] * 1e3)
            for r in rows]


def draw(ax, field, rects, title):
    im = ax.imshow(field, origin="lower", extent=(0, SIDE, 0, SIDE),
                   cmap="coolwarm", vmin=-abs(field).max(), vmax=abs(field).max())
    for x, y, w, h in rects:
        ax.add_patch(plt.Rectangle((x, y), w, h, fill=False, ec="black", lw=0.8))
    ax.set_title(title, fontsize=11)
    ax.set_xlabel("x (mm)")
    ax.set_ylabel("y (mm)")
    return im


def main():
    start = start_boxes()
    for fill in ("air", "epoxy", "underfill", "silicon"):
        g = np.load(RUNS / fill / "grids.npz")
        summary = json.loads((RUNS / fill / "summary.json").read_text(encoding="utf-8"))
        trace = json.loads((RUNS / fill / "trace.json").read_text(encoding="utf-8"))
        at_best = next(r for r in trace if r["step"] == summary["best_step"])
        best = boxes(RUNS / fill / "best_layout.json")
        err = g["best_err"]
        fig, axes = plt.subplots(1, 2, figsize=(11, 5.2))
        span = max(abs(err.min()), abs(err.max()), 1e-12)
        for ax, rects, label in (
                (axes[0], start, "起点材料图"),
                (axes[1], best, f"最优布局，第 {summary['best_step']} 步")):
            im = ax.imshow(err, origin="lower", extent=(0, SIDE, 0, SIDE),
                           cmap="coolwarm", vmin=-span, vmax=span)
            for x, y, w, h in rects:
                ax.add_patch(plt.Rectangle((x, y), w, h, fill=False,
                                           ec="black", lw=0.8))
            ax.set_title(label, fontsize=11)
            ax.set_xlabel("x (mm)")
            ax.set_ylabel("y (mm)")
        fig.colorbar(im, ax=axes, label="frozen − exact (K)", fraction=0.046)
        fig.suptitle(
            f"{NAMES[fill]}，不更新材料图，最优布局上的误差，"
            f"峰值差 {at_best['peak_err_K']:+.2f} K，格最大 {at_best['grid_max_K']:.2f} K",
            fontsize=12)
        fig.savefig(FIG / f"{fill}_误差.png", dpi=140, bbox_inches="tight")
        plt.close(fig)
        print(f"{fill} 写到 {FIG / (fill + '_误差.png')}", flush=True)


if __name__ == "__main__":
    main()
