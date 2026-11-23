#!/usr/bin/env python3
"""Draw HotSpot Layer-0 grid heatmaps for the 16 CPU/GPU/CIM L-BFGS cases."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Rectangle

ROOT = Path("/home/czh/任务仿真产生功耗文件工具")
ATSIM = ROOT / "热应力python工程"
CASES = ATSIM / "cases" / "cpu_gpu_cim_lbfgs"
OUT = ROOT / "热仿真流水线" / "runs" / "hotspot_heatmaps"
BACKUP = ROOT / "成功运行备份" / "三芯片CPU_GPU_CIM_HotSpot热图_2026-10-21"

# sequential blue, light (cool) -> dark (hot)
BLUE_RAMP = [
    "#cde2fb",
    "#9ec5f4",
    "#6da7ec",
    "#3987e5",
    "#2a78d6",
    "#256abf",
    "#1c5cab",
    "#104281",
    "#0d366b",
]
CMAP = LinearSegmentedColormap.from_list("hotspot_T", BLUE_RAMP, N=256)
CHIP_EDGE = {"CPU": "#0b0b0b", "GPU": "#0b0b0b", "CIM": "#0b0b0b"}


def parse_flp(path: Path) -> tuple[list[dict], float, float]:
    chips = []
    for raw in path.read_text(encoding="utf-8", errors="replace").splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) < 5:
            continue
        name, w, h, x, y = parts[0], float(parts[1]), float(parts[2]), float(parts[3]), float(parts[4])
        chips.append({"name": name, "w": w, "h": h, "x": x, "y": y})
    width = max(c["x"] + c["w"] for c in chips)
    height = max(c["y"] + c["h"] for c in chips)
    return chips, width, height


def parse_layer0(path: Path, rows: int = 64, cols: int = 64) -> np.ndarray:
    vals = []
    in_layer = False
    for raw in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if raw.startswith("Layer 0"):
            in_layer = True
            continue
        if in_layer and raw.startswith("Layer"):
            break
        if not in_layer:
            continue
        parts = raw.split()
        if len(parts) >= 2:
            vals.append(float(parts[1]))
    arr = np.array(vals, dtype=float).reshape(rows, cols)
    return arr - 273.15  # C; row 0 = south/bottom if HotSpot y increases with row


def plot_one(mix_id: str, grid_C: np.ndarray, chips: list[dict], width_m: float, height_m: float, path: Path, title: str) -> None:
    fig, ax = plt.subplots(figsize=(7.0, 5.8))
    fig.set_facecolor("#fcfcfb")
    ax.set_facecolor("#fcfcfb")
    extent = [0.0, width_m * 1e3, 0.0, height_m * 1e3]
    im = ax.imshow(
        grid_C,
        origin="lower",
        extent=extent,
        cmap=CMAP,
        aspect="equal",
        interpolation="bilinear",
    )
    cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cb.set_label("T (C)", color="#52514e")
    cb.ax.tick_params(colors="#52514e", labelsize=8)
    for c in chips:
        if c["name"] not in ("CPU", "GPU", "CIM"):
            continue
        ax.add_patch(
            Rectangle(
                (c["x"] * 1e3, c["y"] * 1e3),
                c["w"] * 1e3,
                c["h"] * 1e3,
                fill=False,
                edgecolor="#0b0b0b",
                linewidth=1.0,
            )
        )
        ax.text(
            (c["x"] + 0.5 * c["w"]) * 1e3,
            (c["y"] + 0.5 * c["h"]) * 1e3,
            c["name"],
            ha="center",
            va="center",
            fontsize=8,
            color="#0b0b0b",
        )
    ax.set_xlabel("x (mm)", color="#52514e")
    ax.set_ylabel("y (mm)", color="#52514e")
    ax.set_title(title, loc="left", color="#0b0b0b", fontsize=11)
    ax.tick_params(colors="#52514e", labelsize=8)
    for sp in ax.spines.values():
        sp.set_color("#d9d8d3")
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=140)
    plt.close(fig)


def plot_montage(items: list[tuple[str, np.ndarray, list, float, float]], path: Path, vmin: float, vmax: float) -> None:
    fig, axes = plt.subplots(4, 4, figsize=(14.5, 13.2))
    fig.set_facecolor("#fcfcfb")
    for ax, (mix, grid, chips, w, h) in zip(axes.ravel(), items):
        ax.set_facecolor("#fcfcfb")
        im = ax.imshow(
            grid,
            origin="lower",
            extent=[0, w * 1e3, 0, h * 1e3],
            cmap=CMAP,
            aspect="equal",
            interpolation="bilinear",
            vmin=vmin,
            vmax=vmax,
        )
        for c in chips:
            if c["name"] not in ("CPU", "GPU", "CIM"):
                continue
            ax.add_patch(
                Rectangle(
                    (c["x"] * 1e3, c["y"] * 1e3),
                    c["w"] * 1e3,
                    c["h"] * 1e3,
                    fill=False,
                    edgecolor="#0b0b0b",
                    linewidth=0.6,
                )
            )
        ax.set_title(mix, fontsize=8, color="#0b0b0b")
        ax.set_xticks([])
        ax.set_yticks([])
        for sp in ax.spines.values():
            sp.set_color("#d9d8d3")
    cbar = fig.colorbar(im, ax=axes.ravel().tolist(), fraction=0.02, pad=0.02)
    cbar.set_label("T (C)")
    fig.suptitle("HotSpot Layer 0  ·  16 CPU+GPU+CIM L-BFGS cases", fontsize=13, color="#0b0b0b")
    fig.savefig(path, dpi=140, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    json_files = sorted(p for p in CASES.glob("*.json") if p.name != "index.json")
    items = []
    stats = []
    for jp in json_files:
        spec = json.loads(jp.read_text(encoding="utf-8"))
        mix = spec.get("mix") or jp.stem
        case_name = spec["name"]
        hs_dir = ATSIM / f"hotspot_{jp.stem}"
        flp = hs_dir / f"{case_name}.flp"
        gridp = hs_dir / "outputs" / f"{case_name}.grid.steady"
        if not flp.exists() or not gridp.exists():
            raise SystemExit(f"missing HotSpot files for {mix}: {flp} {gridp}")
        chips, width, height = parse_flp(flp)
        grid = parse_layer0(gridp)
        title = f"{mix}  T {grid.min():.2f}–{grid.max():.2f} C"
        plot_one(mix, grid, chips, width, height, OUT / f"{mix}.png", title)
        items.append((mix, grid, chips, width, height))
        stats.append(
            {
                "mix": mix,
                "Tmin_C": float(grid.min()),
                "Tmax_C": float(grid.max()),
                "Tmean_C": float(grid.mean()),
                "width_mm": width * 1e3,
                "height_mm": height * 1e3,
            }
        )
        print(mix, f"T={grid.min():.2f}–{grid.max():.2f} C")
    vmin = min(s["Tmin_C"] for s in stats)
    vmax = max(s["Tmax_C"] for s in stats)
    plot_montage(items, OUT / "montage_16.png", vmin, vmax)
    (OUT / "heatmap_stats.json").write_text(json.dumps(stats, indent=2) + "\n", encoding="utf-8")
    import shutil

    if BACKUP.exists():
        shutil.rmtree(BACKUP)
    shutil.copytree(OUT, BACKUP)
    print("wrote", OUT / "montage_16.png")


if __name__ == "__main__":
    main()
