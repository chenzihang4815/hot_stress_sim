#!/usr/bin/env python3
"""Resample ATSim die-midplane T onto HotSpot Layer-0 grid.steady and draw official maps.

HotSpot grid convention (temperature_grid.c):
  row 0 = north (high y), last row = south (low y)
  dumped i*cols+j, i from 0..rows-1
Official grid_thermal_map.py uses imshow origin='upper', so row 0 is drawn at the top.
Writing south-first therefore appears as a y-axis flip (top-view vs bottom-view).
"""

from __future__ import annotations

import json
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
from scipy.interpolate import RegularGridInterpolator, griddata

ROOT = Path("/home/czh/任务仿真产生功耗文件工具")
ATSIM_RES = ROOT / "热应力python工程_最终" / "results"
SA = ROOT / "热仿真流水线" / "runs" / "sa_gapped_13chip"
PL = ROOT / "HotSpot-master" / "scripts" / "grid_thermal_map.pl"
PYMAP = ROOT / "HotSpot-master" / "scripts" / "grid_thermal_map.py"
BAK = ROOT / "成功运行备份" / "SA带布线间隙_13芯_2026-10-22"

ROWS = COLS = 64
TAMB = 318.15  # 45 C
DIE_Z = 0.001175


def parse_flp(path: Path):
    chips = []
    for raw in path.read_text(encoding="utf-8", errors="replace").splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) < 5:
            continue
        name, w, h, x, y = parts[0], float(parts[1]), float(parts[2]), float(parts[3]), float(parts[4])
        chips.append((name, x, y, w, h))
    minx = min(c[1] for c in chips)
    miny = min(c[2] for c in chips)
    maxx = max(c[1] + c[3] for c in chips)
    maxy = max(c[2] + c[4] for c in chips)
    return chips, minx, miny, maxx, maxy


def atsim_shift(summary: dict, flp_chips):
    by_flp = {n: (x, y) for n, x, y, w, h in flp_chips}
    shifts = []
    for c in summary["chips"]:
        if c["name"] in by_flp:
            sx = c["x_mm"] * 1e-3 - by_flp[c["name"]][0]
            sy = c["y_mm"] * 1e-3 - by_flp[c["name"]][1]
            shifts.append((sx, sy, c["name"]))
    if not shifts:
        raise SystemExit("cannot match chip names for ATSim↔HotSpot shift")
    sx0, sy0, _ = shifts[0]
    for sx, sy, name in shifts:
        if abs(sx - sx0) > 1e-9 or abs(sy - sy0) > 1e-9:
            raise SystemExit(f"inconsistent chip shift at {name}: {(sx, sy)} vs {(sx0, sy0)}")
    return sx0, sy0


def enlarge(svg: Path) -> None:
    t = svg.read_text(encoding="utf-8")
    t = t.replace('style="font-size:180"', 'style="font-size:420"')
    t = t.replace('style="font-size:250"', 'style="font-size:700"')
    svg.write_text(t, encoding="utf-8")


def cell_edges(minx, miny, maxx, maxy):
    """HotSpot cell edges: row 0 = north (high y), col 0 = west (low x)."""
    W, H = maxx - minx, maxy - miny
    x_e = minx + np.linspace(0.0, 1.0, COLS + 1) * W
    # row 0 north: y decreases as row index increases
    y_e = maxy - np.linspace(0.0, 1.0, ROWS + 1) * H
    return x_e, y_e


def overlap_rect(ax0, ay0, ax1, ay1, bx0, by0, bx1, by1) -> float:
    dx = min(ax1, bx1) - max(ax0, bx0)
    dy = min(ay1, by1) - max(ay0, by0)
    if dx <= 0.0 or dy <= 0.0:
        return 0.0
    return dx * dy


def sample_point_in_overlap(cx, cy, x0, y0, w, h, x_lo, y_lo, x_hi, y_hi):
    """Prefer the cell center if it sits in the chip; else the overlap centroid."""
    if x0 <= cx <= x0 + w and y0 <= cy <= y0 + h:
        return cx, cy
    ox0, oy0 = max(x0, x_lo), max(y0, y_lo)
    ox1, oy1 = min(x0 + w, x_hi), min(y0 + h, y_hi)
    return 0.5 * (ox0 + ox1), 0.5 * (oy0 + oy1)


def sample(
    npz_path: Path,
    flp_path: Path,
    summary: dict,
    out_grid: Path,
    airgap: bool,
) -> tuple[float, float]:
    data = np.load(npz_path)
    T, xc, yc, zc = data["T"], data["xc"], data["yc"], data["zc"]
    k = int(np.argmin(np.abs(zc - DIE_Z)))
    sl = T[:, :, k]
    chips, minx, miny, maxx, maxy = parse_flp(flp_path)
    sx, sy = atsim_shift(summary, chips)
    x_e, y_e = cell_edges(minx, miny, maxx, maxy)

    interp = RegularGridInterpolator((xc, yc), sl, bounds_error=False, fill_value=np.nan)
    solid = sl > (TAMB + 0.5)
    ii, jj = np.where(solid)
    pts = np.column_stack([xc[ii], yc[jj]])
    ts = sl[ii, jj]

    vals = []
    for i in range(ROWS):
        y_hi, y_lo = float(y_e[i]), float(y_e[i + 1])  # north edge, south edge
        cy = 0.5 * (y_lo + y_hi)
        for j in range(COLS):
            x_lo, x_hi = float(x_e[j]), float(x_e[j + 1])
            cx = 0.5 * (x_lo + x_hi)
            best_area = 0.0
            sample_xy = None
            for _n, x0, y0, w, h in chips:
                area = overlap_rect(x_lo, y_lo, x_hi, y_hi, x0, y0, x0 + w, y0 + h)
                if area > best_area:
                    best_area = area
                    sample_xy = sample_point_in_overlap(cx, cy, x0, y0, w, h, x_lo, y_lo, x_hi, y_hi)
            if best_area <= 0.0:
                if airgap:
                    vals.append(TAMB)
                    continue
                x, y = cx + sx, cy + sy
                t = float(interp((x, y)))
                if not np.isfinite(t) or t < TAMB + 0.5:
                    t = float(griddata(pts, ts, [(x, y)], method="nearest")[0]) if pts.size else TAMB
                vals.append(t)
                continue
            px, py = sample_xy
            x, y = px + sx, py + sy
            t = float(interp((x, y)))
            if not np.isfinite(t) or t < TAMB + 0.5:
                t = float(griddata(pts, ts, [(x, y)], method="nearest")[0]) if pts.size else TAMB
            vals.append(t)
    out_grid.write_text("\n".join(f"{n}\t{vals[n]:.2f}" for n in range(len(vals))) + "\n", encoding="utf-8")
    return min(vals), max(vals)


def draw(flp: Path, grid: Path, png: Path, svg: Path, lo=None, hi=None) -> None:
    subprocess.run(
        [sys.executable, str(PYMAP), str(flp), str(grid), "64", "64", str(png)],
        check=True,
    )
    cmd = ["perl", str(PL), str(flp), str(grid), "64", "64"]
    if lo is not None:
        cmd += [f"{lo:.4f}", f"{hi:.4f}"]
    svg.write_text(subprocess.check_output(cmd, text=True), encoding="utf-8")
    enlarge(svg)


def main() -> None:
    jobs = [
        ("init", SA / "layout_init" / "floorplan.flp", ATSIM_RES / "python_sa13_gapped_init"),
        ("best", SA / "layout_best" / "floorplan.flp", ATSIM_RES / "python_sa13_gapped_best"),
    ]
    fill_jobs = [
        ("init", SA / "layout_init" / "floorplan.flp", ATSIM_RES / "python_sa13_gapped_init_sifill"),
        ("best", SA / "layout_best" / "floorplan.flp", ATSIM_RES / "python_sa13_gapped_best_sifill"),
    ]
    SA.mkdir(parents=True, exist_ok=True)

    solid_mm, air_mm, fill_mm = {}, {}, {}
    for tag, flp, resdir in jobs:
        summary = json.loads((resdir / "summary.json").read_text(encoding="utf-8"))
        npz = resdir / "T_field.npz"

        grid = SA / f"atsim_{tag}_layer0.grid.steady"
        tmin, tmax = sample(npz, flp, summary, grid, airgap=False)
        solid_mm[tag] = (tmin, tmax)
        print(tag, "solid K", tmin, tmax, "C", tmin - 273.15, tmax - 273.15)
        draw(flp, grid, SA / f"atsim_{tag}_hotspot.png", SA / f"atsim_{tag}_hotspot.svg")

        grid_a = SA / f"atsim_{tag}_airgap_layer0.grid.steady"
        tmin, tmax = sample(npz, flp, summary, grid_a, airgap=True)
        air_mm[tag] = (tmin, tmax)
        print(tag, "airgap K", tmin, tmax, "C", tmin - 273.15, tmax - 273.15)
        draw(
            flp,
            grid_a,
            SA / f"atsim_{tag}_airgap_hotspot.png",
            SA / f"atsim_{tag}_airgap_hotspot.svg",
        )

    for tag, flp, resdir in fill_jobs:
        summary = json.loads((resdir / "summary.json").read_text(encoding="utf-8"))
        npz = resdir / "T_field.npz"
        grid = SA / f"atsim_{tag}_sifill_layer0.grid.steady"
        tmin, tmax = sample(npz, flp, summary, grid, airgap=False)
        fill_mm[tag] = (tmin, tmax)
        print(tag, "sifill K", tmin, tmax, "C", tmin - 273.15, tmax - 273.15)
        draw(
            flp,
            grid,
            SA / f"atsim_{tag}_sifill_hotspot.png",
            SA / f"atsim_{tag}_sifill_hotspot.svg",
        )

    lo_s = min(v[0] for v in solid_mm.values())
    hi_s = max(v[1] for v in solid_mm.values())
    lo_a = min(v[0] for v in air_mm.values())
    hi_a = max(v[1] for v in air_mm.values())
    lo_f = min(v[0] for v in fill_mm.values())
    hi_f = max(v[1] for v in fill_mm.values())
    print("shared solid", lo_s, hi_s, "air", lo_a, hi_a, "sifill", lo_f, hi_f)
    for tag, flp, _resdir in jobs:
        draw(
            flp,
            SA / f"atsim_{tag}_layer0.grid.steady",
            SA / f"atsim_{tag}_hotspot.png",
            SA / f"atsim_{tag}_hotspot_scaled.svg",
            lo_s,
            hi_s,
        )
        draw(
            flp,
            SA / f"atsim_{tag}_airgap_layer0.grid.steady",
            SA / f"atsim_{tag}_airgap_hotspot.png",
            SA / f"atsim_{tag}_airgap_hotspot_scaled.svg",
            lo_a,
            hi_a,
        )
    for tag, flp, _resdir in fill_jobs:
        draw(
            flp,
            SA / f"atsim_{tag}_sifill_layer0.grid.steady",
            SA / f"atsim_{tag}_sifill_hotspot.png",
            SA / f"atsim_{tag}_sifill_hotspot_scaled.svg",
            lo_f,
            hi_f,
        )

    BAK.mkdir(parents=True, exist_ok=True)
    for p in SA.glob("atsim_*"):
        shutil.copy2(p, BAK / p.name)
    print("copied to", BAK)


if __name__ == "__main__":
    main()
