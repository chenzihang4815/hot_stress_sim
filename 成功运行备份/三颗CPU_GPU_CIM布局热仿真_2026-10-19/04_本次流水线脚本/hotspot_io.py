"""把矩形布局写成 HotSpot 的 .flp / .ptrace / .config 并调用仿真器。"""

from __future__ import annotations

import json
import subprocess
from pathlib import Path

from geometry import Rect, packing_metrics

DEFAULT_CONFIG = """\
-t_chip                 0.00015
-k_chip                 130.0
-p_chip                 1630300
-c_convec               140.4
-r_convec               0.1
-s_sink                 {s_sink:.12g}
-t_sink                 0.0069
-k_sink                 400.0
-p_sink                 3.55e6
-s_spreader             {s_spreader:.12g}
-t_spreader             0.001
-k_spreader             400.0
-p_spreader             3.55e6
-t_interface            2.0e-05
-k_interface            4.0
-p_interface            4.0e6
-ambient                318.15
-init_file              (null)
-init_temp              318.15
-steady_file            (null)
-sampling_intvl         0.01
-base_proc_freq         3e9
-dtm_used               0
-model_type             grid
-leakage_used           0
-leakage_mode           0
-package_model_used     0
-package_config_file    package.config
-block_omit_lateral     0
-grid_rows              {rows}
-grid_cols              {cols}
-grid_layer_file        (null)
-grid_steady_file       (null)
-grid_map_mode          avg
-model_secondary        0
-use_microfluidic_cooling 0
"""


def write_hotspot_inputs(rects: list[Rect], out_dir: Path, rows: int = 64, cols: int = 64) -> dict:
    out_dir.mkdir(parents=True, exist_ok=True)
    metrics = packing_metrics(rects)
    flp = out_dir / "floorplan.flp"
    lines = [
        "# unit-name  width  height  left-x  bottom-y   (meters)",
    ]
    for r in rects:
        lines.append(f"{r.name}\t{r.w:.12g}\t{r.h:.12g}\t{r.x:.12g}\t{r.y:.12g}")
    flp.write_text("\n".join(lines) + "\n", encoding="utf-8")

    ptrace = out_dir / "power.ptrace"
    names = [r.name for r in rects]
    powers = [f"{r.power_W:.12g}" for r in rects]
    ptrace.write_text("\t".join(names) + "\n" + "\t".join(powers) + "\n", encoding="utf-8")

    chip = max(metrics["width_m"], metrics["height_m"])
    cfg = out_dir / "hotspot.config"
    cfg.write_text(
        DEFAULT_CONFIG.format(
            s_sink=max(0.06, chip * 1.2),
            s_spreader=max(0.03, chip * 1.1),
            rows=rows,
            cols=cols,
        ),
        encoding="utf-8",
    )
    manifest = {
        "block_count": len(rects),
        "power_total_W": sum(r.power_W for r in rects),
        "geometry_units": "m",
        "power_units": "W",
        **metrics,
    }
    (out_dir / "layout.json").write_text(
        json.dumps(
            [
                {
                    "name": r.name,
                    "kind": r.kind,
                    "x_m": r.x,
                    "y_m": r.y,
                    "w_m": r.w,
                    "h_m": r.h,
                    "area_mm2": r.area_mm2,
                    "power_W": r.power_W,
                }
                for r in rects
            ],
            indent=2,
            ensure_ascii=False,
        )
        + "\n",
        encoding="utf-8",
    )
    (out_dir / "input_manifest.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    return manifest


def run_hotspot(hotspot_bin: Path, input_dir: Path) -> dict:
    cmd = [
        str(hotspot_bin),
        "-f",
        str(input_dir / "floorplan.flp"),
        "-p",
        str(input_dir / "power.ptrace"),
        "-c",
        str(input_dir / "hotspot.config"),
        "-model_type",
        "grid",
        "-steady_file",
        str(input_dir / "block.steady"),
        "-grid_steady_file",
        str(input_dir / "grid.steady"),
    ]
    completed = subprocess.run(cmd, cwd=input_dir, text=True, capture_output=True, check=False)
    (input_dir / "hotspot.stdout.log").write_text(completed.stdout, encoding="utf-8")
    (input_dir / "hotspot.stderr.log").write_text(completed.stderr, encoding="utf-8")
    result = {
        "command": cmd,
        "returncode": completed.returncode,
        "block_steady": str(input_dir / "block.steady"),
        "grid_steady": str(input_dir / "grid.steady"),
    }
    (input_dir / "command.json").write_text(
        json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    if completed.returncode != 0:
        raise RuntimeError(
            f"HotSpot 失败 rc={completed.returncode}\nSTDOUT:\n{completed.stdout}\nSTDERR:\n{completed.stderr}"
        )
    temps = _read_block_steady(input_dir / "block.steady")
    result["temperatures_K"] = temps
    if temps:
        result["Tmax_K"] = max(temps.values())
        result["Tmin_K"] = min(temps.values())
        result["Tmax_C"] = result["Tmax_K"] - 273.15
        result["Tmin_C"] = result["Tmin_K"] - 273.15
    (input_dir / "thermal_summary.json").write_text(
        json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    return result


def _read_block_steady(path: Path) -> dict[str, float]:
    temps: dict[str, float] = {}
    if not path.exists():
        return temps
    for raw in path.read_text(encoding="utf-8", errors="replace").splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) < 2:
            continue
        try:
            temps[parts[0]] = float(parts[1])
        except ValueError:
            continue
    return temps
