#!/usr/bin/env python3
"""L-BFGS 三颗布局 → HotSpot 与 ATSim 热对比（只比温度）。"""

from __future__ import annotations

import csv
import json
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np

ROOT = Path("/home/czh/任务仿真产生功耗文件工具")
PIPE = ROOT / "热仿真流水线"
SRC = PIPE / "src"
ATSIM = ROOT / "热应力python工程"
HOTSPOT_BIN = ROOT / "HotSpot-master" / "hotspot"
sys.path.insert(0, str(SRC))

from combos import all_mixes  # noqa: E402
from defaults import INITIAL_LAYOUT, INITIAL_LAYOUT_SEED  # noqa: E402
from lbfgs import random_lbfgs_layout  # noqa: E402

OUT = PIPE / "runs" / "hotspot_vs_atsim_16"
BACKUP = ROOT / "成功运行备份" / "三芯片CPU_GPU_CIM_HotSpot_vs_ATSim_2026-10-20"
CASES = ATSIM / "cases" / "cpu_gpu_cim_lbfgs"

PYTHON = "/home/czh/miniconda3/envs/fenicsx/bin/python"


def layout_to_chips(rects) -> list[dict]:
    chips = []
    for r in rects:
        chips.append(
            {
                "name": r.name,
                "x_mm": r.x * 1e3,
                "y_mm": r.y * 1e3,
                "w_mm": r.w * 1e3,
                "h_mm": r.h * 1e3,
                "t_um": 150.0,
                "power_W": r.power_W,
            }
        )
    return chips


def write_atsim_json(mix_id: str, chips: list[dict], power: float) -> Path:
    CASES.mkdir(parents=True, exist_ok=True)
    payload = {
        "name": f"lbfgs_{mix_id}",
        "description": "L-BFGS initial floorplan: CPU + GPU + CIM. Heat-only compare vs HotSpot.",
        "chips": chips,
        "package": {
            "tim_um": 20,
            "t_spreader_mm": 1.0,
            "t_sink_mm": 6.9,
            "t_substrate_mm": 1.0,
            "t_interposer_um": 100.0,
            "r_convec_K_per_W": 0.1,
            "ambient_C": 45,
            "pitch_die_mm": 1.0,
            "pitch_coarse_mm": 4.0,
            "fill_die": False,
        },
        "total_power_W": power,
        "layout_algorithm": INITIAL_LAYOUT,
        "layout_seed": INITIAL_LAYOUT_SEED,
        "mix": mix_id,
    }
    path = CASES / f"lbfgs_{mix_id}.json"
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return path


def run_atsim(json_path: Path, case_name: str, log_stem: str) -> Path:
    cmd = [PYTHON, str(ATSIM / "scripts" / "run_case3.py"), str(json_path)]
    completed = subprocess.run(cmd, cwd=str(ATSIM), text=True, capture_output=True, check=False)
    log_dir = OUT / log_stem
    log_dir.mkdir(parents=True, exist_ok=True)
    (log_dir / "atsim.stdout.log").write_text(completed.stdout, encoding="utf-8")
    (log_dir / "atsim.stderr.log").write_text(completed.stderr, encoding="utf-8")
    if completed.returncode != 0:
        raise RuntimeError(
            f"ATSim failed {json_path}\n{completed.stdout[-2000:]}\n{completed.stderr[-2000:]}"
        )
    return ATSIM / "results" / f"python_{log_stem}"


def run_hotspot_from_atsim(case_name: str, log_stem: str) -> Path:
    hs_dir = ATSIM / f"hotspot_{log_stem}"
    out_dir = hs_dir / "outputs"
    out_dir.mkdir(exist_ok=True)
    cmd = [
        str(HOTSPOT_BIN),
        "-c",
        str(hs_dir / f"{case_name}.config"),
        "-f",
        str(hs_dir / f"{case_name}.flp"),
        "-p",
        str(hs_dir / f"{case_name}.ptrace"),
        "-materials_file",
        str(hs_dir / f"{case_name}.materials"),
        "-model_type",
        "grid",
        "-grid_rows",
        "64",
        "-grid_cols",
        "64",
        "-steady_file",
        str(out_dir / f"{case_name}.steady"),
        "-grid_steady_file",
        str(out_dir / f"{case_name}.grid.steady"),
    ]
    completed = subprocess.run(cmd, cwd=str(hs_dir), text=True, capture_output=True, check=False)
    (OUT / log_stem / "hotspot.stdout.log").write_text(completed.stdout, encoding="utf-8")
    (OUT / log_stem / "hotspot.stderr.log").write_text(completed.stderr, encoding="utf-8")
    if completed.returncode != 0:
        raise RuntimeError(f"HotSpot failed {case_name}\n{completed.stderr}\n{completed.stdout}")
    return out_dir / f"{case_name}.steady"


def parse_steady(path: Path) -> dict[str, float]:
    temps = {}
    for raw in path.read_text(encoding="utf-8", errors="replace").splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) < 2:
            continue
        try:
            temps[parts[0]] = float(parts[1]) - 273.15
        except ValueError:
            continue
    return temps


def chip_stats_from_npz(npz_path: Path, chips: list[dict], t_um: float = 150.0) -> dict:
    data = np.load(npz_path)
    T = data["T"]
    xc, yc, zc = data["xc"], data["yc"], data["zc"]
    z_die0 = 1.0e-3 + 100.0e-6
    z_die1 = z_die0 + t_um * 1e-6
    kz = np.where((zc >= z_die0 - 1e-9) & (zc <= z_die1 + 1e-9))[0]
    out = {}
    for c in chips:
        x0, y0, w, h = c["x_mm"] * 1e-3, c["y_mm"] * 1e-3, c["w_mm"] * 1e-3, c["h_mm"] * 1e-3
        ix = np.where((xc >= x0 - 1e-9) & (xc <= x0 + w + 1e-9))[0]
        jy = np.where((yc >= y0 - 1e-9) & (yc <= y0 + h + 1e-9))[0]
        if ix.size == 0 or jy.size == 0 or kz.size == 0:
            out[c["name"]] = None
            continue
        sub = T[np.ix_(ix, jy, kz)] - 273.15
        out[c["name"]] = {
            "mean_C": float(np.mean(sub)),
            "max_C": float(np.max(sub)),
            "min_C": float(np.min(sub)),
            "ncell": int(sub.size),
        }
    return out


def rel_err(a: float, b: float) -> float:
    if b == 0:
        return float("nan")
    return (a - b) / abs(b) * 100.0


def main() -> None:
    if INITIAL_LAYOUT != "lbfgs":
        raise SystemExit("默认初始布局必须是 lbfgs")
    if not Path(PYTHON).exists():
        raise SystemExit(f"需要带 scipy 的 Python: {PYTHON}")
    if OUT.exists():
        shutil.rmtree(OUT)
    OUT.mkdir(parents=True)

    json_files = sorted(CASES.glob("*.json"))
    json_files = [p for p in json_files if p.name != "index.json"]
    if len(json_files) != 16:
        raise SystemExit(f"期望 16 个三芯片 JSON，实际 {len(json_files)} in {CASES}")
    rows = []
    for json_path in json_files:
        spec = json.loads(json_path.read_text(encoding="utf-8"))
        names = [c["name"] for c in spec["chips"]]
        if names != ["CPU", "GPU", "CIM"]:
            raise SystemExit(f"{json_path.name} 不是 CPU/GPU/CIM 三颗: {names}")
        mix_id = spec.get("mix") or json_path.stem
        print(f"== {mix_id}  lbfgs  {json_path.name}", flush=True)
        name = json_path.stem
        case_out = OUT / name
        case_out.mkdir(parents=True, exist_ok=True)
        shutil.copy2(json_path, case_out / json_path.name)

        case_name = spec["name"]
        atsim_dir = run_atsim(json_path, case_name, name)
        shutil.copy2(atsim_dir / "summary.json", case_out / "atsim_summary.json")
        if (atsim_dir / "probes.csv").exists():
            shutil.copy2(atsim_dir / "probes.csv", case_out / "atsim_probes.csv")
        if (atsim_dir / "die_midplane.png").exists():
            shutil.copy2(atsim_dir / "die_midplane.png", case_out / "atsim_die_midplane.png")

        atsim_sum = json.loads((atsim_dir / "summary.json").read_text(encoding="utf-8"))
        probes = {p["name"]: p["T_C"] for p in atsim_sum["probes"]}
        shifted_chips = atsim_sum["chips"]
        npz = atsim_dir / "T_field.npz"
        field = chip_stats_from_npz(npz, shifted_chips) if npz.exists() else {}

        steady = run_hotspot_from_atsim(case_name, name)
        hs_all = parse_steady(steady)
        hs_chips = {k: v for k, v in hs_all.items() if k in ("CPU", "GPU", "CIM")}
        shutil.copy2(steady, case_out / "hotspot.steady")

        hs_tmax = max(hs_chips.values())
        py_tmax = float(atsim_sum["Tmax_solid_C"])
        row = {
            "mix": mix_id,
            "algorithm": "lbfgs",
            "power_W": sum(float(c["power_W"]) for c in spec["chips"]),
            "HotSpot_Tmax_C": hs_tmax,
            "ATSim_Tmax_C": py_tmax,
            "dTmax_C": hs_tmax - py_tmax,
            "Tmax_rel_pct": rel_err(hs_tmax, py_tmax),
            "ATSim_elapsed_s": atsim_sum.get("elapsed_s"),
        }
        for chip in ("CPU", "GPU", "CIM"):
            hs_t = hs_chips.get(chip)
            ctr = probes.get(f"{chip}_center")
            mean = (field.get(chip) or {}).get("mean_C")
            mx = (field.get(chip) or {}).get("max_C")
            row[f"HS_{chip}_C"] = hs_t
            row[f"ATSim_{chip}_center_C"] = ctr
            row[f"ATSim_{chip}_mean_C"] = mean
            row[f"ATSim_{chip}_max_C"] = mx
            if hs_t is not None and mean is not None:
                row[f"dT_{chip}_mean_C"] = hs_t - mean
                row[f"rel_{chip}_mean_pct"] = rel_err(hs_t, mean)
            if hs_t is not None and ctr is not None:
                row[f"dT_{chip}_center_C"] = hs_t - ctr
                row[f"rel_{chip}_center_pct"] = rel_err(hs_t, ctr)
        rows.append(row)
        (case_out / "compare.json").write_text(
            json.dumps(row, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
        )
        print(
            f"   HS Tmax={hs_tmax:.3f}C  ATSim Tmax={py_tmax:.3f}C  "
            f"dT={row['dTmax_C']:.3f}C ({row['Tmax_rel_pct']:.2f}%)",
            flush=True,
        )

    table_path = OUT / "compare_table.csv"
    fields = list(rows[0].keys())
    with table_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    (OUT / "compare_table.json").write_text(
        json.dumps(rows, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )

    if BACKUP.exists():
        shutil.rmtree(BACKUP)
    shutil.copytree(OUT, BACKUP / "runs")
    shutil.copytree(CASES, BACKUP / "atsim_cases")
    shutil.copy2(table_path, BACKUP / "compare_table.csv")
    shutil.copy2(OUT / "compare_table.json", BACKUP / "compare_table.json")
    print(json.dumps(rows, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
