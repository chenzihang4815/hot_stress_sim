#!/usr/bin/env python3
"""Reuse existing ATSim heat fields, solve elasticity, plot die-midplane stress.

Cases:
  voids  = sa13_gapped_best          (fill_die false)
  sifill = sa13_gapped_best_sifill   (silicon fill of chip bbox)
"""

from __future__ import annotations

import json
import shutil
import sys
from pathlib import Path

import numpy as np

ROOT = Path("/home/czh/任务仿真产生功耗文件工具")
ATSIM = ROOT / "热应力python工程_最终"
SA = ROOT / "热仿真流水线" / "runs" / "sa_gapped_13chip"
BAK = ROOT / "成功运行备份" / "SA带布线间隙_13芯_2026-10-22"

sys.path.insert(0, str(ATSIM))
sys.path.insert(0, str(ATSIM / "scripts"))

from atsim.geometry import load_case, material_field  # noqa: E402
from atsim.fvm3d import build_grid  # noqa: E402
from atsim.elasticity import solve_elasticity  # noqa: E402
from atsim.visualize import die_mid_k  # noqa: E402
from atsim.warpage import package_warpage, save_warpage_csv, warpage_for_json  # noqa: E402
from run_case3_stress import (  # noqa: E402
    plot_field,
    plot_overview,
    plot_probe_bars,
    plot_warpage,
    plot_xz,
    sample_mech,
)

JOBS = [
    {
        "tag": "voids",
        "json": ATSIM / "cases" / "sa_feedback" / "sa13_gapped_best.json",
        "npz": ATSIM / "results" / "python_sa13_gapped_best" / "T_field.npz",
        "out": ATSIM / "results" / "python_sa13_gapped_best_stress",
    },
    {
        "tag": "sifill",
        "json": ATSIM / "cases" / "sa_feedback" / "sa13_gapped_best_sifill.json",
        "npz": ATSIM / "results" / "python_sa13_gapped_best_sifill" / "T_field.npz",
        "out": ATSIM / "results" / "python_sa13_gapped_best_sifill_stress",
    },
]


def thermal_from_npz(case_path: Path, npz_path: Path) -> dict:
    case = load_case(str(case_path))
    grid = build_grid(case)
    data = np.load(npz_path)
    T = data["T"]
    shape = (int(grid["nx"]), int(grid["ny"]), int(grid["nz"]))
    if tuple(T.shape) != shape:
        raise SystemExit(f"{npz_path.name} T{T.shape} != grid {shape}")
    return {
        "T": T,
        "grid": grid,
        "case": case,
        "mat_names": material_field(case, grid),
        "elapsed_s": 0.0,
        "timing": {"reused_npz": str(npz_path), "total_s": 0.0},
        "solver": "fvm3d_reused",
    }


def chip_vm_stats(el) -> list[dict]:
    case = el["case"]
    grid = el["grid"]
    k = die_mid_k({"grid": grid, "case": case, "T": el["thermal"]["T"]})
    xc, yc = grid["xc"], grid["yc"]
    vm = el["vm"][:, :, k]
    solid = grid["solid"][:, :, k]
    ox = getattr(case, "chip_ox", 0.0)
    oy = getattr(case, "chip_oy", 0.0)
    rows = []
    for b in case.blocks:
        x0, y0, w, h = ox + b.x0, oy + b.y0, b.w, b.h
        ix = np.where((xc >= x0 - 1e-9) & (xc <= x0 + w + 1e-9))[0]
        jy = np.where((yc >= y0 - 1e-9) & (yc <= y0 + h + 1e-9))[0]
        if ix.size == 0 or jy.size == 0:
            rows.append({"name": b.name, "n": 0})
            continue
        sub = vm[np.ix_(ix, jy)]
        msk = solid[np.ix_(ix, jy)]
        vals = sub[msk] * 1e-6
        if vals.size == 0:
            rows.append({"name": b.name, "n": 0})
            continue
        rows.append(
            {
                "name": b.name,
                "n": int(vals.size),
                "vm_mean_MPa": float(vals.mean()),
                "vm_max_MPa": float(vals.max()),
                "vm_p95_MPa": float(np.percentile(vals, 95)),
            }
        )
    return rows


def run_one(job: dict) -> dict:
    out = job["out"]
    out.mkdir(parents=True, exist_ok=True)
    print("=" * 64, flush=True)
    print("tag", job["tag"], "json", job["json"].name, flush=True)
    th = thermal_from_npz(job["json"], job["npz"])
    case = th["case"]
    solid = th["grid"]["solid"]
    print(
        "reused T  shape=%s  Tmax=%.3f C  solid=%d" % (
            th["T"].shape,
            float(th["T"][solid].max() - 273.15),
            int(solid.sum()),
        ),
        flush=True,
    )
    el = solve_elasticity(th, T_ref=case.T_amb, verbose=True)
    umax = float(el["umag"][solid].max())
    uzmax = float(np.abs(el["ucz"][solid]).max())
    vmmax = float(el["vm"][solid].max())
    warp = package_warpage(el)
    rows = sample_mech(el)
    chip_vm = chip_vm_stats(el)

    np.savez_compressed(
        out / "stress_die.npz",
        vm=el["vm"],
        sxx=el["sxx"],
        syy=el["syy"],
        ucx=el["ucx"],
        ucy=el["ucy"],
        ucz=el["ucz"],
        umag=el["umag"],
        xc=el["grid"]["xc"],
        yc=el["grid"]["yc"],
        zc=el["grid"]["zc"],
        solid=el["grid"]["solid"],
    )
    save_warpage_csv(el, warp, str(out / "warpage_interposer.csv"))
    import csv

    with (out / "mech_probes.csv").open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    plot_field(
        el["grid"], case, el["vm"], str(out / "vm_die.png"),
        "von Mises die mid  max=%.1f MPa" % (vmmax * 1e-6),
        cmap="inferno", scale=1e-6, unit="MPa",
    )
    plot_field(
        el["grid"], case, el["ucz"], str(out / "uz_die.png"),
        "uz die mid  max=%.2f um" % (uzmax * 1e6),
        cmap="coolwarm", scale=1e6, unit="um",
    )
    plot_field(
        el["grid"], case, el["umag"], str(out / "umag_die.png"),
        "|u| die mid  max=%.2f um" % (umax * 1e6),
        cmap="viridis", scale=1e6, unit="um",
    )
    plot_overview(th, el, str(out / "overview_chips.png"))
    plot_xz(th, str(out / "xz_gpu.png"))
    plot_probe_bars(rows, str(out / "probe_bars.png"))
    plot_warpage(el, warp, str(out / "warpage_interposer.png"))

    summary = {
        "case": case.name,
        "tag": job["tag"],
        "json": str(job["json"]),
        "T_ref_C": case.T_amb - 273.15,
        "Tmax_C": float(th["T"][solid].max() - 273.15),
        "ndof": el["ndof"],
        "u_max_um": umax * 1e6,
        "uz_max_um": uzmax * 1e6,
        "vm_max_MPa": vmmax * 1e-6,
        "warpage": warpage_for_json(warp),
        "elast_s": el["elapsed_s"],
        "elast_timing": el.get("timing"),
        "probes": rows,
        "chip_vm_die_mid": chip_vm,
        "bc": el.get("bc", "substrate_bottom_fixed"),
        "element": el.get("element") or (el.get("timing") or {}).get("element"),
        "coupling": "one-way thermal to mechanical, T reused from npz",
        "chips": [
            {
                "name": b.name,
                "x_mm": b.x0 * 1e3,
                "y_mm": b.y0 * 1e3,
                "w_mm": b.w * 1e3,
                "h_mm": b.h * 1e3,
                "power_W": b.power,
            }
            for b in case.blocks
        ],
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(
        "tag=%s  Tmax=%.2f C  vm_max=%.1f MPa  |u|=%.2f um  elast=%.1fs"
        % (job["tag"], summary["Tmax_C"], summary["vm_max_MPa"], summary["u_max_um"], el["elapsed_s"]),
        flush=True,
    )
    for r in chip_vm:
        if r.get("n"):
            print(
                "  %-16s  vm mean/p95/max = %7.1f / %7.1f / %7.1f MPa"
                % (r["name"], r["vm_mean_MPa"], r["vm_p95_MPa"], r["vm_max_MPa"]),
                flush=True,
            )
    return summary


def main() -> None:
    wanted = None
    if len(sys.argv) > 1:
        wanted = sys.argv[1].strip().lower()
        if wanted not in {j["tag"] for j in JOBS}:
            raise SystemExit("usage: run_gapped_stress_viz.py [voids|sifill]")
    rows = []
    for job in JOBS:
        if wanted and job["tag"] != wanted:
            continue
        rows.append(run_one(job))
    if not rows:
        raise SystemExit("no jobs ran")
    if wanted:
        print("single-tag done:", wanted)
        return
    table = []
    for s in rows:
        table.append(
            {
                "tag": s["tag"],
                "Tmax_C": s["Tmax_C"],
                "vm_max_MPa": s["vm_max_MPa"],
                "u_max_um": s["u_max_um"],
                "uz_max_um": s["uz_max_um"],
                "w_pv_um": (s.get("warpage") or {}).get("w_pv_um"),
                "elast_s": s["elast_s"],
            }
        )
    (SA / "atsim_stress_voids_vs_sifill.json").write_text(
        json.dumps({"cases": table, "details": rows}, indent=2, default=str) + "\n",
        encoding="utf-8",
    )
    SA.mkdir(parents=True, exist_ok=True)
    BAK.mkdir(parents=True, exist_ok=True)
    mapping = {
        "voids": JOBS[0]["out"],
        "sifill": JOBS[1]["out"],
    }
    for tag, src in mapping.items():
        for name in (
            "vm_die.png",
            "uz_die.png",
            "umag_die.png",
            "overview_chips.png",
            "warpage_interposer.png",
            "probe_bars.png",
            "xz_gpu.png",
        ):
            p = src / name
            if p.exists():
                shutil.copy2(p, SA / f"atsim_{tag}_{name}")
                shutil.copy2(p, BAK / f"atsim_{tag}_{name}")
        shutil.copy2(src / "summary.json", SA / f"atsim_{tag}_stress_summary.json")
        shutil.copy2(src / "summary.json", BAK / f"atsim_{tag}_stress_summary.json")
    print("copied plots to", SA)
    print(json.dumps(table, indent=2))


if __name__ == "__main__":
    main()
