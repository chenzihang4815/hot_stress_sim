#!/usr/bin/env python3
"""One-scheme heat-only vs coupled SA. Run serially (12 GB).

Usage:
  python run_sa_compare_scheme.py pkg_4cpu_2gpu_4cim
  python run_sa_compare_scheme.py pkg_5cpu_4gpu_5cim
"""

from __future__ import annotations

import json
import math
import random
import shutil
import sys
import time
from copy import deepcopy
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path("/home/czh/任务仿真产生功耗文件工具")
SRC = ROOT / "热仿真流水线" / "src"
sys.path.insert(0, str(SRC))

import sa_thermal_feedback as sa  # noqa: E402
from atsim_sa_proxy import eval_atsim_proxy, write_proxy_json  # noqa: E402
from geometry import packing_metrics  # noqa: E402
from hotspot_io import run_hotspot, write_hotspot_inputs  # noqa: E402

SCHEMES = ROOT / "热仿真流水线" / "runs" / "sa_compare_schemes"
N_STEPS = 16
SEED = 22
TREF = 25.0
SREF = 100.0  # MPa; worst normal-chip interior von Mises (not center point)
WREF = 4.0
LAMBDA_VM = 1.0
MU_WARP = 0.5


def J_heat(tmax: float) -> float:
    return (tmax - 45.0) / TREF


def J_full(tmax: float, vm: float, wpv: float) -> float:
    return J_heat(tmax) + LAMBDA_VM * (vm / SREF) + MU_WARP * (wpv / WREF)


def vm_of(st: dict) -> float:
    """Worst interior von Mises among normal chiplets (CPU/GPU, CIM ≥ 1 mm²)."""
    v = st.get("vm_normal_chip_max_MPa")
    if v is None:
        v = st.get("vm_center_max_MPa")
    if v is None:
        raise RuntimeError("stress proxy missing chip von Mises")
    return float(v)


def pkg_overrides(scheme: dict) -> dict:
    out = {
        "pitch_die_mm": scheme.get("atsim_pitch_die_mm") or 2.0,
        "pitch_coarse_mm": scheme.get("atsim_pitch_coarse_mm") or 8.0,
    }
    nz = scheme.get("atsim_nz") or {}
    out.update({k: v for k, v in nz.items() if v is not None})
    return out


def eval_stress(chips, tag: str, work: Path, overrides: dict) -> dict:
    jp = work / f"{tag}.json"
    write_proxy_json(chips, jp, f"SA proxy {tag}", fill_die=True, package_overrides=overrides)
    t0 = time.time()
    st = eval_atsim_proxy(jp, verbose=False)
    st["wall_s"] = time.time() - t0
    (work / f"{tag}_stress.json").write_text(
        json.dumps(st, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    return st


def anneal(arm, chips0, intp_w, intp_h, rng, out, overrides, log_md):
    arm_dir = out / arm
    if arm_dir.exists():
        shutil.rmtree(arm_dir)
    arm_dir.mkdir(parents=True)
    jsonl = arm_dir / "steps.jsonl"
    coupled = arm == "coupled"

    tmax0 = sa.eval_hotspot(chips0, arm_dir / "eval_init")
    st0 = eval_stress(chips0, "init", arm_dir, overrides)
    vm0, wpv0 = vm_of(st0), st0["w_pv_um"]
    j0 = J_full(tmax0, vm0, wpv0) if coupled else J_heat(tmax0)
    log_md(
        f"- `{arm}` init  Tmax={tmax0:.3f}C  vm_chipmax={vm0:.2f} MPa "
        f"({st0.get('vm_normal_chip_name')})  "
        f"w_pv={wpv0:.3f} um  J={j0:.4f}  stress={st0['elapsed_s']:.1f}s"
    )

    current = deepcopy(chips0)
    tmax_cur, vm_cur, wpv_cur, j_cur = tmax0, vm0, wpv0, j0
    best = deepcopy(chips0)
    tmax_best, vm_best, wpv_best, j_best = tmax0, vm0, wpv0, j0
    st_best = st0
    T = sa.T0
    accepted = 0
    n_hs = 0
    n_stress = 1
    log = [{
        "step": 0, "op": "init", "accept": True, "Tmax_C": tmax0,
        "vm_chip_max_MPa": vm0, "vm_chip_name": st0.get("vm_normal_chip_name"),
        "w_pv_um": wpv0, "J": j0, "T": T,
    }]
    with jsonl.open("a", encoding="utf-8") as f:
        f.write(json.dumps(log[-1], ensure_ascii=False) + "\n")

    for step in range(1, N_STEPS + 1):
        trial, op = sa.propose(current, intp_w, intp_h, rng)
        rec = {"step": step, "op": op, "T": T, "arm": arm}
        if not sa.legal(trial, intp_w, intp_h):
            rec.update({"accept": False, "reason": "illegal"})
            log.append(rec)
            with jsonl.open("a", encoding="utf-8") as f:
                f.write(json.dumps(rec, ensure_ascii=False) + "\n")
            T *= sa.ALPHA
            continue
        try:
            tmax_new = sa.eval_hotspot(trial, arm_dir / "eval_try")
            n_hs += 1
        except Exception as exc:
            rec.update({"accept": False, "reason": str(exc)[:160]})
            log.append(rec)
            with jsonl.open("a", encoding="utf-8") as f:
                f.write(json.dumps(rec, ensure_ascii=False) + "\n")
            T *= sa.ALPHA
            continue

        st_new = None
        vm_new, wpv_new = vm_cur, wpv_cur
        if coupled:
            st_new = eval_stress(trial, f"try_{step:03d}", arm_dir, overrides)
            n_stress += 1
            vm_new, wpv_new = vm_of(st_new), st_new["w_pv_um"]
        j_new = J_full(tmax_new, vm_new, wpv_new) if coupled else J_heat(tmax_new)
        delta = j_cur - j_new
        ap = 1.0 if delta >= 0 else math.exp(delta / max(T, 1e-9))
        take = rng.random() < ap
        rec.update({
            "Tmax_C": tmax_new, "vm_chip_max_MPa": vm_new,
            "vm_chip_name": None if st_new is None else st_new.get("vm_normal_chip_name"),
            "w_pv_um": wpv_new,
            "J": j_new, "delta_J": delta, "ap": ap, "accept": take, "stressed": st_new is not None,
        })
        if take:
            current = trial
            tmax_cur, vm_cur, wpv_cur, j_cur = tmax_new, vm_new, wpv_new, j_new
            accepted += 1
            if j_new < j_best - 1e-12:
                best = deepcopy(trial)
                tmax_best, vm_best, wpv_best, j_best = tmax_new, vm_new, wpv_new, j_new
                if st_new is not None:
                    st_best = st_new
                rec["new_best"] = True
                log_md(
                    f"  `{arm}` step {step:3d} BEST {op}  Tmax={tmax_new:.3f}C  "
                    f"vm={vm_new:.2f}  wpv={wpv_new:.3f}  J={j_new:.4f}"
                )
            else:
                log_md(f"  `{arm}` step {step:3d} accept {op}  Tmax={tmax_new:.3f}C  J={j_new:.4f}")
        else:
            log_md(
                f"  `{arm}` step {step:3d} reject {op}  Tmax={tmax_new:.3f}C  "
                f"J={j_new:.4f}  ap={ap:.3f}"
            )
        log.append(rec)
        with jsonl.open("a", encoding="utf-8") as f:
            f.write(json.dumps(rec, ensure_ascii=False) + "\n")
        T *= sa.ALPHA

    if not coupled:
        log_md(f"- `{arm}` rescoring best with silicon-fill stress …")
        st_best = eval_stress(best, "best", arm_dir, overrides)
        vm_best, wpv_best = vm_of(st_best), st_best["w_pv_um"]
        n_stress += 1

    write_hotspot_inputs(sa.to_rects(chips0), arm_dir / "layout_init")
    run_hotspot(sa.HS_BIN, arm_dir / "layout_init")
    sa.official_map(arm_dir / "layout_init", arm_dir / "init_hotspot.png")
    write_hotspot_inputs(sa.to_rects(best), arm_dir / "layout_best")
    run_hotspot(sa.HS_BIN, arm_dir / "layout_best")
    sa.official_map(arm_dir / "layout_best", arm_dir / "best_hotspot.png")

    g0, g0gpu = sa.min_clearance(chips0)
    gb, gbgpu = sa.min_clearance(best)
    m0 = packing_metrics(sa.to_rects(chips0))
    mb = packing_metrics(sa.to_rects(best))
    summary = {
        "arm": arm,
        "coupled": coupled,
        "n_steps": N_STEPS,
        "accepted": accepted,
        "n_hotspot": n_hs,
        "n_stress": n_stress,
        "seed": SEED,
        "init": {
            "Tmax_C": tmax0,
            "vm_chip_max_MPa": vm0,
            "vm_chip_name": st0.get("vm_normal_chip_name"),
            "vm_center_mean_MPa": st0.get("vm_center_mean_MPa"),
            "w_pv_um": wpv0,
            "J_heat": J_heat(tmax0),
            "J_full": J_full(tmax0, vm0, wpv0),
            "min_gap_mm": g0 * 1e3,
            "gpu_gap_mm": None if g0gpu is math.inf else g0gpu * 1e3,
            "outline_mm": [m0["width_m"] * 1e3, m0["height_m"] * 1e3],
            "stress": st0,
        },
        "best": {
            "Tmax_C": tmax_best,
            "vm_chip_max_MPa": vm_best,
            "vm_chip_name": st_best.get("vm_normal_chip_name"),
            "vm_center_mean_MPa": st_best.get("vm_center_mean_MPa"),
            "w_pv_um": wpv_best,
            "J_heat": J_heat(tmax_best),
            "J_full": J_full(tmax_best, vm_best, wpv_best),
            "min_gap_mm": gb * 1e3,
            "gpu_gap_mm": None if gbgpu is math.inf else gbgpu * 1e3,
            "outline_mm": [mb["width_m"] * 1e3, mb["height_m"] * 1e3],
            "stress": st_best,
        },
        "dTmax_C": tmax_best - tmax0,
        "d_vm_chip_max_MPa": vm_best - vm0,
        "d_w_pv_um": wpv_best - wpv0,
        "init_layout": sa.snapshot(chips0),
        "best_layout": sa.snapshot(best),
        "log": log,
    }
    (arm_dir / "summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    log_md(
        f"- `{arm}` done  Tmax {tmax0:.3f}→{tmax_best:.3f}C  "
        f"vm {vm0:.2f}→{vm_best:.2f} MPa  wpv {wpv0:.3f}→{wpv_best:.3f} um  "
        f"accepted {accepted}/{N_STEPS}"
    )
    return summary


def main() -> None:
    if len(sys.argv) < 2:
        raise SystemExit("usage: run_sa_compare_scheme.py <scheme_id>")
    sid = sys.argv[1]
    scheme_path = SCHEMES / sid / "scheme.json"
    layout_path = SCHEMES / sid / "layout.json"
    if not scheme_path.exists() or not layout_path.exists():
        raise SystemExit(f"missing scheme files under {SCHEMES / sid}")
    scheme = json.loads(scheme_path.read_text(encoding="utf-8"))
    overrides = pkg_overrides(scheme)
    out = ROOT / "热仿真流水线" / "runs" / f"sa_heat_vs_coupled_{sid}"
    backup = ROOT / "成功运行备份" / f"SA热vs热应力对照_{sid}_2026-10-23"
    t_wall = time.time()
    if out.exists():
        shutil.rmtree(out)
    out.mkdir(parents=True)
    log_path = out / "EXPERIMENT.md"

    def log_md(line: str) -> None:
        with log_path.open("a", encoding="utf-8") as f:
            f.write(line.rstrip() + "\n")
        print(line, flush=True)

    log_path.write_text(
        f"# Heat-only vs heat+stress SA — {sid}\n\n"
        f"- started (UTC): {datetime.now(timezone.utc).isoformat()}\n"
        f"- title: {scheme.get('title')}\n"
        f"- why: {scheme.get('why')}\n"
        f"- n_chiplets: {scheme.get('n_chiplets')}\n"
        f"- seed={SEED}  N_STEPS={N_STEPS}  T0={sa.T0}  alpha={sa.ALPHA}\n"
        f"- gaps: any {sa.GAP_M*1e3:.1f} mm, GPU–GPU {sa.GAP_GPU_M*1e3:.1f} mm\n"
        f"- coupled J = (Tmax-45)/{TREF} + {LAMBDA_VM}*vm_chip_max/{SREF} + {MU_WARP}*w_pv/{WREF}\n"
        f"- vm = worst interior von Mises among normal chiplets (CPU/GPU, CIM area≥1 mm²)\n"
        f"- ATSim package overrides: {json.dumps(overrides)}\n"
        f"- one scheme at a time (12 GB); CG only, no spsolve\n\n"
        "## Run log\n\n",
        encoding="utf-8",
    )

    raw = sa.load_chips(layout_path)
    sa.shift_to_origin(raw)
    m_raw = packing_metrics(sa.to_rects(raw))
    n_gpu = sum(1 for c in raw if c.kind == "GPU")
    intp_w = m_raw["width_m"] + sa.GAP_GPU_M * max(n_gpu, 1) + 0.008
    intp_h = m_raw["height_m"] + sa.GAP_GPU_M * max(n_gpu, 1) + 0.008
    chips0 = sa.spread_to_gap(raw, intp_w, intp_h)
    g0, g0gpu = sa.min_clearance(chips0)
    log_md(
        f"- start outline {packing_metrics(sa.to_rects(chips0))['width_m']*1e3:.2f} x "
        f"{packing_metrics(sa.to_rects(chips0))['height_m']*1e3:.2f} mm  "
        f"min_gap={g0*1e3:.2f} mm  GPU-GPU={g0gpu*1e3:.2f} mm"
    )
    (out / "layout_start.json").write_text(
        json.dumps(sa.snapshot(chips0), indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    (out / "scheme.json").write_text(json.dumps(scheme, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    heat = anneal("heat_only", chips0, intp_w, intp_h, random.Random(SEED), out, overrides, log_md)
    coupled = anneal("coupled", chips0, intp_w, intp_h, random.Random(SEED), out, overrides, log_md)

    table = {
        "scheme": sid,
        "title": scheme.get("title"),
        "n_chiplets": scheme.get("n_chiplets"),
        "n_steps": N_STEPS,
        "seed": SEED,
        "package_overrides": overrides,
        "vm_metric": "worst interior von Mises among CPU/GPU and CIM≥1mm2",
        "heat_only": {
            "Tmax_init_C": heat["init"]["Tmax_C"],
            "Tmax_best_C": heat["best"]["Tmax_C"],
            "vm_chip_max_init_MPa": heat["init"]["vm_chip_max_MPa"],
            "vm_chip_max_best_MPa": heat["best"]["vm_chip_max_MPa"],
            "vm_chip_name_best": heat["best"].get("vm_chip_name"),
            "w_pv_init_um": heat["init"]["w_pv_um"],
            "w_pv_best_um": heat["best"]["w_pv_um"],
            "J_full_best": heat["best"]["J_full"],
            "accepted": heat["accepted"],
        },
        "coupled": {
            "Tmax_init_C": coupled["init"]["Tmax_C"],
            "Tmax_best_C": coupled["best"]["Tmax_C"],
            "vm_chip_max_init_MPa": coupled["init"]["vm_chip_max_MPa"],
            "vm_chip_max_best_MPa": coupled["best"]["vm_chip_max_MPa"],
            "vm_chip_name_best": coupled["best"].get("vm_chip_name"),
            "w_pv_init_um": coupled["init"]["w_pv_um"],
            "w_pv_best_um": coupled["best"]["w_pv_um"],
            "J_full_best": coupled["best"]["J_full"],
            "accepted": coupled["accepted"],
        },
        "delta_coupled_minus_heat": {
            "dTmax_C": coupled["best"]["Tmax_C"] - heat["best"]["Tmax_C"],
            "d_vm_chip_max_MPa": coupled["best"]["vm_chip_max_MPa"] - heat["best"]["vm_chip_max_MPa"],
            "d_w_pv_um": coupled["best"]["w_pv_um"] - heat["best"]["w_pv_um"],
        },
        "wall_s": time.time() - t_wall,
        "finished_utc": datetime.now(timezone.utc).isoformat(),
    }
    (out / "compare_table.json").write_text(
        json.dumps(table, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    log_md("\n## Table\n")
    log_md("| arm | Tmax init→best | vm_chip_max init→best | w_pv init→best | accepted |")
    log_md("|---|---|---|---|---|")
    for name, s in (("heat_only", heat), ("coupled", coupled)):
        log_md(
            f"| {name} | {s['init']['Tmax_C']:.3f}→{s['best']['Tmax_C']:.3f} | "
            f"{s['init']['vm_chip_max_MPa']:.2f}→{s['best']['vm_chip_max_MPa']:.2f} "
            f"({s['best'].get('vm_chip_name')}) | "
            f"{s['init']['w_pv_um']:.3f}→{s['best']['w_pv_um']:.3f} | {s['accepted']}/{N_STEPS} |"
        )
    d = table["delta_coupled_minus_heat"]
    log_md(
        f"\ncoupled − heat_only at best: dTmax={d['dTmax_C']:.3f} C, "
        f"d_vm={d['d_vm_chip_max_MPa']:.2f} MPa, d_wpv={d['d_w_pv_um']:.3f} um"
    )
    log_md(f"\nwall {table['wall_s']:.1f}s  finished {table['finished_utc']}")
    if backup.exists():
        shutil.rmtree(backup)
    shutil.copytree(out, backup)
    print(json.dumps(table, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
