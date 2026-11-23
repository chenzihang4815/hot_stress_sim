#!/usr/bin/env python3
"""TAP-style SA: nudge / swap / jump / rotate.

Default inner cost = HotSpot Tmax (≈1 s). Stress used to be evaluated only
on init vs best (SKIP_STRESS). The thermal-vs-coupled comparison lives in
`sa_thermal_stress_compare.py` (does not overwrite this run directory).

Moves follow TAP-2.5D sim_annealing.py:
  close  = one chiplet by one grid step
  jump   = relocate one chiplet (optional 90° rotate)
  swap   = exchange two chiplets' lower-left if still legal
"""

from __future__ import annotations

import json
import math
import random
import shutil
import subprocess
import sys
from copy import deepcopy
from dataclasses import dataclass
from pathlib import Path

ROOT = Path("/home/czh/任务仿真产生功耗文件工具")
SRC = ROOT / "热仿真流水线" / "src"
sys.path.insert(0, str(SRC))

from geometry import Rect, packing_metrics  # noqa: E402
from hotspot_io import run_hotspot, write_hotspot_inputs  # noqa: E402

HS_BIN = ROOT / "HotSpot-master" / "hotspot"
HS_SCRIPTS = ROOT / "HotSpot-master" / "scripts"
ATSIM = ROOT / "热应力python工程_最终"
PYTHON = "/home/czh/miniconda3/envs/fenicsx/bin/python"
OUT = ROOT / "热仿真流水线" / "runs" / "sa_gapped_13chip"
BACKUP = ROOT / "成功运行备份" / "SA带布线间隙_13芯_2026-10-22"
INIT_LAYOUT = ROOT / "热仿真流水线" / "runs" / "large_package" / "pkg_5cpu_3gpu_5cim" / "layout.json"

GRAN_M = 0.0005  # 0.5 mm 网格
GAP_M = 0.0005  # 任意两颗至少 0.5 mm，给微凸/布线
GAP_GPU_M = 0.002  # GPU–GPU 至少 2 mm
N_STEPS = 120
T0 = 2.0
ALPHA = 0.95
SEED = 22
SKIP_STRESS = True  # 先看间隙热图；应力另评


@dataclass
class Chip:
    name: str
    kind: str
    x: float
    y: float
    w: float
    h: float
    power_W: float
    area_mm2: float
    # CPU 指向 McPAT 官方模板名，其余种类空着。评估按它查漏电表。
    model: str = ""

    def rect(self) -> Rect:
        return Rect(self.name, self.x, self.y, self.w, self.h, self.power_W, self.area_mm2, self.kind)

    def rotate90(self) -> None:
        self.w, self.h = self.h, self.w


def load_chips(path: Path) -> list[Chip]:
    rows = json.loads(path.read_text(encoding="utf-8"))
    return [
        Chip(r["name"], r["kind"], r["x_m"], r["y_m"], r["w_m"], r["h_m"], r["power_W"], r["area_mm2"])
        for r in rows
    ]


def pair_gap(a: Chip, b: Chip) -> float:
    """Clearance between rectangles; 0 if overlapping, negative if overlapping by that much."""
    dx = max(a.x, b.x) - min(a.x + a.w, b.x + b.w)
    dy = max(a.y, b.y) - min(a.y + a.h, b.y + b.h)
    if dx < 0 and dy < 0:
        return -min(-dx, -dy)
    if dx < 0:
        return dy
    if dy < 0:
        return dx
    return math.hypot(dx, dy)


def required_gap(a: Chip, b: Chip) -> float:
    if a.kind == "GPU" and b.kind == "GPU":
        return GAP_GPU_M
    return GAP_M


def min_clearance(chips: list[Chip]) -> tuple[float, float]:
    """Return (min gap any pair, min gap GPU-GPU). Overlap counts as negative."""
    g_all = math.inf
    g_gpu = math.inf
    for i, a in enumerate(chips):
        for b in chips[i + 1 :]:
            g = pair_gap(a, b)
            g_all = min(g_all, g)
            if a.kind == "GPU" and b.kind == "GPU":
                g_gpu = min(g_gpu, g)
    if g_gpu is math.inf:
        g_gpu = math.inf
    return (g_all if g_all is not math.inf else 0.0, g_gpu if g_gpu is not math.inf else math.inf)


def legal(chips: list[Chip], intp_w: float, intp_h: float) -> bool:
    for c in chips:
        if c.x < -1e-12 or c.y < -1e-12:
            return False
        if c.x + c.w > intp_w + 1e-9 or c.y + c.h > intp_h + 1e-9:
            return False
        if c.w <= 0 or c.h <= 0:
            return False
    for i, a in enumerate(chips):
        for b in chips[i + 1 :]:
            if pair_gap(a, b) + 1e-12 < required_gap(a, b):
                return False
    return True


def spread_to_gap(chips: list[Chip], intp_w: float, intp_h: float, rounds: int = 80) -> list[Chip]:
    """Push overlapping / too-close chiplets apart along the shortest axis."""
    cur = deepcopy(chips)
    shift_to_origin(cur)
    for _ in range(rounds):
        moved = False
        for i, a in enumerate(cur):
            for j, b in enumerate(cur[i + 1 :], start=i + 1):
                need = required_gap(a, b)
                g = pair_gap(a, b)
                if g + 1e-12 >= need:
                    continue
                # push along the axis with more overlap (or both if nested)
                ax = (a.x + a.w / 2) - (b.x + b.w / 2)
                ay = (a.y + a.h / 2) - (b.y + b.h / 2)
                push = need - g + GRAN_M
                if abs(ax) >= abs(ay):
                    step = push / 2 if ax != 0 else push / 2
                    if ax >= 0:
                        cur[i].x += step
                        cur[j].x -= step
                    else:
                        cur[i].x -= step
                        cur[j].x += step
                else:
                    step = push / 2
                    if ay >= 0:
                        cur[i].y += step
                        cur[j].y -= step
                    else:
                        cur[i].y -= step
                        cur[j].y += step
                moved = True
        shift_to_origin(cur)
        # snap to grid
        for c in cur:
            c.x = round(c.x / GRAN_M) * GRAN_M
            c.y = round(c.y / GRAN_M) * GRAN_M
        shift_to_origin(cur)
        if not moved and legal(cur, intp_w, intp_h):
            break
    return cur


def shift_to_origin(chips: list[Chip]) -> None:
    minx = min(c.x for c in chips)
    miny = min(c.y for c in chips)
    for c in chips:
        c.x -= minx
        c.y -= miny


def to_rects(chips: list[Chip]) -> list[Rect]:
    return [c.rect() for c in chips]


def eval_hotspot(chips: list[Chip], work: Path, rows: int = 64) -> float:
    if work.exists():
        shutil.rmtree(work)
    work.mkdir(parents=True)
    write_hotspot_inputs(to_rects(chips), work, rows=rows, cols=rows)
    thermal = run_hotspot(HS_BIN, work)
    tmax = thermal.get("Tmax_C")
    if tmax is None:
        raise RuntimeError("HotSpot no Tmax")
    return float(tmax)


def official_map(case_dir: Path, png: Path) -> None:
    all_grid = case_dir / "all.grid.steady"
    shutil.copy2(case_dir / "grid.steady", all_grid)
    subprocess.run(
        [sys.executable, str(HS_SCRIPTS / "split_grid_steady.py"), str(all_grid), "4", "64", "64"],
        cwd=str(case_dir),
        check=True,
        capture_output=True,
        text=True,
    )
    subprocess.run(
        [
            sys.executable,
            str(HS_SCRIPTS / "grid_thermal_map.py"),
            str(case_dir / "floorplan.flp"),
            str(case_dir / "all_layer0.grid.steady"),
            "64",
            "64",
            str(png),
        ],
        cwd=str(case_dir),
        check=True,
        capture_output=True,
        text=True,
    )


def propose(chips: list[Chip], intp_w: float, intp_h: float, rng: random.Random,
            mix: tuple[float, float, float, float] = (0.40, 0.25, 0.20, 0.15)) -> tuple[list[Chip], str]:
    n = len(chips)
    trial = deepcopy(chips)
    gpus = [k for k, c in enumerate(chips) if c.kind == "GPU"]
    r = rng.random()
    # mix = (nudge, jump, swap, rotate). Default 40/25/20/15.
    b_nudge, b_jump, b_swap = mix[0], mix[0] + mix[1], mix[0] + mix[1] + mix[2]
    if r < b_nudge:
        i = rng.choice(gpus) if gpus and rng.random() < 0.55 else rng.randrange(n)
        steps = rng.choice([1, 1, 2, 2, 3, 4])
        dx, dy = rng.choice([(GRAN_M, 0), (-GRAN_M, 0), (0, GRAN_M), (0, -GRAN_M)])
        trial[i].x += dx * steps
        trial[i].y += dy * steps
        op = f"nudge {trial[i].name} ({dx*steps*1e3:.1f},{dy*steps*1e3:.1f})mm"
    elif r < b_jump:
        i = rng.choice(gpus) if gpus and rng.random() < 0.6 else rng.randrange(n)
        max_x = max(0, int((intp_w - trial[i].w) / GRAN_M))
        max_y = max(0, int((intp_h - trial[i].h) / GRAN_M))
        trial[i].x = rng.randint(0, max_x) * GRAN_M
        trial[i].y = rng.randint(0, max_y) * GRAN_M
        if rng.random() < 0.4:
            trial[i].rotate90()
            op = f"jump+rot {trial[i].name}"
        else:
            op = f"jump {trial[i].name}"
    elif r < b_swap:
        i, j = rng.sample(range(n), 2)
        trial[i].x, trial[j].x = trial[j].x, trial[i].x
        trial[i].y, trial[j].y = trial[j].y, trial[i].y
        op = f"swap {trial[i].name}<->{trial[j].name}"
    else:
        i = rng.randrange(n)
        trial[i].rotate90()
        op = f"rotate {trial[i].name}"
    # Coordinates stay relative to the fixed 72 mm frame. Snapping back to the
    # origin would shift every chip, not just the one that moved, and a stored
    # thermal field is only valid on the grid it was computed on.
    return trial, op


def write_atsim_json(chips: list[Chip], path: Path, title: str) -> None:
    payload = {
        "name": path.stem,
        "description": title,
        "chips": [
            {
                "name": c.name,
                "x_mm": c.x * 1e3,
                "y_mm": c.y * 1e3,
                "w_mm": c.w * 1e3,
                "h_mm": c.h * 1e3,
                "t_um": 150.0,
                "power_W": c.power_W,
            }
            for c in chips
        ],
        "package": {
            "tim_um": 20,
            "t_spreader_mm": 1.0,
            "t_sink_mm": 6.9,
            "t_substrate_mm": 1.0,
            "t_interposer_um": 100.0,
            "r_convec_K_per_W": 0.1,
            "ambient_C": 45,
            "pitch_die_mm": 2.0,
            "pitch_coarse_mm": 8.0,
            "nz_die": 3,
            "nz_spreader": 3,
            "nz_sink": 4,
            "nz_tim": 1,
            "nz_substrate": 3,
            "nz_interposer": 2,
            "fill_die": False,
        },
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def run_atsim_stress(json_path: Path) -> dict:
    completed = subprocess.run(
        [PYTHON, str(ATSIM / "scripts" / "run_case3_stress.py"), str(json_path)],
        cwd=str(ATSIM),
        text=True,
        capture_output=True,
        check=False,
    )
    name = json_path.stem
    out = ATSIM / "results" / f"python_{name}_stress"
    log = json_path.parent / f"{name}_atsim.log"
    log.write_text(completed.stdout + "\n--- STDERR ---\n" + completed.stderr, encoding="utf-8")
    if completed.returncode != 0:
        raise RuntimeError(f"ATSim failed {json_path}\n{completed.stderr[-1500:]}")
    return json.loads((out / "summary.json").read_text(encoding="utf-8"))


def snapshot(chips: list[Chip]) -> list[dict]:
    return [
        {
            "name": c.name,
            "kind": c.kind,
            "x_m": c.x,
            "y_m": c.y,
            "w_m": c.w,
            "h_m": c.h,
            "area_mm2": c.area_mm2,
            "power_W": c.power_W,
        }
        for c in chips
    ]


def main() -> None:
    rng = random.Random(SEED)
    raw = load_chips(INIT_LAYOUT)
    shift_to_origin(raw)
    m_raw = packing_metrics(to_rects(raw))
    # Fixed 72 mm chip layer, same outline the thermal solver uses.
    intp_w = intp_h = 0.072
    chips0 = spread_to_gap(raw, intp_w, intp_h)
    if not legal(chips0, intp_w, intp_h):
        print("warn: spread_to_gap did not fully satisfy gaps; SA will keep searching", flush=True)
    m0 = packing_metrics(to_rects(chips0))
    g0, g0gpu = min_clearance(chips0)
    if OUT.exists():
        shutil.rmtree(OUT)
    OUT.mkdir(parents=True)

    print(
        f"gapped init outline {m0['width_m']*1e3:.2f} x {m0['height_m']*1e3:.2f} mm  n={len(chips0)}  "
        f"min_gap={g0*1e3:.2f}mm  GPU-GPU={g0gpu*1e3:.2f}mm  intp={intp_w*1e3:.1f}x{intp_h*1e3:.1f}mm",
        flush=True,
    )
    tmax0 = eval_hotspot(chips0, OUT / "eval_init")
    print(f"init HotSpot Tmax = {tmax0:.3f} C (after forcing routing gaps)", flush=True)

    current = deepcopy(chips0)
    tmax_cur = tmax0
    best = deepcopy(chips0)
    tmax_best = tmax0
    log = [{"step": 0, "op": "init", "accept": True, "Tmax_C": tmax0, "T": T0}]
    T = T0
    accepted = 0
    for step in range(1, N_STEPS + 1):
        trial, op = propose(current, intp_w, intp_h, rng)
        if not legal(trial, intp_w, intp_h):
            log.append({"step": step, "op": op, "accept": False, "reason": "illegal", "T": T})
            T *= ALPHA
            continue
        try:
            tmax_new = eval_hotspot(trial, OUT / "eval_try")
        except Exception as exc:
            log.append({"step": step, "op": op, "accept": False, "reason": str(exc)[:120], "T": T})
            T *= ALPHA
            continue
        delta = tmax_cur - tmax_new
        if delta >= 0:
            ap = 1.0
        else:
            ap = math.exp(delta / max(T, 1e-9))
        take = rng.random() < ap
        rec = {
            "step": step,
            "op": op,
            "Tmax_C": tmax_new,
            "delta_C": delta,
            "ap": ap,
            "accept": take,
            "T": T,
        }
        if take:
            current = trial
            tmax_cur = tmax_new
            accepted += 1
            if tmax_new < tmax_best:
                best = deepcopy(trial)
                tmax_best = tmax_new
                rec["new_best"] = True
                print(f"  step {step:3d} BEST {op}  Tmax={tmax_new:.3f}C", flush=True)
            else:
                print(f"  step {step:3d} accept {op}  Tmax={tmax_new:.3f}C", flush=True)
        else:
            print(f"  step {step:3d} reject {op}  Tmax={tmax_new:.3f}C  ap={ap:.3f}", flush=True)
        log.append(rec)
        T *= ALPHA

    gb, gbgpu = min_clearance(best)
    print(
        f"accepted {accepted}/{N_STEPS}  Tmax {tmax0:.3f} -> {tmax_best:.3f} C  "
        f"min_gap={gb*1e3:.2f}mm GPU-GPU={gbgpu*1e3:.2f}mm",
        flush=True,
    )

    init_dir = OUT / "layout_init"
    best_dir = OUT / "layout_best"
    write_hotspot_inputs(to_rects(chips0), init_dir)
    run_hotspot(HS_BIN, init_dir)
    official_map(init_dir, OUT / "init_hotspot.png")
    write_hotspot_inputs(to_rects(best), best_dir)
    run_hotspot(HS_BIN, best_dir)
    official_map(best_dir, OUT / "best_hotspot.png")

    atsim_cases = ATSIM / "cases" / "sa_feedback"
    write_atsim_json(chips0, atsim_cases / "sa13_gapped_init.json", "13-chip after forcing routing gaps")
    write_atsim_json(best, atsim_cases / "sa13_gapped_best.json", "13-chip SA-best with routing gaps")

    if SKIP_STRESS:
        result = {
            "moves": "TAP-like nudge/swap/jump/rotate with min clearance",
            "gap_any_mm": GAP_M * 1e3,
            "gap_gpu_mm": GAP_GPU_M * 1e3,
            "inner_cost": "HotSpot Tmax",
            "n_steps": N_STEPS,
            "accepted": accepted,
            "seed": SEED,
            "HotSpot_Tmax_init_C": tmax0,
            "HotSpot_Tmax_best_C": tmax_best,
            "dTmax_C": tmax_best - tmax0,
            "init_min_gap_mm": g0 * 1e3,
            "init_gpu_gap_mm": (g0gpu * 1e3) if g0gpu is not math.inf else None,
            "best_min_gap_mm": gb * 1e3,
            "best_gpu_gap_mm": (gbgpu * 1e3) if gbgpu is not math.inf else None,
            "init_layout": snapshot(chips0),
            "best_layout": snapshot(best),
            "log": log,
        }
        (OUT / "sa_result.json").write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
        if BACKUP.exists():
            shutil.rmtree(BACKUP)
        shutil.copytree(OUT, BACKUP)
        print(json.dumps({k: result[k] for k in (
            "HotSpot_Tmax_init_C", "HotSpot_Tmax_best_C", "dTmax_C", "accepted",
            "init_min_gap_mm", "init_gpu_gap_mm", "best_min_gap_mm", "best_gpu_gap_mm",
        )}, indent=2))
        return

    print("ATSim stress on initial ...", flush=True)
    sum_init = run_atsim_stress(atsim_cases / "sa13_gapped_init.json")
    print("ATSim stress on SA-best ...", flush=True)
    sum_best = run_atsim_stress(atsim_cases / "sa13_gapped_best.json")

    def stress_row(s: dict) -> dict:
        centers = {
            p["name"].replace("_center", ""): {"T_C": p["T_C"], "vm_MPa": p["vm_MPa"], "uz_um": p["uz_um"]}
            for p in s["probes"]
            if p["name"].endswith("_center") and not p["name"].startswith(("Sink", "Spreader", "Interposer", "Substrate", "TIM"))
        }
        return {
            "ATSim_Tmax_C": s["Tmax_C"],
            "vm_max_package_MPa": s["vm_max_MPa"],
            "warpage_PV_um": s["warpage"]["w_pv_um"],
            "u_max_um": s["u_max_um"],
            "chip_centers": centers,
        }

    result = {
        "moves": "TAP-like nudge/swap/jump/rotate",
        "inner_cost": "HotSpot Tmax",
        "outer_cost": "ATSim Tmax + interposer warpage PV + chip-center von Mises (evaluated init vs best only)",
        "n_steps": N_STEPS,
        "accepted": accepted,
        "seed": SEED,
        "HotSpot_Tmax_init_C": tmax0,
        "HotSpot_Tmax_best_C": tmax_best,
        "dTmax_C": tmax_best - tmax0,
        "init_stress": stress_row(sum_init),
        "best_stress": stress_row(sum_best),
        "init_layout": snapshot(chips0),
        "best_layout": snapshot(best),
        "log": log,
    }
    # combined score for the user: lower better
    def J(hs_t, st):
        return (
            1.0 * (hs_t - 45.0)
            + 2.0 * st["warpage_PV_um"]
            + 0.02 * max(c["vm_MPa"] for c in st["chip_centers"].values())
        )

    result["J_init"] = J(tmax0, result["init_stress"])
    result["J_best"] = J(tmax_best, result["best_stress"])
    (OUT / "sa_result.json").write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    if BACKUP.exists():
        shutil.rmtree(BACKUP)
    shutil.copytree(OUT, BACKUP)
    print(json.dumps({k: result[k] for k in ("HotSpot_Tmax_init_C", "HotSpot_Tmax_best_C", "dTmax_C", "J_init", "J_best", "accepted")}, indent=2))


if __name__ == "__main__":
    main()
