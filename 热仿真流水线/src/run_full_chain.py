#!/usr/bin/env python3
"""Full chain on area-matched CPU+GPU packages.

  A  architecture feedback: McPAT starts at 350 K, HotSpot Tmax rewrites the
     10 K bin, heat-only SA, repeat until the bin is stable. That layout is L*.
  B  real ATSim stress on L* at the converged bin.
  C  frequency: sweep s, real ATSim each point. Two winners inside
     Tmax<90 C and sigma<200 MPa — lowest CPU leak fraction, and highest
     package throughput. Both are kept.
  D  for each winner: re-apply s on L* as the initial coordinate, heat-only SA
     again, then close the temperature bin once more. Real ATSim on the final
     layout only (12 GB: never inside the SA loop, never beside mcpat).

Leakage is McPAT's own number. Official templates are HP (device_type 0): the
only temperature dependence is I_off, gate leakage is flat, so each bin is the
official 380 K subthreshold scaled by the ITRS HP table with gate fixed. 350 K
is on that table. The small gem5 CPUs are LOP and are not used here — their
area is ~14 mm² against GPUs of 392–471 mm².

Every step appends one row to steps.jsonl and EXPERIMENT.md:
temperature, stress, layout, leakage, ops, power, and the delta from the
previous step.
"""

from __future__ import annotations

import json
import math
import random
import shutil
import sys
import time
from copy import deepcopy

import numpy as np
from datetime import datetime, timezone
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager

ROOT = Path("/home/czh/任务仿真产生功耗文件工具")
SRC = ROOT / "热仿真流水线" / "src"
sys.path.insert(0, str(SRC))

import sa_thermal_feedback as sa  # noqa: E402
from hotspot_io import run_hotspot, write_hotspot_inputs  # noqa: E402
from atsim_sa_proxy import eval_atsim_proxy, write_proxy_json  # noqa: E402
from blocks import Block  # noqa: E402
from geometry import packing_metrics  # noqa: E402
from lbfgs import random_lbfgs_layout  # noqa: E402
from run_sa_compare_scheme import vm_of  # noqa: E402

_font_path = "/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc"
font_manager.fontManager.addfont(_font_path)
_prop = font_manager.FontProperties(fname=_font_path)
plt.rcParams["font.family"] = _prop.get_name()
plt.rcParams["axes.unicode_minus"] = False

HS_BIN = ROOT / "HotSpot-master" / "hotspot"
TPL = ROOT / "成功运行备份" / "McPAT" / "official_templates"
OUT = ROOT / "热仿真流水线" / "runs" / "full_chain_sa120"
BACKUP = ROOT / "成功运行备份" / "全链_SA120步_2026-10-27"

S_LIST = [0.5, 1.0, 1.5, 2.0]
N_SA = 120
SEED = 22
T_LIM_C = 90.0
SIG_LIM_MPA = 200.0
MAX_INNER = 5
MAX_OUTER = 3
T_MIN_K, T_MAX_K = 300, 400
T_START_K = 350

# Official ExampleResults at 380 K. HP: only subthreshold moves with T.
OFFICIAL_380 = {
    "Penryn": {"area_mm2": 91.1348, "dyn_W": 47.9315, "sub_W": 10.7217, "gate_W": 0.588125},
    "T2": {"area_mm2": 273.523, "dyn_W": 50.421, "sub_W": 15.0682, "gate_W": 3.50762},
    "Xeon": {"area_mm2": 421.676, "dyn_W": 66.8485, "sub_W": 35.4132, "gate_W": 1.68429},
}
IOFF_HP = {
    "Penryn": {
        0: 2.8e-7, 10: 3.28e-7, 20: 3.81e-7, 30: 4.39e-7, 40: 5.02e-7, 50: 5.69e-7,
        60: 6.42e-7, 70: 7.2e-7, 80: 8.03e-7, 90: 8.91e-7, 100: 9.84e-7,
    },
    "T2": {
        0: 1.96e-7, 10: 2.29e-7, 20: 2.66e-7, 30: 3.05e-7, 40: 3.49e-7, 50: 3.95e-7,
        60: 4.45e-7, 70: 4.97e-7, 80: 5.48e-7, 90: 5.94e-7, 100: 6.3e-7,
    },
}
IOFF_HP["Xeon"] = IOFF_HP["T2"]

# gem5 task ops. Official templates have no gem5 run, so N_op is the task the
# chip stands in for; it does not change with T or s (IPC held fixed).
OPS = {
    "CPU_matmul": {"task": "x86 192x192 matmul", "N_op": 11675900.0, "t_sim_s": 0.0243},
    "CPU_threads": {"task": "x86 threads", "N_op": 1808180.0, "t_sim_s": 0.004124},
    "CPU_armL2": {"task": "ARM Timing+L2", "N_op": 4994.0, "t_sim_s": 3.5e-05},
    "CPU_atomic": {"task": "ARM Atomic", "N_op": 5020.0, "t_sim_s": 3e-06},
    "CPU_starter": {"task": "ARM starter", "N_op": 5020.0, "t_sim_s": 1e-06},
}

GPU_STATIC_W = {"GPU_3070": 22.5847, "GPU_2060S": 24.4467, "GPU_TITANX": 29.3748}
GPU_CONST_W = {"GPU_3070": 32.3252, "GPU_2060S": 32.3252, "GPU_TITANX": 32.3252}
GPU_AREA = {"GPU_3070": 392.0, "GPU_2060S": 445.0, "GPU_TITANX": 471.0}
GPU_POWER = {"GPU_3070": 66.89, "GPU_2060S": 71.12, "GPU_TITANX": 67.41}

# One of each published CIM. All << 1 mm² except jia, so only jia counts as a
# normal chip for the stress gate.
CIM = [
    ("CIM_jia", 4.281, 0.002773),
    ("CIM_colonnade", 0.2314, 0.0412),
    ("CIM_wang", 0.122, 0.026),
    ("CIM_basic", 0.04672, 0.002625),
    ("CIM_sinangil", 0.002577, 0.0007948),
]

# CPU area within ~2x of the GPU dies (392/445/471). Xeon 422 is the match;
# Penryn 91 is 4–5x smaller and stays out.
SCHEMES = [
    {"id": "pkg_2xeon_2gpu", "label": "2×Xeon+2GPU",
     "cpus": ["Xeon", "Xeon"], "gpus": ["GPU_3070", "GPU_TITANX"], "n_cim": 5},
    {"id": "pkg_2xeon_3gpu", "label": "2×Xeon+3GPU",
     "cpus": ["Xeon", "Xeon"], "gpus": ["GPU_3070", "GPU_2060S", "GPU_TITANX"], "n_cim": 5},
    {"id": "pkg_1xeon_1t2_2gpu", "label": "1Xeon+1T2+2GPU",
     "cpus": ["Xeon", "T2"], "gpus": ["GPU_3070", "GPU_TITANX"], "n_cim": 5},
    {"id": "pkg_2t2_2gpu", "label": "2×T2+2GPU",
     "cpus": ["T2", "T2"], "gpus": ["GPU_2060S", "GPU_TITANX"], "n_cim": 5},
    {"id": "pkg_3xeon_3gpu", "label": "3×Xeon+3GPU",
     "cpus": ["Xeon", "Xeon", "Xeon"], "gpus": ["GPU_3070", "GPU_2060S", "GPU_TITANX"], "n_cim": 5},
]

CPU_SLOTS = ["CPU_matmul", "CPU_threads", "CPU_armL2", "CPU_atomic", "CPU_starter"]

HEAT = "#2a78d6"
COUP = "#eb6834"
AQUA = "#1baf7a"
YEL = "#eda100"
VIO = "#7a5af5"
COLORS = [HEAT, AQUA, YEL, COUP, VIO]
SURFACE = "#fcfcfb"
TEXT = "#0b0b0b"
MUTED = "#52514e"
GRID = "#e8e7e3"


def log_md(path: Path, line: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8") as f:
        f.write(line.rstrip() + "\n")
    print(line, flush=True)


def leak_of(model: str, T_K: int) -> dict:
    off = OFFICIAL_380[model]
    tab = IOFF_HP[model]
    if (T_K - 300) not in tab:
        raise ValueError(f"{model} T={T_K} off the 10 K grid")
    ratio = tab[T_K - 300] / tab[80]
    sub = off["sub_W"] * ratio
    return {
        "model": model, "T_K": T_K, "dyn_W": off["dyn_W"], "sub_W": sub,
        "gate_W": off["gate_W"], "leak_W": sub + off["gate_W"], "Ioff_ratio": ratio,
        "area_mm2": off["area_mm2"],
    }


def round10_k(T_C: float) -> int:
    return max(T_MIN_K, min(T_MAX_K, int(round((T_C + 273.15) / 10.0) * 10)))


def lbfgs_clamped(blocks: list[Block], seed: int = SEED):
    """random_lbfgs_layout, but the log-aspect variable is clamped to ±log(4).

    Past ~9 blocks the unbounded variable overflows math.exp. Clamping keeps
    every chip between 1:4 and 4:1, which is all the packer needs.
    """
    import lbfgs as _lbfgs

    real_decode = _lbfgs._decode

    def clamped(blocks, z):
        z = z.copy()
        for i in range(len(blocks)):
            z[3 * i + 2] = max(-1.386, min(1.386, z[3 * i + 2]))
        return real_decode(blocks, z)

    _lbfgs._decode = clamped
    try:
        return random_lbfgs_layout(blocks, seed=seed)
    finally:
        _lbfgs._decode = real_decode


# One HBM3 stack per package. A stack pays the background once per channel and
# then the traffic of every CPU and every GPU in the scheme, measured back to
# back on HBM3_16Gb_8hi at 6400 Mbps. It does not scale with frequency: the
# traces were captured at the chips' nominal rate.
CHAN_BG_W, CHAN_CPU_W, CHAN_GPU_W, N_CHAN = 0.088, 1.035, 0.844, 16
HBM_AREA_MM2 = 11.0 * 11.0


def hbm_power(n_cpu: int, n_gpu: int) -> float:
    per_chan = CHAN_BG_W + n_cpu * (CHAN_CPU_W - CHAN_BG_W) + n_gpu * (CHAN_GPU_W - CHAN_BG_W)
    return per_chan * N_CHAN


def build_chips(spec: dict) -> tuple[list[sa.Chip], dict[str, str]]:
    blocks = []
    cpu_xml: dict[str, str] = {}
    for i, model in enumerate(spec["cpus"]):
        name = CPU_SLOTS[i]
        off = OFFICIAL_380[model]
        blocks.append(Block(name, "CPU", f"McPAT {model}", off["area_mm2"],
                            off["dyn_W"] + off["sub_W"] + off["gate_W"], model))
        cpu_xml[name] = model
    for g in spec["gpus"]:
        blocks.append(Block(g, "GPU", "AccelWattch", GPU_AREA[g], GPU_POWER[g], g))
    cim_list = spec.get("cims", CIM[: spec["n_cim"]])
    for name, area, power in cim_list:
        blocks.append(Block(name, "CIM", "AccelForge", area, power, name))
    if spec.get("hbm", True):
        blocks.append(Block("HBM3_stack", "HBM", "Ramulator2.1", HBM_AREA_MM2,
                            hbm_power(len(spec["cpus"]), len(spec["gpus"])), "HBM3"))
    rects, meta = lbfgs_clamped(blocks, seed=spec.get("seed", SEED))
    chips = [sa.Chip(r.name, r.kind, r.x, r.y, r.w, r.h, r.power_W, r.area_mm2,
                     model=(cpu_xml.get(r.name, "") if r.kind == "CPU" else "")) for r in rects]
    return chips, cpu_xml


def stash(chips: list[sa.Chip]) -> None:
    for c in chips:
        if c.kind == "GPU":
            c._dyn0 = max(0.0, c.power_W - GPU_STATIC_W[c.name] - GPU_CONST_W[c.name])
            c._static = GPU_STATIC_W[c.name]
            c._const = GPU_CONST_W[c.name]
        else:
            c._p0 = c.power_W  # CIM and HBM: fixed, not scaled by frequency


def apply_power(chips, s, cpu_xml, T_K) -> dict:
    cpu_dyn = cpu_leak = gpu_dyn = gpu_leak = gpu_const = cim = 0.0
    per = {}
    for c in chips:
        if c.kind == "CPU":
            m = leak_of(cpu_xml[c.name], T_K)
            dyn, leak = s * m["dyn_W"], m["leak_W"]
            c.power_W = dyn + leak
            cpu_dyn += dyn
            cpu_leak += leak
            op = OPS.get(c.name, OPS["CPU_matmul"])
            per[c.name] = {"model": cpu_xml[c.name], "dyn_W": dyn, "leak_W": leak,
                           "N_op": op["N_op"], "t_sim_s": op["t_sim_s"]}
        elif c.kind == "GPU":
            dyn = s * c._dyn0
            c.power_W = dyn + c._static + c._const
            gpu_dyn += dyn
            gpu_leak += c._static
            gpu_const += c._const
        else:
            c.power_W = c._p0
            cim += c.power_W
    cpu_tot = cpu_dyn + cpu_leak
    n_op = sum(v["N_op"] for v in per.values())
    return {
        "s": s, "McPAT_T_K": T_K,
        "cpu_dyn_W": cpu_dyn, "cpu_leak_W": cpu_leak,
        "cpu_dyn_over_leak": cpu_dyn / cpu_leak if cpu_leak else None,
        "cpu_leak_frac": cpu_leak / cpu_tot if cpu_tot else None,
        "gpu_dyn_W": gpu_dyn, "gpu_leak_W": gpu_leak,
        "pkg_tot_W": cpu_dyn + cpu_leak + gpu_dyn + gpu_leak + gpu_const + cim,
        "N_op": n_op, "per_cpu": per,
    }


def footprint(chips) -> dict:
    m = packing_metrics(sa.to_rects(chips))
    return {"width_mm": m["width_m"] * 1e3, "height_mm": m["height_m"] * 1e3,
            "deadspace": m.get("deadspace")}


# Search grid. Every candidate is priced from cached per-chip fields on the
# 16 grid: steady conduction is linear, so the full field is the sum of the
# single-chip fields minus (N-1) times ambient, and a chip that did not move
# reuses its stored field. The winner is re-scored once on the 64 grid, which
# is the only temperature reported. On the current dies the 16 grid reads
# 0.2-0.5 C low, below the 1-3 C the search accepts, so the walk is unchanged.
SA_SEARCH_ROWS = 16
SA_FINE_ROWS = 64
FIELD_AMB = 318.15
# Big-jump mix (nudge, jump, swap, rotate). Across 19 layouts it finished
# 0.14 C behind the per-layout best on average, against 0.54 C for the default.
SA_MIX = (0.10, 0.60, 0.20, 0.10)
# Stress re-enters the search only inside this band above the limit. Below it
# the heat-only layout already passes; above it no layout has recovered.
SIG_SEARCH_HI_MPA = 230.0


def field_peak(chips, cache: dict, work: Path) -> float:
    """Peak from cached single-chip fields on the 16 grid.

    The package mesh is fixed, so a field depends only on which chip it is and
    where that chip sits. Chips that did not move are summed from the cache;
    only a chip in a new position costs a solve.
    """
    rows = SA_SEARCH_ROWS
    n = rows * rows
    acc = None
    for i, c in enumerate(chips):
        key = (i, round(c.x, 6), round(c.y, 6), round(c.w, 6), round(c.h, 6))
        if key not in cache:
            solo = [sa.Chip(d.name, d.kind, d.x, d.y, d.w, d.h,
                            c.power_W if k == i else 0.0, d.area_mm2)
                    for k, d in enumerate(chips)]
            slot = work / "eval_try"
            if slot.exists():
                shutil.rmtree(slot)
            write_hotspot_inputs(sa.to_rects(solo), slot, rows=rows, cols=rows)
            run_hotspot(sa.HS_BIN, slot)
            vals = []
            for line in (slot / "grid.steady").read_text().splitlines():
                line = line.strip()
                if not line or line.startswith("Layer"):
                    continue
                vals.append(float(line.split()[-1]))
                if len(vals) == n:
                    break
            if len(vals) != n:
                raise RuntimeError(f"16 网格读到 {len(vals)} 格，期望 {n}")
            cache[key] = np.array(vals)
            shutil.rmtree(slot, ignore_errors=True)
        acc = cache[key] if acc is None else acc + cache[key]
    return float((acc - (len(chips) - 1) * FIELD_AMB).max() - 273.15)


def heat_sa(chips0, work: Path, log) -> tuple[list[sa.Chip], float, dict]:
    if work.exists():
        shutil.rmtree(work)
    work.mkdir(parents=True)
    rng = random.Random(SEED)
    sa.shift_to_origin(chips0)
    # Fixed 72 mm chip layer; chips are placed inside it, not the other way round.
    intp_w = intp_h = 0.072
    chips = sa.spread_to_gap(deepcopy(chips0), intp_w, intp_h)
    cache: dict = {}
    tmax0 = field_peak(chips, cache, work)
    current, best = deepcopy(chips), deepcopy(chips)
    tmax_cur = tmax_best = tmax0
    j_cur = j_best = (tmax0 - 45.0) / 25.0
    T = sa.T0
    accepted = 0
    for step in range(1, N_SA + 1):
        trial, op = sa.propose(current, intp_w, intp_h, rng, SA_MIX)
        if not sa.legal(trial, intp_w, intp_h):
            T *= sa.ALPHA
            continue
        try:
            tmax_new = field_peak(trial, cache, work)
        except Exception:
            T *= sa.ALPHA
            continue
        j_new = (tmax_new - 45.0) / 25.0
        delta = j_cur - j_new
        if rng.random() < (1.0 if delta >= 0 else math.exp(delta / max(T, 1e-9))):
            current, tmax_cur, j_cur = trial, tmax_new, j_new
            accepted += 1
            if j_new < j_best - 1e-12:
                best, tmax_best, j_best = deepcopy(trial), tmax_new, j_new
                log(f"    SA step {step:2d} BEST {op} Tmax={tmax_new:.3f} @{SA_SEARCH_ROWS}")
        T *= sa.ALPHA
    # The search price is a 16-grid underestimate. Only the 64-grid re-score is reported.
    tmax_best = sa.eval_hotspot(best, work / "layout_best", rows=SA_FINE_ROWS)
    return best, tmax_best, {"Tmax_init_C": tmax0, "Tmax_best_C": tmax_best,
                             "accepted": accepted, "field_cache": len(cache)}


def real_stress(chips, work: Path, tag: str) -> dict:
    """One real ATSim solve. Call only outside the SA loop."""
    work.mkdir(parents=True, exist_ok=True)
    jp = work / f"{tag}.json"
    write_proxy_json(chips, jp, tag, fill_die=True,
                     package_overrides={"pitch_die_mm": 2.0, "pitch_coarse_mm": 8.0})
    st = eval_atsim_proxy(jp, verbose=False)
    keep = {k: st.get(k) for k in (
        "Tmax_C", "vm_normal_chip_max_MPa", "vm_normal_chip_name", "vm_normal_chip_p95_MPa",
        "vm_center_max_MPa", "w_pv_um", "u_max_um", "ndof", "elapsed_s", "heat_s", "elast_s")}
    (work / f"{tag}_stress.json").write_text(
        json.dumps(keep, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return keep


class Ledger:
    def __init__(self, slot: Path, md: Path):
        self.slot, self.md = slot, md
        self.rows: list[dict] = []
        self.prev: dict | None = None

    def add(self, stage: str, chips, br: dict, tmax: float, stress: dict | None) -> dict:
        fp = footprint(chips)
        row = {
            "step": len(self.rows), "stage": stage,
            "McPAT_T_K": br["McPAT_T_K"], "s": br["s"],
            "Tmax_C": tmax,
            "sigma_MPa": None if stress is None else stress.get("vm_normal_chip_max_MPa"),
            "sigma_chip": None if stress is None else stress.get("vm_normal_chip_name"),
            "w_pv_um": None if stress is None else stress.get("w_pv_um"),
            "cpu_leak_W": br["cpu_leak_W"], "cpu_dyn_W": br["cpu_dyn_W"],
            "cpu_dyn_over_leak": br["cpu_dyn_over_leak"],
            "cpu_leak_frac": br["cpu_leak_frac"],
            "pkg_tot_W": br["pkg_tot_W"], "N_op": br["N_op"],
            "outline_mm": f"{fp['width_mm']:.1f}x{fp['height_mm']:.1f}",
        }
        if self.prev:
            row["dTmax_C"] = tmax - self.prev["Tmax_C"]
            row["dLeak_W"] = br["cpu_leak_W"] - self.prev["cpu_leak_W"]
            row["dPkg_W"] = br["pkg_tot_W"] - self.prev["pkg_tot_W"]
            if stress and self.prev.get("sigma_MPa") is not None:
                row["dSigma_MPa"] = stress["vm_normal_chip_max_MPa"] - self.prev["sigma_MPa"]
        self.rows.append(row)
        self.prev = row
        with (self.slot / "steps.jsonl").open("a", encoding="utf-8") as f:
            f.write(json.dumps(row, ensure_ascii=False) + "\n")
        sig = "—" if row["sigma_MPa"] is None else f"{row['sigma_MPa']:.1f}"
        dT = "" if "dTmax_C" not in row else f"  ΔT {row['dTmax_C']:+.2f}  Δleak {row['dLeak_W']:+.2f} W"
        log_md(self.md,
               f"| {row['step']} | {stage} | {br['McPAT_T_K']} | {br['s']:g} | {tmax:.2f} | {sig} | "
               f"{br['cpu_leak_W']:.2f} | {br['cpu_dyn_over_leak']:.3f} | {br['cpu_leak_frac']:.3f} | "
               f"{br['pkg_tot_W']:.1f} | {br['N_op']:.3g} | {row['outline_mm']} |{dT}")
        return row


def close_loop(chips, s, cpu_xml, T0, slot, md, ledger, stage_prefix, with_sa: bool) -> tuple:
    """HotSpot ↔ McPAT bin, optional heat-only SA, until the bin is stable."""
    T, br = T0, {}
    for outer in range(MAX_OUTER):
        for i in range(MAX_INNER):
            br = apply_power(chips, s, cpu_xml, T)
            tmax = sa.eval_hotspot(chips, slot / f"{stage_prefix}_o{outer}_i{i}")
            nxt = round10_k(tmax)
            ledger.add(f"{stage_prefix} 内环{i}", chips, br, tmax, None)
            if nxt == T:
                break
            T = nxt
        if not with_sa:
            return chips, T, br, tmax
        log_md(md, f"- {stage_prefix} 外环{outer}: 只热 SA {N_SA} 步，16 网格场叠加，起点为当前布局")
        chips, tmax_sa, info = heat_sa(chips, slot / f"{stage_prefix}_sa{outer}",
                                       lambda m: log_md(md, m))
        br = apply_power(chips, s, cpu_xml, T)
        ledger.add(f"{stage_prefix} SA{outer}", chips, br, tmax_sa, None)
        log_md(md, f"    SA {info['Tmax_init_C']:.2f}→{tmax_sa:.2f} °C  accepted {info['accepted']}/{N_SA}")
        nxt = round10_k(tmax_sa)
        if nxt == T:
            log_md(md, f"- **收敛** {stage_prefix}：温度档 {T} K 与布局一致")
            return chips, T, br, tmax_sa
        T = nxt
    log_md(md, f"- {stage_prefix} 外环用尽，停在 {T} K")
    return chips, T, br, tmax_sa


def theta(s: float) -> float:
    return (s + s + 1.0) / 3.0


def run_scheme(spec, color, md) -> dict:
    sid = spec["id"]
    slot = OUT / sid
    slot.mkdir(parents=True, exist_ok=True)
    ledger = Ledger(slot, md)
    log_md(md, f"## {spec['label']}  (`{sid}`)")
    log_md(md, "")
    log_md(md, "| 步 | 阶段 | 档K | s | Tmax | 应力MPa | CPU漏电W | dyn/leak | 漏电占比 | 封装W | 操作数 | 轮廓mm | 相对上步 |")
    log_md(md, "|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|---|")

    chips, cpu_xml = build_chips(spec)
    stash(chips)
    (slot / "lbfgs_init.json").write_text(
        json.dumps(sa.snapshot(chips), indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    cpu_area = sum(OFFICIAL_380[m]["area_mm2"] for m in spec["cpus"])
    gpu_area = sum(GPU_AREA[g] for g in spec["gpus"])
    log_md(md, "")
    log_md(md, f"- CPU {cpu_area:.0f} mm² / GPU {gpu_area:.0f} mm²，比值 {cpu_area/gpu_area:.2f}")

    # A: architecture feedback from 350 K, SA inside, until the bin is stable.
    log_md(md, "")
    log_md(md, "### A 架构回馈（起点 350 K，只热 SA）")
    chips, T_star, br, tmax = close_loop(
        chips, 1.0, cpu_xml, T_START_K, slot / "A", md, ledger, "A", with_sa=True)
    (slot / "Lstar.json").write_text(
        json.dumps(sa.snapshot(chips), indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    Lstar = deepcopy(chips)

    # B: real stress once, on the converged layout.
    log_md(md, "")
    log_md(md, "### B 收敛布局上的真 ATSim")
    st = real_stress(chips, slot / "B", "Lstar")
    ledger.add("B 真ATSim", chips, br, st["Tmax_C"], st)
    log_md(md, f"- 应力 {st['vm_normal_chip_max_MPa']:.1f} MPa @ {st['vm_normal_chip_name']}，"
               f"翘曲 {st['w_pv_um']:.2f} µm，ndof {st['ndof']}，用时 {st['elapsed_s']:.0f} s")

    # Stress only re-enters the search when the heat-optimal layout lands in the
    # band where a layout change could plausibly pull it back. Far below the
    # limit there is nothing to gain; far above it, no layout has ever recovered.
    sig = st["vm_normal_chip_max_MPa"] or 0.0
    if SIG_LIM_MPA <= sig <= SIG_SEARCH_HI_MPA:
        log_md(md, "")
        log_md(md, f"### B2 应力 {sig:.1f} MPa 落在 {SIG_LIM_MPA:.0f}–{SIG_SEARCH_HI_MPA:.0f}，"
                   f"把应力加进 SA 重跑")
        # The search prices temperature only. The real stress of the layout
        # this produces is computed once, just below.
        chips, T_star, br, tmax = close_loop(
            chips, 1.0, cpu_xml, T_star, slot / "B2", md, ledger, "B2",
            with_sa=True)
        (slot / "Lstar.json").write_text(
            json.dumps(sa.snapshot(chips), indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
        Lstar = deepcopy(chips)
        st = real_stress(chips, slot / "B2", "Lstar")
        ledger.add("B2 真ATSim", chips, br, st["Tmax_C"], st)
        log_md(md, f"- 重跑后应力 {st['vm_normal_chip_max_MPa']:.1f} MPa @ {st['vm_normal_chip_name']}")
    elif sig > SIG_SEARCH_HI_MPA:
        log_md(md, f"- 应力 {sig:.1f} MPa 超过 {SIG_SEARCH_HI_MPA:.0f}，布局搜索拉不回来，跳过应力重跑")

    # C: frequency sweep, real ATSim each point, layout frozen at L*.
    log_md(md, "")
    log_md(md, "### C 选频（L* 冻住，每档真 ATSim）")
    freq = []
    for s in S_LIST:
        T = T_star
        cur = deepcopy(Lstar)
        for i in range(MAX_INNER):
            b = apply_power(cur, s, cpu_xml, T)
            t_hs = sa.eval_hotspot(cur, slot / f"C_s{s}_i{i}")
            nxt = round10_k(t_hs)
            if nxt == T:
                break
            T = nxt
        st_s = real_stress(cur, slot / "C", f"s{s}")
        ok = st_s["Tmax_C"] < T_LIM_C and (st_s["vm_normal_chip_max_MPa"] or 1e9) < SIG_LIM_MPA
        rec = {**b, "Tmax_C": st_s["Tmax_C"], "sigma_MPa": st_s["vm_normal_chip_max_MPa"],
               "sigma_chip": st_s["vm_normal_chip_name"], "w_pv_um": st_s["w_pv_um"],
               "Theta_pkg": theta(s), "in_window": ok}
        freq.append(rec)
        ledger.add(f"C s={s:g}", cur, b, st_s["Tmax_C"], st_s)
        log_md(md, f"- s={s:g} 档 {T} K  Tmax {st_s['Tmax_C']:.2f}  "
                   f"应力 {st_s['vm_normal_chip_max_MPa']:.1f}  Θ {theta(s):.3f}  "
                   f"{'窗内' if ok else '窗外'}")
    feas = [r for r in freq if r["in_window"]]
    winners = {}
    if feas:
        winners["leak"] = min(feas, key=lambda r: r["cpu_leak_frac"])
        winners["theta"] = max(feas, key=lambda r: r["Theta_pkg"])
        log_md(md, f"- 漏电占比赢家 s={winners['leak']['s']:g}（frac {winners['leak']['cpu_leak_frac']:.3f}）")
        log_md(md, f"- 吞吐赢家 s={winners['theta']['s']:g}（Θ {winners['theta']['Theta_pkg']:.3f}）")
    else:
        log_md(md, "- 四个 s 都不在窗内，这一套没有赢家，不做 D。")

    # D: re-SA from L* at the winning s, then close the bin, real ATSim at the end.
    finals = {}
    for rule, w in winners.items():
        if rule == "theta" and winners.get("leak") and w["s"] == winners["leak"]["s"]:
            log_md(md, f"- 吞吐赢家与漏电赢家同为 s={w['s']:g}，D 不重复跑")
            finals[rule] = finals["leak"]
            continue
        log_md(md, "")
        log_md(md, f"### D {rule} 赢家 s={w['s']:g}（初始坐标 = L*，不是 L-BFGS）")
        cur = deepcopy(Lstar)
        cur, T_f, br_f, tmax_f = close_loop(
            cur, w["s"], cpu_xml, w["McPAT_T_K"], slot / f"D_{rule}", md, ledger,
            f"D-{rule}", with_sa=True)
        st_f = real_stress(cur, slot / f"D_{rule}", "final")
        br_f = apply_power(cur, w["s"], cpu_xml, T_f)
        ledger.add(f"D {rule} 最终真ATSim", cur, br_f, st_f["Tmax_C"], st_f)
        ok = st_f["Tmax_C"] < T_LIM_C and (st_f["vm_normal_chip_max_MPa"] or 1e9) < SIG_LIM_MPA
        finals[rule] = {
            "s": w["s"], "McPAT_T_K": T_f, "Tmax_C": st_f["Tmax_C"],
            "sigma_MPa": st_f["vm_normal_chip_max_MPa"], "sigma_chip": st_f["vm_normal_chip_name"],
            "w_pv_um": st_f["w_pv_um"], "cpu_leak_W": br_f["cpu_leak_W"],
            "cpu_dyn_over_leak": br_f["cpu_dyn_over_leak"], "cpu_leak_frac": br_f["cpu_leak_frac"],
            "pkg_tot_W": br_f["pkg_tot_W"], "N_op": br_f["N_op"], "Theta_pkg": theta(w["s"]),
            "in_window": ok,
        }
        (slot / f"D_{rule}_layout.json").write_text(
            json.dumps(sa.snapshot(cur), indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
        log_md(md, f"- 最终 s={w['s']:g} 档 {T_f} K  Tmax {st_f['Tmax_C']:.2f}  "
                   f"应力 {st_f['vm_normal_chip_max_MPa']:.1f}  漏电 {br_f['cpu_leak_W']:.2f} W  "
                   f"{'仍在窗内' if ok else 'SA 后出窗'}")

    (slot / "ledger.json").write_text(
        json.dumps(ledger.rows, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    log_md(md, "")
    return {"id": sid, "label": spec["label"], "color": color,
            "cpu_area_mm2": cpu_area, "gpu_area_mm2": gpu_area,
            "T_star_K": T_star, "freq": freq, "winners": {k: v["s"] for k, v in winners.items()},
            "finals": finals, "steps": ledger.rows}


def main() -> None:
    if OUT.exists():
        shutil.rmtree(OUT)
    OUT.mkdir(parents=True)
    md = OUT / "EXPERIMENT.md"
    log_md(md, "# 全链：350 K 架构回馈 → 真 ATSim 选频 → 从 L* 再 SA → 再收敛")
    log_md(md, "")
    log_md(md, f"- started (UTC): {datetime.now(timezone.utc).isoformat()}")
    log_md(md, "- 只用面积和 GPU 同级的 CPU：Xeon 422 mm²、T2 274 mm²。Penryn 91 mm² 和小核 ~14 mm² 不进这轮。")
    log_md(md, "- 漏电不是现场跑 mcpat。官方模板 HP，温度只改 I_off；每档 = 官方 380 K 亚阈值 × ITRS 表比值，栅漏电固定。350 K 在这张表上，所以起点就是 350 K。")
    log_md(md, f"- 应力是真 ATSim（solve_fvm3d + solve_elasticity），只在 B、每个 s、以及 D 的最终布局上跑。SA {N_SA} 步内只用 HotSpot，避免 12 GB 内存里叠求解器。")
    log_md(md, "- 选频两条都算：窗内（Tmax<90、应力<200）漏电占比最低，和吞吐最高。D 的初始坐标是 A 收敛的 L*，不是 L-BFGS。")
    log_md(md, "- 操作数来自 gem5 任务记录，不随温度和频率变（IPC 固定）。")
    log_md(md, "")

    results = []
    for spec, color in zip(SCHEMES, COLORS):
        results.append(run_scheme(spec, color, md))

    log_md(md, "## 总表")
    log_md(md, "")
    log_md(md, "| 方案 | CPU/GPU面积 | 收敛档 | 判据 | s | 最终Tmax | 应力 | 漏电W | 漏电占比 | Θ | 窗 |")
    log_md(md, "|---|---|---:|---|---:|---:|---:|---:|---:|---:|---|")
    for R in results:
        if not R["finals"]:
            log_md(md, f"| {R['label']} | {R['cpu_area_mm2']:.0f}/{R['gpu_area_mm2']:.0f} | "
                       f"{R['T_star_K']} K | — | — | — | — | — | — | — | 无窗内点 |")
            continue
        for rule, f in R["finals"].items():
            name = "漏电占比" if rule == "leak" else "吞吐"
            log_md(md, f"| {R['label']} | {R['cpu_area_mm2']:.0f}/{R['gpu_area_mm2']:.0f} | {R['T_star_K']} K | "
                       f"{name} | {f['s']:g} | {f['Tmax_C']:.2f} | {f['sigma_MPa']:.1f} | "
                       f"{f['cpu_leak_W']:.2f} | {f['cpu_leak_frac']:.3f} | {f['Theta_pkg']:.3f} | "
                       f"{'是' if f['in_window'] else '否'} |")
    log_md(md, "")

    fig, axes = plt.subplots(1, 2, figsize=(12.4, 4.8))
    fig.set_facecolor(SURFACE)
    for R in results:
        xs = [r["s"] for r in R["freq"]]
        axes[0].plot(xs, [r["Tmax_C"] for r in R["freq"]], "-o", color=R["color"], label=R["label"], lw=2, ms=6)
        axes[1].plot(xs, [r["sigma_MPa"] for r in R["freq"]], "-o", color=R["color"], label=R["label"], lw=2, ms=6)
    axes[0].axhline(T_LIM_C, color="#e34948", ls="--", lw=1)
    axes[1].axhline(SIG_LIM_MPA, color="#e34948", ls="--", lw=1)
    axes[0].set_title("L* 上选频：Tmax", loc="left", color=TEXT, fontsize=12)
    axes[1].set_title("L* 上选频：真 ATSim 应力", loc="left", color=TEXT, fontsize=12)
    for ax, ylab in ((axes[0], "Tmax (°C)"), (axes[1], "最大芯粒应力 (MPa)")):
        ax.set_xlabel("s = f / f0", color=MUTED)
        ax.set_ylabel(ylab, color=MUTED)
        ax.set_facecolor(SURFACE)
        ax.yaxis.grid(True, color=GRID, lw=1)
        ax.set_axisbelow(True)
        ax.tick_params(colors=MUTED, length=0)
        ax.legend(frameon=False, fontsize=8)
        for sp in ax.spines.values():
            sp.set_visible(False)
        ax.spines["bottom"].set_visible(True)
        ax.spines["bottom"].set_color(GRID)
    fig.tight_layout()
    fig.savefig(OUT / "chain_freq.png", dpi=150, bbox_inches="tight", facecolor=SURFACE)
    plt.close(fig)

    (OUT / "result.json").write_text(json.dumps({
        "started_utc": datetime.now(timezone.utc).isoformat(),
        "T_start_K": T_START_K, "s_list": S_LIST, "n_sa": N_SA,
        "results": results,
    }, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    if BACKUP.exists():
        shutil.rmtree(BACKUP)
    shutil.copytree(OUT, BACKUP, ignore=shutil.ignore_patterns("eval_*", "A", "C", "B", "D_leak", "D_theta"))
    # keep the ledgers and final layouts, drop the per-iteration HotSpot grids
    log_md(md, f"日志：`{OUT}`")
    print("WROTE", OUT)


if __name__ == "__main__":
    main()
