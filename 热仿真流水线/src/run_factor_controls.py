#!/usr/bin/env python3
"""七个现存问题各自做控制变量，不改生产退火，也不改默认填充。

绝对温度分三套，互相当成对方：
- python：没有均热片外伸
- HotSpot：只在注明的点上算
- ATSim：phase=stress，对流热阻 0.1，芯片内部 von Mises

功耗模型是现成表，不是把同一个任务重跑 gem5。
"""

from __future__ import annotations

import json
import sys
from copy import deepcopy
from pathlib import Path

import numpy as np

ROOT = Path("/home/czh/任务仿真产生功耗文件工具")
SRC = ROOT / "热仿真流水线" / "src"
sys.path.insert(0, str(SRC))

import run_fill_superposition as hs  # noqa: E402
import run_hw_area_field as hw  # noqa: E402
import run_sa_hw_smooth as sm  # noqa: E402

OUT = ROOT / "热仿真流水线" / "runs" / "factor_controls"
T_AMB = hs.T_AMB
N = hs.N
CW_MM = sm.CW_MM
_FACTOR = {"key": None, "fac": None, "g": None}

# 只抄温度表，不导入 run_full_chain（那个会拉起 fenics）。
OFFICIAL = {
    "Penryn": {"dyn_W": 47.9315, "sub_W": 10.7217, "gate_W": 0.588125, "cores": 2, "f_MHz": 3700},
    "T2": {"dyn_W": 50.421, "sub_W": 15.0682, "gate_W": 3.50762, "cores": 8, "f_MHz": 1400},
    "Xeon": {"dyn_W": 66.8485, "sub_W": 35.4132, "gate_W": 1.68429, "cores": 2, "f_MHz": 3400},
}
IOFF = {
    "Penryn": {40: 5.02e-7, 80: 8.03e-7},
    "T2": {40: 3.49e-7, 80: 5.48e-7},
}
IOFF["Xeon"] = IOFF["T2"]


def leak_W(model: str, T_K: int = 340) -> float:
    off = OFFICIAL[model]
    ratio = IOFF[model][T_K - 300] / IOFF[model][80]
    return off["sub_W"] * ratio + off["gate_W"]


def layers(k_tim: float, k_spread: float):
    return (
        ("chip", 0.00015, None),
        ("tim", 2.0e-05, k_tim),
        ("spreader", 0.001, k_spread),
        ("sink", 0.0069, 400.0),
    )


def solve(k_chip: np.ndarray, power: np.ndarray, k_tim: float, k_spread: float):
    # 只留最近一次分解。冻结场要连解四次，材料图相同，不能各存一份因子。
    key = (float(k_tim), float(k_spread), k_chip.tobytes())
    old = hs.LAYERS
    hs.LAYERS = layers(k_tim, k_spread)
    try:
        if _FACTOR["key"] != key:
            mat, g_amb = hs.build_operator(k_chip)
            from scipy.sparse.linalg import splu
            _FACTOR["key"] = key
            _FACTOR["fac"] = splu(mat.tocsc())
            _FACTOR["g"] = g_amb
        state = hs.solve_state(k_chip, power, factor=_FACTOR["fac"], g_amb=_FACTOR["g"])
    finally:
        hs.LAYERS = old
    return state


def paint(dies, k_fill: float) -> np.ndarray:
    return hw.paint(dies, k_fill)


def zero_but(dies, idx: int) -> list[dict]:
    out = deepcopy(dies)
    for i, d in enumerate(out):
        if i != idx:
            d["power_W"] = 0.0
    return out


def truth(dies, k_fill, k_tim, k_spread):
    return solve(paint(dies, k_fill), hw.power_map(dies), k_tim, k_spread)


def naive_sum(dies, k_fill, k_tim, k_spread) -> np.ndarray:
    acc = np.zeros((N, N))
    for i, d in enumerate(dies):
        one = [d]
        st = solve(paint(one, k_fill), hw.power_map(one), k_tim, k_spread)
        acc += st[0] - T_AMB
    return acc + T_AMB


def frozen_parts(dies, k_fill, k_tim, k_spread) -> list[np.ndarray]:
    k = paint(dies, k_fill)
    parts = []
    for i in range(len(dies)):
        st = solve(k, hw.power_map(zero_but(dies, i)), k_tim, k_spread)
        parts.append(st[0])
    return parts


def sum_parts(parts: list[np.ndarray]) -> np.ndarray:
    acc = np.zeros((N, N))
    for p in parts:
        acc += p - T_AMB
    return acc + T_AMB


def err(a: np.ndarray, b: np.ndarray) -> dict:
    d = a - b
    return {"max_abs_K": float(np.max(np.abs(d))), "rmse_K": float(np.sqrt(np.mean(d ** 2)))}


def tmax(T: np.ndarray) -> float:
    return float(T.max() - 273.15)


def cells_of(dies) -> set:
    s = set()
    for d in dies:
        s.update(hw.cells(d))
    return s


def lateral_fraction(state, k_chip, die) -> dict:
    """CPU 足迹净流出的横向热 / 这颗芯片功耗。向上的热单独记。"""
    t_chip = 1.5e-4
    cc = set(hw.cells(die))
    q_up = 0.0
    q_lat = 0.0
    chip, tim = state[0], state[1]
    for i, j in cc:
        g_up = 1.0 / hs._rz(float(k_chip[i, j]), t_chip)
        q_up += g_up * (chip[i, j] - tim[i, j])
        for di, dj in ((0, 1), (1, 0), (0, -1), (-1, 0)):
            ni, nj = i + di, j + dj
            if not (0 <= ni < N and 0 <= nj < N) or (ni, nj) in cc:
                continue
            r = 0.5 * hs._rx(float(k_chip[i, j]), t_chip) + 0.5 * hs._rx(float(k_chip[ni, nj]), t_chip)
            q_lat += (chip[i, j] - chip[ni, nj]) / r
    p = float(die["power_W"])
    return {"Q_up_W": q_up, "Q_lat_W": q_lat, "lat_frac": q_lat / p if p else None}


def base_dies() -> tuple[list[dict], dict]:
    cat = sm.warehouse()
    dies = sm.start_from(cat, "Penryn", "GPU_3070") if hasattr(sm, "start_from") else None
    if dies is None:
        cpu = next(c for c in cat["CPU"] if c["name"] == "Penryn")
        gpu = next(g for g in cat["GPU"] if g["name"] == "GPU_3070")
        dies = [sm.clonable(cpu), sm.clonable(gpu)] + [sm.clonable(f) for f in cat["fixed"]]
    return dies, cat


def replace_cpu(dies, cat, name: str) -> list[dict]:
    part = next(c for c in cat["CPU"] if c["name"] == name)
    out = deepcopy(dies)
    j0, ib0 = out[0]["j0"], out[0]["ib0"]
    out[0] = sm.clonable(part)
    out[0]["j0"], out[0]["ib0"] = j0, ib0
    return out


def save(report: dict) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "report.json").write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def pack(dies) -> list[dict]:
    return [sm.dump_die(d) for d in dies]


def exp_spread_and_fill(dies) -> dict:
    """问题 1 和 7：谁在横向导热。填充导热和均热片导热分开拧。"""
    cpu = dies[0]
    spread_rows = []
    for k_sp in (400, 200, 80, 40, 20, 10, 4, 2, 0.8):
        st = truth(dies, 0.8, 4.0, k_sp)
        fr = sum_parts(frozen_parts(dies, 0.8, 4.0, k_sp))
        nav = naive_sum(dies, 0.8, 4.0, k_sp)
        flow = lateral_fraction(st, paint(dies, 0.8), cpu)
        spread_rows.append({
            "k_spreader": k_sp, "k_tim": 4.0, "k_fill": 0.8,
            "python_Tmax_C": tmax(st[0]),
            "naive": err(nav, st[0]), "frozen": err(fr, st[0]),
            **flow,
        })
        print(f"  spreader {k_sp:6.1f}  naive {spread_rows[-1]['naive']['max_abs_K']:.4f} K  "
              f"lat {flow['lat_frac']*100:.3f}%  T {spread_rows[-1]['python_Tmax_C']:.1f}", flush=True)
    tim_rows = []
    for k_tim in (4.0, 1.0, 0.4, 0.1):
        st = truth(dies, 0.8, k_tim, 400.0)
        nav = naive_sum(dies, 0.8, k_tim, 400.0)
        flow = lateral_fraction(st, paint(dies, 0.8), cpu)
        tim_rows.append({
            "k_spreader": 400.0, "k_tim": k_tim, "k_fill": 0.8,
            "python_Tmax_C": tmax(st[0]),
            "naive_max_abs_K": err(nav, st[0])["max_abs_K"],
            **flow,
        })
        print(f"  tim {k_tim:4.1f}  naive {tim_rows[-1]['naive_max_abs_K']:.4f} K  "
              f"lat {flow['lat_frac']*100:.3f}%", flush=True)
    grid = []
    for k_sp in (400.0, 40.0, 4.0):
        for k_fill in (130.0, 1.0, 0.8, 0.028):
            st = truth(dies, k_fill, 4.0, k_sp)
            nav = naive_sum(dies, k_fill, 4.0, k_sp)
            grid.append({
                "k_spreader": k_sp, "k_fill": k_fill,
                "python_Tmax_C": tmax(st[0]),
                "naive_max_abs_K": err(nav, st[0])["max_abs_K"],
            })
    return {"spreader_sweep": spread_rows, "tim_sweep": tim_rows, "fill_by_spreader": grid}


def make_moves(dies, cat) -> list[tuple[str, list[dict] | None]]:
    moves = []
    for step, name in ((1, "nudge_1"), (2, "nudge_2"), (4, "nudge_4")):
        trial = deepcopy(dies)
        trial[0]["j0"] += step
        moves.append((name, trial if sm.legal(trial) else None))
    rot = deepcopy(dies)
    rot[0]["w"], rot[0]["h"] = rot[0]["h"], rot[0]["w"]
    moves.append(("rotate", rot if sm.legal(rot) else None))
    jump = None
    for j0, ib0 in ((40, 40), (48, 4), (40, 4), (4, 40)):
        trial = deepcopy(dies)
        trial[0]["j0"], trial[0]["ib0"] = j0, ib0
        if sm.legal(trial):
            jump = trial
            break
    moves.append(("jump", jump))
    swap = deepcopy(dies)
    swap[0]["j0"], swap[1]["j0"] = swap[1]["j0"], swap[0]["j0"]
    swap[0]["ib0"], swap[1]["ib0"] = swap[1]["ib0"], swap[0]["ib0"]
    moves.append(("swap_cpu_gpu", swap if sm.legal(swap) else None))
    moves.append(("hw_T2", replace_cpu(dies, cat, "T2")))
    moves.append(("hw_Xeon", replace_cpu(dies, cat, "Xeon")))
    return moves


def exp_moves(dies, cat) -> list[dict]:
    """问题 2：动作改了多少格材料，忽略材料图的误差有多大。两个均热片导热。"""
    rows = []
    for k_sp in (400.0, 4.0):
        old_parts = frozen_parts(dies, 0.8, 4.0, k_sp)
        old_sum = sum_parts(old_parts)
        st_old = truth(dies, 0.8, 4.0, k_sp)
        # 冻结场相加应当就是真解。
        frozen_check = err(old_sum, st_old[0])["max_abs_K"]
        print(f"  k_sp {k_sp} frozen-check {frozen_check:.2e} K", flush=True)
        for name, trial in make_moves(dies, cat):
            if trial is None:
                rows.append({"move": name, "k_spreader": k_sp, "legal": False})
                print(f"  {name} illegal", flush=True)
                continue
            changed = len(cells_of(dies) ^ cells_of(trial))
            st_new = truth(trial, 0.8, 4.0, k_sp)
            # 只重算被改的那颗，其余沿用旧材料图上的冻结场。
            moved = [i for i, (a, b) in enumerate(zip(dies, trial))
                     if a["name"] != b["name"] or a["j0"] != b["j0"] or a["ib0"] != b["ib0"]
                     or a["w"] != b["w"] or a["h"] != b["h"]]
            k_new = paint(trial, 0.8)
            parts = list(old_parts)
            for i in moved:
                parts[i] = solve(k_new, hw.power_map(zero_but(trial, i)), 4.0, k_sp)[0]
            stale = sum_parts(parts)
            # 材料图不更新，只把功耗换成新的。
            ignore = solve(paint(dies, 0.8), hw.power_map(trial), 4.0, k_sp)[0]
            nav = naive_sum(trial, 0.8, 4.0, k_sp)
            rec = {
                "move": name, "k_spreader": k_sp, "legal": True,
                "cells_changed": changed,
                "hardware": [d["name"] for d in trial],
                "field_moved_K": err(st_new[0], st_old[0])["max_abs_K"],
                "err_stale_cache_K": err(stale, st_new[0])["max_abs_K"],
                "err_ignore_k_K": err(ignore, st_new[0])["max_abs_K"],
                "err_naive_K": err(nav, st_new[0])["max_abs_K"],
                "python_Tmax_C": tmax(st_new[0]),
            }
            rows.append(rec)
            print(f"  {name:12} k {k_sp:5.0f}  cells {changed:4d}  "
                  f"stale {rec['err_stale_cache_K']:.4f}  ignore_k {rec['err_ignore_k_K']:.4f}  "
                  f"naive {rec['err_naive_K']:.4f}", flush=True)
    return rows


def exp_ladder(dies, cat) -> list[dict]:
    """问题 3：同一组布局，python 和 HotSpot 的升降是否同号。"""
    rows = []
    layouts = [("start", dies)]
    for step in (1, 2, 4):
        trial = deepcopy(dies)
        trial[0]["j0"] += step
        if sm.legal(trial):
            layouts.append((f"nudge_{step}", trial))
    layouts.append(("T2", replace_cpu(dies, cat, "T2")))
    layouts.append(("Xeon", replace_cpu(dies, cat, "Xeon")))
    for name, lay in layouts:
        st = truth(lay, 130.0, 4.0, 400.0)
        hs_c = sm.hotspot(lay, OUT / "hotspot" / f"ladder_{name}")
        rec = {
            "step": name,
            "python_Tmax_C": tmax(st[0]),
            "grad_K_per_mm": sm.grad_max(st[0]),
            "hotspot_Tmax_C": hs_c,
            "dies": pack(lay),
            "stress_tag": f"ladder_{name}",
        }
        rows.append(rec)
        print(f"  {name:8} python {rec['python_Tmax_C']:.2f}  HotSpot {hs_c:.2f}  "
              f"grad {rec['grad_K_per_mm']:.3f}", flush=True)
    return rows


def exp_frequency(dies, cat) -> list[dict]:
    """问题 4：几何不动，只改功耗模型。延迟是公式，不是 gem5。"""
    rows = []
    gpu_keep = [deepcopy(d) for d in dies[1:]]
    for cpu_name in ("Penryn", "T2", "Xeon"):
        lay0 = replace_cpu(dies, cat, cpu_name)
        meta = OFFICIAL[cpu_name]
        for s in (0.5, 1.0, 1.5, 2.0):
            lay = deepcopy(lay0)
            p = s * meta["dyn_W"] + leak_W(cpu_name, 340)
            lay[0]["power_W"] = p
            st = truth(lay, 130.0, 4.0, 400.0)
            hs_c = sm.hotspot(lay, OUT / "hotspot" / f"freq_{cpu_name}_s{s:.1f}")
            th_s = meta["f_MHz"] * s / 3700.0
            th_p = meta["cores"] * meta["f_MHz"] * s / 3700.0
            rec = {
                "cpu": cpu_name, "s": s, "power_W": p,
                "python_Tmax_C": tmax(st[0]), "hotspot_Tmax_C": hs_c,
                "grad_K_per_mm": sm.grad_max(st[0]),
                "Theta_serial": th_s, "Theta_parallel": th_p,
                "leak_bin_K": 340,
            }
            rows.append(rec)
            print(f"  {cpu_name:7} s={s:.1f}  P {p:.1f}  HotSpot {hs_c:.2f}  "
                  f"Θser {th_s:.2f} Θpar {th_p:.2f}", flush=True)
        # 仓库里的 380 K 总功耗，当作不随频率变的查表点。
        lay = deepcopy(lay0)
        st = truth(lay, 130.0, 4.0, 400.0)
        hs_c = sm.hotspot(lay, OUT / "hotspot" / f"freq_{cpu_name}_template")
        rows.append({
            "cpu": cpu_name, "s": None, "power_W": lay[0]["power_W"],
            "python_Tmax_C": tmax(st[0]), "hotspot_Tmax_C": hs_c,
            "grad_K_per_mm": sm.grad_max(st[0]),
            "Theta_serial": meta["f_MHz"] / 3700.0,
            "Theta_parallel": meta["cores"] * meta["f_MHz"] / 3700.0,
            "leak_bin_K": 380, "template": True,
        })
    return rows


def exp_power_modes(dies, cat) -> list[dict]:
    """问题 5：面积、总功耗、功率密度拆开。应力在 phase=stress。"""
    pen = replace_cpu(dies, cat, "Penryn")
    dens = pen[0]["power_W"] / pen[0]["grid_area_mm2"]
    p_fixed = pen[0]["power_W"]
    rows = []
    for cpu_name in ("Penryn", "T2", "Xeon"):
        native = replace_cpu(dies, cat, cpu_name)
        modes = {
            "native": native[0]["power_W"],
            "same_power": p_fixed,
            "same_density": dens * native[0]["grid_area_mm2"],
        }
        for mode, power in modes.items():
            lay = deepcopy(native)
            lay[0]["power_W"] = power
            st = truth(lay, 130.0, 4.0, 400.0)
            tag = f"mode_{cpu_name}_{mode}"
            (OUT / "layouts").mkdir(exist_ok=True)
            (OUT / "layouts" / f"{tag}.json").write_text(json.dumps(pack(lay), indent=2) + "\n", encoding="utf-8")
            rec = {
                "cpu": cpu_name, "mode": mode, "tag": tag,
                "power_W": power,
                "grid_area_mm2": lay[0]["grid_area_mm2"],
                "density_W_mm2": power / lay[0]["grid_area_mm2"],
                "python_Tmax_C": tmax(st[0]),
                "grad_K_per_mm": sm.grad_max(st[0]),
                "hotspot_Tmax_C": sm.hotspot(lay, OUT / "hotspot" / tag),
            }
            rows.append(rec)
            print(f"  {cpu_name:7} {mode:12} dens {rec['density_W_mm2']:.3f}  "
                  f"grad {rec['grad_K_per_mm']:.3f}  HotSpot {rec['hotspot_Tmax_C']:.2f}", flush=True)
    return rows


def exp_geometry(cat) -> list[dict]:
    """问题 7：面积变大时，和 CIM 的接触边长并不变。"""
    cim = next(d for d in cat["fixed"] if d["role"] == "CIM")
    gpu = next(g for g in cat["GPU"] if g["name"] == "GPU_3070")
    rows = []
    for cpu in cat["CPU"]:
        gap = (cim["ib0"] - (cpu["ib0"] + cpu["h"])) * CW_MM
        face = min(cpu["w"], cim["w"]) * CW_MM
        gap_g = (gpu["j0"] - (cpu["j0"] + cpu["w"])) * CW_MM
        face_g = min(cpu["w"], gpu["w"]) * CW_MM
        rows.append({
            "cpu": cpu["name"],
            "grid_area_mm2": cpu["grid_area_mm2"],
            "gap_to_cim_mm": gap,
            "facing_cim_mm": face,
            "gap_to_gpu_mm": gap_g,
            "facing_gpu_mm": face_g,
            "cim_w_mm": cim["w"] * CW_MM,
            "cpu_w_mm": cpu["w"] * CW_MM,
        })
    return rows


def thermal() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "layouts").mkdir(exist_ok=True)
    dies, cat = base_dies()
    if not sm.legal(dies):
        raise SystemExit("起点不合法")
    report = {"note": (
        "python 绝对温度不是 HotSpot。naive 是每颗芯片单独成硅再相加。"
        "stale_cache 只重算被动过的那颗，其余沿用旧材料图。"
        "ignore_k 用新功耗、旧材料图。频率实验的漏电固定在 340 K，延迟是公式。"
        "生产退火和默认填充没有改。"
    )}
    print("P1/P7 spreader, tim, fill", flush=True)
    report["conduction"] = exp_spread_and_fill(dies)
    save(report)
    print("P2 moves", flush=True)
    report["moves"] = exp_moves(dies, cat)
    save(report)
    print("P3 ladder", flush=True)
    report["ladder"] = exp_ladder(dies, cat)
    for row in report["ladder"]:
        (OUT / "layouts" / f"{row['stress_tag']}.json").write_text(
            json.dumps(row["dies"], indent=2) + "\n", encoding="utf-8")
    save(report)
    print("P4 frequency", flush=True)
    report["frequency"] = exp_frequency(dies, cat)
    save(report)
    print("P5 power modes", flush=True)
    report["power_modes"] = exp_power_modes(dies, cat)
    save(report)
    report["geometry"] = exp_geometry(cat)
    # 问题 7 的力学填充：只比这一套布局，硅和 organic。organic 用 ATSim 已有材料，不改默认。
    report["fill_stress_tags"] = ["ladder_start", "fill_organic"]
    organic = json.loads((OUT / "layouts" / "ladder_start.json").read_text(encoding="utf-8"))
    (OUT / "layouts" / "fill_organic.json").write_text(json.dumps(organic, indent=2) + "\n", encoding="utf-8")
    save(report)
    print(f"写出 {OUT / 'report.json'}", flush=True)


def run_stress() -> None:
    from atsim_sa_proxy import eval_atsim_proxy, write_proxy_json

    report = json.loads((OUT / "report.json").read_text(encoding="utf-8"))
    tags = [r["stress_tag"] for r in report["ladder"]]
    tags += [r["tag"] for r in report["power_modes"] if r["mode"] != "native"]
    tags += ["fill_organic"]
    (OUT / "stress").mkdir(exist_ok=True)
    seen: dict[str, Path] = {}
    for tag in tags:
        dest = OUT / "stress" / f"{tag}.json"
        dies = json.loads((OUT / "layouts" / f"{tag}.json").read_text(encoding="utf-8"))
        fill = "organic" if tag == "fill_organic" else "silicon"
        sig = fill + json.dumps(dies, sort_keys=True)
        if dest.exists():
            seen.setdefault(sig, dest)
            print(f"skip {tag}", flush=True)
            continue
        if sig in seen:
            keep = json.loads(seen[sig].read_text(encoding="utf-8"))
            keep["tag"] = tag
            keep["copied_from"] = seen[sig].stem
            dest.write_text(json.dumps(keep, indent=2) + "\n", encoding="utf-8")
            print(f"copy {tag} from {seen[sig].stem}", flush=True)
            continue
        rects = sm.to_rects(dies)
        jp = OUT / "stress" / f"{tag}_case.json"
        overrides = {"fill_material": "organic"} if tag == "fill_organic" else None
        write_proxy_json(rects, jp, tag, package_overrides=overrides)
        st = eval_atsim_proxy(jp, verbose=False)
        cpu = next(c for c in st["chip_interiors"] if c["name"] in ("Penryn", "T2", "Xeon"))
        keep = {
            "tag": tag,
            "fill_material": "organic" if tag == "fill_organic" else "silicon",
            "vm_max_MPa": st.get("vm_normal_chip_max_MPa"),
            "vm_name": st.get("vm_normal_chip_name"),
            "cpu_mean_MPa": cpu.get("vm_mean_MPa"),
            "cpu_max_MPa": cpu.get("vm_max_MPa"),
            "cpu_n": cpu.get("n_cells"),
            "atsim_Tmax_C": st.get("Tmax_C"),
            "elapsed_s": st.get("elapsed_s"),
            "ndof": st.get("ndof"),
        }
        dest.write_text(json.dumps(keep, indent=2) + "\n", encoding="utf-8")
        seen[sig] = dest
        print(f"stress {tag}  max {keep['vm_max_MPa']:.1f}  cpuMean {keep['cpu_mean_MPa']:.1f}  "
              f"{keep['elapsed_s']:.0f}s", flush=True)


def _font():
    import matplotlib.pyplot as plt
    from matplotlib import font_manager
    font_manager.fontManager.addfont("/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc")
    plt.rcParams["font.family"] = "Noto Sans CJK JP"
    plt.rcParams["axes.unicode_minus"] = False
    return plt


def _stress(tag: str) -> dict | None:
    p = OUT / "stress" / f"{tag}.json"
    return json.loads(p.read_text(encoding="utf-8")) if p.exists() else None


def plot() -> None:
    plt = _font()
    report = json.loads((OUT / "report.json").read_text(encoding="utf-8"))
    fig_dir = OUT / "figures"
    fig_dir.mkdir(exist_ok=True)
    cond = report["conduction"]

    # 1. 误差和横向热占比随均热片导热
    sp = cond["spreader_sweep"]
    fig, ax = plt.subplots(figsize=(7.2, 4.4), constrained_layout=True)
    xs = [r["k_spreader"] for r in sp]
    ax.plot(xs, [r["naive"]["max_abs_K"] for r in sp], "o-", color="#c4554d", label="朴素叠加最大误差")
    ax.plot(xs, [abs(r["frozen"]["max_abs_K"]) for r in sp], "s--", color="#2a78d6", label="材料图冻结后再叠加")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("均热片导热系数 W/(m·K)")
    ax.set_ylabel("最大误差 / K")
    ax2 = ax.twinx()
    ax2.plot(xs, [r["lat_frac"] * 100 for r in sp], "^-", color="#2f7d4a", label="CPU 横向热占比")
    ax2.set_ylabel("横向流出 / 芯片功耗 / %")
    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, fontsize=8)
    ax.set_title("填充固定为 0.8。铜均热片还在时，叠加误差贴着地。")
    fig.savefig(fig_dir / "p1_spreader.png", dpi=140)
    plt.close(fig)

    # TIM 单独拧，均热片仍是 400
    tim = cond["tim_sweep"]
    fig, ax = plt.subplots(figsize=(6.4, 4.0), constrained_layout=True)
    ax.plot([r["k_tim"] for r in tim], [r["naive_max_abs_K"] for r in tim], "o-", color="#c4554d")
    ax.set_xlabel("TIM 导热系数 W/(m·K)，均热片仍为 400")
    ax.set_ylabel("朴素叠加最大误差 / K")
    ax.set_title("只降低传热界面，不降低均热片")
    fig.savefig(fig_dir / "p7_tim.png", dpi=140)
    plt.close(fig)

    # 填充 × 均热片
    grid = cond["fill_by_spreader"]
    fig, ax = plt.subplots(figsize=(6.6, 4.2), constrained_layout=True)
    for k_sp, color in ((400.0, "#2a78d6"), (40.0, "#c4554d"), (4.0, "#2f7d4a")):
        sub = [r for r in grid if r["k_spreader"] == k_sp]
        ax.plot([r["k_fill"] for r in sub], [r["naive_max_abs_K"] for r in sub], "o-",
                color=color, label=f"均热片 {k_sp:g}")
    ax.set_xscale("log")
    ax.set_xlabel("填充导热系数 W/(m·K)")
    ax.set_ylabel("朴素叠加最大误差 / K")
    ax.legend(fontsize=8)
    ax.set_title("填充要不要紧，取决于均热片还通不通")
    fig.savefig(fig_dir / "p7_fill_x_spreader.png", dpi=140)
    plt.close(fig)

    # 2. 动作
    moves = [r for r in report["moves"] if r.get("legal")]
    names = []
    for r in moves:
        if r["move"] not in names and r["k_spreader"] == 400.0:
            names.append(r["move"])
    fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.4), constrained_layout=True)
    x = np.arange(len(names))
    for ax, key, title in (
        (axes[0], "err_stale_cache_K", "只重算被动的那颗"),
        (axes[1], "err_naive_K", "每颗单独成硅再相加"),
    ):
        for k_sp, color, off in ((400.0, "#2a78d6", -0.15), (4.0, "#c4554d", 0.15)):
            ys = []
            for n in names:
                row = next(r for r in moves if r["move"] == n and r["k_spreader"] == k_sp)
                ys.append(row[key])
            ax.bar(x + off, ys, width=0.3, color=color, label=f"均热片 {k_sp:g}")
        ax.set_xticks(x, names, rotation=30, ha="right", fontsize=8)
        ax.set_ylabel("最大误差 / K")
        ax.set_title(title)
        ax.legend(fontsize=8)
    fig.suptitle("起点都是 Penryn+3070。误差相对这一步的真解。旧材料图上的尖峰另记在表里。", fontsize=11)
    fig.savefig(fig_dir / "p2_moves.png", dpi=140)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.2, 4.4), constrained_layout=True)
    for k_sp, color in ((400.0, "#2a78d6"), (4.0, "#c4554d")):
        sub = [r for r in moves if r["k_spreader"] == k_sp]
        ax.plot([r["cells_changed"] for r in sub], [r["err_naive_K"] for r in sub],
                "o-", color=color, label=f"朴素叠加，均热片 {k_sp:g}")
        ax.plot([r["cells_changed"] for r in sub], [r["err_stale_cache_K"] for r in sub],
                "s--", color=color, label=f"只重算移动芯片，均热片 {k_sp:g}")
        for r in sub:
            ax.annotate(r["move"], (r["cells_changed"], r["err_naive_K"]), fontsize=7, color=color)
    ax.set_xlabel("材料发生变化的格子数")
    ax.set_ylabel("最大误差 / K")
    ax.legend(fontsize=7)
    ax.set_title("格子数本身不决定误差，均热片导热才决定")
    fig.savefig(fig_dir / "p2_cells.png", dpi=140)
    plt.close(fig)

    # 3. 两个求解器的升降
    lad = report["ladder"]
    fig, ax = plt.subplots(figsize=(7.4, 4.4), constrained_layout=True)
    idx = np.arange(len(lad))
    py0 = lad[0]["python_Tmax_C"]
    hs0 = lad[0]["hotspot_Tmax_C"]
    ax.plot(idx, [r["python_Tmax_C"] - py0 for r in lad], "o-", color="#c4554d", label="python 相对起点")
    ax.plot(idx, [r["hotspot_Tmax_C"] - hs0 for r in lad], "s-", color="#2a78d6", label="HotSpot 相对起点")
    atsim = []
    for r in lad:
        st = _stress(r["stress_tag"])
        atsim.append(None if st is None else st["atsim_Tmax_C"])
    if all(v is not None for v in atsim):
        a0 = atsim[0]
        ax.plot(idx, [v - a0 for v in atsim], "^-", color="#2f7d4a", label="ATSim 相对起点")
    ax.axhline(0, color="0.5", lw=0.6)
    ax.set_xticks(idx, [r["step"] for r in lad])
    ax.set_ylabel("ΔT / °C")
    ax.legend(fontsize=8)
    ax.set_title(f"起点 python {py0:.1f} °C，HotSpot {hs0:.1f} °C，不是同一个零点")
    fig.savefig(fig_dir / "p3_ranking.png", dpi=140)
    plt.close(fig)
    # 相邻两步的升降是否同号。绝对温度不放在一起比。
    rank_rows = []
    flips = {"python_vs_hotspot": 0, "python_vs_atsim": 0, "pairs": 0}
    for a, b in zip(lad, lad[1:]):
        dpy = b["python_Tmax_C"] - a["python_Tmax_C"]
        dhs = b["hotspot_Tmax_C"] - a["hotspot_Tmax_C"]
        sa, sb = _stress(a["stress_tag"]), _stress(b["stress_tag"])
        dat = None if sa is None or sb is None else sb["atsim_Tmax_C"] - sa["atsim_Tmax_C"]
        row = {"from": a["step"], "to": b["step"], "d_python_K": dpy, "d_hotspot_K": dhs, "d_atsim_K": dat}
        rank_rows.append(row)
        flips["pairs"] += 1
        if dpy * dhs < 0:
            flips["python_vs_hotspot"] += 1
        if dat is not None and dpy * dat < 0:
            flips["python_vs_atsim"] += 1
    (OUT / "ranking.json").write_text(json.dumps({"steps": rank_rows, "sign_flips": flips}, indent=2) + "\n", encoding="utf-8")

    # 4. 频率
    freq = [r for r in report["frequency"] if not r.get("template")]
    fig, axes = plt.subplots(1, 2, figsize=(10.2, 4.2), constrained_layout=True)
    colors = {"Penryn": "#c4554d", "T2": "#2a78d6", "Xeon": "#2f7d4a"}
    for cpu, color in colors.items():
        sub = [r for r in freq if r["cpu"] == cpu]
        axes[0].plot([r["s"] for r in sub], [r["hotspot_Tmax_C"] for r in sub], "o-", color=color, label=cpu)
        axes[1].plot([r["s"] for r in sub], [r["Theta_serial"] for r in sub], "o-", color=color, label=f"{cpu} 串行")
        axes[1].plot([r["s"] for r in sub], [r["Theta_parallel"] for r in sub], "s--", color=color, label=f"{cpu} 并行")
    axes[0].set_xlabel("动态功耗比例 s（漏电固定 340 K）")
    axes[0].set_ylabel("HotSpot 峰值 / °C")
    axes[0].legend(fontsize=8)
    axes[0].set_title("几何不动，只改功耗模型")
    axes[1].set_xlabel("s")
    axes[1].set_ylabel("相对吞吐（Penryn s=1 单线程 = 1）")
    axes[1].legend(fontsize=6, ncol=2)
    axes[1].set_title("两个延迟定义，赢家不同")
    fig.savefig(fig_dir / "p4_frequency.png", dpi=140)
    plt.close(fig)

    # 5. 应力因素。native 用梯子上的三颗，其余用 mode_ 标签。
    native_tag = {"Penryn": "ladder_start", "T2": "ladder_T2", "Xeon": "ladder_Xeon"}
    pts = []
    for r in report["power_modes"]:
        tag = native_tag[r["cpu"]] if r["mode"] == "native" else r["tag"]
        st = _stress(tag)
        if st is None:
            continue
        pts.append({**r, **st})
    fig, axes = plt.subplots(1, 2, figsize=(10.4, 4.4), constrained_layout=True)
    markers = {"native": "o", "same_power": "s", "same_density": "^"}
    for r in pts:
        c = colors[r["cpu"]]
        axes[0].scatter(r["density_W_mm2"], r["cpu_mean_MPa"], c=c, marker=markers[r["mode"]], s=50)
        axes[1].scatter(r["grad_K_per_mm"], r["cpu_max_MPa"], c=c, marker=markers[r["mode"]], s=50)
    axes[0].set_xlabel("功率密度 W/mm^2")
    axes[0].set_ylabel("CPU 内部平均应力 / MPa")
    axes[0].set_title("平均应力跟着密度，不跟着芯片名字")
    axes[1].set_xlabel("整层最大梯度 / K/mm")
    axes[1].set_ylabel("CPU 内部最大应力 / MPa")
    axes[1].set_title("最大应力不随梯度单调下降")
    # 图例
    for cpu, c in colors.items():
        axes[1].scatter([], [], c=c, marker="o", label=cpu)
    for mode, m in markers.items():
        axes[1].scatter([], [], c="0.3", marker=m, label=mode)
    axes[1].legend(fontsize=7)
    fig.savefig(fig_dir / "p5_stress.png", dpi=140)
    plt.close(fig)

    # 6. 权重。三项可靠性用 native 的 HotSpot、梯度、应力，延迟用 s=1 的公式。
    chips = []
    for cpu in ("Penryn", "T2", "Xeon"):
        lad = next(r for r in report["ladder"] if r["step"] == ("start" if cpu == "Penryn" else cpu))
        fr = next(r for r in freq if r["cpu"] == cpu and r["s"] == 1.0)
        st = _stress(native_tag[cpu])
        if st is None:
            continue
        chips.append({
            "cpu": cpu,
            "T": lad["hotspot_Tmax_C"],
            "grad": lad["grad_K_per_mm"],
            "sig": st["vm_max_MPa"],
            "sig_mean": st["cpu_mean_MPa"],
            "delay_s": 1.0 / fr["Theta_serial"],
            "delay_p": 1.0 / fr["Theta_parallel"],
        })
    if len(chips) == 3:
        def norm(key):
            xs = [c[key] for c in chips]
            lo, hi = min(xs), max(xs)
            return [(x - lo) / (hi - lo) if hi > lo else 0.0 for x in xs]

        keys_raw = ("T", "grad", "sig", "delay_s")
        # 归一之后扫延迟权重。其余权重均分。
        fig, axes = plt.subplots(1, 2, figsize=(10.2, 4.2), constrained_layout=True)
        ws = np.linspace(0, 1, 21)
        for ax, dkey, title in (
            (axes[0], "delay_s", "延迟 = 1/串行吞吐"),
            (axes[1], "delay_p", "延迟 = 1/并行吞吐"),
        ):
            cols = {c["cpu"]: [] for c in chips}
            Tn, Gn, Sn, Dn = norm("T"), norm("grad"), norm("sig"), norm(dkey)
            for w in ws:
                rest = (1.0 - w) / 3.0
                for i, c in enumerate(chips):
                    cols[c["cpu"]].append(rest * Tn[i] + rest * Gn[i] + rest * Sn[i] + w * Dn[i])
            for cpu, color in colors.items():
                ax.plot(ws, cols[cpu], color=color, label=cpu)
            ax.set_xlabel("延迟权重（其余三项平分）")
            ax.set_ylabel("归一化代价，越低越好")
            ax.legend(fontsize=8)
            ax.set_title(title)
        fig.suptitle("应力用的是内部最大值。换成平均应力，赢家会不同，见表。", fontsize=11)
        fig.savefig(fig_dir / "p6_weights.png", dpi=140)
        plt.close(fig)
        # 不归一的原始和
        raw_rows = []
        for c in chips:
            raw_rows.append({
                "cpu": c["cpu"],
                "raw_T_plus_sig_plus_grad": c["T"] + c["sig"] + c["grad"],
                "T": c["T"], "sig": c["sig"], "grad": c["grad"],
                "delay_serial": c["delay_s"], "delay_parallel": c["delay_p"],
            })
        winners = []
        for dkey, skey, label in (
            ("delay_s", "sig", "串行+最大应力"),
            ("delay_p", "sig", "并行+最大应力"),
            ("delay_s", "sig_mean", "串行+平均应力"),
            ("delay_p", "sig_mean", "并行+平均应力"),
        ):
            Tn, Gn, Sn, Dn = norm("T"), norm("grad"), norm(skey), norm(dkey)
            for w in (0.0, 0.25, 0.5, 0.75, 1.0):
                rest = (1.0 - w) / 3.0
                costs = [rest * Tn[i] + rest * Gn[i] + rest * Sn[i] + w * Dn[i] for i in range(3)]
                winners.append({
                    "rule": label, "w_lat": w,
                    "winner": chips[int(np.argmin(costs))]["cpu"],
                    "costs": {chips[i]["cpu"]: costs[i] for i in range(3)},
                })
        (OUT / "weights.json").write_text(
            json.dumps({"chips": chips, "raw": raw_rows, "winners": winners}, indent=2) + "\n", encoding="utf-8")

    # 7. 接触边，和填充材料的应力
    geo = report["geometry"]
    fig, ax = plt.subplots(figsize=(6.4, 4.0), constrained_layout=True)
    xs = np.arange(len(geo))
    ax.plot(xs, [r["facing_cim_mm"] for r in geo], "o-", color="#2a78d6", label="和 CIM 的接触边")
    ax.plot(xs, [r["facing_gpu_mm"] for r in geo], "s-", color="#2f7d4a", label="和 GPU 的接触边")
    ax.plot(xs, [r["gap_to_cim_mm"] for r in geo], "^-", color="#c4554d", label="到 CIM 的缝")
    ax.set_xticks(xs, [f"{r['cpu']}\n{r['grid_area_mm2']:.0f} mm²" for r in geo])
    ax.set_ylabel("mm")
    ax.legend(fontsize=8)
    ax.set_title("CIM 比每颗 CPU 都窄，接触边不随 CPU 变长")
    fig.savefig(fig_dir / "p7_contact.png", dpi=140)
    plt.close(fig)

    sil, org = _stress("ladder_start"), _stress("fill_organic")
    if sil and org:
        fig, ax = plt.subplots(figsize=(5.6, 3.8), constrained_layout=True)
        labels = ["最大应力", "CPU 平均", "ATSim 峰值 °C"]
        a = [sil["vm_max_MPa"], sil["cpu_mean_MPa"], sil["atsim_Tmax_C"]]
        b = [org["vm_max_MPa"], org["cpu_mean_MPa"], org["atsim_Tmax_C"]]
        x = np.arange(3)
        ax.bar(x - 0.15, a, width=0.3, color="#2a78d6", label="填充硅")
        ax.bar(x + 0.15, b, width=0.3, color="#c4554d", label="填充 organic")
        ax.set_xticks(x, labels)
        ax.legend(fontsize=8)
        ax.set_title("同一布局。ATSim 温度不是 HotSpot。")
        fig.savefig(fig_dir / "p7_fill_stress.png", dpi=140)
        plt.close(fig)
    print(f"图在 {fig_dir}", flush=True)


def move_row(dies, trial, name: str, k_sp: float) -> dict:
    old_parts = frozen_parts(dies, 0.8, 4.0, k_sp)
    st_old = truth(dies, 0.8, 4.0, k_sp)
    st_new = truth(trial, 0.8, 4.0, k_sp)
    moved = [i for i, (a, b) in enumerate(zip(dies, trial))
             if a["name"] != b["name"] or a["j0"] != b["j0"] or a["ib0"] != b["ib0"]
             or a["w"] != b["w"] or a["h"] != b["h"]]
    k_new = paint(trial, 0.8)
    parts = list(old_parts)
    for i in moved:
        parts[i] = solve(k_new, hw.power_map(zero_but(trial, i)), 4.0, k_sp)[0]
    ignore = solve(paint(dies, 0.8), hw.power_map(trial), 4.0, k_sp)[0]
    nav = naive_sum(trial, 0.8, 4.0, k_sp)
    d = ignore - st_new[0]
    same = np.abs(paint(dies, 0.8) - k_new) < 1e-6
    rec = {
        "move": name, "k_spreader": k_sp, "legal": True,
        "cells_changed": len(cells_of(dies) ^ cells_of(trial)),
        "hardware": [d["name"] for d in trial],
        "field_moved_K": err(st_new[0], st_old[0])["max_abs_K"],
        "err_stale_cache_K": err(sum_parts(parts), st_new[0])["max_abs_K"],
        "err_ignore_k_K": float(np.max(np.abs(d))),
        "err_ignore_k_rmse_K": float(np.sqrt(np.mean(d ** 2))),
        "err_ignore_k_same_mat_K": float(np.max(np.abs(d[same]))) if same.any() else None,
        "err_naive_K": err(nav, st_new[0])["max_abs_K"],
        "python_Tmax_C": tmax(st_new[0]),
    }
    print(f"  {name:12} k {k_sp:5.0f}  cells {rec['cells_changed']:4d}  "
          f"stale {rec['err_stale_cache_K']:.4f}  naive {rec['err_naive_K']:.4f}  "
          f"spike {rec['err_ignore_k_K']:.1f}  same {rec['err_ignore_k_same_mat_K']:.3f}", flush=True)
    return rec


def extra() -> None:
    """三格挪动，以及旧材料图误差的尖峰和未改区。不重跑已经完成的热学扫描。"""
    report = json.loads((OUT / "report.json").read_text(encoding="utf-8"))
    dies, _cat = base_dies()
    for r in report["moves"]:
        if r["move"] == "nudge_2" and not r.get("legal"):
            r["why"] = "x 投影刚好贴边，间隙函数返回 0；y 方向仍隔 14 格。不是热学结果。"
    trial = deepcopy(dies)
    trial[0]["j0"] += 3
    if not sm.legal(trial):
        raise SystemExit("nudge_3 不合法")
    for k_sp in (400.0, 4.0):
        row = move_row(dies, trial, "nudge_3", k_sp)
        report["moves"] = [r for r in report["moves"] if not (r["move"] == "nudge_3" and r["k_spreader"] == k_sp)]
        # 插在 nudge_1 后面，两套导热各自一组。
        inserted = False
        out = []
        for r in report["moves"]:
            out.append(r)
            if not inserted and r["move"] == "nudge_1" and r["k_spreader"] == k_sp:
                out.append(row)
                inserted = True
        if not inserted:
            out.append(row)
        report["moves"] = out
    # 已有动作补上尖峰的均方根和未改材料区。nudge_3 刚算过。
    for k_sp in (400.0, 4.0):
        for r in report["moves"]:
            if not r.get("legal") or r["k_spreader"] != k_sp or "err_ignore_k_rmse_K" in r:
                continue
            trial = None
            for name, lay in make_moves(dies, _cat):
                if name == r["move"]:
                    trial = lay
            if trial is None:
                continue
            st_new = truth(trial, 0.8, 4.0, k_sp)[0]
            ignore = solve(paint(dies, 0.8), hw.power_map(trial), 4.0, k_sp)[0]
            d = ignore - st_new
            same = np.abs(paint(dies, 0.8) - paint(trial, 0.8)) < 1e-6
            r["err_ignore_k_rmse_K"] = float(np.sqrt(np.mean(d ** 2)))
            r["err_ignore_k_same_mat_K"] = float(np.max(np.abs(d[same]))) if same.any() else None
            print(f"  detail {r['move']} k {k_sp:5.0f}  rmse {r['err_ignore_k_rmse_K']:.3f}  "
                  f"same {r['err_ignore_k_same_mat_K']:.3f}", flush=True)
    if not any(r["step"] == "nudge_3" for r in report["ladder"]):
        lay = deepcopy(dies)
        lay[0]["j0"] += 3
        st = truth(lay, 130.0, 4.0, 400.0)
        hs_c = sm.hotspot(lay, OUT / "hotspot" / "ladder_nudge_3")
        rec = {
            "step": "nudge_3",
            "python_Tmax_C": tmax(st[0]),
            "grad_K_per_mm": sm.grad_max(st[0]),
            "hotspot_Tmax_C": hs_c,
            "dies": pack(lay),
            "stress_tag": "ladder_nudge_3",
        }
        idx = next(i for i, r in enumerate(report["ladder"]) if r["step"] == "nudge_1")
        report["ladder"].insert(idx + 1, rec)
        (OUT / "layouts" / "ladder_nudge_3.json").write_text(
            json.dumps(rec["dies"], indent=2) + "\n", encoding="utf-8")
        print(f"  ladder nudge_3 python {rec['python_Tmax_C']:.2f} HotSpot {hs_c:.2f}", flush=True)
    save(report)


if __name__ == "__main__":
    phase = sys.argv[1] if len(sys.argv) > 1 else "thermal"
    if phase == "thermal":
        thermal()
    elif phase == "extra":
        extra()
    elif phase == "stress":
        run_stress()
    elif phase == "plot":
        plot()
    else:
        raise SystemExit(phase)
