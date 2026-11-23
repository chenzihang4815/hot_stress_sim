#!/usr/bin/env python3
"""退火的选择里加上换硬件，代价同时看峰值温度和温度梯度。

不改 sa_thermal_feedback.py。动作三种，每次只抽一种：
  nudge  挪一格或两格
  rotate 转 90°（原来就有的“改方向”）
  hardware  把最尖的那颗 CPU/GPU 换成仓库里放得下、而且当前代价更低的更大硬件

代价三种，主实验用 both：
  tmax  Python 场的峰值（°C）。不是 HotSpot。
  grad  芯片层最大梯度（K/mm）
  both  (T-45)/(T0-45) + grad/grad0，两项在起点都是 1

HotSpot 峰值和 ATSim 芯片内部 von Mises 只在起点和最优点各算一次。
ATSim 必须用 fenicsx 解释器跑 phase=stress。芯片是硅，填充沿用默认，不改生产配置。
"""

from __future__ import annotations

import json
import random
import sys
from copy import deepcopy
from pathlib import Path

import numpy as np

ROOT = Path("/home/czh/任务仿真产生功耗文件工具")
SRC = ROOT / "热仿真流水线" / "src"
sys.path.insert(0, str(SRC))

import run_fill_superposition as hs  # noqa: E402
import run_sa_hw_smooth as sm  # noqa: E402

OUT = ROOT / "热仿真流水线" / "runs" / "sa_hw_both"
HS_BIN = sm.HS_BIN
N = sm.N
CW_MM = sm.CW_MM
N_STEPS = 30
ALPHA = 0.90
T_AMB_C = 45.0

# 主案例：不同起点。代价都是 both，并各跑一条不许换硬件的对照。
CASES = [
    {"id": "penryn_3070_s7", "cpu": "Penryn", "gpu": "GPU_3070", "seed": 7},
    {"id": "penryn_2060s_s7", "cpu": "Penryn", "gpu": "GPU_2060S", "seed": 7},
    {"id": "penryn_titanx_s7", "cpu": "Penryn", "gpu": "GPU_TITANX", "seed": 7},
    {"id": "penryn_3070_s11", "cpu": "Penryn", "gpu": "GPU_3070", "seed": 11},
    {"id": "penryn_3070_s21", "cpu": "Penryn", "gpu": "GPU_3070", "seed": 21},
    {"id": "t2_3070_s7", "cpu": "T2", "gpu": "GPU_3070", "seed": 7},
    {"id": "t2_titanx_s7", "cpu": "T2", "gpu": "GPU_TITANX", "seed": 7},
    {"id": "xeon_3070_s7", "cpu": "Xeon", "gpu": "GPU_3070", "seed": 7},
]
# 同一个起点、同一串提案，只换代价。用来看“加上梯度”会不会改选择。
COST_ARM = {"id": "penryn_3070_s7", "costs": ("tmax", "grad", "both")}


def score(ev: dict, scales: dict) -> float:
    kind = scales["kind"]
    if kind == "tmax":
        return float(ev["Tmax_C"])
    if kind == "grad":
        return float(ev["grad_max"])
    return (ev["Tmax_C"] - T_AMB_C) / scales["dT"] + ev["grad_max"] / scales["g0"]


def schedule(kind: str, ev0: dict) -> float:
    if kind == "tmax":
        return 1.5
    if kind == "grad":
        return max(0.15, 0.25 * ev0["grad_max"])
    return 0.20


def pick_kind(rng: random.Random, allow_hw: bool) -> str:
    r = rng.random()
    if allow_hw:
        if r < 0.45:
            return "nudge"
        if r < 0.65:
            return "rotate"
        return "hardware"
    if r < 0.70:
        return "nudge"
    return "rotate"


def rotate_move(dies: list[dict], rng: random.Random) -> tuple[list[dict], str]:
    trial = deepcopy(dies)
    i = rng.choice([k for k, d in enumerate(trial) if d["role"] in ("CPU", "GPU")])
    trial[i]["w"], trial[i]["h"] = trial[i]["h"], trial[i]["w"]
    return trial, f"rotate {trial[i]['name']}"


def hardware_choice(dies, info, catalog, scales) -> tuple[list[dict] | None, str, dict | None]:
    """最尖的那颗，换成放得下、而且当前代价最低的更大硬件。左下角不动。"""
    ranked = sorted(info["dies"], key=lambda r: -r["sharp_K_per_mm"])
    for victim in ranked:
        idx = next(i for i, d in enumerate(dies) if d["name"] == victim["name"])
        cur = dies[idx]
        cands = [p for p in catalog[cur["role"]] if p["area_mm2"] > cur["area_mm2"] + 1e-6]
        best = None
        for part in cands:
            trial = deepcopy(dies)
            trial[idx] = sm.clonable(part)
            trial[idx]["j0"] = cur["j0"]
            trial[idx]["ib0"] = cur["ib0"]
            if not sm.legal(trial):
                continue
            ev = sm.evaluate(trial, info["fac"], info["g"])
            j = score(ev, scales)
            if best is None or j < best[0]:
                best = (j, trial, part["name"], ev)
        if best is not None:
            _, trial, name, ev = best
            op = f"swap {cur['name']} -> {name} (sharp {victim['sharp_K_per_mm']:.3f})"
            return trial, op, ev
    return None, "no larger part fits", None


def accept(j_cur: float, j_new: float, temp: float, rng: random.Random) -> tuple[bool, float]:
    """接受随机数每次都抽，代价不同的臂才共用同一条随机序列。"""
    delta = j_cur - j_new
    u = rng.random()
    if delta >= 0:
        return True, 1.0
    ap = float(np.exp(delta / max(temp, 1e-9)))
    return u < ap, ap


def start_from(cat: dict, cpu_name: str, gpu_name: str) -> list[dict]:
    cpu = next(c for c in cat["CPU"] if c["name"] == cpu_name)
    gpu = next(g for g in cat["GPU"] if g["name"] == gpu_name)
    dies = [sm.clonable(cpu), sm.clonable(gpu)] + [sm.clonable(f) for f in cat["fixed"]]
    if not sm.legal(dies):
        raise RuntimeError(f"起点不合法: {cpu_name}+{gpu_name}")
    return dies


def run_sa(start, catalog, fac, g_amb, seed: int, allow_hw: bool, kind: str) -> dict:
    rng_m = random.Random(seed)
    rng_a = random.Random(seed + 7919)
    current = deepcopy(start)
    cur = sm.evaluate(current, fac, g_amb)
    scales = {"kind": kind, "dT": max(cur["Tmax_C"] - T_AMB_C, 1e-6), "g0": max(cur["grad_max"], 1e-9)}
    j_cur = score(cur, scales)
    best = deepcopy(current)
    best_ev = cur
    best_j = j_cur
    temp = schedule(kind, cur)
    log = [{
        "step": 0, "op": "init", "move": "init", "accept": True,
        "J": j_cur, "grad_max": cur["grad_max"], "grad_rms": cur["grad_rms"],
        "Tmax_C": cur["Tmax_C"], "hardware": [d["name"] for d in current],
    }]
    label = f"{'hw' if allow_hw else 'layout'}/{kind}"
    print(
        f"[{label}] seed {seed}  J {j_cur:.3f}  grad {cur['grad_max']:.3f}  "
        f"T {cur['Tmax_C']:.2f}  {[d['name'] for d in current]}",
        flush=True,
    )
    for step in range(1, N_STEPS + 1):
        move = pick_kind(rng_m, allow_hw)
        ev = None
        if move == "nudge":
            trial, op = sm.layout_move(current, rng_m)
            ok = sm.legal(trial)
        elif move == "rotate":
            trial, op = rotate_move(current, rng_m)
            ok = sm.legal(trial)
        else:
            info = {"fac": fac, "g": g_amb, "dies": cur["dies"]}
            trial, op, ev = hardware_choice(current, info, catalog, scales)
            ok = trial is not None
        if not ok:
            log.append({"step": step, "op": op, "move": move, "accept": False, "reason": "illegal"})
            temp *= ALPHA
            continue
        if ev is None:
            ev = sm.evaluate(trial, fac, g_amb)
        j_new = score(ev, scales)
        take, ap = accept(j_cur, j_new, temp, rng_a)
        rec = {
            "step": step, "op": op, "move": move, "accept": take, "ap": ap,
            "J": j_new, "grad_max": ev["grad_max"], "grad_rms": ev["grad_rms"],
            "Tmax_C": ev["Tmax_C"], "hardware": [d["name"] for d in trial], "T": temp,
        }
        if take:
            current = trial
            cur = ev
            j_cur = j_new
            if j_new < best_j - 1e-12:
                best = deepcopy(trial)
                best_ev = ev
                best_j = j_new
                rec["new_best"] = True
                print(
                    f"  {step:02d} BEST {op}  J {j_new:.3f}  grad {ev['grad_max']:.3f}  "
                    f"T {ev['Tmax_C']:.2f}",
                    flush=True,
                )
            else:
                print(f"  {step:02d} take {op}  J {j_new:.3f}  grad {ev['grad_max']:.3f}", flush=True)
        log.append(rec)
        temp *= ALPHA
    return {
        "best": best,
        "best_ev": {k: v for k, v in best_ev.items() if k != "T"},
        "init_ev": {k: v for k, v in sm.evaluate(start, fac, g_amb).items() if k != "T"},
        "T_init": sm.evaluate(start, fac, g_amb)["T"],
        "T_best": best_ev["T"],
        "log": log,
        "scales": scales,
        "final_hardware": [d["name"] for d in current],
    }


def pack(dies: list[dict]) -> list[dict]:
    return [sm.dump_die(d) for d in dies]


def grad_mag(T: np.ndarray) -> np.ndarray:
    gx = np.zeros_like(T)
    gy = np.zeros_like(T)
    gx[:, :-1] = np.abs(np.diff(T, axis=1)) / CW_MM
    gy[:-1, :] = np.abs(np.diff(T, axis=0)) / CW_MM
    return np.hypot(gx, gy)


def thermal() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "fields").mkdir(exist_ok=True)
    (OUT / "layouts").mkdir(exist_ok=True)
    cat = sm.warehouse()
    k = np.full((N, N), hs.K_SI)
    mat, g_amb = hs.build_operator(k)
    from scipy.sparse.linalg import splu
    fac = splu(mat.tocsc())

    endpoints = []
    rows = []

    def add_endpoint(tag: str, dies: list[dict], T: np.ndarray) -> None:
        np.save(OUT / "fields" / f"{tag}.npy", T)
        (OUT / "layouts" / f"{tag}.json").write_text(
            json.dumps(pack(dies), indent=2) + "\n", encoding="utf-8")
        endpoints.append(tag)

    def hotspot_of(tag: str, dies: list[dict]) -> float:
        work = OUT / "hotspot" / tag
        val = sm.hotspot(dies, work)
        return float(val)

    # 同一锚点只换 CPU，不退火。用来看面积和梯度、应力的关系。
    static = []
    for cpu_name in ("Penryn", "T2", "Xeon"):
        dies = start_from(cat, cpu_name, "GPU_3070")
        ev = sm.evaluate(dies, fac, g_amb)
        tag = f"static_{cpu_name}"
        add_endpoint(tag, dies, ev["T"])
        hs_c = hotspot_of(tag, dies)
        rec = {
            "id": tag, "kind": "static", "cpu": cpu_name, "gpu": "GPU_3070",
            "hardware": [d["name"] for d in dies],
            "grad_max": ev["grad_max"], "grad_rms": ev["grad_rms"],
            "python_Tmax_C": ev["Tmax_C"], "hotspot_Tmax_C": hs_c,
            "dies": pack(dies),
        }
        static.append(rec)
        print(
            f"static {cpu_name}: grad {ev['grad_max']:.3f}  "
            f"python {ev['Tmax_C']:.2f}  HotSpot {hs_c:.2f}",
            flush=True,
        )

    seen_init: dict[tuple[str, str], str] = {}
    sa_rows = []
    for case in CASES:
        start = start_from(cat, case["cpu"], case["gpu"])
        key = (case["cpu"], case["gpu"])
        costs = ("both",)
        if case["id"] == COST_ARM["id"]:
            costs = COST_ARM["costs"]
        for allow_hw, mode in ((False, "layout"), (True, "hardware")):
            for kind in costs:
                if mode == "layout" and kind != "both":
                    continue
                res = run_sa(start, cat, fac, g_amb, case["seed"], allow_hw, kind)
                tag_best = f"{case['id']}__{mode}__{kind}__best"
                add_endpoint(tag_best, res["best"], res["T_best"])
                if key not in seen_init:
                    tag_init = f"{case['id']}__init"
                    add_endpoint(tag_init, start, res["T_init"])
                    seen_init[key] = tag_init
                    hs_init = hotspot_of(tag_init, start)
                else:
                    tag_init = seen_init[key]
                    hs_init = next(r["hotspot_init_C"] for r in sa_rows if r["init_tag"] == tag_init)
                hs_best = hotspot_of(tag_best, res["best"])
                row = {
                    "id": case["id"], "mode": mode, "cost": kind,
                    "seed": case["seed"], "cpu": case["cpu"], "gpu": case["gpu"],
                    "init_tag": tag_init, "best_tag": tag_best,
                    "hardware_init": [d["name"] for d in start],
                    "hardware_best": [d["name"] for d in res["best"]],
                    "hardware_final": res["final_hardware"],
                    "grad_init": res["log"][0]["grad_max"],
                    "grad_best": res["best_ev"]["grad_max"],
                    "grad_rms_best": res["best_ev"]["grad_rms"],
                    "python_T_init_C": res["log"][0]["Tmax_C"],
                    "python_T_best_C": res["best_ev"]["Tmax_C"],
                    "hotspot_init_C": hs_init,
                    "hotspot_best_C": hs_best,
                    "J_init": res["log"][0]["J"],
                    "J_best": min(rec["J"] for rec in res["log"] if rec.get("accept")),
                    "swaps": [rec["op"] for rec in res["log"] if rec.get("move") == "hardware" and rec.get("accept")],
                }
                sa_rows.append(row)
                (OUT / "traces" / f"{case['id']}__{mode}__{kind}.json").parent.mkdir(exist_ok=True)
                (OUT / "traces" / f"{case['id']}__{mode}__{kind}.json").write_text(
                    json.dumps(res["log"], indent=2) + "\n", encoding="utf-8")
                print(
                    f"  -> grad {row['grad_init']:.3f} -> {row['grad_best']:.3f}  "
                    f"HotSpot {hs_init:.2f} -> {hs_best:.2f}  best {row['hardware_best']}",
                    flush=True,
                )

    report = {
        "objective": (
            "both = (T_python-45)/(T_init-45) + grad/grad_init。"
            "两项在起点都是 1。tmax 和 grad 是对照，只在 Penryn+3070、种子 7 上跑。"
        ),
        "moves": "nudge 0.45 / rotate 0.20 / hardware 0.35。layout 对照没有 hardware，nudge 0.70 / rotate 0.30。",
        "note": (
            "Python 场是整层硅、没有 spreader 外伸，绝对温度不能当成 HotSpot。"
            "CPU 功耗是各自 McPAT 380 K 总功耗，不是把 matmul 重跑到这三颗核上。"
            "GPU 功耗是同一份 backprop 在三张卡上的平均值。只允许换成面积更大的硬件。"
            "最优是代价最低的那一步，不是退火结束时的状态。"
            "应力在 phase=stress，指标是正常芯片内部 von Mises 最大，不是封装最大。"
        ),
        "static": static,
        "rows": sa_rows,
        "endpoints": endpoints,
    }
    (OUT / "report.json").write_text(
        json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"写出 {OUT / 'report.json'}", flush=True)


def run_stress() -> None:
    """fenicsx 解释器。每个点算完就落盘，已有的跳过。"""
    from atsim_sa_proxy import eval_atsim_proxy, write_proxy_json

    report = json.loads((OUT / "report.json").read_text(encoding="utf-8"))
    tags = list(report["endpoints"])
    (OUT / "stress").mkdir(exist_ok=True)
    for tag in tags:
        dest = OUT / "stress" / f"{tag}.json"
        if dest.exists():
            print(f"skip {tag}", flush=True)
            continue
        dies = json.loads((OUT / "layouts" / f"{tag}.json").read_text(encoding="utf-8"))
        rects = sm.to_rects(dies)
        jp = OUT / "stress" / f"{tag}_case.json"
        write_proxy_json(rects, jp, tag)
        st = eval_atsim_proxy(jp, verbose=False)
        keep = {
            "tag": tag,
            "vm_normal_chip_max_MPa": st.get("vm_normal_chip_max_MPa"),
            "vm_normal_chip_name": st.get("vm_normal_chip_name"),
            "vm_normal_chip_p95_MPa": st.get("vm_normal_chip_p95_MPa"),
            "w_pv_um": st.get("w_pv_um"),
            "atsim_Tmax_C": st.get("Tmax_C"),
            "elapsed_s": st.get("elapsed_s"),
            "ndof": st.get("ndof"),
            "chips": [
                {k: r.get(k) for k in ("name", "normal", "n_cells", "vm_max_MPa", "vm_mean_MPa", "vm_p95_MPa")}
                for r in (st.get("chip_interiors") or [])
            ],
        }
        dest.write_text(json.dumps(keep, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
        print(
            f"stress {tag}  σ={keep['vm_normal_chip_max_MPa']:.2f} MPa"
            f" ({keep['vm_normal_chip_name']})  {keep['elapsed_s']:.1f} s",
            flush=True,
        )


def _font():
    import matplotlib.pyplot as plt
    from matplotlib import font_manager
    font_manager.fontManager.addfont("/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc")
    plt.rcParams["font.family"] = "Noto Sans CJK JP"
    plt.rcParams["axes.unicode_minus"] = False
    return plt


def _outline(ax, dies: list[dict]) -> None:
    from matplotlib.patches import Rectangle
    for d in dies:
        ax.add_patch(Rectangle(
            (d["j0"] * CW_MM, d["ib0"] * CW_MM), d["w"] * CW_MM, d["h"] * CW_MM,
            fill=False, edgecolor="white", lw=0.7,
        ))
        ax.text(
            (d["j0"] + d["w"] / 2) * CW_MM, (d["ib0"] + d["h"] / 2) * CW_MM,
            d["name"].replace("GPU_", ""),
            ha="center", va="center", fontsize=6, color="white",
        )


def _load_dies(tag: str) -> list[dict]:
    return json.loads((OUT / "layouts" / f"{tag}.json").read_text(encoding="utf-8"))


def _load_T(tag: str) -> np.ndarray:
    return np.load(OUT / "fields" / f"{tag}.npy")


def _stress(tag: str) -> dict | None:
    p = OUT / "stress" / f"{tag}.json"
    if not p.exists():
        return None
    return json.loads(p.read_text(encoding="utf-8"))


def plot() -> None:
    plt = _font()
    from matplotlib.patches import Rectangle  # noqa: F401
    report = json.loads((OUT / "report.json").read_text(encoding="utf-8"))
    fig_dir = OUT / "figures"
    fig_dir.mkdir(exist_ok=True)

    # 1. 同一位置换 CPU
    static = report["static"]
    temps = [_load_T(s["id"]) for s in static]
    grads = [grad_mag(T) for T in temps]
    tlo, thi = min(T.min() for T in temps) - 273.15, max(T.max() for T in temps) - 273.15
    glo, ghi = 0.0, max(g.max() for g in grads)
    fig, axes = plt.subplots(2, 3, figsize=(12.6, 8.2), constrained_layout=True)
    for col, s, T, g in zip(range(3), static, temps, grads):
        dies = s["dies"]
        st = _stress(s["id"])
        sig = "" if st is None else f"\n内部 σ {st['vm_normal_chip_max_MPa']:.1f} MPa"
        im0 = axes[0, col].imshow(
            T[::-1] - 273.15, origin="lower", extent=(0, 72, 0, 72),
            vmin=tlo, vmax=thi, cmap="inferno",
        )
        _outline(axes[0, col], dies)
        axes[0, col].set_title(
            f"{s['cpu']}+3070\n梯度 {s['grad_max']:.2f} K/mm  HotSpot {s['hotspot_Tmax_C']:.2f} °C{sig}",
            fontsize=10,
        )
        im1 = axes[1, col].imshow(
            g[::-1], origin="lower", extent=(0, 72, 0, 72),
            vmin=glo, vmax=ghi, cmap="magma",
        )
        _outline(axes[1, col], dies)
        axes[1, col].set_title(f"|梯度| 最大 {g.max():.2f} K/mm", fontsize=10)
        for ax in (axes[0, col], axes[1, col]):
            ax.set_xlim(0, 72)
            ax.set_ylim(0, 72)
            ax.set_xlabel("x / mm")
            ax.set_ylabel("y / mm")
    fig.colorbar(im0, ax=axes[0, :].tolist(), fraction=0.046, pad=0.02, label="°C（快速场，不是 HotSpot）")
    fig.colorbar(im1, ax=axes[1, :].tolist(), fraction=0.046, pad=0.02, label="K/mm")
    fig.suptitle("同一锚点、同一 GPU 任务，只把 CPU 换成更大的一颗", fontsize=13)
    fig.savefig(fig_dir / "static_cpu.png", dpi=140)
    plt.close(fig)

    # 2. 退火前后热图。挑四条 both/hardware。
    want = ["penryn_3070_s7", "penryn_titanx_s7", "t2_3070_s7", "xeon_3070_s7"]
    picked = []
    for cid in want:
        row = next(r for r in report["rows"] if r["id"] == cid and r["mode"] == "hardware" and r["cost"] == "both")
        picked.append(row)
    panels_T = []
    for r in picked:
        panels_T.append(_load_T(r["init_tag"]))
        panels_T.append(_load_T(r["best_tag"]))
    lo = min(T.min() for T in panels_T) - 273.15
    hi = max(T.max() for T in panels_T) - 273.15
    fig, axes = plt.subplots(4, 2, figsize=(8.4, 14.2), constrained_layout=True)
    im = None
    for i, r in enumerate(picked):
        for j, tag, title_bit in (
            (0, r["init_tag"], "起点"),
            (1, r["best_tag"], "最优"),
        ):
            ax = axes[i, j]
            T = _load_T(tag)
            dies = _load_dies(tag)
            im = ax.imshow(T[::-1] - 273.15, origin="lower", extent=(0, 72, 0, 72),
                           vmin=lo, vmax=hi, cmap="inferno")
            _outline(ax, dies)
            st = _stress(tag)
            sig = "" if st is None else f"  σ {st['vm_normal_chip_max_MPa']:.1f}"
            grad = r["grad_init"] if j == 0 else r["grad_best"]
            hs_c = r["hotspot_init_C"] if j == 0 else r["hotspot_best_C"]
            hw = "+".join(n.replace("GPU_", "") for n in (r["hardware_init"] if j == 0 else r["hardware_best"])[:2])
            ax.set_title(f"{title_bit} {hw}\n梯度 {grad:.2f} K/mm  HotSpot {hs_c:.2f} °C{sig}", fontsize=9)
            ax.set_xlim(0, 72)
            ax.set_ylim(0, 72)
    fig.colorbar(im, ax=axes.ravel().tolist(), fraction=0.02, pad=0.01, label="°C（快速场，不是 HotSpot）")
    fig.suptitle("换硬件这一步放进退火之后：起点和代价最优（峰值 + 梯度）", fontsize=13)
    fig.savefig(fig_dir / "sa_heatmaps.png", dpi=140)
    plt.close(fig)

    # 3. 梯度 / HotSpot / 应力 的变化
    both_rows = [r for r in report["rows"] if r["cost"] == "both"]
    labels = []
    for r in both_rows:
        if r["mode"] != "hardware":
            continue
        labels.append(r)
    fig, axes = plt.subplots(1, 3, figsize=(13.2, 4.6), constrained_layout=True)
    names = [r["id"].replace("_", "\n") for r in labels]
    x = np.arange(len(labels))
    w = 0.38

    def pair(mode: str, key_i: str, key_b: str):
        out_i, out_b = [], []
        for r0 in labels:
            r = next(z for z in both_rows if z["id"] == r0["id"] and z["mode"] == mode)
            out_i.append(r[key_i])
            out_b.append(r[key_b])
        return out_i, out_b

    for ax, key_i, key_b, title, ylab in (
        (axes[0], "grad_init", "grad_best", "最大梯度", "K/mm"),
        (axes[1], "hotspot_init_C", "hotspot_best_C", "HotSpot 峰值", "°C"),
    ):
        # init 各案例相同，只画 hardware 的前后，layout 的最优作为对照
        _, hb = pair("hardware", key_i, key_b)
        ib, _ = pair("hardware", key_i, key_b)
        _, lb = pair("layout", key_i, key_b)
        ax.bar(x - w / 2, hb, width=w, color="#2a78d6", label="换硬件后的最优")
        ax.bar(x + w / 2, lb, width=w, color="#c4554d", label="只挪位置/转向")
        ax.plot(x, ib, "k.", label="起点")
        ax.set_xticks(x, names, fontsize=7)
        ax.set_title(title)
        ax.set_ylabel(ylab)
        ax.legend(fontsize=7)
    sig_h, sig_l, sig_i = [], [], []
    for r0 in labels:
        rh = next(z for z in both_rows if z["id"] == r0["id"] and z["mode"] == "hardware")
        rl = next(z for z in both_rows if z["id"] == r0["id"] and z["mode"] == "layout")
        sh, sl, si = _stress(rh["best_tag"]), _stress(rl["best_tag"]), _stress(rh["init_tag"])
        sig_h.append(np.nan if sh is None else sh["vm_normal_chip_max_MPa"])
        sig_l.append(np.nan if sl is None else sl["vm_normal_chip_max_MPa"])
        sig_i.append(np.nan if si is None else si["vm_normal_chip_max_MPa"])
    axes[2].bar(x - w / 2, sig_h, width=w, color="#2a78d6", label="换硬件后的最优")
    axes[2].bar(x + w / 2, sig_l, width=w, color="#c4554d", label="只挪位置/转向")
    axes[2].plot(x, sig_i, "k.", label="起点")
    axes[2].set_xticks(x, names, fontsize=7)
    axes[2].set_title("芯片内部 von Mises")
    axes[2].set_ylabel("MPa")
    axes[2].legend(fontsize=7)
    fig.suptitle("代价 = 峰值温度 + 温度梯度。黑点是起点。", fontsize=12)
    fig.savefig(fig_dir / "bars.png", dpi=140)
    plt.close(fig)

    # 4. 三种代价，同一条提案序列
    fig, axes = plt.subplots(1, 2, figsize=(10.4, 4.2), constrained_layout=True)
    for kind, color in (("tmax", "#c4554d"), ("grad", "#2a78d6"), ("both", "#2f7d4a")):
        path = OUT / "traces" / f"penryn_3070_s7__hardware__{kind}.json"
        log = json.loads(path.read_text(encoding="utf-8"))
        # 已接受状态随步数。被拒的步保持上一个值。
        g, t = [], []
        cg, ct = log[0]["grad_max"], log[0]["Tmax_C"]
        for rec in log:
            if rec.get("accept"):
                cg, ct = rec["grad_max"], rec["Tmax_C"]
            g.append(cg)
            t.append(ct)
        axes[0].plot(g, color=color, label=kind)
        axes[1].plot(t, color=color, label=kind)
    axes[0].set_title("Penryn+3070，种子 7，已接受的梯度")
    axes[0].set_ylabel("K/mm")
    axes[1].set_title("同一条链上的 Python 峰值（不是 HotSpot）")
    axes[1].set_ylabel("°C")
    for ax in axes:
        ax.set_xlabel("步")
        ax.legend(fontsize=8)
    fig.savefig(fig_dir / "cost_traces.png", dpi=140)
    plt.close(fig)

    # 5. 表
    fig, ax = plt.subplots(figsize=(14.5, 6.8))
    ax.axis("off")
    header = ["案例", "动作", "代价", "最优硬件", "梯度初", "梯度优", "HotSpot初", "HotSpot优", "σ初", "σ优", "σ在"]
    cells = []
    for s in static:
        st = _stress(s["id"])
        cells.append([
            f"{s['cpu']}+3070", "不退火", "—", s["cpu"],
            f"{s['grad_max']:.2f}", f"{s['grad_max']:.2f}",
            f"{s['hotspot_Tmax_C']:.2f}", f"{s['hotspot_Tmax_C']:.2f}",
            "—" if st is None else f"{st['vm_normal_chip_max_MPa']:.1f}",
            "—" if st is None else f"{st['vm_normal_chip_max_MPa']:.1f}",
            "" if st is None else str(st["vm_normal_chip_name"]),
        ])
    for r in report["rows"]:
        if not (r["cost"] == "both" or r["id"] == "penryn_3070_s7"):
            continue
        si, sb = _stress(r["init_tag"]), _stress(r["best_tag"])
        hw = "+".join(n.replace("GPU_", "") for n in r["hardware_best"][:2])
        cells.append([
            r["id"], "换硬件" if r["mode"] == "hardware" else "只布局", r["cost"], hw,
            f"{r['grad_init']:.2f}", f"{r['grad_best']:.2f}",
            f"{r['hotspot_init_C']:.2f}", f"{r['hotspot_best_C']:.2f}",
            "—" if si is None else f"{si['vm_normal_chip_max_MPa']:.1f}",
            "—" if sb is None else f"{sb['vm_normal_chip_max_MPa']:.1f}",
            "" if sb is None else str(sb["vm_normal_chip_name"]),
        ])
    table = ax.table(cellText=cells, colLabels=header, loc="center", cellLoc="center")
    table.auto_set_font_size(False)
    table.set_fontsize(7)
    table.scale(1, 1.25)
    ax.set_title("温度梯度、HotSpot 峰值、芯片内部应力。σ 是正常芯片内部 von Mises，不是封装最大。", fontsize=11, pad=12)
    fig.tight_layout()
    fig.savefig(fig_dir / "table.png", dpi=140)
    plt.close(fig)
    print(f"图在 {fig_dir}", flush=True)


if __name__ == "__main__":
    phase = sys.argv[1] if len(sys.argv) > 1 else "thermal"
    if phase == "thermal":
        thermal()
    elif phase == "stress":
        run_stress()
    elif phase == "plot":
        plot()
    else:
        raise SystemExit(f"unknown phase {phase}")
