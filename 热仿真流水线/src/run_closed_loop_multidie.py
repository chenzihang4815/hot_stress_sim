#!/usr/bin/env python3
"""多芯粒闭环的一次案例。不改生产 SA，不改 HotSpot 和 ATSim 的默认封装。

内环：平移、旋转，或把一颗 CPU / GPU 换成库里另一型。填充这一次选整层硅。
代价是 Python 芯粒层的峰值温度和最大梯度。定价用冻结材料图的场叠加：
每颗芯片单独的温升相加。硅是均匀的，算子只分解一次，相加应当只剩舍入误差。

外环：布局冻住。ATSim 算热和应力，对已有门槛（芯粒层峰值 < 90 C，
正常芯片内部 von Mises < 200 MPa）。频率不是全封装一个倍率。每颗 CPU、
每颗 GPU 自己的倍率，步长 0.1，只乘该芯片的动态功耗。CIM 和 HBM 不调。
搜索用同一套温升场按功耗比例相加，再把 Python 温升按起点的 HotSpot
块平均峰值缩放到 90 C 门槛上。缩放出的温度只是搜索用的向导。
认门槛的是 ATSim。HotSpot 另记一块，并出官方 SVG。

HBM 仍用颗数公式。型号变了、颗数没变，瓦数不会变。这条写进报告，不另造公式。
吞吐分列，CPU、GPU、CIM 不加总。
"""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path


ROOT = Path("/home/czh/任务仿真产生功耗文件工具")
SRC = ROOT / "热仿真流水线" / "src"
ATSIM_ROOT = ROOT / "热应力python工程_最终"
FENICS = Path("/home/czh/miniconda3/envs/fenicsx/bin/python")
HS_BIN = ROOT / "HotSpot-master" / "hotspot"
SVG_PL = ROOT / "HotSpot-master" / "scripts" / "grid_thermal_map.pl"
OUT = ROOT / "热仿真流水线" / "runs" / "closed_loop_multidie"
FONT = "/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc"

# 这一次跑的门槛带宽。不是新的项目规定。旧的 200–230 MPa 重退火带不用。
T_GATE = 90.0
S_GATE = 200.0
NEAR_T = 5.0
NEAR_S = 20.0
T_AMB_C = 45.0
S_LO = 0.5
S_HI = 2.5
S_STEP = 0.1


def atsim_entry(case_path: str, out_path: str) -> None:
    """只在 fenicsx 解释器里走这条。系统 Python 不进这里。"""
    import numpy as np

    sys.path.insert(0, str(ATSIM_ROOT))
    sys.path.insert(0, str(SRC))
    from atsim.elasticity import solve_elasticity  # noqa: E402
    from atsim.fvm3d import solve_fvm3d  # noqa: E402
    from atsim.geometry import load_case  # noqa: E402
    import atsim_sa_proxy as proxy  # noqa: E402

    case = load_case(case_path)
    th = solve_fvm3d(case, verbose=False)
    el = solve_elasticity(th, T_ref=case.T_amb, verbose=False)
    grid = th["grid"]
    layer_id = grid["layer_id"]
    die_ids = [
        i for i, layer in enumerate(case.stack)
        if layer.name == "die_fill" or layer.name.startswith("die_")
    ]
    mask = np.isin(layer_id, die_ids) & grid["solid"]
    temperature = th["T"]
    die_vals = temperature[mask] - 273.15
    col = np.zeros(mask.shape[:2], dtype=bool)
    acc = np.zeros(mask.shape[:2], dtype=float)
    cnt = np.zeros(mask.shape[:2], dtype=float)
    for k in range(mask.shape[2]):
        m = mask[:, :, k]
        col |= m
        acc[m] += temperature[:, :, k][m]
        cnt[m] += 1.0
    mean_k = np.divide(acc, cnt, out=np.zeros_like(acc), where=cnt > 0)
    dx = np.asarray(grid["dx"], dtype=float)
    dy = np.asarray(grid["dy"], dtype=float)
    gmax = 0.0
    for i in range(mean_k.shape[0] - 1):
        for j in range(mean_k.shape[1]):
            if col[i, j] and col[i + 1, j]:
                dist = 0.5 * (dx[i] + dx[i + 1])
                gmax = max(gmax, abs(mean_k[i + 1, j] - mean_k[i, j]) / dist * 1e-3)
    for i in range(mean_k.shape[0]):
        for j in range(mean_k.shape[1] - 1):
            if col[i, j] and col[i, j + 1]:
                dist = 0.5 * (dy[j] + dy[j + 1])
                gmax = max(gmax, abs(mean_k[i, j + 1] - mean_k[i, j]) / dist * 1e-3)

    interiors = proxy._chip_interior_vm(el)
    gate = []
    for row in interiors:
        name = row["name"]
        area = row["area_mm2"]
        if name.startswith("HBM") or name.startswith("anchor"):
            keep = False
        elif name.startswith("CIM"):
            keep = area >= 1.0
        else:
            keep = True
        row["gate_chip"] = keep
        if keep and row.get("vm_max_MPa") is not None:
            gate.append(row)
    worst = max(gate, key=lambda r: r["vm_max_MPa"]) if gate else None

    ox = getattr(case, "chip_ox", 0.0)
    oy = getattr(case, "chip_oy", 0.0)
    xc, yc = grid["xc"], grid["yc"]
    per_chip = []
    for block in case.blocks:
        x0, y0 = ox + block.x0, oy + block.y0
        ix = np.where((xc >= x0) & (xc <= x0 + block.w))[0]
        jy = np.where((yc >= y0) & (yc <= y0 + block.h))[0]
        tmax = None
        tmean = None
        if ix.size and jy.size:
            sub = mean_k[np.ix_(ix, jy)]
            m = col[np.ix_(ix, jy)]
            if m.any():
                vals = sub[m] - 273.15
                tmax = float(vals.max())
                tmean = float(vals.mean())
        per_chip.append({"name": block.name, "T_column_max_C": tmax, "T_column_mean_C": tmean})

    out = {
        "die_cell_max_C": float(die_vals.max()),
        "die_cell_min_C": float(die_vals.min()),
        "grad_max_K_per_mm": float(gmax),
        "vm_normal_chip_max_MPa": None if worst is None else float(worst["vm_max_MPa"]),
        "vm_normal_chip_name": None if worst is None else worst["name"],
        "vm_package_max_not_the_gate_MPa": float(el["vm"][grid["solid"]].max() * 1e-6),
        "per_chip_T": per_chip,
        "chip_interiors": interiors,
        "heat_s": float(th["elapsed_s"]),
        "elast_s": float(el["elapsed_s"]),
        "grid": [int(grid["nx"]), int(grid["ny"]), int(grid["nz"])],
        "note": (
            "门槛用 die_cell_max_C 和 vm_normal_chip_max_MPa。"
            "vm_package_max 是整包最大，含边缘奇点，不是门槛。"
            "HBM 不进正常芯片。网格是代理的 2 mm / 8 mm，不是对齐检查用的 1.125 mm。"
            "梯度是芯粒层柱平均温度的相邻差分，单位 K/mm。"
        ),
    }
    Path(out_path).write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(
        f"ATSim die {out['die_cell_max_C']:.2f} C  "
        f"vm {out['vm_normal_chip_max_MPa']:.2f} MPa ({out['vm_normal_chip_name']})  "
        f"grad {out['grad_max_K_per_mm']:.3f}",
        flush=True,
    )


if __name__ == "__main__" and len(sys.argv) > 1 and sys.argv[1] == "--atsim":
    atsim_entry(sys.argv[2], sys.argv[3])
    raise SystemExit(0)


import random
from copy import deepcopy

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import font_manager
from scipy.sparse.linalg import splu

sys.path.insert(0, str(SRC))
import run_fill_superposition as hs  # noqa: E402
import run_hw_area_field as hw  # noqa: E402
import run_sa_hw_smooth as sm  # noqa: E402
from hotspot_io import run_hotspot, write_hotspot_inputs  # noqa: E402

N = hs.N
CW_MM = hs.CW * 1e3
N_INNER = 40
N_FREQ = 24
ALPHA = 0.90
SEED = 11


CPUS = (
    {"part": "Penryn", "role": "CPU", "area_mm2": 91.1348, "dyn_W": 47.9315,
     "leak_W": 10.7217 + 0.588125, "cores": 2, "f_MHz": 3700,
     "theta_serial": 1.0, "theta_parallel": 2.0},
    {"part": "T2", "role": "CPU", "area_mm2": 273.523, "dyn_W": 50.421,
     "leak_W": 15.0682 + 3.50762, "cores": 8, "f_MHz": 1400,
     "theta_serial": 1400 / 3700, "theta_parallel": 8 * 1400 / 3700},
    {"part": "Xeon", "role": "CPU", "area_mm2": 421.676, "dyn_W": 66.8485,
     "leak_W": 35.4132 + 1.68429, "cores": 2, "f_MHz": 3400,
     "theta_serial": 3400 / 3700, "theta_parallel": 2 * 3400 / 3700},
)
GPUS = (
    {"part": "GPU_2060S", "role": "GPU", "area_mm2": 445.0, "dyn_W": 18.410518130135276,
     "leak_W": 20.38510912938833 + 32.3252, "f_MHz": 1905.0,
     "inst_per_s": 29305521412.565357},
    {"part": "GPU_3070", "role": "GPU", "area_mm2": 392.0, "dyn_W": 16.054421410705352,
     "leak_W": 18.50819493837828 + 32.3252, "f_MHz": 1132.0,
     "inst_per_s": 19081322297.51239},
    {"part": "GPU_TITANX", "role": "GPU", "area_mm2": 471.0, "dyn_W": 10.531824961055307,
     "leak_W": 24.55010003977329 + 32.3252, "f_MHz": 1417.0,
     "inst_per_s": 17410487521.12956},
)
CIMS = (
    ("CIM_jia_m64", 273.99848050829274, 0.02051619509144024, 64),
    ("CIM_colonnade_m512", 118.456264372224, 0.2967109271348416, 512),
    ("CIM_colonnade_m256", 59.228132186112, 0.29677526103697566, 256),
    ("CIM_basic_m256", 11.9601949587906, 0.0022112712579585407, 256),
    ("CIM_colonnade_m64", 14.807033046528, 0.29716126444977986, 64),
)


def hbm_power(n_cpu: int, n_gpu: int) -> float:
    return (0.088 + n_cpu * (1.035 - 0.088) + n_gpu * (0.844 - 0.088)) * 16


def part_from(spec: dict) -> dict:
    total = spec["dyn_W"] + spec["leak_W"]
    die = hw.prepare({"name": spec["part"], "area_mm2": spec["area_mm2"], "power_W": total})
    die.update(spec)
    die["s"] = 1.0
    die["power_W"] = total
    return die


def fixed_die(name: str, role: str, area: float, power: float, n_macros: int | None) -> dict:
    die = hw.prepare({"name": name, "area_mm2": area, "power_W": power})
    die.update({
        "role": role, "part": name, "slot": name, "s": 1.0,
        "dyn_W": 0.0, "leak_W": power, "n_macros": n_macros, "power_W": power,
    })
    return die


def catalog() -> dict:
    cpus = [part_from(s) for s in CPUS]
    gpus = [part_from(s) for s in GPUS]
    return {
        "CPU": cpus,
        "GPU": gpus,
        "by": {(d["role"], d["part"]): d for d in cpus + gpus},
    }


def place(cat: dict) -> list[dict]:
    """和 sa_many_macro 同一张货架。槽名带后缀，换成别的型号时后缀不变。"""
    def put(role: str, part: str, slot: str, j0: int, ib0: int) -> dict:
        die = deepcopy(cat["by"][(role, part)])
        die["name"] = f"{part}_{slot}"
        die["slot"] = slot
        die["j0"] = j0
        die["ib0"] = ib0
        die["s"] = 1.0
        return die

    dies = [
        put("CPU", "Xeon", "Ca", 0, 0),
        put("GPU", "GPU_2060S", "Ga", 17, 0),
        put("GPU", "GPU_TITANX", "Gb", 34, 0),
        fixed_die("CIM_jia_m64", "CIM", CIMS[0][1], CIMS[0][2], CIMS[0][3]),
        put("CPU", "Xeon", "Cb", 0, 23),
        put("GPU", "GPU_3070", "Gc", 17, 23),
        fixed_die("HBM3", "HBM", 121.0, 0.0, None),
        fixed_die("CIM_colonnade_m512", "CIM", CIMS[1][1], CIMS[1][2], CIMS[1][3]),
        fixed_die("CIM_colonnade_m256", "CIM", CIMS[2][1], CIMS[2][2], CIMS[2][3]),
        fixed_die("CIM_basic_m256", "CIM", CIMS[3][1], CIMS[3][2], CIMS[3][3]),
        fixed_die("CIM_colonnade_m64", "CIM", CIMS[4][1], CIMS[4][2], CIMS[4][3]),
    ]
    coords = {
        "CIM_jia_m64": (52, 0),
        "HBM3": (32, 23),
        "CIM_colonnade_m512": (41, 23),
        "CIM_colonnade_m256": (0, 46),
        "CIM_basic_m256": (7, 46),
        "CIM_colonnade_m64": (11, 46),
    }
    for die in dies:
        if die["name"] in coords:
            die["j0"], die["ib0"] = coords[die["name"]]
    refresh_hbm(dies)
    return dies


def refresh_hbm(dies: list[dict]) -> float:
    n_cpu = sum(d["role"] == "CPU" for d in dies)
    n_gpu = sum(d["role"] == "GPU" for d in dies)
    power = hbm_power(n_cpu, n_gpu)
    for die in dies:
        if die["role"] == "HBM":
            die["power_W"] = power
            die["leak_W"] = power
            die["n_cpu"] = n_cpu
            die["n_gpu"] = n_gpu
    return power


def apply_power(dies: list[dict]) -> None:
    for die in dies:
        if die["role"] in ("CPU", "GPU"):
            die["power_W"] = die["leak_W"] + die["s"] * die["dyn_W"]


def cache_key(die: dict) -> tuple:
    return (
        die["name"], die["j0"], die["ib0"], die["w"], die["h"],
        round(float(die["power_W"]), 6),
    )


def rise_of(die: dict, fac, g_amb, cache: dict) -> np.ndarray:
    key = cache_key(die)
    if key not in cache:
        state = hw.solve(np.full((N, N), hs.K_SI), hw.power_map([die]), fac, g_amb)
        cache[key] = state[0] - hs.T_AMB
    return cache[key]


def superpose(dies: list[dict], fac, g_amb, cache: dict) -> np.ndarray:
    acc = np.zeros((N, N))
    for die in dies:
        acc += rise_of(die, fac, g_amb, cache)
    return acc + hs.T_AMB


def grad_max(temperature: np.ndarray) -> float:
    gx = np.abs(np.diff(temperature, axis=1)).max() / CW_MM
    gy = np.abs(np.diff(temperature, axis=0)).max() / CW_MM
    return float(max(gx, gy))


def metrics(temperature: np.ndarray) -> dict:
    return {
        "Tmax_C": float(temperature.max() - 273.15),
        "grad_max": grad_max(temperature),
        "T": temperature,
    }


def score(ev: dict, scales: dict) -> float:
    return (ev["Tmax_C"] - T_AMB_C) / scales["dT"] + ev["grad_max"] / scales["g0"]


def nudge(dies: list[dict], rng: random.Random) -> tuple[list[dict], str]:
    trial = deepcopy(dies)
    i = rng.randrange(len(trial))
    step = rng.choice([1, 1, 2])
    dj, di = rng.choice([(1, 0), (-1, 0), (0, 1), (0, -1)])
    trial[i]["j0"] += dj * step
    trial[i]["ib0"] += di * step
    return trial, f"nudge {trial[i]['name']} ({dj * step},{di * step})"


def rotate(dies: list[dict], rng: random.Random) -> tuple[list[dict], str]:
    trial = deepcopy(dies)
    i = rng.randrange(len(trial))
    trial[i]["w"], trial[i]["h"] = trial[i]["h"], trial[i]["w"]
    return trial, f"rotate {trial[i]['name']}"


def swap_hw(dies: list[dict], cat: dict, rng: random.Random) -> tuple[list[dict], str]:
    trial = deepcopy(dies)
    choices = [i for i, d in enumerate(trial) if d["role"] in ("CPU", "GPU")]
    i = rng.choice(choices)
    cur = trial[i]
    options = [p for p in cat[cur["role"]] if p["part"] != cur["part"]]
    src = deepcopy(rng.choice(options))
    src["name"] = f"{src['part']}_{cur['slot']}"
    src["slot"] = cur["slot"]
    src["j0"] = cur["j0"]
    src["ib0"] = cur["ib0"]
    src["s"] = cur["s"]
    src["power_W"] = src["leak_W"] + src["s"] * src["dyn_W"]
    trial[i] = src
    refresh_hbm(trial)
    return trial, f"swap {cur['name']} -> {src['name']}"


def pack(dies: list[dict]) -> list[dict]:
    keys = (
        "name", "role", "part", "slot", "j0", "ib0", "w", "h", "s",
        "power_W", "dyn_W", "leak_W", "area_mm2", "grid_area_mm2", "n_macros",
        "theta_serial", "theta_parallel", "inst_per_s", "cores", "f_MHz",
        "n_cpu", "n_gpu",
    )
    return [{k: d.get(k) for k in keys} for d in dies]


def throughput(die: dict, scale: float) -> dict:
    if die["role"] == "CPU":
        return {
            "Theta_serial": die["theta_serial"] * scale,
            "Theta_parallel": die["theta_parallel"] * scale,
        }
    if die["role"] == "GPU":
        return {"inst_per_s": die["inst_per_s"] * scale}
    return {}


def hotspot_block_max(dies: list[dict], work: Path) -> dict:
    apply_power(dies)
    write_hotspot_inputs(sm.to_rects(dies), work, rows=N, cols=N)
    out = run_hotspot(HS_BIN, work)
    temps = {k: v - 273.15 for k, v in out["temperatures_K"].items() if not k.startswith("anchor")}
    return {"block_max_C": max(temps.values()), "blocks_C": temps}


def layer0_celsius(grid_path: Path) -> list[tuple[int, float]]:
    pairs = []
    grab = False
    for raw in grid_path.read_text(encoding="utf-8", errors="replace").splitlines():
        if raw.startswith("Layer 0"):
            grab = True
            continue
        if grab and raw.startswith("Layer"):
            break
        if not grab:
            continue
        parts = raw.split()
        if len(parts) >= 2:
            pairs.append((int(parts[0]), float(parts[1]) - 273.15))
    if len(pairs) != N * N:
        raise RuntimeError(f"{grid_path} 第 0 层不是 {N*N} 格，是 {len(pairs)}")
    return pairs


def write_layer_file(pairs: list[tuple[int, float]], path: Path) -> None:
    path.write_text("".join(f"{i}\t{t:.4f}\n" for i, t in pairs), encoding="utf-8")


def atsim_case(dies: list[dict], path: Path) -> None:
    apply_power(dies)
    payload = {
        "name": path.stem,
        "description": "closed_loop_multidie",
        "chips": [
            {
                "name": d["name"],
                "x_mm": d["j0"] * CW_MM,
                "y_mm": d["ib0"] * CW_MM,
                "w_mm": d["w"] * CW_MM,
                "h_mm": d["h"] * CW_MM,
                "t_um": 150.0,
                "power_W": d["power_W"],
            }
            for d in dies
        ],
        "package": {
            "tim_um": 20, "t_spreader_mm": 1.0, "t_sink_mm": 6.9,
            "s_spreader_mm": 96.0, "s_sink_mm": 96.0, "frame_mm": 72.0,
            "tim_cover": "frame",
            "t_substrate_mm": 1.0, "t_interposer_um": 100.0,
            "r_convec_K_per_W": 0.0390625, "ambient_C": 45,
            "pitch_die_mm": 2.0, "pitch_coarse_mm": 8.0,
            "nz_die": 3, "nz_spreader": 3, "nz_sink": 4, "nz_tim": 1,
            "nz_substrate": 3, "nz_interposer": 2,
            "fill_die": True, "fill_xy": "chips", "fill_material": "silicon",
        },
    }
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


def run_atsim(dies: list[dict], tag: str) -> dict:
    case_path = OUT / "atsim" / f"{tag}.json"
    out_path = OUT / "atsim" / f"{tag}_result.json"
    case_path.parent.mkdir(parents=True, exist_ok=True)
    atsim_case(dies, case_path)
    log_path = OUT / "atsim" / f"{tag}.log"
    print(f"ATSim {tag} ...", flush=True)
    proc = subprocess.run(
        [str(FENICS), str(Path(__file__).resolve()), "--atsim", str(case_path), str(out_path)],
        text=True, capture_output=True, check=False,
    )
    log_path.write_text(proc.stdout + "\n" + proc.stderr, encoding="utf-8")
    if proc.returncode != 0:
        raise RuntimeError(f"ATSim {tag} 失败 rc={proc.returncode}\n{proc.stderr[-2000:]}")
    return json.loads(out_path.read_text(encoding="utf-8"))


def scales_of(dies: list[dict]) -> dict[str, float]:
    return {d["slot"]: float(d["s"]) for d in dies if d["role"] in ("CPU", "GPU")}


def with_scales(dies: list[dict], scales: dict[str, float]) -> list[dict]:
    trial = deepcopy(dies)
    for die in trial:
        if die["role"] in ("CPU", "GPU"):
            die["s"] = round(float(scales[die["slot"]]), 1)
    apply_power(trial)
    return trial


def python_at(base_rises: dict[str, np.ndarray], base: list[dict], scales: dict[str, float]) -> dict:
    acc = np.zeros((N, N))
    for die in base:
        if die["role"] in ("CPU", "GPU"):
            p0 = die["leak_W"] + die["dyn_W"]
            p = die["leak_W"] + round(float(scales[die["slot"]]), 1) * die["dyn_W"]
        else:
            p0 = die["power_W"]
            p = die["power_W"]
        acc += base_rises[die["name"]] * (p / p0)
    return metrics(acc + hs.T_AMB)


def classify(atsim: dict) -> str:
    t = atsim["die_cell_max_C"]
    s = atsim["vm_normal_chip_max_MPa"]
    if t >= T_GATE or s >= S_GATE:
        return "outside"
    if t >= T_GATE - NEAR_T or s >= S_GATE - NEAR_S:
        return "near"
    return "inside"


def freq_search(base: list[dict], rises: dict[str, np.ndarray], kcal: float, mode: str, rng: random.Random) -> dict:
    """倍率之和最大，同时向导温度留在 90 C 下面。倍率之和不是吞吐。"""
    slots = [d["slot"] for d in base if d["role"] in ("CPU", "GPU")]
    scales = {slot: 1.0 for slot in slots}

    def guide(sc: dict[str, float]) -> dict:
        ev = python_at(rises, base, sc)
        est = T_AMB_C + kcal * (ev["Tmax_C"] - T_AMB_C)
        return {"T_python_C": ev["Tmax_C"], "T_est_C": est, "grad_max": ev["grad_max"], "ok": est < T_GATE}

    def sum_s(sc: dict[str, float]) -> float:
        return float(sum(sc.values()))

    log = []

    def snapshot(step: int, op: str, sc: dict[str, float], take: bool) -> None:
        g = guide(sc)
        log.append({
            "step": step, "op": op, "accept": take, "sum_s": sum_s(sc),
            "scales": dict(sc), **{k: g[k] for k in ("T_python_C", "T_est_C", "grad_max", "ok")},
        })

    # 温度向导已经在门内、但 ATSim 判在门外：出界的是应力。
    # 这里不再把倍率往上加，交给后面的 ATSim 回退往下降。
    if mode == "outside" and guide(scales)["ok"]:
        snapshot(0, "stress_outside_hold", scales, True)
        return {"scales": scales, "sum_s": sum_s(scales), "guide": guide(scales), "log": log, "mode": mode}

    # 门槛外先按 0.1 往下降，每一步降对向导温度最有效的那颗，直到进门或到底。
    if mode == "outside" or not guide(scales)["ok"]:
        for step in range(1, 21):
            if guide(scales)["ok"]:
                break
            best = None
            for slot in slots:
                if scales[slot] <= S_LO + 1e-9:
                    continue
                trial = dict(scales)
                trial[slot] = round(scales[slot] - S_STEP, 1)
                g = guide(trial)
                if best is None or g["T_est_C"] < best[0]:
                    best = (g["T_est_C"], slot, trial)
            if best is None:
                break
            scales = best[2]
            snapshot(step, f"drop {best[1]}", scales, True)
            print(f"  freq drop {best[1]} -> {scales[best[1]]:.1f}  T_est {best[0]:.2f}", flush=True)

    current = dict(scales)
    cur_sum = sum_s(current)
    best = dict(current)
    best_sum = cur_sum
    temp = 0.30
    start_step = len(log)
    snapshot(start_step, "freq_sa_start", current, True)
    for step in range(1, N_FREQ + 1):
        slot = rng.choice(slots)
        direction = rng.choice([S_STEP, -S_STEP])
        trial = dict(current)
        trial[slot] = round(current[slot] + direction, 1)
        if trial[slot] < S_LO - 1e-9 or trial[slot] > S_HI + 1e-9:
            log.append({"step": start_step + step, "op": f"bound {slot}", "accept": False})
            temp *= ALPHA
            continue
        g = guide(trial)
        if not g["ok"]:
            log.append({
                "step": start_step + step, "op": f"{slot} {direction:+.1f}",
                "accept": False, "reason": "over_guide", "T_est_C": g["T_est_C"],
            })
            temp *= ALPHA
            continue
        new_sum = sum_s(trial)
        delta = new_sum - cur_sum
        if delta >= -1e-12:
            take, ap = True, 1.0
        else:
            ap = float(np.exp(delta / max(temp, 1e-9)))
            take = rng.random() < ap
        rec = {
            "step": start_step + step, "op": f"{slot} {direction:+.1f}",
            "accept": take, "ap": ap, "sum_s": new_sum, "scales": trial,
            "T_python_C": g["T_python_C"], "T_est_C": g["T_est_C"], "grad_max": g["grad_max"],
        }
        if take:
            current = trial
            cur_sum = new_sum
            if new_sum > best_sum + 1e-12:
                best = dict(trial)
                best_sum = new_sum
                rec["new_best"] = True
                print(f"  freq {step:02d} BEST sum_s {new_sum:.1f}  T_est {g['T_est_C']:.2f}", flush=True)
        log.append(rec)
        temp *= ALPHA
    return {"scales": best, "sum_s": best_sum, "guide": guide(best), "log": log, "mode": mode}


def style_plots() -> None:
    font_manager.fontManager.addfont(FONT)
    plt.rcParams["font.family"] = "Noto Sans CJK JP"
    plt.rcParams["axes.unicode_minus"] = False


def draw_plots(report: dict, inner_log: list[dict], freq_log: list[dict]) -> None:
    style_plots()
    fig_dir = OUT / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    steps, t_line, g_line = [], [], []
    t_run = inner_log[0]["Tmax_C"]
    g_run = inner_log[0]["grad_max"]
    for row in inner_log:
        if row.get("accept") and "Tmax_C" in row:
            t_run = row["Tmax_C"]
            g_run = row["grad_max"]
        steps.append(row["step"])
        t_line.append(t_run)
        g_line.append(g_run)
    fig, axes = plt.subplots(2, 1, figsize=(8.2, 6.2), sharex=True)
    axes[0].plot(steps, t_line, color="#1f4e79", lw=1.6)
    axes[0].set_ylabel("Python 芯粒层峰值 (C)")
    axes[0].set_title("内环：接受后的温度和梯度", loc="left")
    axes[1].plot(steps, g_line, color="#b85c38", lw=1.6)
    axes[1].set_ylabel("Python 最大梯度 (K/mm)")
    axes[1].set_xlabel("步")
    for ax in axes:
        ax.grid(True, color="#e4e2dc", lw=0.6)
    fig.tight_layout()
    fig.savefig(fig_dir / "inner_trace.png", dpi=140)
    plt.close(fig)

    freq_steps, freq_sum, freq_t = [], [], []
    for row in freq_log:
        if "sum_s" not in row:
            continue
        freq_steps.append(row["step"])
        freq_sum.append(row["sum_s"])
        freq_t.append(row.get("T_est_C"))
    if freq_steps:
        fig, axes = plt.subplots(2, 1, figsize=(8.2, 6.2), sharex=True)
        axes[0].plot(freq_steps, freq_sum, color="#1f4e79", lw=1.6)
        axes[0].set_ylabel("各颗频率倍率之和")
        axes[0].set_title("外环频率搜索。倍率之和不是吞吐", loc="left")
        axes[1].plot(freq_steps, freq_t, color="#b85c38", lw=1.6)
        axes[1].axhline(T_GATE, color="#888780", lw=0.8, ls="--")
        axes[1].set_ylabel("向导温度 (C)")
        axes[1].set_xlabel("步")
        for ax in axes:
            ax.grid(True, color="#e4e2dc", lw=0.6)
        fig.tight_layout()
        fig.savefig(fig_dir / "freq_trace.png", dpi=140)
        plt.close(fig)

    stages = ["起点", "内环最优", "调频后"]
    keys = ("init", "best", "freq")
    fig, axes = plt.subplots(1, 3, figsize=(9.6, 4.2))
    series = (
        ("Python 峰值", [report[k]["python_Tmax_C"] for k in keys], "#1f4e79"),
        ("HotSpot 块平均峰值", [report[k]["hotspot_block_max_C"] for k in keys], "#3d7a5a"),
        ("ATSim 芯粒层峰值", [report[k]["atsim"]["die_cell_max_C"] for k in keys], "#b85c38"),
    )
    x = np.arange(len(stages))
    width = 0.25
    ax = axes[0]
    for i, (label, vals, color) in enumerate(series):
        ax.bar(x + (i - 1) * width, vals, width, label=label, color=color)
    ax.axhline(T_GATE, color="#888780", lw=0.8, ls="--")
    ax.set_xticks(x, stages)
    ax.set_ylabel("温度 (C)")
    ax.set_title("三种求解器，不要混读", loc="left", fontsize=10)
    ax.legend(fontsize=7, frameon=False)
    ax.set_ylim(0, max(v for _, vals, _ in series for v in vals) * 1.15)

    stress = [report[k]["atsim"]["vm_normal_chip_max_MPa"] for k in keys]
    axes[1].bar(stages, stress, color="#1f4e79")
    axes[1].axhline(S_GATE, color="#b85c38", lw=0.8, ls="--")
    axes[1].set_ylabel("正常芯片内部 von Mises (MPa)")
    axes[1].set_title("应力门槛 200 MPa", loc="left", fontsize=10)

    grads = (
        ("Python", [report[k]["python_grad_K_per_mm"] for k in keys], "#1f4e79"),
        ("ATSim", [report[k]["atsim"]["grad_max_K_per_mm"] for k in keys], "#b85c38"),
    )
    for i, (label, vals, color) in enumerate(grads):
        axes[2].bar(x + (i - 0.5) * width, vals, width, label=label, color=color)
    axes[2].set_xticks(x, stages)
    axes[2].set_ylabel("最大梯度 (K/mm)")
    axes[2].set_title("梯度，两套网格", loc="left", fontsize=10)
    axes[2].legend(fontsize=7, frameon=False)
    fig.tight_layout()
    fig.savefig(fig_dir / "stage_compare.png", dpi=140)
    plt.close(fig)

    chips = report["tunable"]
    labels = [c["name_after"] for c in chips]
    x = np.arange(len(labels))
    fig, ax = plt.subplots(figsize=(8.2, 4.2))
    ax.bar(x - 0.18, [c["s_before"] for c in chips], 0.36, label="调频前", color="#c5c1b9")
    ax.bar(x + 0.18, [c["s_after"] for c in chips], 0.36, label="调频后", color="#1f4e79")
    ax.set_xticks(x, labels, rotation=20, ha="right")
    ax.set_ylabel("频率倍率")
    ax.set_title("每颗 CPU / GPU 单独的倍率，步长 0.1", loc="left")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(fig_dir / "freq_scales.png", dpi=140)
    plt.close(fig)

    cpus = [c for c in chips if c["role"] == "CPU"]
    gpus = [c for c in chips if c["role"] == "GPU"]
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.2))
    x = np.arange(len(cpus))
    axes[0].bar(x - 0.18, [c["theta_serial_before"] for c in cpus], 0.36, label="调频前", color="#c5c1b9")
    axes[0].bar(x + 0.18, [c["theta_serial_after"] for c in cpus], 0.36, label="调频后", color="#1f4e79")
    axes[0].set_xticks(x, [c["name_after"] for c in cpus], rotation=15, ha="right")
    axes[0].set_ylabel("串行吞吐，Penryn 单线程 = 1")
    axes[0].set_title("CPU 串行", loc="left")
    axes[0].legend(frameon=False, fontsize=8)
    x = np.arange(len(gpus))
    axes[1].bar(x - 0.18, [c["inst_per_s_before"] / 1e9 for c in gpus], 0.36, label="调频前", color="#c5c1b9")
    axes[1].bar(x + 0.18, [c["inst_per_s_after"] / 1e9 for c in gpus], 0.36, label="调频后", color="#3d7a5a")
    axes[1].set_xticks(x, [c["name_after"] for c in gpus], rotation=15, ha="right")
    axes[1].set_ylabel("GPU 指令每秒 (1e9)")
    axes[1].set_title("GPU，不和 CPU 相加", loc="left")
    axes[1].legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(fig_dir / "throughput.png", dpi=140)
    plt.close(fig)

    if cpus:
        fig, ax = plt.subplots(figsize=(6.4, 4.0))
        x = np.arange(len(cpus))
        ax.bar(x - 0.18, [c["theta_parallel_before"] for c in cpus], 0.36, label="调频前", color="#c5c1b9")
        ax.bar(x + 0.18, [c["theta_parallel_after"] for c in cpus], 0.36, label="调频后", color="#1f4e79")
        ax.set_xticks(x, [c["name_after"] for c in cpus])
        ax.set_ylabel("并行吞吐，Penryn 单线程 = 1")
        ax.set_title("CPU 并行，核数乘主频", loc="left")
        ax.legend(frameon=False)
        fig.tight_layout()
        fig.savefig(fig_dir / "throughput_cpu_parallel.png", dpi=140)
        plt.close(fig)


def write_svgs() -> None:
    maps = []
    for tag in ("init", "best", "freq"):
        pairs = layer0_celsius(OUT / f"hotspot_{tag}" / "grid.steady")
        layer = OUT / f"hotspot_{tag}" / "layer0_C.txt"
        write_layer_file(pairs, layer)
        maps.append((tag, pairs, layer))
    vmin = min(t for _, pairs, _ in maps for _, t in pairs)
    vmax = max(t for _, pairs, _ in maps for _, t in pairs)
    fig_dir = OUT / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)
    for tag, _, layer in maps:
        proc = subprocess.run(
            ["perl", str(SVG_PL), str(OUT / f"hotspot_{tag}" / "floorplan.flp"), str(layer),
             "64", "64", f"{vmin:.4f}", f"{vmax:.4f}"],
            text=True, capture_output=True, check=False,
        )
        if proc.returncode != 0:
            raise RuntimeError(proc.stderr[-2000:])
        (fig_dir / f"hotspot_{tag}.svg").write_text(proc.stdout, encoding="utf-8")
    (fig_dir / "svg_scale.json").write_text(
        json.dumps({"min_C": vmin, "max_C": vmax, "note": "三张 SVG 共用这一档。数值是 HotSpot 第 0 层从 K 换成 C。"}, indent=2)
        + "\n",
        encoding="utf-8",
    )


def tunable_rows(before: list[dict], after: list[dict]) -> list[dict]:
    after_by = {d["slot"]: d for d in after if d["role"] in ("CPU", "GPU")}
    rows = []
    for die in before:
        if die["role"] not in ("CPU", "GPU"):
            continue
        nxt = after_by[die["slot"]]
        row = {
            "slot": die["slot"],
            "role": die["role"],
            "name_before": die["name"],
            "name_after": nxt["name"],
            "s_before": 1.0,
            "s_after": nxt["s"],
            "power_before_W": die["leak_W"] + die["dyn_W"],
            "power_after_W": nxt["power_W"],
            "dyn_W": nxt["dyn_W"],
            "leak_W": nxt["leak_W"],
        }
        b = throughput(die, 1.0)
        a = throughput(nxt, nxt["s"])
        if die["role"] == "CPU":
            row["theta_serial_before"] = b["Theta_serial"]
            row["theta_serial_after"] = a["Theta_serial"]
            row["theta_parallel_before"] = b["Theta_parallel"]
            row["theta_parallel_after"] = a["Theta_parallel"]
        else:
            row["inst_per_s_before"] = b["inst_per_s"]
            row["inst_per_s_after"] = a["inst_per_s"]
        rows.append(row)
    return rows


def stage_block(ev: dict, hotspot: dict, atsim: dict) -> dict:
    return {
        "python_Tmax_C": ev["Tmax_C"],
        "python_grad_K_per_mm": ev["grad_max"],
        "hotspot_block_max_C": hotspot["block_max_C"],
        "atsim": {
            "die_cell_max_C": atsim["die_cell_max_C"],
            "grad_max_K_per_mm": atsim["grad_max_K_per_mm"],
            "vm_normal_chip_max_MPa": atsim["vm_normal_chip_max_MPa"],
            "vm_normal_chip_name": atsim["vm_normal_chip_name"],
        },
    }


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "layouts").mkdir(exist_ok=True)
    (OUT / "fields").mkdir(exist_ok=True)
    cat = catalog()
    start = place(cat)
    if not sm.legal(start):
        raise SystemExit("起点不合法")
    k = np.full((N, N), hs.K_SI)
    mat, g_amb = hs.build_operator(k)
    fac = splu(mat.tocsc())
    cache: dict = {}
    summed = superpose(start, fac, g_amb, cache)
    joint = hw.solve(k, hw.power_map(start), fac, g_amb)[0]
    sup_err = float(np.max(np.abs(summed - joint)))
    print(f"superposition max abs {sup_err:.3e} K", flush=True)
    if sup_err > 0.05:
        raise SystemExit("均匀硅上的场叠加误差过大，停止")

    rng_m = random.Random(SEED)
    rng_a = random.Random(SEED + 7919)
    current = deepcopy(start)
    cur_t = summed
    cur = metrics(cur_t)
    scales = {"dT": max(cur["Tmax_C"] - T_AMB_C, 1e-6), "g0": max(cur["grad_max"], 1e-9)}
    j_cur = score(cur, scales)
    best = deepcopy(current)
    best_ev = cur
    best_j = j_cur
    temp = 0.20
    log = [{
        "step": 0, "op": "init", "accept": True, "J": j_cur,
        "Tmax_C": cur["Tmax_C"], "grad_max": cur["grad_max"],
        "hardware": [d["name"] for d in current],
        "hbm_W": next(d["power_W"] for d in current if d["role"] == "HBM"),
    }]
    print(f"init J {j_cur:.3f}  T {cur['Tmax_C']:.2f}  grad {cur['grad_max']:.3f}", flush=True)
    for step in range(1, N_INNER + 1):
        roll = rng_m.random()
        if roll < 0.60:
            trial, op = nudge(current, rng_m)
            move = "nudge"
        elif roll < 0.80:
            trial, op = rotate(current, rng_m)
            move = "rotate"
        else:
            trial, op = swap_hw(current, cat, rng_m)
            move = "swap"
        if not sm.legal(trial):
            log.append({"step": step, "op": op, "move": move, "accept": False, "reason": "illegal"})
            temp *= ALPHA
            continue
        ev = metrics(superpose(trial, fac, g_amb, cache))
        j_new = score(ev, scales)
        take, ap = sm.accept(j_cur, j_new, temp, rng_a)
        hbm_w = next(d["power_W"] for d in trial if d["role"] == "HBM")
        rec = {
            "step": step, "op": op, "move": move, "accept": take, "ap": ap, "J": j_new,
            "Tmax_C": ev["Tmax_C"], "grad_max": ev["grad_max"], "anneal_T": temp,
            "hardware": [d["name"] for d in trial], "hbm_W": hbm_w,
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
                print(f"  {step:02d} BEST {op}  J {j_new:.3f}  T {ev['Tmax_C']:.2f}", flush=True)
        log.append(rec)
        temp *= ALPHA

    best_sum = superpose(best, fac, g_amb, cache)
    best_joint = hw.solve(k, hw.power_map(best), fac, g_amb)[0]
    best_err = float(np.max(np.abs(best_sum - best_joint)))
    np.save(OUT / "fields" / "init.npy", summed)
    np.save(OUT / "fields" / "best.npy", best_sum)
    (OUT / "layouts" / "init.json").write_text(json.dumps(pack(start), indent=2) + "\n", encoding="utf-8")
    (OUT / "layouts" / "best.json").write_text(json.dumps(pack(best), indent=2) + "\n", encoding="utf-8")
    (OUT / "trace.json").write_text(json.dumps(log, indent=2) + "\n", encoding="utf-8")
    print(f"best superposition max abs {best_err:.3e} K", flush=True)

    hs_init = hotspot_block_max(start, OUT / "hotspot_init")
    hs_best = hotspot_block_max(best, OUT / "hotspot_best")
    kcal = (hs_best["block_max_C"] - T_AMB_C) / max(best_ev["Tmax_C"] - T_AMB_C, 1e-6)
    print(f"HotSpot init {hs_init['block_max_C']:.2f}  best {hs_best['block_max_C']:.2f}  kcal {kcal:.4f}", flush=True)

    at_best = run_atsim(best, "best_s1")
    mode = classify(at_best)
    print(f"outer mode {mode}", flush=True)
    rises = {}
    for die in best:
        rises[die["name"]] = rise_of(die, fac, g_amb, cache)
    freq = freq_search(best, rises, kcal, mode, random.Random(SEED + 17))
    tuned = with_scales(best, freq["scales"])
    # 向导进了门之后，ATSim 仍可能越出门槛。最多再降 4 次，每次重算 ATSim。
    at_freq = None
    walks = []
    for attempt in range(7):
        at_freq = run_atsim(tuned, f"freq_{attempt}")
        walks.append({
            "attempt": attempt,
            "scales": scales_of(tuned),
            "die_cell_max_C": at_freq["die_cell_max_C"],
            "vm_normal_chip_max_MPa": at_freq["vm_normal_chip_max_MPa"],
        })
        if at_freq["die_cell_max_C"] < T_GATE and at_freq["vm_normal_chip_max_MPa"] < S_GATE:
            break
        movable = [d for d in tuned if d["role"] in ("CPU", "GPU") and d["s"] > S_LO + 1e-9]
        if not movable:
            break
        victim = max(movable, key=lambda d: d["dyn_W"])
        victim["s"] = round(victim["s"] - S_STEP, 1)
        apply_power(tuned)
        print(f"  walkback {victim['name']} s {victim['s']:.1f}", flush=True)
    freq_ev = metrics(superpose(tuned, fac, g_amb, cache))
    hs_freq = hotspot_block_max(tuned, OUT / "hotspot_freq")
    at_init = run_atsim(start, "init_s1")
    (OUT / "layouts" / "freq.json").write_text(json.dumps(pack(tuned), indent=2) + "\n", encoding="utf-8")
    (OUT / "freq_trace.json").write_text(json.dumps(freq["log"], indent=2) + "\n", encoding="utf-8")
    np.save(OUT / "fields" / "freq.npy", freq_ev["T"])

    rows = tunable_rows(best, tuned)
    report = {
        "fill": "silicon",
        "k_Si": hs.K_SI,
        "package": "72 mm 芯粒框，96 mm 热沉，r_convec 0.0390625",
        "inner_cost": "(T_python-45)/(T_init-45) + grad/grad_init",
        "inner_steps": N_INNER,
        "seed": SEED,
        "superposition_max_abs_init_K": sup_err,
        "superposition_max_abs_best_K": best_err,
        "superposition_note": "均匀硅、冻结材料图，只开一颗的功耗再相加。这个数是和方法误差，不是 HotSpot 外伸造成的温差。",
        "hbm_formula": "每通道 0.088 + n_cpu*(1.035-0.088) + n_gpu*(0.844-0.088)，再乘 16。只认颗数。",
        "hbm_W_init": next(d["power_W"] for d in start if d["role"] == "HBM"),
        "hbm_W_best": next(d["power_W"] for d in best if d["role"] == "HBM"),
        "kcal_python_rise_to_hotspot_block": kcal,
        "kcal_note": "把 Python 芯粒层峰值的温升，按内环最优那一点的 HotSpot 块平均峰值缩放。只用于频率搜索的向导，不是 HotSpot，也不是 ATSim。",
        "gates": {"T_C": T_GATE, "vm_MPa": S_GATE, "near_T_C": NEAR_T, "near_vm_MPa": NEAR_S},
        "near_note": "附近 = 没越出门槛，但温度距 90 C 不足 5 C，或应力距 200 MPa 不足 20 MPa。这是这一次跑的带宽。",
        "freq_mode": mode,
        "freq_box": {"step": S_STEP, "lo": S_LO, "hi": S_HI},
        "freq_objective": "各颗 CPU/GPU 倍率之和最大，且向导温度 < 90 C。倍率之和不是吞吐。漏电不乘倍率。",
        "walkback": walks,
        "init": stage_block(metrics(summed), hs_init, at_init),
        "best": stage_block(best_ev, hs_best, at_best),
        "freq": stage_block(freq_ev, hs_freq, at_freq),
        "guide_at_freq": freq["guide"],
        "tunable": rows,
        "sum_s_before": float(sum(1.0 for _ in rows)),
        "sum_s_after": float(sum(r["s_after"] for r in rows)),
        "accepted_inner": sum(1 for r in log if r.get("accept") and r["step"] > 0),
        "hardware_init": [d["name"] for d in start],
        "hardware_best": [d["name"] for d in best],
        "note": (
            "没有改生产 SA。Python、HotSpot、ATSim 的温度分列。"
            "ATSim 网格是 2 mm，比对齐检查的 1.125 mm 粗。"
            "正常芯片内部应力才是门槛，整包最大应力不是。"
            "CIM 不换、不调频。HBM 不调频。"
        ),
    }
    (OUT / "report.json").write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    write_svgs()
    draw_plots(report, log, freq["log"])
    print(f"写出 {OUT / 'report.json'}", flush=True)


if __name__ == "__main__":
    main()
