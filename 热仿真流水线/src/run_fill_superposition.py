#!/usr/bin/env python3
"""填充物不再是硅时，场叠加还差多少。

稳态热传导是 div(k grad T) + q = 0。k 在芯粒层逐格取值：芯粒足迹上是硅，
其余是填充物。离散和 HotSpot 非 detailed-3D 网格一致（方形网格上 rx、ry 相同）：

    横向相邻两格     R = rx(k_a)/2 + rx(k_b)/2
    层 n 到层 n+1    R = rz(层 n)          （整层热阻算在下层）
    热沉到环境       R = t_sink/(k A) + r_convec * S_sink^2 / A

k 均匀时，网格内部就是 HotSpot 芯片网格的模板。72 mm 以外的均热片外伸
用 8 个集总节点补上，和 论文证据/填充场叠加/复现/ring_solver.py 同一套：
均热片北南西东各一条，每条上面挂一个热沉节点。这套封装均热片和热沉都是
96 mm，热沉再往外的 4 个节点面积是 0，HotSpot 在那里会除零，所以不建。
外伸条宽 12 mm、长 84 mm。叠加误差仍是在同一个算子里比较。

两种“场叠加”：

- naive：每颗芯片单独求解，此时只有它自己的足迹是硅，别的芯片占的位置仍是
  填充物。这就是现有缓存 (芯片, x, y, w, h) 的做法。k 均匀时填充物也是硅，
  算子不变，相加精确；填充物一变，每张场的算子都和联立解不同。
- frozen：材料图冻结为全部芯片都是硅，只轮流打开功耗。算子不变，相加仍精确。
  用来证明误差来自材料图变了，不是功耗相加本身。
"""

from __future__ import annotations

import json
import shutil
import sys
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve

ROOT = Path("/home/czh/任务仿真产生功耗文件工具")
SRC = ROOT / "热仿真流水线" / "src"
sys.path.insert(0, str(SRC))

from hotspot_io import run_hotspot, write_hotspot_inputs  # noqa: E402
from geometry import Rect  # noqa: E402

HS = ROOT / "HotSpot-master" / "hotspot"
OUT = ROOT / "热仿真流水线" / "runs" / "fill_superposition"

# 与 hotspot_io.py 不传 package 时的封装一致。芯粒层由两个锚点钉成 72 mm。
N = 64
SIDE = 0.072
CW = SIDE / N
T_AMB = 318.15
K_SI = 130.0
# 二氧化硅填充的常规环氧模塑料，文献里常见约 0.8 W/m·K（未增强牌号多在 0.5–1）。
K_FILL = 0.8
# HotSpot package.h 里的空气导热，当作下界，不是推荐的封装填充。
K_AIR = 0.028

LAYERS = (
    # name, thickness_m, k
    ("chip", 0.00015, None),  # k 逐格
    ("tim", 2.0e-05, 4.0),
    ("spreader", 0.001, 400.0),
    ("sink", 0.0069, 400.0),
)
R_CONVEC = 0.1 * (0.06 * 0.06) / (0.096 * 0.096)
S_SINK = 0.096
S_SPREADER = 0.096
NL = len(LAYERS)
N_GRID = NL * N * N
# 外伸节点顺序和 ring_solver.py 一致：均热片北南西东，再是对应的热沉节点。
SP_N, SP_S, SP_E, SP_W = 0, 1, 2, 3
SINK_C_N, SINK_C_S, SINK_C_E, SINK_C_W = 4, 5, 6, 7
N_EXTRA = 8
N_UNK = N_GRID + N_EXTRA


def _overhang() -> dict:
    """一边的外伸条。均热片和热沉同宽，只有均热片这一条有面积。"""
    span = (S_SPREADER - SIDE) / 2.0
    length = (S_SPREADER + SIDE) / 2.0
    t_sp = LAYERS[2][1]
    k_cu = LAYERS[2][2]
    g_from_edge = 1.0 / ((span / 2.0) / (k_cu * length * t_sp))
    g_up = 1.0 / (t_sp / (k_cu * span * length))
    area = span * length
    t_sink = LAYERS[3][1]
    g_to_air = 1.0 / (t_sink / (k_cu * area) + R_CONVEC * (S_SINK * S_SINK) / area)
    return {"g_edge": g_from_edge / N, "g_up": g_up, "g_to_air": g_to_air}


_OVER = _overhang()
G_EDGE = _OVER["g_edge"]
G_UP = _OVER["g_up"]
G_TO_AIR = _OVER["g_to_air"]

# 网格对齐的三颗芯片。j 是 x（列），ib 是从下往上的行，和 HotSpot 的 y 一致。
# 间隙 4.5 mm，三颗互不重叠。
DIES = (
    {"name": "cpu", "j0": 16, "j1": 32, "ib0": 24, "ib1": 40, "power_W": 80.0},
    {"name": "gpu", "j0": 36, "j1": 48, "ib0": 28, "ib1": 40, "power_W": 40.0},
    {"name": "cim", "j0": 24, "j1": 32, "ib0": 12, "ib1": 20, "power_W": 15.0},
)


def rows_of(die: dict) -> tuple[int, int]:
    """HotSpot 行号：第 0 行在 y 最大处。返回半开区间 [i0, i1)。"""
    i_top = N - die["ib1"]
    i_bot = N - die["ib0"]
    return i_top, i_bot


def k_map(present: list[dict], k_fill: float) -> np.ndarray:
    k = np.full((N, N), k_fill, dtype=float)
    for d in present:
        i0, i1 = rows_of(d)
        k[i0:i1, d["j0"]:d["j1"]] = K_SI
    return k


def power_map(powered: list[dict]) -> np.ndarray:
    p = np.zeros((N, N), dtype=float)
    for d in powered:
        i0, i1 = rows_of(d)
        ncell = (i1 - i0) * (d["j1"] - d["j0"])
        p[i0:i1, d["j0"]:d["j1"]] = d["power_W"] / ncell
    return p


def _rx(k: float, thickness: float) -> float:
    """HotSpot getr(k, cw, ch*t)。方格上和 ry 相同。"""
    return CW / (k * CW * thickness)


def _rz(k: float, thickness: float) -> float:
    return thickness / (k * CW * CW)


def idx_of(layer: int, i: int, j: int) -> int:
    return (layer * N + i) * N + j


def build_operator(k_chip: np.ndarray) -> tuple[sparse.csr_matrix, float]:
    """导热矩阵，以及热沉每格到环境的电导。功耗不在矩阵里。

    矩阵比网格多 8 个外伸节点。solve_state 会把这 8 个未知数去掉，
    返回的仍是 (NL, N, N)。
    """
    nunk = N_UNK
    acc = [[] for _ in range(nunk)]

    def add(a: int, b: int, g: float) -> None:
        acc[a].append((b, -g))
        acc[b].append((a, -g))
        acc[a].append((a, g))
        acc[b].append((b, g))

    t_chip = LAYERS[0][1]
    for i in range(N):
        for j in range(N):
            if j + 1 < N:
                r = 0.5 * _rx(k_chip[i, j], t_chip) + 0.5 * _rx(k_chip[i, j + 1], t_chip)
                add(idx_of(0, i, j), idx_of(0, i, j + 1), 1.0 / r)
            if i + 1 < N:
                r = 0.5 * _rx(k_chip[i, j], t_chip) + 0.5 * _rx(k_chip[i + 1, j], t_chip)
                add(idx_of(0, i, j), idx_of(0, i + 1, j), 1.0 / r)

    for layer in range(1, NL):
        t = LAYERS[layer][1]
        k = LAYERS[layer][2]
        g = 1.0 / _rx(k, t)
        for i in range(N):
            for j in range(N):
                if j + 1 < N:
                    add(idx_of(layer, i, j), idx_of(layer, i, j + 1), g)
                if i + 1 < N:
                    add(idx_of(layer, i, j), idx_of(layer, i + 1, j), g)

    for i in range(N):
        for j in range(N):
            g = 1.0 / _rz(k_chip[i, j], t_chip)
            add(idx_of(0, i, j), idx_of(1, i, j), g)
    for layer in range(1, NL - 1):
        t = LAYERS[layer][1]
        k = LAYERS[layer][2]
        g = 1.0 / _rz(k, t)
        for i in range(N):
            for j in range(N):
                add(idx_of(layer, i, j), idx_of(layer + 1, i, j), g)

    t_sink = LAYERS[3][1]
    k_sink = LAYERS[3][2]
    area = CW * CW
    r_sink = _rz(k_sink, t_sink) + R_CONVEC * (S_SINK * S_SINK) / area
    g_amb = 1.0 / r_sink
    for i in range(N):
        for j in range(N):
            a = idx_of(3, i, j)
            acc[a].append((a, g_amb))

    # 均热片边界格接到外伸条，外伸条再向上到热沉节点，热沉节点对空气。
    # 角上的格子同时属于两条边，会接到两个外伸节点。
    base = N_GRID
    for cells, node in (
        ([idx_of(2, 0, j) for j in range(N)], SP_N),
        ([idx_of(2, N - 1, j) for j in range(N)], SP_S),
        ([idx_of(2, i, 0) for i in range(N)], SP_W),
        ([idx_of(2, i, N - 1) for i in range(N)], SP_E),
    ):
        for cell in cells:
            add(cell, base + node, G_EDGE)
    for a, b in (
        (SP_N, SINK_C_N), (SP_S, SINK_C_S),
        (SP_W, SINK_C_W), (SP_E, SINK_C_E),
    ):
        add(base + a, base + b, G_UP)
    for node in (SINK_C_N, SINK_C_S, SINK_C_W, SINK_C_E):
        a = base + node
        acc[a].append((a, G_TO_AIR))

    rows, cols, data = [], [], []
    for a, terms in enumerate(acc):
        merged: dict[int, float] = {}
        for b, v in terms:
            merged[b] = merged.get(b, 0.0) + v
        for b, v in merged.items():
            rows.append(a)
            cols.append(b)
            data.append(v)
    mat = sparse.csr_matrix((data, (rows, cols)), shape=(nunk, nunk))
    return mat, g_amb


def solve_state(k_chip: np.ndarray, power: np.ndarray,
                 factor=None, g_amb: float | None = None) -> np.ndarray:
    """四层温度，形状 (NL, N, N)，单位 K。

    factor 是 scipy.sparse.linalg.splu(build_operator(...))。传入时必须同时
    给出那个算子的 g_amb，矩阵里已经含对环境的电导，这里只补右端项。
    """
    if factor is None:
        mat, g_amb = build_operator(k_chip)
    if g_amb is None:
        raise ValueError("g_amb is required")
    rhs = np.zeros(N_UNK)
    for i in range(N):
        for j in range(N):
            rhs[idx_of(3, i, j)] = g_amb * T_AMB
            rhs[idx_of(0, i, j)] += power[i, j]
    for node in (SINK_C_N, SINK_C_S, SINK_C_W, SINK_C_E):
        rhs[N_GRID + node] = G_TO_AIR * T_AMB
    temp = spsolve(mat, rhs) if factor is None else factor.solve(rhs)
    return temp[:N_GRID].reshape(NL, N, N)


def solve(k_chip: np.ndarray, power: np.ndarray) -> np.ndarray:
    """返回芯粒层温度，形状 (N, N)，单位 K。第 0 行是 y 最大的一行。"""
    return solve_state(k_chip, power)[0]


def superpose(k_fill: float, mode: str) -> tuple[np.ndarray, np.ndarray]:
    """mode = naive 或 frozen。返回 (联立真解, 叠加场)，都是芯粒层 K。"""
    truth = solve(k_map(list(DIES), k_fill), power_map(list(DIES)))
    parts = []
    for n in range(len(DIES)):
        present = list(DIES) if mode == "frozen" else [DIES[n]]
        powered = [DIES[n]]
        field = solve(k_map(present, k_fill), power_map(powered))
        parts.append(field - T_AMB)
    summed = np.sum(parts, axis=0) + T_AMB
    return truth, summed


def stats(truth: np.ndarray, summed: np.ndarray) -> dict:
    err = summed - truth
    return {
        "Tmax_truth_C": float(truth.max() - 273.15),
        "Tmax_sum_C": float(summed.max() - 273.15),
        "dTmax_K": float(summed.max() - truth.max()),
        "max_abs_K": float(np.abs(err).max()),
        "rmse_K": float(np.sqrt(np.mean(err ** 2))),
        "rise_truth_K": float(truth.max() - T_AMB),
    }


def to_rects() -> list[Rect]:
    rects = []
    for d in DIES:
        rects.append(
            Rect(
                name=d["name"],
                x=d["j0"] * CW,
                y=d["ib0"] * CW,
                w=(d["j1"] - d["j0"]) * CW,
                h=(d["ib1"] - d["ib0"]) * CW,
                power_W=d["power_W"],
                area_mm2=(d["j1"] - d["j0"]) * (d["ib1"] - d["ib0"]) * (CW * 1e3) ** 2,
                kind=d["name"],
            )
        )
    return rects


def hotspot_uniform(work: Path) -> np.ndarray:
    if work.exists():
        shutil.rmtree(work)
    write_hotspot_inputs(to_rects(), work, rows=N, cols=N)
    run_hotspot(HS, work)
    vals = []
    for line in (work / "grid.steady").read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("Layer"):
            continue
        vals.append(float(line.split()[-1]))
        if len(vals) == N * N:
            break
    return np.array(vals).reshape(N, N)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    report: dict = {
        "grid": N,
        "side_m": SIDE,
        "k_silicon": K_SI,
        "k_fill_epoxy": K_FILL,
        "k_fill_air": K_AIR,
        "ambient_K": T_AMB,
        "r_convec": R_CONVEC,
        "dies": [],
        "note": (
            "绝对温度不含 HotSpot 的 spreader 外伸环，只和同一次离散里的真解比。"
            "naive 是每颗芯片单独求解（其余位置是填充物）再把温升相加。"
            "frozen 是材料图固定为全部芯片都是硅，只分摊功耗。"
        ),
    }
    for d in DIES:
        i0, i1 = rows_of(d)
        report["dies"].append(
            {
                "name": d["name"],
                "x_mm": d["j0"] * CW * 1e3,
                "y_mm": d["ib0"] * CW * 1e3,
                "w_mm": (d["j1"] - d["j0"]) * CW * 1e3,
                "h_mm": (d["ib1"] - d["ib0"]) * CW * 1e3,
                "power_W": d["power_W"],
                "rows": [i0, i1],
                "cols": [d["j0"], d["j1"]],
            }
        )

    cases = []
    for name, k_fill, mode in (
        ("silicon_naive", K_SI, "naive"),
        ("epoxy_naive", K_FILL, "naive"),
        ("epoxy_frozen", K_FILL, "frozen"),
        ("air_naive", K_AIR, "naive"),
    ):
        truth, summed = superpose(k_fill, mode)
        st = stats(truth, summed)
        st["case"] = name
        st["k_fill"] = k_fill
        st["mode"] = mode
        cases.append(st)
        np.save(OUT / f"{name}_truth.npy", truth)
        np.save(OUT / f"{name}_sum.npy", summed)
        print(
            f"{name:16s} 真解峰值 {st['Tmax_truth_C']:.3f} °C  "
            f"叠加 {st['Tmax_sum_C']:.3f} °C  "
            f"峰值差 {st['dTmax_K']:+.4f} K  "
            f"逐格最大 {st['max_abs_K']:.4f} K",
            flush=True,
        )
    report["cases"] = cases

    print("HotSpot 整层硅对照 …", flush=True)
    hs = hotspot_uniform(OUT / "hotspot_silicon")
    py = np.load(OUT / "silicon_naive_truth.npy")
    report["hotspot_check"] = {
        "Tmax_hotspot_C": float(hs.max() - 273.15),
        "Tmax_python_C": float(py.max() - 273.15),
        "dTmax_python_minus_hotspot_K": float(py.max() - hs.max()),
        "max_abs_K": float(np.abs(py - hs).max()),
        "rmse_K": float(np.sqrt(np.mean((py - hs) ** 2))),
    }
    np.save(OUT / "hotspot_silicon.npy", hs)
    h = report["hotspot_check"]
    print(
        f"HotSpot {h['Tmax_hotspot_C']:.3f} °C  "
        f"本求解器 {h['Tmax_python_C']:.3f} °C  "
        f"逐格最大 {h['max_abs_K']:.3f} K",
        flush=True,
    )

    (OUT / "report.json").write_text(
        json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    _plot(report)
    print(f"写出 {OUT / 'report.json'}", flush=True)


def _plot(report: dict) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    truth = np.load(OUT / "epoxy_naive_truth.npy") - 273.15
    summed = np.load(OUT / "epoxy_naive_sum.npy") - 273.15
    err = summed - truth
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.8))
    extent = [0, SIDE * 1e3, 0, SIDE * 1e3]
    vmax = max(truth.max(), summed.max())
    vmin = min(truth.min(), summed.min())
    # 数组第 0 行是 y 最大处，origin=lower 前先翻过来，避免热图上下颠倒。
    truth, summed, err = truth[::-1], summed[::-1], err[::-1]
    for ax, data, title, cmap, v0, v1 in (
        (axes[0], truth, "epoxy truth", "inferno", vmin, vmax),
        (axes[1], summed, "naive superposition", "inferno", vmin, vmax),
        (axes[2], err, "sum − truth (K)", "coolwarm", None, None),
    ):
        kw = {"origin": "lower", "extent": extent, "cmap": cmap}
        if v0 is not None:
            kw["vmin"] = v0
            kw["vmax"] = v1
        im = ax.imshow(data, **kw)
        ax.set_title(title)
        ax.set_xlabel("x (mm)")
        ax.set_ylabel("y (mm)")
        fig.colorbar(im, ax=ax, fraction=0.046)
    fig.tight_layout()
    fig.savefig(OUT / "epoxy_naive_error.png", dpi=120)
    plt.close(fig)


if __name__ == "__main__":
    main()
