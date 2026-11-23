#!/usr/bin/env python3
"""芯片变大变小时，温度场错在哪，以及怎么补。

上一份实验已经说明：在这套封装里，芯粒层填充物换成环氧或空气，温度几乎不动。
横向热走的是 1 mm 铜散热片，不走 0.15 mm 的芯粒层。所以面积变化里，
“填充被硅替换”可以忽略。真正变的是功耗图：多少瓦、从哪些格子送进散热片。

同一块芯片、同一个任务，面积和性能可以一起变。这里用两条明确的功耗假设
把“可能变”的两端摆出来，不把它们说成 McPAT 的结果：

- 等功耗：总功率保持 80 W，大芯片只是把同样的瓦数铺开。
- 等密度：功耗密度不变，总功率随面积成正比。运行时更大的芯片有更多电路在翻转。

基准是正中 16×16 格（18 mm）的 80 W。算子用整层硅，只分解一次。
面积或功耗变了之后，温升的修正是

    T_new = T_old + M^{-1} (p_new - p_old)

右端项只在功耗变过的格子上非零。材料图没有改，别的芯片已经缓存的场不用重算。
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy.sparse.linalg import splu

ROOT = Path("/home/czh/任务仿真产生功耗文件工具")
SRC = ROOT / "热仿真流水线" / "src"
sys.path.insert(0, str(SRC))

import run_fill_superposition as hs  # noqa: E402

OUT = ROOT / "热仿真流水线" / "runs" / "area_correction"
SIDES = (8, 12, 16, 20, 24, 32)
SIDE0 = 16
P0 = 80.0
T_CHIP = hs.LAYERS[0][1]


def origin(side: int) -> int:
    return (hs.N - side) // 2


def power_of(side: int, law: str) -> np.ndarray:
    p = np.zeros((hs.N, hs.N))
    a = origin(side)
    if law == "fixed_W":
        total = P0
        p[a : a + side, a : a + side] = total / side**2
    elif law == "prop_area":
        total = P0 * (side / SIDE0) ** 2
        p[a : a + side, a : a + side] = total / side**2
    elif law == "same_map":
        # 功耗图钉在基准 16 格上。更大的芯片只是外圈从填充变成硅，不额外发热。
        b = origin(SIDE0)
        p[b : b + SIDE0, b : b + SIDE0] = P0 / SIDE0**2
    else:
        raise ValueError(law)
    return p


def upward_watts(state: np.ndarray, k_chip: np.ndarray) -> np.ndarray:
    """每个芯粒格向上进入 TIM 的功率，W。稳态下全场之和等于芯片总功耗。"""
    g = k_chip * hs.CW * hs.CW / T_CHIP
    return g * (state[0] - state[1])


def pack(state: np.ndarray, power: np.ndarray, k_chip: np.ndarray, side: int) -> dict:
    chip = state[0]
    flux = upward_watts(state, k_chip)
    a = origin(side)
    inside = flux[a : a + side, a : a + side].sum()
    total = float(flux.sum())
    cell = hs.CW * hs.CW
    return {
        "power_W": float(power.sum()),
        "Tmax_C": float(chip.max() - 273.15),
        "rise_K": float(chip.max() - hs.T_AMB),
        "flux_sum_W": total,
        "flux_inside_die_W": float(inside),
        "frac_into_spreader_under_die": float(inside / total) if total else 0.0,
        "peak_flux_W_per_mm2": float(flux.max() / cell / 1e6),
    }


def field_err(approx: np.ndarray, truth: np.ndarray) -> dict:
    d = approx - truth
    return {
        "dTmax_K": float(approx.max() - truth.max()),
        "max_abs_K": float(np.abs(d).max()),
        "rmse_K": float(np.sqrt(np.mean(d**2))),
    }


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    k_si = np.full((hs.N, hs.N), hs.K_SI)
    t0 = time.perf_counter()
    mat, g_amb = hs.build_operator(k_si)
    factor = splu(mat.tocsc())
    t_factor = time.perf_counter() - t0

    ref_p = power_of(SIDE0, "fixed_W")
    t0 = time.perf_counter()
    ref = hs.solve_state(k_si, ref_p, factor=factor, g_amb=g_amb)
    t_apply = time.perf_counter() - t0

    rows = []
    saved = {}

    def apply(power: np.ndarray) -> np.ndarray:
        return hs.solve_state(k_si, power, factor=factor, g_amb=g_amb)

    def correct(power: np.ndarray) -> np.ndarray:
        # 右端项只有功耗差，没有环境项。和基准场相加就是新场。
        # 两场的环境项相同，相减后外伸节点和热沉格的环境项都是 0。
        rhs = np.zeros(hs.N_UNK)
        rhs[: hs.N * hs.N] = (power - ref_p).ravel()
        delta = factor.solve(rhs)[: hs.N_GRID].reshape(hs.NL, hs.N, hs.N)
        return ref + delta

    for law in ("fixed_W", "prop_area"):
        for side in SIDES:
            power = power_of(side, law)
            truth = apply(power)
            info = pack(truth, power, k_si, side)
            info.update({
                "law": law,
                "side_cells": side,
                "side_mm": side * hs.CW * 1e3,
                "reuse": field_err(ref[0], truth[0]),
                "scale_power": field_err(
                    hs.T_AMB + (ref[0] - hs.T_AMB) * (power.sum() / P0),
                    truth[0],
                ),
                "delta_p": field_err(correct(power)[0], truth[0]),
            })
            rows.append(info)
            print(
                f"{law:10s} {info['side_mm']:5.2f} mm  P {info['power_W']:7.1f} W  "
                f"Tmax {info['Tmax_C']:.2f} °C  "
                f"直达散热片 {info['frac_into_spreader_under_die']:.3f}  "
                f"原场不动 {info['reuse']['max_abs_K']:.2f} K  "
                f"只缩放 {info['scale_power']['max_abs_K']:.2f} K  "
                f"功耗差修正 {info['delta_p']['max_abs_K']:.2e} K",
                flush=True,
            )
            if side in (8, 16, 32):
                saved[f"{law}_{side}"] = truth[0]

    # 外圈变成硅、但功耗图不动。这是上一份实验的面积版。
    ring = []
    for side in (16, 24, 32):
        power = power_of(side, "same_map")
        a = origin(side)
        k_ring = np.full((hs.N, hs.N), hs.K_FILL)
        k_ring[a : a + side, a : a + side] = hs.K_SI
        hetero = hs.solve_state(k_ring, power)[0]
        uniform = apply(power)[0]
        err = field_err(uniform, hetero)
        # uniform 是“外圈仍按硅算”的场。和真的环氧外圈比，差多少。
        # 功耗图和基准相同，所以 uniform 应该就等于基准场。
        same = field_err(uniform, ref[0])
        ring.append({
            "side_mm": side * hs.CW * 1e3,
            "silicon_vs_epoxy_ring_max_abs_K": err["max_abs_K"],
            "silicon_vs_epoxy_ring_dTmax_K": err["dTmax_K"],
            "power_map_unchanged_vs_reference_max_abs_K": same["max_abs_K"],
        })
        print(
            f"外圈 {ring[-1]['side_mm']:.2f} mm  "
            f"硅/环氧最大差 {err['max_abs_K']:.4f} K  "
            f"功耗图没变时相对基准 {same['max_abs_K']:.2e} K",
            flush=True,
        )

    report = {
        "package_side_mm": hs.SIDE * 1e3,
        "grid": hs.N,
        "reference": {"side_mm": SIDE0 * hs.CW * 1e3, "power_W": P0,
                      "Tmax_C": float(ref[0].max() - 273.15)},
        "laws": {
            "fixed_W": "总功率固定 80 W，密度随面积下降",
            "prop_area": "功耗密度不变，总功率与面积成正比",
        },
        "timing_s": {"factor": t_factor, "one_apply": t_apply},
        "cases": rows,
        "ring_material": ring,
        "note": (
            "绝对温度不含 HotSpot 的 spreader 外伸环。"
            "delta_p 修正用的是同一个已分解的整层硅算子。"
        ),
    }
    (OUT / "report.json").write_text(
        json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    for name, field in saved.items():
        np.save(OUT / f"{name}.npy", field)
    np.save(OUT / "reference.npy", ref[0])
    _plot(rows)
    print(f"写出 {OUT / 'report.json'}", flush=True)


def _plot(rows: list[dict]) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(9.2, 3.6))
    for law, color in (("fixed_W", "C0"), ("prop_area", "C3")):
        pts = [r for r in rows if r["law"] == law]
        x = [r["side_mm"] for r in pts]
        axes[0].plot(x, [r["Tmax_C"] for r in pts], "o-", color=color, label=law)
        axes[1].plot(x, [r["reuse"]["max_abs_K"] for r in pts], "o-", color=color, label=f"{law} reuse")
        axes[1].plot(
            x, [r["scale_power"]["max_abs_K"] for r in pts], "s--", color=color, label=f"{law} scale"
        )
    axes[0].set_xlabel("die side (mm)")
    axes[0].set_ylabel("Tmax (°C)")
    axes[0].legend()
    axes[1].set_xlabel("die side (mm)")
    axes[1].set_ylabel("max |ΔT| vs truth (K)")
    axes[1].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT / "area_tmax_and_error.png", dpi=120)
    plt.close(fig)


if __name__ == "__main__":
    main()
