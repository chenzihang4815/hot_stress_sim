"""L-BFGS 连续松弛布局（用户写的 L-FBGS 按 L-BFGS 实现）。

决策变量：每个模块 (x, y, log(aspect))。面积固定。
目标 = 轮廓面积 + 方形惩罚 + 重叠惩罚。
不用 scipy：numpy 实现 limited-memory BFGS 两循环。
"""

from __future__ import annotations

import math

import numpy as np

from blocks import Block
from geometry import Rect, normalize_origin, packing_metrics, size_from_area


def _decode(blocks: list[Block], z: np.ndarray) -> list[Rect]:
    n = len(blocks)
    xywh = []
    for i, b in enumerate(blocks):
        x, y, la = z[3 * i : 3 * i + 3]
        w, h = size_from_area(b.area_mm2, float(math.exp(la)))
        xywh.append((float(x), float(y), w, h))
    rects = [
        Rect(b.name, x, y, w, h, b.power_W, b.area_mm2, b.kind)
        for b, (x, y, w, h) in zip(blocks, xywh)
    ]
    return rects


def _objective(blocks: list[Block], z: np.ndarray) -> float:
    rects = _decode(blocks, z)
    m = packing_metrics(rects)
    # 把模块尽量推到原点附近，避免整体漂移
    origin = sum(max(0.0, -r.x) + max(0.0, -r.y) for r in rects)
    aspect_pen = abs(math.log((m["width_m"] + 1e-18) / (m["height_m"] + 1e-18)))
    return (
        m["outline_mm2"]
        + 0.2 * m["outline_mm2"] * aspect_pen
        + 5e10 * m["overlap_m2"]
        + 1e3 * origin
    )


def _grad(blocks: list[Block], z: np.ndarray, eps: float = 1e-7) -> np.ndarray:
    g = np.zeros_like(z)
    f0 = _objective(blocks, z)
    for i in range(z.size):
        z2 = z.copy()
        step = eps if abs(z[i]) < 1 else eps * abs(z[i])
        z2[i] += step
        g[i] = (_objective(blocks, z2) - f0) / step
    return g


def _lbfgs_direction(grad: np.ndarray, mem_s: list[np.ndarray], mem_y: list[np.ndarray]) -> np.ndarray:
    q = grad.copy()
    rho = []
    alpha = []
    for s, y in zip(reversed(mem_s), reversed(mem_y)):
        r = 1.0 / (float(y @ s) + 1e-18)
        rho.append(r)
        a = r * float(s @ q)
        alpha.append(a)
        q = q - a * y
    if mem_y:
        y = mem_y[-1]
        s = mem_s[-1]
        gamma = float(s @ y) / (float(y @ y) + 1e-18)
        rvec = gamma * q
    else:
        rvec = q
    for s, y, r, a in zip(mem_s, mem_y, reversed(rho), reversed(alpha)):
        b = r * float(y @ rvec)
        rvec = rvec + s * (a - b)
    return -rvec


def random_lbfgs_layout(
    blocks: list[Block],
    seed: int = 2,
    max_iter: int = 40,
    mem: int = 8,
) -> tuple[list[Rect], dict]:
    rng = np.random.default_rng(seed)
    n = len(blocks)
    # 先用随机方形网格当初值，再连续松弛
    order = rng.permutation(n)
    z = np.zeros(3 * n)
    cursor_x = 0.0
    row_h = 0.0
    row_y = 0.0
    row_w = 0.0
    cols = max(1, int(math.ceil(math.sqrt(n))))
    for k, idx in enumerate(order):
        if k and k % cols == 0:
            cursor_x = 0.0
            row_y += row_h + 5e-5
            row_h = 0.0
        la = float(rng.uniform(-0.4, 0.4))
        w, h = size_from_area(blocks[idx].area_mm2, math.exp(la))
        z[3 * idx] = cursor_x
        z[3 * idx + 1] = row_y
        z[3 * idx + 2] = la
        cursor_x += w + 5e-5
        row_h = max(row_h, h)
        row_w = max(row_w, cursor_x)

    mem_s: list[np.ndarray] = []
    mem_y: list[np.ndarray] = []
    f = _objective(blocks, z)
    g = _grad(blocks, z)
    history = [f]

    for _ in range(max_iter):
        d = _lbfgs_direction(g, mem_s, mem_y)
        # 回溯线搜索
        step = 1.0
        accepted = False
        for _ls in range(12):
            z_new = z + step * d
            f_new = _objective(blocks, z_new)
            if f_new < f - 1e-6 * step * float(g @ d) or f_new < f:
                g_new = _grad(blocks, z_new)
                s = z_new - z
                y = g_new - g
                if float(y @ s) > 1e-16:
                    mem_s.append(s)
                    mem_y.append(y)
                    if len(mem_s) > mem:
                        mem_s.pop(0)
                        mem_y.pop(0)
                z, f, g = z_new, f_new, g_new
                accepted = True
                break
            step *= 0.5
        history.append(f)
        if not accepted:
            break
        if abs(history[-1] - history[-2]) < 1e-6 * (1 + abs(history[-1])):
            break

    rects = normalize_origin(_decode(blocks, z))
    # 若仍重叠，把模块按 x 排序后水平摊开（保证 HotSpot 能跑）
    m = packing_metrics(rects)
    if m["overlap_m2"] > 1e-18:
        rects = _sweep_apart(rects)
        m = packing_metrics(rects)
    meta = m
    meta.update(
        {
            "algorithm": "lbfgs",
            "seed": seed,
            "max_iter": max_iter,
            "cost": f,
            "history_tail": history[-5:],
        }
    )
    return rects, meta


def _sweep_apart(rects: list[Rect]) -> list[Rect]:
    ordered = sorted(rects, key=lambda r: (r.x, r.y))
    x = 0.0
    out = []
    for r in ordered:
        out.append(Rect(r.name, x, 0.0, r.w, r.h, r.power_W, r.area_mm2, r.kind))
        x += r.w + 2e-5
    return out
