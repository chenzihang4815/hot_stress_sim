"""B*-tree 随机初始布局（Chang et al., DAC 2000）。

每个模块是树节点：左孩子 = 右边邻接，右孩子 = 上方邻接。
先随机生成树和长宽比，再确定性装箱，保证无重叠。
"""

from __future__ import annotations

import math
import random
from dataclasses import dataclass

from blocks import Block
from geometry import Rect, normalize_origin, packing_metrics, size_from_area


@dataclass
class Node:
    idx: int
    left: int | None = None
    right: int | None = None


def _random_tree(n: int, rng: random.Random) -> list[Node]:
    nodes = [Node(i) for i in range(n)]
    if n == 1:
        return nodes
    order = list(range(n))
    rng.shuffle(order)
    root = order[0]
    parent_of = {root: None}
    frontier = [root]
    for idx in order[1:]:
        while True:
            p = rng.choice(frontier)
            slots = []
            if nodes[p].left is None:
                slots.append("L")
            if nodes[p].right is None:
                slots.append("R")
            if not slots:
                frontier.remove(p)
                continue
            side = rng.choice(slots)
            if side == "L":
                nodes[p].left = idx
            else:
                nodes[p].right = idx
            parent_of[idx] = p
            frontier.append(idx)
            if nodes[p].left is not None and nodes[p].right is not None:
                frontier.remove(p)
            break
    # 把根放到 0 号逻辑：返回 nodes 以及 root
    nodes[0].__dict__["root"] = root  # 轻量标记
    for nobj in nodes:
        nobj.__dict__["root"] = root
    return nodes


def _pack(blocks: list[Block], nodes: list[Node], aspects: list[float]) -> list[Rect]:
    n = len(blocks)
    sizes = [size_from_area(b.area_mm2, a) for b, a in zip(blocks, aspects)]
    x = [0.0] * n
    y = [0.0] * n
    contour: list[tuple[float, float, float]] = []  # (x1, x2, height)

    def contour_height(x1: float, x2: float) -> float:
        h = 0.0
        for a, b, hh in contour:
            if b <= x1 or a >= x2:
                continue
            h = max(h, hh)
        return h

    def raise_contour(x1: float, x2: float, new_h: float) -> None:
        nonlocal contour
        updated: list[tuple[float, float, float]] = []
        for a, b, hh in contour:
            if b <= x1 or a >= x2:
                updated.append((a, b, hh))
                continue
            if a < x1:
                updated.append((a, x1, hh))
            if b > x2:
                updated.append((x2, b, hh))
        updated.append((x1, x2, new_h))
        updated.sort()
        merged: list[tuple[float, float, float]] = []
        for seg in updated:
            if merged and abs(merged[-1][1] - seg[0]) < 1e-15 and abs(merged[-1][2] - seg[2]) < 1e-15:
                merged[-1] = (merged[-1][0], seg[1], merged[-1][2])
            else:
                merged.append(seg)
        contour = merged

    def dfs(idx: int | None, parent: int | None, as_left: bool) -> None:
        if idx is None:
            return
        w, h = sizes[idx]
        if parent is None:
            x[idx] = 0.0
        elif as_left:
            x[idx] = x[parent] + sizes[parent][0]
        else:
            x[idx] = x[parent]
        y[idx] = contour_height(x[idx], x[idx] + w)
        raise_contour(x[idx], x[idx] + w, y[idx] + h)
        dfs(nodes[idx].left, idx, True)
        dfs(nodes[idx].right, idx, False)

    root = nodes[0].__dict__["root"]
    dfs(root, None, True)
    rects = [
        Rect(b.name, x[i], y[i], sizes[i][0], sizes[i][1], b.power_W, b.area_mm2, b.kind)
        for i, b in enumerate(blocks)
    ]
    return normalize_origin(rects)


def random_bstree_layout(
    blocks: list[Block],
    seed: int = 0,
    n_candidates: int = 48,
) -> tuple[list[Rect], dict]:
    rng = random.Random(seed)
    best = None
    best_cost = math.inf
    best_meta = None
    for k in range(n_candidates):
        nodes = _random_tree(len(blocks), rng)
        aspects = [math.exp(rng.uniform(-0.7, 0.7)) for _ in blocks]
        rects = _pack(blocks, nodes, aspects)
        m = packing_metrics(rects)
        # 优先无重叠（B* 本身无重叠），再压轮廓面积和接近方形
        aspect_pen = abs(math.log((m["width_m"] + 1e-18) / (m["height_m"] + 1e-18)))
        cost = m["outline_mm2"] * (1.0 + 0.15 * aspect_pen) + 1e9 * m["overlap_m2"]
        if cost < best_cost:
            best_cost = cost
            best = rects
            best_meta = {"candidate": k, "cost": cost, **m}
    assert best is not None and best_meta is not None
    best_meta["algorithm"] = "bstree"
    best_meta["seed"] = seed
    best_meta["n_candidates"] = n_candidates
    return best, best_meta
