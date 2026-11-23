"""MCTS 随机初始布局。

把「先放哪个模块、放在已有轮廓的哪条边」当成树搜索。
随机 rollout 补全其余模块，用轮廓面积 + 重叠惩罚当回报。
这是 floorplan 用的 Monte Carlo Tree Search，不是博弈棋类。
"""

from __future__ import annotations

import math
import random
from dataclasses import dataclass, field

from blocks import Block
from geometry import Rect, normalize_origin, packing_metrics, size_from_area


@dataclass
class Place:
    idx: int
    x: float
    y: float
    w: float
    h: float


def _legal_anchors(placed: list[Place]) -> list[tuple[float, float]]:
    if not placed:
        return [(0.0, 0.0)]
    pts = {(0.0, 0.0)}
    for p in placed:
        pts.add((p.x + p.w, p.y))
        pts.add((p.x, p.y + p.h))
    return sorted(pts)


def _overlaps(x: float, y: float, w: float, h: float, placed: list[Place]) -> bool:
    x2, y2 = x + w, y + h
    for p in placed:
        if x2 <= p.x or p.x + p.w <= x or y2 <= p.y or p.y + p.h <= y:
            continue
        return True
    return False


def _place_one(
    idx: int,
    area_mm2: float,
    placed: list[Place],
    rng: random.Random,
) -> Place:
    aspect = math.exp(rng.uniform(-0.7, 0.7))
    w, h = size_from_area(area_mm2, aspect)
    anchors = _legal_anchors(placed)
    rng.shuffle(anchors)
    for x, y in anchors:
        if not _overlaps(x, y, w, h, placed):
            return Place(idx, x, y, w, h)
    # 实在贴不上就放到当前轮廓右侧
    max_x = max((p.x + p.w for p in placed), default=0.0)
    return Place(idx, max_x, 0.0, w, h)


def _complete(order_prefix: list[int], remaining: list[int], blocks: list[Block], rng: random.Random) -> list[Place]:
    placed: list[Place] = []
    for idx in order_prefix + remaining:
        placed.append(_place_one(idx, blocks[idx].area_mm2, placed, rng))
    return placed


def _rects(blocks: list[Block], placed: list[Place]) -> list[Rect]:
    by_idx = {p.idx: p for p in placed}
    rects = [
        Rect(b.name, by_idx[i].x, by_idx[i].y, by_idx[i].w, by_idx[i].h, b.power_W, b.area_mm2, b.kind)
        for i, b in enumerate(blocks)
    ]
    return normalize_origin(rects)


def _cost(blocks: list[Block], placed: list[Place]) -> float:
    m = packing_metrics(_rects(blocks, placed))
    aspect_pen = abs(math.log((m["width_m"] + 1e-18) / (m["height_m"] + 1e-18)))
    return m["outline_mm2"] * (1.0 + 0.15 * aspect_pen) + 1e9 * m["overlap_m2"]


@dataclass
class Node:
    prefix: tuple[int, ...]
    remaining: tuple[int, ...]
    parent: "Node | None" = None
    children: dict[int, "Node"] = field(default_factory=dict)
    untried: list[int] = field(default_factory=list)
    visits: int = 0
    value: float = 0.0  # 越小越好，存 -cost 的均值

    def uct(self, parent_visits: int, c: float = 1.25) -> float:
        if self.visits == 0:
            return 1e9
        return self.value / self.visits + c * math.sqrt(math.log(parent_visits + 1) / self.visits)


def random_mcts_layout(
    blocks: list[Block],
    seed: int = 1,
    n_sims: int = 80,
) -> tuple[list[Rect], dict]:
    rng = random.Random(seed)
    n = len(blocks)
    root = Node(prefix=(), remaining=tuple(range(n)), untried=list(range(n)))
    best_placed = None
    best_cost = math.inf

    def rollout(node: Node) -> float:
        remaining = list(node.remaining)
        rng.shuffle(remaining)
        placed = _complete(list(node.prefix), remaining, blocks, rng)
        cost = _cost(blocks, placed)
        nonlocal best_placed, best_cost
        if cost < best_cost:
            best_cost = cost
            best_placed = placed
        return -cost

    for _ in range(n_sims):
        node = root
        # select
        while not node.untried and node.children:
            node = max(node.children.values(), key=lambda ch: ch.uct(node.visits))
        # expand
        if node.untried:
            choice = node.untried.pop(rng.randrange(len(node.untried)))
            new_prefix = node.prefix + (choice,)
            new_rem = tuple(x for x in node.remaining if x != choice)
            child = Node(prefix=new_prefix, remaining=new_rem, parent=node, untried=list(new_rem))
            node.children[choice] = child
            node = child
        reward = rollout(node)
        while node is not None:
            node.visits += 1
            node.value += reward
            node = node.parent

    assert best_placed is not None
    rects = _rects(blocks, best_placed)
    meta = packing_metrics(rects)
    meta.update(
        {
            "algorithm": "mcts",
            "seed": seed,
            "n_sims": n_sims,
            "cost": best_cost,
        }
    )
    return rects, meta
