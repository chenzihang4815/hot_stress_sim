"""无重叠矩形几何：面积固定，长宽比可变。"""

from __future__ import annotations

import math
from dataclasses import dataclass

from blocks import Block


@dataclass
class Rect:
    name: str
    x: float  # m, 左下角
    y: float
    w: float
    h: float
    power_W: float
    area_mm2: float
    kind: str

    @property
    def x2(self) -> float:
        return self.x + self.w

    @property
    def y2(self) -> float:
        return self.y + self.h


def size_from_area(area_mm2: float, aspect: float) -> tuple[float, float]:
    """aspect = width/height。返回米。"""
    area_m2 = area_mm2 * 1e-6
    aspect = min(max(aspect, 0.25), 4.0)
    w = math.sqrt(area_m2 * aspect)
    h = area_m2 / w
    return w, h


def overlap_area(a: Rect, b: Rect) -> float:
    dx = min(a.x2, b.x2) - max(a.x, b.x)
    dy = min(a.y2, b.y2) - max(a.y, b.y)
    if dx <= 0 or dy <= 0:
        return 0.0
    return dx * dy


def packing_metrics(rects: list[Rect]) -> dict[str, float]:
    if not rects:
        return {
            "width_m": 0.0,
            "height_m": 0.0,
            "outline_mm2": 0.0,
            "block_area_mm2": 0.0,
            "deadspace_ratio": 0.0,
            "overlap_m2": 0.0,
            "overlap_count": 0.0,
        }
    min_x = min(r.x for r in rects)
    min_y = min(r.y for r in rects)
    width = max(r.x2 for r in rects) - min_x
    height = max(r.y2 for r in rects) - min_y
    outline = width * height * 1e6
    block = sum(r.area_mm2 for r in rects)
    overlap = 0.0
    n_ov = 0
    for i, a in enumerate(rects):
        for b in rects[i + 1 :]:
            oa = overlap_area(a, b)
            if oa > 1e-18:
                overlap += oa
                n_ov += 1
    return {
        "width_m": width,
        "height_m": height,
        "outline_mm2": outline,
        "block_area_mm2": block,
        "deadspace_ratio": max(0.0, 1.0 - block / outline) if outline > 0 else 0.0,
        "overlap_m2": overlap,
        "overlap_count": float(n_ov),
    }


def normalize_origin(rects: list[Rect]) -> list[Rect]:
    min_x = min(r.x for r in rects)
    min_y = min(r.y for r in rects)
    out = []
    for r in rects:
        out.append(
            Rect(r.name, r.x - min_x, r.y - min_y, r.w, r.h, r.power_W, r.area_mm2, r.kind)
        )
    return out


def rects_from_xywh(blocks: list[Block], xywh: list[tuple[float, float, float, float]]) -> list[Rect]:
    rects = [
        Rect(b.name, x, y, w, h, b.power_W, b.area_mm2, b.kind)
        for b, (x, y, w, h) in zip(blocks, xywh)
    ]
    return normalize_origin(rects)
