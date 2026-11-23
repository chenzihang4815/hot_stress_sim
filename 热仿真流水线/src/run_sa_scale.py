#!/usr/bin/env python3
"""Shrink one die until a coarse search stops matching a fine one.

In each of the six layouts the smallest die that actually dissipates power is
the variable. It is scaled about its own centre through 1/2, 1/4 and 1/8 while
every other die keeps its size and place, and its power follows its area so the
density stays put. The same 120-step big-jump search then runs twice, priced
entirely on the 64 grid and entirely on the 32 grid, and both winners are
re-scored on the 64 grid. The scale where the two finals first differ is the
critical size for that layout.

One HotSpot at a time. Resumable from results.jsonl.
"""

from __future__ import annotations

import importlib.util
import json
import sys
import time
from pathlib import Path

ROOT = Path("/home/czh/任务仿真产生功耗文件工具")
SRC = ROOT / "热仿真流水线" / "src"
sys.path.insert(0, str(SRC))

import sa_thermal_feedback as sa  # noqa: E402

spec = importlib.util.spec_from_file_location("m", SRC / "run_sa_multicase.py")
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)

OUT = ROOT / "热仿真流水线" / "runs" / "superposition" / "scale"
SCALES = [0.5, 0.25, 0.125]


def variable_of(chips):
    """The powered die with the shortest side. Ties go to the hotter one."""
    live = [c for c in chips if c.power_W > 0.1]
    return min(live, key=lambda c: (min(c.w, c.h), -c.power_W)).name


def shrink_one(chips, name, scale):
    """Scale one die about its centre. Everyone else is untouched."""
    out = []
    for c in chips:
        if c.name != name:
            out.append(sa.Chip(c.name, c.kind, c.x, c.y, c.w, c.h, c.power_W, c.area_mm2))
            continue
        w, h = c.w * scale, c.h * scale
        out.append(sa.Chip(c.name, c.kind, c.x + (c.w - w) / 2, c.y + (c.h - h) / 2,
                           w, h, c.power_W * scale * scale, c.area_mm2 * scale * scale))
    return out


def grids_main() -> None:
    """Run the search priced entirely on the 16, 32 and 64 grids.

    Winners are re-scored on the 64 grid, so a match means the coarse grid
    ranked every move the same way, not merely that its temperatures were close.
    """
    OUT.mkdir(parents=True, exist_ok=True)
    m.f.OUT = OUT / "_work"
    log = OUT / "grids.jsonl"
    done = set()
    if log.exists():
        for line in log.read_text().splitlines():
            r = json.loads(line)
            done.add((r["label"], r["rows"]))
    for label, path in m.CASES:
        chips = m.load(path)
        for rows in (64, 32, 16):
            if (label, rows) in done:
                continue
            saved = m.field_temp
            def priced(cs, cache, grid, _rows=rows):
                return saved(cs, cache, _rows)
            m.field_temp = priced
            t0 = time.time()
            try:
                _, final, secs = m.search(chips, m.JUMP, m.JUMP, 1.0, False)
            finally:
                m.field_temp = saved
            rec = {"label": label, "chips": len(chips), "rows": rows,
                   "final_64_C": final, "secs": secs, "wall_s": round(time.time() - t0, 1)}
            with log.open("a") as fh:
                fh.write(json.dumps(rec, ensure_ascii=False) + "\n")
            print(json.dumps(rec, ensure_ascii=False), flush=True)


def conv_main() -> None:
    """Per-step best temperature for the four move mixes, priced on the 16 grid."""
    OUT.mkdir(parents=True, exist_ok=True)
    m.f.OUT = OUT / "_work"
    log = OUT / "conv.jsonl"
    done = set()
    if log.exists():
        for line in log.read_text().splitlines():
            r = json.loads(line)
            done.add((r["label"], r["mix"]))
    saved = m.field_temp
    m.field_temp = lambda cs, cache, grid: saved(cs, cache, 16)
    extra = [
        ("xeon1_9芯", "runs/full_chain/pkg_1xeon_1t2_2gpu/Lstar.json"),
        ("xeon2g3_10芯", "runs/full_chain/pkg_2xeon_3gpu/Lstar.json"),
        ("cmp_3芯", "runs/sa_compare_schemes/pkg_1cpu_1gpu_1cim/layout.json"),
        ("cmp_5芯", "runs/sa_compare_schemes/pkg_2cpu_1gpu_2cim/layout.json"),
        ("cmp_6芯", "runs/sa_compare_schemes/pkg_2cpu_2gpu_2cim/layout.json"),
        ("cmp_7芯", "runs/sa_compare_schemes/pkg_3cpu_1gpu_3cim/layout.json"),
        ("cmp_10芯", "runs/sa_compare_schemes/pkg_4cpu_2gpu_4cim/layout.json"),
        ("cmp_11芯", "runs/sa_compare_schemes/pkg_4cpu_3gpu_4cim/layout.json"),
        ("cmp_12芯", "runs/sa_compare_schemes/pkg_4cpu_3gpu_5cim/layout.json"),
        ("cmp_13芯", "runs/sa_compare_schemes/pkg_5cpu_3gpu_5cim/layout.json"),
        ("cmp_14芯", "runs/sa_compare_schemes/pkg_5cpu_4gpu_5cim/layout.json"),
        ("cmp_15芯", "runs/sa_compare_schemes/pkg_5cpu_5gpu_5cim/layout.json"),
        ("penryn_13芯", "runs/sa_compare_schemes/pkg_5penryn_3gpu_5cim/layout.json"),
    ]
    try:
        for label, path in list(m.CASES) + extra:
            chips = m.load(path)
            for name, mix in m.MIXES:
                if (label, name) in done:
                    continue
                t0 = time.time()
                best_at, final, secs = m.search(chips, mix, mix, 1.0, True)
                rec = {"label": label, "chips": len(chips), "mix": name, "rows": 16,
                       "final_64_C": final, "secs": secs, "wall_s": round(time.time() - t0, 1),
                       "best_C": [round(v, 3) for v in best_at]}
                with log.open("a") as fh:
                    fh.write(json.dumps(rec, ensure_ascii=False) + "\n")
                print(json.dumps({k: v for k, v in rec.items() if k != "best_C"},
                                 ensure_ascii=False), flush=True)
    finally:
        m.field_temp = saved


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    m.f.OUT = OUT / "_work"
    log = OUT / "results.jsonl"
    done = set()
    if log.exists():
        for line in log.read_text().splitlines():
            r = json.loads(line)
            done.add((r["label"], r["scale"], r["coarse_frac"]))
    for label, path in m.CASES:
        base = m.load(path)
        name = variable_of(base)
        for scale in SCALES:
            chips = shrink_one(base, name, scale)
            if not sa.legal(chips, m.INTP, m.INTP):
                chips = sa.spread_to_gap(chips, m.INTP, m.INTP)
            if not sa.legal(chips, m.INTP, m.INTP):
                print(f"{label} 的 {name} 缩到 {scale} 后摊不开，跳过", flush=True)
                continue
            var = next(c for c in chips if c.name == name)
            side = min(var.w, var.h) * 1e3
            for frac in (0.0, 1.0):
                if (label, scale, frac) in done:
                    continue
                t0 = time.time()
                _, final, secs = m.search(chips, m.JUMP, m.JUMP, frac, False)
                rec = {"label": label, "chips": len(chips), "variable": name,
                       "scale": scale, "side_mm": round(side, 2),
                       "power_W": round(var.power_W, 3), "coarse_frac": frac,
                       "final_64_C": final, "secs": secs,
                       "wall_s": round(time.time() - t0, 1)}
                with log.open("a") as fh:
                    fh.write(json.dumps(rec, ensure_ascii=False) + "\n")
                print(json.dumps(rec, ensure_ascii=False), flush=True)


if __name__ == "__main__":
    import sys
    if "--grids" in sys.argv:
        grids_main()
    elif "--conv" in sys.argv:
        conv_main()
    else:
        main()
