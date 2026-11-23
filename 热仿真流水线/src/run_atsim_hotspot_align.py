#!/usr/bin/env python3
"""只算热，把 ATSim 的封装改成和 HotSpot 同一套之后，温度还差多少。

不跑弹性。芯片是 sa_many_macro 的起点。HotSpot 的格子温度已经在
hotspot_init/grid.steady 里，这里不再重跑 HotSpot。
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path("/home/czh/任务仿真产生功耗文件工具")
ATSIM = ROOT / "热应力python工程_最终"
sys.path.insert(0, str(ATSIM))

from atsim.fvm3d import solve_fvm3d  # noqa: E402
from atsim.geometry import load_case  # noqa: E402

PIPE = ROOT / "热仿真流水线"
OUT = PIPE / "runs" / "atsim_hotspot_align"
LAYOUT = PIPE / "runs" / "sa_many_macro" / "layouts" / "init.json"
GRID = PIPE / "runs" / "sa_many_macro" / "hotspot_init" / "grid.steady"
CW_MM = 0.072 / 64 * 1e3


def chips() -> list[dict]:
    rows = json.loads(LAYOUT.read_text(encoding="utf-8"))
    out = []
    for d in rows:
        out.append({
            "name": d["name"],
            "x_mm": d["j0"] * CW_MM,
            "y_mm": d["ib0"] * CW_MM,
            "w_mm": d["w"] * CW_MM,
            "h_mm": d["h"] * CW_MM,
            "t_um": 150.0,
            "power_W": d["power_W"],
        })
    return out


def base_pkg() -> dict:
    return {
        "tim_um": 20, "t_spreader_mm": 1.0, "t_sink_mm": 6.9,
        "t_substrate_mm": 1.0, "t_interposer_um": 100.0,
        "ambient_C": 45,
        "pitch_die_mm": 1.125, "pitch_coarse_mm": 3.0,
        "nz_die": 3, "nz_spreader": 3, "nz_sink": 4, "nz_tim": 1,
        "nz_substrate": 3, "nz_interposer": 2,
        "fill_die": True, "fill_material": "silicon",
    }


def packages() -> dict[str, dict]:
    old = base_pkg()
    old.update(r_convec_K_per_W=0.1, fill_xy="chips", tim_cover="chips")
    aligned = base_pkg()
    aligned.update(
        s_spreader_mm=96.0, s_sink_mm=96.0, frame_mm=72.0,
        tim_cover="frame", r_convec_K_per_W=0.0390625, fill_xy="chips",
    )
    no_under = dict(aligned)
    no_under.update(t_substrate_mm=0.0, t_interposer_um=0.0)
    return {"legacy_margin": old, "aligned_with_substrate": aligned, "aligned_no_underlayer": no_under}


def hotspot_layer0() -> dict:
    text = GRID.read_text().split()
    vals = []
    i = 0
    while i < len(text):
        if text[i] == "Layer" and text[i + 1].startswith("0"):
            i += 2
            while i < len(text) and text[i] != "Layer":
                vals.append(float(text[i + 1]))
                i += 2
            break
        i += 1
    a = np.array(vals) - 273.15
    return {"n": int(a.size), "min_C": float(a.min()), "mean_C": float(a.mean()), "max_C": float(a.max())}


def die_stats(th) -> dict:
    grid = th["grid"]
    case = th["case"]
    layer_id = grid["layer_id"]
    die_ids = [i for i, L in enumerate(case.stack) if L.name == "die_fill" or L.name.startswith("die_")]
    mask = np.isin(layer_id, die_ids) & grid["solid"]
    T = th["T"]
    vals = T[mask] - 273.15
    # 每一列先把芯片层厚度平均，再按格子面积平均。HotSpot 芯片层是一整层一个节点。
    col = np.zeros(mask.shape[:2], dtype=bool)
    acc = np.zeros(mask.shape[:2], dtype=float)
    cnt = np.zeros(mask.shape[:2], dtype=float)
    for k in range(mask.shape[2]):
        m = mask[:, :, k]
        col |= m
        acc[m] += T[:, :, k][m]
        cnt[m] += 1.0
    area = grid["dx"][:, None] * grid["dy"][None, :]
    mean_col = np.divide(acc, cnt, out=np.zeros_like(acc), where=cnt > 0)
    w = area[col]
    mean = float(np.sum(mean_col[col] * w) / np.sum(w) - 273.15)
    meta = case.raw.get("_layout_meta", {})
    return {
        "cell_max_C": float(vals.max()),
        "cell_min_C": float(vals.min()),
        "column_mean_C": mean,
        "n_cells": int(vals.size),
        "grid": [int(grid["nx"]), int(grid["ny"]), int(grid["nz"])],
        "layout": {
            "spreader_mm": meta.get("spreader_mm"),
            "sink_mm": meta.get("sink_mm"),
            "frame_mm": meta.get("frame_mm"),
            "h_W_m2K": meta.get("h_W_m2K"),
            "t_substrate_mm": meta.get("t_substrate_mm"),
            "t_interposer_um": meta.get("t_interposer_um"),
        },
        "elapsed_s": th["elapsed_s"],
    }


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    hs = hotspot_layer0()
    rows = {}
    for name, pkg in packages().items():
        path = OUT / f"{name}.json"
        payload = {"name": name, "chips": chips(), "package": pkg}
        path.write_text(json.dumps(payload) + "\n", encoding="utf-8")
        case = load_case(str(path))
        th = solve_fvm3d(case, verbose=True)
        st = die_stats(th)
        st["minus_hotspot_max_K"] = st["cell_max_C"] - hs["max_C"]
        st["minus_hotspot_mean_K"] = st["column_mean_C"] - hs["mean_C"]
        rows[name] = st
        print(
            f"{name}: die max {st['cell_max_C']:.2f}  mean {st['column_mean_C']:.2f}  "
            f"vs HS max {st['minus_hotspot_max_K']:+.2f}  mean {st['minus_hotspot_mean_K']:+.2f}",
            flush=True,
        )
    report = {
        "hotspot_layer0_C": hs,
        "note": (
            "HotSpot 芯片层 72 mm，spreader 和热沉 96 mm，r_convec=0.0390625。"
            "legacy 是旧代理：r_convec=0.1，热沉跟着外框加边。"
            "aligned 用 96 mm 和 0.0390625，芯片框 72 mm，基板和中介层还在。"
            "aligned_no_underlayer 把基板和中介层去掉，才是 HotSpot 的导热叠层。"
            "网格 pitch 1.125 mm，比代理应力用的 2 mm 细。没有算应力。"
            "Python 场和 ATSim 的温度不要混。"
        ),
        "cases": rows,
    }
    (OUT / "report.json").write_text(
        json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"写出 {OUT / 'report.json'}", flush=True)


if __name__ == "__main__":
    main()
