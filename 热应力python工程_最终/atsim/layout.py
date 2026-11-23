# -*- coding: utf-8 -*-
"""Expand a chips[] layout into a full package Case dict.

User-facing input (mm / um / W) — everything else is filled in:

    {
      "name": "my_3chip",
      "chips": [
        {"name": "CPU", "x_mm": 5, "y_mm": 5, "w_mm": 8, "h_mm": 8,
         "t_um": 150, "power_W": 20}
      ]
    }

Optional package{} overrides TIM / spreader / sink / substrate / interposer.
Existing case1 JSON (no chips[]) is left unchanged. Case 2 writes stack[] itself.
"""

from __future__ import print_function

import copy
import math


DEFAULT_MATERIALS = {
    "silicon": {
        "kappa0_W_mK": 130.0,
        "T_ref_K": 300.0,
        "alpha": 1.5,
        "rho_kg_m3": 2330.0,
        "cp_J_kgK": 700.0,
        "E_Pa": 1.7e11,
        "nu": 0.28,
        "alpha_cte_1_K": 2.6e-6,
    },
    "tim": {
        "kappa0_W_mK": 4.0,
        "T_ref_K": 300.0,
        "alpha": 0.0,
        "rho_kg_m3": 2500.0,
        "cp_J_kgK": 1600.0,
        "E_Pa": 5.0e9,
        "nu": 0.40,
        "alpha_cte_1_K": 4.0e-5,
    },
    "copper": {
        "kappa0_W_mK": 400.0,
        "T_ref_K": 300.0,
        "alpha": 0.0,
        "rho_kg_m3": 8960.0,
        "cp_J_kgK": 385.0,
        "E_Pa": 1.1e11,
        "nu": 0.35,
        "alpha_cte_1_K": 1.7e-5,
    },
    "organic": {
        "kappa0_W_mK": 0.5,
        "T_ref_K": 300.0,
        "alpha": 0.0,
        "rho_kg_m3": 1800.0,
        "cp_J_kgK": 1000.0,
        "E_Pa": 2.5e10,
        "nu": 0.30,
        "alpha_cte_1_K": 1.5e-5,
    },
}

DEFAULT_PACKAGE = {
    "tim_um": 20.0,
    "t_spreader_mm": 1.0,
    "t_sink_mm": 6.9,
    "s_spreader_mm": None,  # None: bbox + 2*spreader_margin (not HotSpot 30 mm)
    "s_sink_mm": None,      # None: spreader + 2*sink_margin (not HotSpot 60 mm)
    "spreader_margin_mm": 11.0,
    "sink_margin_mm": 15.0,
    "r_convec_K_per_W": 0.1,
    "ambient_C": 45.0,
    "pitch_die_mm": 0.5,
    "pitch_coarse_mm": 3.0,
    "nz_die": 4,
    "nz_tim": 2,
    "nz_spreader": 4,
    "nz_sink": 6,
    "nz_substrate": 4,
    "nz_interposer": 3,
    "t_substrate_mm": 1.0,
    "t_interposer_um": 100.0,
    "interposer_margin_mm": 0.0,
    "k_chip": 130.0,
    "k_tim": 4.0,
    "k_spreader": 400.0,
    "k_sink": 400.0,
    "fill_die": False,
    "fill_material": "silicon",
    "fill_xy": "spreader",
    "tim_cover": "chips",
}


def _f(d, *keys, **kw):
    default = kw.get("default", None)
    for k in keys:
        if k in d and d[k] is not None:
            return float(d[k])
    return default


def parse_chip(raw, idx):
    """Return a chip dict in meters, origin = bottom-left of package."""
    name = raw.get("name", "chip%d" % idx)
    w = _f(raw, "w_mm", "width_mm")
    h = _f(raw, "h_mm", "height_mm")
    x = _f(raw, "x_mm", "x0_mm")
    y = _f(raw, "y_mm", "y0_mm")
    if w is None:
        w = _f(raw, "w_m", "width_m")
        if w is not None:
            w *= 1e3
    if h is None:
        h = _f(raw, "h_m", "height_m")
        if h is not None:
            h *= 1e3
    if x is None:
        x = _f(raw, "x0_m", "x_m")
        if x is not None:
            x *= 1e3
    if y is None:
        y = _f(raw, "y0_m", "y_m")
        if y is not None:
            y *= 1e3
    t_um = _f(raw, "t_um", "thickness_um")
    if t_um is None:
        t_mm = _f(raw, "t_mm", "thickness_mm")
        if t_mm is not None:
            t_um = t_mm * 1e3
        else:
            t_m = _f(raw, "t_m", "thickness_m")
            t_um = (t_m * 1e6) if t_m is not None else 150.0
    power = _f(raw, "power_W", "power", default=0.0)
    if None in (w, h, x, y):
        raise ValueError(
            "chip %s needs x_mm, y_mm, w_mm, h_mm (or *_m equivalents)" % name
        )
    if w <= 0 or h <= 0:
        raise ValueError("chip %s has non-positive size" % name)
    return {
        "name": str(name),
        "x0": x * 1e-3,
        "y0": y * 1e-3,
        "w": w * 1e-3,
        "h": h * 1e-3,
        "t": t_um * 1e-6,
        "power": power,
        "material": raw.get("material", "silicon"),
    }


def _overlap(a, b, tol=1e-12):
    ox = min(a["x0"] + a["w"], b["x0"] + b["w"]) - max(a["x0"], b["x0"])
    oy = min(a["y0"] + a["h"], b["y0"] + b["h"]) - max(a["y0"], b["y0"])
    return ox > tol and oy > tol


def _unique_sorted(vals, tol=1e-12):
    s = sorted(vals)
    out = [s[0]]
    for v in s[1:]:
        if abs(v - out[-1]) > tol:
            out.append(v)
    return out


def _regions_from_cuts(cuts, chips, pitch_fine, pitch_coarse):
    """Split [cuts[0], cuts[-1]] into regions; refine intervals that hit a chip."""
    regions = []
    for a, b in zip(cuts[:-1], cuts[1:]):
        span = b - a
        if span <= 1e-15:
            continue
        hits = False
        for c in chips:
            if c["x0"] < b - 1e-15 and c["x0"] + c["w"] > a + 1e-15:
                hits = True
                break
            # y-version passed via chips with x0/w rebound by caller
        n = int(math.ceil(span / (pitch_fine if hits else pitch_coarse)))
        n = max(n, 2 if hits else 1)
        regions.append({"x0_m": a, "x1_m": b, "n": n})
    return regions


def _xy_regions(cuts, chips, key0, keyw, pitch_fine, pitch_coarse, k0, k1):
    regions = []
    for a, b in zip(cuts[:-1], cuts[1:]):
        span = b - a
        if span <= 1e-15:
            continue
        hits = False
        for c in chips:
            c0 = c[key0]
            c1 = c[key0] + c[keyw]
            if c0 < b - 1e-15 and c1 > a + 1e-15:
                hits = True
                break
        n = int(math.ceil(span / (pitch_fine if hits else pitch_coarse)))
        n = max(n, 2 if hits else 1)
        regions.append({k0: a, k1: b, "n": n})
    return regions


def _hotspot_fillers(chips):
    """Tile the chip bounding box; cells not covered by a chip become 0 W fillers.

    HotSpot floorplans should fill a rectangle without holes or overlaps.
    """
    xs = []
    ys = []
    for c in chips:
        xs.extend([c["x0"], c["x0"] + c["w"]])
        ys.extend([c["y0"], c["y0"] + c["h"]])
    xs = _unique_sorted(xs)
    ys = _unique_sorted(ys)
    ox, oy = xs[0], ys[0]
    fillers = []
    n = 0
    for i in range(len(xs) - 1):
        for j in range(len(ys) - 1):
            x0, x1 = xs[i], xs[i + 1]
            y0, y1 = ys[j], ys[j + 1]
            cx = 0.5 * (x0 + x1)
            cy = 0.5 * (y0 + y1)
            covered = False
            for c in chips:
                if (c["x0"] - 1e-12) <= cx <= (c["x0"] + c["w"] + 1e-12) and (
                    c["y0"] - 1e-12
                ) <= cy <= (c["y0"] + c["h"] + 1e-12):
                    covered = True
                    break
            if not covered:
                n += 1
                fillers.append(
                    {
                        "name": "FILL_%d" % n,
                        "x0": x0 - ox,
                        "y0": y0 - oy,
                        "w": x1 - x0,
                        "h": y1 - y0,
                        "power": 0.0,
                    }
                )
    hs_chips = []
    for c in chips:
        hs_chips.append(
            {
                "name": c["name"],
                "x0": c["x0"] - ox,
                "y0": c["y0"] - oy,
                "w": c["w"],
                "h": c["h"],
                "power": c["power"],
            }
        )
    return hs_chips, fillers, ox, oy


def expand_chips_layout(d):
    """Return a new dict with domain/stack/blocks/probes/grid filled in."""
    out = copy.deepcopy(d)
    pkg_in = dict(DEFAULT_PACKAGE)
    pkg_in.update(d.get("package") or {})
    chips = [parse_chip(c, i) for i, c in enumerate(d["chips"])]
    if not chips:
        raise ValueError("chips[] is empty")
    for i in range(len(chips)):
        for j in range(i + 1, len(chips)):
            if _overlap(chips[i], chips[j]):
                raise ValueError(
                    "chips overlap: %s vs %s" % (chips[i]["name"], chips[j]["name"])
                )

    xmin = min(c["x0"] for c in chips)
    ymin = min(c["y0"] for c in chips)
    xmax = max(c["x0"] + c["w"] for c in chips)
    ymax = max(c["y0"] + c["h"] for c in chips)
    t_die_max = max(c["t"] for c in chips)
    t_sub = float(pkg_in.get("t_substrate_mm", 0.0) or 0.0) * 1e-3
    t_ip = float(pkg_in.get("t_interposer_um", 0.0) or 0.0) * 1e-6
    t_tim = float(pkg_in["tim_um"]) * 1e-6
    t_sp = float(pkg_in["t_spreader_mm"]) * 1e-3
    t_sink = float(pkg_in["t_sink_mm"]) * 1e-3
    z_sub0 = 0.0
    z_sub1 = t_sub
    z_ip0 = z_sub1
    z_ip1 = z_ip0 + t_ip
    z_die0 = z_ip1
    z_die1 = z_die0 + t_die_max
    z_tim_top = z_die1 + t_tim
    z_sp_top = z_tim_top + t_sp
    z_sink_top = z_sp_top + t_sink

    bbox_w = xmax - xmin
    bbox_h = ymax - ymin
    margin_sp = float(pkg_in["spreader_margin_mm"]) * 1e-3
    margin_sk = float(pkg_in["sink_margin_mm"]) * 1e-3
    # Package grows from the chip-cluster bounding box, then the cluster
    # is shifted so it sits in the center. User x_mm/y_mm are layout-relative.
    sp_w = bbox_w + 2.0 * margin_sp
    sp_h = bbox_h + 2.0 * margin_sp
    sink_w = sp_w + 2.0 * margin_sk
    sink_h = sp_h + 2.0 * margin_sk
    if pkg_in.get("s_spreader_mm"):
        side = float(pkg_in["s_spreader_mm"]) * 1e-3
        sp_w = sp_h = max(side, bbox_w, bbox_h)
    if pkg_in.get("s_sink_mm"):
        side = float(pkg_in["s_sink_mm"]) * 1e-3
        sink_w = sink_h = max(side, sp_w, sp_h)

    # frame_mm is optional. Chip x/y are then positions inside a fixed square
    # (HotSpot's floorplan), and that square is centered on the sink. Absent,
    # the cluster bbox is centered, which is the old behavior.
    frame_m = None
    if pkg_in.get("frame_mm"):
        frame_m = float(pkg_in["frame_mm"]) * 1e-3
        for c in chips:
            if (
                c["x0"] < -1e-6
                or c["y0"] < -1e-6
                or c["x0"] + c["w"] > frame_m + 1e-6
                or c["y0"] + c["h"] > frame_m + 1e-6
            ):
                raise ValueError(
                    "chip %s sits outside the %.3f mm frame"
                    % (c["name"], frame_m * 1e3)
                )
        shift_x = 0.5 * (sink_w - frame_m)
        shift_y = 0.5 * (sink_h - frame_m)
    else:
        shift_x = 0.5 * (sink_w - bbox_w) - xmin
        shift_y = 0.5 * (sink_h - bbox_h) - ymin
    for c in chips:
        c["x0"] += shift_x
        c["y0"] += shift_y
    xmin += shift_x
    ymin += shift_y
    xmax += shift_x
    ymax += shift_y
    frame_x0 = shift_x if frame_m is not None else None
    frame_y0 = shift_y if frame_m is not None else None
    sp_x0 = xmin - margin_sp
    sp_y0 = ymin - margin_sp
    if pkg_in.get("s_spreader_mm"):
        sp_x0 = 0.5 * (sink_w - sp_w)
        sp_y0 = 0.5 * (sink_h - sp_h)
    ip_margin = float(pkg_in.get("interposer_margin_mm", 0.0) or 0.0) * 1e-3
    ip_x0 = xmin - ip_margin
    ip_y0 = ymin - ip_margin
    ip_w = bbox_w + 2.0 * ip_margin
    ip_h = bbox_h + 2.0 * ip_margin

    r_convec = float(pkg_in["r_convec_K_per_W"])
    h = 1.0 / (r_convec * sink_w * sink_h)
    amb_c = float(pkg_in.get("ambient_C", 45.0))
    if "ambient_K" in d:
        amb_k = float(d["ambient_K"])
    else:
        amb_k = amb_c + 273.15

    pitch_f = float(pkg_in["pitch_die_mm"]) * 1e-3
    pitch_c = float(pkg_in["pitch_coarse_mm"]) * 1e-3
    x_cuts = [0.0, sink_w, sp_x0, sp_x0 + sp_w]
    y_cuts = [0.0, sink_h, sp_y0, sp_y0 + sp_h]
    if frame_m is not None:
        x_cuts.extend([frame_x0, frame_x0 + frame_m])
        y_cuts.extend([frame_y0, frame_y0 + frame_m])
    if t_ip > 0.0:
        x_cuts.extend([ip_x0, ip_x0 + ip_w])
        y_cuts.extend([ip_y0, ip_y0 + ip_h])
    for c in chips:
        x_cuts.extend([c["x0"], c["x0"] + c["w"]])
        y_cuts.extend([c["y0"], c["y0"] + c["h"]])
    x_cuts = _unique_sorted(x_cuts)
    y_cuts = _unique_sorted(y_cuts)
    x_regions = _xy_regions(x_cuts, chips, "x0", "w", pitch_f, pitch_c, "x0_m", "x1_m")
    y_regions = _xy_regions(y_cuts, chips, "y0", "h", pitch_f, pitch_c, "y0_m", "y1_m")

    nz_die = int(pkg_in["nz_die"])
    nz_tim = int(pkg_in["nz_tim"])
    stack = []
    if t_sub > 0.0:
        stack.append(
            {
                "name": "substrate",
                "material": "organic",
                "x0_m": xmin,
                "y0_m": ymin,
                "w_m": bbox_w,
                "h_m": bbox_h,
                "z0_m": z_sub0,
                "z1_m": z_sub1,
                "nz": int(pkg_in.get("nz_substrate", 4)),
                "has_power": False,
            }
        )
    if t_ip > 0.0:
        stack.append(
            {
                "name": "interposer",
                "material": "silicon",
                "x0_m": ip_x0,
                "y0_m": ip_y0,
                "w_m": ip_w,
                "h_m": ip_h,
                "z0_m": z_ip0,
                "z1_m": z_ip1,
                "nz": int(pkg_in.get("nz_interposer", 3)),
                "has_power": False,
            }
        )
    fill_die = bool(pkg_in.get("fill_die", False))
    fill_mat = pkg_in.get("fill_material", "silicon")
    fill_xy = pkg_in.get("fill_xy", "spreader")
    if fill_die:
        if frame_m is not None:
            fx0, fy0, fw, fh = frame_x0, frame_y0, frame_m, frame_m
        elif fill_xy == "sink":
            fx0, fy0, fw, fh = 0.0, 0.0, sink_w, sink_h
        elif fill_xy == "chips":
            fx0, fy0, fw, fh = xmin, ymin, bbox_w, bbox_h
        else:
            fx0, fy0, fw, fh = sp_x0, sp_y0, sp_w, sp_h
        stack.append(
            {
                "name": "die_fill",
                "material": fill_mat,
                "x0_m": fx0,
                "y0_m": fy0,
                "w_m": fw,
                "h_m": fh,
                "z0_m": z_die0,
                "z1_m": z_die1,
                "nz": nz_die,
                "has_power": False,
            }
        )
    for c in chips:
        stack.append(
            {
                "name": "die_" + c["name"],
                "material": c["material"],
                "x0_m": c["x0"],
                "y0_m": c["y0"],
                "w_m": c["w"],
                "h_m": c["h"],
                "z0_m": z_die0,
                "z1_m": z_die0 + c["t"],
                "nz": nz_die,
                "has_power": True,
            }
        )
    # Substrate, interposer and TIM follow the chip-layer bbox by default.
    tim_cover = pkg_in.get("tim_cover", "chips")
    if tim_cover in ("none", "per_chip"):
        for c in chips:
            stack.append(
                {
                    "name": "tim_" + c["name"],
                    "material": "tim",
                    "x0_m": c["x0"],
                    "y0_m": c["y0"],
                    "w_m": c["w"],
                    "h_m": c["h"],
                    "z0_m": z_die0 + c["t"],
                    "z1_m": z_tim_top,
                    "nz": nz_tim,
                    "has_power": False,
                }
            )
    else:
        if tim_cover == "sink":
            tx0, ty0, tw, th = 0.0, 0.0, sink_w, sink_h
        elif tim_cover == "spreader":
            tx0, ty0, tw, th = sp_x0, sp_y0, sp_w, sp_h
        elif tim_cover == "frame":
            if frame_m is None:
                raise ValueError("tim_cover frame needs package.frame_mm")
            tx0, ty0, tw, th = frame_x0, frame_y0, frame_m, frame_m
        else:
            tx0, ty0, tw, th = xmin, ymin, bbox_w, bbox_h
        stack.append(
            {
                "name": "tim_blanket",
                "material": "tim",
                "x0_m": tx0,
                "y0_m": ty0,
                "w_m": tw,
                "h_m": th,
                "z0_m": z_die1,
                "z1_m": z_tim_top,
                "nz": nz_tim,
                "has_power": False,
            }
        )
        for c in chips:
            if c["t"] < t_die_max - 1e-12:
                stack.append(
                    {
                        "name": "tim_gap_" + c["name"],
                        "material": "tim",
                        "x0_m": c["x0"],
                        "y0_m": c["y0"],
                        "w_m": c["w"],
                        "h_m": c["h"],
                        "z0_m": z_die0 + c["t"],
                        "z1_m": z_die1,
                        "nz": max(1, nz_tim),
                        "has_power": False,
                    }
                )
    stack.append(
        {
            "name": "spreader",
            "material": "copper",
            "x0_m": sp_x0,
            "y0_m": sp_y0,
            "w_m": sp_w,
            "h_m": sp_h,
            "z0_m": z_tim_top,
            "z1_m": z_sp_top,
            "nz": int(pkg_in["nz_spreader"]),
            "has_power": False,
        }
    )
    stack.append(
        {
            "name": "sink",
            "material": "copper",
            "x0_m": 0.0,
            "y0_m": 0.0,
            "w_m": sink_w,
            "h_m": sink_h,
            "z0_m": z_sp_top,
            "z1_m": z_sink_top,
            "nz": int(pkg_in["nz_sink"]),
            "has_power": False,
        }
    )

    blocks = []
    probes = []
    for c in chips:
        blocks.append(
            {
                "name": c["name"],
                "x0_m": c["x0"],
                "y0_m": c["y0"],
                "w_m": c["w"],
                "h_m": c["h"],
                "z0_m": z_die0,
                "z1_m": z_die0 + c["t"],
                "power_W": c["power"],
                "material": c["material"],
            }
        )
        probes.append(
            {
                "name": c["name"] + "_center",
                "x_m": c["x0"] + 0.5 * c["w"],
                "y_m": c["y0"] + 0.5 * c["h"],
                "z_m": z_die0 + 0.5 * c["t"],
            }
        )
    probes.append(
        {
            "name": "Sink_top_center",
            "x_m": 0.5 * sink_w,
            "y_m": 0.5 * sink_h,
            "z_m": z_sink_top,
        }
    )
    probes.append(
        {
            "name": "Spreader_center",
            "x_m": sp_x0 + 0.5 * sp_w,
            "y_m": sp_y0 + 0.5 * sp_h,
            "z_m": 0.5 * (z_tim_top + z_sp_top),
        }
    )
    if t_ip > 0.0:
        probes.append(
            {
                "name": "Interposer_center",
                "x_m": ip_x0 + 0.5 * ip_w,
                "y_m": ip_y0 + 0.5 * ip_h,
                "z_m": 0.5 * (z_ip0 + z_ip1),
            }
        )
    if t_sub > 0.0:
        probes.append(
            {
                "name": "Substrate_center",
                "x_m": xmin + 0.5 * bbox_w,
                "y_m": ymin + 0.5 * bbox_h,
                "z_m": 0.5 * (z_sub0 + z_sub1),
            }
        )
    if tim_cover in ("none", "per_chip") and chips:
        c0 = chips[0]
        probes.append(
            {
                "name": "TIM_mid",
                "x_m": c0["x0"] + 0.5 * c0["w"],
                "y_m": c0["y0"] + 0.5 * c0["h"],
                "z_m": 0.5 * (z_die1 + z_tim_top),
            }
        )
    else:
        probes.append(
            {
                "name": "TIM_mid",
                "x_m": xmin + 0.5 * bbox_w,
                "y_m": ymin + 0.5 * bbox_h,
                "z_m": 0.5 * (z_die1 + z_tim_top),
            }
        )

    hs_chips, fillers, hs_ox, hs_oy = _hotspot_fillers(chips)
    hotspot_units = hs_chips + fillers

    materials = copy.deepcopy(DEFAULT_MATERIALS)
    materials.update(d.get("materials") or {})

    out["ambient_K"] = amb_k
    out["convection"] = {
        "face": "sink_top",
        "h_W_m2K": h,
        "r_convec_K_per_W": r_convec,
    }
    out["domain"] = {"width_m": sink_w, "height_m": sink_h}
    out["chip_origin_m"] = [0.0, 0.0]
    out["grid"] = {"x_regions": x_regions, "y_regions": y_regions}
    out["stack"] = stack
    out["blocks"] = blocks
    if not d.get("probes"):
        out["probes"] = probes
    out["material_default"] = d.get("material_default", "silicon")
    out["materials"] = materials
    out.setdefault(
        "solver",
        {
            "type": "fvm3d",
            "nonlinear": False,
            "leakage": False,
            "kirchhoff": False,
            "tolerance_K": 0.05,
            "max_iter": 20,
        },
    )
    out["hotspot"] = {
        "r_convec_K_per_W": r_convec,
        "t_chip_m": t_die_max,
        "k_chip": float(pkg_in["k_chip"]),
        "t_interface_m": t_tim,
        "k_interface": float(pkg_in["k_tim"]),
        "s_spreader_m": max(sp_w, sp_h),
        "t_spreader_m": t_sp,
        "k_spreader": float(pkg_in["k_spreader"]),
        "s_sink_m": max(sink_w, sink_h),
        "t_sink_m": t_sink,
        "k_sink": float(pkg_in["k_sink"]),
        "grid_rows": 64,
        "grid_cols": 64,
        "units": hotspot_units,
        "origin_m": [hs_ox, hs_oy],
    }
    out["_layout_meta"] = {
        "n_chips": len(chips),
        "bbox_mm": [bbox_w * 1e3, bbox_h * 1e3],
        "spreader_mm": [sp_w * 1e3, sp_h * 1e3],
        "sink_mm": [sink_w * 1e3, sink_h * 1e3],
        "frame_mm": None if frame_m is None else frame_m * 1e3,
        "shift_mm": [shift_x * 1e3, shift_y * 1e3],
        "h_W_m2K": h,
        "total_power_W": sum(c["power"] for c in chips),
        "fill_die": fill_die,
        "tim_cover": tim_cover,
        "t_substrate_mm": t_sub * 1e3,
        "t_interposer_um": t_ip * 1e6,
        "interposer_mm": [ip_w * 1e3, ip_h * 1e3] if t_ip > 0.0 else [0.0, 0.0],
    }
    out.setdefault("mechanics", {"bc": "substrate_bottom_fixed"})
    return out
