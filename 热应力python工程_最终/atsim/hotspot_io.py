# -*- coding: utf-8 -*-
"""Emit HotSpot floorplan / power / config for the shared Case JSON."""

from __future__ import print_function

import os


def _write(path, text):
    parent = os.path.dirname(path)
    if parent and not os.path.isdir(parent):
        os.makedirs(parent)
    with open(path, "w") as f:
        f.write(text)


HOTSPOT_CONFIG_TEMPLATE = """# Auto-generated from cases/case1_2d_ic.json
# Chip-level comparison uses Python/COMSOL (equivalent-h convection).
# HotSpot always includes TIM + spreader + sink; use it for trend check only.

-t_chip				{t_chip}
-k_chip				{k_chip}
-p_chip				1630300

-c_convec			140.4
-r_convec			{r_convec}
-s_sink				{s_sink}
-t_sink				{t_sink}
-k_sink				{k_sink}
-p_sink				3.55e6

-s_spreader			{s_spreader}
-t_spreader			{t_spreader}
-k_spreader			{k_spreader}
-p_spreader			3.55e6

-t_interface		{t_interface}
-k_interface		{k_interface}
-p_interface		4.0e6

-model_secondary	0
-ambient			{ambient}
-init_temp			{ambient}
-sampling_intvl		0.01
-base_proc_freq		3e+09
-dtm_used			0
-model_type			grid
-leakage_used		0
-package_model_used	0
-block_omit_lateral	0
-grid_rows			{grid_rows}
-grid_cols			{grid_cols}
-grid_map_mode		avg
-use_microfluidic_cooling 0
"""

MATERIALS_FILE = """# material name / type / k / volumetric heat capacity
silicon
solid
130.0
1630300

copper
solid
400.0
3.45e6
"""


def write_hotspot_case(case, out_dir):
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)

    hs = case.hotspot or {}
    flp_lines = [
        "# name  width  height  left-x  bottom-y",
        "# generated from %s" % (case.path or case.name),
        "",
    ]
    names = []
    powers = []
    units = hs.get("units")
    if units:
        for u in units:
            flp_lines.append(
                "{0}\t{1:.10e}\t{2:.10e}\t{3:.10e}\t{4:.10e}".format(
                    u["name"], u["w"], u["h"], u["x0"], u["y0"]
                )
            )
            names.append(u["name"])
            powers.append(float(u.get("power", 0.0)))
    else:
        for b in case.blocks:
            flp_lines.append(
                "{0}\t{1:.10e}\t{2:.10e}\t{3:.10e}\t{4:.10e}".format(
                    b.name, b.w, b.h, b.x0, b.y0
                )
            )
            names.append(b.name)
            powers.append(b.power)
    flp_path = os.path.join(out_dir, case.name + ".flp")
    _write(flp_path, "\n".join(flp_lines) + "\n")

    ptrace_path = os.path.join(out_dir, case.name + ".ptrace")
    _write(
        ptrace_path,
        "\t".join(names)
        + "\n"
        + "\t".join("{0:.8f}".format(p) for p in powers)
        + "\n",
    )

    cfg = HOTSPOT_CONFIG_TEMPLATE.format(
        t_chip=hs.get("t_chip_m", getattr(case, "chip_Lz", case.Lz)),
        k_chip=hs.get("k_chip", case.material().kappa0),
        r_convec=hs.get("r_convec_K_per_W", 0.1),
        s_sink=hs.get("s_sink_m", 0.06),
        t_sink=hs.get("t_sink_m", 0.0069),
        k_sink=hs.get("k_sink", 400.0),
        s_spreader=hs.get("s_spreader_m", 0.03),
        t_spreader=hs.get("t_spreader_m", 0.001),
        k_spreader=hs.get("k_spreader", 400.0),
        t_interface=hs.get("t_interface_m", 2e-5),
        k_interface=hs.get("k_interface", 4.0),
        ambient=case.T_amb,
        grid_rows=hs.get("grid_rows", 64),
        grid_cols=hs.get("grid_cols", 64),
    )
    cfg_path = os.path.join(out_dir, case.name + ".config")
    _write(cfg_path, cfg)
    mat_path = os.path.join(out_dir, case.name + ".materials")
    _write(mat_path, MATERIALS_FILE)

    rows = int(hs.get("grid_rows", 64))
    cols = int(hs.get("grid_cols", 64))
    # Windows-friendly .bat
    bat = (
        "@echo off\n"
        "set ROOT=%~dp0..\HotSpot-master\\HotSpot-master\n"
        "if not exist outputs mkdir outputs\n"
        "if exist \"%ROOT%\\hotspot.exe\" (\n"
        "  set HS=%ROOT%\\hotspot.exe\n"
        ") else if exist \"%ROOT%\\hotspot\" (\n"
        "  set HS=%ROOT%\\hotspot\n"
        ") else (\n"
        "  echo [error] compile HotSpot first: cd HotSpot-master\\HotSpot-master && make\n"
        "  exit /b 1\n"
        ")\n"
        "\"%HS%\" -c {name}.config -f {name}.flp -p {name}.ptrace "
        "-materials_file {name}.materials -model_type grid "
        "-grid_rows {rows} -grid_cols {cols} "
        "-steady_file outputs/{name}.steady "
        "-grid_steady_file outputs/{name}.grid.steady\n"
    ).format(name=case.name, rows=rows, cols=cols)
    _write(os.path.join(out_dir, "run_hotspot.bat"), bat)

    sh = (
        "#!/usr/bin/env bash\n"
        "set -e\n"
        "ROOT=\"$(cd \"$(dirname \"$0\")/../HotSpot-master/HotSpot-master\" && pwd)\"\n"
        "mkdir -p outputs\n"
        "HS=\"$ROOT/hotspot\"\n"
        "if [ ! -x \"$HS\" ]; then echo \"compile HotSpot: cd $ROOT && make\"; exit 1; fi\n"
        "\"$HS\" -c {name}.config -f {name}.flp -p {name}.ptrace "
        "-materials_file {name}.materials -model_type grid "
        "-grid_rows {rows} -grid_cols {cols} "
        "-steady_file outputs/{name}.steady "
        "-grid_steady_file outputs/{name}.grid.steady\n"
    ).format(name=case.name, rows=rows, cols=cols)
    _write(os.path.join(out_dir, "run_hotspot.sh"), sh)

    return {
        "flp": flp_path,
        "ptrace": ptrace_path,
        "config": cfg_path,
        "materials": mat_path,
    }
