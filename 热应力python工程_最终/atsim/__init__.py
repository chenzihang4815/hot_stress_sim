# -*- coding: utf-8 -*-
"""ATSim: chip-level steady thermal solver (Chapter-2 style)."""

from .geometry import load_case, CaseConfig
from .materials import Material
from .fvm3d import solve_fvm3d
from .global_local import solve_atsim3d
from .hotspot_io import write_hotspot_case
from .visualize import plot_midplane, save_probes_csv
from .compare import compare_midplanes, load_comsol_csv

__all__ = [
    "load_case",
    "CaseConfig",
    "Material",
    "solve_fvm3d",
    "solve_atsim3d",
    "write_hotspot_case",
    "plot_midplane",
    "save_probes_csv",
    "compare_midplanes",
    "load_comsol_csv",
]
