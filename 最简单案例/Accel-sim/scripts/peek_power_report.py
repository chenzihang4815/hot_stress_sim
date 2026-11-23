#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Parse AccelWattch accelwattch_power_report.log (real or fake)."""
from __future__ import print_function
import sys

SHOW = [
    "kernel_name",
    "kernel_avg_power",
    "gpu_avg_INTP,",
    "gpu_avg_FPUP,",
    "gpu_avg_RFP,",
    "gpu_avg_L2CP,",
    "gpu_avg_NOCP,",
    "gpu_avg_DRAMP,",
    "gpu_avg_IDLE_COREP,",
    "gpu_avg_CONSTP",
    "gpu_avg_STATICP",
]


def open_text(path):
    try:
        return open(path, "r", encoding="utf-8", errors="ignore")
    except TypeError:
        return open(path, "r")


def main():
    if len(sys.argv) < 2:
        print("usage: python peek_power_report.py <accelwattch_power_report.log>")
        sys.exit(1)
    path = sys.argv[1]
    kv = {}
    for line in open_text(path):
        line = line.strip().replace(" ", "")
        if "=" not in line:
            continue
        k, v = line.split("=", 1)
        kv[k] = v
    print("file:", path)
    print("%-22s %s" % ("field", "value"))
    print("-" * 40)
    for k in SHOW:
        print("%-22s %s" % (k, kv.get(k, "(missing)")))
    print()
    print("kernel_avg_power = chip-total average Watts for this kernel.")
    print("INTP/FPUP/RFP    = units inside SMs (no XY).")
    print("L2CP/NOCP/DRAMP  = uncore / memory. You still need a floorplan for thermal.")


if __name__ == "__main__":
    main()
