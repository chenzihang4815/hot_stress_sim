#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Parse McPAT text output and print the Processor-level power/area block."""
from __future__ import print_function
import re
import sys

KEYS = [
    "Area",
    "Peak Power",
    "Total Leakage",
    "Peak Dynamic",
    "Subthreshold Leakage",
    "Gate Leakage",
    "Runtime Dynamic",
]


def parse(path):
    try:
        text = open(path, "r", encoding="utf-8", errors="ignore").read()
    except TypeError:
        text = open(path, "r").read()
    out = {}
    for key in KEYS:
        mm = re.search(r"%s\s*=\s*([0-9.eE+-]+)" % re.escape(key), text)
        if mm:
            out[key] = float(mm.group(1))
    return out


def main():
    if len(sys.argv) < 2:
        print("usage: python read_mcpat_output.py <mcpat_stdout.txt>")
        sys.exit(1)
    vals = parse(sys.argv[1])
    if not vals:
        print("no Processor block found in", sys.argv[1])
        sys.exit(2)
    print("file:", sys.argv[1])
    for k in KEYS:
        if k in vals:
            unit = "mm^2" if k == "Area" else "W"
            print("  %-22s %10.6g %s" % (k, vals[k], unit))
    if "Runtime Dynamic" in vals:
        print()
        print("Runtime Dynamic = this workload's average dynamic power.")
        print("Peak Power      = architectural worst case (can exceed TDP).")
        print("Leakage depends on XML <param name=\"temperature\"> (Kelvin).")


if __name__ == "__main__":
    main()
