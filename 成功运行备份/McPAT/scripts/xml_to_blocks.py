#!/usr/bin/env python3
"""Parse McPAT stdout into a CSV of named blocks (indent, area, power)."""
from __future__ import print_function
import re
import sys

KEYS = [
    "Area",
    "Peak Dynamic",
    "Subthreshold Leakage",
    "Gate Leakage",
    "Runtime Dynamic",
    "Peak Power",
    "Total Leakage",
]


def parse(path):
    lines = open(path, encoding="utf-8", errors="ignore").read().splitlines()
    stack = []
    blocks = []
    kv_re = re.compile(r"^(\s+)(.+?)\s*=\s*([0-9.eE+-]+)")
    header_re = re.compile(r"^(\s*)(.+?):\s*(.*)$")
    for line in lines:
        if not line.strip() or set(line.strip()) <= set("*"):
            continue
        km = kv_re.match(line)
        if km:
            ind = len(km.group(1))
            k = km.group(2).strip()
            v = float(km.group(3))
            for i in range(len(stack) - 1, -1, -1):
                if stack[i][0] < ind:
                    if k in KEYS:
                        stack[i][2][k] = v
                    break
            continue
        hm = header_re.match(line)
        if hm and "=" not in line:
            ind = len(hm.group(1))
            name = hm.group(2).strip()
            extra = hm.group(3).strip()
            if extra:
                name = name + " " + extra
            while stack and stack[-1][0] >= ind:
                blocks.append(stack.pop())
            stack.append([ind, name, {}])
    while stack:
        blocks.append(stack.pop())
    return blocks


def main():
    if len(sys.argv) < 2:
        print("usage: xml_to_blocks.py <mcpat_stdout.txt> [out.csv]")
        sys.exit(1)
    blocks = parse(sys.argv[1])
    out = sys.argv[2] if len(sys.argv) > 2 else sys.stdout
    fh = open(out, "w") if out is not sys.stdout else sys.stdout
    fh.write(
        "indent,block,area_mm2,peak_dynamic_W,subthreshold_leakage_W,"
        "gate_leakage_W,runtime_dynamic_W,peak_power_W,total_leakage_W\n"
    )
    n = 0
    for ind, name, v in blocks:
        if not v:
            continue
        n += 1
        fh.write(
            "%d,%s,%s,%s,%s,%s,%s,%s,%s\n"
            % (
                ind,
                name.replace(",", " "),
                v.get("Area", ""),
                v.get("Peak Dynamic", ""),
                v.get("Subthreshold Leakage", ""),
                v.get("Gate Leakage", ""),
                v.get("Runtime Dynamic", ""),
                v.get("Peak Power", ""),
                v.get("Total Leakage", ""),
            )
        )
    if fh is not sys.stdout:
        fh.close()
        print("wrote", n, "blocks ->", out)


if __name__ == "__main__":
    main()
