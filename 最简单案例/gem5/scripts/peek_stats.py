#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Print a few gem5 stats that McPAT actually cares about."""
from __future__ import print_function
import sys

WANT = [
    "simSeconds",
    "simInsts",
    "system.cpu.numCycles",
    "system.cpu.idleCycles",
    "system.cpu.ipc",
    "system.cpu.committedInsts",
    "system.cpu.numInsts",
    "system.cpu.numIntInsts",
    "system.cpu.numFpInsts",
    "system.cpu.numLoadInsts",
    "system.cpu.numStoreInsts",
    "system.cpu.numBranches",
    "system.cpu.icache.overallAccesses::total",
    "system.cpu.icache.overallMisses::total",
    "system.cpu.dcache.ReadReq.accesses::total",
    "system.cpu.dcache.WriteReq.accesses::total",
    "system.l2.overallAccesses::total",
    "system.mem_ctrl.readReqs",
]


def load_stats(path):
    d = {}
    for line in open(path, "r"):
        line = line.strip()
        if not line or line.startswith("-") or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) >= 2:
            try:
                d[parts[0]] = parts[1]
            except Exception:
                pass
    return d


def main():
    if len(sys.argv) < 2:
        print("usage: python peek_stats.py <stats.txt>")
        sys.exit(1)
    d = load_stats(sys.argv[1])
    print("file:", sys.argv[1])
    print("%-45s %s" % ("stat", "value"))
    print("-" * 60)
    for k in WANT:
        print("%-45s %s" % (k, d.get(k, "(missing)")))
    print()
    print("These counts become McPAT <stat> via gem5_stats_to_mcpat.py")


if __name__ == "__main__":
    main()
