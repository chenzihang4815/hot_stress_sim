#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Fill McPAT XML <stat> fields from a gem5 stats.txt.

Does NOT convert gem5 -> Accel-Sim. Accel-Sim is a GPU tracer/simulator.

Usage (Python 2.7 / 3.6+):
  python gem5_stats_to_mcpat.py \\
      --template ../../McPAT/input/ARM_A9_2000.xml \\
      --stats    ../sample_m5out/stats.txt \\
      --out      ../sample_m5out/from_gem5.xml

Only <stat> values listed in MAP are overwritten. All <param> stay as in the
template (core count, tech node, cache geometry). Missing gem5 stats are
left unchanged and printed as warnings.

This is a teaching mapper for Timing/Atomic + classic caches. An O3 CPU has
many more McPAT fields (ROB, rename, IQ); those stay at template defaults
unless you extend MAP.
"""
from __future__ import print_function
import argparse
import os
import re
import sys

# (xml_component_id, xml_stat_name) -> list of gem5 stat names to try in order
MAP = {
    ("system", "total_cycles"): [
        "system.cpu.numCycles",
        "system.cpu0.numCycles",
    ],
    ("system", "busy_cycles"): [
        "system.cpu.numCycles",
        "system.cpu0.numCycles",
    ],
    ("system", "idle_cycles"): [
        "system.cpu.idleCycles",
        "system.cpu.numIdleCycles",
        "system.cpu0.idleCycles",
    ],
    ("system.core0", "total_cycles"): [
        "system.cpu.numCycles",
        "system.cpu0.numCycles",
    ],
    ("system.core0", "busy_cycles"): [
        "system.cpu.numCycles",
        "system.cpu0.numCycles",
    ],
    ("system.core0", "idle_cycles"): [
        "system.cpu.idleCycles",
        "system.cpu0.idleCycles",
    ],
    ("system.core0", "total_instructions"): [
        "system.cpu.committedInsts",
        "system.cpu.numInsts",
        "simInsts",
        "system.cpu.commitStats0.numInsts",
    ],
    ("system.core0", "committed_instructions"): [
        "system.cpu.committedInsts",
        "system.cpu.numInsts",
        "simInsts",
    ],
    ("system.core0", "int_instructions"): [
        "system.cpu.numIntInsts",
        "system.cpu.commitStats0.numIntInsts",
    ],
    ("system.core0", "committed_int_instructions"): [
        "system.cpu.numIntInsts",
        "system.cpu.commitStats0.numIntInsts",
    ],
    ("system.core0", "fp_instructions"): [
        "system.cpu.numFpInsts",
        "system.cpu.commitStats0.numFpInsts",
    ],
    ("system.core0", "committed_fp_instructions"): [
        "system.cpu.numFpInsts",
        "system.cpu.commitStats0.numFpInsts",
    ],
    ("system.core0", "load_instructions"): [
        "system.cpu.numLoadInsts",
        "system.cpu.commitStats0.numLoadInsts",
    ],
    ("system.core0", "store_instructions"): [
        "system.cpu.numStoreInsts",
        "system.cpu.commitStats0.numStoreInsts",
    ],
    ("system.core0", "branch_instructions"): [
        "system.cpu.numBranches",
        "system.cpu.commitStats0.numBranches",
        "system.cpu.branchPred.lookups",
    ],
    ("system.core0", "branch_mispredictions"): [
        "system.cpu.numBranchMispred",
        "system.cpu.branchPred.mispredicted",
        "system.cpu.branchPred.condIncorrect",
    ],
    ("system.core0", "pipeline_duty_cycle"): [
        "system.cpu.ipc",
    ],
    ("system.core0", "ialu_accesses"): [
        "system.cpu.numIntAluAccesses",
        "system.cpu.numIntInsts",
    ],
    ("system.core0", "fpu_accesses"): [
        "system.cpu.numFpAluAccesses",
        "system.cpu.numFpInsts",
    ],
    ("system.core0", "function_calls"): [
        "system.cpu.numCalls",
    ],
    ("system.core0.icache", "read_accesses"): [
        "system.cpu.icache.overallAccesses::total",
        "system.cpu.icache.demandAccesses::total",
        "system.cpu.icache.ReadReq.accesses::total",
    ],
    ("system.core0.icache", "read_misses"): [
        "system.cpu.icache.overallMisses::total",
        "system.cpu.icache.demandMisses::total",
        "system.cpu.icache.ReadReq.misses::total",
    ],
    ("system.core0.dcache", "read_accesses"): [
        "system.cpu.dcache.ReadReq.accesses::total",
        "system.cpu.dcache.overallAccesses::total",
    ],
    ("system.core0.dcache", "write_accesses"): [
        "system.cpu.dcache.WriteReq.accesses::total",
    ],
    ("system.core0.dcache", "read_misses"): [
        "system.cpu.dcache.ReadReq.misses::total",
        "system.cpu.dcache.overallMisses::total",
    ],
    ("system.core0.dcache", "write_misses"): [
        "system.cpu.dcache.WriteReq.misses::total",
    ],
    ("system.core0.itlb", "total_accesses"): [
        "system.cpu.itb.rdAccesses",
        "system.cpu.mmu.itb.rdAccesses",
    ],
    ("system.core0.itlb", "total_misses"): [
        "system.cpu.itb.rdMisses",
        "system.cpu.mmu.itb.rdMisses",
    ],
    ("system.core0.dtlb", "total_accesses"): [
        "system.cpu.dtb.rdAccesses",
        "system.cpu.mmu.dtb.rdAccesses",
    ],
    ("system.core0.dtlb", "total_misses"): [
        "system.cpu.dtb.rdMisses",
        "system.cpu.mmu.dtb.rdMisses",
    ],
    ("system.core0.BTB", "read_accesses"): [
        "system.cpu.numBranches",
        "system.cpu.branchPred.BTBLookups",
    ],
    ("system.L20", "read_accesses"): [
        "system.l2.ReadReq.accesses::total",
        "system.l2.overallAccesses::total",
        "system.cpu.l2cache.ReadReq.accesses::total",
    ],
    ("system.L20", "write_accesses"): [
        "system.l2.WriteReq.accesses::total",
        "system.cpu.l2cache.WriteReq.accesses::total",
    ],
    ("system.L20", "read_misses"): [
        "system.l2.ReadReq.misses::total",
        "system.l2.overallMisses::total",
    ],
    ("system.L20", "write_misses"): [
        "system.l2.WriteReq.misses::total",
    ],
    ("system.mc", "memory_reads"): [
        "system.mem_ctrl.readReqs",
        "system.mem_ctrl.numReads",
    ],
    ("system.mc", "memory_writes"): [
        "system.mem_ctrl.writeReqs",
        "system.mem_ctrl.numWrites",
    ],
}


def load_stats(path):
    d = {}
    for line in open(path, "r"):
        line = line.strip()
        if not line or line.startswith("-"):
            continue
        # name  value  # comment
        parts = line.split()
        if len(parts) < 2:
            continue
        name, val = parts[0], parts[1]
        d[name] = val
    return d


def lookup(stats, names):
    for n in names:
        if n in stats:
            return stats[n], n
    return None, None


def replace_stat_in_component(xml, component_id, stat_name, new_value):
    """Replace <stat name="stat_name" value="..."/> inside component id=..."""
    # Find the component opening tag
    open_re = re.compile(
        r'(<component\s+id="%s"[^>]*>)' % re.escape(component_id),
        re.I,
    )
    m = open_re.search(xml)
    if not m:
        return xml, False
    start = m.end()
    # Find matching close by scanning nested <component>
    depth = 1
    pos = start
    while depth > 0:
        nxt_open = xml.find("<component", pos)
        nxt_close = xml.find("</component>", pos)
        if nxt_close < 0:
            return xml, False
        if nxt_open >= 0 and nxt_open < nxt_close:
            depth += 1
            pos = nxt_open + 10
        else:
            depth -= 1
            if depth == 0:
                end = nxt_close
                break
            pos = nxt_close + 12
    body = xml[start:end]
    stat_re = re.compile(
        r'(<stat\s+name="%s"\s+value=")([^"]*)("/>)' % re.escape(stat_name)
    )
    new_body, n = stat_re.subn(r"\g<1>%s\g<3>" % new_value, body, count=1)
    if n == 0:
        return xml, False
    return xml[:start] + new_body + xml[end:], True


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--template", required=True, help="McPAT XML template")
    ap.add_argument("--stats", required=True, help="gem5 stats.txt")
    ap.add_argument("--out", required=True, help="output XML")
    args = ap.parse_args()

    stats = load_stats(args.stats)
    xml = open(args.template, "r").read()

    filled, missing, skipped = [], [], []
    for (cid, sname), gem5_names in MAP.items():
        val, src = lookup(stats, gem5_names)
        if val is None:
            missing.append((cid, sname, gem5_names[0]))
            continue
        xml2, ok = replace_stat_in_component(xml, cid, sname, val)
        if not ok:
            skipped.append((cid, sname, val))
        else:
            xml = xml2
            filled.append((cid, sname, src, val))

    # memory_accesses = reads + writes if both present
    r, _ = lookup(stats, MAP[("system.mc", "memory_reads")])
    w, _ = lookup(stats, MAP[("system.mc", "memory_writes")])
    if r is not None and w is not None:
        try:
            total = str(int(float(r) + float(w)))
            xml2, ok = replace_stat_in_component(
                xml, "system.mc", "memory_accesses", total
            )
            if ok:
                xml = xml2
                filled.append(("system.mc", "memory_accesses", "reads+writes", total))
        except ValueError:
            pass

    out_dir = os.path.dirname(args.out)
    if out_dir and not os.path.isdir(out_dir):
        os.makedirs(out_dir)
    open(args.out, "w").write(xml)

    print("template :", args.template)
    print("stats    :", args.stats)
    print("wrote    :", args.out)
    print("filled   :", len(filled))
    for cid, sname, src, val in filled:
        print("  %-22s %-32s <- %s = %s" % (cid, sname, src, val))
    if missing:
        print("kept template (gem5 stat missing):", len(missing))
        for cid, sname, g in missing:
            print("  %-22s %-32s (tried %s)" % (cid, sname, g))
    if skipped:
        print("component/stat not in template:", len(skipped))
        for cid, sname, val in skipped:
            print("  %-22s %-32s (had %s)" % (cid, sname, val))
    print()
    print("Next (Linux):")
    print("  ./mcpat -infile %s -print_level 5 > mcpat_out.txt" % args.out)


if __name__ == "__main__":
    main()
