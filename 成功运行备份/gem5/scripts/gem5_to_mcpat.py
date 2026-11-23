#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Build a *new* McPAT XML from one gem5 run (config.ini + stats.txt).

This is not "edit the official ARM_A9 template". It writes a McPAT 1.0
input whose cores/caches/clocks/activity come from gem5.

Intended for: SE + TimingSimpleCPU/AtomicSimpleCPU + classic L1[/L2].
Not a full O3 / Ruby / FS translator.

Usage:
  python3 gem5_to_mcpat.py --config m5out/config.ini --stats m5out/stats.txt \\
      --out mcpat.xml
  ./mcpat -infile mcpat.xml -print_level 5 -opt_for_clk 0

McPAT still needs process assumptions gem5 does not provide:
  --tech-nm (default 40) --temperature-k (default 350) --device-type (default 2=LOP)
"""
from __future__ import print_function

import argparse
import os
import sys
from collections import defaultdict


def parse_ini(path):
    sec = None
    cfg = defaultdict(dict)
    for raw in open(path):
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        if line.startswith("[") and line.endswith("]"):
            sec = line[1:-1]
            continue
        if "=" in line and sec is not None:
            k, v = line.split("=", 1)
            cfg[sec][k.strip()] = v.strip()
    return cfg


def parse_stats(path):
    d = {}
    for raw in open(path):
        line = raw.strip()
        if not line or line.startswith("-") or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) >= 2:
            d[parts[0]] = parts[1]
    return d


def fnum(d, *names, default=0.0):
    for n in names:
        if n in d:
            try:
                return float(d[n])
            except ValueError:
                continue
    return float(default)


def inum(d, *names, default=0):
    return int(fnum(d, *names, default=default))


def cpu_mhz(cfg, stats):
    sim_freq = fnum(stats, "simFreq", default=1e12)  # ticks / simulated second
    period = fnum(cfg.get("system.cpu_clk_domain", {}), "clock", default=0)
    if period <= 0:
        period = fnum(cfg.get("system.cpu_cluster.clk_domain", {}), "clock", default=0)
    if period <= 0:
        period = fnum(cfg.get("system.clk_domain", {}), "clock", default=1000)
    hz = sim_freq / period
    return max(1, int(round(hz / 1e6)))


def count_cpus(cfg):
    n = 0
    if "system.cpu" in cfg:
        n = 1
    i = 0
    while "system.cpu%d" % i in cfg:
        n = max(n, i + 1)
        i += 1
    cluster = cfg.get("system.cpu_cluster", {})
    if "cpus" in str(cluster):
        pass
    return max(1, n)


def cpu_is_o3(cfg):
    t = cfg.get("system.cpu", {}).get("type", "")
    return "O3" in t or "o3" in t


def cache_cfg(section, line, default_size, default_assoc, default_lat):
    size = inum(section, "size", default=default_size)
    assoc = inum(section, "assoc", default=default_assoc)
    lat = inum(section, "data_latency", "tag_latency", default=default_lat)
    mshr = inum(section, "mshrs", default=4)
    wb = inum(section, "write_buffers", default=4)
    return size, line, assoc, lat, mshr, wb


def first_stat(stats, names, default=0):
    for n in names:
        if n in stats:
            try:
                return float(stats[n])
            except ValueError:
                continue
    return float(default)


def emit_xml(cfg, stats, args):
    line = inum(cfg.get("system", {}), "cache_line_size", default=64)
    mhz = cpu_mhz(cfg, stats)
    ncores = count_cpus(cfg)
    o3 = cpu_is_o3(cfg)
    machine_type = 0 if o3 else 1
    vdd = fnum(cfg.get("system.cpu_voltage_domain", {}), "voltage", default=1.0)
    if vdd <= 0:
        vdd = fnum(cfg.get("system.voltage_domain", {}), "voltage", default=1.0)

    ic = cfg.get("system.cpu.icache", cfg.get("system.cpu_cluster.cpus.icache", {}))
    dc = cfg.get("system.cpu.dcache", cfg.get("system.cpu_cluster.cpus.dcache", {}))
    l2 = cfg.get("system.l2", cfg.get("system.cpu_cluster.l2", {}))
    has_l2 = bool(l2.get("size"))
    itb = cfg.get(
        "system.cpu.mmu.itb",
        cfg.get(
            "system.cpu.itb",
            cfg.get("system.cpu_cluster.cpus.mmu.itb", {}),
        ),
    )
    dtb = cfg.get(
        "system.cpu.mmu.dtb",
        cfg.get(
            "system.cpu.dtb",
            cfg.get("system.cpu_cluster.cpus.mmu.dtb", {}),
        ),
    )

    ic_size, _, ic_assoc, ic_lat, ic_mshr, _ = cache_cfg(ic, line, 32768, 2, 2)
    dc_size, _, dc_assoc, dc_lat, dc_mshr, dc_wb = cache_cfg(dc, line, 65536, 2, 2)
    if has_l2:
        l2_size, _, l2_assoc, l2_lat, l2_mshr, l2_wb = cache_cfg(l2, line, 2097152, 8, 20)
    else:
        l2_size = l2_assoc = l2_lat = l2_mshr = l2_wb = 0

    tlb_i = inum(itb, "size", default=64)
    tlb_d = inum(dtb, "size", default=64)

    P = "system.cpu_cluster.cpus"
    cycles = first_stat(
        stats,
        [
            "system.cpu.numCycles",
            "system.cpu0.numCycles",
            P + ".numCycles",
        ],
        0,
    )
    insts = first_stat(
        stats,
        [
            "system.cpu.commitStats0.numInsts",
            "system.cpu.committedInsts",
            P + ".commitStats0.numInsts",
            "simInsts",
        ],
        0,
    )
    int_i = first_stat(
        stats,
        [
            "system.cpu.commitStats0.numIntInsts",
            "system.cpu.numIntInsts",
            P + ".commitStats0.numIntInsts",
        ],
        insts,
    )
    fp_i = first_stat(
        stats,
        [
            "system.cpu.commitStats0.numFpInsts",
            "system.cpu.numFpInsts",
            P + ".commitStats0.numFpInsts",
        ],
        0,
    )
    ld_i = first_stat(
        stats,
        [
            "system.cpu.commitStats0.numLoadInsts",
            "system.cpu.numLoadInsts",
            P + ".commitStats0.numLoadInsts",
        ],
        0,
    )
    st_i = first_stat(
        stats,
        [
            "system.cpu.commitStats0.numStoreInsts",
            "system.cpu.numStoreInsts",
            P + ".commitStats0.numStoreInsts",
        ],
        0,
    )
    ipc = first_stat(
        stats,
        [
            "system.cpu.ipc",
            "system.cpu.commitStats0.ipc",
            P + ".ipc",
            P + ".commitStats0.ipc",
        ],
        0,
    )
    if ipc > 1:
        ipc = 1.0
    alu = first_stat(
        stats,
        [
            "system.cpu.issuedInstType_0::IntAlu",
            "system.cpu.numIntAluAccesses",
            P + ".issuedInstType_0::IntAlu",
            P + ".commitStats0.committedInstType::IntAlu",
        ],
        int_i,
    )
    mul = first_stat(
        stats,
        [
            "system.cpu.issuedInstType_0::IntMult",
            "system.cpu.issuedInstType_0::IntDiv",
            P + ".issuedInstType_0::IntMult",
        ],
        0,
    )
    fpu = first_stat(
        stats,
        ["system.cpu.numFpAluAccesses", P + ".numFpAluAccesses"],
        fp_i,
    )
    br = first_stat(
        stats,
        [
            "system.cpu.numBranches",
            "system.cpu.commitStats0.numBranches",
            "system.cpu.branchPred.lookups",
        ],
        0,
    )
    brm = first_stat(
        stats,
        [
            "system.cpu.numBranchMispred",
            "system.cpu.branchPred.condIncorrect",
        ],
        0,
    )

    ic_a = first_stat(
        stats,
        [
            "system.cpu.icache.overallAccesses::total",
            "system.cpu.icache.demandAccesses::total",
        ],
        insts,
    )
    ic_m = first_stat(
        stats,
        [
            "system.cpu.icache.overallMisses::total",
            "system.cpu.icache.demandMisses::total",
        ],
        0,
    )
    dc_r = first_stat(
        stats,
        [
            "system.cpu.dcache.ReadReq.accesses::total",
            "system.cpu.dcache.overallAccesses::total",
        ],
        ld_i,
    )
    dc_w = first_stat(stats, ["system.cpu.dcache.WriteReq.accesses::total"], st_i)
    dc_rm = first_stat(
        stats,
        [
            "system.cpu.dcache.ReadReq.misses::total",
            "system.cpu.dcache.overallMisses::total",
        ],
        0,
    )
    dc_wm = first_stat(stats, ["system.cpu.dcache.WriteReq.misses::total"], 0)
    l2_a = first_stat(
        stats,
        ["system.l2.overallAccesses::total", "system.l2.demandAccesses::total"],
        0,
    )
    l2_m = first_stat(
        stats,
        ["system.l2.overallMisses::total", "system.l2.demandMisses::total"],
        0,
    )
    mem_r = first_stat(
        stats,
        [
            "system.mem_ctrls.readReqs",
            "system.mem_ctrl.readReqs",
            "system.mem_ctrls.dram.numReads::total",
        ],
        l2_m,
    )
    mem_w = first_stat(
        stats,
        [
            "system.mem_ctrls.writeReqs",
            "system.mem_ctrl.writeReqs",
            "system.mem_ctrls.dram.writeBursts",
        ],
        0,
    )

    idle = first_stat(
        stats,
        [
            "system.cpu.idleCycles",
            "system.cpu.numIdleCycles",
            P + ".idleCycles",
        ],
        0,
    )
    busy = max(0.0, cycles - idle)
    duty = min(1.0, max(ipc, 0.05)) if cycles else 0.5

    notes = {
        "ncores": ncores,
        "mhz": mhz,
        "o3": o3,
        "has_l2": has_l2,
        "ic_size": ic_size,
        "dc_size": dc_size,
        "l2_size": l2_size,
        "cycles": int(cycles),
        "insts": int(insts),
        "tech": args.tech_nm,
        "temp": args.temperature_k,
    }

    def p(name, val, indent=2):
        return "%s<param name=\"%s\" value=\"%s\"/>\n" % ("\t" * indent, name, val)

    def s(name, val, indent=2):
        if isinstance(val, float) and val == int(val):
            val = int(val)
        return "%s<stat name=\"%s\" value=\"%s\"/>\n" % ("\t" * indent, name, val)

    ic_config = "%d,%d,%d,1,1,%d,32,0" % (ic_size, line, ic_assoc, ic_lat)
    dc_config = "%d,%d,%d,1,1,%d,32,1" % (dc_size, line, dc_assoc, dc_lat)
    l2_config = "%d,%d,%d,1,8,%d,32,1" % (l2_size or 2097152, line, l2_assoc or 8, l2_lat or 20)

    xml = []
    xml.append('<?xml version="1.0" ?>\n')
    xml.append("<!-- gem5_to_mcpat.py\n")
    xml.append("     cores=%s clock=%sMHz o3=%s L2=%s I$=%s D$=%s L2size=%s\n" % (
        ncores, mhz, o3, has_l2, ic_size, dc_size, l2_size))
    xml.append("     cycles=%s insts=%s tech=%snm T=%sK (tech/T not from gem5)\n" % (
        int(cycles), int(insts), args.tech_nm, args.temperature_k))
    xml.append("-->\n")
    xml.append('<component id="root" name="root">\n')
    xml.append('\t<component id="system" name="system">\n')
    xml.append(p("number_of_cores", ncores))
    xml.append(p("number_of_L1Directories", 0))
    xml.append(p("number_of_L2Directories", 0))
    xml.append(p("number_of_L2s", 1 if has_l2 else 0))
    xml.append(p("Private_L2", 0))
    xml.append(p("number_of_L3s", 0))
    xml.append(p("number_of_NoCs", 1))
    xml.append(p("homogeneous_cores", 1))
    xml.append(p("homogeneous_L2s", 1))
    xml.append(p("homogeneous_L1Directories", 1))
    xml.append(p("homogeneous_L2Directories", 1))
    xml.append(p("homogeneous_L3s", 1))
    xml.append(p("homogeneous_ccs", 1))
    xml.append(p("homogeneous_NoCs", 1))
    xml.append(p("core_tech_node", args.tech_nm))
    xml.append(p("target_core_clockrate", mhz))
    xml.append(p("temperature", args.temperature_k))
    xml.append(p("number_cache_levels", 2 if has_l2 else 1))
    xml.append(p("interconnect_projection_type", 1))
    xml.append(p("device_type", args.device_type))
    xml.append(p("longer_channel_device", 1))
    xml.append(p("power_gating", 0))
    xml.append(p("Embedded", 1 if args.embedded else 0))
    xml.append(p("opt_clockrate", 0))
    xml.append(p("machine_bits", args.machine_bits))
    xml.append(p("virtual_address_width", args.machine_bits))
    xml.append(p("physical_address_width", args.machine_bits))
    xml.append(p("virtual_memory_page_size", 4096))
    xml.append(s("total_cycles", cycles))
    xml.append(s("idle_cycles", idle))
    xml.append(s("busy_cycles", busy))

    xml.append('\t\t<component id="system.core0" name="core0">\n')
    xml.append(p("clock_rate", mhz, 3))
    xml.append(p("vdd", vdd, 3))
    xml.append(p("opt_local", 0, 3))
    xml.append(p("instruction_length", 32, 3))
    xml.append(p("opcode_width", 7, 3))
    xml.append(p("x86", 1 if args.machine_bits == 64 and "X86" in cfg.get("system.cpu", {}).get("type", "") else 0, 3))
    xml.append(p("micro_opcode_width", 8, 3))
    xml.append(p("machine_type", machine_type, 3))
    xml.append(p("number_hardware_threads", inum(cfg.get("system.cpu", {}), "numThreads", default=1), 3))
    width = 4 if o3 else 1
    xml.append(p("fetch_width", width, 3))
    xml.append(p("number_instruction_fetch_ports", 1, 3))
    xml.append(p("decode_width", width, 3))
    xml.append(p("issue_width", width, 3))
    xml.append(p("peak_issue_width", width, 3))
    xml.append(p("commit_width", width, 3))
    xml.append(p("fp_issue_width", 1, 3))
    xml.append(p("prediction_width", 0 if not o3 else 1, 3))
    xml.append(p("pipelines_per_core", "1,1", 3))
    xml.append(p("pipeline_depth", "5,5" if not o3 else "10,12", 3))
    xml.append(p("ALU_per_core", 1 if not o3 else 3, 3))
    xml.append(p("MUL_per_core", 1, 3))
    xml.append(p("FPU_per_core", 1, 3))
    xml.append(p("instruction_buffer_size", 16, 3))
    xml.append(p("decoded_stream_buffer_size", 16, 3))
    xml.append(p("instruction_window_scheme", 0, 3))
    xml.append(p("instruction_window_size", 8 if not o3 else 32, 3))
    xml.append(p("fp_instruction_window_size", 8 if not o3 else 16, 3))
    xml.append(p("ROB_size", 0 if not o3 else 128, 3))
    xml.append(p("archi_Regs_IRF_size", 16 if args.machine_bits == 32 else 32, 3))
    xml.append(p("archi_Regs_FRF_size", 32, 3))
    xml.append(p("phy_Regs_IRF_size", 16 if not o3 else 128, 3))
    xml.append(p("phy_Regs_FRF_size", 32 if not o3 else 128, 3))
    xml.append(p("rename_scheme", 0, 3))
    xml.append(p("register_windows_size", 0, 3))
    xml.append(p("LSU_order", "inorder", 3))
    xml.append(p("store_buffer_size", 4, 3))
    xml.append(p("load_buffer_size", 0 if not o3 else 16, 3))
    xml.append(p("memory_ports", 1, 3))
    xml.append(p("RAS_size", 8, 3))
    xml.append(s("total_instructions", insts, 3))
    xml.append(s("int_instructions", int_i, 3))
    xml.append(s("fp_instructions", fp_i, 3))
    xml.append(s("branch_instructions", br, 3))
    xml.append(s("branch_mispredictions", brm, 3))
    xml.append(s("load_instructions", ld_i, 3))
    xml.append(s("store_instructions", st_i, 3))
    xml.append(s("committed_instructions", insts, 3))
    xml.append(s("committed_int_instructions", int_i, 3))
    xml.append(s("committed_fp_instructions", fp_i, 3))
    xml.append(s("pipeline_duty_cycle", ipc if ipc > 0 else duty, 3))
    xml.append(s("total_cycles", cycles, 3))
    xml.append(s("idle_cycles", idle, 3))
    xml.append(s("busy_cycles", busy, 3))
    xml.append(s("ROB_reads", insts if o3 else 0, 3))
    xml.append(s("ROB_writes", insts if o3 else 0, 3))
    xml.append(s("rename_reads", 0 if not o3 else insts * 2, 3))
    xml.append(s("rename_writes", 0 if not o3 else insts, 3))
    xml.append(s("fp_rename_reads", 0, 3))
    xml.append(s("fp_rename_writes", 0, 3))
    xml.append(s("inst_window_reads", insts, 3))
    xml.append(s("inst_window_writes", insts, 3))
    xml.append(s("inst_window_wakeup_accesses", 0 if not o3 else insts, 3))
    xml.append(s("fp_inst_window_reads", 0, 3))
    xml.append(s("fp_inst_window_writes", 0, 3))
    xml.append(s("fp_inst_window_wakeup_accesses", 0, 3))
    xml.append(s("int_regfile_reads", int_i, 3))
    xml.append(s("float_regfile_reads", fp_i, 3))
    xml.append(s("int_regfile_writes", int_i, 3))
    xml.append(s("float_regfile_writes", fp_i, 3))
    xml.append(s("function_calls", 0, 3))
    xml.append(s("context_switches", 0, 3))
    xml.append(s("ialu_accesses", alu, 3))
    xml.append(s("fpu_accesses", fpu, 3))
    xml.append(s("mul_accesses", mul, 3))
    xml.append(s("cdb_alu_accesses", alu, 3))
    xml.append(s("cdb_mul_accesses", mul, 3))
    xml.append(s("cdb_fpu_accesses", fpu, 3))
    xml.append(s("IFU_duty_cycle", duty, 3))
    xml.append(s("BR_duty_cycle", min(1.0, duty * 0.2), 3))
    xml.append(s("LSU_duty_cycle", duty, 3))
    xml.append(s("MemManU_I_duty_cycle", duty, 3))
    xml.append(s("MemManU_D_duty_cycle", duty, 3))
    xml.append(s("ALU_duty_cycle", duty, 3))
    xml.append(s("MUL_duty_cycle", 0.1, 3))
    xml.append(s("FPU_duty_cycle", 0.0 if fp_i == 0 else duty, 3))
    xml.append(s("ALU_cdb_duty_cycle", duty, 3))
    xml.append(s("MUL_cdb_duty_cycle", 0.1, 3))
    xml.append(s("FPU_cdb_duty_cycle", 0.0 if fp_i == 0 else duty, 3))
    xml.append(p("number_of_BPT", 1, 3))
    xml.append('\t\t\t<component id="system.core0.predictor" name="PBT">\n')
    xml.append(p("local_predictor_size", "10,3", 4))
    xml.append(p("local_predictor_entries", 4, 4))
    xml.append(p("global_predictor_entries", 256, 4))
    xml.append(p("global_predictor_bits", 2, 4))
    xml.append(p("chooser_predictor_entries", 256, 4))
    xml.append(p("chooser_predictor_bits", 2, 4))
    xml.append("\t\t\t</component>\n")
    xml.append('\t\t\t<component id="system.core0.itlb" name="itlb">\n')
    xml.append(p("number_entries", tlb_i, 4))
    xml.append(s("total_accesses", ic_a, 4))
    xml.append(s("total_misses", 0, 4))
    xml.append(s("conflicts", 0, 4))
    xml.append("\t\t\t</component>\n")
    xml.append('\t\t\t<component id="system.core0.icache" name="icache">\n')
    xml.append(p("icache_config", ic_config, 4))
    xml.append(p("buffer_sizes", "%d,4,4,0" % ic_mshr, 4))
    xml.append(s("read_accesses", ic_a, 4))
    xml.append(s("read_misses", ic_m, 4))
    xml.append(s("conflicts", 0, 4))
    xml.append("\t\t\t</component>\n")
    xml.append('\t\t\t<component id="system.core0.dtlb" name="dtlb">\n')
    xml.append(p("number_entries", tlb_d, 4))
    xml.append(s("total_accesses", dc_r + dc_w, 4))
    xml.append(s("total_misses", 0, 4))
    xml.append(s("conflicts", 0, 4))
    xml.append("\t\t\t</component>\n")
    xml.append('\t\t\t<component id="system.core0.dcache" name="dcache">\n')
    xml.append(p("dcache_config", dc_config, 4))
    xml.append(p("buffer_sizes", "%d,4,4,%d" % (dc_mshr, dc_wb), 4))
    xml.append(s("read_accesses", dc_r, 4))
    xml.append(s("write_accesses", dc_w, 4))
    xml.append(s("read_misses", dc_rm, 4))
    xml.append(s("write_misses", dc_wm, 4))
    xml.append(s("conflicts", 0, 4))
    xml.append("\t\t\t</component>\n")
    xml.append(p("number_of_BTB", 1, 3))
    xml.append('\t\t\t<component id="system.core0.BTB" name="BTB">\n')
    xml.append(p("BTB_config", "1024,4,2,1,1,1", 4))
    xml.append(s("read_accesses", br, 4))
    xml.append(s("write_accesses", 0, 4))
    xml.append("\t\t\t</component>\n")
    xml.append("\t\t</component>\n")

    # unused directory placeholders (number_of_* = 0)
    xml.append('\t\t<component id="system.L1Directory0" name="L1Directory0">\n')
    xml.append(p("Directory_type", 0, 3))
    xml.append(p("Dir_config", "2048,1,0,1,4,4,8", 3))
    xml.append(p("buffer_sizes", "8,8,8,8", 3))
    xml.append(p("clockrate", mhz, 3))
    xml.append(p("ports", "1,1,1", 3))
    xml.append(p("device_type", args.device_type, 3))
    xml.append(s("read_accesses", 0, 3))
    xml.append(s("write_accesses", 0, 3))
    xml.append(s("read_misses", 0, 3))
    xml.append(s("write_misses", 0, 3))
    xml.append(s("conflicts", 0, 3))
    xml.append(s("duty_cycle", 0.1, 3))
    xml.append("\t\t</component>\n")
    xml.append('\t\t<component id="system.L2Directory0" name="L2Directory0">\n')
    xml.append(p("Directory_type", 1, 3))
    xml.append(p("Dir_config", "1048576,16,16,1,2,100", 3))
    xml.append(p("buffer_sizes", "8,8,8,8", 3))
    xml.append(p("clockrate", mhz, 3))
    xml.append(p("ports", "1,1,1", 3))
    xml.append(p("device_type", args.device_type, 3))
    xml.append(s("read_accesses", 0, 3))
    xml.append(s("write_accesses", 0, 3))
    xml.append(s("read_misses", 0, 3))
    xml.append(s("write_misses", 0, 3))
    xml.append(s("conflicts", 0, 3))
    xml.append(s("duty_cycle", 0.1, 3))
    xml.append("\t\t</component>\n")

    xml.append('\t\t<component id="system.L20" name="L20">\n')
    xml.append(p("L2_config", l2_config, 3))
    xml.append(p("buffer_sizes", "%d,16,16,%d" % (l2_mshr or 16, l2_wb or 8), 3))
    xml.append(p("clockrate", mhz, 3))
    xml.append(p("ports", "1,1,1", 3))
    xml.append(p("device_type", args.device_type, 3))
    xml.append(s("read_accesses", l2_a if has_l2 else 0, 3))
    xml.append(s("write_accesses", 0, 3))
    xml.append(s("read_misses", l2_m if has_l2 else 0, 3))
    xml.append(s("write_misses", 0, 3))
    xml.append(s("conflicts", 0, 3))
    xml.append(s("duty_cycle", 1.0, 3))
    xml.append("\t\t</component>\n")

    xml.append('\t\t<component id="system.L30" name="L30">\n')
    xml.append(p("L3_config", "16777216,64,16,16,16,100,1", 3))
    xml.append(p("clockrate", 800, 3))
    xml.append(p("ports", "1,1,1", 3))
    xml.append(p("device_type", args.device_type, 3))
    xml.append(p("buffer_sizes", "16,16,16,16", 3))
    xml.append(s("read_accesses", 0, 3))
    xml.append(s("write_accesses", 0, 3))
    xml.append(s("read_misses", 0, 3))
    xml.append(s("write_misses", 0, 3))
    xml.append(s("conflicts", 0, 3))
    xml.append(s("duty_cycle", 1.0, 3))
    xml.append("\t\t</component>\n")

    xml.append('\t\t<component id="system.NoC0" name="noc0">\n')
    xml.append(p("clockrate", max(1, mhz // 2), 3))
    xml.append(p("type", 0, 3))
    xml.append(p("horizontal_nodes", 1, 3))
    xml.append(p("vertical_nodes", 1, 3))
    xml.append(p("has_global_link", 0, 3))
    xml.append(p("link_throughput", 1, 3))
    xml.append(p("link_latency", 1, 3))
    xml.append(p("input_ports", 1, 3))
    xml.append(p("output_ports", 1, 3))
    xml.append(p("flit_bits", 64, 3))
    xml.append(p("chip_coverage", 1, 3))
    xml.append(p("link_routing_over_percentage", 0.5, 3))
    xml.append(s("total_accesses", l2_a if has_l2 else dc_r + dc_w, 3))
    xml.append(s("duty_cycle", 0.2, 3))
    xml.append("\t\t</component>\n")

    xml.append('\t\t<component id="system.mc" name="mc">\n')
    xml.append(p("type", 1, 3))
    xml.append(p("mc_clock", 400, 3))
    xml.append(p("peak_transfer_rate", 6400, 3))
    xml.append(p("block_size", line, 3))
    xml.append(p("number_mcs", 1, 3))
    xml.append(p("memory_channels_per_mc", 1, 3))
    xml.append(p("number_ranks", 1, 3))
    xml.append(p("req_window_size_per_channel", 32, 3))
    xml.append(p("IO_buffer_size_per_channel", 32, 3))
    xml.append(p("databus_width", 64, 3))
    xml.append(p("addressbus_width", args.machine_bits, 3))
    xml.append(s("memory_accesses", mem_r + mem_w, 3))
    xml.append(s("memory_reads", mem_r, 3))
    xml.append(s("memory_writes", mem_w, 3))
    xml.append(p("withPHY", 1, 3))
    xml.append("\t\t</component>\n")

    for extra in (
        (
            "niu",
            [p("type", 1, 3), p("clockrate", 350, 3), p("number_units", 0, 3),
             s("duty_cycle", 1.0, 3), s("total_load_perc", 0.0, 3)],
        ),
        (
            "pcie",
            [p("type", 1, 3), p("withPHY", 1, 3), p("clockrate", 350, 3),
             p("number_units", 0, 3), p("num_channels", 8, 3),
             s("duty_cycle", 1.0, 3), s("total_load_perc", 0.0, 3)],
        ),
        (
            "flashc",
            [p("number_flashcs", 0, 3), p("type", 1, 3), p("withPHY", 1, 3),
             p("peak_transfer_rate", 200, 3), s("duty_cycle", 1.0, 3),
             s("total_load_perc", 0.0, 3)],
        ),
    ):
        xml.append('\t\t<component id="system.%s" name="%s">\n' % (extra[0], extra[0]))
        xml.extend(extra[1])
        xml.append("\t\t</component>\n")

    xml.append("\t</component>\n</component>\n")
    return "".join(xml), notes


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", required=True, help="gem5 config.ini")
    ap.add_argument("--stats", required=True, help="gem5 stats.txt")
    ap.add_argument("--out", required=True, help="McPAT XML to write")
    ap.add_argument("--tech-nm", type=int, default=40)
    ap.add_argument("--temperature-k", type=int, default=350)
    ap.add_argument("--device-type", type=int, default=2, help="0=HP 1=LSTP 2=LOP")
    ap.add_argument("--machine-bits", type=int, default=32)
    ap.add_argument("--embedded", type=int, default=1)
    args = ap.parse_args()

    cfg = parse_ini(args.config)
    stats = parse_stats(args.stats)
    xml, notes = emit_xml(cfg, stats, args)
    out_dir = os.path.dirname(args.out)
    if out_dir and not os.path.isdir(out_dir):
        os.makedirs(out_dir)
    open(args.out, "w").write(xml)
    print("wrote", args.out)
    print("  cores=%(ncores)s  %(mhz)s MHz  O3=%(o3)s  L2=%(has_l2)s" % notes)
    print("  I$=%(ic_size)s B  D$=%(dc_size)s B  L2=%(l2_size)s B" % notes)
    print("  cycles=%(cycles)s  insts=%(insts)s" % notes)
    print("  tech=%(tech)s nm  T=%(temp)s K  (not from gem5)" % notes)
    print("next:")
    print("  mcpat -infile %s -print_level 5 -opt_for_clk 0" % args.out)


if __name__ == "__main__":
    main()
