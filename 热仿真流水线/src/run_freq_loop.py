#!/usr/bin/env python3
"""Frequency loop after SA. Layout L* is frozen.

CPU and GPU share scale s = f/f0; CIM analog power and throughput do not scale.
Window: Tmax < 90 C and worst-normal-chip interior von Mises < 200 MPa.
s_min = 0.5: if even that exceeds T or sigma, abandon the whole scheme.
Among feasible s, the winner is the highest package throughput
    Theta_pkg = (1/3) Theta_CPU/ref + (1/3) Theta_GPU/ref + (1/3) Theta_CIM/ref
with IPC held fixed so CPU/GPU relative throughput = s, CIM analog = 1.

Does not overwrite SA run directories. Log: EXPERIMENT.md in the out dir.
"""

from __future__ import annotations

import json
import shutil
import sys
import time
from copy import deepcopy
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path("/home/czh/任务仿真产生功耗文件工具")
SRC = ROOT / "热仿真流水线" / "src"
sys.path.insert(0, str(SRC))

import sa_thermal_feedback as sa  # noqa: E402
from atsim_sa_proxy import eval_atsim_proxy, write_proxy_json  # noqa: E402
from hotspot_io import run_hotspot, write_hotspot_inputs  # noqa: E402
from run_sa_compare_scheme import pkg_overrides, vm_of  # noqa: E402

T_LIM_C = 90.0
SIG_LIM_MPA = 200.0
S_MIN = 0.5
S_MAX_TRY = 4.0

# Processor-level McPAT Total Leakage (W). Dynamic = layout power − leak.
CPU_LEAK_W = {
    "CPU_matmul": 0.78389,
    "CPU_threads": 0.78389,
    "CPU_armL2": 0.78389,
    "CPU_atomic": 0.6197,
    "CPU_starter": 0.6197,
}
# AccelWattch CONST + STATIC (W), unscaled with s. GPU_3070b is a second 3070.
GPU_STATIC_W = {
    "GPU_3070": 32.3252 + 22.5847,  # 54.9099
    "GPU_3070b": 32.3252 + 22.5847,
    "GPU_2060S": 32.3252 + 24.4467,  # 56.7719
    "GPU_TITANX": 32.3252 + 29.3748,  # 61.7000
}


def log_md(path: Path, line: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8") as f:
        f.write(line.rstrip() + "\n")
    print(line, flush=True)


def split_power(c: sa.Chip) -> tuple[float, float]:
    """Return (P_dyn_at_s1, P_unscaled)."""
    if c.kind == "CPU":
        leak = CPU_LEAK_W.get(c.name, 0.0)
        leak = min(leak, 0.95 * c.power_W)
        return (c.power_W - leak, leak)
    if c.kind == "GPU":
        st = GPU_STATIC_W.get(c.name, 0.4 * c.power_W)
        st = min(st, 0.95 * c.power_W)
        return (c.power_W - st, st)
    return (0.0, c.power_W)  # CIM analog: latency_s fixed


def scale_chips(chips: list[sa.Chip], s: float) -> list[sa.Chip]:
    out = deepcopy(chips)
    for c in out:
        dyn, uns = split_power(c)
        c.power_W = uns + s * dyn
    return out


def theta_pkg(s: float) -> dict:
    t_cpu = s
    t_gpu = s
    t_cim = 1.0
    pkg = (t_cpu + t_gpu + t_cim) / 3.0
    return {
        "Theta_CPU_rel": t_cpu,
        "Theta_GPU_rel": t_gpu,
        "Theta_CIM_rel": t_cim,
        "Theta_pkg": pkg,
        "note": "CPU/GPU IPC held fixed so relative throughput = s; CIM analog latency_s unchanged",
    }


def eval_s(chips0, s: float, work: Path, overrides: dict) -> dict:
    scaled = scale_chips(chips0, s)
    slot = work / f"s_{s:.4f}".replace(".", "p")
    if slot.exists():
        shutil.rmtree(slot)
    slot.mkdir(parents=True)
    t0 = time.time()
    tmax = sa.eval_hotspot(scaled, slot / "hotspot")
    write_hotspot_inputs(sa.to_rects(scaled), slot / "layout")
    run_hotspot(sa.HS_BIN, slot / "layout")
    jp = slot / "atsim.json"
    write_proxy_json(
        scaled, jp, f"freq s={s}", fill_die=True, package_overrides=overrides
    )
    st = eval_atsim_proxy(jp, verbose=False)
    sig = vm_of(st)
    th = theta_pkg(s)
    feasible = (tmax < T_LIM_C) and (sig < SIG_LIM_MPA)
    rec = {
        "s": s,
        "Tmax_C": tmax,
        "sigma_chip_max_MPa": sig,
        "sigma_chip_name": st.get("vm_normal_chip_name"),
        "w_pv_um": st.get("w_pv_um"),
        "ATSim_Tmax_C": st.get("Tmax_C"),
        "feasible": feasible,
        "fail_T": tmax >= T_LIM_C,
        "fail_sigma": sig >= SIG_LIM_MPA,
        "elapsed_s": time.time() - t0,
        "power_W": {c.name: c.power_W for c in scaled},
        "stress": {k: st[k] for k in (
            "vm_normal_chip_max_MPa", "vm_normal_chip_name", "w_pv_um", "Tmax_C", "elapsed_s"
        ) if k in st},
        **th,
    }
    (slot / "eval.json").write_text(json.dumps(rec, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return rec


def midpoint(a: float, b: float) -> float:
    return 0.5 * (a + b)


def main() -> None:
    sid = "pkg_2cpu_1gpu_2cim"
    if len(sys.argv) > 1:
        sid = sys.argv[1]
    sa_dir = ROOT / "热仿真流水线" / "runs" / f"sa_heat_vs_coupled_{sid}"
    scheme_path = sa_dir / "scheme.json"
    coupled = json.loads((sa_dir / "coupled" / "summary.json").read_text(encoding="utf-8"))
    chips0 = [
        sa.Chip(r["name"], r["kind"], r["x_m"], r["y_m"], r["w_m"], r["h_m"], r["power_W"], r["area_mm2"])
        for r in coupled["best_layout"]
    ]
    scheme = json.loads(scheme_path.read_text(encoding="utf-8")) if scheme_path.exists() else {"id": sid}
    overrides = pkg_overrides(scheme)

    out = ROOT / "热仿真流水线" / "runs" / f"freq_loop_{sid}"
    backup = ROOT / "成功运行备份" / f"频率环_{sid}_2026-10-23"
    if out.exists():
        shutil.rmtree(out)
    out.mkdir(parents=True)
    logp = out / "EXPERIMENT.md"
    logp.write_text(
        f"# Frequency loop — {sid}\n\n"
        f"- started (UTC): {datetime.now(timezone.utc).isoformat()}\n"
        f"- L* = coupled SA best (seed=22, 16 steps), layout frozen\n"
        f"- T_lim = {T_LIM_C} C, sigma_lim = {SIG_LIM_MPA} MPa (worst normal-chip interior von Mises)\n"
        f"- s_min = {S_MIN}; CPU & GPU scale together; CIM analog does not scale\n"
        f"- Theta_pkg = (Theta_CPU_rel + Theta_GPU_rel + Theta_CIM_rel)/3, refs at s=1\n"
        f"- **winner = highest Theta_pkg among feasible s only**\n"
        f"- s=0.5 still over T_lim or sigma_lim → abandon entire scheme\n"
        f"- one eval at a time (12 GB)\n\n"
        "## Run log\n\n",
        encoding="utf-8",
    )
    (out / "layout_Lstar.json").write_text(
        json.dumps(sa.snapshot(chips0), indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    log_md(logp, f"- chips: {[c.name for c in chips0]}")
    log_md(logp, f"- L* from `{sa_dir / 'coupled'}`")

    history: list[dict] = []

    def go(s: float, why: str) -> dict:
        log_md(logp, f"- eval s={s:.4f}  ({why})")
        rec = eval_s(chips0, s, out, overrides)
        history.append({**rec, "why": why})
        flag = "OK" if rec["feasible"] else (
            "FAIL_T" if rec["fail_T"] else "FAIL_SIG" if rec["fail_sigma"] else "FAIL"
        )
        log_md(
            logp,
            f"    Tmax={rec['Tmax_C']:.3f}C  sigma={rec['sigma_chip_max_MPa']:.2f} MPa "
            f"({rec['sigma_chip_name']})  Theta={rec['Theta_pkg']:.4f}  {flag}  "
            f"{rec['elapsed_s']:.1f}s"
        )
        return rec

    r_min = go(S_MIN, "lowest frequency; abandon scheme if still over window")
    abandoned = False
    if not r_min["feasible"]:
        abandoned = True
        log_md(
            logp,
            f"- **ABANDON scheme**: s={S_MIN} already "
            f"{'Tmax>=90' if r_min['fail_T'] else 'sigma>=200'}. No frequency is usable.",
        )

    winner = None
    if not abandoned:
        r1 = go(1.0, "nominal frequency (throughput reference)")
        if r1["feasible"]:
            winner = r1
            log_md(logp, "- s=1 in window → try up-clock for throughput")
            r2 = go(2.0, "first up-clock f*2")
            if r2["feasible"]:
                winner = r2
                r4 = go(S_MAX_TRY, "still in window, try s=4")
                if r4["feasible"]:
                    winner = r4
                else:
                    lo, hi = 2.0, S_MAX_TRY
                    for i in range(4):
                        mid = midpoint(lo, hi)
                        rm = go(mid, f"bisect ({lo:.3f},{hi:.3f}) after s=4 fail")
                        if rm["feasible"]:
                            lo = mid
                            if rm["Theta_pkg"] > winner["Theta_pkg"]:
                                winner = rm
                        else:
                            hi = mid
            else:
                lo, hi = 1.0, 2.0
                for i in range(4):
                    mid = midpoint(lo, hi)
                    rm = go(mid, f"bisect ({lo:.3f},{hi:.3f}) after s=2 fail")
                    if rm["feasible"]:
                        lo = mid
                        if rm["Theta_pkg"] > winner["Theta_pkg"]:
                            winner = rm
                    else:
                        hi = mid
        else:
            log_md(logp, "- s=1 out of window → only down-clock in [0.5, 1)")
            lo, hi = S_MIN, 1.0
            winner = r_min if r_min["feasible"] else None
            for i in range(4):
                mid = midpoint(lo, hi)
                rm = go(mid, f"bisect down ({lo:.3f},{hi:.3f})")
                if rm["feasible"]:
                    lo = mid
                    if winner is None or rm["Theta_pkg"] > winner["Theta_pkg"]:
                        winner = rm
                else:
                    hi = mid

    result = {
        "scheme": sid,
        "T_lim_C": T_LIM_C,
        "sigma_lim_MPa": SIG_LIM_MPA,
        "s_min": S_MIN,
        "abandoned": abandoned,
        "winner_rule": "max Theta_pkg among feasible s; thermal/stress are gates only",
        "winner": None if winner is None else {
            k: winner[k] for k in (
                "s", "Tmax_C", "sigma_chip_max_MPa", "sigma_chip_name",
                "Theta_pkg", "Theta_CPU_rel", "Theta_GPU_rel", "Theta_CIM_rel",
                "feasible", "w_pv_um",
            )
        },
        "history": history,
        "finished_utc": datetime.now(timezone.utc).isoformat(),
    }
    (out / "result.json").write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    log_md(logp, "\n## Outcome\n")
    if abandoned or winner is None:
        log_md(logp, "- no feasible frequency; scheme discarded")
    else:
        log_md(
            logp,
            f"- **winner s={winner['s']:.4f}  Theta_pkg={winner['Theta_pkg']:.4f}** "
            f"(CPU=GPU={winner['s']:.4f}, CIM=1)  "
            f"Tmax={winner['Tmax_C']:.3f}C  sigma={winner['sigma_chip_max_MPa']:.2f} MPa",
        )
        log_md(logp, f"- vs s=1: throughput ×{winner['Theta_pkg']:.3f}  (1.0 is nominal)")
    if backup.exists():
        shutil.rmtree(backup)
    shutil.copytree(out, backup)
    print(json.dumps(result["winner"], indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
