#!/usr/bin/env python3
"""Frequency loop that minimizes leakage *fraction* in the T/sigma window.

Throughput is recorded but is NOT the winner rule (ops are not T-modeled).
CPU_matmul uses McPAT from gem5 192x192 matmul (heavier than hello/64x64).
Leakage scales with Tmax after a first HotSpot pass:
    P_leak(T) = P_leak(Tref) * exp(B*(1/Tref - 1/T)), Tref=350 K (McPAT).

Winner = lowest package leak_frac among feasible s.
"""

from __future__ import annotations

import json
import math
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
TREF_K = 350.0  # McPAT XML temperature
LEAK_B = 7000.0  # K; compact subthreshold

# McPAT Total Leakage at Tref. CPU_matmul leakage same 0.78389 after 192x192.
CPU_LEAK_W = {
    "CPU_matmul": 0.78389,
    "CPU_threads": 0.78389,
    "CPU_armL2": 0.78389,
    "CPU_atomic": 0.6197,
    "CPU_starter": 0.6197,
}
# Heavier CPU task: 192x192 gem5+McPAT Runtime Dynamic.
CPU_DYN_OVERRIDE_W = {
    "CPU_matmul": 0.152807,
}
GPU_STATIC_W = {
    "GPU_3070": 22.5847,  # STATICP only (true-ish leakage); CONSTP is not leak
    "GPU_3070b": 22.5847,
    "GPU_2060S": 24.4467,
    "GPU_TITANX": 29.3748,
}
GPU_CONST_W = {
    "GPU_3070": 32.3252,
    "GPU_3070b": 32.3252,
    "GPU_2060S": 32.3252,
    "GPU_TITANX": 32.3252,
}


def log_md(path: Path, line: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8") as f:
        f.write(line.rstrip() + "\n")
    print(line, flush=True)


def leak_scale(T_C: float) -> float:
    T = T_C + 273.15
    T = max(T, 250.0)
    return math.exp(LEAK_B * (1.0 / TREF_K - 1.0 / T))


def parts(c: sa.Chip) -> dict:
    if c.kind == "CPU":
        leak = CPU_LEAK_W.get(c.name, 0.0)
        dyn = CPU_DYN_OVERRIDE_W.get(c.name)
        if dyn is None:
            dyn = max(0.0, c.power_W - leak)
        return {"dyn_s1": dyn, "leak_ref": leak, "const": 0.0}
    if c.kind == "GPU":
        st = GPU_STATIC_W.get(c.name, 0.0)
        cst = GPU_CONST_W.get(c.name, 0.0)
        dyn = max(0.0, c.power_W - st - cst)
        return {"dyn_s1": dyn, "leak_ref": st, "const": cst}
    return {"dyn_s1": 0.0, "leak_ref": 0.0, "const": c.power_W}


def apply_power(chips: list[sa.Chip], s: float, T_C: float | None) -> tuple[list[sa.Chip], dict]:
    sc = leak_scale(T_C) if T_C is not None else 1.0
    out = deepcopy(chips)
    tot_dyn = tot_leak = tot_const = 0.0
    cpu_dyn = cpu_leak = 0.0
    for c in out:
        p = parts(c)
        dyn = s * p["dyn_s1"]
        leak = p["leak_ref"] * sc
        cst = p["const"]
        c.power_W = dyn + leak + cst
        tot_dyn += dyn
        tot_leak += leak
        tot_const += cst
        if c.kind == "CPU":
            cpu_dyn += dyn
            cpu_leak += leak
    tot = tot_dyn + tot_leak + tot_const
    br = {
        "P_dyn_W": tot_dyn,
        "P_leak_W": tot_leak,
        "P_const_W": tot_const,
        "P_tot_W": tot,
        "leak_frac_pkg": tot_leak / tot if tot else None,
        "leak_frac_cpu": cpu_leak / (cpu_dyn + cpu_leak) if (cpu_dyn + cpu_leak) else None,
        "leak_T_scale": sc,
        "T_for_leak_C": T_C,
    }
    return out, br


def eval_s(chips0, s: float, work: Path, overrides: dict) -> dict:
    slot = work / f"s_{s:.4f}".replace(".", "p")
    if slot.exists():
        shutil.rmtree(slot)
    slot.mkdir(parents=True)
    t0 = time.time()
    # pass 1: leak at McPAT Tref
    scaled, _ = apply_power(chips0, s, T_C=None)
    tmax1 = sa.eval_hotspot(scaled, slot / "hotspot1")
    # pass 2: leak(T), then T/sigma for the window
    scaled, br = apply_power(chips0, s, T_C=tmax1)
    tmax = sa.eval_hotspot(scaled, slot / "hotspot")
    write_hotspot_inputs(sa.to_rects(scaled), slot / "layout")
    run_hotspot(sa.HS_BIN, slot / "layout")
    jp = slot / "atsim.json"
    write_proxy_json(scaled, jp, f"leakage s={s}", fill_die=True, package_overrides=overrides)
    st = eval_atsim_proxy(jp, verbose=False)
    sig = vm_of(st)
    feasible = (tmax < T_LIM_C) and (sig < SIG_LIM_MPA)
    rec = {
        "s": s,
        "Tmax_C": tmax,
        "Tmax_pass1_C": tmax1,
        "sigma_chip_max_MPa": sig,
        "sigma_chip_name": st.get("vm_normal_chip_name"),
        "w_pv_um": st.get("w_pv_um"),
        "feasible": feasible,
        "fail_T": tmax >= T_LIM_C,
        "fail_sigma": sig >= SIG_LIM_MPA,
        "elapsed_s": time.time() - t0,
        "power_W": {c.name: c.power_W for c in scaled},
        "Theta_pkg_record_only": (s + s + 1.0) / 3.0,
        **br,
    }
    (slot / "eval.json").write_text(json.dumps(rec, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return rec


def midpoint(a: float, b: float) -> float:
    return 0.5 * (a + b)


def better(a: dict, b: dict) -> bool:
    """True if a is a better winner than b (lower leak_frac_pkg)."""
    fa, fb = a.get("leak_frac_pkg"), b.get("leak_frac_pkg")
    if fa is None:
        return False
    if fb is None:
        return True
    return fa < fb - 1e-12


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

    out = ROOT / "热仿真流水线" / "runs" / f"leakage_loop_{sid}"
    backup = ROOT / "成功运行备份" / f"漏电占比环_{sid}_2026-10-25"
    if out.exists():
        shutil.rmtree(out)
    out.mkdir(parents=True)
    logp = out / "EXPERIMENT.md"
    logp.write_text(
        f"# Leakage-fraction loop — {sid}\n\n"
        f"- started (UTC): {datetime.now(timezone.utc).isoformat()}\n"
        f"- L* frozen from coupled SA best\n"
        f"- CPU_matmul power from gem5 192×192 matmul + McPAT: "
        f"Runtime Dynamic 0.152807 W, Leakage 0.78389 W (was 64×64 ~0.00026 W dyn)\n"
        f"- T_lim={T_LIM_C} C, sigma_lim={SIG_LIM_MPA} MPa, s_min={S_MIN}\n"
        f"- leak(T)=leak(350K)*exp({LEAK_B}*(1/350-1/T)); HotSpot twice per s\n"
        f"- GPU STATICP as leak; CONSTP as unscaled non-leak\n"
        f"- **winner = lowest leak_frac_pkg among feasible s** "
        f"(throughput recorded only, not the objective)\n"
        f"- s=0.5 still over window → abandon scheme\n\n"
        "## Run log\n\n",
        encoding="utf-8",
    )
    (out / "layout_Lstar.json").write_text(
        json.dumps(sa.snapshot(chips0), indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    log_md(logp, f"- chips: {[c.name for c in chips0]}")

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
            f"({rec['sigma_chip_name']})  leak_frac_pkg={rec['leak_frac_pkg']:.4f}  "
            f"leak_frac_cpu={rec['leak_frac_cpu']:.4f}  {flag}  {rec['elapsed_s']:.1f}s"
        )
        return rec

    r_min = go(S_MIN, "lowest frequency; abandon if still over window")
    abandoned = not r_min["feasible"]
    if abandoned:
        log_md(logp, "- **ABANDON scheme** at s_min")

    winner = None
    if not abandoned:
        r1 = go(1.0, "nominal")
        if r1["feasible"]:
            winner = r1
            log_md(logp, "- s=1 in window → search for lower leak_frac")
            r2 = go(2.0, "s=2")
            cand = [r2]
            if r2["feasible"]:
                cand.append(go(S_MAX_TRY, "s=4"))
                lo, hi = (2.0, S_MAX_TRY) if not cand[-1]["feasible"] else (S_MAX_TRY, S_MAX_TRY)
            else:
                lo, hi = 1.0, 2.0
            if lo < hi - 1e-6:
                for _ in range(4):
                    mid = midpoint(lo, hi)
                    rm = go(mid, f"bisect ({lo:.3f},{hi:.3f})")
                    cand.append(rm)
                    if rm["feasible"]:
                        lo = mid
                    else:
                        hi = mid
            for rm in cand:
                if rm["feasible"] and better(rm, winner):
                    winner = rm
        else:
            log_md(logp, "- s=1 out of window → down-clock only")
            winner = r_min
            lo, hi = S_MIN, 1.0
            for _ in range(4):
                mid = midpoint(lo, hi)
                rm = go(mid, f"bisect down ({lo:.3f},{hi:.3f})")
                if rm["feasible"]:
                    lo = mid
                    if better(rm, winner):
                        winner = rm
                else:
                    hi = mid

    result = {
        "scheme": sid,
        "objective": "min leak_frac_pkg among feasible s",
        "cpu_task": "x86 Timing+L1/L2 192x192 matmul (CPU_matmul)",
        "abandoned": abandoned,
        "winner": None if winner is None else {
            k: winner[k] for k in (
                "s", "Tmax_C", "sigma_chip_max_MPa", "sigma_chip_name",
                "leak_frac_pkg", "leak_frac_cpu", "P_dyn_W", "P_leak_W", "P_tot_W",
                "Theta_pkg_record_only",
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
            f"- **winner s={winner['s']:.4f}  leak_frac_pkg={winner['leak_frac_pkg']:.4f}  "
            f"leak_frac_cpu={winner['leak_frac_cpu']:.4f}**  "
            f"Tmax={winner['Tmax_C']:.3f}C  sigma={winner['sigma_chip_max_MPa']:.2f} MPa",
        )
        log_md(logp, "- throughput is recorded only and was not used to pick the winner")
    if backup.exists():
        shutil.rmtree(backup)
    shutil.copytree(out, backup)
    print(json.dumps(result["winner"], indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
