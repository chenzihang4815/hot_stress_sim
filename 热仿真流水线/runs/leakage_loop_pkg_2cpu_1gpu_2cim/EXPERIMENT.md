# Leakage-fraction loop — pkg_2cpu_1gpu_2cim

- started (UTC): 2026-10-25T04:18:29.176782+00:00
- L* frozen from coupled SA best
- CPU_matmul power from gem5 192×192 matmul + McPAT: Runtime Dynamic 0.152807 W, Leakage 0.78389 W (was 64×64 ~0.00026 W dyn)
- T_lim=90.0 C, sigma_lim=200.0 MPa, s_min=0.5
- leak(T)=leak(350K)*exp(7000.0*(1/350-1/T)); HotSpot twice per s
- GPU STATICP as leak; CONSTP as unscaled non-leak
- **winner = lowest leak_frac_pkg among feasible s** (throughput recorded only, not the objective)
- s=0.5 still over window → abandon scheme

## Run log

- chips: ['CPU_matmul', 'CPU_threads', 'GPU_3070', 'CIM_jia', 'CIM_colonnade']
- eval s=0.5000  (lowest frequency; abandon if still over window)
    Tmax=51.840C  sigma=22.58 MPa (CIM_jia)  leak_frac_pkg=0.1385  leak_frac_cpu=0.7121  OK  23.8s
- eval s=1.0000  (nominal)
    Tmax=52.830C  sigma=25.88 MPa (CIM_jia)  leak_frac_pkg=0.1284  leak_frac_cpu=0.5678  OK  21.8s
- s=1 in window → search for lower leak_frac
- eval s=2.0000  (s=2)
    Tmax=54.830C  sigma=32.51 MPa (CIM_jia)  leak_frac_pkg=0.1153  leak_frac_cpu=0.4258  OK  22.8s
- eval s=4.0000  (s=4)
    Tmax=58.880C  sigma=45.92 MPa (CIM_jia)  leak_frac_pkg=0.1035  leak_frac_cpu=0.3198  OK  22.9s

## Outcome
- **winner s=4.0000  leak_frac_pkg=0.1035  leak_frac_cpu=0.3198**  Tmax=58.880C  sigma=45.92 MPa
- throughput is recorded only and was not used to pick the winner
