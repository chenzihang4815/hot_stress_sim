# Frequency loop — pkg_2cpu_1gpu_2cim

- started (UTC): 2026-10-23T10:55:35.550082+00:00
- L* = coupled SA best (seed=22, 16 steps), layout frozen
- T_lim = 90.0 C, sigma_lim = 200.0 MPa (worst normal-chip interior von Mises)
- s_min = 0.5; CPU & GPU scale together; CIM analog does not scale
- Theta_pkg = (Theta_CPU_rel + Theta_GPU_rel + Theta_CIM_rel)/3, refs at s=1
- **winner = highest Theta_pkg among feasible s only**
- s=0.5 still over T_lim or sigma_lim → abandon entire scheme
- one eval at a time (12 GB)

## Run log

- chips: ['CPU_matmul', 'CPU_threads', 'GPU_3070', 'CIM_jia', 'CIM_colonnade']
- L* from `/home/czh/任务仿真产生功耗文件工具/热仿真流水线/runs/sa_heat_vs_coupled_pkg_2cpu_1gpu_2cim/coupled`
- eval s=0.5000  (lowest frequency; abandon scheme if still over window)
    Tmax=54.550C  sigma=31.62 MPa (CIM_jia)  Theta=0.6667  OK  22.8s
- eval s=1.0000  (nominal frequency (throughput reference))
    Tmax=55.490C  sigma=34.74 MPa (CIM_jia)  Theta=1.0000  OK  20.3s
- s=1 in window → try up-clock for throughput
- eval s=2.0000  (first up-clock f*2)
    Tmax=57.370C  sigma=40.97 MPa (CIM_jia)  Theta=1.6667  OK  20.5s
- eval s=4.0000  (still in window, try s=4)
    Tmax=61.140C  sigma=53.45 MPa (CIM_jia)  Theta=3.0000  OK  22.9s

## Outcome
- **winner s=4.0000  Theta_pkg=3.0000** (CPU=GPU=4.0000, CIM=1)  Tmax=61.140C  sigma=53.45 MPa
- vs s=1: throughput ×3.000  (1.0 is nominal)
