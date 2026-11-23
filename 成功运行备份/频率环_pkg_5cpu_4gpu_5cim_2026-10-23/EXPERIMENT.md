# Frequency loop — pkg_5cpu_4gpu_5cim

- started (UTC): 2026-10-24T11:30:13.625279+00:00
- L* = coupled SA best (seed=22, 16 steps), layout frozen
- T_lim = 90.0 C, sigma_lim = 200.0 MPa (worst normal-chip interior von Mises)
- s_min = 0.5; CPU & GPU scale together; CIM analog does not scale
- Theta_pkg = (Theta_CPU_rel + Theta_GPU_rel + Theta_CIM_rel)/3, refs at s=1
- **winner = highest Theta_pkg among feasible s only**
- s=0.5 still over T_lim or sigma_lim → abandon entire scheme
- one eval at a time (12 GB)

## Run log

- chips: ['CPU_matmul', 'CPU_threads', 'CPU_armL2', 'CPU_atomic', 'CPU_starter', 'GPU_3070', 'GPU_2060S', 'GPU_TITANX', 'GPU_3070b', 'CIM_jia', 'CIM_colonnade', 'CIM_wang', 'CIM_basic', 'CIM_sinangil']
- L* from `/home/czh/任务仿真产生功耗文件工具/热仿真流水线/runs/sa_heat_vs_coupled_pkg_5cpu_4gpu_5cim/coupled`
- eval s=0.5000  (lowest frequency; abandon scheme if still over window)
    Tmax=73.610C  sigma=173.31 MPa (GPU_3070b)  Theta=0.6667  OK  108.7s
- eval s=1.0000  (nominal frequency (throughput reference))
    Tmax=76.310C  sigma=189.07 MPa (GPU_3070b)  Theta=1.0000  OK  107.3s
- s=1 in window → try up-clock for throughput
- eval s=2.0000  (first up-clock f*2)
    Tmax=82.510C  sigma=220.59 MPa (GPU_3070b)  Theta=1.6667  FAIL_SIG  106.4s
- eval s=1.5000  (bisect (1.000,2.000) after s=2 fail)
    Tmax=79.410C  sigma=204.83 MPa (GPU_3070b)  Theta=1.3333  FAIL_SIG  106.8s
- eval s=1.2500  (bisect (1.000,1.500) after s=2 fail)
    Tmax=77.860C  sigma=196.95 MPa (GPU_3070b)  Theta=1.1667  OK  105.8s
- eval s=1.3750  (bisect (1.250,1.500) after s=2 fail)
    Tmax=78.630C  sigma=200.89 MPa (GPU_3070b)  Theta=1.2500  FAIL_SIG  105.4s
- eval s=1.3125  (bisect (1.250,1.375) after s=2 fail)
    Tmax=78.240C  sigma=198.92 MPa (GPU_3070b)  Theta=1.2083  OK  105.7s

## Outcome
- **winner s=1.3125  Theta_pkg=1.2083** (CPU=GPU=1.3125, CIM=1)  Tmax=78.240C  sigma=198.92 MPa
- vs s=1: throughput ×1.208  (1.0 is nominal)
