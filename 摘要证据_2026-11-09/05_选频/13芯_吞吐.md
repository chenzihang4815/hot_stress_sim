# Frequency loop — pkg_5cpu_3gpu_5cim

- started (UTC): 2026-10-24T11:05:34.722772+00:00
- L* = coupled SA best (seed=22, 16 steps), layout frozen
- T_lim = 90.0 C, sigma_lim = 200.0 MPa (worst normal-chip interior von Mises)
- s_min = 0.5; CPU & GPU scale together; CIM analog does not scale
- Theta_pkg = (Theta_CPU_rel + Theta_GPU_rel + Theta_CIM_rel)/3, refs at s=1
- **winner = highest Theta_pkg among feasible s only**
- s=0.5 still over T_lim or sigma_lim → abandon entire scheme
- one eval at a time (12 GB)

## Run log

- chips: ['CPU_matmul', 'CPU_threads', 'CPU_armL2', 'CPU_atomic', 'CPU_starter', 'GPU_3070', 'GPU_2060S', 'GPU_TITANX', 'CIM_jia', 'CIM_colonnade', 'CIM_wang', 'CIM_basic', 'CIM_sinangil']
- L* from `/home/czh/任务仿真产生功耗文件工具/热仿真流水线/runs/sa_heat_vs_coupled_pkg_5cpu_3gpu_5cim/coupled`
- eval s=0.5000  (lowest frequency; abandon scheme if still over window)
    Tmax=68.060C  sigma=139.95 MPa (GPU_2060S)  Theta=0.6667  OK  229.7s
- eval s=1.0000  (nominal frequency (throughput reference))
    Tmax=70.040C  sigma=152.40 MPa (GPU_2060S)  Theta=1.0000  OK  151.9s
- s=1 in window → try up-clock for throughput
- eval s=2.0000  (first up-clock f*2)
    Tmax=74.860C  sigma=177.30 MPa (GPU_2060S)  Theta=1.6667  OK  150.7s
- eval s=4.0000  (still in window, try s=4)
    Tmax=84.520C  sigma=227.09 MPa (GPU_2060S)  Theta=3.0000  FAIL_SIG  148.8s
- eval s=3.0000  (bisect (2.000,4.000) after s=4 fail)
    Tmax=79.690C  sigma=202.20 MPa (GPU_2060S)  Theta=2.3333  FAIL_SIG  149.4s
- eval s=2.5000  (bisect (2.000,3.000) after s=4 fail)
    Tmax=77.280C  sigma=189.75 MPa (GPU_2060S)  Theta=2.0000  OK  148.3s
- eval s=2.7500  (bisect (2.500,3.000) after s=4 fail)
    Tmax=78.480C  sigma=195.97 MPa (GPU_2060S)  Theta=2.1667  OK  147.9s
- eval s=2.8750  (bisect (2.750,3.000) after s=4 fail)
    Tmax=79.090C  sigma=199.08 MPa (GPU_2060S)  Theta=2.2500  OK  148.1s

## Outcome
- **winner s=2.8750  Theta_pkg=2.2500** (CPU=GPU=2.8750, CIM=1)  Tmax=79.090C  sigma=199.08 MPa
- vs s=1: throughput ×2.250  (1.0 is nominal)
