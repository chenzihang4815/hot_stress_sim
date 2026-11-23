# Heat-only vs heat+stress SA — pkg_3cpu_1gpu_3cim

- started (UTC): 2026-10-23T08:18:48.505690+00:00
- title: 3 CPU + 1 GPU + 3 CIM（7 颗）
- why: 三颗 CPU + TitanX + 三 CIM
- n_chiplets: 7
- seed=22  N_STEPS=16  T0=2.0  alpha=0.95
- gaps: any 0.5 mm, GPU–GPU 2.0 mm
- coupled J = (Tmax-45)/25.0 + 1.0*vm_chip_max/100.0 + 0.5*w_pv/4.0
- vm = worst interior von Mises among normal chiplets (CPU/GPU, CIM area≥1 mm²)
- ATSim package overrides: {"pitch_die_mm": 2.0, "pitch_coarse_mm": 8.0}
- one scheme at a time (12 GB); CG only, no spsolve

## Run log

- start outline 25.78 x 34.25 mm  min_gap=0.59 mm  GPU-GPU=inf mm
- `heat_only` init  Tmax=55.300C  vm_chipmax=36.68 MPa (CIM_jia)  w_pv=1.086 um  J=0.4120  stress=38.2s
  `heat_only` step   1 accept rotate CPU_threads  Tmax=55.300C  J=0.4120
  `heat_only` step   2 accept jump GPU_TITANX  Tmax=55.310C  J=0.4124
  `heat_only` step   3 accept rotate CPU_starter  Tmax=55.310C  J=0.4124
  `heat_only` step   7 accept rotate CIM_wang  Tmax=55.310C  J=0.4124
  `heat_only` step   8 accept rotate CPU_matmul  Tmax=55.320C  J=0.4128
  `heat_only` step   9 accept jump GPU_TITANX  Tmax=55.320C  J=0.4128
  `heat_only` step  10 accept swap CPU_threads<->CPU_starter  Tmax=55.320C  J=0.4128
  `heat_only` step  12 accept nudge GPU_TITANX (0.0,1.0)mm  Tmax=55.320C  J=0.4128
  `heat_only` step  16 accept nudge GPU_TITANX (0.0,-0.5)mm  Tmax=55.340C  J=0.4136
- `heat_only` rescoring best with silicon-fill stress …
- `heat_only` done  Tmax 55.300→55.300C  vm 36.68→36.68 MPa  wpv 1.086→1.086 um  accepted 9/16
- `coupled` init  Tmax=55.300C  vm_chipmax=36.68 MPa (CIM_jia)  w_pv=1.086 um  J=0.9146  stress=37.5s
  `coupled` step   1 BEST rotate CPU_threads  Tmax=55.300C  vm=36.20  wpv=1.041  J=0.9042
  `coupled` step   2 accept jump GPU_TITANX  Tmax=55.310C  J=0.9517
  `coupled` step   3 accept rotate CPU_starter  Tmax=55.310C  J=0.9512
  `coupled` step   7 accept rotate CIM_wang  Tmax=55.310C  J=0.9498
  `coupled` step   8 accept rotate CPU_matmul  Tmax=55.320C  J=0.9503
  `coupled` step   9 accept jump GPU_TITANX  Tmax=55.320C  J=0.9514
  `coupled` step  10 accept swap CPU_threads<->CPU_starter  Tmax=55.320C  J=0.9466
  `coupled` step  12 accept nudge GPU_TITANX (0.0,1.0)mm  Tmax=55.320C  J=0.9558
  `coupled` step  16 accept nudge GPU_TITANX (0.0,-0.5)mm  Tmax=55.340C  J=0.9503
- `coupled` done  Tmax 55.300→55.300C  vm 36.68→36.20 MPa  wpv 1.086→1.041 um  accepted 9/16

## Table
| arm | Tmax init→best | vm_chip_max init→best | w_pv init→best | accepted |
|---|---|---|---|---|
| heat_only | 55.300→55.300 | 36.68→36.68 (CIM_jia) | 1.086→1.086 | 9/16 |
| coupled | 55.300→55.300 | 36.68→36.20 (CIM_jia) | 1.086→1.041 | 9/16 |

coupled − heat_only at best: dTmax=0.000 C, d_vm=-0.48 MPa, d_wpv=-0.045 um

wall 434.7s  finished 2026-10-23T08:26:03.250249+00:00
