# Heat-only vs heat+stress SA — pkg_2cpu_1gpu_1cim

- started (UTC): 2026-10-25T12:54:09.097200+00:00
- title: 2 CPU + 1 GPU + 1 CIM（4 颗）
- why: 两颗 x86 CPU + 一颗 GPU
- n_chiplets: 4
- seed=22  N_STEPS=16  T0=2.0  alpha=0.95
- gaps: any 0.5 mm, GPU–GPU 2.0 mm
- coupled J = (Tmax-45)/25.0 + 1.0*vm_chip_max/100.0 + 0.5*w_pv/4.0
- vm = worst interior von Mises among normal chiplets (CPU/GPU, CIM area≥1 mm²)
- ATSim package overrides: {"pitch_die_mm": 2.0, "pitch_coarse_mm": 8.0}
- one scheme at a time (12 GB); CG only, no spsolve

## Run log

- start outline 21.50 x 26.87 mm  min_gap=0.81 mm  GPU-GPU=inf mm
- `heat_only` init  Tmax=55.590C  vm_chipmax=32.15 MPa (GPU_3070)  w_pv=0.738 um  J=0.4236  stress=13.9s
  `heat_only` step   2 accept nudge CIM_jia (0.5,0.0)mm  Tmax=55.590C  J=0.4236
  `heat_only` step   3 BEST rotate CPU_matmul  Tmax=55.580C  vm=32.15  wpv=0.738  J=0.4232
  `heat_only` step   4 BEST rotate GPU_3070  Tmax=55.500C  vm=32.15  wpv=0.738  J=0.4200
  `heat_only` step   8 accept nudge GPU_3070 (0.0,1.0)mm  Tmax=55.520C  J=0.4208
  `heat_only` step  10 accept nudge GPU_3070 (1.0,0.0)mm  Tmax=55.520C  J=0.4208
  `heat_only` step  11 accept rotate CIM_jia  Tmax=55.520C  J=0.4208
  `heat_only` step  13 BEST nudge GPU_3070 (0.0,-1.5)mm  Tmax=55.490C  vm=32.15  wpv=0.738  J=0.4196
  `heat_only` step  14 accept rotate CIM_jia  Tmax=55.490C  J=0.4196
  `heat_only` step  15 accept rotate GPU_3070  Tmax=55.560C  J=0.4224
  `heat_only` step  16 accept swap CPU_matmul<->CIM_jia  Tmax=55.570C  J=0.4228
- `heat_only` rescoring best with silicon-fill stress …
- `heat_only` done  Tmax 55.590→55.490C  vm 32.15→35.63 MPa  wpv 0.738→0.970 um  accepted 10/16
- `coupled` init  Tmax=55.590C  vm_chipmax=32.15 MPa (GPU_3070)  w_pv=0.738 um  J=0.8373  stress=13.3s
  `coupled` step   2 accept nudge CIM_jia (0.5,0.0)mm  Tmax=55.590C  J=0.8468
  `coupled` step   3 accept rotate CPU_matmul  Tmax=55.580C  J=0.8527
  `coupled` step   4 accept rotate GPU_3070  Tmax=55.500C  J=0.8812
  `coupled` step   8 accept nudge GPU_3070 (0.0,1.0)mm  Tmax=55.520C  J=0.8789
  `coupled` step  10 accept nudge GPU_3070 (1.0,0.0)mm  Tmax=55.520C  J=0.8840
  `coupled` step  11 accept rotate CIM_jia  Tmax=55.520C  J=0.8845
  `coupled` step  13 accept nudge GPU_3070 (0.0,-1.5)mm  Tmax=55.490C  J=0.8971
  `coupled` step  14 accept rotate CIM_jia  Tmax=55.490C  J=0.8972
  `coupled` step  15 accept rotate GPU_3070  Tmax=55.560C  J=0.8731
  `coupled` step  16 accept swap CPU_matmul<->CIM_jia  Tmax=55.570C  J=0.8617
- `coupled` done  Tmax 55.590→55.590C  vm 32.15→32.15 MPa  wpv 0.738→0.738 um  accepted 10/16

## Table
| arm | Tmax init→best | vm_chip_max init→best | w_pv init→best | accepted |
|---|---|---|---|---|
| heat_only | 55.590→55.490 | 32.15→35.63 (GPU_3070) | 0.738→0.970 | 10/16 |
| coupled | 55.590→55.590 | 32.15→32.15 (GPU_3070) | 0.738→0.738 | 10/16 |

coupled − heat_only at best: dTmax=0.100 C, d_vm=-3.48 MPa, d_wpv=-0.232 um

wall 238.1s  finished 2026-10-25T12:58:07.246161+00:00
