# Heat-only vs heat+stress SA — pkg_5cpu_3gpu_5cim

- started (UTC): 2026-10-23T08:38:34.221979+00:00
- title: 5 CPU + 3 GPU + 5 CIM（13 颗）
- why: 已有大封装；对照实验的基准规模
- n_chiplets: 13
- seed=22  N_STEPS=16  T0=2.0  alpha=0.95
- gaps: any 0.5 mm, GPU–GPU 2.0 mm
- coupled J = (Tmax-45)/25.0 + 1.0*vm_chip_max/100.0 + 0.5*w_pv/4.0
- vm = worst interior von Mises among normal chiplets (CPU/GPU, CIM area≥1 mm²)
- ATSim package overrides: {"pitch_die_mm": 2.0, "pitch_coarse_mm": 8.0}
- one scheme at a time (12 GB); CG only, no spsolve

## Run log

- start outline 57.12 x 40.55 mm  min_gap=0.54 mm  GPU-GPU=2.62 mm
- `heat_only` init  Tmax=70.990C  vm_chipmax=150.16 MPa (GPU_2060S)  w_pv=4.457 um  J=1.0396  stress=163.4s
  `heat_only` step   1 accept rotate CPU_armL2  Tmax=70.990C  J=1.0396
  `heat_only` step   4 accept rotate CIM_wang  Tmax=70.990C  J=1.0396
  `heat_only` step   6 BEST jump GPU_3070  Tmax=70.040C  vm=150.16  wpv=4.457  J=1.0016
  `heat_only` step   8 accept nudge GPU_2060S (-1.0,0.0)mm  Tmax=70.620C  J=1.0248
  `heat_only` step   9 accept nudge CIM_wang (-1.5,0.0)mm  Tmax=70.630C  J=1.0252
  `heat_only` step  10 accept nudge GPU_2060S (0.0,0.5)mm  Tmax=70.330C  J=1.0132
  `heat_only` step  11 accept nudge CIM_jia (0.0,-1.0)mm  Tmax=70.330C  J=1.0132
  `heat_only` step  12 accept rotate CPU_armL2  Tmax=70.330C  J=1.0132
  `heat_only` step  15 accept nudge GPU_3070 (0.0,-0.5)mm  Tmax=70.780C  J=1.0312
  `heat_only` step  16 accept nudge CPU_atomic (-1.5,0.0)mm  Tmax=70.730C  J=1.0292
- `heat_only` rescoring best with silicon-fill stress …
- `heat_only` done  Tmax 70.990→70.040C  vm 150.16→152.40 MPa  wpv 4.457→4.479 um  accepted 10/16
- `coupled` init  Tmax=70.990C  vm_chipmax=150.16 MPa (GPU_2060S)  w_pv=4.457 um  J=3.0983  stress=148.1s
  `coupled` step   1 BEST rotate CPU_armL2  Tmax=70.990C  vm=150.16  wpv=4.456  J=3.0982
  `coupled` step   4 accept rotate CIM_wang  Tmax=70.990C  J=3.0982
  `coupled` step   6 BEST jump GPU_3070  Tmax=70.040C  vm=152.40  wpv=4.479  J=3.0855
  `coupled` step   8 accept nudge GPU_2060S (-1.0,0.0)mm  Tmax=70.620C  J=3.1245
  `coupled` step   9 accept nudge CIM_wang (-1.5,0.0)mm  Tmax=70.630C  J=3.1257
  `coupled` step  10 accept nudge GPU_2060S (0.0,0.5)mm  Tmax=70.330C  J=3.1124
  `coupled` step  11 accept nudge CIM_jia (0.0,-1.0)mm  Tmax=70.330C  J=3.1124
  `coupled` step  12 accept rotate CPU_armL2  Tmax=70.330C  J=3.1124
  `coupled` step  15 accept nudge GPU_3070 (0.0,-0.5)mm  Tmax=70.780C  J=3.1250
  `coupled` step  16 accept nudge CPU_atomic (-1.5,0.0)mm  Tmax=70.730C  J=3.0961
- `coupled` done  Tmax 70.990→70.040C  vm 150.16→152.40 MPa  wpv 4.457→4.479 um  accepted 10/16

## Table
| arm | Tmax init→best | vm_chip_max init→best | w_pv init→best | accepted |
|---|---|---|---|---|
| heat_only | 70.990→70.040 | 150.16→152.40 (GPU_2060S) | 4.457→4.479 | 10/16 |
| coupled | 70.990→70.040 | 150.16→152.40 (GPU_2060S) | 4.457→4.479 | 10/16 |

coupled − heat_only at best: dTmax=0.000 C, d_vm=0.00 MPa, d_wpv=0.000 um

wall 2020.8s  finished 2026-10-23T09:12:15.052208+00:00
