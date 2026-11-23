# Heat-only vs heat+stress SA — pkg_5cpu_4gpu_5cim

- started (UTC): 2026-10-23T09:12:16.684083+00:00
- title: 5 CPU + 4 GPU + 5 CIM（14 颗）
- why: 多一颗同规格 3070，看同类 GPU 间距是否被应力拉开
- n_chiplets: 14
- seed=22  N_STEPS=16  T0=2.0  alpha=0.95
- gaps: any 0.5 mm, GPU–GPU 2.0 mm
- coupled J = (Tmax-45)/25.0 + 1.0*vm_chip_max/100.0 + 0.5*w_pv/4.0
- vm = worst interior von Mises among normal chiplets (CPU/GPU, CIM area≥1 mm²)
- ATSim package overrides: {"pitch_die_mm": 2.5, "pitch_coarse_mm": 10.0, "nz_die": 2, "nz_spreader": 2, "nz_sink": 3, "nz_tim": 1, "nz_substrate": 2, "nz_interposer": 2}
- one scheme at a time (12 GB); CG only, no spsolve

## Run log

- start outline 49.17 x 51.21 mm  min_gap=0.57 mm  GPU-GPU=2.25 mm
- `heat_only` init  Tmax=76.290C  vm_chipmax=189.78 MPa (GPU_3070b)  w_pv=5.499 um  J=1.2516  stress=106.0s
  `heat_only` step   2 accept nudge GPU_TITANX (0.5,0.0)mm  Tmax=76.310C  J=1.2524
  `heat_only` step   3 accept rotate CIM_basic  Tmax=76.310C  J=1.2524
  `heat_only` step   5 accept rotate CIM_colonnade  Tmax=76.310C  J=1.2524
  `heat_only` step   6 accept jump CIM_jia  Tmax=77.810C  J=1.3124
  `heat_only` step   7 accept nudge CIM_jia (1.5,0.0)mm  Tmax=77.380C  J=1.2952
  `heat_only` step  10 accept nudge GPU_2060S (0.0,1.0)mm  Tmax=77.220C  J=1.2888
  `heat_only` step  13 accept nudge GPU_3070 (0.0,-1.0)mm  Tmax=77.240C  J=1.2896
  `heat_only` step  15 accept nudge GPU_3070b (-0.5,0.0)mm  Tmax=77.300C  J=1.2920
- `heat_only` rescoring best with silicon-fill stress …
- `heat_only` done  Tmax 76.290→76.290C  vm 189.78→189.78 MPa  wpv 5.499→5.499 um  accepted 8/16
- `coupled` init  Tmax=76.290C  vm_chipmax=189.78 MPa (GPU_3070b)  w_pv=5.499 um  J=3.8368  stress=104.1s
  `coupled` step   2 accept nudge GPU_TITANX (0.5,0.0)mm  Tmax=76.310C  J=3.8436
  `coupled` step   3 accept rotate CIM_basic  Tmax=76.310C  J=3.8436
  `coupled` step   5 BEST rotate CIM_colonnade  Tmax=76.310C  vm=189.07  wpv=5.543  J=3.8360
  `coupled` step   6 accept jump CIM_jia  Tmax=77.810C  J=4.0598
  `coupled` step   7 accept nudge CIM_jia (1.5,0.0)mm  Tmax=77.380C  J=4.0682
  `coupled` step  10 accept nudge GPU_2060S (0.0,1.0)mm  Tmax=77.220C  J=4.0557
  `coupled` step  13 accept nudge GPU_3070 (0.0,-1.0)mm  Tmax=77.240C  J=4.0585
  `coupled` step  15 accept nudge GPU_3070b (-0.5,0.0)mm  Tmax=77.300C  J=4.0813
- `coupled` done  Tmax 76.290→76.310C  vm 189.78→189.07 MPa  wpv 5.499→5.543 um  accepted 8/16

## Table
| arm | Tmax init→best | vm_chip_max init→best | w_pv init→best | accepted |
|---|---|---|---|---|
| heat_only | 76.290→76.290 | 189.78→189.78 (GPU_3070b) | 5.499→5.499 | 8/16 |
| coupled | 76.290→76.310 | 189.78→189.07 (GPU_3070b) | 5.499→5.543 | 8/16 |

coupled − heat_only at best: dTmax=0.020 C, d_vm=-0.71 MPa, d_wpv=0.044 um

wall 1158.2s  finished 2026-10-23T09:31:34.831230+00:00
