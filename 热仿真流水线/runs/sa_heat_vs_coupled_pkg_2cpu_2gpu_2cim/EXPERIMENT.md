# Heat-only vs heat+stress SA — pkg_2cpu_2gpu_2cim

- started (UTC): 2026-10-23T08:13:13.816231+00:00
- title: 2 CPU + 2 GPU + 2 CIM（6 颗）
- why: 两颗不同 GPU，看间距是否被应力拉开
- n_chiplets: 6
- seed=22  N_STEPS=16  T0=2.0  alpha=0.95
- gaps: any 0.5 mm, GPU–GPU 2.0 mm
- coupled J = (Tmax-45)/25.0 + 1.0*vm_chip_max/100.0 + 0.5*w_pv/4.0
- vm = worst interior von Mises among normal chiplets (CPU/GPU, CIM area≥1 mm²)
- ATSim package overrides: {"pitch_die_mm": 2.0, "pitch_coarse_mm": 8.0}
- one scheme at a time (12 GB); CG only, no spsolve

## Run log

- start outline 27.94 x 47.16 mm  min_gap=0.92 mm  GPU-GPU=2.05 mm
- `heat_only` init  Tmax=63.650C  vm_chipmax=74.33 MPa (GPU_2060S)  w_pv=2.418 um  J=0.7460  stress=32.0s
  `heat_only` step   5 accept rotate CIM_wang  Tmax=63.650C  J=0.7460
  `heat_only` step   8 accept jump+rot GPU_2060S  Tmax=63.750C  J=0.7500
  `heat_only` step  10 accept nudge GPU_3070 (0.0,1.0)mm  Tmax=63.760C  J=0.7504
  `heat_only` step  12 accept nudge CIM_jia (0.0,-1.5)mm  Tmax=63.760C  J=0.7504
  `heat_only` step  13 BEST rotate GPU_2060S  Tmax=63.570C  vm=74.33  wpv=2.418  J=0.7428
  `heat_only` step  14 accept rotate GPU_3070  Tmax=63.910C  J=0.7564
- `heat_only` rescoring best with silicon-fill stress …
- `heat_only` done  Tmax 63.650→63.570C  vm 74.33→76.75 MPa  wpv 2.418→2.334 um  accepted 6/16
- `coupled` init  Tmax=63.650C  vm_chipmax=74.33 MPa (GPU_2060S)  w_pv=2.418 um  J=1.7916  stress=31.0s
  `coupled` step   5 BEST rotate CIM_wang  Tmax=63.650C  vm=74.33  wpv=2.402  J=1.7896
  `coupled` step   8 accept jump+rot GPU_2060S  Tmax=63.750C  J=1.8374
  `coupled` step  10 accept nudge GPU_3070 (0.0,1.0)mm  Tmax=63.760C  J=1.8432
  `coupled` step  12 accept nudge CIM_jia (0.0,-1.5)mm  Tmax=63.760C  J=1.8432
  `coupled` step  13 accept rotate GPU_2060S  Tmax=63.570C  J=1.8020
  `coupled` step  14 accept rotate GPU_3070  Tmax=63.910C  J=1.8043
- `coupled` done  Tmax 63.650→63.650C  vm 74.33→74.33 MPa  wpv 2.418→2.402 um  accepted 6/16

## Table
| arm | Tmax init→best | vm_chip_max init→best | w_pv init→best | accepted |
|---|---|---|---|---|
| heat_only | 63.650→63.570 | 74.33→76.75 (GPU_2060S) | 2.418→2.334 | 6/16 |
| coupled | 63.650→63.650 | 74.33→74.33 (GPU_2060S) | 2.418→2.402 | 6/16 |

coupled − heat_only at best: dTmax=0.080 C, d_vm=-2.42 MPa, d_wpv=0.068 um

wall 333.2s  finished 2026-10-23T08:18:46.944454+00:00
