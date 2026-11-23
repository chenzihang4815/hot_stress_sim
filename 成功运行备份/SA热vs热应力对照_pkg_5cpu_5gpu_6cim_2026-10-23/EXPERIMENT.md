# Heat-only vs heat+stress SA — pkg_5cpu_5gpu_6cim

- started (UTC): 2026-10-25T14:01:14.132477+00:00
- title: 5 CPU + 5 GPU + 6 CIM（16 颗）
- why: 15 芯再加第二颗 jia
- n_chiplets: 16
- seed=22  N_STEPS=16  T0=2.0  alpha=0.95
- gaps: any 0.5 mm, GPU–GPU 2.0 mm
- coupled J = (Tmax-45)/25.0 + 1.0*vm_chip_max/100.0 + 0.5*w_pv/4.0
- vm = worst interior von Mises among normal chiplets (CPU/GPU, CIM area≥1 mm²)
- ATSim package overrides: {"pitch_die_mm": 2.5, "pitch_coarse_mm": 10.0, "nz_die": 2, "nz_spreader": 2, "nz_sink": 3, "nz_tim": 1, "nz_substrate": 2, "nz_interposer": 2}
- one scheme at a time (12 GB); CG only, no spsolve

## Run log

- start outline 47.10 x 86.51 mm  min_gap=0.52 mm  GPU-GPU=2.36 mm
- `heat_only` init  Tmax=87.870C  vm_chipmax=256.48 MPa (GPU_TITANX)  w_pv=7.308 um  J=1.7148  stress=133.4s
  `heat_only` step   1 accept rotate CPU_starter  Tmax=87.870C  J=1.7148
  `heat_only` step   3 accept nudge GPU_3070 (-1.0,0.0)mm  Tmax=87.970C  J=1.7188
  `heat_only` step   4 accept swap CIM_basic<->CIM_colonnade  Tmax=87.970C  J=1.7188
  `heat_only` step   8 accept rotate CPU_starter  Tmax=87.970C  J=1.7188
  `heat_only` step   9 accept nudge CIM_basic (0.0,0.5)mm  Tmax=87.970C  J=1.7188
  `heat_only` step  14 accept nudge GPU_3070 (0.0,-1.0)mm  Tmax=88.050C  J=1.7220
  `heat_only` step  16 accept nudge GPU_3070b (-0.5,0.0)mm  Tmax=87.910C  J=1.7164
- `heat_only` rescoring best with silicon-fill stress …
- `heat_only` done  Tmax 87.870→87.870C  vm 256.48→256.48 MPa  wpv 7.308→7.308 um  accepted 7/16
- `coupled` init  Tmax=87.870C  vm_chipmax=256.48 MPa (GPU_TITANX)  w_pv=7.308 um  J=5.1931  stress=140.7s
  `coupled` step   1 accept rotate CPU_starter  Tmax=87.870C  J=5.1931
  `coupled` step   3 accept nudge GPU_3070 (-1.0,0.0)mm  Tmax=87.970C  J=5.1942
  `coupled` step   4 accept swap CIM_basic<->CIM_colonnade  Tmax=87.970C  J=5.1943
  `coupled` step   8 accept rotate CPU_starter  Tmax=87.970C  J=5.1943
  `coupled` step   9 accept nudge CIM_basic (0.0,0.5)mm  Tmax=87.970C  J=5.1943
  `coupled` step  14 accept nudge GPU_3070 (0.0,-1.0)mm  Tmax=88.050C  J=5.2020
  `coupled` step  16 accept nudge GPU_3070b (-0.5,0.0)mm  Tmax=87.910C  J=5.1980
- `coupled` done  Tmax 87.870→87.870C  vm 256.48→256.48 MPa  wpv 7.308→7.308 um  accepted 7/16

## Table
| arm | Tmax init→best | vm_chip_max init→best | w_pv init→best | accepted |
|---|---|---|---|---|
| heat_only | 87.870→87.870 | 256.48→256.48 (GPU_TITANX) | 7.308→7.308 | 7/16 |
| coupled | 87.870→87.870 | 256.48→256.48 (GPU_TITANX) | 7.308→7.308 | 7/16 |

coupled − heat_only at best: dTmax=0.000 C, d_vm=0.00 MPa, d_wpv=0.000 um

wall 1560.5s  finished 2026-10-25T14:27:14.659856+00:00
