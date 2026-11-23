# Heat-only vs heat+stress SA — pkg_1cpu_1gpu_1cim

- started (UTC): 2026-10-25T12:51:29.803252+00:00
- title: 1 CPU + 1 GPU + 1 CIM（3 颗）
- why: 最小异构：matmul + 3070 + jia
- n_chiplets: 3
- seed=22  N_STEPS=16  T0=2.0  alpha=0.95
- gaps: any 0.5 mm, GPU–GPU 2.0 mm
- coupled J = (Tmax-45)/25.0 + 1.0*vm_chip_max/100.0 + 0.5*w_pv/4.0
- vm = worst interior von Mises among normal chiplets (CPU/GPU, CIM area≥1 mm²)
- ATSim package overrides: {"pitch_die_mm": 2.0, "pitch_coarse_mm": 8.0}
- one scheme at a time (12 GB); CG only, no spsolve

## Run log

- start outline 26.92 x 21.56 mm  min_gap=1.11 mm  GPU-GPU=inf mm
- `heat_only` init  Tmax=55.470C  vm_chipmax=31.21 MPa (GPU_3070)  w_pv=0.716 um  J=0.4188  stress=17.2s
  `heat_only` step   2 BEST nudge GPU_3070 (0.5,0.0)mm  Tmax=55.460C  vm=31.21  wpv=0.716  J=0.4184
  `heat_only` step   3 BEST rotate CPU_matmul  Tmax=55.450C  vm=31.21  wpv=0.716  J=0.4180
  `heat_only` step   4 accept rotate GPU_3070  Tmax=55.490C  J=0.4196
  `heat_only` step   5 accept nudge CIM_jia (1.0,0.0)mm  Tmax=55.490C  J=0.4196
  `heat_only` step  10 accept nudge GPU_3070 (1.0,0.0)mm  Tmax=55.470C  J=0.4188
  `heat_only` step  11 BEST rotate GPU_3070  Tmax=55.440C  vm=31.21  wpv=0.716  J=0.4176
  `heat_only` step  12 BEST jump GPU_3070  Tmax=55.410C  vm=31.21  wpv=0.716  J=0.4164
  `heat_only` step  15 accept nudge GPU_3070 (-0.5,0.0)mm  Tmax=55.410C  J=0.4164
- `heat_only` rescoring best with silicon-fill stress …
- `heat_only` done  Tmax 55.470→55.410C  vm 31.21→32.49 MPa  wpv 0.716→0.877 um  accepted 8/16
- `coupled` init  Tmax=55.470C  vm_chipmax=31.21 MPa (GPU_3070)  w_pv=0.716 um  J=0.8205  stress=11.6s
  `coupled` step   2 accept nudge GPU_3070 (0.5,0.0)mm  Tmax=55.460C  J=0.8237
  `coupled` step   3 BEST rotate CPU_matmul  Tmax=55.450C  vm=31.00  wpv=0.694  J=0.8148
  `coupled` step   4 BEST rotate GPU_3070  Tmax=55.490C  vm=30.37  wpv=0.635  J=0.8027
  `coupled` step   5 BEST nudge CIM_jia (1.0,0.0)mm  Tmax=55.490C  vm=30.39  wpv=0.599  J=0.7984
  `coupled` step  10 accept nudge GPU_3070 (1.0,0.0)mm  Tmax=55.470C  J=0.8058
  `coupled` step  11 accept rotate GPU_3070  Tmax=55.440C  J=0.8160
  `coupled` step  12 accept jump GPU_3070  Tmax=55.410C  J=0.8509
  `coupled` step  15 accept nudge GPU_3070 (-0.5,0.0)mm  Tmax=55.410C  J=0.8476
- `coupled` done  Tmax 55.470→55.490C  vm 31.21→30.39 MPa  wpv 0.716→0.599 um  accepted 8/16

## Table
| arm | Tmax init→best | vm_chip_max init→best | w_pv init→best | accepted |
|---|---|---|---|---|
| heat_only | 55.470→55.410 | 31.21→32.49 (GPU_3070) | 0.716→0.877 | 8/16 |
| coupled | 55.470→55.490 | 31.21→30.39 (GPU_3070) | 0.716→0.599 | 8/16 |

coupled − heat_only at best: dTmax=0.080 C, d_vm=-2.10 MPa, d_wpv=-0.278 um

wall 157.7s  finished 2026-10-25T12:54:07.486968+00:00
