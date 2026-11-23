# Heat-only vs heat+stress SA — pkg_4cpu_3gpu_5cim

- started (UTC): 2026-10-25T13:18:41.494649+00:00
- title: 4 CPU + 3 GPU + 5 CIM（12 颗）
- why: 11 芯再加 sinangil
- n_chiplets: 12
- seed=22  N_STEPS=16  T0=2.0  alpha=0.95
- gaps: any 0.5 mm, GPU–GPU 2.0 mm
- coupled J = (Tmax-45)/25.0 + 1.0*vm_chip_max/100.0 + 0.5*w_pv/4.0
- vm = worst interior von Mises among normal chiplets (CPU/GPU, CIM area≥1 mm²)
- ATSim package overrides: {"pitch_die_mm": 2.0, "pitch_coarse_mm": 8.0}
- one scheme at a time (12 GB); CG only, no spsolve

## Run log

- start outline 59.12 x 44.98 mm  min_gap=0.59 mm  GPU-GPU=2.62 mm
- `heat_only` init  Tmax=70.540C  vm_chipmax=150.22 MPa (GPU_TITANX)  w_pv=4.241 um  J=1.0216  stress=123.6s
  `heat_only` step   1 accept rotate CPU_armL2  Tmax=70.540C  J=1.0216
  `heat_only` step   3 accept nudge GPU_3070 (-1.0,0.0)mm  Tmax=70.570C  J=1.0228
  `heat_only` step   6 accept swap CIM_wang<->CIM_colonnade  Tmax=70.570C  J=1.0228
  `heat_only` step   9 accept nudge GPU_3070 (-2.0,0.0)mm  Tmax=70.670C  J=1.0268
  `heat_only` step  10 accept nudge GPU_TITANX (0.0,0.5)mm  Tmax=70.750C  J=1.0300
  `heat_only` step  12 accept nudge CPU_armL2 (0.0,1.0)mm  Tmax=70.750C  J=1.0300
  `heat_only` step  15 accept nudge GPU_2060S (0.0,-0.5)mm  Tmax=70.760C  J=1.0304
  `heat_only` step  16 accept nudge CPU_starter (-1.5,0.0)mm  Tmax=70.770C  J=1.0308
- `heat_only` rescoring best with silicon-fill stress …
- `heat_only` done  Tmax 70.540→70.540C  vm 150.22→150.22 MPa  wpv 4.241→4.241 um  accepted 8/16
- `coupled` init  Tmax=70.540C  vm_chipmax=150.22 MPa (GPU_TITANX)  w_pv=4.241 um  J=3.0539  stress=124.1s
  `coupled` step   1 BEST rotate CPU_armL2  Tmax=70.540C  vm=150.22  wpv=4.241  J=3.0539
  `coupled` step   3 BEST nudge GPU_3070 (-1.0,0.0)mm  Tmax=70.570C  vm=149.47  wpv=4.253  J=3.0492
  `coupled` step   6 accept swap CIM_wang<->CIM_colonnade  Tmax=70.570C  J=3.0492
  `coupled` step   9 BEST nudge GPU_3070 (-2.0,0.0)mm  Tmax=70.670C  vm=121.94  wpv=3.797  J=2.7208
  `coupled` step  10 accept nudge GPU_TITANX (0.0,0.5)mm  Tmax=70.750C  J=2.7771
  `coupled` step  12 accept nudge CPU_armL2 (0.0,1.0)mm  Tmax=70.750C  J=2.7770
  `coupled` step  15 accept nudge GPU_2060S (0.0,-0.5)mm  Tmax=70.760C  J=2.8283
  `coupled` step  16 accept nudge CPU_starter (-1.5,0.0)mm  Tmax=70.770C  J=2.8428
- `coupled` done  Tmax 70.540→70.670C  vm 150.22→121.94 MPa  wpv 4.241→3.797 um  accepted 8/16

## Table
| arm | Tmax init→best | vm_chip_max init→best | w_pv init→best | accepted |
|---|---|---|---|---|
| heat_only | 70.540→70.540 | 150.22→150.22 (GPU_TITANX) | 4.241→4.241 | 8/16 |
| coupled | 70.540→70.670 | 150.22→121.94 (GPU_TITANX) | 4.241→3.797 | 8/16 |

coupled − heat_only at best: dTmax=0.130 C, d_vm=-28.28 MPa, d_wpv=-0.444 um

wall 1492.3s  finished 2026-10-25T13:43:33.639051+00:00
