# Heat-only vs heat+stress SA — pkg_5cpu_5gpu_5cim

- started (UTC): 2026-10-25T13:43:36.929177+00:00
- title: 5 CPU + 5 GPU + 5 CIM（15 颗）
- why: 14 芯再加第三颗 3070
- n_chiplets: 15
- seed=22  N_STEPS=16  T0=2.0  alpha=0.95
- gaps: any 0.5 mm, GPU–GPU 2.0 mm
- coupled J = (Tmax-45)/25.0 + 1.0*vm_chip_max/100.0 + 0.5*w_pv/4.0
- vm = worst interior von Mises among normal chiplets (CPU/GPU, CIM area≥1 mm²)
- ATSim package overrides: {"pitch_die_mm": 2.5, "pitch_coarse_mm": 10.0, "nz_die": 2, "nz_spreader": 2, "nz_sink": 3, "nz_tim": 1, "nz_substrate": 2, "nz_interposer": 2}
- one scheme at a time (12 GB); CG only, no spsolve

## Run log

- start outline 49.17 x 70.00 mm  min_gap=0.57 mm  GPU-GPU=2.25 mm
- `heat_only` init  Tmax=85.470C  vm_chipmax=253.00 MPa (GPU_3070b)  w_pv=7.325 um  J=1.6188  stress=127.4s
  `heat_only` step   2 BEST nudge GPU_TITANX (0.5,0.0)mm  Tmax=85.420C  vm=253.00  wpv=7.325  J=1.6168
  `heat_only` step   3 accept rotate CIM_wang  Tmax=85.420C  J=1.6168
  `heat_only` step   5 accept rotate CIM_jia  Tmax=85.420C  J=1.6168
  `heat_only` step   7 accept nudge GPU_3070c (1.5,0.0)mm  Tmax=85.480C  J=1.6192
  `heat_only` step  16 accept nudge CPU_armL2 (0.0,-0.5)mm  Tmax=85.480C  J=1.6192
- `heat_only` rescoring best with silicon-fill stress …
- `heat_only` done  Tmax 85.470→85.420C  vm 253.00→253.52 MPa  wpv 7.325→7.344 um  accepted 5/16
- `coupled` init  Tmax=85.470C  vm_chipmax=253.00 MPa (GPU_3070b)  w_pv=7.325 um  J=5.0644  stress=119.2s
  `coupled` step   2 accept nudge GPU_TITANX (0.5,0.0)mm  Tmax=85.420C  J=5.0700
  `coupled` step   3 accept rotate CIM_wang  Tmax=85.420C  J=5.0700
  `coupled` step   5 accept rotate CIM_jia  Tmax=85.420C  J=5.0701
  `coupled` step   7 accept nudge GPU_3070c (1.5,0.0)mm  Tmax=85.480C  J=5.0702
  `coupled` step  16 accept nudge CPU_armL2 (0.0,-0.5)mm  Tmax=85.480C  J=5.1308
- `coupled` done  Tmax 85.470→85.470C  vm 253.00→253.00 MPa  wpv 7.325→7.325 um  accepted 5/16

## Table
| arm | Tmax init→best | vm_chip_max init→best | w_pv init→best | accepted |
|---|---|---|---|---|
| heat_only | 85.470→85.420 | 253.00→253.52 (GPU_3070b) | 7.325→7.344 | 5/16 |
| coupled | 85.470→85.470 | 253.00→253.00 (GPU_3070b) | 7.325→7.325 | 5/16 |

coupled − heat_only at best: dTmax=0.050 C, d_vm=-0.52 MPa, d_wpv=-0.020 um

wall 1055.3s  finished 2026-10-25T14:01:12.254048+00:00
