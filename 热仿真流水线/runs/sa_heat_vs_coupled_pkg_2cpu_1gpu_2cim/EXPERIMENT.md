# Heat-only vs heat+stress SA — pkg_2cpu_1gpu_2cim

- started (UTC): 2026-10-23T08:07:27.715873+00:00
- title: 2 CPU + 1 GPU + 2 CIM（5 颗）
- why: 小规模基线：两颗 CPU、一颗 3070、jia/colonnade
- n_chiplets: 5
- seed=22  N_STEPS=16  T0=2.0  alpha=0.95
- gaps: any 0.5 mm, GPU–GPU 2.0 mm
- coupled J = (Tmax-45)/25.0 + 1.0*vm_chip_max/100.0 + 0.5*w_pv/4.0
- vm = worst interior von Mises among normal chiplets (CPU/GPU, CIM area≥1 mm²)
- ATSim package overrides: {"pitch_die_mm": 2.0, "pitch_coarse_mm": 8.0}
- one scheme at a time (12 GB); CG only, no spsolve

## Run log

- start outline 22.97 x 29.53 mm  min_gap=0.64 mm  GPU-GPU=inf mm
- `heat_only` init  Tmax=55.490C  vm_chipmax=34.74 MPa (CIM_jia)  w_pv=0.744 um  J=0.4196  stress=18.4s
  `heat_only` step   1 accept rotate CPU_threads  Tmax=55.500C  J=0.4200
  `heat_only` step   2 BEST jump GPU_3070  Tmax=55.480C  vm=34.74  wpv=0.744  J=0.4192
  `heat_only` step   4 accept rotate CPU_threads  Tmax=55.480C  J=0.4192
  `heat_only` step   5 accept swap CIM_jia<->CPU_matmul  Tmax=55.490C  J=0.4196
  `heat_only` step   6 accept nudge GPU_3070 (-0.5,0.0)mm  Tmax=55.490C  J=0.4196
  `heat_only` step   7 accept jump CIM_colonnade  Tmax=55.480C  J=0.4192
  `heat_only` step   9 accept rotate CIM_colonnade  Tmax=55.480C  J=0.4192
  `heat_only` step  10 accept swap CPU_threads<->CIM_jia  Tmax=55.480C  J=0.4192
  `heat_only` step  11 accept jump GPU_3070  Tmax=55.520C  J=0.4208
  `heat_only` step  12 accept rotate CPU_threads  Tmax=55.520C  J=0.4208
  `heat_only` step  13 accept nudge GPU_3070 (0.0,-1.5)mm  Tmax=55.510C  J=0.4204
  `heat_only` step  14 accept rotate CIM_jia  Tmax=55.510C  J=0.4204
  `heat_only` step  15 accept rotate GPU_3070  Tmax=55.490C  J=0.4196
  `heat_only` step  16 accept swap CPU_matmul<->CPU_threads  Tmax=55.490C  J=0.4196
- `heat_only` rescoring best with silicon-fill stress …
- `heat_only` done  Tmax 55.490→55.480C  vm 34.74→36.92 MPa  wpv 0.744→0.933 um  accepted 14/16
- `coupled` init  Tmax=55.490C  vm_chipmax=34.74 MPa (CIM_jia)  w_pv=0.744 um  J=0.8600  stress=17.8s
  `coupled` step   1 accept rotate CPU_threads  Tmax=55.500C  J=0.8618
  `coupled` step   2 accept jump GPU_3070  Tmax=55.480C  J=0.9050
  `coupled` step   4 accept rotate CPU_threads  Tmax=55.480C  J=0.9007
  `coupled` step   5 accept swap CIM_jia<->CPU_matmul  Tmax=55.490C  J=0.9031
  `coupled` step   6 accept nudge GPU_3070 (-0.5,0.0)mm  Tmax=55.490C  J=0.9016
  `coupled` step   7 accept jump CIM_colonnade  Tmax=55.480C  J=0.9006
  `coupled` step   9 accept rotate CIM_colonnade  Tmax=55.480C  J=0.9006
  `coupled` step  10 accept swap CPU_threads<->CIM_jia  Tmax=55.480C  J=0.8968
  `coupled` step  11 accept jump GPU_3070  Tmax=55.520C  J=0.8856
  `coupled` step  12 accept rotate CPU_threads  Tmax=55.520C  J=0.8856
  `coupled` step  13 accept nudge GPU_3070 (0.0,-1.5)mm  Tmax=55.510C  J=0.8977
  `coupled` step  14 accept rotate CIM_jia  Tmax=55.510C  J=0.9142
  `coupled` step  15 accept rotate GPU_3070  Tmax=55.490C  J=0.8771
  `coupled` step  16 accept swap CPU_matmul<->CPU_threads  Tmax=55.490C  J=0.8773
- `coupled` done  Tmax 55.490→55.490C  vm 34.74→34.74 MPa  wpv 0.744→0.744 um  accepted 14/16

## Table
| arm | Tmax init→best | vm_chip_max init→best | w_pv init→best | accepted |
|---|---|---|---|---|
| heat_only | 55.490→55.480 | 34.74→36.92 (GPU_3070) | 0.744→0.933 | 14/16 |
| coupled | 55.490→55.490 | 34.74→34.74 (CIM_jia) | 0.744→0.744 | 14/16 |

coupled − heat_only at best: dTmax=0.010 C, d_vm=-2.19 MPa, d_wpv=-0.189 um

wall 344.5s  finished 2026-10-23T08:13:12.193826+00:00
