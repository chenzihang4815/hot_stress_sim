"""给 python 求解器补上 HotSpot 的均热片外伸环，和 HotSpot 对同一套硅填充。

外伸不是一整圈网格。HotSpot 网格模型在 64x64 之外加 12 个节点：
均热片北南西东 4 个，热沉内圈 4 个，热沉外圈 4 个。
这套封装均热片和热沉同为 96 mm，热沉外圈的面积是 0，
populate_package_R 的 r_hs_per 在这里会除零，所以外圈的 4 个节点不建，
只留 8 个。外伸条按真实几何：宽 12 mm，长 84 mm。
边上每一格分到 g_from_edge/N，外伸条向上到热沉节点，热沉节点到环境。

层堆叠和 run_fill_superposition.py 相同：芯片、TIM、均热片、热沉。
没有中介层，没有冷边。封装 96 mm 对 72 mm，对流热阻 0.0390625 K/W。
"""
import json
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve

N = 64
SIDE = 0.072
CW = SIDE / N
T_AMB = 318.15
K_SI = 130.0
S_SPREADER = 0.096
S_SINK = 0.096
T_SPREADER = 0.001
T_SINK = 0.0069
K_CU = 400.0
R_CONVEC = 0.0390625
T_CHIP = 0.00015
T_TIM = 2.0e-05
K_TIM = 4.0

# 和 fill_superposition 的 DIES 同一套，行号已经换成 HotSpot 的（第 0 行在 y 最大处）。
DIES = (
    {"name": "cpu", "j0": 16, "j1": 32, "i0": 24, "i1": 40, "power_W": 80.0},
    {"name": "gpu", "j0": 36, "j1": 48, "i0": 24, "i1": 36, "power_W": 40.0},
    {"name": "cim", "j0": 24, "j1": 32, "i0": 44, "i1": 52, "power_W": 15.0},
)
SP_N, SP_S, SP_E, SP_W = 0, 1, 2, 3
SINK_C_N, SINK_C_S, SINK_C_E, SINK_C_W = 4, 5, 6, 7
SINK_N, SINK_S, SINK_E, SINK_W = 8, 9, 10, 11
N_EXTRA = 8  # 热沉和均热片同宽，没有外圈的 4 个节点


def getr(k, thickness, area):
    return thickness / (k * area)


def overhang(s_outer, s_inner, k, thickness):
    """一边的外伸条。均热片和热沉同宽时这条面积是 0，HotSpot 的公式会除零，这里跳过。"""
    span = (s_outer - s_inner) / 2.0
    if span <= 0.0:
        return None
    length = (s_outer + s_inner) / 2.0
    return {
        "g_from_edge": 1.0 / getr(k, span / 2.0, length * thickness),
        "g_up": 1.0 / getr(k, thickness, span * length),
        "area": span * length,
    }


def package():
    width = height = SIDE
    k_sp, k_hs = K_CU, K_CU
    p = {}
    p["r_sp1_x"] = getr(k_sp, (S_SPREADER - width) / 4.0,
                        (S_SPREADER + 3 * height) / 4.0 * T_SPREADER)
    p["r_sp1_y"] = getr(k_sp, (S_SPREADER - height) / 4.0,
                        (S_SPREADER + 3 * width) / 4.0 * T_SPREADER)
    p["r_hs1_x"] = getr(k_hs, (S_SPREADER - width) / 4.0,
                        (S_SPREADER + 3 * height) / 4.0 * T_SINK)
    p["r_hs1_y"] = getr(k_hs, (S_SPREADER - height) / 4.0,
                        (S_SPREADER + 3 * width) / 4.0 * T_SINK)
    p["r_hs2_x"] = getr(k_hs, (S_SPREADER - width) / 4.0,
                        (3 * S_SPREADER + height) / 4.0 * T_SINK)
    p["r_hs2_y"] = getr(k_hs, (S_SPREADER - height) / 4.0,
                        (3 * S_SPREADER + width) / 4.0 * T_SINK)
    p["r_hs"] = getr(k_hs, (S_SINK - S_SPREADER) / 4.0,
                     (S_SINK + 3 * S_SPREADER) / 4.0 * T_SINK)
    p["r_sp_per_x"] = getr(k_sp, T_SPREADER,
                           (S_SPREADER + height) * (S_SPREADER - width) / 4.0)
    p["r_sp_per_y"] = getr(k_sp, T_SPREADER,
                           (S_SPREADER + width) * (S_SPREADER - height) / 4.0)
    p["r_hs_c_per_x"] = getr(k_hs, T_SINK,
                             (S_SPREADER + height) * (S_SPREADER - width) / 4.0)
    p["r_hs_c_per_y"] = getr(k_hs, T_SINK,
                             (S_SPREADER + width) * (S_SPREADER - height) / 4.0)
    p["r_hs_per"] = getr(k_hs, T_SINK,
                         (S_SINK * S_SINK - S_SPREADER * S_SPREADER) / 4.0)
    p["r_amb_c_per_x"] = R_CONVEC * (S_SINK * S_SINK) / (
        (S_SPREADER + height) * (S_SPREADER - width) / 4.0)
    p["r_amb_c_per_y"] = R_CONVEC * (S_SINK * S_SINK) / (
        (S_SPREADER + width) * (S_SPREADER - height) / 4.0)
    p["r_amb_per"] = R_CONVEC * (S_SINK * S_SINK) / (
        (S_SINK * S_SINK - S_SPREADER * S_SPREADER) / 4.0)
    return p


def rx(k, thickness):
    return CW / (k * CW * thickness)


def rz(k, thickness):
    return thickness / (k * CW * CW)


def build(k_chip, with_ring):
    """四层网格，外加 8 个外伸节点。返回矩阵和热沉每格到环境的电导。"""
    nl = 4
    ngrid = nl * N * N
    nunk = ngrid + (N_EXTRA if with_ring else 0)
    acc = [[] for _ in range(nunk)]

    def add(a, b, g):
        acc[a].append((b, -g))
        acc[b].append((a, -g))
        acc[a].append((a, g))
        acc[b].append((b, g))

    def gid(layer, i, j):
        return (layer * N + i) * N + j

    t_chip = T_CHIP
    for i in range(N):
        for j in range(N):
            if j + 1 < N:
                r = 0.5 * rx(k_chip[i, j], t_chip) + 0.5 * rx(k_chip[i, j + 1], t_chip)
                add(gid(0, i, j), gid(0, i, j + 1), 1.0 / r)
            if i + 1 < N:
                r = 0.5 * rx(k_chip[i, j], t_chip) + 0.5 * rx(k_chip[i + 1, j], t_chip)
                add(gid(0, i, j), gid(0, i + 1, j), 1.0 / r)
    for layer, (t, k) in ((1, (T_TIM, K_TIM)), (2, (T_SPREADER, K_CU)),
                          (3, (T_SINK, K_CU))):
        g = 1.0 / rx(k, t)
        for i in range(N):
            for j in range(N):
                if j + 1 < N:
                    add(gid(layer, i, j), gid(layer, i, j + 1), g)
                if i + 1 < N:
                    add(gid(layer, i, j), gid(layer, i + 1, j), g)
    for i in range(N):
        for j in range(N):
            add(gid(0, i, j), gid(1, i, j), 1.0 / rz(k_chip[i, j], t_chip))
    for layer, (t, k) in ((1, (T_TIM, K_TIM)), (2, (T_SPREADER, K_CU))):
        g = 1.0 / rz(k, t)
        for i in range(N):
            for j in range(N):
                add(gid(layer, i, j), gid(layer + 1, i, j), g)

    area = CW * CW
    g_amb = 1.0 / (rz(K_CU, T_SINK) + R_CONVEC * (S_SINK * S_SINK) / area)
    for i in range(N):
        for j in range(N):
            a = gid(3, i, j)
            acc[a].append((a, g_amb))

    ring_g = {}
    if with_ring:
        base = ngrid
        ry = rx(K_CU, T_SPREADER)
        rx_hs = rx(K_CU, T_SINK)
        sp_over = overhang(S_SPREADER, SIDE, K_CU, T_SPREADER)
        hs_over = overhang(S_SINK, S_SPREADER, K_CU, T_SINK)
        ring_g = {"sp": sp_over, "hs": hs_over}

        def link(edge_cells, extra, r_edge, r_ring, n_cells):
            g = 1.0 / (r_edge / 2.0 + n_cells * r_ring)
            for cell in edge_cells:
                add(cell, base + extra, g)

        sp_north = [gid(2, 0, j) for j in range(N)]
        sp_south = [gid(2, N - 1, j) for j in range(N)]
        sp_west = [gid(2, i, 0) for i in range(N)]
        sp_east = [gid(2, i, N - 1) for i in range(N)]
        hs_north = [gid(3, 0, j) for j in range(N)]
        hs_south = [gid(3, N - 1, j) for j in range(N)]
        hs_west = [gid(3, i, 0) for i in range(N)]
        hs_east = [gid(3, i, N - 1) for i in range(N)]
        g_edge = sp_over["g_from_edge"] / N
        for cells, node in ((sp_north, SP_N), (sp_south, SP_S),
                            (sp_west, SP_W), (sp_east, SP_E)):
            for cell in cells:
                add(cell, base + node, g_edge)
        g_up = sp_over["g_up"]
        for a, b in ((SP_N, SINK_C_N), (SP_S, SINK_C_S), (SP_W, SINK_C_W), (SP_E, SINK_C_E)):
            add(base + a, base + b, g_up)
        r_metal = T_SINK / (K_CU * sp_over["area"])
        r_air = R_CONVEC * (S_SINK * S_SINK) / sp_over["area"]
        g_to_air = 1.0 / (r_metal + r_air)
        for node in (SINK_C_N, SINK_C_S, SINK_C_W, SINK_C_E):
            a = base + node
            acc[a].append((a, g_to_air))
        ring_g["to_air"] = g_to_air

    rows, cols, data = [], [], []
    for a, terms in enumerate(acc):
        merged = {}
        for b, v in terms:
            merged[b] = merged.get(b, 0.0) + v
        for b, v in merged.items():
            rows.append(a)
            cols.append(b)
            data.append(v)
    mat = sparse.csr_matrix((data, (rows, cols)), shape=(nunk, nunk))
    return mat, g_amb, ring_g


def solve(with_ring):
    k = np.full((N, N), K_SI)
    power = np.zeros((N, N))
    for d in DIES:
        ncell = (d["i1"] - d["i0"]) * (d["j1"] - d["j0"])
        power[d["i0"]:d["i1"], d["j0"]:d["j1"]] = d["power_W"] / ncell
    mat, g_amb, ring_g = build(k, with_ring)
    ngrid = 4 * N * N
    rhs = np.zeros(mat.shape[0])
    rhs[3 * N * N:4 * N * N] = g_amb * T_AMB
    rhs[:N * N] = power.ravel()
    if with_ring:
        for node in (SINK_C_N, SINK_C_S, SINK_C_W, SINK_C_E):
            rhs[ngrid + node] = ring_g["to_air"] * T_AMB
    x = spsolve(mat, rhs)
    return x[:N * N].reshape(N, N)


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import font_manager

    base = Path("/home/czh/任务仿真产生功耗文件工具/热仿真流水线/runs/fill_superposition")
    hs = np.load(base / "hotspot_silicon.npy")
    plain = solve(False)
    ring = solve(True)

    def stats(name, field):
        diff = field - hs
        return {
            "name": name,
            "peak_C": float(field.max() - 273.15),
            "min_C": float(field.min() - 273.15),
            "peak_minus_hotspot_K": float(field.max() - hs.max()),
            "mean_minus_hotspot_K": float(diff.mean()),
            "max_abs_K": float(np.max(np.abs(diff))),
            "rmse_K": float(np.sqrt(np.mean(diff ** 2))),
        }

    rows = [stats("不加外伸环", plain), stats("加 8 个外伸节点", ring)]
    rows.append({
        "name": "HotSpot",
        "peak_C": float(hs.max() - 273.15),
        "min_C": float(hs.min() - 273.15),
    })
    out_dir = Path("/home/czh/任务仿真产生功耗文件工具/论文证据/填充场叠加/复现/spreader_ring")
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "compare.json").write_text(
        json.dumps({"dies": list(DIES), "rows": rows}, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8")
    np.savez(out_dir / "fields.npz", hotspot=hs, plain=plain, ring=ring)
    for r in rows:
        print(f"{r['name']:16} 峰值 {r['peak_C']:.3f}  最低 {r['min_C']:.3f}"
              + (f"  峰值差 {r['peak_minus_hotspot_K']:+.3f}  "
                 f"逐格最大 {r['max_abs_K']:.3f}  均方根 {r['rmse_K']:.3f}"
                 if "max_abs_K" in r else ""), flush=True)

    font_manager.fontManager.addfont("/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc")
    plt.rcParams["font.family"] = "Noto Sans CJK JP"
    plt.rcParams["axes.unicode_minus"] = False
    diff = ring - hs
    fig, axes = plt.subplots(1, 3, figsize=(13.2, 4.6))
    panels = (
        (hs - 273.15, f"HotSpot\n峰值 {hs.max()-273.15:.2f} °C", "inferno", None),
        (ring - 273.15, f"python 加外伸环\n峰值 {ring.max()-273.15:.2f} °C", "inferno", None),
        (diff, f"python 减 HotSpot\n最大 {np.max(np.abs(diff)):.3f} K", "coolwarm",
         float(np.max(np.abs(diff)))),
    )
    lo = min(float((hs - 273.15).min()), float((ring - 273.15).min()))
    hi = max(float(hs.max()), float(ring.max())) - 273.15
    for ax, (fld, title, cmap, lim) in zip(axes, panels):
        kw = {"vmin": -lim, "vmax": lim} if lim else {"vmin": lo, "vmax": hi}
        im = ax.imshow(fld[::-1], origin="lower", extent=[0, 72, 0, 72], cmap=cmap, **kw)
        ax.set_title(title, fontsize=11)
        ax.set_xlabel("mm")
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    axes[0].set_ylabel("mm")
    fig.suptitle("同一布局、同一功耗、整层硅。python 补上 8 个外伸节点", fontsize=12)
    fig.tight_layout()
    fig_path = Path("/home/czh/任务仿真产生功耗文件工具/论文证据/填充场叠加/图/python对hotspot_加外伸环.png")
    fig.savefig(fig_path, dpi=140)
    print("写到", fig_path, flush=True)


if __name__ == "__main__":
    main()
