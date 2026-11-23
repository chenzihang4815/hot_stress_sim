"""四层网格加上均热片外伸环。

层堆叠和 run_fill_superposition.py 相同：芯片 150 um、TIM 20 um k=4、
铜均热片 1 mm、铜热沉 6.9 mm。没有中介层，没有冷边。网格到 72 mm 截止，
72 mm 以外的散热由 8 个额外节点补上：均热片北南西东各一条外伸，
每条上面再挂一个热沉节点。

这套封装均热片和热沉同为 96 mm，热沉外圈的面积是 0。HotSpot 的
populate_package_R 在这里会除零，所以外圈的 4 个节点不建。
外伸条按真实几何：宽 (0.096-0.072)/2 = 12 mm，长 (0.096+0.072)/2 = 84 mm。
边上每一格分到 g_from_edge/N，外伸条向上到热沉节点是 g_up，
热沉节点到环境是金属热阻加对流热阻。

行号是 HotSpot 的，第 0 行在 y 最大处。
"""
import numpy as np
from scipy import sparse

N = 64
SIDE = 0.072
CW = SIDE / N
T_AMB = 318.15
K_SI = 130.0
K_FILL = 0.8
S_SPREADER = 0.096
S_SINK = 0.096
T_SPREADER = 0.001
T_SINK = 0.0069
K_CU = 400.0
R_CONVEC = 0.1 * (0.06 * 0.06) / (0.096 * 0.096)
T_CHIP = 0.00015
T_TIM = 2.0e-05
K_TIM = 4.0
NL = 4

SP_N, SP_S, SP_E, SP_W = 0, 1, 2, 3
SINK_C_N, SINK_C_S, SINK_C_E, SINK_C_W = 4, 5, 6, 7
N_EXTRA = 8
N_GRID = NL * N * N
N_UNK = N_GRID + N_EXTRA


def getr(k, thickness, area):
    return thickness / (k * area)


def _overhang():
    """一边的外伸条。均热片和热沉同宽，只有均热片这一条有面积。"""
    span = (S_SPREADER - SIDE) / 2.0
    length = (S_SPREADER + SIDE) / 2.0
    return {
        "g_from_edge": 1.0 / getr(K_CU, span / 2.0, length * T_SPREADER),
        "g_up": 1.0 / getr(K_CU, T_SPREADER, span * length),
        "area": span * length,
    }


OVER = _overhang()
G_EDGE = OVER["g_from_edge"] / N
G_UP = OVER["g_up"]
G_TO_AIR = 1.0 / (T_SINK / (K_CU * OVER["area"])
                  + R_CONVEC * (S_SINK * S_SINK) / OVER["area"])


def _rx(k, thickness):
    return CW / (k * CW * thickness)


def _rz(k, thickness):
    return thickness / (k * CW * CW)


def _gid(layer, i, j):
    return (layer * N + i) * N + j


def build(k_chip):
    """返回 (矩阵, 热沉每格到环境的电导)。矩阵含 8 个外伸节点。"""
    acc = [[] for _ in range(N_UNK)]

    def add(a, b, g):
        acc[a].append((b, -g))
        acc[b].append((a, -g))
        acc[a].append((a, g))
        acc[b].append((b, g))

    for i in range(N):
        for j in range(N):
            if j + 1 < N:
                r = 0.5 * _rx(k_chip[i, j], T_CHIP) + 0.5 * _rx(k_chip[i, j + 1], T_CHIP)
                add(_gid(0, i, j), _gid(0, i, j + 1), 1.0 / r)
            if i + 1 < N:
                r = 0.5 * _rx(k_chip[i, j], T_CHIP) + 0.5 * _rx(k_chip[i + 1, j], T_CHIP)
                add(_gid(0, i, j), _gid(0, i + 1, j), 1.0 / r)
    for layer, (t, k) in ((1, (T_TIM, K_TIM)), (2, (T_SPREADER, K_CU)),
                          (3, (T_SINK, K_CU))):
        g = 1.0 / _rx(k, t)
        for i in range(N):
            for j in range(N):
                if j + 1 < N:
                    add(_gid(layer, i, j), _gid(layer, i, j + 1), g)
                if i + 1 < N:
                    add(_gid(layer, i, j), _gid(layer, i + 1, j), g)
    for i in range(N):
        for j in range(N):
            add(_gid(0, i, j), _gid(1, i, j), 1.0 / _rz(k_chip[i, j], T_CHIP))
    for layer, (t, k) in ((1, (T_TIM, K_TIM)), (2, (T_SPREADER, K_CU))):
        g = 1.0 / _rz(k, t)
        for i in range(N):
            for j in range(N):
                add(_gid(layer, i, j), _gid(layer + 1, i, j), g)

    g_amb = 1.0 / (_rz(K_CU, T_SINK) + R_CONVEC * (S_SINK * S_SINK) / (CW * CW))
    for i in range(N):
        for j in range(N):
            a = _gid(3, i, j)
            acc[a].append((a, g_amb))

    base = N_GRID
    for cells, node in (
        ([_gid(2, 0, j) for j in range(N)], SP_N),
        ([_gid(2, N - 1, j) for j in range(N)], SP_S),
        ([_gid(2, i, 0) for i in range(N)], SP_W),
        ([_gid(2, i, N - 1) for i in range(N)], SP_E),
    ):
        for cell in cells:
            add(cell, base + node, G_EDGE)
    for a, b in ((SP_N, SINK_C_N), (SP_S, SINK_C_S),
                 (SP_W, SINK_C_W), (SP_E, SINK_C_E)):
        add(base + a, base + b, G_UP)
    for node in (SINK_C_N, SINK_C_S, SINK_C_W, SINK_C_E):
        a = base + node
        acc[a].append((a, G_TO_AIR))

    rows, cols, data = [], [], []
    for a, terms in enumerate(acc):
        merged = {}
        for b, v in terms:
            merged[b] = merged.get(b, 0.0) + v
        for b, v in merged.items():
            rows.append(a)
            cols.append(b)
            data.append(v)
    mat = sparse.csr_matrix((data, (rows, cols)), shape=(N_UNK, N_UNK))
    return mat, g_amb


def rhs_of(power, g_amb):
    """power 是芯粒层，第 0 行在 y 最大处。"""
    rhs = np.zeros(N_UNK)
    rhs[:N * N] = np.asarray(power).ravel()
    rhs[3 * N * N:4 * N * N] = g_amb * T_AMB
    for node in (SINK_C_N, SINK_C_S, SINK_C_W, SINK_C_E):
        rhs[N_GRID + node] = G_TO_AIR * T_AMB
    return rhs


def chip_layer(x):
    return np.asarray(x[:N * N]).reshape(N, N).copy()
