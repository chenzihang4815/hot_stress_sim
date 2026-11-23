"""一个设计点的全部可调量。

链路里每个被冻住的常数都收在这里，带它现在的取值和它属于哪一层。评估函数
只读 Design，不再读散落在各脚本里的模块级常量。各层的划分按量改变的是什么：

- package：热封装，全局标量，对所有芯粒一视同仁，场叠加照用。
- materials：结构和材料，决定应力约束线的位置，只在外层场景扫。
- operating：每颗芯粒的频率和电压，低维，场叠加改缓存键即可。
- architecture：芯片和任务的选择，等于换一颗芯片，不进搜索内循环。

布局不在这里，它是搜索的输出。任务分配和时序调度改变功耗的空间或时间
分布，场叠加失效，还没有对应的字段。
"""

from __future__ import annotations

from dataclasses import dataclass, field, asdict


@dataclass
class Package:
    """热封装。默认值就是 hotspot_io.DEFAULT_CONFIG 里写死的那一套。"""

    # 芯粒层的合法区域，退火在这个框里走。
    die_mm: float = 72.0
    spreader_mm: float = 96.0
    sink_mm: float = 96.0
    # 0.1 K/W 是 HotSpot 对 60 mm 热沉的默认值，按面积缩到现在的热沉。
    r_convec_K_per_W: float = 0.0390625
    # 界面层，也就是 TIM。普通硅脂约 4，液态金属和铟焊能到几十。
    tim_um: float = 20.0
    k_tim_W_per_mK: float = 4.0
    t_die_um: float = 150.0
    k_die_W_per_mK: float = 130.0
    t_spreader_mm: float = 1.0
    t_sink_mm: float = 6.9
    k_spreader_W_per_mK: float = 400.0
    k_sink_W_per_mK: float = 400.0
    ambient_C: float = 45.0
    # HotSpot 自带的微流道冷却，默认关。
    microfluidic: bool = False


@dataclass
class Materials:
    """结构和材料。改动会整体平移应力约束，不进每一步的搜索。"""

    # 芯粒之间的空隙填什么。silicon 是现在的做法，void 是挖空。
    fill: str = "silicon"
    t_interposer_um: float = 100.0
    t_substrate_mm: float = 1.0
    # 任意两颗的最小间距，GPU 之间更宽，怕两颗高功耗贴在一起。
    gap_mm: float = 0.5
    gap_gpu_mm: float = 2.0


@dataclass
class ChipOp:
    """一颗芯粒的工作点。s 是相对标称频率的倍率，v 是相对标称电压的倍率。"""

    s: float = 1.0
    v: float = 1.0


@dataclass
class Architecture:
    """芯片与任务的选择。改这里等于换芯片，不是调工作点。"""

    cpus: list = field(default_factory=list)
    gpus: list = field(default_factory=list)
    n_cim: int = 0
    hbm: bool = True
    # 空闲部件的漏电关掉多少。McPAT 的 power_gating 现在是关的。
    power_gating: float = 0.0


@dataclass
class Design:
    package: Package = field(default_factory=Package)
    materials: Materials = field(default_factory=Materials)
    operating: dict = field(default_factory=dict)
    architecture: Architecture = field(default_factory=Architecture)
    seed: int = 22

    def to_dict(self) -> dict:
        return asdict(self)
