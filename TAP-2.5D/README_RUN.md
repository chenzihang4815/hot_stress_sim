# TAP-2.5D 运行指南

## 📋 快速开始

### 1. 将文件复制到Linux系统

```bash
# 方法1：使用scp（从本地Windows到Linux）
scp -r "E:\南航研究生总文档\韩老师的保研后任务\HotSpot代码运行笔记\3D-IC_Workspace\TAP-2.5D" username@linux-server:/path/to/destination/

# 方法2：如果已经在Windows上通过WSL访问
# 直接在WSL中访问路径：
cd /mnt/e/南航研究生总文档/韩老师的保研后任务/HotSpot代码运行笔记/3D-IC_Workspace/TAP-2.5D
```

### 2. 进入目录并添加执行权限

```bash
cd TAP-2.5D

# 添加执行权限
chmod +x check_env.sh
chmod +x run_tap.sh
chmod +x quick_test.py
chmod +x util/hotspot
chmod +x util/*.pl
```

### 3. 检查环境

```bash
./check_env.sh
```

**预期输出：**
```
======================================
TAP-2.5D 环境依赖检查
======================================

1. 检查Python版本...
Python 3.8.10
✅ Python已安装

2. 检查CPLEX优化器...
CPLEX版本: 12.10
✅ CPLEX已安装

3. 检查HotSpot工具...
✅ HotSpot可执行文件存在

4. 检查Perl（用于生成热图）...
This is perl 5, version 30
```

---

## 🚀 运行方式

### 方式A：快速测试（推荐新手）

**不需要CPLEX，适合快速验证环境**

```bash
python3 quick_test.py
```

**输出示例：**
```
==================================================
TAP-2.5D 简化测试（仅热仿真）
==================================================

步骤1: 读取配置文件...
  ✅ 芯粒数量: 6
  ✅ 中介层尺寸: 50 mm

步骤2: 生成初始布局...
  芯粒坐标:
    Chiplet 0: (25.00, 25.00) mm
    Chiplet 1: (15.00, 35.00) mm
    ...

步骤3: 生成HotSpot输入文件...
  ✅ 已生成: ./output/ascend910/test_initial_L4_ChipLayer.flp

步骤4: 运行HotSpot热仿真...
  ✅ 热仿真完成
  最高温度: 85.32 °C
  耗时: 2.15 秒

步骤5: 估算线长（曼哈顿距离）...
  芯粒 0 -> 1: 距离=20.00mm, 带宽=200, 加权=4000.00
  ...
  ✅ 平均加权线长: 3250.50 mm
```

---

### 方式B：完整优化（需要CPLEX）

**运行完整的模拟退火优化**

```bash
./run_tap.sh
```

**或直接运行：**

```bash
python3 sim_annealing.py -c configs/sys_ascend910_local.cfg
```

**运行时间：**
- 小规模（6芯粒）：10-30分钟
- 中等规模（10芯粒）：30-60分钟
- 大规模（20+芯粒）：1-3小时

---

## 📊 查看结果

### 1. 最终结果

```bash
cat ./output/ascend910/output.txt
```

**输出格式：**
```
125                    # 优化步数
78.45                  # 最终温度 (°C)
8.32                   # 平均线长 (mm)
[25.3, 15.2, ...]      # 芯粒X坐标
[30.1, 20.5, ...]      # 芯粒Y坐标
```

### 2. 优化日志

```bash
cat ./output/ascend910/log.txt
```

### 3. 热图（如果有Perl）

```bash
# 查看SVG热图
ls ./output/ascend910/*.svg

# 在本地查看（使用浏览器）
# 或转换成PDF
convert output.svg output.pdf
```

---

## ⚠️ 常见问题

### 问题1：HotSpot执行权限错误

```bash
bash: ./util/hotspot: Permission denied
```

**解决：**
```bash
chmod +x ./util/hotspot
```

---

### 问题2：找不到CPLEX

```python
ModuleNotFoundError: No module named 'cplex'
```

**解决方案A：安装CPLEX（商业软件）**
```bash
# 从IBM官网下载
# https://www.ibm.com/products/ilog-cplex-optimization-studio

# 安装Python API
cd /path/to/cplex/python/3.8/
pip install cplex
```

**解决方案B：使用免费替代方案（修改routing.py）**
```bash
# 可以用开源优化器替代（如PuLP + CBC）
# 需要修改 routing.py 中的 CPLEX 调用
```

**解决方案C：跳过走线优化**
```bash
# 使用 quick_test.py 仅测试热仿真
python3 quick_test.py
```

---

### 问题3：HotSpot运行失败

```
Error: Cannot open floorplan file
```

**解决：**
```bash
# 检查文件是否存在
ls ./output/ascend910/*.flp

# 检查路径配置
cat configs/sys_ascend910_local.cfg | grep path
```

---

### 问题4：Python版本过低

```bash
Python 3.5.1
```

**TAP-2.5D需要Python 3.6+**

**解决：**
```bash
# Ubuntu/Debian
sudo apt-get update
sudo apt-get install python3.8

# CentOS/RHEL
sudo yum install python38

# 使用conda
conda create -n tap python=3.8
conda activate tap
```

---

## 🔧 自定义配置

### 创建你自己的测试案例

1. **复制配置文件**
```bash
cp configs/sys_ascend910_local.cfg configs/my_test.cfg
```

2. **编辑配置文件**
```bash
vim configs/my_test.cfg
```

3. **修改关键参数**
```ini
[general]
path = ./output/my_test/           # 输出目录
decay = 0.85                        # 退火衰减（越大收敛越慢）

[chiplets]
chiplet_count = 4                   # 芯粒数量
widths = 10, 8, 6, 6                # 宽度 (mm)
heights = 10, 8, 6, 6               # 高度 (mm)
powers = 50, 30, 20, 20             # 功耗 (W)

connections = 0, 100, 50, 0;       # 连接矩阵
              100, 0, 0, 0;
              50, 0, 0, 100;
              0, 0, 100, 0
```

4. **运行**
```bash
python3 sim_annealing.py -c configs/my_test.cfg
```

---

## 📈 进阶使用

### 覆盖命令行参数

```bash
# 不修改配置文件，直接指定参数
python3 sim_annealing.py \
    -c configs/sys_ascend910_local.cfg \
    -d ./output/custom/ \
    -g intp_size=45 \
    -g decay=0.9
```

### 查看实时进度

```bash
# 运行并在另一个终端监控日志
tail -f ./output/ascend910/step.txt
```

### 仅评估给定布局（不优化）

修改配置文件：
```ini
[general]
initial_placement = given  # 使用给定坐标
```

然后指定坐标：
```ini
[chiplets]
x = 25.0, 15.0, 35.0, ...
y = 25.0, 35.0, 15.0, ...
```

运行：
```bash
python3 config.py -c configs/my_test.cfg
```

---

## 📚 参考资料

- TAP-2.5D论文: https://github.com/bu-icsg/TAP-2.5D
- HotSpot手册: https://github.com/uvahotspot/hotspot
- CPLEX文档: https://www.ibm.com/docs/en/icos/22.1.1

---

## 💡 提示

1. **首次运行建议使用quick_test.py**验证环境
2. **优化时间较长**，建议在后台运行：
   ```bash
   nohup ./run_tap.sh > run.log 2>&1 &
   ```
3. **如果内存不足**，减小`chiplet_count`或增加`granularity`
4. **HotSpot结果不可靠时**，检查：
   - 芯粒是否重叠
   - 功耗是否合理
   - 材料参数是否正确

---

## 🎯 下一步

运行成功后，你可以：

1. **分析结果**
   - 对比不同布局的温度
   - 绘制热图
   - 分析线长分布

2. **修改算法**
   - 修改`sim_annealing.py`中的接受准则
   - 改进邻域搜索策略
   - 优化初始布局

3. **添加自己的芯粒系统**
   - 创建新的配置文件
   - 定义芯粒尺寸、功耗、连接关系

祝实验顺利！🎉
