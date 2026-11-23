#!/bin/bash
# TAP-2.5D 运行脚本

# 设置脚本在遇到错误时退出
set -e

echo "======================================"
echo "TAP-2.5D 芯粒布局优化工具"
echo "======================================"
echo ""

# 检查环境
echo "步骤1: 检查环境依赖..."
bash check_env.sh
echo ""

# 创建输出目录
echo "步骤2: 创建输出目录..."
mkdir -p ./output/ascend910
echo "✅ 输出目录已创建"
echo ""

# 运行TAP-2.5D
echo "步骤3: 开始运行优化..."
echo "----------------------------------------"
echo "配置文件: configs/sys_ascend910_local.cfg"
echo "输出目录: ./output/ascend910/"
echo "----------------------------------------"
echo ""

# 检查是否有CPLEX
python3 -c "import cplex" 2>/dev/null
if [ $? -eq 0 ]; then
    echo "✅ 检测到CPLEX，运行完整优化（热+线长）..."
    python3 sim_annealing.py -c configs/sys_ascend910_local.cfg
else
    echo "⚠️  未检测到CPLEX，运行简化版本（仅热优化）..."
    echo "   注意：走线部分将被跳过"
    python3 sim_annealing_no_cplex.py -c configs/sys_ascend910_local.cfg || \
    python3 sim_annealing.py -c configs/sys_ascend910_local.cfg 2>&1 | grep -v "cplex"
fi

echo ""
echo "======================================"
echo "运行完成！"
echo "======================================"
echo ""
echo "结果文件位置:"
echo "  - 布局坐标: ./output/ascend910/output.txt"
echo "  - 优化日志: ./output/ascend910/log.txt"
echo "  - 热图文件: ./output/ascend910/*.svg (如果有Perl)"
echo ""
echo "查看结果:"
echo "  cat ./output/ascend910/output.txt"
echo ""
