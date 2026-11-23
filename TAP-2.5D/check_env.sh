#!/bin/bash
# TAP-2.5D 环境检查脚本

echo "======================================"
echo "TAP-2.5D 环境依赖检查"
echo "======================================"
echo ""

# 1. 检查Python版本
echo "1. 检查Python版本..."
python3 --version
if [ $? -ne 0 ]; then
    echo "❌ 未找到python3，请先安装Python 3.6+"
    exit 1
fi
echo "✅ Python已安装"
echo ""

# 2. 检查CPLEX
echo "2. 检查CPLEX优化器..."
python3 -c "import cplex; print('CPLEX版本:', cplex.__version__)" 2>/dev/null
if [ $? -ne 0 ]; then
    echo "⚠️  未找到CPLEX Python API"
    echo "   走线优化需要CPLEX，但你仍然可以运行热仿真部分"
    echo "   如需完整功能，请安装IBM CPLEX:"
    echo "   - 下载: https://www.ibm.com/products/ilog-cplex-optimization-studio"
    echo "   - 或使用教育版: https://www.ibm.com/academic/home"
else
    echo "✅ CPLEX已安装"
fi
echo ""

# 3. 检查HotSpot
echo "3. 检查HotSpot工具..."
if [ -f "./util/hotspot" ]; then
    if [ -x "./util/hotspot" ]; then
        echo "✅ HotSpot可执行文件存在"
    else
        echo "⚠️  HotSpot文件没有执行权限，正在添加..."
        chmod +x ./util/hotspot
    fi
else
    echo "❌ 未找到util/hotspot"
    exit 1
fi
echo ""

# 4. 检查Perl脚本
echo "4. 检查Perl（用于生成热图）..."
perl --version | head -1
if [ $? -ne 0 ]; then
    echo "⚠️  未找到Perl，无法生成热图SVG"
fi
echo ""

echo "======================================"
echo "环境检查完成！"
echo "======================================"
