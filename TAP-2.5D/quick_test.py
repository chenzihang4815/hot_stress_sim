#!/usr/bin/env python3
"""
TAP-2.5D 简化测试版本
不需要CPLEX，仅测试热仿真和布局生成
"""

import config
import os
import sys

def main():
    print("=" * 50)
    print("TAP-2.5D 简化测试（仅热仿真）")
    print("=" * 50)
    print()

    # 读取配置
    print("步骤1: 读取配置文件...")
    system = config.read_config()
    print(f"  ✅ 芯粒数量: {system.chiplet_count}")
    print(f"  ✅ 中介层尺寸: {system.intp_size} mm")
    print()

    # 生成初始布局
    print("步骤2: 生成初始布局...")
    print(f"  ✅ 初始布局方式: bstree")
    print(f"  芯粒坐标:")
    for i in range(system.chiplet_count):
        print(f"    Chiplet {i}: ({system.x[i]:.2f}, {system.y[i]:.2f}) mm")
    print()

    # 生成HotSpot输入文件
    filename = "test_initial"
    print(f"步骤3: 生成HotSpot输入文件...")
    system.gen_flp(filename)
    system.gen_ptrace(filename)
    print(f"  ✅ 已生成: {system.path}{filename}_L4_ChipLayer.flp")
    print(f"  ✅ 已生成: {system.path}{filename}.ptrace")
    print()

    # 运行热仿真
    print("步骤4: 运行HotSpot热仿真...")
    import time
    start_time = time.time()
    try:
        temp = system.run_hotspot(filename)
        elapsed = time.time() - start_time
        print(f"  ✅ 热仿真完成")
        print(f"  最高温度: {temp:.2f} °C")
        print(f"  耗时: {elapsed:.2f} 秒")
    except Exception as e:
        print(f"  ❌ HotSpot运行失败: {e}")
        return 1
    finally:
        system.clean_hotspot(filename)
    print()

    # 简化的线长估算（曼哈顿距离）
    print("步骤5: 估算线长（曼哈顿距离）...")
    total_length = 0
    count = 0
    for i in range(system.chiplet_count):
        for j in range(i + 1, system.chiplet_count):
            if system.connection_matrix[i][j] > 0:
                dist = abs(system.x[i] - system.x[j]) + abs(system.y[i] - system.y[j])
                bandwidth = system.connection_matrix[i][j]
                total_length += dist * bandwidth
                count += 1
                print(f"  芯粒 {i} -> {j}: 距离={dist:.2f}mm, 带宽={bandwidth}, 加权={dist*bandwidth:.2f}")

    avg_length = total_length / count if count > 0 else 0
    print(f"  ✅ 平均加权线长: {avg_length:.2f} mm")
    print()

    # 保存结果
    output_file = system.path + "test_result.txt"
    with open(output_file, 'w') as f:
        f.write(f"温度: {temp:.2f} C\n")
        f.write(f"平均线长: {avg_length:.2f} mm\n")
        f.write(f"芯粒坐标:\n")
        for i in range(system.chiplet_count):
            f.write(f"  Chiplet {i}: ({system.x[i]:.2f}, {system.y[i]:.2f})\n")

    print("=" * 50)
    print("测试完成！")
    print("=" * 50)
    print(f"结果已保存到: {output_file}")
    print()
    print("下一步:")
    print("  1. 查看结果: cat output/ascend910/test_result.txt")
    print("  2. 如果有CPLEX，运行完整优化: bash run_tap.sh")
    print()

    return 0

if __name__ == "__main__":
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("\n\n用户中断")
        sys.exit(1)
    except Exception as e:
        print(f"\n❌ 错误: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
