// core/sim/Integrator.h
// PD 时间步：预测 → (局部步 → 散射 → 全局步)^K → 速度更新 → 阻尼。
//
// 与 docs/plan.md §2.4(3) 的热路径逐条对应，刻意保持"一个函数一条流程"的直白写法，
// 便于对照文档与逐项插桩。
#pragma once

#include "core/sim/Scene.h"

namespace pd {

/// 执行一个子步（不负责子步循环与帧计时）。
/// 返回实际使用的 PD 迭代数。
int stepOnce(SimContext& ctx);

/// 执行一个完整帧：substepsPerFrame 个子步。
void stepFrame(SimContext& ctx);

/// 把预测、局部步、散射、组装、求解等各阶段的累计耗时清零。
void resetStageTimes(SimContext& ctx);

/// 当前代理能量（用于单调性与有界性检查）：
///   Σ_c (κ_c/2)(‖x_a-x_b‖ - l_c)²
Scalar surrogateEnergy(const SimContext& ctx);

/// 当前总能量 = 动能 + 弹性势能 + 重力势能（供能量有界性测试使用）。
Scalar totalEnergy(const SimContext& ctx);

/// 最大速度绝对值（诊断用）。
Scalar maxSpeed(const SimContext& ctx);

/// 当前位形的**非线性残差**（单位 m/s²）：最大顶点上的不平衡力除以该点质量。
///
/// 不动点处应为 0。它衡量的是"离本子步的解还有多远"，而不是"这一步动了多少"：
/// 高刚度下 PD 外层迭代收敛极慢（收缩因子 q = κ/(κ + m/h²)），相邻迭代位移
/// 很小但离正确解仍很远，此时只能用残差判断收敛质量。
/// 要求 ctx.targets 是当前位置的投影（stepOnce 每个迭代都会更新它）。
Scalar nonlinearResidual(const SimContext& ctx);

}  // namespace pd
