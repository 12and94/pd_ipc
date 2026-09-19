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

}  // namespace pd
