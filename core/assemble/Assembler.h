// core/assemble/Assembler.h
// 把"约束贡献 + 惯性项 + pinned 行"组装成全局系统 L x = b。
//
// 关键约定（docs/plan.md §2.2.1 / §2.3）：
//   1. L 的数值只依赖 拓扑 / 刚度 κ / 质量 m / 子步长 h / pin 掩码，
//      **与位形无关** —— 所以正常运行时只需要分解一次，每帧只做回代。
//   2. 与 pinned 顶点相关的约束既不进 L 也不进 b；pinned 顶点所在的行被整行覆盖为
//      (对角 1, 右端 pin 位置)，这样解出的 pinned 位置**严格等于**目标位置。
//   3. 右端的散射量是 2 κ_c d_c（系数 2 不可省，见文档说明）。
#pragma once

#include <Eigen/Sparse>

#include "core/energy/DistanceTerm.h"
#include "core/math/Vec3.h"
#include "core/mesh/Mesh.h"
#include "core/solver/IGlobalSolver.h"

namespace pd {

/// 计算触发重新分解的判定戳。只做精确值比较，不做容差比较。
SolverStamp computeStamp(const Mesh& mesh, Scalar dt, Scalar damping, uint64_t topologyId);

/// 组装 L 的数值部分：M_γ/h² + Σ_c κ_c A_cᵀ A_c，并覆盖 pinned 行。
/// 其中 M_γ = M/(1-k_d) 是"带阻尼的有效质量"（见 Integrator.cpp 第 5 步说明）。
/// 只在 stamp 变化时调用（正常运行时每帧不调用）。
void assembleLeftHandSide(const Mesh& mesh, Scalar dt, Scalar damping,
                          Eigen::SparseMatrix<Scalar>& L);

/// 组装右端的惯性部分：b = (M_γ/h²)·x̂。
/// 重力**不在这里**：它只出现在预测位置 x̂ = x + hv + h²g 里（只出现一次）。
void assembleInertialRhs(const Mesh& mesh, const std::vector<Vec3>& predicted,
                         Scalar dt, Scalar damping, Eigen::VectorXd& b);

/// 覆盖 pinned 行：b 的对应分量置为 pin 位置。
void applyPinRhs(const Mesh& mesh, Eigen::VectorXd& b);

/// 便捷封装：把 Eigen::VectorXd 解包成 Vec3 数组。
void unpackPositions(const Eigen::VectorXd& x, std::vector<Vec3>& out);
void packPositions(const std::vector<Vec3>& positions, Eigen::VectorXd& out);

}  // namespace pd
