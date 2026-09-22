// core/assemble/Assembler.h
// 把"约束贡献 + 惯性项 + pinned 行"组装成全局系统 L x = b。
//
// 关键约定（docs/plan.md §2.2.1 / §2.3）：
//   1. L 的数值只依赖 拓扑 / 刚度 κ / 质量 m / 子步长 h / pin **掩码**，
//      **与位形无关**、也与 pin 的**位置**无关 —— 所以正常运行时只需要分解一次，
//      每帧只做回代；拖拽把手只改右端 b，不会触发重分解。阻尼同样不进 L。
//   2. **只有两端都 pinned 的约束才整条跳过**；单端 pinned 时，自由端的对角 +κ_c 仍然
//      要进 L（丢了这份刚度系统会明显偏软），散射时还要给自由端补上被消去的那一项
//      b_free += κ_c·x̄_pin（见 DistanceTerm.cpp 里的推导）。pinned 顶点所在的行被整行
//      覆盖为 (对角 1, 右端 pin 位置)，这样解出的 pinned 位置**严格等于**目标位置。
//   3. 右端的散射量是 **κ_c d_c**（b_a += κ_c d_c、b_b -= κ_c d_c），与矩阵块的 κ_c 同源；
//      **没有系数 2**（早期文档曾写 2κ_c d_c，是错的，已全局更正）。
#pragma once

#include <vector>

#include <Eigen/Sparse>

#include "core/energy/DistanceTerm.h"
#include "core/math/Vec3.h"
#include "core/mesh/Mesh.h"
#include "core/solver/IGlobalSolver.h"

namespace pd {

/// 计算触发重新分解的判定戳。只做精确值比较，不做容差比较。
SolverStamp computeStamp(const Mesh& mesh, Scalar dt, Scalar damping, uint64_t topologyId);

/// 计算"稀疏结构戳"：决定 L 的**非零结构**（而不是数值）的量。
///
/// 目前包含：拓扑标识 + 边表（端点对）+ pin **掩码**。
/// 为什么 pin 掩码会改变结构：`DistanceTerm::assembleMatrix` 只在**两端都自由**时写耦合块
/// `-κI`，所以给一个原本自由的顶点加 pin 会让那一行的非对角块消失。
/// 结构变了就必须重新做符号分解（analyzePattern），否则 factorize 一个结构不同的矩阵后，
/// solve 会访问未定义数据（见 IGlobalSolver.h 的说明；本项目因符号结构与真实矩阵不一致
/// 踩过一次访问违例）。
uint64_t computeStructureStamp(const Mesh& mesh, uint64_t topologyId);

/// 组装 L 的数值部分：M/h² + Σ_c κ_c A_cᵀ A_c + Σ_s k_s·w_p w_q·I₃，并覆盖 pinned 行。
/// 惯性项是**未缩放**的 M/h²；阻尼不进入 L（它只影响速度更新，见 Integrator.cpp 第 5 步）。
/// 只在 stamp 变化时调用（正常运行时每帧不调用）。
///
/// 为什么弯曲的右端是一个**常向量**：整条 `A_c` 是常数矩阵、投影 `p_c ≡ 0`（流形是线性子空间），
/// 所以每迭代唯一会变的那部分（散射）**根本不存在**。右端唯一多出来的项来自
/// "stencil 里有 pinned 顶点"时的消元补偿（推导见 .cpp），而 pinned 位置是**常数**
/// （把手不动时）⇒ 它可以在装配期一次算完，不进每迭代热路径、不需要任何额外 barrier。
void assembleLeftHandSide(const Mesh& mesh, Scalar dt, Scalar damping,
                          Eigen::SparseMatrix<Scalar>& L,
                          std::vector<Scalar>* bendingRhs = nullptr);

/// 只算弯曲约束的常向量右端（尺寸 3N）。**等价于**跑一遍 assembleLeftHandSide 的右端部分，
/// 代价只有 O(Σ_s 3)（不动矩阵、不建三元组）。
///
/// 为什么要单独有这个入口：`assembleLeftHandSide` 只在**数值 stamp 变化**时调用，而它的
/// 右端部分依赖 `mesh.pinPositions` —— 拖拽把手**只改 b、不改 L**（docs/plan.md §2.2.1 的
/// 明确约定，`stampChangesOnlyWhenLeftHandSideValuesChange` 把它钉住了）。于是"pin 位置变了
/// 但 stamp 没变"时，早先算好的常向量就过期了。所以这个函数每子步都要跑一次（O(B)，可忽略），
/// 而矩阵那一趟仍然只在 stamp 变化时跑。
void assembleBendingRhs(const Mesh& mesh, std::vector<Scalar>& out);

/// 组装右端的惯性部分：b = (M/h²)·x̂（**只有这一项**）。
/// 重力**不在这里**：它只出现在预测位置 x̂ = x + hv + h²g 里（只出现一次，
/// 再补一项 -Mg 就是重复计入）。
void assembleInertialRhs(const Mesh& mesh, const std::vector<Vec3>& predicted,
                         Scalar dt, Scalar damping, Eigen::VectorXd& b);

/// 覆盖 pinned 行：b 的对应分量置为 pin 位置。
void applyPinRhs(const Mesh& mesh, Eigen::VectorXd& b);

/// 只覆盖 pinned 行的**第 c 个分量**（配合全局步按 x/y/z 分量并行，见
/// IGlobalSolver::parallelComponents）。`components == 1` 时等价于 applyPinRhs。
///
/// 为什么可以拆：pin 的三个分量互不相干（各写各的 3v+c），而 L 的每一块都是标量 × I₃
/// ⇒ 分量之间内存完全不相交，可以由不同线程并发写同一向量 b 的不同位置。
void applyPinRhsComponent(const Mesh& mesh, Eigen::VectorXd& b, int c, int components);

/// 便捷封装：把 Eigen::VectorXd 解包成 Vec3 数组。
void unpackPositions(const Eigen::VectorXd& x, std::vector<Vec3>& out);
void packPositions(const std::vector<Vec3>& positions, Eigen::VectorXd& out);

}  // namespace pd
