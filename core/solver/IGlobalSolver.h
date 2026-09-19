// core/solver/IGlobalSolver.h
// 全局步求解器的抽象。三条实现路线（CPU 直接法 / GPU Chebyshev / GPU 局部+CPU 全局）
// 的差异全部集中在这里（见 docs/plan.md §1.1 的架构约定）。
//
// 三个操作粒度必须分开（docs/plan.md §2.4(2)）：
//   analyze   —— 符号分解与排序，只依赖稀疏结构（拓扑）=> 通常只做一次
//   factorize —— 数值分解，依赖矩阵数值（拓扑/刚度/h/pin）=> 只在 stamp 变化时做
//   solve     —— 前代 + 回代 => 每次全局步都做，是热路径
#pragma once

#include <vector>

#include <Eigen/Sparse>

#include "core/math/Vec3.h"
#include "core/solver/SparsePattern.h"

namespace pd {

/// 触发重新数值分解的条件集合。
/// L 的数值只依赖这些量（与位形无关），因此可以用"逐字段精确比较"判断是否需要重分解：
/// 相比浮点容差比较，精确比较不会因为"几乎没变"而误判，语义也更明确。
struct SolverStamp {
  uint64_t topologyId = 0;      ///< 拓扑标识（顶点数、边数、结构版本）
  uint64_t stiffnessBits = 0;   ///< 刚度集合的位模式哈希（Σ κ_c 的精确位）
  uint64_t massBits = 0;        ///< 质量集合的位模式哈希
  uint64_t dtBits = 0;          ///< 子步长 h 的位模式（h 进入 M/h²）
  uint64_t pinBits = 0;         ///< pin 掩码哈希
  uint64_t dampingBits = 0;     ///< 阻尼系数的位模式（阻尼进入有效质量 M/(1-k_d)）

  bool operator==(const SolverStamp& o) const {
    return topologyId == o.topologyId && stiffnessBits == o.stiffnessBits &&
           massBits == o.massBits && dtBits == o.dtBits && pinBits == o.pinBits &&
           dampingBits == o.dampingBits;
  }
  bool operator!=(const SolverStamp& o) const { return !(*this == o); }
};

class IGlobalSolver {
 public:
  virtual ~IGlobalSolver() = default;

  /// 符号分解：给定自由度数与**真实矩阵的结构**（而非理想化的块模式）。
  ///
  /// 为什么必须用真实矩阵：理想块模式（对角块 + 每条约束两个耦合块）会把
  /// "结构上恒为零"的位置也算进去（例如某顶点只与 pinned 顶点相连时，它的
  /// 非对角块会被跳过）。Eigen 的 SimplicialLDLT 在符号分解阶段就把这些位置
  /// 计入消去树，随后 factorize 一个结构不同的矩阵会让它在 solve 阶段访问
  /// 未定义的数据（实测表现为访问违例）。
  /// 用真实矩阵的结构做符号分解可确保两者严格一致。
  virtual void analyze(int n, const Eigen::SparseMatrix<Scalar>& patternMatrix) = 0;

  /// 数值分解：给定已组装好的 L。
  virtual void factorize(const Eigen::SparseMatrix<Scalar>& L) = 0;

  /// 是否已经完成符号分解（决定是否需要先调 analyze）。
  virtual bool analyzed() const = 0;

  /// 求解 L x = b。支持多右端（Eigen 的 solve(B)）以便将来批处理。
  virtual void solve(const Eigen::VectorXd& b, Eigen::VectorXd& x) = 0;

  /// 诊断：符号分解 / 数值分解 / 回代的累计次数与最近一次耗时（秒）。
  struct Stats {
    int analyzeCalls = 0;
    int factorizeCalls = 0;
    int solveCalls = 0;
    double lastAnalyzeSeconds = 0.0;
    double lastFactorizeSeconds = 0.0;
    double lastSolveSeconds = 0.0;
    double totalFactorizeSeconds = 0.0;
    double totalSolveSeconds = 0.0;
    long long nnz = 0;
  };
  virtual const Stats& stats() const = 0;
};

}  // namespace pd
