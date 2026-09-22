// core/solver/EigenDirectSolver.h
// 稀疏直接法后端（CPU 版本唯一的直接法实现）：分解用 Eigen 的 SimplicialLDLT，
// **排序 / 置换 / 前代 / 回代由本文件自己实现**（2026-09-22 起，见 docs/perf.md §8）。
//
// 为什么是 LDLᵀ 而不是 LLᵀ(Cholesky)：矩阵数值跨多个量级（刚度 κ 与 M/h² 相差很大），
// LDLᵀ 在数值上更稳；且 Eigen 的 SimplicialLDLT 支持 analyzePattern / factorize 分离调用，
// 正好对应本项目的"符号分解只做一次、数值分解只在 stamp 变化时做"的用法。
// 实测还多一条理由：同样因子结构下 LDLT 的回代比 LLT 快 24–31 %（docs/solver-feasibility.md §5）。
//
// 为什么不用 `SimplicialLDLT::solve()`：它每次调用都要做两次置换、重算 n 次除法
// （`m_diag.asDiagonal().inverse()`）、并对因子的**转置视图**做回代。自己实现这三步后
// 回代快 10–20 %（回代占单子步 76 % ⇒ 整体 −8…−18 %），而**数值结果与它逐位相同**
// （累加顺序没变）。详见 docs/perf.md §8。
//
// 诚实说明（docs/plan.md D12）：稀疏三角求解是**单线程**的（我们的实现也是），
// 因此全局步的并行加速无从谈起；全局步的快来自"分解复用"而不是并行。
// 并行回代的可行性也量过了：层调度原型只值 0.81×/1.45×/1.39×（40×40/100×100/200×200），
// 不值得做 —— 见 docs/solver-feasibility.md §3。
#pragma once

#include <memory>

#include <Eigen/Sparse>

#include "core/solver/IGlobalSolver.h"

namespace pd {

class EigenDirectSolver final : public IGlobalSolver {
 public:
  EigenDirectSolver();
  ~EigenDirectSolver() override;

  void analyze(int n, const Eigen::SparseMatrix<Scalar>& patternMatrix) override;
  void factorize(const Eigen::SparseMatrix<Scalar>& L) override;
  void solve(const Eigen::VectorXd& b, Eigen::VectorXd& x) override;
  const Stats& stats() const override { return stats_; }

  /// 是否已完成符号分解。
  bool analyzed() const override { return analyzed_; }

 private:
  struct Impl;
  std::unique_ptr<Impl> impl_;
  bool analyzed_ = false;
  Stats stats_;
};

}  // namespace pd
