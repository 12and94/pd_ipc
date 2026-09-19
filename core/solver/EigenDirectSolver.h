// core/solver/EigenDirectSolver.h
// 基于 Eigen SimplicialLDLT 的稀疏直接法后端（CPU 版本唯一的直接法实现）。
//
// 为什么是 LDLᵀ 而不是 LLᵀ(Cholesky)：矩阵数值跨多个量级（刚度 κ 与 M/h² 相差很大），
// LDLᵀ 在数值上更稳；且 Eigen 的 SimplicialLDLT 支持 analyzePattern / factorize 分离调用，
// 正好对应本项目的"符号分解只做一次、数值分解只在 stamp 变化时做"的用法。
//
// 诚实说明（docs/plan.md D12）：Eigen 的稀疏三角求解是**单线程**的，
// 因此全局步的并行加速无从谈起；全局步的快来自"分解复用"而不是并行。
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
