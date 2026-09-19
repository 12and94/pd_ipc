// core/solver/EigenDirectSolver.cpp
#include "core/solver/EigenDirectSolver.h"

#include <chrono>

namespace pd {

namespace {
using Clock = std::chrono::steady_clock;
double secondsSince(const Clock::time_point& t0) {
  return std::chrono::duration<double>(Clock::now() - t0).count();
}
}  // namespace

struct EigenDirectSolver::Impl {
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<Scalar>, Eigen::Lower, Eigen::AMDOrdering<int>> llt;
};

EigenDirectSolver::EigenDirectSolver() : impl_(std::make_unique<Impl>()) {}
EigenDirectSolver::~EigenDirectSolver() = default;

void EigenDirectSolver::analyze(int n, const Eigen::SparseMatrix<Scalar>& patternMatrix) {
  const auto t0 = Clock::now();

  // 符号分解只依赖稀疏结构：直接用真实矩阵的结构（见接口注释里的原因）。
  // 结构相同（拓扑不变）时可以永久复用，因此本函数正常只被调用一次。
  impl_->llt.analyzePattern(patternMatrix);

  analyzed_ = true;
  stats_.analyzeCalls += 1;
  stats_.lastAnalyzeSeconds = secondsSince(t0);
  stats_.nnz = static_cast<long long>(patternMatrix.nonZeros());
}

void EigenDirectSolver::factorize(const Eigen::SparseMatrix<Scalar>& L) {
  const auto t0 = Clock::now();
  impl_->llt.factorize(L);
  stats_.factorizeCalls += 1;
  stats_.lastFactorizeSeconds = secondsSince(t0);
  stats_.totalFactorizeSeconds += stats_.lastFactorizeSeconds;
}

void EigenDirectSolver::solve(const Eigen::VectorXd& b, Eigen::VectorXd& x) {
  const auto t0 = Clock::now();
  x = impl_->llt.solve(b);
  stats_.solveCalls += 1;
  stats_.lastSolveSeconds = secondsSince(t0);
  stats_.totalSolveSeconds += stats_.lastSolveSeconds;
}

}  // namespace pd
