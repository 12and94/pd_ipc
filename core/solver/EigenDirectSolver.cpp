// core/solver/EigenDirectSolver.cpp
//
// LDLᵀ 稀疏直接法 + **自己实现的前代/回代**。
//
// 为什么不用 `SimplicialLDLT::solve()`（2026-09-22 的测量结论，见 docs/solver-feasibility.md）：
// 在同一轮会话内交替测量、各取最小值，手写的回代比 Eigen 的 solve 快：40×40 0.0749 vs 0.0968 ms、
// 100×100 0.8227 vs 0.9809、200×200 5.3961 vs 6.4449。回代占单子步 76 %，所以这一项直接决定整体。
// 差距来自 Eigen 那条通用路径的开销，都是可以省掉的：
//   · 每次 solve 都要 `dest = P·b` 与 `dest = P⁻¹·dest`（两次置换）—— 本文件自己管置换，同样要做，
//     但省掉了 PermutationMatrix 表达式的中间对象；
//   · `dest = m_diag.asDiagonal().inverse() * dest` —— 每次 solve 重新做 **n 次除法**；
//     这里在数值分解后预存 1/D（`dinv`），每次只乘；
//   · 回代走的是**存储因子的转置视图**（`matrixU()`）。本文件直接扫 CSC 的列（作用相同、访存更直）。
// 分解本身仍然交给 Eigen（符号/数值 LDLᵀ），语义、数值稳健性、pin 处理方式都不变。
//
// ---- 三个实现要点（都踩过，别改回去）----
//   1) 因子的存储方向决定了两个三角求解怎么写：
//        · 回代 Lᵀx = z 要的是 x[j] -= L(i,j)·x[i]（i>j）—— `L(i,j)` 正好在**第 j 列**（CSC），
//          所以直接扫列、按 i>j 筛，天然无冲突；
//        · 前代 L·y = bp 要的是 y[j] -= L(j,i)·y[i]（i<j）—— 按"扫第 j 列、按 i<j 筛"写会**一条都
//          命不中**（第 j 列装的是 L(k,j)，k>j），等于把前代写成恒等式 y = bp。
//          第一版就是这么错的，症状是"解差 5 %、残差 5.9 而 Eigen 是 1e-12"。
//          正确写法是 **scatter 形式**：扫第 j 列时把 y[j]（此刻已是终值）推给所有 i>j。
//   2) 两种形式的**累加顺序完全相同**（每个节点的各项都按"对方的列号递增"逐个相减），
//      因此结果与 Eigen 的 solve **逐位一致**：实测 40×40 与 200×200 两次运行的全部物理输出
//      （收敛统计 / 末次残差 / 应变 / 能量）与本改动前的二进制逐字相同。
//      这是本改造最重要的性质 —— 它不改变任何物理结果，只省时间。
//   3) 前代就地做在 `bp` 上（scatter 到 i>j，j 递增时 bp[j] 已终值），回代也就地做（j 递减时
//      bp[i]（i>j）已终值）。整个 solve 只用一块 n 维缓冲，且**不建因子的第二份副本**
//      （试过"另存一份 CSR 给前代按行 gather"，在真实管线的冷数据条件下反而更慢：
//      多读 0.8 MB（40×40）/ 39 MB（200×200）因子数据，把收益从 −23 % 压到 −10 %）。
#include "core/solver/EigenDirectSolver.h"

#include <chrono>
#include <cstdio>
#include <cstdlib>

namespace pd {

namespace {
using Clock = std::chrono::steady_clock;
double secondsSince(const Clock::time_point& t0) {
  return std::chrono::duration<double>(Clock::now() - t0).count();
}
}  // namespace

struct EigenDirectSolver::Impl {
  using CholMatrix = Eigen::SparseMatrix<Scalar, Eigen::ColMajor, int>;
  using Solver = Eigen::SimplicialLDLT<CholMatrix, Eigen::Lower, Eigen::NaturalOrdering<int>>;

  int n = 0;
  /// 置换：分解的是 ap = P·A·Pᵀ（下三角存储）。语义与 Eigen 内部一致（见 analyze 的注释）。
  Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic, int> P;
  Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic, int> Pinv;
  CholMatrix ap;                        ///< 置换后的矩阵（符号分解与数值分解都用它）
  Solver solver;                        ///< 分解本身仍是 Eigen 的 SimplicialLDLT
  const CholMatrix* factor = nullptr;   ///< 指向 solver 内部的因子（下三角 CSC，单位对角）
  Eigen::VectorXd dinv;                 ///< 1/D，预计算：省掉每次 solve 的 n 次除法
  Eigen::VectorXd work;                 ///< 复用缓冲（置换后的右端 → 解）
  bool ready = false;
};

EigenDirectSolver::EigenDirectSolver() : impl_(std::make_unique<Impl>()) {}
EigenDirectSolver::~EigenDirectSolver() = default;

void EigenDirectSolver::analyze(int n, const Eigen::SparseMatrix<Scalar>& patternMatrix) {
  const auto t0 = Clock::now();

  // 1) 自己算 AMD 置换。注意 Eigen 的约定：ordering 给出的是**逆**置换
  //    （见 Eigen/src/SparseCholesky/SimplicialCholesky.h 的 `ordering()`），
  //    内部随即取 `m_P = m_Pinv.inverse()`，这里照抄，保证分解与它完全一致
  //    （实测因子 nnz 逐位相同，例如 40×40 都是 66,027）。
  impl_->n = n;
  const Impl::CholMatrix sym = patternMatrix.selfadjointView<Eigen::Lower>();
  Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic, int> pinv;
  Eigen::AMDOrdering<int> ordering;
  ordering(sym, pinv);
  if (pinv.size() > 0) {
    impl_->P = pinv.inverse();
    impl_->Pinv = pinv;
  } else {
    impl_->P.resize(n);
    impl_->Pinv.resize(n);
  }

  // 2) 置换后的矩阵交给 NaturalOrdering 的 LDLᵀ 做符号分解（不再二次排序）。
  //    用真实矩阵的结构（不是理想块模式）—— 与 IGlobalSolver::analyze 的说明一致。
  impl_->ap.resize(n, n);
  impl_->ap.selfadjointView<Eigen::Lower>() =
      sym.selfadjointView<Eigen::Lower>().twistedBy(impl_->P);
  impl_->solver.analyzePattern(impl_->ap);

  impl_->ready = false;
  analyzed_ = true;
  stats_.analyzeCalls += 1;
  stats_.lastAnalyzeSeconds = secondsSince(t0);
  stats_.nnz = static_cast<long long>(patternMatrix.nonZeros());
}

void EigenDirectSolver::factorize(const Eigen::SparseMatrix<Scalar>& L) {
  const auto t0 = Clock::now();

  // 数值分解：结构与 analyzePattern 时**必须一致**（调用方保证）。仍然只填下三角。
  impl_->ap.selfadjointView<Eigen::Lower>() = L.selfadjointView<Eigen::Lower>().twistedBy(impl_->P);
  impl_->solver.factorize(impl_->ap);

  if (impl_->solver.info() == Eigen::Success) {
    impl_->factor = &impl_->solver.matrixL().nestedExpression();
    // 预存 1/D：Eigen 每次 solve 都要重算 n 次除法，这里一劳永逸。
    impl_->dinv = impl_->solver.vectorD().cwiseInverse();
    impl_->work.resize(impl_->n);
    impl_->ready = true;
  } else {
    impl_->ready = false;
    std::fprintf(stderr,
                 "[EigenDirectSolver] 数值分解失败（LDLᵀ info != Success）——矩阵可能不是正定的，"
                 "或者结构与 analyzePattern 时不一致。\n");
    std::fflush(stderr);
  }

  stats_.factorizeCalls += 1;
  stats_.lastFactorizeSeconds = secondsSince(t0);
  stats_.totalFactorizeSeconds += stats_.lastFactorizeSeconds;
  stats_.nnz = static_cast<long long>(L.nonZeros());
}

void EigenDirectSolver::solve(const Eigen::VectorXd& b, Eigen::VectorXd& x) {
  const auto t0 = Clock::now();
  const int n = impl_->n;
  const char* kPre = "[EigenDirectSolver::solve] 违反前置条件";

  if (!impl_->ready) {
    std::fprintf(stderr, "%s：尚未完成数值分解（或分解失败）。\n", kPre);
    std::fflush(stderr);
    std::abort();
  }
  if (b.size() != n) {
    std::fprintf(stderr, "%s：右端维度不匹配（b %lld，期望 %d）。\n", kPre,
                 static_cast<long long>(b.size()), n);
    std::fflush(stderr);
    std::abort();
  }
  // 输出向量按需要重新分配：原来的 `x = llt.solve(b)` 就会替你 resize，
  // 所以这里也必须允许传进来一个空/尺寸不对的 x（测试里就有这种用法，别收紧这条契约）。
  if (x.size() != n) x.resize(n);

  // 1) 置换到分解顺序：work = P·b
  impl_->work.noalias() = impl_->P * b;
  Scalar* w = impl_->work.data();

  // 2) 前代 L·y = work（单位下三角，scatter 形式，见文件头要点 1）
  for (int j = 0; j < n; ++j) {
    const Scalar yj = w[j];  // j 递增 ⇒ 此刻 w[j] 已累加完毕
    for (Impl::CholMatrix::InnerIterator it(*impl_->factor, j); it; ++it) {
      const int i = static_cast<int>(it.row());
      if (i > j) {
        w[i] -= it.value() * yj;
      }
    }
  }

  // 3) 对角缩放 z = y / D（乘预存的 1/D）
  for (int j = 0; j < n; ++j) {
    w[j] *= impl_->dinv[static_cast<Eigen::Index>(j)];
  }

  // 4) 回代 Lᵀ·x = z（扫 CSC 第 j 列取 i>j 的项；j 递减 ⇒ w[i]（i>j）已终值）
  for (int j = n - 1; j >= 0; --j) {
    Scalar s = w[j];
    for (Impl::CholMatrix::InnerIterator it(*impl_->factor, j); it; ++it) {
      const int i = static_cast<int>(it.row());
      if (i > j) {
        s -= it.value() * w[i];
      }
    }
    w[j] = s;
  }

  // 5) 置换回原顺序：x = P⁻¹·y
  x.noalias() = impl_->Pinv * impl_->work;

  stats_.solveCalls += 1;
  stats_.lastSolveSeconds = secondsSince(t0);
  stats_.totalSolveSeconds += stats_.lastSolveSeconds;
}

}  // namespace pd
