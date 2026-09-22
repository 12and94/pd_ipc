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

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

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

  // ---- x/y/z 三分量（L = Ã ⊗ I₃）的拆分信息 ----
  // 见 IGlobalSolver::parallelComponents 的完整说明。填不成 3 分量时 components 保持 1，
  // 调用方自动退回整趟求解。
  int components = 1;
  std::vector<int> origOfPerm;            ///< 置换后列号 → 原索引
  std::vector<std::vector<int>> compAsc;  ///< 每个分量的置换后列号（升序）
  std::vector<std::vector<int>> compDesc; ///< 同上（降序；回代按降序扫列）

  /// 在数值分解之后分析"能否按 3 个分量拆分"，并把列号表建好（一次性）。
  void buildComponents();
};

void EigenDirectSolver::Impl::buildComponents() {
  components = 1;
  origOfPerm.clear();
  compAsc.clear();
  compDesc.clear();
  if (n <= 0 || n % 3 != 0 || factor == nullptr) return;

  // 置换后列号 → 原索引：把 P 作用在 label 向量上一次拿到。
  // （Eigen 的置换索引语义容易记错，用数值探针是唯一不会写错的办法。）
  Eigen::VectorXd labels(n);
  for (int i = 0; i < n; ++i) labels[i] = static_cast<double>(i);
  const Eigen::VectorXd mapped = P * labels;
  origOfPerm.assign(static_cast<std::size_t>(n), 0);
  for (int k = 0; k < n; ++k) {
    origOfPerm[static_cast<std::size_t>(k)] = static_cast<int>(std::lround(mapped[k]));
  }
  {  // 自检：确实是一个置换（排序后应恰好是 0..n-1）
    std::vector<int> sorted = origOfPerm;
    std::sort(sorted.begin(), sorted.end());
    for (int i = 0; i < n; ++i) {
      if (sorted[static_cast<std::size_t>(i)] != i) return;
    }
  }

  // 组建分量列号表（分量 c = 原索引 ≡ c (mod 3) 的那些列），并检查每个分量等大。
  compAsc.assign(3, {});
  for (int k = 0; k < n; ++k) {
    const int c = origOfPerm[static_cast<std::size_t>(k)] % 3;
    compAsc[static_cast<std::size_t>(c)].push_back(k);
  }
  for (int c = 0; c < 3; ++c) {
    if (static_cast<int>(compAsc[static_cast<std::size_t>(c)].size()) != n / 3) return;
  }

  // 关键检查：填充不能跨分量。原因：L = Ã ⊗ I₃ ⇒ 它的图是 3 个互不相连的连通分量，
  // 而消元只会"在某个分量的节点之间"生成填充 ⇒ 因子必然按分量块对角。
  // 一旦将来引入各向异性/耦合 3×3 块的能量项，这里会命不中，components 自动保持 1。
  for (int j = 0; j < n; ++j) {
    const int cj = origOfPerm[static_cast<std::size_t>(j)] % 3;
    for (CholMatrix::InnerIterator it(*factor, j); it; ++it) {
      const int i = static_cast<int>(it.row());
      if (i == j) continue;
      if (origOfPerm[static_cast<std::size_t>(i)] % 3 != cj) return;
    }
  }

  compDesc = compAsc;
  for (auto& v : compDesc) std::reverse(v.begin(), v.end());
  components = 3;
}

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
    // 数值分解之后才知道因子结构 ⇒ 在这里做一次"能否按 x/y/z 三分量拆分"的分析
    // （一次性 O(nnz)；失败就保持 1 个分量，调用方退回整趟求解）。
    impl_->buildComponents();
  } else {
    impl_->ready = false;
    impl_->components = 1;
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

int EigenDirectSolver::parallelComponents() const { return impl_->components; }

void EigenDirectSolver::solveComponent(int c, const Eigen::VectorXd& b, Eigen::VectorXd& x) {
  // 结构不满足（不是 Ã ⊗ I₃）时不拆分：转回整趟求解，行为与改动前完全一致。
  if (impl_->components <= 1) {
    if (c != 0) {
      std::fprintf(stderr, "[EigenDirectSolver::solveComponent] 分量号越界：c=%d，但只有 1 个分量。\n",
                   c);
      std::fflush(stderr);
      std::abort();
    }
    solve(b, x);
    return;
  }

  const auto t0 = Clock::now();
  const int n = impl_->n;
  const char* kPre = "[EigenDirectSolver::solveComponent] 违反前置条件";

  if (!impl_->ready) {
    std::fprintf(stderr, "%s：尚未完成数值分解（或分解失败）。\n", kPre);
    std::fflush(stderr);
    std::abort();
  }
  if (c < 0 || c >= impl_->components) {
    std::fprintf(stderr, "%s：分量号越界（c=%d，共 %d 个分量）。\n", kPre, c, impl_->components);
    std::fflush(stderr);
    std::abort();
  }
  if (b.size() != n) {
    std::fprintf(stderr, "%s：右端维度不匹配（b %lld，期望 %d）。\n", kPre,
                 static_cast<long long>(b.size()), n);
    std::fflush(stderr);
    std::abort();
  }
  // 拆分路径**刻意不 resize**：多线程下并发 resize 同一个向量是竞争。
  // 调用方（stepOnce）在进并行区域之前就把 xSolution 备成了长度 n。
  if (x.size() != n) {
    std::fprintf(stderr, "%s：拆分求解要求 x 已经是长度 %d（实际 %lld）—— 请调用方先备好。\n", kPre,
                 n, static_cast<long long>(x.size()));
    std::fflush(stderr);
    std::abort();
  }

  const std::vector<int>& ac = impl_->compAsc[static_cast<std::size_t>(c)];
  const std::vector<int>& de = impl_->compDesc[static_cast<std::size_t>(c)];
  const int* o = impl_->origOfPerm.data();
  double* w = impl_->work.data();
  const Impl::CholMatrix& F = *impl_->factor;

  // 下面的三趟循环与 solve() 里逐行对应，唯一区别是"只遍历本分量的列"。
  // 因为分量之间内存完全不相交、且列顺序（升序 / 降序）与整趟一致，所以
  // 每个分量的算术序列与整趟**逐位相同** —— 这是本改造最重要的性质
  // （tests/primitives 里有"三分量之和 == 整趟"的逐位断言；pd_solvecomp 也钉住了它）。
  for (int k : ac) w[k] = b[o[k]];  // work = P·b（只取本分量）
  for (int j : ac) {                // 前代 L·y = work（scatter 形式，见文件头要点 1）
    const Scalar yj = w[j];
    for (Impl::CholMatrix::InnerIterator it(F, j); it; ++it) {
      const int i = static_cast<int>(it.row());
      if (i > j) {
        w[i] -= it.value() * yj;
      }
    }
  }
  for (int j : ac) {  // 对角缩放 z = y / D
    w[j] *= impl_->dinv[static_cast<Eigen::Index>(j)];
  }
  for (int j : de) {  // 回代 Lᵀ·x = z（扫 CSC 第 j 列取 i>j 的项）
    Scalar s = w[j];
    for (Impl::CholMatrix::InnerIterator it(F, j); it; ++it) {
      const int i = static_cast<int>(it.row());
      if (i > j) {
        s -= it.value() * w[i];
      }
    }
    w[j] = s;
  }
  for (int k : ac) x[o[k]] = w[k];  // x = P⁻¹·y（只写本分量）

  // 统计口径：`solveCalls` 只在分量 0 的那次调用上自增 ⇒ 仍然等于"全局步次数"，
  // 与改动前（每次 solve 一次）以及 tests/bench 的口径一致。
  // `totalSolveSeconds` 累加每个分量自己的耗时 ⇒ 仍近似等于"整趟串行的工作量"，
  // 因此"累计回代耗时 / 回代次数"这个比值仍与历史可比。
  // 注意 `lastSolveSeconds` 在拆分时是**单个分量的切片**（约为整趟的 1/3）。
  if (c == 0) stats_.solveCalls += 1;
  stats_.lastSolveSeconds = secondsSince(t0);
  stats_.totalSolveSeconds += stats_.lastSolveSeconds;
}

}  // namespace pd
