// _verify/solve_audit.cpp
// 全局步（Eigen 稀疏三角回代）的可行性审计 —— **只测量，不改生产代码**。
//
// 背景：Phase 1–3 之后，40×40 工况里"全局回代"占单子步的 76 %（docs/perf.md §5），
// 而它是**单线程**的（Eigen 的 SimplicialLDLT 不做并行回代），所以它是 Amdahl 天花板。
// 要判断值不值得动它，本工具按顺序回答四个可测的问题：
//
//   ① 这一趟是什么性质的开销？把"等效流量 ÷ 实测耗时"算成 GB/s，再与同规模的
//      **流式读写**实测带宽对照。远低于流式带宽 ⇒ 延迟受限（依赖链上的 load→use），
//      那么"少搬数据"没用，只有"缩短关键路径 / 提高访存并行度"才有用。
//   ② 因式结构里还剩多少并行度？从 L 的结构算**层集**（前代/回代各一遍）：
//      关键路径 = 层数 D，理想并行度 = 总工作量 / D。层数决定"至少要同步几次"。
//   ③ 换排序/换分解能拿到什么？AMD（现状）/ Natural / COLAMD，以及 LLT vs LDLT
//      （看 nnz(因子)、层结构、实测单次回代）。
//   ④ **决定性实验**：真的写一个层调度的并行前代/回代原型（自建置换，vcomp barrier 与
//      自旋 barrier 两种同步各测一遍），与 Eigen 的串行回代正面对比 —— 并且先用
//      "与 Eigen 解逐元素对齐"证明原型是对的。原型只在本工具里，不进生产。
//
// **2026-09-22 更新（本工具结论的后续）**：这个 Amdahl 天花板**已经被绕过去了**，但走的不是
// 本工具探索的"三角求解内部并行"那条路 —— 那条路（层调度）实测仍然不值得做。正确的问题是
// "这个系统里是不是本来就有独立的右端项"：本项目 `L = Ã ⊗ I₃`（每块都是标量 × I₃），
// 图是 3 个互不相连的连通分量 ⇒ 按 x/y/z 拆成 3 条**零同步**的独立链，实测单次回代快 2.0–2.8×、
// 主工况 −42.3 %、且物理输出逐位不变。见 `docs/perf.md` §9 与工具 `pd_solvecomp`
// （`_verify/solve_components.cpp`）。
//
// 用法：
//   pd_solveaudit [--grid 40 40] [--stiffness 1e4] [--solves 200] [--threads 4]
//                 [--order amd|natural|colamd] [--no-proto] [--proto-threads 1,2,4,8]
//
// 输出全部是中文的可复核数字，结论段落由数字直接推出（见最后"判定"一节）。

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <thread>
#include <vector>

#include <Eigen/Sparse>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "core/sim/Scene.h"
#include "core/assemble/Assembler.h"

using namespace pd;

namespace {

using Clock = std::chrono::steady_clock;
using SpMat = Eigen::SparseMatrix<Scalar>;

double secondsSince(const Clock::time_point& t0) {
  return std::chrono::duration<double>(Clock::now() - t0).count();
}

// ---------------------------------------------------------------- 层集分析
// 从**因子**（已按排序置换过的 L，下三角）算前代/回代的层集与每列工作量。
//
// 前代 L y = b：y[j] 依赖 y[i]（i < j 且 L(j,i) ≠ 0）。
//   按 i 递增扫：每条 (i → j>i) 的边把 level[j] 抬到 level[i]+1。
//   因为所有边都指向更大的下标，扫到 i 时 level[i] 已最终确定。
// 回代 Lᵀ z = y：z[j] 依赖 z[i]（i > j 且 L(i,j) ≠ 0）—— 同一批非零、方向相反，
//   所以按下标**递减**扫，取 level[i]+1 的最大值。
// 两者顺序执行（回代要用前代的结果），所以关键路径 = D_F + D_B。
struct LevelInfo {
  int levelsF = 0;
  int levelsB = 0;
  long long totalEntries = 0;
  long long upperEntries = 0;  // >0 说明因子按满阵存（正常情况下为 0）
  int maxWidthF = 0;
  int maxWidthB = 0;
  std::vector<long long> workF;          // 每层工作量（= 该层各列的"非对角非零"个数之和）
  std::vector<long long> workB;
  std::vector<long long> colNnz;         // 每列的非对角非零个数（= 该列的计算量）
  std::vector<std::vector<int>> nodesF;  // 每层的列号（供原型用）
  std::vector<std::vector<int>> nodesB;
};

LevelInfo computeLevels(const SpMat& F) {
  const int n = static_cast<int>(F.rows());
  LevelInfo out;
  std::vector<int> levelF(static_cast<std::size_t>(n), 0);
  std::vector<int> levelB(static_cast<std::size_t>(n), 0);
  out.colNnz.assign(static_cast<std::size_t>(n), 0);

  for (int i = 0; i < n; ++i) {
    const int lv = levelF[static_cast<std::size_t>(i)];
    for (SpMat::InnerIterator it(F, i); it; ++it) {
      const int j = static_cast<int>(it.row());
      out.totalEntries += 1;
      if (j == i) continue;  // 对角（LDLT 的单位对角）不算工作
      out.colNnz[static_cast<std::size_t>(i)] += 1;
      if (j > i) {
        int& target = levelF[static_cast<std::size_t>(j)];
        if (target < lv + 1) target = lv + 1;
      } else {
        out.upperEntries += 1;
      }
    }
  }
  for (int j = n - 1; j >= 0; --j) {
    int lv = 0;
    for (SpMat::InnerIterator it(F, j); it; ++it) {
      const int i = static_cast<int>(it.row());
      if (i > j) lv = std::max(lv, levelB[static_cast<std::size_t>(i)] + 1);
    }
    levelB[static_cast<std::size_t>(j)] = lv;
  }

  for (int v = 0; v < n; ++v) {
    out.levelsF = std::max(out.levelsF, levelF[static_cast<std::size_t>(v)] + 1);
    out.levelsB = std::max(out.levelsB, levelB[static_cast<std::size_t>(v)] + 1);
  }
  out.workF.assign(static_cast<std::size_t>(out.levelsF), 0);
  out.workB.assign(static_cast<std::size_t>(out.levelsB), 0);
  out.nodesF.resize(static_cast<std::size_t>(out.levelsF));
  out.nodesB.resize(static_cast<std::size_t>(out.levelsB));
  for (int v = 0; v < n; ++v) {
    const std::size_t lf = static_cast<std::size_t>(levelF[static_cast<std::size_t>(v)]);
    const std::size_t lb = static_cast<std::size_t>(levelB[static_cast<std::size_t>(v)]);
    out.workF[lf] += out.colNnz[static_cast<std::size_t>(v)];
    out.workB[lb] += out.colNnz[static_cast<std::size_t>(v)];
    out.nodesF[lf].push_back(v);
    out.nodesB[lb].push_back(v);
  }
  for (const auto w : out.workF) out.maxWidthF = static_cast<int>(std::max<long long>(out.maxWidthF, w));
  for (const auto w : out.workB) out.maxWidthB = static_cast<int>(std::max<long long>(out.maxWidthB, w));
  return out;
}

struct WorkStats {
  long long total = 0;
  long long max = 0;
  double mean = 0.0;
  double median = 0.0;
};

WorkStats summarize(const std::vector<long long>& work) {
  WorkStats s;
  if (work.empty()) return s;
  std::vector<long long> sorted = work;
  std::sort(sorted.begin(), sorted.end());
  for (const auto w : work) { s.total += w; s.max = std::max(s.max, w); }
  s.mean = static_cast<double>(s.total) / static_cast<double>(work.size());
  s.median = static_cast<double>(sorted[sorted.size() / 2]);
  return s;
}

// ---------------------------------------------------------------- barrier 实测
double measureOmpBarrierUs(int threads, int reps) {
#ifdef _OPENMP
  const auto t0 = Clock::now();
#pragma omp parallel num_threads(threads)
  {
    for (int r = 0; r < reps; ++r) {
#pragma omp barrier
    }
  }
  return secondsSince(t0) / reps * 1e6;
#else
  (void)threads;
  (void)reps;
  return 0.0;
#endif
}

/// 自旋 barrier（sense-reversing）：层调度每层同步一次，同步原语本身的成本直接决定成败，
/// 所以这里两种都测（vcomp 的 `omp barrier` vs 自己写的纯自旋）。
///
/// **纯自旋在高线程数下会卡死**（实测 18 线程 / 200×200 就一直不返回）：本机是 8P(含HT)+4E
/// 的异构布局，自旋线程会把正在干活的线程挤掉，形成活锁。所以这里带 yield 回退 ——
/// 这也说明"用一个便宜的自旋 barrier 替掉 omp barrier"并不是免费的午餐。
struct SpinBarrier {
  std::atomic<int> count{0};
  std::atomic<int> gen{0};

  void wait(int n) {
    if (n <= 1) return;
    const int g = gen.load(std::memory_order_acquire);
    if (count.fetch_add(1, std::memory_order_acq_rel) == n - 1) {
      count.store(0, std::memory_order_relaxed);
      gen.fetch_add(1, std::memory_order_release);
    } else {
      int spins = 0;
      while (gen.load(std::memory_order_acquire) == g) {
        if (++spins < 512) {
#if defined(_MSC_VER)
          _mm_pause();
#endif
        } else {
          std::this_thread::yield();  // 回退：避免活锁（见上面的说明）
        }
      }
    }
  }
};

double measureSpinBarrierUs(int threads, int reps) {
#ifdef _OPENMP
  SpinBarrier bar;
  const auto t0 = Clock::now();
#pragma omp parallel num_threads(threads)
  {
    const int nt = omp_get_num_threads();
    for (int r = 0; r < reps; ++r) bar.wait(nt);
  }
  return secondsSince(t0) / reps * 1e6;
#else
  (void)threads;
  (void)reps;
  return 0.0;
#endif
}

double measureStreamGBs(std::size_t n) {
  if (n < 1024) n = 1024;
  std::vector<double> a(n, 1.5), b(n, 2.5);
  const int reps = 20;
  for (std::size_t i = 0; i < n; ++i) b[i] = a[i] * 2.0 + b[i];
  const auto t0 = Clock::now();
  for (int r = 0; r < reps; ++r) {
    for (std::size_t i = 0; i < n; ++i) b[i] = a[i] * 2.0 + b[i];
  }
  const double sec = secondsSince(t0);
  return static_cast<double>(n) * 24.0 * reps / sec / 1e9;
}

// ---------------------------------------------------------------- Eigen 各配置
struct CaseResult {
  std::string name;
  bool ok = false;
  int n = 0;
  long long nnzMatrix = 0;
  long long nnzFactor = 0;
  double analyzeMs = 0.0;
  double factorMs = 0.0;
  double solveMinMs = 0.0;
  double solveMedMs = 0.0;
  LevelInfo levels;
};

template <typename SolverType>
CaseResult runCase(const std::string& name, const SpMat& L, int solves) {
  CaseResult r;
  r.name = name;
  r.n = static_cast<int>(L.rows());
  r.nnzMatrix = L.nonZeros();

  SolverType solver;
  const auto tAnalyze = Clock::now();
  solver.analyzePattern(L);
  r.analyzeMs = secondsSince(tAnalyze) * 1e3;
  const auto tFactor = Clock::now();
  solver.factorize(L);
  r.factorMs = secondsSince(tFactor) * 1e3;
  if (solver.info() != Eigen::Success) {
    std::printf("  %-10s 分解失败（info != Success）\n", name.c_str());
    return r;
  }
  r.ok = true;
  const auto& F = solver.matrixL().nestedExpression();
  r.nnzFactor = F.nonZeros();
  r.levels = computeLevels(F);

  Eigen::VectorXd b = Eigen::VectorXd::NullaryExpr(r.n, []() { return 1.0; });
  Eigen::VectorXd x;
  std::vector<double> samples;
  samples.reserve(static_cast<std::size_t>(solves));
  x = solver.solve(b);  // 预热（与生产代码同一写法）
  for (int s = 0; s < solves; ++s) {
    b[static_cast<int>(s % r.n)] += 0.25;
    const auto t0 = Clock::now();
    x = solver.solve(b);
    samples.push_back(secondsSince(t0) * 1e3);
  }
  std::sort(samples.begin(), samples.end());
  r.solveMinMs = samples.front();
  r.solveMedMs = samples[samples.size() / 2];
  return r;
}

void printCase(const CaseResult& r, int threads) {
  if (!r.ok) return;
  const WorkStats wf = summarize(r.levels.workF);
  const WorkStats wb = summarize(r.levels.workB);
  const int D = r.levels.levelsF + r.levels.levelsB;
  const long long totalWork = wf.total + wb.total;
  std::printf("\n--- %s ---\n", r.name.c_str());
  std::printf("  自由度 %d   组装矩阵 nnz %lld   因子 L nnz %lld（放大 %.2fx）\n", r.n, r.nnzMatrix,
              r.nnzFactor,
              static_cast<double>(r.nnzFactor) / static_cast<double>(std::max<long long>(1, r.nnzMatrix)));
  std::printf("  符号分解 %.2f ms   数值分解 %.2f ms   单次回代 min %.4f ms / 中位 %.4f ms\n",
              r.analyzeMs, r.factorMs, r.solveMinMs, r.solveMedMs);
  if (r.levels.upperEntries > 0) {
    std::printf("  注意：因子按满阵存储（上三角项 %lld），层集只按下三角算。\n", r.levels.upperEntries);
  }
  std::printf("  前代：层数 %d   每层工作量 中位 %.0f / 均值 %.0f / 最大 %d\n", r.levels.levelsF,
              wf.median, wf.mean, static_cast<int>(wf.max));
  std::printf("  回代：层数 %d   每层工作量 中位 %.0f / 均值 %.0f / 最大 %d\n", r.levels.levelsB,
              wb.median, wb.mean, static_cast<int>(wb.max));
  std::printf("  关键路径 %d 层；总工作量 %lld 项 ⇒ 理想并行度上限 %.1fx\n", D, totalWork,
              static_cast<double>(totalWork) / static_cast<double>(std::max(1, D)));
  const double bytes =
      static_cast<double>(r.nnzFactor) * 8.0 * 2.0 + static_cast<double>(r.n) * 8.0 * 4.0;
  std::printf("  单次回代等效流量 %.2f MB ⇒ 达成带宽 %.2f GB/s\n", bytes / 1e6,
              bytes / (r.solveMinMs * 1e-3) / 1e9);
  (void)threads;
}

// ---------------------------------------------------------------- ④ 原型：层调度并行回代
//
// 为什么要自建置换：Eigen 的 m_P 是 protected，拿不到就没办法自己做前代/回代。
// 这里复刻 Eigen 内部的流程（见 Eigen/src/SparseCholesky/SimplicialCholesky.h 的 ordering()）：
//     C  = A.selfadjointView<Lower>()            // 对称化
//     ordering(C, Pinv)                          // 注意 ordering 给的是**逆**置换
//     P  = Pinv.inverse()
//     Ap = P A Pᵀ                                // 只填下三角
// 然后把 Ap 交给 NaturalOrdering 的 LDLT（不再二次排序），于是**置换与因子都在自己手里**。
//
// 原型只做两件事（与 Eigen 的 solve 对应）：
//     前代 L y = b_perm        回代 Lᵀ x = z     中间 z = y / D
// 每层一次同步，层内按"工作量前缀和"手工静态分块（确定性 + 均衡），
// 同步用两种实现各测一遍：vcomp 的 `omp barrier` 与自旋 barrier。
struct ProtoFactor {
  bool ok = false;
  int n = 0;
  Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic, int> P;
  Eigen::SimplicialLDLT<SpMat, Eigen::Lower, Eigen::NaturalOrdering<int>> solver;
  const SpMat* F = nullptr;
  Eigen::VectorXd D;
  LevelInfo levels;
  double amdMs = 0.0;        // 自建 AMD 排序耗时
  double analyzeMs = 0.0;    // 符号分解
  double factorMs = 0.0;     // 数值分解
  double csrBuildMs = 0.0;   // 因子 CSR 副本的构建（一次性）
  // 因子的 **CSR（行）** 形式：前代必须用它（原因见 substForward）。
  std::vector<int> rowStart;
  std::vector<int> rowCol;
  std::vector<double> rowVal;
};

// 注意：不能用"返回 ProtoFactor"的写法 —— Eigen 的求解器不可拷贝/移动，
/// 具名返回值的 NRVO 也不是强制的（C++17 只保证纯右值省略），所以填引用参数。
bool buildProto(const SpMat& L, ProtoFactor& pf) {
  pf.n = static_cast<int>(L.rows());
  const SpMat C = L.selfadjointView<Eigen::Lower>();
  Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic, int> pinv;
  Eigen::AMDOrdering<int> ordering;
  const auto tAmd = Clock::now();
  ordering(C, pinv);
  pf.amdMs = secondsSince(tAmd) * 1e3;
  if (pinv.size() > 0) {
    pf.P = pinv.inverse();
  } else {
    pf.P.resize(pf.n);
  }
  SpMat Ap(pf.n, pf.n);
  Ap.selfadjointView<Eigen::Lower>() = C.selfadjointView<Eigen::Lower>().twistedBy(pf.P);
  const auto tAnalyze = Clock::now();
  pf.solver.analyzePattern(Ap);
  pf.analyzeMs = secondsSince(tAnalyze) * 1e3;
  const auto tFactor = Clock::now();
  pf.solver.factorize(Ap);
  pf.factorMs = secondsSince(tFactor) * 1e3;
  if (pf.solver.info() != Eigen::Success) return false;
  pf.F = &pf.solver.matrixL().nestedExpression();
  pf.D = pf.solver.vectorD();
  pf.levels = computeLevels(*pf.F);

  // 因子的 CSR 形式（严格下三角：行 j 存的是列 i<j 的 L(j,i)）。
  // 为什么非要有它：前代 y[j] = rhs[j] − Σ_{i<j} L(j,i)·y[i] 需要**按行**取数据，
  // 而因子是 CSC（按列存的）—— "第 j 列"里装的是 L(k,j)（k>j），不是我们要的 L(j,i)。
  // 用 CSC 直接写前代，`i < j` 一条都命不中，等于 y = rhs（这个坑真踩过：
  // 症状是"解出来差 5 %、残差 5.9 而 Eigen 是 1e-12"）。
  {
    const SpMat& F = *pf.F;
    const auto tCsr = Clock::now();
    std::vector<int> rowCount(static_cast<std::size_t>(pf.n), 0);
    for (int j = 0; j < F.outerSize(); ++j) {
      for (SpMat::InnerIterator it(F, j); it; ++it) rowCount[static_cast<std::size_t>(it.row())] += 1;
    }
    pf.rowStart.assign(static_cast<std::size_t>(pf.n) + 1, 0);
    for (int i = 0; i < pf.n; ++i) {
      pf.rowStart[static_cast<std::size_t>(i) + 1] =
          pf.rowStart[static_cast<std::size_t>(i)] + rowCount[static_cast<std::size_t>(i)];
    }
    const int total = pf.rowStart[static_cast<std::size_t>(pf.n)];
    pf.rowCol.assign(static_cast<std::size_t>(total), 0);
    pf.rowVal.assign(static_cast<std::size_t>(total), 0.0);
    std::vector<int> cursor(pf.rowStart.begin(), pf.rowStart.end() - 1);
    for (int j = 0; j < F.outerSize(); ++j) {
      for (SpMat::InnerIterator it(F, j); it; ++it) {
        const int r = static_cast<int>(it.row());
        const int p = cursor[static_cast<std::size_t>(r)]++;
        pf.rowCol[static_cast<std::size_t>(p)] = j;
        pf.rowVal[static_cast<std::size_t>(p)] = it.value();
      }
    }
    pf.csrBuildMs = secondsSince(tCsr) * 1e3;
  }
  pf.ok = true;
  return true;
}

/// 层内按工作量前缀和做确定性静态分块：返回 [begin,end) 下标范围。
void chunkRange(const std::vector<long long>& prefix, int nThreads, int tid, std::size_t& begin,
                std::size_t& end) {
  const std::size_t nNodes = prefix.size() - 1;
  const long long total = prefix.back();
  begin = 0;
  end = nNodes;
  if (nThreads <= 1 || total == 0) {
    if (tid != 0) { begin = end = 0; }
    return;
  }
  const long long lo = total * tid / nThreads;
  const long long hi = total * (tid + 1) / nThreads;
  begin = static_cast<std::size_t>(
      std::lower_bound(prefix.begin(), prefix.end(), lo) - prefix.begin());
  end = static_cast<std::size_t>(std::upper_bound(prefix.begin(), prefix.end(), hi) - prefix.begin());
  if (begin > nNodes) begin = nNodes;
  if (end > nNodes) end = nNodes;
  if (end < begin) end = begin;
}

struct LevelPlan {
  std::vector<std::vector<int>> nodes;
  std::vector<std::vector<long long>> prefix;  // 每层的工作量前缀和（含 0 与总和）
};

LevelPlan makePlan(const LevelInfo& lv, const std::vector<std::vector<int>>& nodesPerLevel) {
  LevelPlan p;
  p.nodes = nodesPerLevel;
  p.prefix.resize(p.nodes.size());
  for (std::size_t l = 0; l < p.nodes.size(); ++l) {
    auto& pre = p.prefix[l];
    pre.assign(p.nodes[l].size() + 1, 0);
    for (std::size_t k = 0; k < p.nodes[l].size(); ++k) {
      pre[k + 1] = pre[k] + lv.colNnz[static_cast<std::size_t>(p.nodes[l][k])];
    }
  }
  return p;
}

/// 前代 L y = rhs（L 是单位下三角，只存非对角项）：y[j] = rhs[j] − Σ_{i<j} L(j,i)·y[i]。
///
/// **必须按行（CSR）走**：这样每个节点只写 y[j]、只读 y[i]（i<j，属于更早的层），
/// 层内不冲突。若按列（CSC）走，第 j 列里是 L(k,j)（k>j），要写成
/// "算完 i 后把贡献推给 j" 就变成了 scatter（同层的多个 i 可能推同一个 j ⇒ 要归约/原子）。
///
/// 回代 Lᵀ x = rhs：x[j] = rhs[j] − Σ_{i>j} L(i,j)·x[i] —— 这里 L(i,j) 正好在**第 j 列**，
/// 所以回代用 CSC 反而天然不冲突。**同一个三角求解，两个方向各要一种存储**，这是
/// 并行稀疏三角求解最初的坑之一。
void substForward(const ProtoFactor& pf, const LevelPlan& plan, const double* rhs, double* out,
                  int threads, bool useOmpBarrier, bool indexOrder1T) {
  auto runLevel = [&](int tid, int nt, std::size_t l) {
    std::size_t b = 0, e = 0;
    chunkRange(plan.prefix[l], nt, tid, b, e);
    const auto& nodes = plan.nodes[l];
    for (std::size_t k = b; k < e; ++k) {
      const int j = nodes[k];
      double s = rhs[j];
      for (int p = pf.rowStart[static_cast<std::size_t>(j)];
           p < pf.rowStart[static_cast<std::size_t>(j) + 1]; ++p) {
        s -= pf.rowVal[static_cast<std::size_t>(p)] * out[pf.rowCol[static_cast<std::size_t>(p)]];
      }
      out[j] = s;
    }
  };
  const std::size_t nLev = plan.nodes.size();
  if (threads <= 1) {
    // 单线程：**按索引序**走，而不是层序。置换空间里"索引序"本身就是合法的消元顺序
    // （层集只是它的一个分组），而层序会让访问顺序在索引空间里跳来跳去、伤缓存 ——
    // 这个变体就是把"层遍历本身的局部性代价"单独量出来。
    if (indexOrder1T) {
      for (int j = 0; j < pf.n; ++j) {
        double s = rhs[j];
        for (int p = pf.rowStart[static_cast<std::size_t>(j)];
             p < pf.rowStart[static_cast<std::size_t>(j) + 1]; ++p) {
          s -= pf.rowVal[static_cast<std::size_t>(p)] * out[pf.rowCol[static_cast<std::size_t>(p)]];
        }
        out[j] = s;
      }
      return;
    }
    for (std::size_t l = 0; l < nLev; ++l) runLevel(0, 1, l);
    return;
  }
#ifdef _OPENMP
  SpinBarrier spin;
#pragma omp parallel num_threads(threads)
  {
    const int tid = omp_get_thread_num();
    const int nt = omp_get_num_threads();
    for (std::size_t l = 0; l < nLev; ++l) {
      runLevel(tid, nt, l);
      if (useOmpBarrier) {
#pragma omp barrier
      } else {
        spin.wait(nt);
      }
    }
  }
#else
  (void)useOmpBarrier;
  for (std::size_t l = 0; l < nLev; ++l) runLevel(0, 1, l);
#endif
}

/// 回代 Lᵀ x = rhs：x[j] = rhs[j] − Σ_{i>j} L(i,j)·x[i]（CSC：第 j 列的 i>j 项）。
void substBackward(const ProtoFactor& pf, const LevelPlan& plan, const double* rhs, double* out,
                   int threads, bool useOmpBarrier, bool indexOrder1T) {
  const SpMat& F = *pf.F;
  auto runLevel = [&](int tid, int nt, std::size_t l) {
    std::size_t b = 0, e = 0;
    chunkRange(plan.prefix[l], nt, tid, b, e);
    const auto& nodes = plan.nodes[l];
    for (std::size_t k = b; k < e; ++k) {
      const int j = nodes[k];
      double s = rhs[j];
      for (SpMat::InnerIterator it(F, j); it; ++it) {
        const int i = static_cast<int>(it.row());
        if (i > j) s -= it.value() * out[i];
      }
      out[j] = s;
    }
  };
  const std::size_t nLev = plan.nodes.size();
  if (threads <= 1) {
    // 单线程：按索引序（**逆序**，回代要从后往前），同样是为了局部性。
    if (indexOrder1T) {
      for (int j = pf.n - 1; j >= 0; --j) {
        double s = rhs[j];
        for (SpMat::InnerIterator it(F, j); it; ++it) {
          const int i = static_cast<int>(it.row());
          if (i > j) s -= it.value() * out[i];
        }
        out[j] = s;
      }
      return;
    }
    for (std::size_t l = 0; l < nLev; ++l) runLevel(0, 1, l);
    return;
  }
#ifdef _OPENMP
  SpinBarrier spin;
#pragma omp parallel num_threads(threads)
  {
    const int tid = omp_get_thread_num();
    const int nt = omp_get_num_threads();
    for (std::size_t l = 0; l < nLev; ++l) {
      runLevel(tid, nt, l);
      if (useOmpBarrier) {
#pragma omp barrier
      } else {
        spin.wait(nt);
      }
    }
  }
#else
  (void)useOmpBarrier;
  for (std::size_t l = 0; l < nLev; ++l) runLevel(0, 1, l);
#endif
}

struct ProtoTiming {
  int threads = 0;
  double ompMinMs = 0.0;
  double spinMinMs = 0.0;
};

/// 置换空间内的两步三角求解（不含置换本身）：L y = bp → z = y/D → Lᵀ x = z。
void protoPermutedSolve(const ProtoFactor& pf, const Eigen::VectorXd& bp, Eigen::VectorXd& xp,
                        int threads, bool useOmpBarrier, bool indexOrder1T, const LevelPlan& planF,
                        const LevelPlan& planB) {
  std::vector<double> y(static_cast<std::size_t>(pf.n), 0.0);
  substForward(pf, planF, bp.data(), y.data(), threads, useOmpBarrier, indexOrder1T);
  for (int i = 0; i < pf.n; ++i) y[static_cast<std::size_t>(i)] /= pf.D[i];  // z = y / D
  std::vector<double> xs(static_cast<std::size_t>(pf.n), 0.0);
  substBackward(pf, planB, y.data(), xs.data(), threads, useOmpBarrier, indexOrder1T);
  xp = Eigen::Map<Eigen::VectorXd>(xs.data(), pf.n);
}

/// 原型求解一次（含置换与 D 缩放），返回耗时（ms）。
double protoSolve(const ProtoFactor& pf, const Eigen::VectorXd& b, Eigen::VectorXd& x, int threads,
                  bool useOmpBarrier, const LevelPlan& planF, const LevelPlan& planB,
                  bool indexOrder1T = false) {
  const auto t0 = Clock::now();
  const Eigen::VectorXd bp = pf.P * b;
  Eigen::VectorXd xp;
  protoPermutedSolve(pf, bp, xp, threads, useOmpBarrier, indexOrder1T, planF, planB);
  x = pf.P.inverse() * xp;
  return secondsSince(t0) * 1e3;
}

void runProto(const SpMat& L, const Eigen::VectorXd& bRef, const Eigen::VectorXd& xRef,
              const CaseResult& eigenCase, const std::vector<int>& threadList) {
  using LdltAmd = Eigen::SimplicialLDLT<SpMat, Eigen::Lower, Eigen::AMDOrdering<int>>;
  ProtoFactor pf;
  if (!buildProto(L, pf)) {
    std::printf("\n④ 原型未能建立（分解失败）。\n");
    return;
  }
  std::printf("\n=== ④ 原型：层调度并行回代（自建置换 + 自己的前代/回代）===\n");
  std::printf("  因子 L nnz %lld（Eigen 的 AMD 配置：%lld，两者应接近）\n",
              static_cast<long long>(pf.F->nonZeros()), eigenCase.nnzFactor);
  std::printf("  一次性成本：自建 AMD %.2f ms + 符号分解 %.2f ms + 数值分解 %.2f ms + 因子 CSR 副本 "
              "%.2f ms = %.2f ms\n    （Eigen 的 analyzePattern + factorize = %.2f + %.2f = %.2f ms）\n",
              pf.amdMs, pf.analyzeMs, pf.factorMs, pf.csrBuildMs,
              pf.amdMs + pf.analyzeMs + pf.factorMs + pf.csrBuildMs, eigenCase.analyzeMs,
              eigenCase.factorMs, eigenCase.analyzeMs + eigenCase.factorMs);
  std::printf("  层数 前代 %d / 回代 %d；总工作量 %lld 项\n", pf.levels.levelsF, pf.levels.levelsB,
              summarize(pf.levels.workF).total + summarize(pf.levels.workB).total);

  const LevelPlan planF = makePlan(pf.levels, pf.levels.nodesF);
  const LevelPlan planB = makePlan(pf.levels, pf.levels.nodesB);

  // 正确性自证分两步，出问题才能定位到"三角求解实现"还是"置换应用"：
  //   (a) 置换空间内：与"直接用 Eigen 解置换后的矩阵 Ap"对比（同一 b_perm）；
  //   (b) 原空间：与 Eigen 的 AMD 解对比。
  Eigen::VectorXd bp = pf.P * bRef;
  Eigen::VectorXd xpMine;
  protoPermutedSolve(pf, bp, xpMine, 1, true, false, planF, planB);
  Eigen::VectorXd xpRef = pf.solver.solve(bp);  // 注意 pf.solver 就建在 Ap 上、Natural 序（无置换）
  const double diffPerm = (xpMine - xpRef).cwiseAbs().maxCoeff();
  const double scalePerm = xpRef.cwiseAbs().maxCoeff();
  std::printf("  正确性(a) 置换空间内 原型 vs Eigen(同一 Ap)：最大差 %.3e（相对 %.3e）⇒ %s\n", diffPerm,
              diffPerm / std::max(1e-300, scalePerm),
              (diffPerm / std::max(1e-300, scalePerm) < 1e-9) ? "一致" : "**不一致**");

  // 诊断（只在中小规模做）：谁是对的？注意**不能**用"只存了下三角"的矩阵去做乘法
  // （那样乘积不是对称矩阵），所以这里用真正的乘积 P L Pᵀ（L 由组装器填满两个三角）。
  // 大网格上 Eigen 的通用稀疏乘法会非常慢（实测 200×200 会让本工具超时），所以设个上限。
  if (pf.n <= 30000) {
    const SpMat ApFull = pf.P * L * pf.P.inverse();
    const SpMat ApFullAlt = pf.P.inverse() * L * pf.P;
    const double rMine = (ApFull * xpMine - bp).cwiseAbs().maxCoeff();
    const double rRef = (ApFull * xpRef - bp).cwiseAbs().maxCoeff();
    const double rMineAlt = (ApFullAlt * xpMine - bp).cwiseAbs().maxCoeff();
    const double rRefAlt = (ApFullAlt * xpRef - bp).cwiseAbs().maxCoeff();
    int zeroDiag = 0;
    for (int i = 0; i < L.rows(); ++i) {
      if (L.coeff(i, i) == 0.0) zeroDiag += 1;
    }
    std::printf("  诊断：‖P L Pᵀ·x − b_perm‖∞  原型 %.3e / Eigen %.3e   ⇒ 原型%s\n", rMine, rRef,
                (rMine <= rRef * 1.001) ? "至少不差" : "**更差，实现有问题**");
    std::printf("        若把方向当成 Pᵀ L P：原型 %.3e / Eigen %.3e（用来判断置换方向）\n", rMineAlt,
                rRefAlt);
    std::printf("        零对角 %d 个；有限性 原型=%d Eigen=%d\n", zeroDiag,
                static_cast<int>(xpMine.allFinite()), static_cast<int>(xpRef.allFinite()));
  } else {
    std::printf("  诊断：规模较大（自由度 %d > 30000），跳过置换乘积诊断（Eigen 通用稀疏乘法太慢）。\n",
                pf.n);
  }

  Eigen::VectorXd xSer;
  protoSolve(pf, bRef, xSer, 1, true, planF, planB);
  const double scale = xRef.cwiseAbs().maxCoeff();
  const double diff = (xSer - xRef).cwiseAbs().maxCoeff();
  std::printf("  正确性(b) 原空间  原型 vs Eigen(AMD)：最大差 %.3e（相对 %.3e）⇒ %s\n", diff,
              diff / std::max(1e-300, scale), (diff / scale < 1e-9) ? "一致" : "**不一致**");

  const int solves = 60;

  // ---------- 计时（纪律见 docs/perf.md §0）----------
  //
  // 两条必须遵守的规矩，这里都踩过：
  //   ① **交替测量**：A/B 必须背靠背轮换、各取最小值，否则机器漂移会被算进差值；
  //   ② **串行变体必须排在所有多线程区域之前**：vcomp 的线程池在区域结束后还会自旋一阵，
  //      紧跟 18 线程跑"单线程"变体会被邻居抢 CPU —— 实测同一个索引序变体出现过
  //      5.28 ms 与 10.29 ms 两个读数（差一倍），全是这个原因。
  {
    std::printf("\n  【单线程三方对照】在跑任何多线程区域**之前**、交替测量（各 %d 次取最小）：\n", solves);
    LdltAmd prod;  // 与生产同配置：LDLT + AMD（含 Eigen 自己的内部置换）
    prod.analyzePattern(L);
    prod.factorize(L);
    Eigen::VectorXd xEig;
    std::vector<double> tEig, tIdx, tLvl;
    tEig.reserve(solves); tIdx.reserve(solves); tLvl.reserve(solves);
    Eigen::VectorXd xTmp;
    xEig = prod.solve(bRef);  // 预热
    for (int rep = 0; rep < solves; ++rep) {
      const auto t0 = Clock::now();
      xEig = prod.solve(bRef);
      tEig.push_back(secondsSince(t0) * 1e3);
      tIdx.push_back(protoSolve(pf, bRef, xTmp, 1, false, planF, planB, true));
      tLvl.push_back(protoSolve(pf, bRef, xTmp, 1, false, planF, planB, false));
    }
    std::sort(tEig.begin(), tEig.end());
    std::sort(tIdx.begin(), tIdx.end());
    std::sort(tLvl.begin(), tLvl.end());
    std::printf("    Eigen solve（含置换）   %.4f ms\n", tEig.front());
    std::printf("    手写回代 · 索引序(1T)   %.4f ms   ⇒ 相对 Eigen %.2fx\n", tIdx.front(),
                tEig.front() / tIdx.front());
    std::printf("    手写回代 · 层序(1T)     %.4f ms   ⇒ 层遍历本身慢 %.0f %%\n", tLvl.front(),
                100.0 * (tLvl.front() / tIdx.front() - 1.0));
    std::fflush(stdout);
  }

  // ---------- 并行对照（此时线程池已经热了，所以 Eigen 也在同一个循环里交替测）----------
  std::printf("\n  【并行对照】每个 rep 内交替测 Eigen 与各线程数（各取最小值）：\n");
  std::printf("  %-8s %-14s %-14s %-14s %-10s %s\n", "线程", "vcomp barrier", "自旋 barrier", "串行(Eigen)",
              "加速(vcomp)", "加速(自旋)");
  {
    LdltAmd prod;
    prod.analyzePattern(L);
    prod.factorize(L);
    Eigen::VectorXd xEig;
    std::vector<double> tEig;
    std::vector<std::vector<double>> tOmp(threadList.size());
    std::vector<std::vector<double>> tSpin(threadList.size());
    Eigen::VectorXd xTmp;
    xEig = prod.solve(bRef);
    for (int rep = 0; rep < solves; ++rep) {
      const auto t0 = Clock::now();
      xEig = prod.solve(bRef);
      tEig.push_back(secondsSince(t0) * 1e3);
      for (std::size_t k = 0; k < threadList.size(); ++k) {
        const int t = threadList[k];
        tOmp[k].push_back(protoSolve(pf, bRef, xTmp, t, true, planF, planB));
        tSpin[k].push_back(protoSolve(pf, bRef, xTmp, t, false, planF, planB));
      }
    }
    std::sort(tEig.begin(), tEig.end());
    const double eigenMs = tEig.front();
    for (std::size_t k = 0; k < threadList.size(); ++k) {
      std::sort(tOmp[k].begin(), tOmp[k].end());
      std::sort(tSpin[k].begin(), tSpin[k].end());
      std::printf("  %-8d %-14.4f %-14.4f %-14.4f %-10.2fx %.2fx\n", threadList[k], tOmp[k].front(),
                  tSpin[k].front(), eigenMs, eigenMs / tOmp[k].front(), eigenMs / tSpin[k].front());
      std::fflush(stdout);
    }
  }
  std::printf("\n  同步成本实测（µs/次，空载；真实负载下会高得多，见 docs/solver-feasibility.md §4）：\n"
              "  %-8s %-16s %s\n", "线程", "vcomp `omp barrier`", "自旋 barrier");
  for (int t : threadList) {
    std::printf("  %-8d %-16.2f %.2f\n", t, measureOmpBarrierUs(t, 50000), measureSpinBarrierUs(t, 50000));
  }
}

void printUsage() {
  std::printf(
      "用法：pd_solveaudit [--grid N M] [--stiffness K] [--solves S] [--threads T]\n"
      "                    [--order amd|natural|colamd] [--no-proto] [--proto-threads 1,2,4,8]\n");
}

}  // namespace

int main(int argc, char** argv) {
  int nx = 40;
  int ny = 40;
  double stiffness = 1.0e4;
  int solves = 200;
  int threads = 4;
  std::string onlyOrder;
  bool proto = true;
  std::vector<int> protoThreads{1, 2, 4, 8};

  for (int i = 1; i < argc; ++i) {
    const std::string a = argv[i];
    const bool hasNext = (i + 1 < argc);
    if (a == "--grid" && i + 2 < argc) {
      nx = std::atoi(argv[++i]);
      ny = std::atoi(argv[++i]);
    } else if (a == "--stiffness" && hasNext) {
      stiffness = std::atof(argv[++i]);
    } else if (a == "--solves" && hasNext) {
      solves = std::atoi(argv[++i]);
    } else if (a == "--threads" && hasNext) {
      threads = std::atoi(argv[++i]);
    } else if (a == "--order" && hasNext) {
      onlyOrder = argv[++i];
    } else if (a == "--no-proto") {
      proto = false;
    } else if (a == "--proto-threads" && hasNext) {
      protoThreads.clear();
      std::string list = argv[++i];
      std::size_t pos = 0;
      while (pos <= list.size()) {
        const std::size_t comma = list.find(',', pos);
        const std::string item = list.substr(pos, comma == std::string::npos ? std::string::npos : comma - pos);
        if (!item.empty()) protoThreads.push_back(std::atoi(item.c_str()));
        if (comma == std::string::npos) break;
        pos = comma + 1;
      }
    } else if (a == "--help" || a == "-h") {
      printUsage();
      return 0;
    } else {
      std::printf("未知参数：%s\n", a.c_str());
      printUsage();
      return 2;
    }
  }

  std::printf("=== 全局步求解器可行性审计 ===\n");
  std::printf("网格 %dx%d   刚度 %.3g   每例回代次数 %d   线程（barrier 实测用）%d\n", nx, ny,
              stiffness, solves, threads);

  SceneConfig cfg;
  cfg.gridNx = nx;
  cfg.gridNy = ny;
  cfg.stiffness = stiffness;
  cfg.maxIterations = 1;  // 本工具只关心 L 与回代，不跑迭代
  SimContext ctx = makeScene(cfg);

  SpMat L;
  assembleLeftHandSide(ctx.mesh, cfg.dt, cfg.velocityDamping, L);
  std::printf("顶点 %d  约束 %d  自由度 %d  组装矩阵 nnz %lld\n", ctx.mesh.vertexCount(),
              ctx.mesh.edgeCount(), 3 * ctx.mesh.vertexCount(),
              static_cast<long long>(L.nonZeros()));

  const double streamGBs = measureStreamGBs(static_cast<std::size_t>(L.nonZeros()));
  std::printf("参考：同规模（%.0f 万项）流式读写实测 %.2f GB/s\n",
              static_cast<double>(L.nonZeros()) / 1e4, streamGBs);

  using LdltAmd = Eigen::SimplicialLDLT<SpMat, Eigen::Lower, Eigen::AMDOrdering<int>>;
  using LdltNat = Eigen::SimplicialLDLT<SpMat, Eigen::Lower, Eigen::NaturalOrdering<int>>;
  using LdltColamd = Eigen::SimplicialLDLT<SpMat, Eigen::Lower, Eigen::COLAMDOrdering<int>>;
  using LltAmd = Eigen::SimplicialLLT<SpMat, Eigen::Lower, Eigen::AMDOrdering<int>>;

  const CaseResult rAmd = runCase<LdltAmd>("LDLT + AMD（现状）", L, solves);
  printCase(rAmd, threads);

  const bool wantAll = onlyOrder.empty();
  if (wantAll || onlyOrder == "natural") {
    printCase(runCase<LdltNat>("LDLT + Natural", L, std::max(20, solves / 2)), threads);
  }
  if (wantAll || onlyOrder == "colamd") {
    printCase(runCase<LdltColamd>("LDLT + COLAMD", L, std::max(20, solves / 2)), threads);
  }
  if (wantAll) {
    printCase(runCase<LltAmd>("LLT + AMD", L, std::max(20, solves / 2)), threads);
  }

  if (proto && rAmd.ok) {
    Eigen::VectorXd b = Eigen::VectorXd::NullaryExpr(rAmd.n, []() { return 1.0; });
    Eigen::VectorXd xRef;
    {
      LdltAmd s;
      s.analyzePattern(L);
      s.factorize(L);
      xRef = s.solve(b);
    }
    runProto(L, b, xRef, rAmd, protoThreads);
  }

  // ---------------------------------------------------------------- 判定
  std::printf("\n=== 判定 ===\n");
  if (rAmd.ok) {
    const double bytes =
        static_cast<double>(rAmd.nnzFactor) * 8.0 * 2.0 + static_cast<double>(rAmd.n) * 8.0 * 4.0;
    const double gbs = bytes / (rAmd.solveMinMs * 1e-3) / 1e9;
    const double ratio = gbs / std::max(1e-9, streamGBs);
    const int D = rAmd.levels.levelsF + rAmd.levels.levelsB;
    const long long totalWork = summarize(rAmd.levels.workF).total + summarize(rAmd.levels.workB).total;
    std::printf("  ① 性质：达成带宽 %.2f GB/s，是流式参考的 %.0f%% ⇒ %s\n", gbs, ratio * 100.0,
                ratio < 0.5 ? "**延迟受限**（依赖链上 load→use，不是搬不动数据）"
                            : "接近带宽受限（搬数据本身就是成本）");
    std::printf("  ② 结构并行度：关键路径 %d 层，理想并行度上限 %.1fx\n", D,
                static_cast<double>(totalWork) / static_cast<double>(std::max(1, D)));
    std::printf("  ③ 层数 ⇒ 至少 %d 次同步；本机实测 vcomp %.2f µs/次、自旋 %.2f µs/次（4 线程）\n", D,
                measureOmpBarrierUs(4, 50000), measureSpinBarrierUs(4, 50000));
    std::printf("     同步总成本：vcomp %.3f ms、自旋 %.3f ms；串行回代 %.4f ms\n",
                D * measureOmpBarrierUs(4, 50000) / 1e3, D * measureSpinBarrierUs(4, 50000) / 1e3,
                rAmd.solveMinMs);
  }
  return 0;
}
