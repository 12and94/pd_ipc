// _verify/solve_components.cpp
//
// 目的：验证并测量「按 x/y/z 三分量拆分全局步回代」的可行性与收益。
//
// 背景（2026-09-22 新发现，此前未被记录）：
//   L 的组装保证**每一块都是标量 × I₃**：
//     · 惯性项 (m_v/h²) 出现在 (3v+d, 3v+d) —— d = 0,1,2 同一个值；
//     · 距离约束块 ±κ 出现在 (3a+d, 3b+d) —— d 不变（对角块 (3a+d,3a+d) 同理）；
//     · pin 行 (3v+d, 3v+d) = 1。
//   （见 core/assemble/Assembler.cpp L93-117 与 core/energy/DistanceTerm.cpp L275-290。）
//   即 L = Ã ⊗ I₃ —— 3n 个方程其实是 3 个**互不相连**的标量系统（x/y/z 完全解耦）。
//   消元的填充不可能跨越连通分量 ⇒ 因子是 3 份结构同构的标量因子，
//   前代 / 对角缩放 / 回代 / 两次置换都能按分量切开，而且每个分量的**算术序列与现状
//   逐位相同**（列顺序、i>j 筛选、累加顺序都不变）⇒ 零同步、零竞争、结果逐位不变。
//
//   为什么它重要：solver-feasibility.md §2 实测回代的**达成带宽只有流式读写的 9–13 %**，
//   也就是说瓶颈是"未决访存请求数（MLP）不足"，而不是搬运能力。3 条互不相干的依赖链
//   分给 3 条线程，恰好把 MLP 提高约 3 倍，而**不需要三角求解内部的任何并行化**
//   （层调度那条路已被实测否掉：层内工作量与同步成本同阶）。
//
// 本工具做三件事（全部只测量，不改生产代码）：
//   ① 结构自检：把 P 作用在 label 向量上得到"置换后列号 → 原索引"的映射，逐列断言
//      "第 j 列里所有 i≠j 的项都与 j 同分量"（违例数必须为 0），并报告每个分量的
//      列数与因子 nnz（应当三份相等）。
//   ② 正确性：V1（按分量拆开做）与 V0（现状写法做整趟）的解**逐位比较**。
//   ③ 微基准：在**同一个并行区域内**逐轮交替测量
//        A = omp master 里跑整趟（等价于现状的 omp single 形状）
//        B = omp for 三个分量（生产要用的形状；计时含"等最慢者"的墙钟）
//        C = omp master 里按分量顺序串行跑三趟（只有循环结构变了，没有并行）
//        N = 空任务（量化这一对 barrier 的固定成本）
//      各取最小值。A vs C 分离"循环结构收益"，B vs A 给出"并行收益"。
//
// 用法：
//   pd_solvecomp [--grid 40 40] [--stiffness 1e4] [--reps 200] [--threads 4]

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

#include <Eigen/Sparse>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "core/assemble/Assembler.h"
#include "core/sim/Scene.h"

using namespace pd;

namespace {

using Clock = std::chrono::steady_clock;
using SpMat = Eigen::SparseMatrix<Scalar>;
using Solver = Eigen::SimplicialLDLT<SpMat, Eigen::Lower, Eigen::NaturalOrdering<int>>;

double secondsSince(const Clock::time_point& t0) {
  return std::chrono::duration<double>(Clock::now() - t0).count();
}

constexpr int kComponents = 3;

struct FactorSet {
  bool ok = false;
  int n = 0;
  Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic, int> P;
  Solver solver;
  const SpMat* F = nullptr;                 ///< 因子（下三角 CSC，单位对角）
  Eigen::VectorXd dinv;                     ///< 1/D（与生产代码一样预存）
  std::vector<int> origOfPerm;              ///< 置换后列号 k → 原索引
  std::vector<std::vector<int>> asc;        ///< 每个分量的列号（升序）
  std::vector<std::vector<int>> desc;       ///< 每个分量的列号（降序）
  std::vector<long long> nnzByComp;         ///< 每个分量占据的因子非零个数
  long long crossLinks = 0;                 ///< 跨分量的 i≠j 项（必须为 0）
};

/// 与生产代码 EigenDirectSolver::analyze/factorize 同构：自己调 AMD 拿置换，
/// 用 NaturalOrdering 分解置换后的矩阵，只取矩阵真实结构。
bool buildFactor(const SpMat& L, FactorSet& fs) {
  fs.n = static_cast<int>(L.rows());
  const SpMat sym = L.selfadjointView<Eigen::Lower>();

  Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic, int> pinv;
  Eigen::AMDOrdering<int> ordering;
  ordering(sym, pinv);
  if (pinv.size() > 0) {
    fs.P = pinv.inverse();
  } else {
    fs.P.resize(fs.n);
  }

  SpMat ap(fs.n, fs.n);
  ap.selfadjointView<Eigen::Lower>() = sym.selfadjointView<Eigen::Lower>().twistedBy(fs.P);
  fs.solver.analyzePattern(ap);
  fs.solver.factorize(ap);
  if (fs.solver.info() != Eigen::Success) return false;

  fs.F = &fs.solver.matrixL().nestedExpression();
  fs.dinv = fs.solver.vectorD().cwiseInverse();

  // ---- "置换后列号 → 原索引" 的映射：把 P 作用在 label 向量上一次拿到 ----
  // Eigen 的索引语义容易记错，用数值探针是唯一不会写错的办法。
  Eigen::VectorXd labels(fs.n);
  for (int i = 0; i < fs.n; ++i) labels[i] = static_cast<double>(i);
  const Eigen::VectorXd mapped = fs.P * labels;
  fs.origOfPerm.assign(static_cast<std::size_t>(fs.n), 0);
  for (int k = 0; k < fs.n; ++k) {
    fs.origOfPerm[static_cast<std::size_t>(k)] = static_cast<int>(std::lround(mapped[k]));
  }

  // 自检：确实是一个置换（排序后等于 0..n-1）
  {
    std::vector<int> sorted = fs.origOfPerm;
    std::sort(sorted.begin(), sorted.end());
    for (int i = 0; i < fs.n; ++i) {
      if (sorted[static_cast<std::size_t>(i)] != i) return false;
    }
  }

  // ---- 结构性检查 + 分量列号表 ----
  fs.asc.assign(kComponents, {});
  fs.desc.assign(kComponents, {});
  fs.nnzByComp.assign(kComponents, 0);
  for (int k = 0; k < fs.n; ++k) {
    const int c = fs.origOfPerm[static_cast<std::size_t>(k)] % kComponents;
    fs.asc[static_cast<std::size_t>(c)].push_back(k);
  }
  for (int c = 0; c < kComponents; ++c) {
    fs.desc[static_cast<std::size_t>(c)] = fs.asc[static_cast<std::size_t>(c)];
    std::reverse(fs.desc[static_cast<std::size_t>(c)].begin(), fs.desc[static_cast<std::size_t>(c)].end());
  }
  for (int j = 0; j < fs.n; ++j) {
    const int cj = fs.origOfPerm[static_cast<std::size_t>(j)] % kComponents;
    for (SpMat::InnerIterator it(*fs.F, j); it; ++it) {
      const int i = static_cast<int>(it.row());
      fs.nnzByComp[static_cast<std::size_t>(cj)] += 1;
      if (i == j) continue;
      if (fs.origOfPerm[static_cast<std::size_t>(i)] % kComponents != cj) fs.crossLinks += 1;
    }
  }
  fs.ok = true;
  return true;
}

/// V0：与生产代码 EigenDirectSolver::solve 逐行等价的写法（单线程跑完整 3n 系统）。
/// w 是长度为 n 的复用缓冲（调用方提供）。
inline void solveFull(const FactorSet& fs, double* w, const double* b, double* x) {
  const int n = fs.n;
  const int* o = fs.origOfPerm.data();
  const SpMat& F = *fs.F;

  for (int k = 0; k < n; ++k) w[k] = b[o[k]];              // work = P·b
  for (int j = 0; j < n; ++j) {                            // 前代（scatter 形式）
    const double yj = w[j];
    for (SpMat::InnerIterator it(F, j); it; ++it) {
      const int i = static_cast<int>(it.row());
      if (i > j) w[i] -= it.value() * yj;
    }
  }
  for (int j = 0; j < n; ++j) w[j] *= fs.dinv[j];          // z = y / D
  for (int j = n - 1; j >= 0; --j) {                       // 回代（扫 CSC 列）
    double s = w[j];
    for (SpMat::InnerIterator it(F, j); it; ++it) {
      const int i = static_cast<int>(it.row());
      if (i > j) s -= it.value() * w[i];
    }
    w[j] = s;
  }
  for (int k = 0; k < n; ++k) x[o[k]] = w[k];              // x = P⁻¹·y
}

/// V1：只处理第 c 个分量的列（其余与 V0 完全相同 ⇒ 每个分量的算术逐位不变）。
/// 分量之间的内存完全不相交 ⇒ 3 条线程并行时零竞争、零同步。
inline void solveOneComponent(const FactorSet& fs, int c, double* w, const double* b, double* x) {
  const int* o = fs.origOfPerm.data();
  const SpMat& F = *fs.F;
  const std::vector<int>& ac = fs.asc[static_cast<std::size_t>(c)];
  const std::vector<int>& de = fs.desc[static_cast<std::size_t>(c)];

  for (int k : ac) w[k] = b[o[k]];
  for (int j : ac) {
    const double yj = w[j];
    for (SpMat::InnerIterator it(F, j); it; ++it) {
      const int i = static_cast<int>(it.row());
      if (i > j) w[i] -= it.value() * yj;
    }
  }
  for (int j : ac) w[j] *= fs.dinv[j];
  for (int j : de) {
    double s = w[j];
    for (SpMat::InnerIterator it(F, j); it; ++it) {
      const int i = static_cast<int>(it.row());
      if (i > j) s -= it.value() * w[i];
    }
    w[j] = s;
  }
  for (int k : ac) x[o[k]] = w[k];
}

struct BenchResult {
  double serialMinMs = 0.0;   ///< 区域外纯串行（无池参与）
  double aMinMs = 0.0;        ///< master 整趟
  double bMinMs = 0.0;        ///< omp for 三分量
  double cMinMs = 0.0;        ///< master 串行跑三分量
  double nullMinMs = 0.0;     ///< 空任务（一对 barrier 的固定成本）
  int reps = 0;
  int threads = 0;
};

/// 微基准。串行对照先跑（必须在任何多线程区域之前 —— vcomp 的池在区域结束后
/// 还会自旋一阵，紧跟其后测"串行"会被邻居抢 CPU，这是 docs/solver-feasibility.md §1
/// 记录过的测量事故）。
BenchResult runBench(const FactorSet& fs, const std::vector<double>& b, std::vector<double>& x,
                     int reps, int threads) {
  BenchResult r;
  r.reps = reps;
  r.threads = threads;
  // 复用缓冲：3 条线程各写自己的分量 ⇒ 内存不相交，共用一块不会竞争（也与生产一致）。
  std::vector<double> w(static_cast<std::size_t>(fs.n), 0.0);

  r.serialMinMs = 1e30;
  for (int i = 0; i < reps; ++i) {
    std::fill(x.begin(), x.end(), 0.0);
    const auto t0 = Clock::now();
    solveFull(fs, w.data(), b.data(), x.data());
    r.serialMinMs = std::min(r.serialMinMs, secondsSince(t0) * 1e3);
  }

  double aMin = 1e30;
  double bMin = 1e30;
  double cMin = 1e30;
  double nMin = 1e30;

#ifdef _OPENMP
#pragma omp parallel num_threads(threads)
#endif
  {
#ifdef _OPENMP
    const int tid = omp_get_thread_num();
#else
    const int tid = 0;
#endif
    for (int rep = 0; rep < reps; ++rep) {
      // ---- N：空的 omp for（3 个迭代、零工作）—— 量化 worksharing + barrier 的固定成本 ----
      {
        const auto tn0 = Clock::now();
#ifdef _OPENMP
#pragma omp for schedule(static, 1)
#endif
        for (int c = 0; c < kComponents; ++c) {
          // 故意空转：只量构造本身
        }
        const double t = secondsSince(tn0) * 1e3;
        if (tid == 0) nMin = std::min(nMin, t);
      }

      // ---- A：现状形状（单线程整趟）----
      {
        double t = 0.0;
#ifdef _OPENMP
#pragma omp master
#endif
        {
          const auto t0 = Clock::now();
          solveFull(fs, w.data(), b.data(), x.data());
          t = secondsSince(t0) * 1e3;
        }
#ifdef _OPENMP
#pragma omp barrier
#endif
        if (tid == 0) aMin = std::min(aMin, t);
      }

      // ---- B：生产要用的形状（3 个分量分给线程；计时含等最慢者）----
      {
        const auto tb0 = Clock::now();
#ifdef _OPENMP
#pragma omp for schedule(static, 1)
#endif
        for (int c = 0; c < kComponents; ++c) {
          solveOneComponent(fs, c, w.data(), b.data(), x.data());
        }
        const double t = secondsSince(tb0) * 1e3;
        if (tid == 0) bMin = std::min(bMin, t);
      }

      // ---- C：同样按分量拆，但串行跑三趟（只改循环结构，无并行）----
      {
        double t = 0.0;
#ifdef _OPENMP
#pragma omp master
#endif
        {
          const auto t0 = Clock::now();
          for (int c = 0; c < kComponents; ++c) solveOneComponent(fs, c, w.data(), b.data(), x.data());
          t = secondsSince(t0) * 1e3;
        }
#ifdef _OPENMP
#pragma omp barrier
#endif
        if (tid == 0) cMin = std::min(cMin, t);
      }
    }
  }

  r.aMinMs = aMin;
  r.bMinMs = bMin;
  r.cMinMs = cMin;
  r.nullMinMs = nMin;
  return r;
}

void printUsage() {
  std::printf("用法：pd_solvecomp [--grid N M] [--stiffness K] [--reps R] [--threads T]\n");
}

}  // namespace

int main(int argc, char** argv) {
  int nx = 40;
  int ny = 40;
  double stiffness = 1.0e4;
  int reps = 200;
  int threads = 4;

  for (int i = 1; i < argc; ++i) {
    const std::string a = argv[i];
    const bool hasNext = (i + 1 < argc);
    if (a == "--grid" && i + 2 < argc) {
      nx = std::atoi(argv[++i]);
      ny = std::atoi(argv[++i]);
    } else if (a == "--stiffness" && hasNext) {
      stiffness = std::atof(argv[++i]);
    } else if (a == "--reps" && hasNext) {
      reps = std::atoi(argv[++i]);
    } else if (a == "--threads" && hasNext) {
      threads = std::atoi(argv[++i]);
    } else if (a == "--help" || a == "-h") {
      printUsage();
      return 0;
    } else {
      std::printf("未知参数：%s\n", a.c_str());
      printUsage();
      return 2;
    }
  }
  if (reps < 1) reps = 1;
  if (threads < 1) threads = 1;

  std::printf("=== 全局步回代：按 x/y/z 三分量拆分的可行性与收益 ===\n");
  std::printf("网格 %dx%d  刚度 %.3g  每档次数 %d  区域内线程数 %d\n", nx, ny, stiffness, reps,
              threads);

  SceneConfig cfg;
  cfg.gridNx = nx;
  cfg.gridNy = ny;
  cfg.stiffness = stiffness;
  cfg.maxIterations = 1;
  SimContext ctx = makeScene(cfg);

  SpMat L;
  assembleLeftHandSide(ctx.mesh, cfg.dt, cfg.velocityDamping, L);
  const int n = 3 * ctx.mesh.vertexCount();
  std::printf("顶点 %d  约束 %d  自由度 %d  L 的 nnz %lld\n", ctx.mesh.vertexCount(), ctx.mesh.edgeCount(),
              n, static_cast<long long>(L.nonZeros()));

  FactorSet fs;
  if (!buildFactor(L, fs)) {
    std::printf("**分解失败或置换自检不通过** —— 中止。\n");
    return 1;
  }
  std::printf("因子 nnz %lld\n", static_cast<long long>(fs.F->nonZeros()));

  // ---- ① 结构自检 ----
  std::printf("\n--- ① 结构自检：因子是否按 x/y/z 分量块对角 ---\n");
  for (int c = 0; c < kComponents; ++c) {
    std::printf("  分量 %d：列数 %zu（期望 %d）  因子 nnz %lld\n", c,
                fs.asc[static_cast<std::size_t>(c)].size(), n / kComponents,
                fs.nnzByComp[static_cast<std::size_t>(c)]);
  }
  std::printf("  跨分量的 i≠j 项（填充跨分量）：%lld  ⇒ %s\n", fs.crossLinks,
              fs.crossLinks == 0 ? "**块对角成立，3 条链完全独立**" : "**不成立：不能按分量拆**");
  if (fs.crossLinks != 0) {
    std::printf("  结论：本项目当前的 L 不是 Ã ⊗ I₃ 结构，拆分方案不适用。\n");
    return 1;
  }

  // ---- ② 正确性：V1（三分量） vs V0（整趟）逐位比较 ----
  std::vector<double> b(static_cast<std::size_t>(n), 0.0);
  {
    // 用真实右端的形状：b = (M/h²)·x̂（与 assembleInertialRhs 一致），再覆盖 pin 行
    std::vector<Vec3> predicted = ctx.mesh.positions;
    for (std::size_t i = 0; i < predicted.size(); ++i) {
      predicted[i].x += 1e-3 * static_cast<double>(i % 17);
    }
    Eigen::VectorXd bVec;
    assembleInertialRhs(ctx.mesh, predicted, cfg.dt, cfg.velocityDamping, bVec);
    applyPinRhs(ctx.mesh, bVec);
    for (int i = 0; i < n; ++i) b[static_cast<std::size_t>(i)] = bVec[i];
  }

  std::vector<double> x0(static_cast<std::size_t>(n), 0.0);
  std::vector<double> x1(static_cast<std::size_t>(n), 0.0);
  std::vector<double> wCheck(static_cast<std::size_t>(n), 0.0);
  solveFull(fs, wCheck.data(), b.data(), x0.data());
  for (int c = 0; c < kComponents; ++c) {
    solveOneComponent(fs, c, wCheck.data(), b.data(), x1.data());
  }

  std::printf("\n--- ② 正确性：按分量拆 vs 整趟 ---\n");
  const bool bitwise = std::memcmp(x0.data(), x1.data(), sizeof(double) * static_cast<std::size_t>(n)) == 0;
  double maxAbs = 0.0;
  for (int i = 0; i < n; ++i) {
    maxAbs = std::max(maxAbs, std::fabs(x0[static_cast<std::size_t>(i)] - x1[static_cast<std::size_t>(i)]));
  }
  std::printf("  最大逐元素差 %.3e  ⇒ %s\n", maxAbs,
              bitwise ? "**逐位相同**（预期的性质：每个分量的算术序列没变）" : "**不是逐位相同，需查**");

  // ---- ③ 微基准 ----
  std::printf("\n--- ③ 微基准（同一并行区域内逐轮交替，取最小值）---\n");
  const BenchResult r = runBench(fs, b, x1, reps, threads);
  const double nnzTotal = static_cast<double>(fs.F->nonZeros());
  const double bytes = nnzTotal * 12.0 * 2.0;  // 值 8B + 行号 4B，前代+回代各读一遍
  const auto gbs = [bytes](double ms) { return ms > 0.0 ? bytes / (ms * 1e-3) / 1e9 : 0.0; };

  std::printf("  区域外纯串行（无池）      %8.4f ms   %6.1f GB/s\n", r.serialMinMs, gbs(r.serialMinMs));
  std::printf("  A master 整趟（现状形状）  %8.4f ms   相对串行 %.3fx\n", r.aMinMs, r.aMinMs / r.serialMinMs);
  std::printf("  C master 串行跑三分量      %8.4f ms   相对整趟 %.3fx\n", r.cMinMs, r.cMinMs / r.aMinMs);
  std::printf("  B omp for 三分量（%d 线程）%8.4f ms   相对整趟 %.3fx   %6.1f GB/s\n", threads,
              r.bMinMs, r.bMinMs / r.aMinMs, gbs(r.bMinMs));
  std::printf("  N 空 omp for（3 迭代）    %8.4f ms   （占 B 的 %.1f%%）\n", r.nullMinMs,
              r.nullMinMs / r.bMinMs * 100.0);
  std::printf("\n  判读：① 若 C≈A，说明按分量拆本身不改变单线程速度（只是换了循环结构）；\n");
  std::printf("        ② B 相对 A 的倍数就是这 3 条独立链带来的并行收益（上界 3.00x）；\n");
  std::printf("        ③ N 是这一对同步的固定成本，剩下的才是有效工作。\n");
  return bitwise ? 0 : 3;
}
