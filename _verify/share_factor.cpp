// _verify/share_factor.cpp
//
// 目的：实测「三分量因子只存一份（三条线程共用）」vs「各存一份（现状）」值不值得做。
// 只测量、不改生产代码。
//
// 背景（`docs/open-issues.md` §1.3 的机会）：
//   L = Ã ⊗ I₃ ⇒ 因子的三个分量块**结构同构、数值完全相同**（`pd_solvecomp` 已经钉住"每分量
//   因子 nnz 相等、填充不跨分量"）。生产代码现在是"一份 CSC 因子矩阵里三块各存一份"，
//   三条线程各读自己那 1/3。若只存一份：
//     · 内存足迹降到 1/3（200×200：39 MB → 13 MB，按每非零 12 B 估）；
//     · 三条线程读**同一份**数据，第二次起的访问落在 L3（13 MB < 本机 25 MB L3），
//       而现状是 39 MB 三块不同地址、L3 装不下 ⇒ **DRAM 流量可降到约 1/3**。
//   这一条正好打在瓶颈上：回代占单子步 64–74 %，且实测是**延迟受限**（达成带宽只有流式读写的
//   9–13 %），所以"少搬数据 + 让数据住进 L3"是有可能真赚的。
//
// 本工具做四件事：
//   ① 前提自检：把参考分量压成**紧凑单份**（列指针 + 分量内相对行号 + 值），再逐列核对
//      另外两个分量与它"结构同构 + 数值逐位相同"（含对角 1/D）—— 这是"只存一份"成立的前提。
//   ② 正确性：单份求解（B）与各份求解（A）的解**逐位比较**（算术序列按构造相同）。
//   ③ 微基准：在同一并行区域内逐轮交替测 A / B，各取最小值；并且**冷热两档**都测
//      （冷 = 每轮之间清一块 64 MB 缓冲，把因子从 L3 里赶出去）。
//      冷档必须做：热数据下的紧循环会高估 L3 命中的好处 —— 2026-09-22 那次"另存一份 CSR"
//      就是被这个坑骗过（微基准 +23~29 %，真实管线只剩 −8~−10 %，见 docs/perf.md §8.4）。
//   ④ 一次性成本：紧凑单份的构建耗时（O(nnz) 一遍），用来判断它能不能被 12000 次回代摊掉。
//
// 用法：
//   pd_sharefactor [--grid 40 40] [--stiffness 1e4] [--reps 40] [--threads 4] [--flush-mb 64]

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

/// 逐位相等（浮点不能用 == 之外的东西判断"位级相同"）
inline bool bitwiseEqual(double a, double b) {
  return std::memcmp(&a, &b, sizeof(double)) == 0;
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
  std::vector<int> rankInComp;              ///< 置换后列号 k → 它在自己分量里的序号
  std::vector<std::vector<int>> asc;        ///< 每个分量的列号（升序）
  std::vector<std::vector<int>> desc;       ///< 同上（降序）
  long long crossLinks = 0;                 ///< 跨分量的 i≠j 项（必须为 0）
};

/// 与生产代码 EigenDirectSolver::analyze/factorize 同构（自己调 AMD + NaturalOrdering LDLᵀ）。
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

  // 置换后列号 → 原索引（数值探针，别凭 Eigen 的索引语义猜）
  Eigen::VectorXd labels(fs.n);
  for (int i = 0; i < fs.n; ++i) labels[i] = static_cast<double>(i);
  const Eigen::VectorXd mapped = fs.P * labels;
  fs.origOfPerm.assign(static_cast<std::size_t>(fs.n), 0);
  for (int k = 0; k < fs.n; ++k) {
    fs.origOfPerm[static_cast<std::size_t>(k)] = static_cast<int>(std::lround(mapped[k]));
  }
  {  // 自检：确实是一个置换
    std::vector<int> sorted = fs.origOfPerm;
    std::sort(sorted.begin(), sorted.end());
    for (int i = 0; i < fs.n; ++i) {
      if (sorted[static_cast<std::size_t>(i)] != i) return false;
    }
  }

  fs.asc.assign(kComponents, {});
  fs.desc.assign(kComponents, {});   // 必须先开好 3 个槽位：下面按分量下标赋值
  fs.rankInComp.assign(static_cast<std::size_t>(fs.n), -1);
  for (int k = 0; k < fs.n; ++k) {
    const int c = fs.origOfPerm[static_cast<std::size_t>(k)] % kComponents;
    fs.rankInComp[static_cast<std::size_t>(k)] =
        static_cast<int>(fs.asc[static_cast<std::size_t>(c)].size());
    fs.asc[static_cast<std::size_t>(c)].push_back(k);
  }
  for (int c = 0; c < kComponents; ++c) {
    if (static_cast<int>(fs.asc[static_cast<std::size_t>(c)].size()) != fs.n / kComponents) return false;
    fs.desc[static_cast<std::size_t>(c)] = fs.asc[static_cast<std::size_t>(c)];
    std::reverse(fs.desc[static_cast<std::size_t>(c)].begin(), fs.desc[static_cast<std::size_t>(c)].end());
  }
  for (int j = 0; j < fs.n; ++j) {
    const int cj = fs.origOfPerm[static_cast<std::size_t>(j)] % kComponents;
    for (SpMat::InnerIterator it(*fs.F, j); it; ++it) {
      const int i = static_cast<int>(it.row());
      if (i == j) continue;
      if (fs.origOfPerm[static_cast<std::size_t>(i)] % kComponents != cj) fs.crossLinks += 1;
    }
  }
  fs.ok = true;
  return true;
}

/// 紧凑单份因子：只有参考分量一份，行号存"分量内相对序号"。
struct Compact {
  int m = 0;                        ///< 每个分量的列数（= n/3）
  std::vector<int> ptr;             ///< m+1，列指针
  std::vector<int> rowRel;          ///< 非零项的分量内相对行号
  std::vector<double> val;          ///< 非零项的值
  std::vector<double> dinv;         ///< 每列的对角倒数（参考分量的那 m 个）
  long long nnz = 0;
  int mismatchedColumns = 0;        ///< 与参考分量对不上的列数（必须为 0）
  int mismatchedDiag = 0;           ///< 对角 1/D 对不上的列数（必须为 0）
  double buildMs = 0.0;             ///< 构建耗时（一次性）
};

/// 把参考分量（c=0）压成紧凑单份，并逐列核对另外两个分量与它"同构 + 逐位相同"。
Compact buildCompact(const FactorSet& fs) {
  Compact cf;
  const auto t0 = Clock::now();
  cf.m = fs.n / kComponents;
  cf.ptr.assign(static_cast<std::size_t>(cf.m) + 1, 0);
  const SpMat& F = *fs.F;
  const std::vector<int>& ref = fs.asc[0];

  // 1) 参考分量：列 k ←→ 置换后列号 ref[k]；行号存分量内相对序号 rankInComp[i]
  for (int k = 0; k < cf.m; ++k) {
    cf.ptr[static_cast<std::size_t>(k)] = static_cast<int>(cf.nnz);
    for (SpMat::InnerIterator it(F, ref[static_cast<std::size_t>(k)]); it; ++it) {
      cf.rowRel.push_back(fs.rankInComp[static_cast<std::size_t>(it.row())]);
      cf.val.push_back(it.value());
      cf.nnz += 1;
    }
    cf.dinv.push_back(fs.dinv[ref[static_cast<std::size_t>(k)]]);
  }
  cf.ptr[static_cast<std::size_t>(cf.m)] = static_cast<int>(cf.nnz);
  cf.buildMs = secondsSince(t0) * 1e3;

  // 2) 逐列核对另外两个分量：结构（相对行号序列）与数值都必须与参考完全一致
  for (int c = 1; c < kComponents; ++c) {
    for (int k = 0; k < cf.m; ++k) {
      const int j = fs.asc[static_cast<std::size_t>(c)][static_cast<std::size_t>(k)];
      int e = cf.ptr[static_cast<std::size_t>(k)];
      bool bad = false;
      for (SpMat::InnerIterator it(F, j); it; ++it) {
        if (e >= cf.ptr[static_cast<std::size_t>(k) + 1]) { bad = true; break; }
        const int rel = fs.rankInComp[static_cast<std::size_t>(it.row())];
        if (rel != cf.rowRel[static_cast<std::size_t>(e)] ||
            !bitwiseEqual(it.value(), cf.val[static_cast<std::size_t>(e)])) {
          bad = true;
          break;
        }
        ++e;
      }
      if (bad || e != cf.ptr[static_cast<std::size_t>(k) + 1]) cf.mismatchedColumns += 1;
      if (!bitwiseEqual(fs.dinv[j], cf.dinv[static_cast<std::size_t>(k)])) cf.mismatchedDiag += 1;
    }
  }
  return cf;
}

/// A：现状形状 —— 每条线程扫**自己分量**的那些列（生产 `solveComponent` 逐行等价）。
inline void solveOwnSlice(const FactorSet& fs, int c, double* w, const double* b, double* x) {
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

/// B：只存一份 —— 三条线程读**同一份**紧凑因子（各自带自己的列映射与右端）。
/// 算术序列与 A 逐项对应（列升序/降序、项序、筛选条件都不变）⇒ 结果必然逐位相同。
inline void solveSharedSlice(const FactorSet& fs, const Compact& cf, int c, double* w, const double* b,
                             double* x) {
  const int* o = fs.origOfPerm.data();
  const std::vector<int>& ac = fs.asc[static_cast<std::size_t>(c)];
  const int m = cf.m;

  for (int k = 0; k < m; ++k) w[ac[static_cast<std::size_t>(k)]] = b[o[ac[static_cast<std::size_t>(k)]]];
  for (int k = 0; k < m; ++k) {
    const int j = ac[static_cast<std::size_t>(k)];
    const double yj = w[j];
    for (int e = cf.ptr[static_cast<std::size_t>(k)]; e < cf.ptr[static_cast<std::size_t>(k) + 1]; ++e) {
      const int rel = cf.rowRel[static_cast<std::size_t>(e)];
      if (rel > k) w[ac[static_cast<std::size_t>(rel)]] -= cf.val[static_cast<std::size_t>(e)] * yj;
    }
  }
  for (int k = 0; k < m; ++k) {
    w[ac[static_cast<std::size_t>(k)]] *= cf.dinv[static_cast<std::size_t>(k)];
  }
  for (int k = m - 1; k >= 0; --k) {
    const int j = ac[static_cast<std::size_t>(k)];
    double s = w[j];
    for (int e = cf.ptr[static_cast<std::size_t>(k)]; e < cf.ptr[static_cast<std::size_t>(k) + 1]; ++e) {
      const int rel = cf.rowRel[static_cast<std::size_t>(e)];
      if (rel > k) s -= cf.val[static_cast<std::size_t>(e)] * w[ac[static_cast<std::size_t>(rel)]];
    }
    w[j] = s;
  }
  for (int k = 0; k < m; ++k) x[o[ac[static_cast<std::size_t>(k)]]] = w[ac[static_cast<std::size_t>(k)]];
}

/// 清缓存：在一块远大于 L3 的缓冲上写一遍，把因子从 L3 里赶出去（冷档专用）。
inline void flushCache(std::vector<unsigned char>& buf) {
  std::memset(buf.data(), 0, buf.size());
}

struct Bench {
  double aHot = 0, bHot = 0;      ///< 1 线程（区域外，无池）
  double aHotP = 0, bHotP = 0;    ///< 区域内 P 线程（omp for）
  double aCold = 0, bCold = 0;    ///< 1 线程 + 每轮清缓存
  double aColdP = 0, bColdP = 0;  ///< 区域内 P 线程 + 每轮清缓存
  double nullMs = 0;              ///< 空 omp for（3 迭代）的固定成本
};

Bench runBench(const FactorSet& fs, const Compact& cf, const std::vector<double>& b,
               std::vector<double>& x, int reps, int threads, std::vector<unsigned char>& flushBuf) {
  Bench r;
  std::vector<double> w(static_cast<std::size_t>(fs.n), 0.0);
  const double inf = 1e30;
  r.aHot = r.bHot = r.aHotP = r.bHotP = r.aCold = r.bCold = r.aColdP = r.bColdP = r.nullMs = inf;

  // ---- 1 线程（区域外；必须排在任何多线程区域之前 —— 池会自旋抢核，见 solver-feasibility §1）----
  for (int i = 0; i < reps; ++i) {
    {
      const auto t0 = Clock::now();
      for (int c = 0; c < kComponents; ++c) solveOwnSlice(fs, c, w.data(), b.data(), x.data());
      r.aHot = std::min(r.aHot, secondsSince(t0) * 1e3);
    }
    {
      const auto t0 = Clock::now();
      for (int c = 0; c < kComponents; ++c) solveSharedSlice(fs, cf, c, w.data(), b.data(), x.data());
      r.bHot = std::min(r.bHot, secondsSince(t0) * 1e3);
    }
    if (!flushBuf.empty()) {
      flushCache(flushBuf);
      const auto t1 = Clock::now();
      for (int c = 0; c < kComponents; ++c) solveOwnSlice(fs, c, w.data(), b.data(), x.data());
      r.aCold = std::min(r.aCold, secondsSince(t1) * 1e3);
      flushCache(flushBuf);
      const auto t2 = Clock::now();
      for (int c = 0; c < kComponents; ++c) solveSharedSlice(fs, cf, c, w.data(), b.data(), x.data());
      r.bCold = std::min(r.bCold, secondsSince(t2) * 1e3);
    }
  }

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
      {
        const auto t0 = Clock::now();
#ifdef _OPENMP
#pragma omp for schedule(static, 1)
#endif
        for (int c = 0; c < kComponents; ++c) solveOwnSlice(fs, c, w.data(), b.data(), x.data());
        const double t = secondsSince(t0) * 1e3;
        if (tid == 0) r.aHotP = std::min(r.aHotP, t);
      }
      {
        const auto t0 = Clock::now();
#ifdef _OPENMP
#pragma omp for schedule(static, 1)
#endif
        for (int c = 0; c < kComponents; ++c) solveSharedSlice(fs, cf, c, w.data(), b.data(), x.data());
        const double t = secondsSince(t0) * 1e3;
        if (tid == 0) r.bHotP = std::min(r.bHotP, t);
      }
      if (!flushBuf.empty()) {
        // 冷档：清缓存放在计时之外，只测解本身
#ifdef _OPENMP
#pragma omp master
#endif
        flushCache(flushBuf);
#ifdef _OPENMP
#pragma omp barrier
#endif
        {
          const auto t0 = Clock::now();
#ifdef _OPENMP
#pragma omp for schedule(static, 1)
#endif
          for (int c = 0; c < kComponents; ++c) solveOwnSlice(fs, c, w.data(), b.data(), x.data());
          const double t = secondsSince(t0) * 1e3;
          if (tid == 0) r.aColdP = std::min(r.aColdP, t);
        }
#ifdef _OPENMP
#pragma omp master
#endif
        flushCache(flushBuf);
#ifdef _OPENMP
#pragma omp barrier
#endif
        {
          const auto t0 = Clock::now();
#ifdef _OPENMP
#pragma omp for schedule(static, 1)
#endif
          for (int c = 0; c < kComponents; ++c) solveSharedSlice(fs, cf, c, w.data(), b.data(), x.data());
          const double t = secondsSince(t0) * 1e3;
          if (tid == 0) r.bColdP = std::min(r.bColdP, t);
        }
      }
      {
        const auto t0 = Clock::now();
#ifdef _OPENMP
#pragma omp for schedule(static, 1)
#endif
        for (int c = 0; c < kComponents; ++c) { /* 空转：量 worksharing + barrier 的固定成本 */
        }
        const double t = secondsSince(t0) * 1e3;
        if (tid == 0) r.nullMs = std::min(r.nullMs, t);
      }
    }
  }
  return r;
}

void printUsage() {
  std::printf(
      "用法：pd_sharefactor [--grid N M] [--stiffness K] [--reps R] [--threads T] [--flush-mb F]\n"
      "  --flush-mb 0 表示不做冷档（不推荐：热气会高估 L3 命中的好处）\n");
}

}  // namespace

int main(int argc, char** argv) {
  // 关掉 stdout 缓冲：这个工具要跑大内存，万一中途崩了，至少要看得见它走到哪一步
  // （printf 的行缓冲/块缓冲在 abort/段错误时会把已打印的内容丢掉）。
  std::setvbuf(stdout, nullptr, _IONBF, 0);
  int nx = 40;
  int ny = 40;
  double stiffness = 1.0e4;
  int reps = 40;
  int threads = 4;
  int flushMb = 64;

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
    } else if (a == "--flush-mb" && hasNext) {
      flushMb = std::atoi(argv[++i]);
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

  std::printf("=== 三分量因子：只存一份（三条线程共用） vs 各存一份（现状） ===\n");
  std::printf("网格 %dx%d  刚度 %.3g  每档次数 %d  区域内线程数 %d  冷档清缓存 %d MB\n", nx, ny,
              stiffness, reps, threads, flushMb);

  SceneConfig cfg;
  cfg.gridNx = nx;
  cfg.gridNy = ny;
  cfg.stiffness = stiffness;
  cfg.maxIterations = 1;
  SimContext ctx = makeScene(cfg);

  SpMat L;
  assembleLeftHandSide(ctx.mesh, cfg.dt, cfg.velocityDamping, L);
  const int n = 3 * ctx.mesh.vertexCount();
  std::printf("顶点 %d  约束 %d  自由度 %d  L 的 nnz %lld\n", ctx.mesh.vertexCount(),
              ctx.mesh.edgeCount(), n, static_cast<long long>(L.nonZeros()));

  FactorSet fs;
  if (!buildFactor(L, fs)) {
    std::printf("**分解失败或置换自检不通过** —— 中止。\n");
    return 1;
  }
  const long long factorNnz = static_cast<long long>(fs.F->nonZeros());
  std::printf("因子 nnz %lld（整份，三块各一份）\n", factorNnz);
  if (fs.crossLinks != 0) {
    std::printf("  跨分量的 i≠j 项 %lld ≠ 0 ⇒ 不是 Ã ⊗ I₃ 结构，本实验不适用。\n", fs.crossLinks);
    return 1;
  }

  // ---- ① 前提自检 ----
  std::printf("\n--- ① 前提自检：三个分量切片是否同构且数值逐位相同 ---\n");
  const Compact cf = buildCompact(fs);
  const double wholeMb = static_cast<double>(factorNnz) * 12.0 / (1024.0 * 1024.0);
  const double oneMb = static_cast<double>(cf.nnz) * 12.0 / (1024.0 * 1024.0);
  std::printf("  每分量列数 m = %d（期望 %d）\n", cf.m, n / kComponents);
  std::printf("  因子 nnz：整份 %lld（≈%.1f MB）  单份 %lld（≈%.1f MB，占 %.1f%%）\n", factorNnz,
              wholeMb, cf.nnz, oneMb, 100.0 * static_cast<double>(cf.nnz) / static_cast<double>(factorNnz));
  std::printf("  逐列核对分量 1/2 与参考分量：(相对行号, 值) 不一致的列 = %d；对角 1/D 不一致 = %d  ⇒ %s\n",
              cf.mismatchedColumns, cf.mismatchedDiag,
              (cf.mismatchedColumns == 0 && cf.mismatchedDiag == 0)
                  ? "**三份完全一致 ⇒ 只存一份成立**"
                  : "**三份并不一致 ⇒ 只存一份不成立，下面的数字只作参考**");
  std::printf("  紧凑单份的构建（一次性，O(nnz)）：%.3f ms\n", cf.buildMs);

  // ---- ② 正确性：B 与 A 逐位比较 ----
  std::vector<double> b(static_cast<std::size_t>(n), 0.0);
  {
    std::vector<Vec3> predicted = ctx.mesh.positions;
    for (std::size_t i = 0; i < predicted.size(); ++i) {
      predicted[i].x += 1e-3 * static_cast<double>(i % 17);
    }
    Eigen::VectorXd bVec;
    assembleInertialRhs(ctx.mesh, predicted, cfg.dt, cfg.velocityDamping, bVec);
    applyPinRhs(ctx.mesh, bVec);
    for (int i = 0; i < n; ++i) b[static_cast<std::size_t>(i)] = bVec[i];
  }
  std::vector<double> xA(static_cast<std::size_t>(n), 0.0);
  std::vector<double> xB(static_cast<std::size_t>(n), 0.0);
  std::vector<double> wCheck(static_cast<std::size_t>(n), 0.0);
  for (int c = 0; c < kComponents; ++c) solveOwnSlice(fs, c, wCheck.data(), b.data(), xA.data());
  for (int c = 0; c < kComponents; ++c) solveSharedSlice(fs, cf, c, wCheck.data(), b.data(), xB.data());
  const bool bitwise =
      std::memcmp(xA.data(), xB.data(), sizeof(double) * static_cast<std::size_t>(n)) == 0;
  double maxAbs = 0.0;
  for (int i = 0; i < n; ++i) {
    maxAbs = std::max(maxAbs, std::fabs(xA[static_cast<std::size_t>(i)] - xB[static_cast<std::size_t>(i)]));
  }
  std::printf("\n--- ② 正确性：单份求解 vs 各份求解 ---\n");
  std::printf("  最大逐元素差 %.3e  ⇒ %s\n", maxAbs,
              bitwise ? "**逐位相同**（算术序列按构造一致）" : "**不是逐位相同，需查**");

  // ---- ③ 微基准 ----
  std::printf("\n--- ③ 微基准（同一并行区域内逐轮交替，取最小值）---\n");
  std::vector<unsigned char> flushBuf;
  if (flushMb > 0) flushBuf.assign(static_cast<std::size_t>(flushMb) * 1024 * 1024, 0);
  const Bench r = runBench(fs, cf, b, xB, reps, threads, flushBuf);

  const double bytesA = static_cast<double>(factorNnz) * 12.0 * 2.0;  // 值 8B + 行号 4B，前代+回代
  const auto gbs = [bytesA](double ms) { return ms > 0.0 ? bytesA / (ms * 1e-3) / 1e9 : 0.0; };
  const auto row = [](const char* label, double a, double bb, const char* unit) {
    std::printf("    %-22s %9.4f %9.4f   %5.3fx%s\n", label, a, bb, bb / a, unit);
  };
  std::printf("                           各存一份(A) 只存一份(B)   B/A\n");
  row("1 线程 · 热数据", r.aHot, r.bHot, "");
  row("1 线程 · 冷数据", r.aCold, r.bCold, "   ← 隔离内存足迹/L3 效应（这里没有并行）");
  row("并行 · 热数据", r.aHotP, r.bHotP, "");
  row("并行 · 冷数据", r.aColdP, r.bColdP, "   ← 生产工况的形状");
  std::printf("  （空 omp for 固定成本 %.4f ms；A 的达成带宽 %.1f GB/s，B %.1f GB/s）\n", r.nullMs,
              gbs(r.aColdP), gbs(r.bColdP));

  // ---- ④ 判读 ----
  std::printf("\n--- ④ 判读 ---\n");
  const double hotGain = 100.0 * (r.aHotP - r.bHotP) / r.aHotP;
  const double coldGain = 100.0 * (r.aColdP - r.bColdP) / r.aColdP;
  const double serialColdGain = flushBuf.empty() ? 0.0 : 100.0 * (r.aCold - r.bCold) / r.aCold;
  std::printf("  并行档：热数据 B 快 %.1f%%，冷数据 B 快 %.1f%%\n", hotGain, coldGain);
  if (!flushBuf.empty()) std::printf("  单线程冷档：B 快 %.1f%%（这一档只有内存足迹/L3 效应，不含并行）\n", serialColdGain);
  std::printf(
      "  判读规则：**只看冷档**。热档的紧循环会让因子住在 L3 里，高估「只存一份」的收益\n"
      "  （2026-09-22 的 CSR 变体就是被这一点骗过：微基准 +23~29 %%，真实管线只剩 −8~−10 %%，\n"
      "  见 docs/perf.md §8.4）。冷档收益若在噪声（<5 %%）内，结论就是「不做」。\n");
  return bitwise ? 0 : 3;
}
