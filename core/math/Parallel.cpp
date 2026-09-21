// core/math/Parallel.cpp
#include "core/math/Parallel.h"

namespace pd {

namespace {
// 默认线程数。
//
// **2026-09-21 标定**：原来是 `hardware_concurrency() - 2`（本机 = 18），理由是
// "20 个逻辑核里含 4 个 E-core，留点余量"。实测下来这个默认值是**每个工况下最差的**：
//
//   工况                        1T      2T      4T      8T     18T
//   40×40 / iters40 / 300 子步  6.06    5.94    6.01    6.56   8.34 ms   ← 18T 慢 39 %
//   100×100 / iters10           14.79   14.19   13.69   14.69  15.60 ms  ← 18T 慢 14 %
//   200×200 / iters10           108.2   105.9   102.5   103.3  103.8 ms
//   60×60 / iters2（交互档）     0.84    0.89    0.90    1.02   1.23 ms   ← 18T 慢 47 %
//
// 原因是三件事一起作用：① 每次并行调用的同步成本随参与者增多而上升
// （本机 8 P-core + 4 E-core 异构，barrier 要等最慢的那条线程，线程还会跨核迁移）；
// ② 串行段（Eigen 回代，占 50–70 %）在有池在旁边时自己变慢；
// ③ 每线程分到的活太小，摊不平同步成本。所以**最优线程数很小，且与规模关系不大**。
//
// 4 是这几个工况上共同的最优点附近（2 与 4 在噪声内，4 在大规模上更稳）。
// 这是**机器相关参数**，换机器或换数量级（≫10 万条约束）都要重新标定 ——
// 与刚度的处理方式一致（见 docs/perf.md §6 的标定记录）。
// 运行时仍可用 `setNumThreads()` 覆盖（bench 的 `--threads`、查看器的 `;` `'`）。
constexpr int kDefaultThreadCap = 4;

int g_numThreads = [] {
  const unsigned hw = std::thread::hardware_concurrency();
  if (hw == 0) return 1;
  const int headroom = static_cast<int>(hw > 4 ? hw - 2 : hw);
  return std::max(1, std::min(headroom, kDefaultThreadCap));
}();
}  // namespace

int numThreads() { return g_numThreads; }

void setNumThreads(int n) {
  if (n < 1) n = 1;
  g_numThreads = n;
#ifdef _OPENMP
  omp_set_num_threads(n);
#endif
}

int threadId() {
#ifdef _OPENMP
  return omp_in_parallel() ? omp_get_thread_num() : 0;
#else
  return 0;
#endif
}

int regionThreadCount() {
#ifdef _OPENMP
  return omp_in_parallel() ? omp_get_num_threads() : 1;
#else
  return 1;
#endif
}

}  // namespace pd
