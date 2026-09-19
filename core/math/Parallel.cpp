// core/math/Parallel.cpp
#include "core/math/Parallel.h"

namespace pd {

namespace {
// 默认线程数：留出余量，且不假设所有逻辑核等价
// （本机 20 个逻辑核含 4 个 E-core，见 docs/environment.md）。
int g_numThreads = [] {
  const unsigned hw = std::thread::hardware_concurrency();
  if (hw == 0) return 1;
  return static_cast<int>(std::max(1u, hw > 4 ? hw - 2 : hw));
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

}  // namespace pd
