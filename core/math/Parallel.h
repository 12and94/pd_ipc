// core/math/Parallel.h
// 极简并行原语。
//
// 设计说明（对应 docs/plan.md 决策 D7）：
//   首期**不使用浮点原子加**做归约，而是"每线程私有缓冲 + 固定顺序二次归约"。
//   这样并行结果与线程数无关，可以做逐位回归测试（determinism 测试依赖这一点），
//   而且这套语义可以原样映射到 GPU（规避对 float atomic 扩展的依赖）。
#pragma once

#include <algorithm>
#include <cstddef>
#include <thread>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace pd {

/// 全局可配置的线程数。
/// 注意：本机 20 个逻辑核里含 4 个 E-core（见 docs/environment.md），
/// 因此默认不取 hardware_concurrency 的全部，也不假设线性加速。
int numThreads();
void setNumThreads(int n);

/// 当前线程在**并行区域内**的编号。
/// 不在区域内（或未启用 OpenMP）时返回 0 —— 这样"区域内代码"在串行构建下
/// 也能原样跑通，不需要 `#ifdef` 分支（"一个子步一个区域"方案的基础）。
int threadId();

/// 当前**并行区域内**的实际线程数。
/// 与 numThreads() 的区别：后者是"请求值"，本函数是"区域内真值"
/// （`if` 子句把区域串行化时它是 1）。线程私有缓冲必须按这个值分配，
/// 否则串行化的区域会白白多清零/多归约 (P-1)×dim 个槽位。
int regionThreadCount();

/// 按块划分的并行 for。fn(begin, end) 处理半开区间 [begin, end)。
/// 调度方式为动态（schedule(dynamic)），因为局部步各元素耗时不完全一致。
template <typename Fn>
void parallelFor(std::size_t count, Fn&& fn) {
  const int nThreads = numThreads();
  if (nThreads <= 1 || count < 256) {
    fn(std::size_t{0}, count);
    return;
  }
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 64) num_threads(nThreads)
  for (long long i = 0; i < static_cast<long long>(count); ++i) {
    fn(static_cast<std::size_t>(i), static_cast<std::size_t>(i) + 1);
  }
#else
  // 无 OpenMP 时的串行回退：保证同一份代码在任何构建下都能跑。
  fn(std::size_t{0}, count);
#endif
}

/// 线程私有累加缓冲的辅助类：每个线程拿到一块独立缓冲，最后按固定顺序合并。
/// 用途：并行散射（scatter）到全局右端向量，避免原子竞争。
template <typename T>
class ThreadLocalBuffers {
 public:
  ThreadLocalBuffers() = default;

  /// 分配 nThreads 块，每块 size 个元素，并清零。
  void reset(std::size_t size) {
    const int n = std::max(1, numThreads());
    buffers_.assign(static_cast<std::size_t>(n), std::vector<T>(size, T{}));
  }

  /// 需要扩容时按当前线程数重建（内容清零）。
  void ensure(std::size_t size) {
    if (buffers_.empty() || buffers_.size() != static_cast<std::size_t>(std::max(1, numThreads())) ||
        buffers_[0].size() != size) {
      reset(size);
    } else {
      for (auto& b : buffers_) std::fill(b.begin(), b.end(), T{});
    }
  }

  /// 按**区域内实际线程数**准备：见 regionThreadCount() 的说明。
  void ensureFor(int nThreads, std::size_t size) {
    const std::size_t n = static_cast<std::size_t>(std::max(1, nThreads));
    if (buffers_.size() != n || buffers_.empty() || buffers_[0].size() != size) {
      buffers_.assign(n, std::vector<T>(size, T{}));
    } else {
      for (auto& b : buffers_) std::fill(b.begin(), b.end(), T{});
    }
  }

  int threadCount() const { return static_cast<int>(buffers_.size()); }
  std::vector<T>& operator[](int t) { return buffers_[static_cast<std::size_t>(t)]; }
  const std::vector<T>& operator[](int t) const { return buffers_[static_cast<std::size_t>(t)]; }

  /// 按线程编号升序合并到 out（固定顺序 => 结果可复现）。
  void reduceInto(std::vector<T>& out) const {
    std::fill(out.begin(), out.end(), T{});
    for (const auto& b : buffers_) {
      for (std::size_t i = 0; i < out.size(); ++i) out[i] += b[i];
    }
  }

 private:
  std::vector<std::vector<T>> buffers_;
};

}  // namespace pd
