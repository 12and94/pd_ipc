// core/energy/DistanceTerm.cpp
#include "core/energy/DistanceTerm.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>

#include "core/math/Parallel.h"

namespace pd {

namespace {

/// 并行阈值：边数小于它时，本地循环比"开一个并行区域"更划算。
/// 依据：实测每次区域开关在 18 线程下要 10–25 µs，而一次投影/散射的
/// 有效工作在几千条边时只有几微秒（见 docs/perf.md）。与改造前 `parallelFor`
/// 里的 `count < 256` 是同一个思路，只是那里的门槛没有按"区域开关的实测成本"校准过。
constexpr int kMinParallelEdges = 256;

/// 前置条件检查（快速失败）：一旦有 pinned 顶点，`pinPositions` 必须按顶点数填好，
/// 因为散射要读 `pinPositions[pin]` 做消元补偿。手工搭 Mesh 的代码极易只写 pinned
/// 而漏掉它，越界读会表现为随机崩溃（0xC0000005），所以这里直接查出来并停下。
/// 快路径只是一次 size 比较，开销可忽略。
void checkPinContract(const Mesh& mesh) {
  if (mesh.pinPositions.size() >= mesh.positions.size()) return;
  bool anyPinned = false;
  for (int v = 0; v < mesh.vertexCount() && !anyPinned; ++v) anyPinned = mesh.isPinned(v);
  if (!anyPinned) return;
  std::fprintf(stderr,
               "[DistanceTerm::scatter] 违反前置条件：mesh 存在 pinned 顶点，但 "
               "pinPositions 只有 %zu 项（顶点数 %zu）。\n"
               "  pin 消元补偿需要 b_free += κ·q 里那个 q = pinPositions[pin 顶点]。\n"
               "  正常流程由 makeScene / refreshPinPositions 填充；手工搭 Mesh 时请显式设置，\n"
               "  或在设置 pinned 后把 pinPositions 同步为当前位置。\n",
               mesh.pinPositions.size(), mesh.positions.size());
  std::fflush(stderr);
  std::abort();
}

/// 前置条件检查（快速失败）：约束着色必须已经建好。
/// 漏了它会让"按颜色分组"扫到空区间，结果是右端悄悄变成只有基值 —— 一个不会崩溃、
/// 但物理完全错的静默错误，所以这里直接查出来并停下（与 pinPositions 契约同样的风格）。
/// 正常路径：`Mesh::buildSparsityPattern()` 会建好；手搭 Mesh 的调用方由
/// `Mesh::ensureConstraintColoring()` 兜底（两条散射入口都在进区域之前调用它）。
void checkColoringContract(const Mesh& mesh) {
  const ConstraintColoring& coloring = mesh.constraintColoring();
  if (!coloring.built) {
    std::fprintf(stderr,
                 "[DistanceTerm::scatter] 违反前置条件：约束着色尚未构建（网格 %d 条约束）。\n"
                 "  请调用 Mesh::buildSparsityPattern() 或 Mesh::ensureConstraintColoring()，\n"
                 "  并且必须在进入并行区域**之前**调用。\n",
                 mesh.edgeCount());
    std::fflush(stderr);
    std::abort();
  }
  // 便宜的陈旧检查（O(1)）：着色与当前边表规模必须一致。
  // 这条能挡住"改了拓扑却没重建"的最常见形式（makeScene 建好的着色被后续替换边表复用），
  // 把它变成一句明确的报错，而不是 scatterEdge 里的越界读（实测会以 0xC0000005 崩掉）。
  if (static_cast<int>(coloring.colorOf.size()) != mesh.edgeCount()) {
    std::fprintf(stderr,
                 "[DistanceTerm::scatter] 违反前置条件：约束着色已过期 —— 着色对应 %zu 条约束，"
                 "当前网格有 %d 条。\n  改了拓扑（增删边）之后必须重新调用 "
                 "Mesh::buildSparsityPattern()。\n",
                 coloring.colorOf.size(), mesh.edgeCount());
    std::fflush(stderr);
    std::abort();
  }
  // 关联表（CSR）的规模校验**不在这里**做。
  //
  // 原因是一条实测（很反直觉，务必保留）：把"vertexStart/vertexEdges 规模对不对"这种
  // O(1) 校验加进本函数，会让 MSVC 重新安排 scatter 循环的代码布局 ——
  // 40×40 / 300 子步 / 40 迭代（12,000 次散射调用）实测散射阶段
  // **100–104 ms → 128–133 ms（+28 %）**，而"多出来的指令"本身只有几条、
  // 与增量完全不相称（1T 总时间 +4 %，4T 看不到差别）。
  // 结论：本函数在**每次迭代**都被调用，任何改动都可能被编译器的布局变化放大。
  // 关联表校验因此挪到 `stepOnce` 里**一次/子步**的位置
  //（`checkAdjacencyContract`，见 Integrator.cpp），只由它的唯一使用者（残差的 gather）负责。
}

/// 投影的核心（一条边）：Π_c(A_c x) = ℓ·unit(x_a − x_b)。
///
/// ---- 约束方向（全项目唯一约定）----
/// 方向**只由边的端点顺序 (a,b) 决定，没有任何例外**：
///   u_c   := unit(x_a - x_b)
///   A_c x  = x_a - x_b = r·u_c        （r = ‖x_a - x_b‖）
///   d_c    = Π_c(A_c x) = ℓ·u_c       （只把长度改成 ℓ，方向就是 u_c）
/// 不允许按"哪一端是 pinned"去翻转方向：投影只是把长度定长化，它必须与矩阵块
/// κ·A_cᵀA_c 来自**同一个** A_c。谁支撑谁（pinned / 自由）由组装阶段处理：
///   · 矩阵块：对角 +κI、耦合 -κI，只跳过"两端都 pin"的约束（assembleMatrix）
///   · 散射  ：b_a -= κ d_c、b_b += κ d_c（scatterEdge）
///   · pinned 行：由 applyPinRhs 整行覆盖（Assembler.cpp）
inline Vec3 projectEdgeValue(const Mesh& mesh, const Edge& e) {
  const Vec3 raw = mesh.positions[static_cast<std::size_t>(e.a)] -
                   mesh.positions[static_cast<std::size_t>(e.b)];
  const Scalar len = length(raw);
  if (len <= DistanceTerm::kMinLength) {
    return Vec3{0, 0, 0};  // 退化保护：两点重合时方向未定义
  }
  return raw * (e.restLength / len);
}

/// 分段两趟路径用的薄封装：把投影写进 targets。
inline void projectEdge(const Mesh& mesh, std::vector<Vec3>& targets, std::size_t c) {
  targets[c] = projectEdgeValue(mesh, mesh.edges[c]);
}

/// 散射的核心（一条边，累加到某个缓冲上）。
///
/// 散射项 κ_c d_c：与对角块 +κ_c I、耦合块 -κ_c I 配对。
/// 配错符号会让弹簧力整体反号，稳态会跑到 ℓ + m g/κ（而不是正确的 ℓ - m g/κ）。
///
/// ---- pin 消元的补偿项（关键）----
/// pinned 端那一侧的行被 applyPinRhs 覆盖为 (对角 1, 右端 q)，等价于把该自由度消去。
/// 消去后必须把原来耦合项 -κ·x_pin 的贡献移到自由端的右端：b_free += κ·q。
/// 漏掉它会导致自由端缺了 pinned 点的位置偏移 —— 表现为自由点被拉向原点，
/// 平移系统后行为改变（无重力、已在静止长度时都会自行变形）。
/// 以单弹簧为例，正确的全局步是
///     (m/h² + κ)·x = (m/h²)·x̂ + κ·d_c + κ·q
inline void scatterEdgeValue(const Mesh& mesh, const Edge& e, const Vec3& d, Scalar* out) {
  const bool aPinned = mesh.isPinned(e.a);
  const bool bPinned = mesh.isPinned(e.b);
  if (aPinned && bPinned) return;
  const Vec3 contribution = d * e.stiffness;
  const std::size_t ia = static_cast<std::size_t>(e.a) * 3;
  const std::size_t ib = static_cast<std::size_t>(e.b) * 3;

  // 标准 PD 散射：b_a += κ d_c，b_b -= κ d_c。
  if (!aPinned) {
    out[ia + 0] += contribution.x;
    out[ia + 1] += contribution.y;
    out[ia + 2] += contribution.z;
  }
  if (!bPinned) {
    out[ib + 0] -= contribution.x;
    out[ib + 1] -= contribution.y;
    out[ib + 2] -= contribution.z;
  }
  // pin 消元补偿：b_free += κ·q
  if (aPinned && !bPinned) {
    const Vec3& q = mesh.pinPositions[static_cast<std::size_t>(e.a)];
    out[ib + 0] += e.stiffness * q.x;
    out[ib + 1] += e.stiffness * q.y;
    out[ib + 2] += e.stiffness * q.z;
  } else if (bPinned && !aPinned) {
    const Vec3& q = mesh.pinPositions[static_cast<std::size_t>(e.b)];
    out[ia + 0] += e.stiffness * q.x;
    out[ia + 1] += e.stiffness * q.y;
    out[ia + 2] += e.stiffness * q.z;
  }
}

/// 分段两趟路径（project 之后再 scatter）用的薄封装：从物化好的 targets 里取投影。
/// 融合路径不走这里 —— 两条路共用上面的 `scatterEdgeValue`，这正是"逐位相同"的来源。
inline void scatterEdge(const Mesh& mesh, const std::vector<Vec3>& targets, std::size_t c,
                        Scalar* out) {
  scatterEdgeValue(mesh, mesh.edges[c], targets[c], out);
}

}  // namespace

Vec3 DistanceTerm::projectOne(const Vec3& pa, const Vec3& pb, Scalar restLength) {
  const Vec3 d = pa - pb;
  const Scalar len = length(d);
  if (len <= kMinLength) {
    // 退化保护：两点重合时方向未定义。返回零向量，由调用方（或测试）处理。
    return Vec3{0, 0, 0};
  }
  // Π_c(d) = ℓ · d/‖d‖ —— 投影结果的长度**恰好**等于静止长度，这是定义本身。
  return d * (restLength / len);
}

void DistanceTerm::projectInRegion(const Mesh& mesh, std::vector<Vec3>& targets) {
  const int ne = mesh.edgeCount();
  targets.resize(static_cast<std::size_t>(ne));
  // nowait：调用方（stepOnce 的迭代循环）紧接着的同步点会覆盖它 ——
  // 散射的入口 single 带隐式 barrier。独立入口 project() 不使用 nowait。
#ifdef _OPENMP
#pragma omp for schedule(dynamic, 64) nowait
#endif
  for (long long c = 0; c < static_cast<long long>(ne); ++c) {
    projectEdge(mesh, targets, static_cast<std::size_t>(c));
  }
}

void DistanceTerm::project(const Mesh& mesh, std::vector<Vec3>& targets) {
  const int ne = mesh.edgeCount();
  targets.resize(static_cast<std::size_t>(ne));
  if (numThreads() <= 1 || ne < kMinParallelEdges) {
    for (long long c = 0; c < static_cast<long long>(ne); ++c) {
      projectEdge(mesh, targets, static_cast<std::size_t>(c));
    }
    return;
  }
#ifdef _OPENMP
#pragma omp parallel num_threads(numThreads())
  { projectInRegion(mesh, targets); }
#endif
}

void DistanceTerm::scatterIntoInRegion(const Mesh& mesh, const std::vector<Vec3>& targets,
                                       const Scalar* base, Scalar* b, std::size_t dim,
                                       const Scalar* extraRhs) {
  // 语义：b = base + extraRhs + Σ_c κ_c (A_cᵀ d_c)。base == b 表示"累加到现有值"
  // （独立入口 scatterInto 复用这条实现时用），此时要求 b 已是基值。
  const bool accumulateOnly = (base == b);

  // 前置检查各查一次即可（违反时 abort，不属于热路径；读-only，无竞争）。
  if (threadId() == 0) {
    checkPinContract(mesh);
    checkColoringContract(mesh);
  }

  const ConstraintColoring& coloring = mesh.constraintColoring();

  if (!accumulateOnly) {
    // 种子趟：b = base + extraRhs。取代了旧方案的"每线程一份全维缓冲清零 + 归约"（O(P·dim)），
    // 这一趟只有 O(dim) —— 这就是 Phase 2 省下来的主要流量。
    //
    // `extraRhs` 就是**线性（中点）弯曲的常向量**（见 Assembler.h 的推导）：它与位形无关，
    // 所以在唯一一处"b 被整个覆盖"的地方顺手加上，每迭代自动重新加一遍，**零额外 barrier**。
    // 不要在热循环里判 `extraRhs` 是否为空以外的任何条件 —— 按 docs/perf.md §7.3 的教训，
    // 每次迭代都跑的代码里加任何一点东西都可能被编译器的布局变化放大。
#ifdef _OPENMP
#pragma omp for
#endif
    for (long long i = 0; i < static_cast<long long>(dim); ++i) {
      b[i] = base[i] + (extraRhs ? extraRhs[i] : Scalar{0});
    }
  }

  // 逐色累加：同色边两两不共享顶点 ⇒ 并发写 b **没有任何冲突**，不需要私有缓冲、
  // 不需要原子加、也不需要归约。
  //
  // 求和顺序 = 颜色序 0,1,…,C-1，而同一顶点在每一色里**最多被写一次**
  // ⇒ 结果与线程数、分块、调度完全无关（位级可复现，见 tests/primitives 的跨线程断言）。
  // 这与旧方案（每线程按分块累加、再按线程号归约）不同：那时"哪些边先相加"取决于
  // schedule(dynamic) 的分块分配，位级结果随运行漂移。
  //
  // 代价：每个颜色轮是一次 `omp for`，因此每次散射有 C 个 barrier。
  // 换来的是 O(P·dim) → O(dim + E) 的流量，实测在大网格上划算（docs/perf.md）。
  //
  // 调度用默认的 static（连续分块）而不是 dynamic：同色内已按顶点序排列
  // （见 buildConstraintColoring），连续分块能让每个线程只写**自己那一段顶点**，
  // 把 cache line 争抢降到区段边界；而且 static 的分块分配是确定的，与运行无关。
  for (int c = 0; c < coloring.colorCount; ++c) {
    const long long begin = static_cast<long long>(coloring.begin(c));
    const long long end = static_cast<long long>(coloring.end(c));
    if (begin == end) continue;  // 空色类：不浪费一次 barrier（所有线程看到同样的值）
#ifdef _OPENMP
#pragma omp for
#endif
    for (long long k = begin; k < end; ++k) {
      scatterEdge(mesh, targets, coloring.order[static_cast<std::size_t>(k)], b);
    }
  }
}

void DistanceTerm::scatterInto(const Mesh& mesh, const std::vector<Vec3>& targets, Scalar* b,
                               std::size_t dim, const Scalar* extraRhs) {
  const int ne = mesh.edgeCount();
  // 懒构建兜底：**必须在进区域之前**（见 Mesh::ensureConstraintColoring 的说明）。
  mesh.ensureConstraintColoring();
  // 与 stepOnce 同一套做法：小网格 / 单线程时用 `if` 子句把区域**串行化**，
  // 而不是另写一条"直接按边序累加"的串行分支 ——
  // 求和顺序全项目只有一条（颜色序），这样 1 线程与 18 线程的结果才能逐位相同。
  const bool useParallel = (numThreads() > 1) && (ne >= kMinParallelEdges);
#ifdef _OPENMP
#pragma omp parallel num_threads(numThreads()) if (useParallel)
#endif
  {
    // base == b ⇒ 累加到现有值，与改造前 parallel 分支的语义相同
    // （此时 extraRhs 被忽略 —— 调用方必须先把它加进 b 的基值里）。
    scatterIntoInRegion(mesh, targets, b, b, dim, extraRhs);
  }
}

void DistanceTerm::projectAndScatterIntoInRegion(const Mesh& mesh, const Scalar* base, Scalar* b,
                                                std::size_t dim, const Scalar* extraRhs) {
  // 与 scatterIntoInRegion 同一套契约与前置检查（各查一次，O(1)，非热路径）。
  const bool accumulateOnly = (base == b);
  if (threadId() == 0) {
    checkPinContract(mesh);
    checkColoringContract(mesh);
  }

  const ConstraintColoring& coloring = mesh.constraintColoring();

  if (!accumulateOnly) {
    // 种子趟：b = base + extraRhs（extraRhs = 弯曲的常向量，见 Assembler.h）。
    // 保留不动 —— 它只有 O(dim) 流量与 1 个 barrier，
    // 而"把首色并进基数写"需要在每条边的热循环里判"该顶点是否首次被触碰"，
    // 按 docs/perf.md §7.3 的教训（热循环里加一点点东西都可能被布局放大）不划算。
#ifdef _OPENMP
#pragma omp for
#endif
    for (long long i = 0; i < static_cast<long long>(dim); ++i) {
      b[i] = base[i] + (extraRhs ? extraRhs[i] : Scalar{0});
    }
  }

  // 逐色融合：同色内的边两两不共享顶点 ⇒ 就地算投影、直接累加**没有任何冲突**。
  // 求和顺序与两趟版本完全一致（颜色序；同色内每个顶点最多被写一次）⇒ 结果逐位相同。
  for (int c = 0; c < coloring.colorCount; ++c) {
    const long long begin = static_cast<long long>(coloring.begin(c));
    const long long end = static_cast<long long>(coloring.end(c));
    if (begin == end) continue;  // 空色类：不浪费一次 barrier
#ifdef _OPENMP
#pragma omp for
#endif
    for (long long k = begin; k < end; ++k) {
      const Edge& e = mesh.edges[coloring.order[static_cast<std::size_t>(k)]];
      scatterEdgeValue(mesh, e, projectEdgeValue(mesh, e), b);
    }
  }
}

void DistanceTerm::assembleMatrix(const Mesh& mesh, const std::vector<Vec3>&,
                                  std::vector<Eigen::Triplet<Scalar>>& triplets) {
  // 约束块（标准 PD 形式，与位形无关 —— 这是"预分解复用"的前提）：
  //   对角块  +κ_c I      （顶点 a、b 各一份）
  //   耦合块  -κ_c I      （仅当两端都自由）
  //
  // 与投影方向的对应关系：矩阵块是 κ·A_cᵀA_c，其中 A_c 的取法见 project() 顶部
  // 那一段"约束方向（全项目唯一约定）"。A_cᵀA_c 本身与方向无关（翻转 A_c 不变），
  // 所以本函数**不需要**也不会引入任何与 pin 有关的方向分支 —— 这一点与
  // scatterInto 不同（那里用的是 d_c，与方向有关）。
  // 两处必须同源：投影定方向 → 散射用该方向 → 矩阵提供同源的刚度。
  //
  // 只跳过"两端都被 pin"的约束；只要还有自由端点，那个端点的对角贡献就必须进矩阵，
  // 否则该顶点丢掉这份刚度、系统明显偏软。
  for (const auto& e : mesh.edges) {
    const bool aPinned = mesh.isPinned(e.a);
    const bool bPinned = mesh.isPinned(e.b);
    if (aPinned && bPinned) continue;

    const int ia = e.a * 3;
    const int ib = e.b * 3;
    for (int d = 0; d < 3; ++d) {
      if (!aPinned) triplets.emplace_back(ia + d, ia + d, e.stiffness);
      if (!bPinned) triplets.emplace_back(ib + d, ib + d, e.stiffness);
      if (!aPinned && !bPinned) {
        triplets.emplace_back(ia + d, ib + d, -e.stiffness);
        triplets.emplace_back(ib + d, ia + d, -e.stiffness);
      }
    }
  }
}

}  // namespace pd
