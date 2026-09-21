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
inline void projectEdge(const Mesh& mesh, std::vector<Vec3>& targets, std::size_t c) {
  const Edge& e = mesh.edges[c];
  const Vec3 raw = mesh.positions[static_cast<std::size_t>(e.a)] -
                   mesh.positions[static_cast<std::size_t>(e.b)];
  const Scalar len = length(raw);
  if (len <= DistanceTerm::kMinLength) {
    targets[c] = Vec3{0, 0, 0};  // 退化保护：两点重合时方向未定义
  } else {
    targets[c] = raw * (e.restLength / len);
  }
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
inline void scatterEdge(const Mesh& mesh, const std::vector<Vec3>& targets, std::size_t c, Scalar* out) {
  const Edge& e = mesh.edges[c];
  const bool aPinned = mesh.isPinned(e.a);
  const bool bPinned = mesh.isPinned(e.b);
  if (aPinned && bPinned) return;
  const Vec3 contribution = targets[c] * e.stiffness;
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

/// 每线程一份全维缓冲：散射的加法是"同一顶点被多条边累加"，并行会竞争，
/// 而决策 D7 明确不用浮点原子加（要保证结果与线程数无关、且能映射到 GPU），
/// 于是用"线程私有缓冲 + 固定顺序归约"。
/// 注意它的代价是 O(P·dim)：清零与归约都按"线程数 × 自由度数"走，
/// 与边数无关 —— 这正是 docs/parallel-refactor.md Phase 2（图着色）要消除的东西。
ThreadLocalBuffers<Scalar>& scatterBuffers() {
  static ThreadLocalBuffers<Scalar> buffers;
  return buffers;
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
                                       const Scalar* base, Scalar* b, std::size_t dim) {
  // 语义：b = base + Σ_c κ_c (A_cᵀ d_c)。base == b 表示"累加到现有值"
  // （独立入口 scatterInto 复用这条实现时用），此时要求 b 已是基值。
  const bool accumulateOnly = (base == b);
  const int nt = regionThreadCount();

  if (nt <= 1) {
    // 单线程路径：基值先写、再按**边序号顺序**逐条累加 —— 与独立入口的串行分支逐位相同。
    // 保留这条路径有两个理由：
    //   ① 小网格上开区域本来就亏（见 kMinParallelEdges）；
    //   ② 它是"1 线程下新旧实现位级一致"的对照基准（Phase 1 的验收之一）。
    checkPinContract(mesh);
    if (!accumulateOnly) {
      for (std::size_t i = 0; i < dim; ++i) b[i] = base[i];
    }
    const int ne = mesh.edgeCount();
    for (long long c = 0; c < static_cast<long long>(ne); ++c) {
      scatterEdge(mesh, targets, static_cast<std::size_t>(c), b);
    }
    return;
  }

  // 多线程：清零 → 各线程累加自己那份 → 按线程号升序归约。
  // 求和规则与改造前完全一致（每线程按边序累加自己的分块；归约按线程号升序），
  // 只是把"拷基值"融进了归约那一趟，省掉整次 O(3N) 的 b = bBase 拷贝。
  ThreadLocalBuffers<Scalar>& buffers = scatterBuffers();
#ifdef _OPENMP
#pragma omp single
#endif
  {
    checkPinContract(mesh);
    buffers.ensureFor(nt, dim);
  }

  const int ne = mesh.edgeCount();
  const std::size_t chunk = 64;
#ifdef _OPENMP
#pragma omp for schedule(dynamic, 64)
#endif
  for (long long begin = 0; begin < static_cast<long long>(ne); begin += static_cast<long long>(chunk)) {
    const long long end =
        std::min<long long>(begin + static_cast<long long>(chunk), static_cast<long long>(ne));
    Scalar* local = buffers[threadId()].data();
    for (long long c = begin; c < end; ++c) scatterEdge(mesh, targets, static_cast<std::size_t>(c), local);
  }

  // 归约趟：每个输出元素独立，且对同一元素仍按 t = 0..P-1 的顺序相加 ⇒ 与改造前逐位相同。
  // 归约本身也并行（改造前是单线程跑的 P×dim 循环）。
#ifdef _OPENMP
#pragma omp for
#endif
  for (long long i = 0; i < static_cast<long long>(dim); ++i) {
    Scalar partial = 0.0;
    for (int t = 0; t < nt; ++t) partial += buffers[t][static_cast<std::size_t>(i)];
    b[i] = accumulateOnly ? (b[i] + partial) : (base[i] + partial);
  }
}

void DistanceTerm::scatterInto(const Mesh& mesh, const std::vector<Vec3>& targets, Scalar* b,
                               std::size_t dim) {
  const int ne = mesh.edgeCount();
  if (numThreads() <= 1 || ne < kMinParallelEdges) {
    // 串行分支：直接按边序累加进 b（不经过线程私有缓冲）。
    checkPinContract(mesh);
    for (long long c = 0; c < static_cast<long long>(ne); ++c) {
      scatterEdge(mesh, targets, static_cast<std::size_t>(c), b);
    }
    return;
  }
#ifdef _OPENMP
#pragma omp parallel num_threads(numThreads())
  {
    // base == b ⇒ 累加到现有值，与改造前的 parallel 分支语义相同（含同样的归约顺序）。
    scatterIntoInRegion(mesh, targets, b, b, dim);
  }
#endif
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
