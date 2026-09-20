// core/energy/DistanceTerm.cpp
#include "core/energy/DistanceTerm.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>

#include "core/math/Parallel.h"

namespace pd {

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

void DistanceTerm::project(const Mesh& mesh, std::vector<Vec3>& targets) {
  const int ne = mesh.edgeCount();
  targets.resize(static_cast<std::size_t>(ne));
  parallelFor(static_cast<std::size_t>(ne), [&](std::size_t begin, std::size_t end) {
    for (std::size_t c = begin; c < end; ++c) {
      const Edge& e = mesh.edges[c];

      // ---- 约束方向（全项目唯一约定）----
      //
      // 方向**只由边的端点顺序 (a,b) 决定，没有任何例外**：
      //
      //   u_c   := unit(x_a - x_b)
      //   A_c x  = x_a - x_b = r·u_c        （r = ‖x_a - x_b‖）
      //   d_c    = Π_c(A_c x) = ℓ·u_c       （只把长度改成 ℓ，方向就是 u_c）
      //
      // 不允许按"哪一端是 pinned"去翻转方向：投影只是把长度定长化，它必须与矩阵块
      // κ·A_cᵀA_c 来自**同一个** A_c。谁支撑谁（pinned / 自由）由组装阶段处理：
      //   · 矩阵块：对角 +κI、耦合 -κI，只跳过"两端都 pin"的约束（assembleMatrix）
      //   · 散射  ：b_a -= κ d_c、b_b += κ d_c（scatterInto）
      //   · pinned 行：由 applyPinRhs 整行覆盖（Assembler.cpp）
      const Vec3 raw = mesh.positions[static_cast<std::size_t>(e.a)] -
                       mesh.positions[static_cast<std::size_t>(e.b)];
      const Scalar len = length(raw);
      if (len <= kMinLength) {
        targets[c] = Vec3{0, 0, 0};
      } else {
        targets[c] = raw * (e.restLength / len);
      }
    }
  });
}

void DistanceTerm::scatterInto(const Mesh& mesh, const std::vector<Vec3>& targets, Scalar* b,
                               std::size_t dim) {
  const int ne = mesh.edgeCount();

  // ---- 前置条件检查（快速失败）----
  // 散射在两端只有一端被 pin 时要读 mesh.pinPositions[*] 做消元补偿。
  // 该缓冲的填充由调用方负责（Scene 的 makeScene / refreshPinPositions）。
  // 手工搭 Mesh 的调用方（测试、_verify 程序）很容易只写 pinned 而漏掉 pinPositions；
  // 那时越界读会表现为随机崩溃（0xC0000005），极难定位，所以这里直接查出来并停下。
  if (mesh.pinPositions.size() < mesh.positions.size()) {
    bool anyPinned = false;
    for (int v = 0; v < mesh.vertexCount() && !anyPinned; ++v) anyPinned = mesh.isPinned(v);
    if (anyPinned) {
      std::fprintf(stderr,
                   "[DistanceTerm::scatterInto] 违反前置条件：mesh 存在 pinned 顶点，但 "
                   "pinPositions 只有 %zu 项（顶点数 %zu）。\n"
                   "  pin 消元补偿需要 b_free += κ·q 里那个 q = pinPositions[pin 顶点]。\n"
                   "  正常流程由 makeScene / refreshPinPositions 填充；手工搭 Mesh 时请显式设置，\n"
                   "  或在设置 pinned 后把 pinPositions 同步为当前位置。\n",
                   mesh.pinPositions.size(), mesh.positions.size());
      std::fflush(stderr);
      std::abort();
    }
  }

  // 每线程私有缓冲 + 固定顺序归约：结果与线程数无关（可复现）。
  // 首期刻意不用浮点原子加，见 docs/plan.md 决策 D6。
  static ThreadLocalBuffers<Scalar> buffers;
  buffers.ensure(dim);

  const int nThreads = buffers.threadCount();

  if (nThreads <= 1 || ne < 256) {
    for (int c = 0; c < ne; ++c) {
      const Edge& e = mesh.edges[static_cast<std::size_t>(c)];
      // 散射项 κ_c d_c：与对角块 +κ_c I、耦合块 -κ_c I 配对。
      // 配错符号会让弹簧力整体反号（稳态跑到 ℓ - m g/κ）。
      const Vec3 contribution = targets[static_cast<std::size_t>(c)] * e.stiffness;

      const bool aPinned = mesh.isPinned(e.a);
      const bool bPinned = mesh.isPinned(e.b);
      if (aPinned && bPinned) continue;
      const std::size_t ia = static_cast<std::size_t>(e.a) * 3;
      const std::size_t ib = static_cast<std::size_t>(e.b) * 3;

      // 标准 PD 散射：b_a += κ d_c，b_b -= κ d_c。
      if (!aPinned) {
        b[ia + 0] += contribution.x;
        b[ia + 1] += contribution.y;
        b[ia + 2] += contribution.z;
      }
      if (!bPinned) {
        b[ib + 0] -= contribution.x;
        b[ib + 1] -= contribution.y;
        b[ib + 2] -= contribution.z;
      }

      // ---- pin 消元的补偿项（关键）----
      // pinned 端那一侧的行被 applyPinRhs 覆盖为 (对角 1, 右端 q)，等价于把该自由度消去。
      // 消去后必须把原来耦合项 -κ·x_pin 的贡献移到自由端的右端：b_free += κ·q。
      // 漏掉它会导致自由端缺了 pinned 点的位置偏移 —— 表现为自由点被拉向原点，
      // 平移系统后行为改变（无重力、已在静止长度时都会自行变形）。
      // 以单弹簧为例，正确的全局步是
      //     (m/h² + κ)·x = (m/h²)·x̂ + κ·d_c + κ·q
      if (aPinned && !bPinned) {
        const Vec3& q = mesh.pinPositions[static_cast<std::size_t>(e.a)];
        b[ib + 0] += e.stiffness * q.x;
        b[ib + 1] += e.stiffness * q.y;
        b[ib + 2] += e.stiffness * q.z;
      } else if (bPinned && !aPinned) {
        const Vec3& q = mesh.pinPositions[static_cast<std::size_t>(e.b)];
        b[ia + 0] += e.stiffness * q.x;
        b[ia + 1] += e.stiffness * q.y;
        b[ia + 2] += e.stiffness * q.z;
      }
    }
    return;
  }

#ifdef _OPENMP
#pragma omp parallel num_threads(nThreads)
  {
    const int tid = omp_get_thread_num();
    std::vector<Scalar>& local = buffers[tid];
#pragma omp for schedule(dynamic, 64)
    for (long long c = 0; c < ne; ++c) {
      const Edge& e = mesh.edges[static_cast<std::size_t>(c)];
      const bool aPinned = mesh.isPinned(e.a);
      const bool bPinned = mesh.isPinned(e.b);
      if (aPinned && bPinned) continue;
      // 标准 PD 散射：b_a += κ d_c，b_b -= κ d_c
      const Vec3 contribution = targets[static_cast<std::size_t>(c)] * e.stiffness;
      const std::size_t ia = static_cast<std::size_t>(e.a) * 3;
      const std::size_t ib = static_cast<std::size_t>(e.b) * 3;
      if (!aPinned) {
        local[ia + 0] += contribution.x;
        local[ia + 1] += contribution.y;
        local[ia + 2] += contribution.z;
      }
      if (!bPinned) {
        local[ib + 0] -= contribution.x;
        local[ib + 1] -= contribution.y;
        local[ib + 2] -= contribution.z;
      }
      // pin 消元补偿：b_free += κ·q（详见串行分支处的推导）
      if (aPinned && !bPinned) {
        const Vec3& q = mesh.pinPositions[static_cast<std::size_t>(e.a)];
        local[ib + 0] += e.stiffness * q.x;
        local[ib + 1] += e.stiffness * q.y;
        local[ib + 2] += e.stiffness * q.z;
      } else if (bPinned && !aPinned) {
        const Vec3& q = mesh.pinPositions[static_cast<std::size_t>(e.b)];
        local[ia + 0] += e.stiffness * q.x;
        local[ia + 1] += e.stiffness * q.y;
        local[ia + 2] += e.stiffness * q.z;
      }
    }
  }
  // 固定顺序（按线程号升序）归约。
  // 注：ThreadLocalBuffers 内部按线程号顺序累加，因此结果与线程数无关。
  static std::vector<Scalar> merged;
  merged.assign(dim, Scalar{0});
  {
    const int n = buffers.threadCount();
    for (int t = 0; t < n; ++t) {
      const std::vector<Scalar>& local = buffers[t];
      for (std::size_t i = 0; i < dim; ++i) merged[i] += local[i];
    }
  }
  // 注意：merged 是全量累加（含零），等价于按线程号顺序求和。
  // 由于只有被当前线程写过的位置非零，逐位结果与串行一致到浮点求和顺序。
  for (std::size_t i = 0; i < dim; ++i) b[i] += merged[i];
  return;
#else
  for (int c = 0; c < ne; ++c) {
    const Edge& e = mesh.edges[static_cast<std::size_t>(c)];
    const bool aPinned = mesh.isPinned(e.a);
    const bool bPinned = mesh.isPinned(e.b);
    if (aPinned && bPinned) continue;
    // 散射项 κ_c d_c：与对角块 +κ_c I、耦合块 -κ_c I 配对。
    // 配错符号会让弹簧力整体反号（稳态跑到 ℓ - m g/κ）。
    const Vec3 contribution = targets[static_cast<std::size_t>(c)] * e.stiffness;
    const std::size_t ia = static_cast<std::size_t>(e.a) * 3;
    const std::size_t ib = static_cast<std::size_t>(e.b) * 3;
    if (!aPinned) {
      b[ia + 0] += contribution.x;
      b[ia + 1] += contribution.y;
      b[ia + 2] += contribution.z;
    }
    if (!bPinned) {
      b[ib + 0] -= contribution.x;
      b[ib + 1] -= contribution.y;
      b[ib + 2] -= contribution.z;
    }
    // pin 消元补偿：b_free += κ·q（详见串行分支处的推导）
    if (aPinned && !bPinned) {
      const Vec3& q = mesh.pinPositions[static_cast<std::size_t>(e.a)];
      b[ib + 0] += e.stiffness * q.x;
      b[ib + 1] += e.stiffness * q.y;
      b[ib + 2] += e.stiffness * q.z;
    } else if (bPinned && !aPinned) {
      const Vec3& q = mesh.pinPositions[static_cast<std::size_t>(e.b)];
      b[ia + 0] += e.stiffness * q.x;
      b[ia + 1] += e.stiffness * q.y;
      b[ia + 2] += e.stiffness * q.z;
    }
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
