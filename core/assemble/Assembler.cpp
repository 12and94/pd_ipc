// core/assemble/Assembler.cpp
#include "core/assemble/Assembler.h"

#include <bit>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "core/math/Parallel.h"

namespace pd {

namespace {

/// 把 double 的位模式取出来，用于 SolverStamp 的精确比较。
uint64_t bitsOf(Scalar v) {
  return static_cast<uint64_t>(std::bit_cast<uint64_t>(v));
}

/// FNV-1a 风格的 64 位混合，用于把一串精确位模式折成一个哈希。
/// 目的只是"变了就能看出来"，不比密码学强度。
uint64_t mix(uint64_t h, uint64_t v) {
  h ^= v + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
  return h;
}

}  // namespace

SolverStamp computeStamp(const Mesh& mesh, Scalar dt, Scalar damping, uint64_t topologyId) {
  // damping 参数保留在签名里（调用方语义完整），但它不影响数值戳 —— 见 computeStamp 尾部说明。
  SolverStamp s;
  s.topologyId = topologyId;
  s.dtBits = bitsOf(dt);

  uint64_t stiffnessHash = 0xcbf29ce484222325ULL;
  uint64_t massHash = 0xcbf29ce484222325ULL;
  for (const auto& e : mesh.edges) stiffnessHash = mix(stiffnessHash, bitsOf(e.stiffness));
  // **弯曲刚度也要进数值戳**：弯曲块 (p,q) = k·w_p·w_q·I₃ 直接乘在 k 上，
  // 因此 k 变了 L 就变了、必须重分解。漏掉它的表现是"调弯曲刚度好像没反应"——
  // 因为分解出来的因子还是旧的（这个坑与"把不改变 L 的量放进 stamp"是镜像的：
  // 那一边是白做一次分解，这一边是**该做而不做**，后者是错的物理结果）。
  for (const auto& bend : mesh.bends) stiffnessHash = mix(stiffnessHash, bitsOf(bend.stiffness));
  for (Scalar m : mesh.masses) massHash = mix(massHash, bitsOf(m));

  // pin **掩码**：只有"哪些顶点被 pin"会改变 L。
  // 被 pin 的行被覆盖为 (对角 1, 右端 q)，对角值与 q 无关；约束块只看端点是否被 pin。
  // 所以 pinPositions 的数值**不进 stamp**——否则拖拽把手（只改 b）会触发一次
  // 毫无必要的数值重分解，而 docs/plan.md §2.2.1 明确要求"移动把手不触发重分解"。
  // 注意这里不能对 pinned 顶点取"位置相关"的哈希：掩码是拓扑性质的量。
  uint64_t pinMask = 0xcbf29ce484222325ULL;
  for (std::size_t v = 0; v < mesh.pinned.size(); ++v) {
    pinMask = mix(pinMask, mesh.pinned[v] ? 1ULL : 0ULL);
  }

  s.stiffnessBits = stiffnessHash;
  s.massBits = massHash;
  s.pinBits = pinMask;

  // 阻尼不进 L：L 的惯性项是未缩放的 M/h²，阻尼只体现在速度更新
  // v_{n+1} = (1-k_d)(x_{n+1}-x_n)/h 上（见 Integrator.cpp 第 5 步）。
  // 因此改阻尼只需要重算右端，不该触发重分解。dampingBits 仍按位记录，
  // 供上层诊断/回归对比使用，但不参与"是否需要重分解"的判定。
  s.dampingBits = bitsOf(damping);
  return s;
}

uint64_t computeStructureStamp(const Mesh& mesh, uint64_t topologyId) {
  // 顺序敏感：同一个 mesh 每次调用都得到同一个值（edges / bends 与 pinned 都不会在
  // 帧内被重排），所以可以直接逐项混合。
  uint64_t h = 0xcbf29ce484222325ULL;
  h = mix(h, topologyId);
  h = mix(h, static_cast<uint64_t>(mesh.vertexCount()));
  h = mix(h, static_cast<uint64_t>(mesh.edgeCount()));
  for (const auto& e : mesh.edges) {
    h = mix(h, static_cast<uint64_t>(static_cast<uint32_t>(e.a)));
    h = mix(h, static_cast<uint64_t>(static_cast<uint32_t>(e.b)));
  }
  // 弯曲 stencil 的顶点三元组：它们引入新的 3×3 块（可能落在原先完全为空的位置），
  // 所以**结构变了必须重新符号分解**。只哈希顶点号、不哈希刚度 —— 刚度改变的是
  // 数值（走 computeStamp），稀疏结构一个字都不变，不该触发 analyzePattern。
  h = mix(h, static_cast<uint64_t>(mesh.bendCount()));
  for (const auto& s : mesh.bends) {
    h = mix(h, static_cast<uint64_t>(static_cast<uint32_t>(s.a)));
    h = mix(h, static_cast<uint64_t>(static_cast<uint32_t>(s.b)));
    h = mix(h, static_cast<uint64_t>(static_cast<uint32_t>(s.c)));
  }
  // pin 掩码：它决定哪些行被覆盖、哪些耦合块被跳过 —— 会改变稀疏结构。
  for (std::size_t v = 0; v < mesh.pinned.size(); ++v) {
    h = mix(h, mesh.pinned[v] ? 1ULL : 0ULL);
  }
  return h;
}

void assembleBendingRhs(const Mesh& mesh, std::vector<Scalar>& out) {
  const std::size_t dim = static_cast<std::size_t>(3 * mesh.vertexCount());
  out.assign(dim, Scalar{0});
  if (mesh.bends.empty()) return;

  // 右端 = -k·w_p·C，其中 C = Σ_{t ∈ stencil 内 pinned 顶点} w_t · q_t（**常向量**）。
  //
  // 推导（与距离约束的 pin 消元补偿同源，只是"一端"变成"两端/中间"）：
  //   stencil 的能量梯度对自由顶点 p 是  [L·x]_p = k·w_p·Σ_q w_q x_q。
  //   把 x_q = q_q（pinned）的那部分移到右端，得到
  //       b_p  +=  -k·w_p·C,        C = Σ_{q pinned} w_q·q_q。
  //   **符号必须是负**：p 的右端多出来的是"被消去的那一侧已经替它承担了多少"。
  //   若写成 +k·w_p·C，自由点会被推着**朝 pin 位移的方向跑**（而不是被它带着走），
  //   平直构型下的不动点会从"三点共线"变成"偏离"—— 一条容易被忽略的方向错。
  //   自检：单侧 pin 的 stencil（a 被 pin、b/c 自由）在共线位置上，b 的平衡是
  //   "b 正好落在 q_a 与 x_c 的中点"，把上面的式子代进去恰好给出这个解。
  //
  // 全自由 stencil ⇒ C = 0 ⇒ 本条 stencil 对右端零贡献（绝大多数情况）。
  for (const auto& s : mesh.bends) {
    if (s.stiffness == 0.0) continue;  // 刚度 0 = 这条 stencil 等效不存在（也不进 L）
    const int vertex[3] = {s.a, s.b, s.c};
    const Scalar weight[3] = {1.0, -2.0, 1.0};

    Vec3 c{0, 0, 0};
    bool anyPinned = false;
    for (int k = 0; k < 3; ++k) {
      if (!mesh.isPinned(vertex[k])) continue;
      anyPinned = true;
      const Vec3& q = mesh.pinPositions[static_cast<std::size_t>(vertex[k])];
      c += q * weight[k];
    }
    if (!anyPinned) continue;

    for (int k = 0; k < 3; ++k) {
      const int p = vertex[k];
      if (mesh.isPinned(p)) continue;  // pinned 行会被 applyPinRhs 整行覆盖，写了也白写
      const Vec3 contribution = c * (-s.stiffness * weight[k]);
      Scalar* slot = out.data() + static_cast<std::size_t>(p) * 3;
      slot[0] += contribution.x;
      slot[1] += contribution.y;
      slot[2] += contribution.z;
    }
  }
}

void assembleLeftHandSide(const Mesh& mesh, Scalar dt, Scalar /*damping*/,
                          Eigen::SparseMatrix<Scalar>& L, std::vector<Scalar>* bendingRhs) {
  const int n = mesh.vertexCount();
  const int dim = 3 * n;
  // 惯性项用**未缩放**的 M/h²。damping 参数在这里刻意不使用：
  // 阻尼只体现在 Integrator 的速度更新 v_{n+1} = (1-k_d)(x_{n+1}-x_n)/h 上，
  // 不改变 L 的数值 —— 这也正是"改阻尼不需要重分解"的原因（见 computeStamp）。
  const Scalar invDt2 = 1.0 / (dt * dt);

  std::vector<Eigen::Triplet<Scalar>> triplets;
  triplets.reserve(static_cast<std::size_t>(n) * 3 + mesh.edgeCount() * 12 +
                   mesh.bendCount() * 27 + 3);

  // (1) 惯性项 M/h²（对角），与位形无关。
  for (int v = 0; v < n; ++v) {
    const Scalar d = mesh.masses[static_cast<std::size_t>(v)] * invDt2;
    triplets.emplace_back(v * 3 + 0, v * 3 + 0, d);
    triplets.emplace_back(v * 3 + 1, v * 3 + 1, d);
    triplets.emplace_back(v * 3 + 2, v * 3 + 2, d);
  }

  // (2) 约束贡献（常数块）。
  DistanceTerm::assembleMatrix(mesh, {}, triplets);

  // (2b) 线性（中点）弯曲约束的常数块。
  //
  // 能量 E_s = (k/2)‖A_s x‖²，A_s x = x_a - 2 x_b + x_c ⇒ 梯度对顶点 p 是
  //   [∇E_s]_p = k·w_p·(Σ_q w_q x_q),   w = {1, -2, 1}
  // 于是 Hessian 的块就是  k·w_p·w_q·I₃ —— **仍然是标量 × I₃**，
  // 所以 L = Ã ⊗ I₃ 与"按 x/y/z 三分量并行回代"都保住（这是选这个弯曲形式的全部理由）。
  //
  // ---- 规则：p、q **两个都自由**的块才进 L；任何一端被 pin 的块一律不进 ----
  //
  // 这是"pinned 行整行覆盖"这套消元法的**硬前提**，不是可选的口味问题。
  // 被 pin 的行被覆盖为 (对角 1, 右端 q) 等价于把该自由度消去；这套写法只在
  // **pinned 列在所有自由行里都是 0** 时才自洽 —— 否则后面 Lᵀ 回代会用那个耦合
  // 把 x_pin 从 (b_pin = q) 上"顶"回一个别的值，pinned 顶点就不再严格等于把手位置。
  // 实测症状（本次先写错成"只跳过 pinned 行、保留 pinned 列"，很值得记）：
  //   · 三顶点小链、只 pin 两端：自由端点的小块是 [[1, 0, -2k], [0, 18400, -2k], [-2k, 0, 1]]，
  //     它不是半正定的（行列式 18400 - 4×10⁶ < 0），**求解器完全不报错**，
  //     而是给出一个谁都不满足的解（x_0 = -0.924 而 pin 行要求 x_0 = 0）；
  //   · 线性残差 |Lx-b|∞/|b|∞ = 8.9e-2（应该到 1e-12），而 Eigen 自带的
  //     SimplicialLDLT 给出**完全相同**的错误解 ⇒ 说明错在系统、不在求解器。
  //
  // 与距离约束的写法一致：`DistanceTerm::assembleMatrix` 也是"只有两端都自由才写耦合块"。
  //
  // **pinned 那边的贡献去了右端**：把 `x_q = q_q`（pinned）那部分移到右端，就是
  // `assembleBendingRhs` 里的常向量 `-k·w_p·C`（C = Σ_{pinned q} w_q·q_q）。
  // 两边是同一个式子的两半，必须成对出现：漏了右端 ⇒ 自由点丢掉 pin 的位置信息
  // （与距离约束漏 κ·q 补偿同源，README §4.2）；漏了这里的"不写" ⇒ pinned 顶点不严格。
  //
  // 写进 pinned **行**的东西同样不用管（下面 (3) 会把它补成 1）。
  for (const auto& s : mesh.bends) {
    if (s.stiffness == 0.0) continue;
    const int vertex[3] = {s.a, s.b, s.c};
    const Scalar weight[3] = {1.0, -2.0, 1.0};
    if (mesh.isPinned(s.a) && mesh.isPinned(s.b) && mesh.isPinned(s.c)) continue;
    for (int p = 0; p < 3; ++p) {
      if (mesh.isPinned(vertex[p])) continue;  // pinned 行会被整行覆盖：只算自由端的行
      const int ip = vertex[p] * 3;
      for (int q = 0; q < 3; ++q) {
        if (mesh.isPinned(vertex[q])) continue;  // pinned 列必须为 0（见上面的推导）
        const Scalar value = s.stiffness * weight[p] * weight[q];
        const int iq = vertex[q] * 3;
        for (int d = 0; d < 3; ++d) {
          triplets.emplace_back(ip + d, iq + d, value);
        }
      }
    }
  }

  // (3) pinned 行整行覆盖为对角 1。注意必须放在最后：
  //     setFromTriplets 会对同一位置求和，若先放 1.0 会被之前的项污染。
  //
  // "该行对角上已有多少"必须与 (1)/(2)/(2b) 真正写进去的数**逐个一致**。
  // 弯曲这一项容易配错：本次两个方向都踩了一遍，两次都表现为"pinned 顶点不够严格"
  // （少减 ⇒ 对角变成 1+k·w²；多减 ⇒ 变成 1-k·w² 且右端被"补 1-existing"一起带偏，
  // 实测 pinned 对角 -999、线性残差 0.217）。所以判据写成**共用的小函数**，
  // 而不是在两处各写一遍条件。
  //
  // `bendDiagonalAt(v)` 的定义就是"顶点 v 这一行会不会被 (2b) 写、写多少"：
  //   · v 必须是自由点（(2b) 只写自由行）；且
  //   · stencil 里至少有一个自由顶点（三个全 pin 时 (2b) 整条跳过）。
  // 满足时该行对角拿到 k·w_v²。
  //
  // **结论：v 是 pinned ⇒ 这一项恒为 0**（pinned 行 (2b) 一个字都没写），
  // 所以下面的求和实际只在"v 自由"时非零；写成一个函数只是为了让条件和 (2b) 同步、
  // 以后再加能量项时不必回来数第二遍。
  auto bendDiagonalAt = [&mesh](int v) -> Scalar {
    Scalar sum = 0.0;
    if (mesh.isPinned(v)) return sum;
    for (const auto& s : mesh.bends) {
      if (s.stiffness == 0.0) continue;
      if (!(s.a == v || s.b == v || s.c == v)) continue;
      if (mesh.isPinned(s.a) && mesh.isPinned(s.b) && mesh.isPinned(s.c)) continue;
      const int wv = (s.a == v) ? 1 : ((s.b == v) ? -2 : 1);
      sum += s.stiffness * static_cast<Scalar>(wv) * static_cast<Scalar>(wv);
    }
    return sum;
  };

  for (int v = 0; v < n; ++v) {
    if (!mesh.isPinned(v)) continue;
    // 该行对角上目前已有的贡献：惯性项 + 所有未跳过的约束（对角各 +κ_c）+ 弯曲。
    Scalar existing = mesh.masses[static_cast<std::size_t>(v)] * invDt2;
    for (const auto& e : mesh.edges) {
      if (mesh.isPinned(e.a) || mesh.isPinned(e.b)) continue;
      if (e.a == v || e.b == v) existing += e.stiffness;
    }
    existing += bendDiagonalAt(v);
    // 补上 (1 - existing)，使该对角最终恰好等于 1；非对角项因约束被跳过而保持为零。
    triplets.emplace_back(v * 3 + 0, v * 3 + 0, 1.0 - existing);
    triplets.emplace_back(v * 3 + 1, v * 3 + 1, 1.0 - existing);
    triplets.emplace_back(v * 3 + 2, v * 3 + 2, 1.0 - existing);
  }

  Eigen::SparseMatrix<Scalar> fresh(dim, dim);
  fresh.setFromTriplets(triplets.begin(), triplets.end());
  fresh.makeCompressed();

  // 调试开关：PD_DEBUG_ASSEMBLE=1 时打印三元组数量与几个关键位置。
  if (std::getenv("PD_DEBUG_ASSEMBLE") != nullptr) {
    std::fprintf(stderr, "[assemble] dim=%d triplets=%zu edges=%d bends=%d nnz=%lld\n", dim,
                 triplets.size(), mesh.edgeCount(), mesh.bendCount(), (long long)fresh.nonZeros());
    for (const auto& tr : triplets) {
      if (tr.row() == 4 || tr.col() == 4) {
        std::fprintf(stderr, "  triplet (%d,%d) = %.10g\n", tr.row(), tr.col(), tr.value());
      }
    }
  }

  // 替换数值。直接赋值，不要写成 `L = cond ? std::move(fresh) : fresh;`
  // —— 条件表达式两个分支类型不同（右值 vs 左值）会退化为右值，
  //    结果是从已被移走的对象拷贝构造，矩阵内容变成未定义（实测为全零）。
  L = fresh;
  L.makeCompressed();

  // 常向量右端（只有调用方要的时候才算；不算时连一次遍历都不发生）。
  if (bendingRhs != nullptr) assembleBendingRhs(mesh, *bendingRhs);
}

void assembleInertialRhs(const Mesh& mesh, const std::vector<Vec3>& predicted, Scalar dt,
                         Scalar /*damping*/, Eigen::VectorXd& b) {
  const int n = mesh.vertexCount();
  // damping 参数同样不在这里使用：有效质量不做任何缩放，惯性项就是 M/h²
  // （阻尼由 Integrator 的速度更新承担，见上）。
  const Scalar invDt2 = 1.0 / (dt * dt);
  b.resize(3 * n);

  // 右端 = (M/h²)·x̂      —— 只有这一项。
  //
  // **重力不在这里出现**：它已经包含在预测位置 x̂ = x + h v + h² g 里面，
  // 而且只能出现这一次。若在右端再补一项 -M g，等于把重力算了两遍，
  // 静止条件会从 κ(x-ℓ) = m g 变成 κ(x-ℓ) = 0，平衡位置随之退到 x = ℓ。
  // （本项目在排查符号问题的过程中犯过这个错，特此记录。）
  //
  // 标准 PD 的增量势能：g(x) = 1/(2h²)‖x-x̂‖²_M + Ψ(x)，
  //   ∇g = (M/h²)(x - x̂) + ∇Ψ = 0
  //   ⇒ L x = (M/h²) x̂ + Σ_c κ_c A_cᵀ d_c
  // 重力势能 -M g·x 的作用已经通过 x̂ 里的 +h²g 体现（这正是 x̂ 的来历）。
  for (int v = 0; v < n; ++v) {
    const std::size_t i = static_cast<std::size_t>(v);
    const Scalar mOverH2 = mesh.masses[i] * invDt2;
    const Vec3& p = predicted[i];
    b[v * 3 + 0] = mOverH2 * p.x;
    b[v * 3 + 1] = mOverH2 * p.y;
    b[v * 3 + 2] = mOverH2 * p.z;
  }
}

void applyPinRhs(const Mesh& mesh, Eigen::VectorXd& b) {
  for (int v = 0; v < mesh.vertexCount(); ++v) {
    if (!mesh.isPinned(v)) continue;
    const Vec3& p = mesh.pinPositions[static_cast<std::size_t>(v)];
    b[v * 3 + 0] = p.x;
    b[v * 3 + 1] = p.y;
    b[v * 3 + 2] = p.z;
  }
}

void applyPinRhsComponent(const Mesh& mesh, Eigen::VectorXd& b, int c, int components) {
  // components == 1：不分量，三个分量一起写（与 applyPinRhs 等价）。
  // components == 3：只写 3v+c —— 不同 c 写的下标互不相交 ⇒ 可并发。
  if (components <= 1) {
    applyPinRhs(mesh, b);
    return;
  }
  for (int v = 0; v < mesh.vertexCount(); ++v) {
    if (!mesh.isPinned(v)) continue;
    const Vec3& p = mesh.pinPositions[static_cast<std::size_t>(v)];
    const Scalar d[3] = {p.x, p.y, p.z};
    b[v * 3 + c] = d[c];
  }
}

void unpackPositions(const Eigen::VectorXd& x, std::vector<Vec3>& out) {
  const int n = static_cast<int>(x.size() / 3);
  out.resize(static_cast<std::size_t>(n));
  for (int v = 0; v < n; ++v) {
    out[static_cast<std::size_t>(v)] = Vec3{x[v * 3 + 0], x[v * 3 + 1], x[v * 3 + 2]};
  }
}

void packPositions(const std::vector<Vec3>& positions, Eigen::VectorXd& out) {
  const int n = static_cast<int>(positions.size());
  out.resize(3 * n);
  for (int v = 0; v < n; ++v) {
    const Vec3& p = positions[static_cast<std::size_t>(v)];
    out[v * 3 + 0] = p.x;
    out[v * 3 + 1] = p.y;
    out[v * 3 + 2] = p.z;
  }
}

}  // namespace pd
