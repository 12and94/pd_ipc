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
  SolverStamp s;
  s.topologyId = topologyId;
  s.dtBits = bitsOf(dt);
  s.dampingBits = bitsOf(damping);  // 阻尼进入有效质量，故也是失效条件之一

  uint64_t stiffnessHash = 0xcbf29ce484222325ULL;
  uint64_t massHash = 0xcbf29ce484222325ULL;
  for (const auto& e : mesh.edges) stiffnessHash = mix(stiffnessHash, bitsOf(e.stiffness));
  for (Scalar m : mesh.masses) massHash = mix(massHash, bitsOf(m));

  uint64_t pinHash = 0xcbf29ce484222325ULL;
  for (std::size_t v = 0; v < mesh.pinned.size(); ++v) {
    if (mesh.pinned[v]) {
      pinHash = mix(pinHash, static_cast<uint64_t>(v));
      pinHash = mix(pinHash, bitsOf(mesh.pinPositions[v].x));
      pinHash = mix(pinHash, bitsOf(mesh.pinPositions[v].y));
      pinHash = mix(pinHash, bitsOf(mesh.pinPositions[v].z));
    }
  }

  s.stiffnessBits = stiffnessHash;
  s.massBits = massHash;
  s.pinBits = pinHash;
  return s;
}

void assembleLeftHandSide(const Mesh& mesh, Scalar dt, Scalar damping,
                          Eigen::SparseMatrix<Scalar>& L) {
  const int n = mesh.vertexCount();
  const int dim = 3 * n;
  // 带阻尼的有效质量：M_γ = M/(1-k_d)。
  // 与 Integrator 里 v_{n+1} = (1-k_d)(x_{n+1}-x_n)/h 配对，
  // 才等价于"带集中阻尼的隐式欧拉"（详见 Integrator.cpp 第 5 步的说明）。
  // 惯性项用 M/h²：阻尼由 Integrator 里 v_{n+1}=(1-k_d)[v_n+Δx/h] 的 v_n 项体现，
  // 不需要（也不应该）在这里缩放有效质量。
  const Scalar invDt2 = 1.0 / (dt * dt);

  std::vector<Eigen::Triplet<Scalar>> triplets;
  triplets.reserve(static_cast<std::size_t>(n) * 3 + mesh.edgeCount() * 12 + 3);

  // (1) 惯性项 M/h²（对角），与位形无关。
  for (int v = 0; v < n; ++v) {
    const Scalar d = mesh.masses[static_cast<std::size_t>(v)] * invDt2;
    triplets.emplace_back(v * 3 + 0, v * 3 + 0, d);
    triplets.emplace_back(v * 3 + 1, v * 3 + 1, d);
    triplets.emplace_back(v * 3 + 2, v * 3 + 2, d);
  }

  // (2) 约束贡献（常数块）。
  DistanceTerm::assembleMatrix(mesh, {}, triplets);

  // (3) pinned 行整行覆盖为对角 1。注意必须放在最后：
  //     setFromTriplets 会对同一位置求和，若先放 1.0 会被之前的项污染。
  for (int v = 0; v < n; ++v) {
    if (!mesh.isPinned(v)) continue;
    // 该行对角上目前已有的贡献：惯性项 + 所有未跳过的约束（对角各 +κ_c）。
    Scalar existing = mesh.masses[static_cast<std::size_t>(v)] * invDt2;
    for (const auto& e : mesh.edges) {
      if (mesh.isPinned(e.a) || mesh.isPinned(e.b)) continue;
      if (e.a == v || e.b == v) existing += e.stiffness;
    }
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
    std::fprintf(stderr, "[assemble] dim=%d triplets=%zu edges=%d nnz=%lld\n", dim, triplets.size(),
                 mesh.edgeCount(), (long long)fresh.nonZeros());
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
}

void assembleInertialRhs(const Mesh& mesh, const std::vector<Vec3>& predicted, Scalar dt,
                         Scalar damping, Eigen::VectorXd& b) {
  const int n = mesh.vertexCount();
  // 带阻尼的有效质量：M_γ = M/(1-k_d)。
  // 与 Integrator 里 v_{n+1} = (1-k_d)(x_{n+1}-x_n)/h 配对，
  // 才等价于"带集中阻尼的隐式欧拉"（详见 Integrator.cpp 第 5 步的说明）。
  // 惯性项用 M/h²：阻尼由 Integrator 里 v_{n+1}=(1-k_d)[v_n+Δx/h] 的 v_n 项体现，
  // 不需要（也不应该）在这里缩放有效质量。
  const Scalar invDt2 = 1.0 / (dt * dt);
  b.resize(3 * n);

  // 右端 = (M_γ/h²)·x̂      —— 只有这一项。
  //
  // **重力不在这里出现**：它已经包含在预测位置 x̂ = x + h v + h² g 里面，
  // 而且只能出现这一次。若在右端再补一项 -M g，等于把重力算了两遍，
  // 静止条件会从 κ(x-ℓ) = m g 变成 κ(x-ℓ) = 0，平衡位置随之退到 x = ℓ。
  // （本项目在排查符号问题的过程中犯过这个错，特此记录。）
  //
  // 标准 PD 的增量势能：g(x) = 1/(2h²)‖x-x̂‖²_{M_γ} + Ψ(x)，
  //   ∇g = (M_γ/h²)(x - x̂) + ∇Ψ = 0
  //   ⇒ L x = (M_γ/h²) x̂ + Σ_c κ_c A_cᵀ d_c
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
