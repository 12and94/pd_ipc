// _verify/iteration_spectrum.cpp
// 判据 A：**PD 迭代算子的谱半径 vs 实测收缩率** —— 一个不依赖任何外部文献的"实现正确性"判据。
//
// 背景（见 docs/pd-convergence.md §4.4）：跟文献比做不到闭环（他们的 κ 量纲、pin 方式、
// 网格结构、误差归一化口径都缺）。那就换个问法：**我们实测的收缩率，是否正好等于 PD 应有的那个？**
//
// 推导（线性化，围绕当前位形）：PD 的一次迭代 = 局部投影 + 全局求解
//     局部：d_c = Π(A_c x)（距离约束：d = ℓ·unit(A_c x)）
//     全局：A x⁺ = M x̂/h² + Σ_c κ_c A_cᵀ d_c ，其中 A = M/h² + Σ_c κ_c A_cᵀA_c
//   距离约束的投影在 |A_c x| ≈ ℓ 处线性化为 P_c = I − n_c n_cᵀ（n_c = unit(A_c x)）
//   ⇒ 误差 e = x − x* 满足 e⁺ = T e，
//     **T = I − A⁻¹(M/h² + N)**，  **N = Σ_c κ_c A_cᵀ (n_c n_cᵀ) A_c**，  **B := I − T = A⁻¹(M/h² + N)**
//   （惯性项确实是未缩放的 m/h² —— 见 Assembler.cpp：damping 刻意不进 L。）
//
// 于是要算的只有两件事，而它们都只用到**已经很便宜的 A⁻¹**（预分解复用）：
//   · ρ_plain = ρ(T)         —— 幂迭代（x ← T x）
//   · λ_max(B)               —— 幂迭代（x ← B x）；且 λ_min(B) = 1 − ρ(T)（T = I − B 的特征值一一对应）
//   · Chebyshev 半迭代的理论收缩率（最优多项式在 [λ_min, λ_max] 上压零）：
//        ρ_cheb = (√λ_max − √λ_min) / (√λ_max + √λ_min)
//     ⇒ 直接回答"上 Chebyshev 能到多少"，不用先实现。
//
// 同时用**真实 integrator** 测同一状态下的收缩率（stepOnce 8 次 vs 16 次，残差之比的 1/8 次方），
// 与 ρ(T) 对照：**两者一致 ⇒ 实现按 PD 应有的速率收敛**（若实测显著大于 ρ(T)，才说明实现有问题）。
//
// 用法：pd_itspectrum [--grid 100] [--steps 600] [--iters 400] [--stiffness 2305] [--dt 0.008333]
//                    [--pin top|corners] [--settle-iters 40]
// 说明：这是**诊断工具**，不进生产、不进门禁的计时项。
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <random>
#include <string>
#include <vector>

#include <Eigen/Sparse>

#include "core/assemble/Assembler.h"
#include "core/sim/Integrator.h"
#include "core/sim/Scene.h"

using namespace pd;

namespace {

using SpMat = Eigen::SparseMatrix<Scalar>;

struct EdgeGeom {
  int a = 0;
  int b = 0;
  Scalar kappa = 0.0;
  Vec3 n{};  // unit(x_a − x_b)，退化边为 0
};

void usage() {
  std::printf(
      "用法: pd_itspectrum [--grid N] [--steps N] [--iters N] [--settle-iters N]\n"
      "                    [--stiffness K] [--dt H] [--damping D] [--pin top|corners]\n"
      "  --grid N         网格边长（默认 100）\n"
      "  --steps N        先跑多少子步到准静态（默认 600；与 pd_restrace 同口径可比）\n"
      "  --iters N        幂迭代次数（默认 400）\n"
      "  --settle-iters N 稳态阶段每子步迭代数（默认 40）\n"
      "输出：λ_max(B)、λ_min(B)、ρ(T)（预测）与实测收缩率、以及 Chebyshev 的理论收缩率与所需迭代数。\n");
}

}  // namespace

int main(int argc, char** argv) {
  int grid = 100;
  int steps = 600;
  int powerIters = 400;
  int settleIters = 40;
  Scalar stiffness = 2305.0;
  Scalar dt = 1.0 / 120.0;
  Scalar damping = 0.02;
  bool pinTop = true;

  for (int i = 1; i < argc; ++i) {
    const std::string a = argv[i];
    auto next = [&](Scalar& out) { if (i + 1 < argc) out = std::atof(argv[++i]); };
    auto nextInt = [&](int& out) { if (i + 1 < argc) out = std::atoi(argv[++i]); };
    if (a == "--grid") nextInt(grid);
    else if (a == "--steps") nextInt(steps);
    else if (a == "--iters") nextInt(powerIters);
    else if (a == "--settle-iters") nextInt(settleIters);
    else if (a == "--stiffness") next(stiffness);
    else if (a == "--dt") next(dt);
    else if (a == "--damping") next(damping);
    else if (a == "--pin" && i + 1 < argc) pinTop = (std::strcmp(argv[++i], "top") == 0);
    else if (a == "--help" || a == "-h") { usage(); return 0; }
  }

  SceneConfig cfg;
  cfg.gridNx = grid;
  cfg.gridNy = grid;
  cfg.gridSpacing = 0.02;
  cfg.dt = dt;
  cfg.stiffness = stiffness;
  cfg.velocityDamping = damping;
  cfg.maxIterations = settleIters;
  cfg.relTolerance = 0.0;
  cfg.absTolerance = 0.0;
  cfg.residualTolerance = 1e-30;

  SimContext ctx = makeScene(cfg);
  if (ctx.mesh.bendCount() > 0) {
    std::fprintf(stderr, "⚠ 本工具算的是**距离约束**的投影线性化；请用无弯曲场景。\n");
    return 3;
  }
  if (pinTop) {
    for (auto& s : ctx.mesh.pinned) s = 0;
    for (int x = 0; x < grid; ++x) ctx.mesh.pinned[static_cast<std::size_t>((grid - 1) * grid + x)] = 1;
  }
  refreshPinPositions(ctx);
  ensureBuffers(ctx);

  for (int s = 0; s < steps; ++s) stepOnce(ctx);
  const Mesh& m = ctx.mesh;
  const int n = m.vertexCount();
  const int dim = 3 * n;

  // A = M/h² + ΣκAᵀA（与生产同一函数 ⇒ 矩阵定义不可能不一致）；求解器已在该 stamp 上分解过
  SpMat A;
  assembleLeftHandSide(m, dt, damping, A);
  if (!ctx.solver->analyzed()) ctx.solver->analyze(dim, A);
  ctx.solver->factorize(A);

  // 每条边的单位方向（当前位形）与权重
  std::vector<EdgeGeom> eg(m.edges.size());
  for (std::size_t c = 0; c < m.edges.size(); ++c) {
    const Edge& e = m.edges[c];
    const Vec3 rel = m.positions[static_cast<std::size_t>(e.a)] - m.positions[static_cast<std::size_t>(e.b)];
    const Scalar len = length(rel);
    eg[c].a = e.a;
    eg[c].b = e.b;
    eg[c].kappa = e.stiffness;
    eg[c].n = (len > 0.0) ? rel * (Scalar{1} / len) : Vec3{};
  }

  const Scalar invH2 = Scalar{1} / (dt * dt);
  std::vector<int> freeIdx;      // 自由顶点（幂迭代只在自由自由度上做，pin 行恒 0）
  for (int v = 0; v < n; ++v) {
    if (!m.isPinned(v)) freeIdx.push_back(v);
  }

  Eigen::VectorXd rhs(dim), sol(dim);

  // ---- 自检：本工具**全部**建立在「A⁻¹ 就是 solve()」这个假设上 ----
  // 幂迭代第一版给出 ρ(T) ≈ 9194（理论上不可能：B 的广义特征值必 ≤ 1）⇒ 先验一遍这个假设：
  // 随机 v → w = A·v → solve(w) 是否回到 v。
  {
    Eigen::VectorXd v = Eigen::VectorXd::Zero(dim);
    for (const int idx : freeIdx) {
      const std::size_t i = static_cast<std::size_t>(idx) * 3;
      for (int d = 0; d < 3; ++d) v[i + d] = 0.1 * static_cast<double>((idx * 7 + d * 13) % 17) - 0.8;
    }
    const Eigen::VectorXd wFull = A * v;                                   // 若只存下三角，这里会不等
    const Eigen::VectorXd wSym = A.selfadjointView<Eigen::Lower>() * v;
    ctx.solver->solve(wFull, sol);
    const Scalar errFull = (sol - v).norm() / std::fmax(1e-300, v.norm());
    ctx.solver->solve(wSym, sol);
    const Scalar errSym = (sol - v).norm() / std::fmax(1e-300, v.norm());
    std::printf("自检 solve()：|solve(A·v) − v|/|v| = %.3e（满阵乘）/ %.3e（selfadjointView<Lower>）\n",
                errFull, errSym);
    if (std::fmin(errFull, errSym) > 1e-6) {
      std::fprintf(stderr,
                   "⚠ solve() 没有还原 A·v ⇒「A⁻¹ = solve()」不成立，下面的谱数值都不可信。\n");
      return 4;
    }
    std::printf("      ⇒「A⁻¹ = solve()」成立（与%s存储一致）\n", (errSym < errFull) ? "下三角" : "满阵");
  }

  // r = (M/h²)x + N x
  auto applyInertiaPlusN = [&](const Eigen::VectorXd& x) {
    rhs.setZero();
    for (const int v : freeIdx) {
      const Scalar d = m.masses[static_cast<std::size_t>(v)] * invH2;
      const std::size_t i = static_cast<std::size_t>(v) * 3;
      rhs[i + 0] = d * x[i + 0];
      rhs[i + 1] = d * x[i + 1];
      rhs[i + 2] = d * x[i + 2];
    }
    for (const EdgeGeom& g : eg) {
      if (g.kappa == 0.0 || (g.n.x == 0.0 && g.n.y == 0.0 && g.n.z == 0.0)) continue;
      const std::size_t ia = static_cast<std::size_t>(g.a) * 3;
      const std::size_t ib = static_cast<std::size_t>(g.b) * 3;
      const Vec3 rel{x[ia + 0] - x[ib + 0], x[ia + 1] - x[ib + 1], x[ia + 2] - x[ib + 2]};
      const Scalar dot = g.n.x * rel.x + g.n.y * rel.y + g.n.z * rel.z;
      const Vec3 s = g.n * (g.kappa * dot);
      rhs[ia + 0] += s.x; rhs[ia + 1] += s.y; rhs[ia + 2] += s.z;
      rhs[ib + 0] -= s.x; rhs[ib + 1] -= s.y; rhs[ib + 2] -= s.z;
    }
    // ⚠️ **pin 行的右端必须清零**（踩过的坑，务必保留）：A 的 pin 行被组装成**单位行**，
    // 所以 (A⁻¹) 在 pin 行上是 1；而分子 (M/h²+N) 在 pin 行上有 ~κ 量级的耦合项
    // （来自与自由顶点相连的边）。若不清零，A⁻¹(M/h²+N) 的 pin 行元素 ~κ，
    // 幂迭代立刻给出 ρ ≈ 9194 这种**理论上不可能**的值（B 的广义特征值必 ≤ 1）。
    // 物理上：pin 自由度不是未知量，扰动恒为 0 ⇒ 右端在 pin 行就是 0
    //（真实循环里那里被 applyPinRhs 写成 pin 位置，解出来仍是 pin 位置，扰动为 0）。
    for (int v = 0; v < n; ++v) {
      if (!m.isPinned(v)) continue;
      const std::size_t i = static_cast<std::size_t>(v) * 3;
      rhs[i + 0] = rhs[i + 1] = rhs[i + 2] = 0.0;
    }
  };

  // T x = x − A⁻¹(M/h²x + N x)     B x = A⁻¹(M/h²x + N x)
  auto applyT = [&](const Eigen::VectorXd& x, Eigen::VectorXd& out) {
    applyInertiaPlusN(x);
    ctx.solver->solve(rhs, sol);
    out = x - sol;
  };
  auto applyB = [&](const Eigen::VectorXd& x, Eigen::VectorXd& out) {
    applyInertiaPlusN(x);
    ctx.solver->solve(rhs, sol);
    out = sol;
  };

  // 幂迭代（同一个随机起点 + 固定种子，便于复现）
  std::mt19937 rng(12345);
  auto powerIterate = [&](auto apply, Eigen::VectorXd& x) {
    std::vector<Scalar> ratios;
    Eigen::VectorXd y(dim);
    for (int it = 0; it < powerIters; ++it) {
      apply(x, y);
      const Scalar nrm = y.norm();
      if (nrm <= 0.0) break;
      if (it >= powerIters - 20) ratios.push_back(nrm / x.norm());
      x = y / nrm;
    }
    if (ratios.size() >= 2) {
      const Scalar last = ratios.back();
      const Scalar prev = ratios[ratios.size() - 2];
      if (std::fabs(last - prev) > 1e-4 * std::fmax(1.0, std::fabs(last))) {
        std::printf("      （提示：最后两次估计 %.6f / %.6f 还没稳住 ⇒ 再把 --iters 调大）\n", prev, last);
      }
    }
    return ratios.empty() ? Scalar{0} : ratios.back();
  };

  std::printf("=== PD 迭代算子的谱 vs 实测收缩率（_verify/iteration_spectrum.cpp）===\n");
  std::printf("网格 %dx%d（%d 顶点）  κ=%.6g  h=%.6g  阻尼 %.3g  pin=%s  稳态 %d 子步 × %d 迭代\n",
              grid, grid, n, stiffness, dt, damping, pinTop ? "整条上边" : "顶边两角", steps, settleIters);
  std::printf("自检：A 的 nnz = %lld；自由顶点 %d/%d\n", static_cast<long long>(A.nonZeros()),
              static_cast<int>(freeIdx.size()), n);

  Eigen::VectorXd x = Eigen::VectorXd::Zero(dim);
  for (const int v : freeIdx) {
    const std::size_t i = static_cast<std::size_t>(v) * 3;
    for (int d = 0; d < 3; ++d) x[i + d] = static_cast<double>(rng() % 2000) / 1000.0 - 1.0;
  }
  const Eigen::VectorXd x0 = x;
  const Scalar rhoPlain = powerIterate(applyT, x);
  std::printf("\n  幂迭代 T（%d 次）：ρ(T) = %.6f   ← PD 的**预测**收缩率\n", powerIters, rhoPlain);

  x = x0;
  const Scalar lamMaxB = powerIterate(applyB, x);
  const Scalar lamMinB = Scalar{1} - rhoPlain;   // T = I − B ⇒ 特征值一一对应
  std::printf("  幂迭代 B（%d 次）：λ_max(B) = %.6f   λ_min(B) = 1 − ρ(T) = %.6f\n", powerIters, lamMaxB,
              lamMinB);
  if (lamMaxB > Scalar{1} + 1e-6) {
    std::printf("  ⚠ λ_max(B) > 1：线性化谱不在 (0,1] 内 ⇒ 下面 Chebyshev 的界只能当参考\n");
  }

  const Scalar sq = std::sqrt(lamMinB / lamMaxB);
  const Scalar rhoCheb = (1.0 - sq) / (1.0 + sq);
  std::printf("  Chebyshev 理论收缩率 ρ_cheb = (1−√(λmin/λmax))/(1+√(λmin/λmax)) = **%.6f**\n", rhoCheb);
  auto itersTo = [](Scalar rate, Scalar thr) {
    if (rate <= 0.0 || rate >= 1.0) return -1;
    return static_cast<int>(std::ceil(std::log(thr) / std::log(rate)));
  };
  for (const Scalar thr : {1e-2, 1e-3, 1e-4}) {
    const int kp = itersTo(rhoPlain, thr);
    const int kc = itersTo(rhoCheb, thr);
    if (kp > 0 && kc > 0) {
      std::printf("      到相对 %.0e：plain ≈ %d 次，Chebyshev ≈ %d 次（%.1f× 更少）\n", thr, kp, kc,
                  static_cast<double>(kp) / static_cast<double>(kc));
    } else {
      std::printf("      到相对 %.0e：plain ≈ %s，Chebyshev ≈ %s\n", thr,
                  kp > 0 ? std::to_string(kp).c_str() : "达不到", kc > 0 ? std::to_string(kc).c_str() : "达不到");
    }
  }

  // ---- 实测对照：真实 integrator 在同一状态下的收缩率 ----
  // 用非线性残差（max 范数）即可 —— 衰减**速率**与范数无关
  const std::vector<Vec3> pos = m.positions;
  const std::vector<Vec3> vel = m.velocities;
  auto residualAfter = [&](int k) {
    ctx.mesh.positions = pos;
    ctx.mesh.velocities = vel;
    ctx.config.maxIterations = k;
    stepOnce(ctx);
    return nonlinearResidual(ctx);
  };
  const Scalar r8 = residualAfter(8);
  const Scalar r16 = residualAfter(16);
  const Scalar measured = std::pow(r16 / r8, 1.0 / 8.0);
  std::printf("\n  实测（同一状态，真实 integrator）：r(8)=%.6g → r(16)=%.6g ⇒ **%.6f / 次**\n", r8, r16,
              measured);
  const Scalar ratio = measured / rhoPlain;
  std::printf("  实测 / 预测 = %.4f ⇒ %s\n", ratio,
              (ratio < 1.15)
                  ? "**实测不超过理论预测 ⇒ 实现按 PD 应有的速率收敛（正常）**"
                  : "**实测明显大于预测 ⇒ 需要查实现（不是「PD 本来就这样」能解释的）**");
  std::printf("\n口径提醒：这是**线性化**谱（围绕当前位形），只在局部有效；\n"
              "残差的实测值受非线性/相位影响，所以比较的是**速率**而不是绝对值。\n");
  return 0;
}
