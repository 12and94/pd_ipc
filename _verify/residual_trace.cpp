// _verify/residual_trace.cpp
// 把 Wang 2015（SIGGRAPH Asia，Chebyshev 加速 PD）Fig. 9 那种
// "**单个子步内：相对残差 vs 迭代次数**"曲线，在**我们自己的场景**上量出来。
//
// 为什么要它：公开数据显示未加速的 PD 用**直接法**时，桌布（10K 顶点）约 **8–10 次**迭代把相对残差
// 压到 1e-2（Fig. 9(b)）。我们此前只有"末次残差 + --iters 扫描"，且锚点是**第 1 次迭代后的 r1**，
// 而论文的纵轴是**相对该子步起始 r0** —— 两个口径不能直接比。本工具把口径对齐，并给出 r_k/r0 曲线。
//
// 口径（对齐论文）：`e^(k) = ∇E(q^k)`，即"当前迭代位形 q^k 处的不平衡力"，逐次除以**本子步的 r0**。
//   · r0    —— **手工复现"解全局步之前"的那半次迭代**（沿用 _verify/residual_audit.cpp 的公开 API 组合）：
//              x̂ = x + h·v + h²·g  →  局部投影  →  装配右端  →  此处测残差。
//   · r_k≥1 —— **用真实 integrator**：恢复准静态快照 → `maxIterations = k` → `stepOnce` →
//              对当前位形重新投影 → `nonlinearResidual`（= ∇E(q^k)）。
//   ⇒ r_k 序列来自生产循环本身，不是复刻；两者都作用在同一份准静态状态上，因此可比。
//
// 范数口径提醒：本项目 `nonlinearResidual` 是 **max 范数**（最大顶点的 |不平衡力|/m），
// 论文画的是 **2 范数**。同一迭代法的**衰减速率**可比，**绝对水平**不可比（max ≥ RMS）。
//
// 用法：
//   pd_restrace --grid 100 --iters 64 [--steps 600] [--stiffness 2305] [--dt 0.008333]
//               [--damping 0.02] [--pin top|corners]
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

#include "core/assemble/Assembler.h"
#include "core/energy/DistanceTerm.h"
#include "core/sim/Integrator.h"
#include "core/sim/Scene.h"

using namespace pd;

namespace {

struct Snapshot {
  std::vector<Vec3> positions;
  std::vector<Vec3> velocities;
};

Snapshot take(const SimContext& ctx) {
  return Snapshot{ctx.mesh.positions, ctx.mesh.velocities};
}

void restore(SimContext& ctx, const Snapshot& s) {
  ctx.mesh.positions = s.positions;
  ctx.mesh.velocities = s.velocities;
}

/// 手工复现"解全局步之前"的半个迭代，返回此刻的残差 r0。
/// 序列与 stepOnce 的前半段一致：预测 → 局部投影 → 装配 L 与 b → 覆盖 pin 行 → 测残差。
/// （本函数写 ctx.predicted / targets / L / b，**不改** m.positions | velocities。）
Scalar residualAtPredictor(SimContext& ctx) {
  const Mesh& m = ctx.mesh;
  const Scalar h = ctx.config.dt;
  for (int v = 0; v < m.vertexCount(); ++v) {
    const std::size_t i = static_cast<std::size_t>(v);
    ctx.predicted[i] = m.positions[i] + m.velocities[i] * h + ctx.config.gravity * (h * h);
  }
  DistanceTerm::project(m, ctx.targets);
  assembleLeftHandSide(m, h, ctx.config.velocityDamping, ctx.L);
  if (!ctx.solver->analyzed()) ctx.solver->analyze(3 * m.vertexCount(), ctx.L);
  ctx.solver->factorize(ctx.L);
  assembleInertialRhs(m, ctx.predicted, h, ctx.config.velocityDamping, ctx.b);
  DistanceTerm::scatterInto(m, ctx.targets, ctx.b.data(), static_cast<std::size_t>(ctx.b.size()));
  applyPinRhs(m, ctx.b);
  return nonlinearResidual(ctx);
}

// ---------------------------------------------------------------- 两种范数
//
// **为什么要自己复算**：生产的 `nonlinearResidual` 只返回 **max 范数**，而 Wang 2015 的图用的是
// **2 范数**。实测发现 max 范数被**单个顶点主导**（60×60 稳态下最坏顶点一直是 0 号——悬挂布料的
// 左下角，它的残差 ~1·|g| 且几乎不降），于是 max 范数的曲线会"假平台"，与论文不可比。
// 因此这里按 **nonlinearResidual 的同式**复算每顶点残差，同时给出 max 与 2 范数，
// 并**把复算的 max 与生产函数逐位互校**（一致才继续，见 main）——复算因此是被验证过的，不是另写一套。
struct Norms {
  Scalar maxv = 0.0;
  Scalar l2 = 0.0;
  int count = 0;
};

Norms residualNorms(const SimContext& ctx) {
  const Mesh& m = ctx.mesh;
  const int n = m.vertexCount();
  const Scalar invH2 = 1.0 / (ctx.config.dt * ctx.config.dt);
  static std::vector<Scalar> force;
  force.assign(static_cast<std::size_t>(n) * 3, Scalar{0});
  // 距离约束力：与 nonlinearResidual 逐字同形（就地重算投影 d_c，不用 ctx.targets）
  for (std::size_t c = 0; c < m.edges.size(); ++c) {
    const Edge& e = m.edges[c];
    const Vec3 rel = m.positions[static_cast<std::size_t>(e.a)] - m.positions[static_cast<std::size_t>(e.b)];
    const Scalar len = length(rel);
    const Vec3 d = (len <= DistanceTerm::kMinLength) ? Vec3{} : rel * (e.restLength / len);
    const Vec3 contrib = rel * e.stiffness - d * e.stiffness;
    const std::size_t ia = static_cast<std::size_t>(e.a) * 3;
    const std::size_t ib = static_cast<std::size_t>(e.b) * 3;
    force[ia + 0] += contrib.x; force[ia + 1] += contrib.y; force[ia + 2] += contrib.z;
    force[ib + 0] -= contrib.x; force[ib + 1] -= contrib.y; force[ib + 2] -= contrib.z;
  }
  Norms out;
  Scalar sum2 = 0.0;
  for (int v = 0; v < n; ++v) {
    if (m.isPinned(v)) continue;   // 与 nonlinearResidual 一致：pin 行不是未知量
    const std::size_t i = static_cast<std::size_t>(v);
    const Scalar mass = m.masses[i];
    if (mass <= 0.0) continue;
    const Vec3 inertia = (m.positions[i] - ctx.predicted[i]) * (mass * invH2);
    const Vec3 total{inertia.x + force[i * 3 + 0], inertia.y + force[i * 3 + 1],
                     inertia.z + force[i * 3 + 2]};
    const Scalar r = length(total) / mass;
    if (r > out.maxv) out.maxv = r;
    sum2 += r * r;
    out.count += 1;
  }
  out.l2 = std::sqrt(sum2);
  return out;
}

void usage() {
  std::printf(
      "用法: pd_restrace [--grid N] [--iters K] [--steps N] [--settle-iters N]\n"
      "                  [--stiffness K] [--dt H] [--damping D] [--pin top|corners]\n"
      "  --grid N        网格边长（默认 100 ⇒ 10000 顶点，与 Wang 2015 的桌布同量级）\n"
      "  --iters K       子步内最大迭代数（曲线画到 K，默认 64）\n"
      "  --steps N       先跑多少子步达到准静态（默认 600）\n"
      "  --settle-iters N 稳态阶段每子步迭代数（默认 40 = 查看器实时档）\n"
      "  --pin top       钉住整条上边（默认；会真正静止下来；会先清掉 makeScene 的默认两角 pin）\n"
      "  --pin corners   只钉顶边两角（makeScene 默认；布料会一直摆）\n"
      "输出：k、r_k（max 与 2 范数两种）、r_k/r0、逐步收缩率；以及到 1e-1…1e-5 需要的迭代数。\n");
}

}  // namespace

int main(int argc, char** argv) {
  int grid = 100;
  int iters = 64;
  int steps = 600;
  int settleIters = 40;   // 稳态阶段每子步跑几次迭代（取查看器"实时档"的 40，别把稳态跑得过收敛）
  Scalar stiffness = 2305.0;
  Scalar spacing = 0.02;
  Scalar density = 1.0;
  Scalar softMass = 1.0e5;  // 软 pin 的附加质量（对齐 Wang 2015 的 fixed[i]=100000）
  Scalar dt = 1.0 / 120.0;
  Scalar damping = 0.02;
  int pinMode = 1;          // 0 = 两角（硬）；1 = 整条上边（硬）；2 = 两角（软：附加质量）
  bool pinTop = true;       // 兼容旧写法 --pin top|corners

  for (int i = 1; i < argc; ++i) {
    const std::string a = argv[i];
    auto next = [&](Scalar& out) { if (i + 1 < argc) out = std::atof(argv[++i]); };
    auto nextInt = [&](int& out) { if (i + 1 < argc) out = std::atoi(argv[++i]); };
    if (a == "--grid") nextInt(grid);
    else if (a == "--iters") nextInt(iters);
    else if (a == "--steps") nextInt(steps);
    else if (a == "--settle-iters") nextInt(settleIters);
    else if (a == "--stiffness") next(stiffness);
    else if (a == "--spacing") next(spacing);
    else if (a == "--density") next(density);
    else if (a == "--soft-mass") next(softMass);
    else if (a == "--dt") next(dt);
    else if (a == "--damping") next(damping);
    else if (a == "--pin" && i + 1 < argc) {
      const char* v = argv[++i];
      if (std::strcmp(v, "top") == 0) { pinMode = 1; pinTop = true; }
      else if (std::strcmp(v, "corners") == 0) { pinMode = 0; pinTop = false; }
      else if (std::strcmp(v, "soft-corners") == 0) { pinMode = 2; pinTop = false; }
      else { std::fprintf(stderr, "--pin 只认 top|corners|soft-corners\n"); return 2; }
    }
    else if (a == "--help" || a == "-h") { usage(); return 0; }
  }
  (void)pinTop;

  SceneConfig cfg;
  cfg.gridNx = grid;
  cfg.gridNy = grid;
  cfg.gridSpacing = spacing;
  cfg.density = density;
  cfg.dt = dt;
  cfg.stiffness = stiffness;
  cfg.velocityDamping = damping;
  cfg.maxIterations = settleIters;
  cfg.relTolerance = 0.0;       // 关位移判据
  cfg.absTolerance = 0.0;
  cfg.residualTolerance = 1e-30;  // 残差判据"开着但不可达" ⇒ 一定跑满 maxIterations，且每迭代都算残差

  SimContext ctx = makeScene(cfg);
  if (ctx.mesh.bendCount() > 0) {
    std::fprintf(stderr, "⚠ 本工具复算的残差只含距离约束；请用无弯曲场景（默认即无弯曲）。\n");
    return 3;
  }
  // pin 的三种口径（为了和 Wang 2015 的代码对齐）：
  //   hard-top      整条上边，硬固定（我们默认的诊断口径）
  //   hard-corners  顶边两角，硬固定（我们 makeScene 的默认）
  //   soft-corners  **他们的做法**：只钉两角，且不写 pinned，而是给该顶点**加大质量**
  //                 （他们 `fixed[i]=100000` 进的是 c=(M[i]+fixed[i])/h²；他们还在预测里
  //                  直接跳过这些顶点，我们让预测动一点点 g·h²/1e5 ≈ 1e-8 m/子步，可忽略）
  if (pinMode == 1) {
    for (auto& s : ctx.mesh.pinned) s = 0;
    for (int x = 0; x < grid; ++x) ctx.mesh.pinned[static_cast<std::size_t>((grid - 1) * grid + x)] = 1;
  } else if (pinMode == 2) {
    // 软 pin：清掉硬标记，改用附加质量
    const int c0 = (grid - 1) * grid + 0;
    const int c1 = (grid - 1) * grid + (grid - 1);
    for (auto& s : ctx.mesh.pinned) s = 0;
    ctx.mesh.masses[static_cast<std::size_t>(c0)] += softMass;
    ctx.mesh.masses[static_cast<std::size_t>(c1)] += softMass;
  }
  refreshPinPositions(ctx);
  ensureBuffers(ctx);

  // 跑到准静态
  for (int s = 0; s < steps; ++s) stepOnce(ctx);
  const Scalar vmax = maxSpeed(ctx);

  const char* pinName = (pinMode == 1) ? "整条上边（硬）"
                        : (pinMode == 2) ? "顶边两角（软：附加质量）"
                                         : "顶边两角（硬）";
  std::printf("=== 单子步内的 残差-迭代 曲线（对齐 Wang 2015 Fig.9 的口径）===\n");
  std::printf("网格 %dx%d（%d 顶点）  间距 %.4g  密度 %.4g  κ=%.6g  h=%.6g  阻尼 %.3g  pin=%s  每子步 %d 次迭代\n",
              grid, grid, grid * grid, spacing, density, stiffness, dt, damping, pinName, settleIters);
  std::printf("先跑 %d 子步到准静态：|v|∞ = %.4g m/s（越接近 0 越接近静止）\n", steps, vmax);

  const Snapshot snap = take(ctx);
  restore(ctx, snap);
  const Scalar r0max = residualAtPredictor(ctx);
  restore(ctx, snap);
  const Norms n0 = residualNorms(ctx);

  // ---- E1：论文口径的 r0 锚点（Wang 2015 Algorithm 1: q^(0) ← s_t，即**预测子位形**）----
  // 在预测子位形处评 ∇E：把 positions 与 predicted 都设成 x̂ ⇒ 惯性项恒为 0，只剩弹性梯度。
  // 注意：上面那个 r0max("预测子处")其实是在**当前位形**上算的（只换了新预测子的惯性项），
  // 不是这个量；这里才是与论文 e(0)=∇E(s_t) 对应的锚点。
  const std::vector<Vec3> posKeep = ctx.mesh.positions;
  const std::vector<Vec3> predKeep = ctx.predicted;
  for (int v = 0; v < ctx.mesh.vertexCount(); ++v) {
    const std::size_t i = static_cast<std::size_t>(v);
    const Vec3 xhat = posKeep[i] + ctx.mesh.velocities[i] * dt + ctx.config.gravity * (dt * dt);
    ctx.mesh.positions[i] = xhat;
    ctx.predicted[i] = xhat;
  }
  const Norms nPred = residualNorms(ctx);
  ctx.mesh.positions = posKeep;
  ctx.predicted = predKeep;
  std::printf("E1 论文锚点 e0 = ∇E(预测子位形)：max %.6g / 2范数 %.6g ｜ 当前位形 r0：max %.6g / 2范数 %.6g\n",
              nPred.maxv, nPred.l2, n0.maxv, n0.l2);

  // 自校验：复算的 max 必须与生产函数一致（否则本工具的 2 范数不可信）
  const Scalar prod0 = nonlinearResidual(ctx);
  const Scalar diff0 = std::fabs(n0.maxv - prod0) / std::fmax(1.0, prod0);
  std::printf("自校验 r0：复算 max = %.9g，生产 nonlinearResidual = %.9g（相对差 %.2e）⇒ %s\n",
              n0.maxv, prod0, diff0, diff0 < 1e-12 ? "一致 ✓" : "**不一致，工具不可信**");
  if (diff0 >= 1e-12) return 4;
  // 手工路径（预测子处、解全局步之前）与"当前位形"两条 r0 口径都给出来，便于与论文对照
  std::printf("r0 两种口径：预测子处（解全局步之前）= %.6g；当前位形处 = %.6g（max） / %.6g（2 范数）\n",
              r0max, n0.maxv, n0.l2);

  std::vector<Scalar> rmax(static_cast<std::size_t>(iters) + 1, 0.0);
  std::vector<Scalar> rl2(static_cast<std::size_t>(iters) + 1, 0.0);
  rmax[0] = n0.maxv;
  rl2[0] = n0.l2;
  for (int k = 1; k <= iters; ++k) {
    restore(ctx, snap);
    ctx.config.maxIterations = k;
    stepOnce(ctx);
    const Norms nn = residualNorms(ctx);   // 与生产同式；生产值在下面交叉核对
    const Scalar prod = nonlinearResidual(ctx);
    const Scalar d = std::fabs(nn.maxv - prod) / std::fmax(1.0, prod);
    if (d >= 1e-12) {
      std::fprintf(stderr, "⚠ k=%d 复算 max 与生产不一致（%.9g vs %.9g，相对差 %.2e）\n", k, nn.maxv, prod, d);
      return 4;
    }
    rmax[static_cast<std::size_t>(k)] = nn.maxv;
    rl2[static_cast<std::size_t>(k)] = nn.l2;
  }

  std::printf("\n  k      r_k max         max/r0      r_k 2范数       2范数/r0    逐步收缩(2范数)\n");
  for (int k = 0; k <= iters; ++k) {
    const std::size_t u = static_cast<std::size_t>(k);
    const Scalar relM = rmax[u] / rmax[0];
    const Scalar relL = rl2[u] / rl2[0];
    if (k == 0) {
      std::printf("  %2d   %12.6g   %10.4g  %12.6g   %10.4g          —\n", k, rmax[u], relM, rl2[u], relL);
    } else if (k <= 12 || k % 4 == 0 || k == iters) {
      const Scalar rate = rl2[u] / rl2[u - 1];
      std::printf("  %2d   %12.6g   %10.4g  %12.6g   %10.4g      %8.4f\n", k, rmax[u], relM, rl2[u], relL,
                  rate);
    }
  }

  std::printf("\n到相对门槛需要的迭代数（**按 2 范数**，与论文同口径；r0 = %.6g）：\n", rl2[0]);
  for (const Scalar thr : {1e-1, 1e-2, 1e-3, 1e-4, 1e-5}) {
    int hit = -1;
    for (int k = 1; k <= iters; ++k) {
      if (rl2[static_cast<std::size_t>(k)] <= thr * rl2[0]) { hit = k; break; }
    }
    if (hit > 0) std::printf("  %8.0e  →  %3d 次\n", thr, hit);
    else std::printf("  %8.0e  →  %d 次内未达到\n", thr, iters);
  }
  std::printf("（同口径的 max 范数仅供对照：max 范数会被单个最坏顶点主导，见本文件顶部说明。）\n");
  std::printf("\n【E1】以**论文锚点** e0 = %.6g（2 范数）为分母的迭代数：\n", nPred.l2);
  for (const Scalar thr : {1e-1, 1e-2, 1e-3, 1e-4}) {
    int hit = -1;
    for (int k = 1; k <= iters; ++k) {
      if (rl2[static_cast<std::size_t>(k)] <= thr * nPred.l2) { hit = k; break; }
    }
    if (hit > 0) std::printf("  %8.0e  →  %3d 次\n", thr, hit);
    else std::printf("  %8.0e  →  %d 次内未达到\n", thr, iters);
  }  const Scalar geo = std::pow(rl2[static_cast<std::size_t>(iters)] / rl2[1], 1.0 / (iters - 1));
  std::printf("2 范数几何平均收缩因子（1…K）：%.6f  ← 论文就是用连续两次残差之比估 ρ 的\n", geo);
  return 0;
}
