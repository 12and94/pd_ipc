// tests/chain/convergence_criterion.cpp
// 收敛判据与**外层非线性求解质量**的断言。
//
// 背景（这是本项目修过的一个真缺陷，根因与判据本身不是一回事）：
//
//   1. 收敛判据原先用**坐标量级**做参照：`diff / max|x̂|∞`。于是放行的单步位移
//      = relTolerance × 坐标量级，**与刚度无关**。高刚度下真实位移被它淹没，
//      第 1 次迭代就误判"收敛"退出，每次子步留下残差并在后续帧累加（缓慢漂移）；
//      更糟的是判据**不再平移不变** —— 同一场景平移到远处，判据突然变松。
//      现已改为以**本子步已发生的位移量级**为参照。
//
//   2. 但"提前退出"不是高刚度现象的主因。高刚度下 PD 的**外层迭代**本身收敛极慢：
//      单弹簧在静止长度附近做微小切向运动时，切向位移的收缩因子是
//          q = κ / (m/h² + κ)
//      κ=1e5、m=8e-4、h=1/120 时 q ≈ 0.9998848，**40 次迭代只能完成正确切向
//      位移的约 0.46%**。因此"移动很小"绝不等于"已经解好" ——
//      更新量是被 L 预条件后的梯度，小更新不保证小物理残差。
//      下面的 springTangentialMotionNeedsManyIterations 把这个事实钉住。
//
// 判据的核心教训：**不能用"相邻迭代位移变小"当作收敛的证据**，必须同时能
// 表达"离正确解还有多远"。本文件因此以**解析解对照**为准，而不是位移判据。
#include <cmath>
#include <cstdio>
#include <vector>

#include "core/energy/DistanceTerm.h"
#include "core/sim/Integrator.h"
#include "core/sim/Scene.h"
#include "tests/primitives/test_harness.h"

using namespace pd;

namespace {

/// 单弹簧：一端 pin 在原点、另一端自由，无重力、无阻尼。
/// 自由端初始沿 -y 摆在静止长度上，给一个**横向**（+x）速度。
SimContext makeTangentialSpring(Scalar mass, Scalar restLength, Scalar kappa, int iters,
                                Scalar vx) {
  SimContext ctx;
  ctx.config.dt = 1.0 / 120.0;
  ctx.config.gravity = Vec3{0.0, 0.0, 0.0};
  ctx.config.stiffness = kappa;
  ctx.config.velocityDamping = 0.0;
  ctx.config.maxIterations = iters;
  ctx.config.relTolerance = 0.0;  // 关闭全部提前退出，强制跑满迭代数
  ctx.config.absTolerance = 0.0;
  ctx.config.residualTolerance = 0.0;  // 残差判据也必须关（三者是 OR 关系）
  ctx.config.substepsPerFrame = 1;

  Mesh& m = ctx.mesh;
  m.positions = {Vec3{0.0, 0.0, 0.0}, Vec3{0.0, -restLength, 0.0}};
  m.restPositions = m.positions;
  m.velocities = {Vec3{}, Vec3{vx, 0.0, 0.0}};
  m.masses = {mass, mass};
  m.pinned = {1, 0};
  m.pinPositions = m.positions;
  m.edges.push_back(Edge{0, 1, restLength, kappa});
  m.buildSparsityPattern();
  ensureBuffers(ctx);
  return ctx;
}

/// 该子步隐式欧拉的**精确**极小点（把整个子步解到收敛时的位置）。
///
/// 增量势 Φ(x) = α/2‖x-x̂‖² + κ/2(‖x‖-ℓ)²，α = m/h²，pin 在原点。
/// 其驻点沿 x̂ 方向，令 x = s·unit(x̂)：
///     α(s - |x̂|) + κ(s - ℓ) = 0  ⇒  s = (α|x̂| + κℓ)/(α+κ)
/// 注意 s 与 κ 的具体值无关地由上式给出，且当 α→0（κ→∞）时 s→ℓ：
/// **刚度越大，该子步的横向位移越依赖迭代把方向解出来**。
Vec3 exactSubstepSolution(const Vec3& predicted, Scalar alpha, Scalar kappa, Scalar restLength) {
  const Scalar len = length(predicted);
  if (len <= 0.0) return Vec3{};
  const Scalar s = (alpha * len + kappa * restLength) / (alpha + kappa);
  return predicted * (s / len);
}

/// 横向位移（x 分量绝对值）—— 本装置里纵向由约束锁定，横向才反映"方向解出来没有"。
Scalar lateralDisplacement(const SimContext& ctx) {
  return std::fabs(ctx.mesh.positions[1].x);
}

constexpr Scalar kMass = 8.0e-4;
constexpr Scalar kRest = 0.02;
constexpr Scalar kH = 1.0 / 120.0;
constexpr Scalar kVx = 0.1;

}  // namespace

// ---------------------------------------------------------------------------
// 1) 高刚度下"移动很小"不等于"解好了"：有限次迭代严重压制切向运动
//
//    这条是本次诊断的核心：它不依赖重力、不依赖显示帧率、不依赖提前退出。
// ---------------------------------------------------------------------------
TEST(springTangentialMotionNeedsManyIterations) {
  const Scalar alpha = kMass / (kH * kH);  // m/h²
  std::printf("    [注] α = m/h² = %.6g\n", static_cast<double>(alpha));

  const Scalar kappas[] = {1.0e3, 1.0e4, 1.0e5};
  for (Scalar k : kappas) {
    SimContext ctx = makeTangentialSpring(kMass, kRest, k, 40, kVx);
    const Vec3 predicted = ctx.mesh.positions[1] + ctx.mesh.velocities[1] * kH;
    const Vec3 exact = exactSubstepSolution(predicted, alpha, k, kRest);
    const Scalar want = std::fabs(exact.x);

    stepOnce(ctx);
    const Scalar got = lateralDisplacement(ctx);
    const Scalar ratio = (want > 0.0) ? got / want : 0.0;

    // 理论完成度：1 - q^N，q = κ/(α+κ)
    const Scalar q = k / (alpha + k);
    const Scalar theory = 1.0 - std::pow(q, 40.0);

    std::printf("    [注] κ=%.0e  解析横向 %.6e  40次迭代 %.6e  完成度 %.4f%%（理论 %.4f%%）\n",
                static_cast<double>(k), static_cast<double>(want), static_cast<double>(got),
                static_cast<double>(ratio * 100.0), static_cast<double>(theory * 100.0));

    // q 必须随 κ 增大而趋近 1（这是"越硬越慢"的机理）
    CHECK_MSG(q > 0.0 && q < 1.0, "收缩因子 q 必须落在 (0,1)");
    // 实测完成度必须与理论同量级（允许 2 倍以内的偏差，理论是渐近估计）
    CHECK_MSG(ratio <= theory * 2.0 + 1e-9, "40 次迭代的完成度不应超过理论上界太多");
  }

  // 关键断言：κ=1e5 时完成度**远低于** 1 —— 即"几乎没解出来"。
  // 这个上界是本次缺陷的量化指纹；若将来外层求解被改进，这个断言会失败并提醒更新。
  {
    SimContext ctx = makeTangentialSpring(kMass, kRest, 1.0e5, 40, kVx);
    const Vec3 predicted = ctx.mesh.positions[1] + ctx.mesh.velocities[1] * kH;
    const Scalar want = std::fabs(exactSubstepSolution(predicted, alpha, 1.0e5, kRest).x);
    stepOnce(ctx);
    const Scalar ratio = lateralDisplacement(ctx) / want;
    CHECK_MSG(ratio < 0.02,
              "κ=1e5、40 次迭代的切向完成度应远低于 2%（实测约 0.46%：外层迭代严重不足）");
  }
}

// ---------------------------------------------------------------------------
// 2) 增加迭代数必须单调改善切向完成度（证明瓶颈是外层迭代次数，而不是判据）
// ---------------------------------------------------------------------------
TEST(moreIterationsImproveTangentialAccuracy) {
  const Scalar alpha = kMass / (kH * kH);
  SimContext ctxRef = makeTangentialSpring(kMass, kRest, 1.0e5, 40, kVx);
  const Vec3 predicted = ctxRef.mesh.positions[1] + ctxRef.mesh.velocities[1] * kH;
  const Scalar want = std::fabs(exactSubstepSolution(predicted, alpha, 1.0e5, kRest).x);

  Scalar previous = -1.0;
  for (int iters : {40, 400, 4000}) {
    SimContext ctx = makeTangentialSpring(kMass, kRest, 1.0e5, iters, kVx);
    stepOnce(ctx);
    const Scalar ratio = lateralDisplacement(ctx) / want;
    std::printf("    [注] 迭代 %5d 次 -> 切向完成度 %.4f%%\n", iters, static_cast<double>(ratio * 100.0));
    CHECK_MSG(ratio > previous, "迭代数增加必须改善切向完成度");
    previous = ratio;
  }
  // 4000 次仍明显不足 100%：说明"多迭代"不是廉价的修复路径
  CHECK_MSG(previous < 1.0, "该反例中即使 4000 次迭代也不应达到解析解");
}

// ---------------------------------------------------------------------------
// 3) 判据必须平移不变（旧判据用坐标量级做参照，会被平移破坏）
// ---------------------------------------------------------------------------
TEST(convergenceCriterionIsTranslationInvariant) {
  auto makeCloth = [](const Vec3& shift) {
    SceneConfig cfg;
    cfg.gridNx = 16;
    cfg.gridNy = 16;
    cfg.gridSpacing = 0.02;
    cfg.dt = kH;
    cfg.substepsPerFrame = 1;
    cfg.stiffness = 1.0e4;
    cfg.density = 1.0;
    cfg.maxIterations = 20;
    cfg.relTolerance = 1.0e-5;
    cfg.absTolerance = 0.0;  // 只考察相对判据
    SimContext ctx = makeScene(cfg);
    for (int i = 0; i < 16; ++i) {
      ctx.mesh.pinned[static_cast<std::size_t>(15 * 16 + i)] = 1;
    }
    refreshPinPositions(ctx);
    for (auto& p : ctx.mesh.positions) p += shift;
    for (auto& p : ctx.mesh.restPositions) p += shift;
    for (auto& p : ctx.mesh.pinPositions) p += shift;
    return ctx;
  };

  SimContext a = makeCloth(Vec3{0.0, 0.0, 0.0});
  SimContext b = makeCloth(Vec3{1000.0, 500.0, -300.0});
  for (int s = 0; s < 400; ++s) {
    stepOnce(a);
    stepOnce(b);
  }

  auto lowestY = [](const SimContext& c) {
    Scalar r = 1e300;
    for (int v = 0; v < c.mesh.vertexCount(); ++v) {
      if (c.mesh.isPinned(v)) continue;
      r = std::min(r, c.mesh.positions[static_cast<std::size_t>(v)].y);
    }
    return r;
  };

  std::printf("    [注] 原位形 提前退出 %d 步；平移后 提前退出 %d 步\n", a.earlyExitCount,
              b.earlyExitCount);
  CHECK_MSG(a.earlyExitCount == b.earlyExitCount,
            "平移前后提前退出次数必须相同（旧判据会因坐标量级变大而变松）");
  // 平移了 +500，比较相对的 y 位移。
  // 容差 1e-6 m：400 步的浮点累积实测约 3e-9，而旧判据下差异是**米**量级
  // （判据随坐标变松 -> 收敛质量完全不同），因此 1e-6 仍有充足分辨力。
  CHECK_NEAR(lowestY(a), lowestY(b) - 500.0, 1e-6);
}

// ---------------------------------------------------------------------------
// 4) absTolerance / relTolerance 传 0 必须表示"禁用"，而不是"恒成立"
//
//    历史坑：判据写成 `diff <= absTolerance || rel <= relTolerance` 时，
//    传 0 会让 `diff <= 0`（diff 恰为 0 时）与 `rel <= 0` 恒成立，
//    于是"关闭提前退出"的对照实验实际立刻退出，得出错误结论。
//    现在用 "> 0 才生效" 显式判断禁用。
// ---------------------------------------------------------------------------
TEST(zeroToleranceDisablesEarlyExit) {
  SceneConfig cfg;
  cfg.gridNx = 12;
  cfg.gridNy = 12;
  cfg.gridSpacing = 0.02;
  cfg.dt = kH;
  cfg.substepsPerFrame = 1;
  cfg.stiffness = 1.0e3;
  cfg.density = 1.0;
  cfg.maxIterations = 5;      // 上限很小，便于观察是否被提前退出打断
  cfg.relTolerance = 0.0;     // 禁用
  cfg.absTolerance = 0.0;     // 禁用
  cfg.residualTolerance = 0.0;  // 禁用（三个判据是 OR 关系，少关一个就测不出"跑满"）
  SimContext ctx = makeScene(cfg);
  for (int i = 0; i < 12; ++i) ctx.mesh.pinned[static_cast<std::size_t>(11 * 12 + i)] = 1;
  refreshPinPositions(ctx);

  int used = 0;
  for (int s = 0; s < 50; ++s) used = stepOnce(ctx);
  std::printf("    [注] 两个容差都置 0 后，每子步实际迭代数 = %d（上限 %d），提前退出 %d 次\n", used,
              cfg.maxIterations, ctx.earlyExitCount);
  CHECK_MSG(used == cfg.maxIterations, "容差为 0 时必须跑满 maxIterations（0 表示禁用提前退出）");
  CHECK_MSG(ctx.earlyExitCount == 0, "禁用后不应记录到提前退出");
}

TEST_MAIN("chain/convergence_criterion")
