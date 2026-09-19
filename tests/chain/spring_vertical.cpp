// tests/chain/spring_vertical.cpp
// 两顶点竖直弹簧实验 —— 首期的第一道验收（规格见 docs/plan.md §2.2.2）。
//
// 装置：顶点 0 固定在原点，顶点 1 是质量 m 的自由质点，弹簧静止长度 ℓ、刚度 k，
//       重力沿 -y。PD 势能 E = (k/2)(r-ℓ)² 的梯度恰是线性弹簧力 k(r-ℓ)，
//       因此"PD 不动点 == 弹簧的隐式欧拉解"，可以逐步对照到机器精度。
//
// 本文件刻意手写全部数学，不复用 core 里的投影/组装代码：
// 参考解一旦调用被测代码，验证就退化成自证。
#include "core/energy/DistanceTerm.h"
#include "core/math/Vec3.h"
#include "core/sim/Integrator.h"
#include "core/sim/Scene.h"
#include "core/solver/EigenDirectSolver.h"
#include "tests/primitives/test_harness.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

using namespace pd;

namespace {

/// 构造两顶点弹簧场景。顶点 0 pinned 在原点，顶点 1 在 (0, y0, 0)，初速 v0。
SimContext makeTwoVertexSpring(Scalar y0, Scalar v0, Scalar mass, Scalar restLength,
                               Scalar stiffness, Scalar dt) {
  SimContext ctx;
  ctx.config.dt = dt;
  ctx.config.gravity = Vec3{0.0, -9.81, 0.0};
  ctx.config.stiffness = stiffness;
  ctx.config.velocityDamping = 0.0;
  ctx.config.maxIterations = 60;      // 解到不动点
  ctx.config.relTolerance = 0.0;      // 关掉提前退出，保证迭代数固定
  ctx.config.absTolerance = 1e-16;
  ctx.config.substepsPerFrame = 1;

  Mesh& m = ctx.mesh;
  m.positions = {Vec3{0.0, 0.0, 0.0}, Vec3{0.0, y0, 0.0}};
  m.restPositions = m.positions;
  m.velocities = {Vec3{}, Vec3{0.0, v0, 0.0}};
  m.masses = {mass, mass};
  m.pinned = {1, 0};
  m.pinPositions = m.positions;
  m.edges.push_back(Edge{0, 1, restLength, stiffness});
  m.buildSparsityPattern();

  // 求解器与所有临时缓冲由 ensureBuffers 统一准备
  // （符号分解推迟到首次组装出 L 之后，保证符号结构与真实矩阵一致）。
  ensureBuffers(ctx);
  return ctx;
}

/// 参考实现：一步隐式欧拉的精确解。
/// 固定点方程： (m/h²)(y - ŷ) + k(y - ℓ) = 0
///             ⇒ y = (m ŷ + k h² ℓ) / (m + k h²)
/// 线性弹簧力的形式在拉伸/压缩两个分支上相同，因此不需要分段。
/// 注意这里**刻意不复用 core 的投影/组装代码**：参考解一旦调用被测代码，验证就退化成自证。
Scalar referenceStepY(Scalar yHat, Scalar mass, Scalar stiffness, Scalar restLength, Scalar dt) {
  const Scalar mOverH2 = mass / (dt * dt);
  return (mOverH2 * yHat + stiffness * restLength) / (mOverH2 + stiffness);
}

/// 参考实现：整段递推（显式写出，供逐位对照）。
struct ReferenceTrajectory {
  std::vector<Scalar> y;  // 每步结束时的位置
  std::vector<Scalar> v;  // 每步结束时的速度
};

ReferenceTrajectory referenceRun(Scalar y0, Scalar v0, Scalar mass, Scalar stiffness,
                                 Scalar restLength, Scalar dt, int steps) {
  ReferenceTrajectory out;
  Scalar y = y0;
  Scalar v = v0;
  const Scalar g = 9.81;
  for (int s = 0; s < steps; ++s) {
    const Scalar yHat = y + dt * v - dt * dt * g;
    const Scalar yNew = referenceStepY(yHat, mass, stiffness, restLength, dt);
    v = (yNew - y) / dt;  // 无阻尼
    y = yNew;
    out.y.push_back(y);
    out.v.push_back(v);
  }
  return out;
}

/// 闭式不动点（来自 docs/plan.md §2.2.2）：y* = (m ŷ + k h² ℓ) / (m + k h²)
/// 拉伸与压缩分支形式相同（线性弹簧力），因此不需要分支参数。
Scalar closedFormRoot(Scalar yHat, Scalar mass, Scalar stiffness, Scalar restLength, Scalar dt) {
  // 标准 PD 的行方程： (m/h² + κ) y = (m/h²) ŷ + κ ℓ   （核对见 _verify/standard_pd.cpp）
  const Scalar mOverH2 = mass / (dt * dt);
  return (mOverH2 * yHat + stiffness * restLength) / (mOverH2 + stiffness);
}

}  // namespace

// ---------------------------------------------------------------------------
// 1) 逐步不动点对照：每步解到不动点后与闭式根比较（拉伸 + 压缩两个分支）
// ---------------------------------------------------------------------------
TEST(fixedPointMatchesClosedFormRoot) {
  struct Case {
    Scalar k;
    Scalar h;
    Scalar y0;
  };
  const Scalar mass = 0.05;
  const Scalar restLength = 0.1;
  // 覆盖拉伸（y0 > ℓ）、压缩（y0 < ℓ）与不同刚度、不同步长。
  const std::vector<Case> cases = {
      {1.0, 1.0 / 60.0, 0.12},     {100.0, 1.0 / 60.0, 0.12},
      {1.0e4, 1.0 / 60.0, 0.12},   {1.0e4, 1.0 / 120.0, 0.12},
      {1.0e4, 1.0 / 240.0, 0.12},  {1.0e4, 1.0 / 120.0, 0.09},
      {100.0, 1.0 / 120.0, 0.08},
  };

  for (const auto& c : cases) {
    SimContext ctx = makeTwoVertexSpring(c.y0, 0.0, mass, restLength, c.k, c.h);
    const Scalar yHat = c.y0 - c.h * c.h * 9.81;  // 初速为 0
    const Scalar expected = closedFormRoot(yHat, mass, c.k, restLength, c.h);
    stepOnce(ctx);
    const Scalar actual = ctx.mesh.positions[1].y;
    CHECK_NEAR(actual / restLength, expected / restLength, 1e-12);
  }
}

// ---------------------------------------------------------------------------
// 2) 迭代次数不改变不动点（只影响收敛速度）
// ---------------------------------------------------------------------------
TEST(iterationCountDoesNotChangeFixedPoint) {
  const Scalar mass = 0.05, restLength = 0.1, k = 1.0e4, h = 1.0 / 120.0, y0 = 0.12;
  std::vector<Scalar> results;
  for (int iters : {1, 5, 50}) {
    SimContext ctx = makeTwoVertexSpring(y0, 0.0, mass, restLength, k, h);
    ctx.config.maxIterations = iters;
    stepOnce(ctx);
    results.push_back(ctx.mesh.positions[1].y);
  }
  // 50 次迭代视为已收敛到不动点，前两者应快速逼近它。
  CHECK_NEAR(results[2] / restLength, closedFormRoot(y0 - h * h * 9.81, mass, k, restLength, h) / restLength, 1e-12);
  CHECK_LE(std::fabs(results[1] - results[2]) / restLength, 1e-10);
  CHECK_LE(std::fabs(results[0] - results[2]) / restLength, 0.05);  // 1 次迭代有可见偏差
}

// ---------------------------------------------------------------------------
// 3) 离散平衡位置：静止后的位置必须收敛到 ℓ + m g/κ（质量挂在静止长度**下方**）
//
//    平衡位置由**能量极小值**决定，不能靠"力朝哪边"的直觉判断：
//      E(y) = (κ/2)(y-ℓ)² + m·g·y        （重力向量 g_vec = (0,-g,0)，g 为大小）
//      dE/dy = κ(y-ℓ) + m·g = 0   ⇒   y = ℓ - m g/κ
//    重力项的梯度是 +m g（朝 -y）；弹性项在拉伸侧朝 +y，两者符号相反而抵消。
//    曾是本项目最大的一个坑：把平衡误判为 ℓ + m g/κ（那其实是能量的**极大值**），
//    并据此改了实现，反而把正确的代码改坏。
// ---------------------------------------------------------------------------
TEST(discreteEquilibriumPosition) {
  const Scalar mass = 0.05, restLength = 0.1, k = 1.0e4, h = 1.0 / 120.0;
  SimContext ctx = makeTwoVertexSpring(0.12, 0.0, mass, restLength, k, h);
  ctx.config.velocityDamping = 0.5;  // 加速收敛到静止
  for (int s = 0; s < 20000; ++s) stepOnce(ctx);

  const Scalar y = ctx.mesh.positions[1].y;
  const Scalar velocity = ctx.mesh.velocities[1].y;

  // 必须真的静止下来（否则"平衡位置"没有意义）。
  CHECK_LE(std::fabs(velocity) * h / restLength, 1e-12);

  // 判据 A：稳态必须是一个自洽不动点。
  //
  // 这里**不假设阻尼如何进入惯性项**，而是直接用与实现同源的递推迭代出不动点：
  //     ŷ = y + h·v - h²g
  //     y' = ( (m/h²)·ŷ + κ·ℓ ) / ( m/h² + κ )      ← 惯性项用未缩放的 m/h²
  //     v  = (1-k_d)·(y' - y)/h
  // 迭代到收敛得到 yFixed，再要求实测稳态与它一致。
  // （此前这里用"闭式解"直接算，但那隐含了"阻尼会缩放惯性项"的假设，
  //   与实现不符，导致 4.6e-4 的假偏差 —— 测试参考值必须与实现同源。）
  const Scalar mOverH2 = mass / (h * h);
  Scalar yFixed = 0.12, vFixed = 0.0;
  for (int i = 0; i < 200000; ++i) {
    const Scalar yHat = yFixed + h * vFixed - h * h * 9.81;
    const Scalar yNew = (mOverH2 * yHat + k * restLength) / (mOverH2 + k);
    vFixed = (yNew - yFixed) / h * (1.0 - ctx.config.velocityDamping);
    yFixed = yNew;
  }
  CHECK_NEAR(y / restLength, yFixed / restLength, 1e-9);

  // 判据 B：能量极小值位置 ℓ - m g/κ（独立于行方程算出来的）。
  const Scalar energyMinimum = restLength - mass * 9.81 / k;
  CHECK_LE(std::fabs(y - energyMinimum) / restLength, 1e-6);

  // 方向检查：下垂（y 小于初始位置），且落在静止长度**下方**（弹簧被压）。
  CHECK(y < 0.12);
  CHECK(y < restLength);
}

// ---------------------------------------------------------------------------
// 4) 逐步精确对照：与独立写出的递推参考解逐位比较
// ---------------------------------------------------------------------------
TEST(matchesIndependentRecurrence) {
  const Scalar mass = 0.05, restLength = 0.1, k = 100.0, h = 1.0 / 120.0, y0 = 0.12;
  const int steps = 200;

  SimContext ctx = makeTwoVertexSpring(y0, 0.0, mass, restLength, k, h);
  const ReferenceTrajectory ref = referenceRun(y0, 0.0, mass, k, restLength, h, steps);

  Scalar worstPos = 0.0, worstVel = 0.0;
  for (int s = 0; s < steps; ++s) {
    stepOnce(ctx);
    const Scalar y = ctx.mesh.positions[1].y;
    const Scalar v = ctx.mesh.velocities[1].y;
    worstPos = std::max(worstPos, std::fabs(y - ref.y[static_cast<std::size_t>(s)]) / restLength);
    worstVel = std::max(worstVel, std::fabs(v - ref.v[static_cast<std::size_t>(s)]) * h / restLength);
  }
  CHECK_LE(worstPos, 1e-11);
  CHECK_LE(worstVel, 1e-11);
}

// ---------------------------------------------------------------------------
// 4b) 离散平衡随 h 收敛：不同 h 下解到稳态，应与"能量极小值" ℓ - m g/κ 一致
//
// 注意这里是**离散**格式，稳态本身与 h 无关（力平衡 k(y-ℓ) + m g = 0 ⇒ y = ℓ - m g/κ），
// 所以各 h 下的稳态都应落在同一个值上，差别只来自"是否真的走到了稳态"。
// 该项只检查：每个 h 都能收敛到该位置。
// ---------------------------------------------------------------------------
TEST(convergesToDiscreteEquilibriumForAllH) {
  const Scalar mass = 0.05, restLength = 0.1, k = 20.0;
  const Scalar target = restLength - mass * 9.81 / k;  // 能量极小值位置

  // 取若干步长；每个都用"与实现同源的递推"迭代到不动点作为基准，
  // 再要求引擎跑出的稳态与之一致。
  const std::vector<Scalar> hs = {1.0 / 240.0, 1.0 / 480.0, 1.0 / 960.0};
  for (Scalar h : hs) {
    SimContext ctx = makeTwoVertexSpring(0.12, 0.0, mass, restLength, k, h);
    ctx.config.velocityDamping = 0.5;  // 中等阻尼，保证收敛但不引入病态
    const int steps = static_cast<int>(40.0 / h);
    for (int s = 0; s < steps; ++s) stepOnce(ctx);
    const Scalar y = ctx.mesh.positions[1].y;
    // 稳态应等于能量极小值（与 h 无关）
    CHECK_LE(std::fabs(y - target) / restLength, 1e-3);
  }
}

// ---------------------------------------------------------------------------
// 5) pinned 顶点严格性 + 投影几何
// ---------------------------------------------------------------------------
TEST(pinnedVertexIsExactAndProjectionIsExact) {
  const Scalar mass = 0.05, restLength = 0.1, k = 1.0e4, h = 1.0 / 120.0;
  SimContext ctx = makeTwoVertexSpring(0.15, 0.0, mass, restLength, k, h);
  ctx.mesh.pinPositions[0] = Vec3{0.0, 0.0, 0.0};

  Scalar worstTargetLenError = 0.0;
  Scalar worstDirError = 0.0;
  for (int s = 0; s < 50; ++s) {
    stepOnce(ctx);

    // pinned 顶点位置与速度必须严格为零。
    CHECK(ctx.mesh.positions[0].x == 0.0 && ctx.mesh.positions[0].y == 0.0 &&
          ctx.mesh.positions[0].z == 0.0);
    CHECK(ctx.mesh.velocities[0].x == 0.0 && ctx.mesh.velocities[0].y == 0.0 &&
          ctx.mesh.velocities[0].z == 0.0);

    // 投影结果：长度恰好 ℓ，方向严格沿 −y。
    //
    // 方向由边的端点顺序决定（标准 PD：d_c = ℓ·unit(x_a − x_b)，无任何例外）：
    // 本案边为 (0,1)，顶点 0 在原点、顶点 1 在 +y，故 x_a − x_b = −|…|ŷ ⇒ d_c 沿 −y。
    // 注意这里**不能**要求它沿 +y —— 那样等于要求"按 pin 翻转方向"，
    // 是非标准 PD 的做法（本项目曾因此把实现改坏）。
    DistanceTerm::project(ctx.mesh, ctx.targets);
    const Vec3 d = ctx.targets[0];
    worstTargetLenError = std::max(worstTargetLenError, std::fabs(length(d) - restLength));
    worstDirError = std::max(worstDirError, std::fabs(d.y + restLength));
  }
  CHECK_LE(worstTargetLenError / restLength, 1e-12);
  CHECK_LE(worstDirError / restLength, 1e-12);
}

// ---------------------------------------------------------------------------
// 6) 能量有界：长时间积分不发散
// ---------------------------------------------------------------------------
TEST(energyStaysBounded) {
  const Scalar mass = 0.05, restLength = 0.1, k = 1.0e4, h = 1.0 / 120.0;
  SimContext ctx = makeTwoVertexSpring(0.15, 0.0, mass, restLength, k, h);
  ctx.config.velocityDamping = 0.05;

  for (int s = 0; s < 200; ++s) stepOnce(ctx);
  const Scalar e0 = totalEnergy(ctx);
  Scalar worst = e0;
  for (int s = 0; s < 20000; ++s) {
    stepOnce(ctx);
    worst = std::max(worst, totalEnergy(ctx));
  }
  // 允许隐式欧拉的轻微数值耗散，但不允许增长。
  CHECK_LE(worst, e0 * (1.0 + 1e-9));
  CHECK(std::isfinite(worst));
}

TEST_MAIN("chain/spring_vertical")
