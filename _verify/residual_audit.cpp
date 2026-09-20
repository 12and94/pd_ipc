// _verify/residual_audit.cpp
// 校验 nonlinearResidual 公式本身 —— 用三个**有已知答案**的位形。
//
// 背景：第一版把散射项 κd_c 当成了弹性力，实测残差 7.07e6 m/s²（而应变 1e-4 的
// 弹簧其力只有 0.2 N 量级），差了 4 个数量级。真实弹性力是
//     κ(r−ℓ)·unit = κA x − κ d_c
// 必须**两项都有**：x = x̂、v = 0 时 κAᵀA x 与 κAᵀd 各自都是 O(κℓ)，
// 只有两者之差才是 O(κ(r−ℓ))。
//
// **测试设计要点（第一版搞错过）**：不能在 stepOnce 之后测。项目自己记录过
// "单弹簧（两顶点、单约束）一次 PD 迭代即精确到达不动点"（docs/plan.md §2.2.2），
// 所以 stepOnce 之后残差必然是 0，测不出公式对错。必须**手工复现一次 PD 迭代**，
// 并在**解全局步之前**测残差 —— 那才是"离不动点还有多远"。
//
// 说明：诊断工具，不是验收程序。
#include <cmath>
#include <cstdio>

#include "core/assemble/Assembler.h"
#include "core/energy/DistanceTerm.h"
#include "core/mesh/Mesh.h"
#include "core/sim/Integrator.h"
#include "core/sim/Scene.h"
#include "core/solver/EigenDirectSolver.h"

using namespace pd;

namespace {

int gFail = 0;

void report(const char* name, Scalar got, Scalar want, Scalar tol) {
  const bool ok = std::fabs(got - want) <= tol;
  if (!ok) ++gFail;
  std::printf("  %-54s 实测 %.6e  期望 %.6e  %s\n", name, static_cast<double>(got),
              static_cast<double>(want), ok ? "OK" : "NG");
}

/// 单弹簧竖直：pin 在上（原点），自由端在 y = -len 处。
SimContext makeVerticalSpring(Scalar len, Scalar mass, Scalar rest, Scalar kappa, Scalar h,
                              Scalar g) {
  SimContext ctx;
  ctx.config.dt = h;
  ctx.config.gravity = Vec3{0.0, -g, 0.0};
  ctx.config.velocityDamping = 0.0;
  ctx.config.stiffness = kappa;
  ctx.config.maxIterations = 1;
  ctx.config.relTolerance = 0.0;
  ctx.config.absTolerance = 0.0;
  ctx.config.residualTolerance = 0.0;
  Mesh& m = ctx.mesh;
  m.positions = {Vec3{0.0, 0.0, 0.0}, Vec3{0.0, -len, 0.0}};
  m.restPositions = m.positions;
  m.velocities = {Vec3{}, Vec3{}};
  m.masses = {mass, mass};
  m.pinned = {1, 0};
  m.pinPositions = m.positions;
  m.edges.push_back(Edge{0, 1, rest, kappa});
  m.buildSparsityPattern();
  ensureBuffers(ctx);
  return ctx;
}

/// 手工复现"一次 PD 迭代"到**解全局步之前**，返回此时的最大残差。
///
/// 步骤与 stepOnce 的 4a–4c 完全一致：
///   预测 x̂ -> 局部步 d=Π(Ax) -> b=(M/h²)x̂+散射 -> 覆盖 pin 行 -> 此处测残差
/// 注意：残差正是 (★) 式的左端，而 (★★) 式的右端 b 与 L x 之差就是它。
Scalar residualAfterFirstHalfIteration(SimContext& ctx) {
  const Mesh& m = ctx.mesh;
  const Scalar h = ctx.config.dt;
  const Scalar n = static_cast<Scalar>(m.vertexCount());

  // 预测
  for (int v = 0; v < m.vertexCount(); ++v) {
    const std::size_t i = static_cast<std::size_t>(v);
    ctx.predicted[i] = m.positions[i] + m.velocities[i] * h + ctx.config.gravity * (h * h);
  }
  (void)n;

  // 局部步
  DistanceTerm::project(m, ctx.targets);

  // 组装 L 与右端（首轮必然要组装）
  assembleLeftHandSide(m, h, ctx.config.velocityDamping, ctx.L);
  if (!ctx.solver->analyzed()) ctx.solver->analyze(3 * m.vertexCount(), ctx.L);
  ctx.solver->factorize(ctx.L);
  assembleInertialRhs(m, ctx.predicted, h, ctx.config.velocityDamping, ctx.b);
  DistanceTerm::scatterInto(m, ctx.targets, ctx.b.data(), static_cast<std::size_t>(ctx.b.size()));
  applyPinRhs(m, ctx.b);

  // 解全局步**之前**测残差 —— 这就是"离不动点还有多远"
  return nonlinearResidual(ctx);
}

}  // namespace

int main() {
  const Scalar h = 1.0 / 120.0, g = 9.81;
  const Scalar mass = 8.0e-4, rest = 0.02, kappa = 1.0e4;
  std::printf("=== nonlinearResidual 公式自校验（h=%.6g g=%.2f m=%.3g l=%.3g k=%.3g）===\n", h, g,
              mass, rest, kappa);
  std::printf("    （在全局步**之前**测，此时位形尚未被求解修正）\n\n");

  // ---- A. 处于解析平衡长度 ℓ + mg/κ ----
  //     位形已是本子步的不动点（v=0 且力平衡），残差必须 ≈ 0。
  {
    const Scalar eqLen = rest + mass * g / kappa;
    SimContext ctx = makeVerticalSpring(eqLen, mass, rest, kappa, h, g);
    report("A 解析平衡长度 ℓ+mg/κ：残差应为 0", residualAfterFirstHalfIteration(ctx), 0.0, 1e-6);
  }

  // ---- B. 处于静止长度 ℓ（比平衡短 mg/κ）：弹簧无力，不平衡力就是重力 ----
  {
    SimContext ctx = makeVerticalSpring(rest, mass, rest, kappa, h, g);
    report("B 静止长度 ℓ（弹簧无力）：残差应为 g", residualAfterFirstHalfIteration(ctx), g, 1e-3);
  }

  // ---- C. 拉长 10%：弹性力朝上 = κ·0.1ℓ/m，与重力反向，残差为其差 ----
  {
    const Scalar len = rest * 1.1;
    SimContext ctx = makeVerticalSpring(len, mass, rest, kappa, h, g);
    const Scalar springAcc = kappa * (len - rest) / mass;
    report("C 拉长 10%：残差应为 |κ·0.1ℓ/m − g|", residualAfterFirstHalfIteration(ctx),
           std::fabs(springAcc - g), 1e-3);
  }

  std::printf("\n  %s（失败 %d 项）\n", gFail == 0 ? "全部正确" : "存在错误", gFail);
  return gFail == 0 ? 0 : 1;
}
