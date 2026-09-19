// 测量：0.09995095 是在第几次 PD 迭代 / 第几个时间步出现的
#include <cstdio>

#include "core/sim/Integrator.h"
#include "core/sim/Scene.h"
#include "core/solver/EigenDirectSolver.h"

using namespace pd;

namespace {
SimContext spring(Scalar y0, Scalar kd, int maxIters) {
  SimContext ctx;
  ctx.config.dt = 1.0 / 120.0;
  ctx.config.gravity = Vec3{0.0, -9.81, 0.0};
  ctx.config.stiffness = 1.0e4;
  ctx.config.velocityDamping = kd;
  ctx.config.maxIterations = maxIters;
  ctx.config.relTolerance = 0.0;
  ctx.config.absTolerance = 1e-18;
  Mesh& m = ctx.mesh;
  m.positions = {Vec3{0.0, 0.0, 0.0}, Vec3{0.0, y0, 0.0}};
  m.restPositions = m.positions;
  m.velocities = {Vec3{}, Vec3{}};
  m.masses = {0.05, 0.05};
  m.pinned = {1, 0};
  m.pinPositions = m.positions;
  m.edges.push_back(Edge{0, 1, 0.1, 1.0e4});
  m.buildSparsityPattern();
  ensureBuffers(ctx);
  return ctx;
}
}  // namespace

int main() {
  const Scalar tail = 0.09995095;
  std::printf("问题：%.8f 是几次迭代 / 几个时间步的结果？\n\n", tail);

  // A. 时间步扫描：maxIterations=40（每个子步迭代到不动点）
  std::printf("A. 每个子步迭代 40 次，看第几个时间步落到 0.09995095\n");
  {
    SimContext ctx = spring(0.12, 0.0, 40);
    int hitStep = -1;
    for (int s = 1; s <= 40; ++s) {
      stepOnce(ctx);
      const Scalar y = ctx.mesh.positions[1].y;
      if (s <= 6 || s == 40) {
        std::printf("   时间步 %2d: y = %.10f   迭代 %d 次\n", s, y, ctx.iterationsUsed);
      }
      if (hitStep < 0 && std::fabs(y - tail) < 5e-9) hitStep = s;
    }
    std::printf("   首次命中 0.09995095 于第 %d 个时间步\n\n", hitStep);
  }

  // B. 迭代次数扫描：单个时间步，限制迭代数 1..40，看解随迭代数怎么变
  std::printf("B. 只做 1 个时间步，限制 PD 迭代次数 1..40（起点 y=0.12）\n");
  for (int iters : {1, 2, 3, 5, 10, 20, 40}) {
    SimContext ctx = spring(0.12, 0.0, iters);
    stepOnce(ctx);
    std::printf("   迭代上限 %2d 次 -> 实际用 %2d 次, y = %.10f\n", iters, ctx.iterationsUsed,
                ctx.mesh.positions[1].y);
  }
  std::printf("\n   参考：正确平衡 0.10004905，代码稳态 0.09995095\n");
  return 0;
}
