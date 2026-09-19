// _verify/steady.cpp —— 只看多步稳态：从首步正确 → 稳态错位，问题出在哪一步
#include <cstdio>

#include "core/sim/Integrator.h"
#include "core/sim/Scene.h"
#include "core/solver/EigenDirectSolver.h"

using namespace pd;

namespace {
SimContext makeSpring(Scalar y0, Scalar h, Scalar m, Scalar ell, Scalar k, Scalar kd) {
  SimContext ctx;
  ctx.config.dt = h;
  ctx.config.gravity = Vec3{0.0, -9.81, 0.0};
  ctx.config.stiffness = k;
  ctx.config.velocityDamping = kd;
  ctx.config.maxIterations = 60;
  ctx.config.relTolerance = 0.0;
  ctx.config.absTolerance = 1e-18;
  Mesh& mesh = ctx.mesh;
  mesh.positions = {Vec3{0.0, 0.0, 0.0}, Vec3{0.0, y0, 0.0}};
  mesh.restPositions = mesh.positions;
  mesh.velocities = {Vec3{}, Vec3{}};
  mesh.masses = {m, m};
  mesh.pinned = {1, 0};
  mesh.pinPositions = mesh.positions;
  mesh.edges.push_back(Edge{0, 1, ell, k});
  mesh.buildSparsityPattern();
  ensureBuffers(ctx);
  return ctx;
}
}  // namespace

int main() {
  const Scalar h = 1.0 / 120.0, m = 0.05, ell = 0.1, k = 1.0e4, g = 9.81;
  const Scalar correct = ell + m * g / k;

  for (Scalar kd : {0.0, 0.2, 0.5}) {
    // 从"候选平衡位置"出发：若它是真平衡，就应保持不动
    SimContext ctx = makeSpring(correct, h, m, ell, k, kd);
    std::printf("阻尼 k_d = %.1f : 从 y = ℓ+mg/κ = %.12g 出发\n", kd, correct);
    Scalar drift = 0.0;
    for (int s = 0; s < 4000; ++s) {
      stepOnce(ctx);
      drift = std::max(drift, std::fabs(ctx.mesh.positions[1].y - correct));
    }
    const Scalar y = ctx.mesh.positions[1].y;
    std::printf("   4000 步后 y = %.12g   最大漂移 = %.3g   末速 = %+.3g  %s\n\n", y, drift,
                ctx.mesh.velocities[1].y, drift < 1e-9 * ell ? "[真平衡，保持不动]" : "[会漂走]");
  }

  for (Scalar kd : {0.0, 0.5}) {
    SimContext ctx = makeSpring(0.12, h, m, ell, k, kd);
    std::printf("阻尼 k_d = %.1f : 从 y = 0.12 出发\n", kd);
    for (int s = 0; s < 20000; ++s) stepOnce(ctx);
    const Scalar y = ctx.mesh.positions[1].y;
    std::printf("   → 收敛到 y = %.12g   正确值 = %.12g   差 = %+.3g  %s\n\n", y, correct,
                y - correct, std::fabs(y - correct) / ell < 1e-6 ? "[正确]" : "[错误]");
  }
  return 0;
}
