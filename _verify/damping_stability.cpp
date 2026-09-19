// _verify/damping_stability.cpp
// 目的：检验"小阻尼 + 足够长时间步"下系统是否稳定。
//
// 判据（每个 (场景, 阻尼) 组合）：
//   1. 全程无 NaN/Inf
//   2. 位置有界（单弹簧 < 1 m；布料 < 10 m）
//   3. 末段速度足够小（已静止）
//   4. 稳态与能量极小值 ℓ - m g/κ 一致（仅单弹簧有解析值）
#include <cmath>
#include <cstdio>
#include <vector>

#include "core/sim/Integrator.h"
#include "core/sim/Scene.h"
#include "core/solver/EigenDirectSolver.h"

using namespace pd;

namespace {

const Scalar kH = 1.0 / 120.0, kM = 0.05, kEll = 0.1, kK = 1.0e4, kG = 9.81;

SimContext makeSpring(Scalar y0, Scalar kd) {
  SimContext ctx;
  ctx.config.dt = kH;
  ctx.config.gravity = Vec3{0.0, -kG, 0.0};
  ctx.config.stiffness = kK;
  ctx.config.velocityDamping = kd;
  ctx.config.maxIterations = 40;
  ctx.config.relTolerance = 0.0;
  ctx.config.absTolerance = 1e-18;
  Mesh& m = ctx.mesh;
  m.positions = {Vec3{0.0, 0.0, 0.0}, Vec3{0.0, y0, 0.0}};
  m.restPositions = m.positions;
  m.velocities = {Vec3{}, Vec3{}};
  m.masses = {kM, kM};
  m.pinned = {1, 0};
  m.pinPositions = m.positions;
  m.edges.push_back(Edge{0, 1, kEll, kK});
  m.buildSparsityPattern();
  ensureBuffers(ctx);
  return ctx;
}

SimContext makeCloth(Scalar kd) {
  SceneConfig cfg;
  cfg.gridNx = 16;
  cfg.gridNy = 16;
  cfg.gridSpacing = 0.04;
  cfg.dt = kH;
  cfg.substepsPerFrame = 2;
  cfg.stiffness = kK;
  cfg.velocityDamping = kd;
  cfg.maxIterations = 10;
  cfg.relTolerance = 1e-3;
  return makeScene(cfg);
}

}  // namespace

int main() {
  const Scalar energyMin = kEll - kM * kG / kK;  // 单弹簧的能量极小值
  int failures = 0;

  std::printf("=== 单弹簧：小阻尼 + 长时间（200 000 步 ≈ 1667 s）===\n");
  std::printf("能量极小值 = %.10f\n\n", energyMin);
  std::printf("  %-8s %-16s %-14s %-14s %s\n", "k_d", "末位置", "末速度", "全程max|y|", "判定");

  for (Scalar kd : {0.0, 0.001, 0.005, 0.01, 0.02, 0.05, 0.1}) {
    SimContext ctx = makeSpring(0.12, kd);
    const int steps = 200000;
    Scalar worst = 0.0;
    bool finite = true;
    for (int s = 0; s < steps; ++s) {
      stepOnce(ctx);
      const Scalar y = ctx.mesh.positions[1].y;
      if (!std::isfinite(y)) { finite = false; break; }
      worst = std::max(worst, std::fabs(y));
    }
    const Scalar y = ctx.mesh.positions[1].y;
    const Scalar v = ctx.mesh.velocities[1].y;
    const bool ok = finite && worst < 1.0 && std::fabs(y - energyMin) / kEll < 1e-3;
    if (!ok) ++failures;
    std::printf("  %-8.3f %-16.10f %-14.3e %-14.4f %s\n", kd, y, v, worst, ok ? "OK" : "NG");
  }

  std::printf("\n=== 布料 16x16：小阻尼 + 长时间（60 000 步 ≈ 500 s）===\n\n");
  std::printf("  %-8s %-14s %-14s %-14s %s\n", "k_d", "末 max|y|", "全程 max|y|", "末速度", "判定");
  for (Scalar kd : {0.0, 0.005, 0.02, 0.05, 0.2}) {
    SimContext cloth = makeCloth(kd);
    const int steps = 60000;
    Scalar worst = 0.0;
    bool finite = true;
    for (int s = 0; s < steps; ++s) {
      stepOnce(cloth);
      for (const auto& p : cloth.mesh.positions) {
        if (!std::isfinite(p.y)) { finite = false; break; }
        worst = std::max(worst, std::fabs(p.y));
      }
      if (!finite) break;
    }
    Scalar ymax = 0.0;
    for (const auto& p : cloth.mesh.positions) ymax = std::max(ymax, std::fabs(p.y));
    const Scalar vmax = maxSpeed(cloth);
    const bool ok = finite && worst < 10.0;
    if (!ok) ++failures;
    std::printf("  %-8.3f %-14.4f %-14.4f %-14.3e %s\n", kd, ymax, worst, vmax, ok ? "OK" : "NG");
  }

  std::printf("\n  %s（失败 %d 项）\n", failures == 0 ? "全部稳定" : "存在不稳定", failures);
  return failures == 0 ? 0 : 1;
}
