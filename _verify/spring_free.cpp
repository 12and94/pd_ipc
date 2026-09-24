#include <cmath>
#include <cstdio>
#include "core/sim/Integrator.h"
#include "core/sim/Scene.h"
using namespace pd;
int main() {
  // 弹簧水平放在 x 轴上，左端 pin，右端自由且给一个横向偏移
  SceneConfig cfg;
  cfg.dt = 1.0/120.0;
  cfg.gravity = Vec3{0,0,0};      // 无重力
  cfg.stiffness = 1.0e4;
  cfg.velocityDamping = 0.0;
  cfg.maxIterations = 40;
  cfg.relTolerance = 0.0;
  cfg.absTolerance = 1e-16;
  cfg.substepsPerFrame = 1;
  SimContext ctx = makeScene(cfg);
  Mesh& m = ctx.mesh;
  m.positions = {Vec3{0,0,0}, Vec3{0.25,0,0}};   // 从 0.25 开始（ℓ 设为 0.1）
  m.restPositions = m.positions;
  m.velocities = {Vec3{}, Vec3{}};
  m.masses = {0.05, 0.05};
  m.pinned = {1, 0};
  m.pinPositions = m.positions;
  m.edges.clear();
  m.edges.push_back(Edge{0, 1, 0.1, 1.0e4});
  m.buildSparsityPattern();
  ensureBuffers(ctx);
  std::printf("无重力，单弹簧 ℓ=0.1，从 r=0.25 释放（无阻尼）\n");
  std::printf("  %-8s %-14s %-14s\n", "步", "长度r", "速度");
  for (int s = 1; s <= 12; ++s) {
    stepOnce(ctx);
    const double r = length(m.positions[1] - m.positions[0]);
    std::printf("  %-8d %-14.6f %-14.6f\n", s, r, m.velocities[1].x);
  }
  for (int s = 13; s <= 2000; ++s) stepOnce(ctx);
  std::printf("  2000 步后 r = %.6f（应回到 ℓ=0.1 附近）\n", length(m.positions[1]-m.positions[0]));
  return 0;
}
