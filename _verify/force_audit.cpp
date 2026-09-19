// _verify/force_audit.cpp
// 量一条典型竖直约束上的实际力，与重力对比，找 1000 倍差距的来源。
#include <cmath>
#include <cstdio>

#include "core/energy/DistanceTerm.h"
#include "core/sim/Integrator.h"
#include "core/sim/Scene.h"
#include "core/solver/IGlobalSolver.h"

using namespace pd;

int main() {
  SceneConfig cfg;
  cfg.gridNx = 8;
  cfg.gridNy = 8;
  cfg.gridSpacing = 0.05;
  cfg.dt = 1.0 / 120.0;
  cfg.stiffness = 1.0e4;
  cfg.density = 1.0;
  cfg.substepsPerFrame = 1;
  cfg.maxIterations = 40;
  cfg.relTolerance = 0.0;
  cfg.absTolerance = 1e-16;

  SimContext ctx = makeScene(cfg);
  Mesh& m = ctx.mesh;
  for (int i = 0; i < 8; ++i) m.pinned[static_cast<std::size_t>(7 * 8 + i)] = 1;
  refreshPinPositions(ctx);
  ensureBuffers(ctx);

  // 找一个自由顶点：第 4 行第 4 列（j=3, i=3）=> 索引 3*8+3 = 27
  const int v = 27;
  const double mass = m.masses[static_cast<std::size_t>(v)];
  const double weight = mass * 9.81;
  std::printf("自由顶点 %d: 质量 %.6g, 重力 %.6g N\n", v, mass, weight);
  std::printf("该顶点入射约束数: ");
  int cnt = 0;
  for (const auto& e : m.edges) if (e.a == v || e.b == v) ++cnt;
  std::printf("%d\n\n", cnt);

  // 一步之后，量"约束贡献的力"与"惯性项"
  stepOnce(ctx);

  // 用投影 + 散射单独算出约束对该顶点的贡献
  std::vector<Vec3> targets;
  DistanceTerm::project(m, targets);
  std::vector<Scalar> b(3 * m.vertexCount(), 0.0);
  DistanceTerm::scatterInto(m, targets, b.data(), b.size());
  const Vec3 scatter{b[v * 3 + 0], b[v * 3 + 1], b[v * 3 + 2]};

  const double mOverH2 = mass * (1.0 / (cfg.dt * cfg.dt));
  const Vec3 xHat = m.positions[static_cast<std::size_t>(v)] +
                    m.velocities[static_cast<std::size_t>(v)] * cfg.dt +
                    cfg.gravity * (cfg.dt * cfg.dt);
  const Vec3 inertial = xHat * mOverH2;

  std::printf("一步之后（顶点 %d 的 y 分量）：\n", v);
  std::printf("  惯性项 (m/h²)·ŷ = %.6g × %.6g = %.6g\n", mOverH2, xHat.y, inertial.y);
  std::printf("  约束散射项        = %.6g\n", scatter.y);
  std::printf("  重力等效           = -m·g = %.6g\n", -weight);
  std::printf("  散射/重力 = %.4g\n", std::fabs(scatter.y) / std::max(1e-12, weight));
  std::printf("\n");
  std::printf("  该顶点入射的约束刚度合计 = %.6g（%d 条 × κ）\n", cnt * m.edges[0].stiffness, cnt);
  std::printf("  若应变 ε，单条边力 = κ·ε·ℓ = %.6g × ε × %.4g\n", m.edges[0].stiffness,
              m.edges[0].restLength);
  std::printf("  平衡要求 Σ 边力 ≈ 重力 %.6g\n", weight);
  std::printf("  => 理论 ε ≈ %.6g\n", weight / (cnt * m.edges[0].stiffness * m.edges[0].restLength));
  return 0;
}
