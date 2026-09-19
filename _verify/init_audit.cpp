// _verify/init_audit.cpp
// 量初始状态：质量、静止长度、真实边长、每步应变。
#include <cmath>
#include <cstdio>

#include "core/mesh/Mesh.h"
#include "core/sim/Integrator.h"
#include "core/sim/Scene.h"

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
  const Mesh& m = ctx.mesh;

  std::printf("=== 初始状态审计（8x8，间距 0.05，density 1.0）===\n\n");

  // 质量
  double mMin = 1e30, mMax = 0, mSum = 0;
  for (double mm : m.masses) { mMin = std::min(mMin, mm); mMax = std::max(mMax, mm); mSum += mm; }
  std::printf("顶点数 %d, 约束数 %d\n", m.vertexCount(), m.edgeCount());
  std::printf("质量: 单顶点 [%.6g, %.6g], 总质量 %.6g\n", mMin, mMax, mSum);
  std::printf("  => 每顶点重力 = m*g = %.6g N\n\n", mMin * 9.81);

  // 静止长度 vs 当前长度
  double rlMin = 1e30, rlMax = 0, rlSum = 0, curSum = 0;
  for (const auto& e : m.edges) {
    rlMin = std::min(rlMin, e.restLength);
    rlMax = std::max(rlMax, e.restLength);
    rlSum += e.restLength;
    curSum += length(m.positions[static_cast<std::size_t>(e.a)] -
                     m.positions[static_cast<std::size_t>(e.b)]);
  }
  const double rlMean = rlSum / m.edgeCount();
  const double curMean = curSum / m.edgeCount();
  std::printf("静止长度: [%.6g, %.6g], 均值 %.6g\n", rlMin, rlMax, rlMean);
  std::printf("初始当前长度均值 %.6g\n", curMean);
  std::printf("  => 初始应变 = %.6g  (应为 0)\n\n", std::fabs(curMean - rlMean) / rlMean);

  // 刚度与缩放
  double kMin = 1e30, kMax = 0;
  for (const auto& e : m.edges) { kMin = std::min(kMin, e.stiffness); kMax = std::max(kMax, e.stiffness); }
  std::printf("约束刚度: [%.6g, %.6g]  (配置值 %.6g)\n", kMin, kMax, cfg.stiffness);
  std::printf("  => 配置 κ 被缩放系数 %.6g 改变\n\n", kMax / cfg.stiffness);

  // 理论应变估计：一条竖直边承载多少力
  const double perVertexWeight = mMin * 9.81;
  std::printf("理论估计: 每顶点重力 %.4g N；若一条边承担约 1 个顶点的重量，\n", perVertexWeight);
  std::printf("          应变 ≈ F/(κ·ℓ) = %.4g / (%.4g × %.4g) = %.6g\n", perVertexWeight, kMax,
              rlMean, perVertexWeight / (kMax * rlMean));

  // 前 3 步的实际数据
  std::printf("\n=== 前 3 步 ===\n");
  for (int s = 0; s < 3; ++s) {
    stepOnce(ctx);
    double cs = 0;
    for (const auto& e : m.edges) {
      cs += length(m.positions[static_cast<std::size_t>(e.a)] -
                   m.positions[static_cast<std::size_t>(e.b)]) / e.restLength;
    }
    std::printf("  第 %d 步: 迭代 %d, 长度/静止长度 均值 = %.6f  (应变 %.4f)\n", s + 1,
                ctx.iterationsUsed, cs / m.edgeCount(), cs / m.edgeCount() - 1.0);
  }
  return 0;
}
