// _verify/hang_test.cpp
// 最小定位：钉住顶边一行 + 重力，只改 κ，看
//   (a) 底部下沉量   (b) 应变   (c) 是否随时间持续增长（说明约束没在起作用）
#include <cmath>
#include <cstdio>

#include "core/sim/Integrator.h"
#include "core/sim/Scene.h"

using namespace pd;

namespace {
SimContext hang(double k, int n, double spacing, Scalar kd, int pinMode) {  // 0=两角 1=整行 2=单点
  SceneConfig cfg;
  cfg.gridNx = n;
  cfg.gridNy = n;
  cfg.gridSpacing = spacing;
  cfg.dt = 1.0 / 120.0;
  cfg.substepsPerFrame = 1;
  cfg.stiffness = k;
  cfg.velocityDamping = kd;
  cfg.maxIterations = 40;   // 每步解到不动点，排除内迭代不足
  cfg.relTolerance = 0.0;
  cfg.absTolerance = 1e-16;
  SimContext ctx = makeScene(cfg);
  if (pinMode == 2) {
    ctx.mesh.pinned[static_cast<std::size_t>((n - 1) * n + n / 2)] = 1;   // 顶边中点
  } else if (pinMode == 1) {
    for (int i = 0; i < n; ++i) ctx.mesh.pinned[static_cast<std::size_t>((n - 1) * n + i)] = 1;
  } else {
    ctx.mesh.pinned[static_cast<std::size_t>((n - 1) * n)] = 1;           // 两角
    ctx.mesh.pinned[static_cast<std::size_t>((n - 1) * n + n - 1)] = 1;
  }
  refreshPinPositions(ctx);
  return ctx;
}

double bottomY(const SimContext& c) {
  double lowest = 1e30;
  for (const auto& p : c.mesh.positions) lowest = std::min(lowest, p.y);
  return lowest;
}
}  // namespace

int main() {
  const int n = 8;
  const double spacing = 0.05;
  const double ell = spacing;                 // 静止长度
  const double hangLength = (n - 1) * ell;    // 0.35

  std::printf("=== 钉住顶边 + 重力：8x8 网格，间距 %.3f，悬挂总长 %.3f ===\n", spacing, hangLength);
  std::printf("无阻尼，每步迭代到不动点（40 次上限）\n\n");
  std::printf("  %-9s %-13s %-13s %-13s %-13s %s\n", "κ", "200步底部y", "800步底部y", "应变均值",
              "应变最大", "判定");

  for (double k : {1e3, 1e4, 1e5, 1e6, 1e7}) {
    SimContext ctx = hang(k, n, spacing, 0.0, 1);
    for (int s = 0; s < 200; ++s) stepOnce(ctx);
    const double y200 = bottomY(ctx);
    const double s200 = ctx.mesh.meanRelativeStrain();
    for (int s = 0; s < 600; ++s) stepOnce(ctx);
    const double y800 = bottomY(ctx);
    const double s800 = ctx.mesh.meanRelativeStrain();
    const bool grows = (y800 < y200 - 1e-6);  // 还在继续下沉
    std::printf("  %-9.0e %-13.5f %-13.5f %-13.5f %-13.5f %s\n", k, y200, y800,
                ctx.mesh.meanRelativeStrain(), ctx.mesh.maxRelativeStrain(),
                grows ? "仍在增长" : "已静止");
    (void)s200;
    (void)s800;
  }

  std::printf("\n  参考：顶边在 y=0（网格铺在 y=0 平面），悬挂长度 %.3f\n", hangLength);
  std::printf("  期望：κ 越大 → 下沉越少、应变越小；且 800 步时应已静止。\n");
  std::printf("\n=== 同时检查：每步实际用了多少次 PD 迭代 ===\n");
  {
    SimContext c = hang(1e4, n, spacing, 0.0, 1);
    for (int s = 0; s < 5; ++s) {
      stepOnce(c);
      std::printf("  第 %d 步：迭代 %d 次，底部 y = %.6f，应变均值 %.6f\n", s + 1,
                  c.iterationsUsed, bottomY(c), c.mesh.meanRelativeStrain());
    }
  }
  return 0;
}
