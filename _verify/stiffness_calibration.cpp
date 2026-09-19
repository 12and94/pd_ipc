// _verify/stiffness_calibration.cpp
// 刚度标定：找一组 κ 使布料应变在小/中/大网格上都落在合理区间。
// 注意顶点质量随网格变细而减小，而 (m/h²) 与 κ 的比值决定等效软硬，
// 因此固定 κ 会在细网格上偏软 —— 这正是扫描要量化的事。
#include <cmath>
#include <cstdio>
#include <vector>

#include "core/sim/Integrator.h"
#include "core/sim/Scene.h"

using namespace pd;

int main() {
  std::printf("=== 刚度扫描：应变(均值/最大) vs κ 与网格 ===\n\n");
  const std::vector<std::pair<int, int>> grids = {{20, 20}, {40, 40}, {60, 60}};
  const std::vector<double> ks = {1e4, 1e5, 5e5, 1e6, 5e6, 1e7};

  std::printf("  %-9s %-12s %-12s %-12s %s\n", "网格", "κ", "应变均值", "应变最大", "判定");
  for (const auto& g : grids) {
    for (double k : ks) {
      SceneConfig cfg;
      cfg.gridNx = g.first;
      cfg.gridNy = g.second;
      cfg.gridSpacing = 0.02;
      cfg.dt = 1.0 / 120.0;
      cfg.substepsPerFrame = 2;
      cfg.stiffness = k;
      cfg.maxIterations = 10;
      cfg.relTolerance = 1e-3;
      SimContext ctx = makeScene(cfg);
      for (int s = 0; s < 600; ++s) stepOnce(ctx);
      const double mean = ctx.mesh.meanRelativeStrain();
      const double mx = ctx.mesh.maxRelativeStrain();
      const bool stable = std::isfinite(mean) && std::isfinite(mx);
      const bool good = stable && mean < 0.05;
      std::printf("  %-9s %-12.0e %-12.4f %-12.4f %s\n",
                  (std::to_string(g.first) + "x" + std::to_string(g.second)).c_str(), k, mean, mx,
                  !stable ? "不稳定" : (good ? "良好" : "偏软"));
    }
    std::printf("\n");
  }
  std::printf("  说明：目标是平均应变 < 5%%。\n");
  return 0;
}
