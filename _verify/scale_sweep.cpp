// _verify/scale_sweep.cpp
// 网格规模扫描：稳定性、应变、性能、分解复用是否随规模正确工作。
#include <chrono>
#include <cmath>
#include <cstdio>
#include <vector>

#include "core/sim/Integrator.h"
#include "core/sim/Scene.h"

using namespace pd;

namespace {
struct Row {
  int nx, ny, vertices, edges;
  double stepMs;
  double strainMean, strainMax;
  int factorizations;
  double solveMsTotal;
  double factorMsTotal;
};
}  // namespace

int main() {
  std::printf("=== 网格规模扫描（刚度 1e4，dt=1/120，子步 2，PD 迭代 10）===\n\n");
  std::printf("  %-10s %-9s %-9s %-11s %-11s %-11s %-8s %s\n", "网格", "顶点", "约束", "单子步ms",
              "应变均值", "应变最大", "分解次数", "判定");

  std::vector<Row> rows;
  int failures = 0;
  for (auto dims : {std::pair<int, int>{20, 20}, {40, 40}, {60, 60}, {100, 100}}) {
    SceneConfig cfg;
    cfg.gridNx = dims.first;
    cfg.gridNy = dims.second;
    cfg.gridSpacing = 0.02;
    cfg.dt = 1.0 / 120.0;
    cfg.substepsPerFrame = 2;
    cfg.stiffness = 1.0e4;
    cfg.maxIterations = 10;
    cfg.relTolerance = 1e-3;

    SimContext ctx = makeScene(cfg);
    const int steps = 400;

    // 预热
    for (int s = 0; s < 20; ++s) stepOnce(ctx);

    const auto t0 = std::chrono::steady_clock::now();
    bool finite = true;
    for (int s = 0; s < steps; ++s) {
      stepOnce(ctx);
      for (const auto& p : ctx.mesh.positions) {
        if (!std::isfinite(p.y)) { finite = false; break; }
      }
      if (!finite) break;
    }
    const double ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0)
                          .count() / steps;

    Row r;
    r.nx = dims.first; r.ny = dims.second;
    r.vertices = ctx.mesh.vertexCount();
    r.edges = ctx.mesh.edgeCount();
    r.stepMs = ms;
    r.strainMean = ctx.mesh.meanRelativeStrain();
    r.strainMax = ctx.mesh.maxRelativeStrain();
    r.factorizations = ctx.factorizeCount;
    r.solveMsTotal = ctx.solver->stats().totalSolveSeconds * 1e3;
    r.factorMsTotal = ctx.solver->stats().totalFactorizeSeconds * 1e3;
    rows.push_back(r);

    const bool ok = finite && r.factorizations == 1;
    if (!ok) ++failures;
    std::printf("  %-10s %-9d %-9d %-11.3f %-11.4f %-11.4f %-8d %s\n",
                (std::to_string(dims.first) + "x" + std::to_string(dims.second)).c_str(),
                r.vertices, r.edges, r.stepMs, r.strainMean, r.strainMax, r.factorizations,
                ok ? "OK" : "NG");
  }

  std::printf("\n=== 分解复用（每次数值分解的代价 vs 每步回代）===\n\n");
  std::printf("  %-10s %-14s %-14s %-14s %s\n", "网格", "数值分解总ms", "回代总ms", "单次回代ms", "说明");
  for (const auto& r : rows) {
    const double perSolve = r.solveMsTotal / 400.0;
    std::printf("  %-10s %-14.3f %-14.1f %-14.4f %s\n",
                (std::to_string(r.nx) + "x" + std::to_string(r.ny)).c_str(), r.factorMsTotal,
                r.solveMsTotal, perSolve,
                r.factorizations == 1 ? "只分解一次，全程复用" : "重复分解(异常)");
  }

  std::printf("\n  %s（失败 %d 项）\n", failures == 0 ? "全部稳定且复用正常" : "存在异常", failures);
  return failures == 0 ? 0 : 1;
}
