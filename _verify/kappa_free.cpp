// _verify/kappa_free.cpp
// 决定性实验：无 pin、无重力的自由漂浮布料，观察 κ 是否影响应变。
//   有 pin  → 应变与 κ 无关（锁死）
//   无 pin  → 若应变随 κ 明显下降，则问题定位在 pinned 边界处理
#include <cmath>
#include <cstdio>

#include "core/sim/Integrator.h"
#include "core/sim/Scene.h"

using namespace pd;

namespace {
SimContext makeFreeCloth(double k) {
  SceneConfig cfg;
  cfg.gridNx = 6;
  cfg.gridNy = 6;
  cfg.gridSpacing = 0.05;
  cfg.dt = 1.0 / 120.0;
  cfg.substepsPerFrame = 1;
  cfg.stiffness = k;
  cfg.maxIterations = 10;
  cfg.relTolerance = 1e-3;
  cfg.gravity = Vec3{0.0, 0.0, 0.0};  // 无重力
  SimContext ctx = makeScene(cfg);
  for (auto& p : ctx.mesh.pinned) p = 0;  // 去掉 pin
  return ctx;
}
}  // namespace

int main() {
  std::printf("=== 无 pin、无重力：给一个初始拉伸，看 κ 是否控制应变 ===\n");
  std::printf("（初始位形：把网格沿 x 方向拉伸 20%%）\n\n");
  std::printf("  %-10s %-16s %-16s %s\n", "κ", "第1步后均值", "400步后均值", "判定");

  for (double k : {1e3, 1e4, 1e5, 1e6, 1e7}) {
    SimContext ctx = makeFreeCloth(k);
    // 沿 x 拉伸 20%
    for (int v = 0; v < ctx.mesh.vertexCount(); ++v) {
      ctx.mesh.positions[static_cast<std::size_t>(v)].x *= 1.2;
    }
    stepOnce(ctx);
    const double after1 = ctx.mesh.meanRelativeStrain();
    for (int s = 0; s < 400; ++s) stepOnce(ctx);
    const double after400 = ctx.mesh.meanRelativeStrain();
    std::printf("  %-10.0e %-16.6f %-16.6f %s\n", k, after1, after400,
                after400 < 0.01 ? "收敛回原长" : "未恢复");
  }
  std::printf("\n  期望：κ 越大，400 步后应变越接近 0（弹簧把网格拉回原长）。\n");
  return 0;
}
