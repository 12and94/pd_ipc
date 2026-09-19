// _verify/kappa_effect.cpp
// 对症下药：κ 到底有没有进入求解？用小网格直接看
//   1) 矩阵对角是否随 κ 变
//   2) 一步解是否随 κ 变
//   3) 长时稳态是否随 κ 变
#include <cmath>
#include <cstdio>

#include "core/assemble/Assembler.h"
#include "core/sim/Integrator.h"
#include "core/sim/Scene.h"

using namespace pd;

namespace {
SimContext cloth(double k, int n) {
  SceneConfig cfg;
  cfg.gridNx = n;
  cfg.gridNy = n;
  cfg.gridSpacing = 0.05;
  cfg.dt = 1.0 / 120.0;
  cfg.substepsPerFrame = 1;
  cfg.stiffness = k;
  cfg.maxIterations = 10;
  cfg.relTolerance = 1e-3;
  return makeScene(cfg);
}
}  // namespace

int main() {
  std::printf("=== κ 的影响（6x6 网格，便于观察）===\n\n");
  std::printf("  %-10s %-16s %-16s %-16s\n", "κ", "一步后最大应变", "400步后均值", "400步后最大");

  for (double k : {1e3, 1e4, 1e5, 1e6}) {
    // 一步
    SimContext c1 = cloth(k, 6);
    stepOnce(c1);
    const double s1 = c1.mesh.maxRelativeStrain();

    // 400 步
    SimContext c2 = cloth(k, 6);
    for (int s = 0; s < 400; ++s) stepOnce(c2);
    std::printf("  %-10.0e %-16.6f %-16.6f %-16.6f\n", k, s1, c2.mesh.meanRelativeStrain(),
                c2.mesh.maxRelativeStrain());
  }

  std::printf("\n=== 矩阵对角是否随 κ 变（取一个自由顶点）===\n");
  for (double k : {1e3, 1e6}) {
    SimContext c = cloth(k, 6);
    Eigen::SparseMatrix<Scalar> L;
    assembleLeftHandSide(c.mesh, c.config.dt, 0.0, L);
    // 找一个内部自由顶点（索引靠近中间）
    const int n = c.mesh.vertexCount();
    const int v = n / 2;
    std::printf("  κ=%.0e  顶点%d 有效刚度(对角项) = %.6g   （m/h² 应 = %.6g）\n", k, v,
                L.coeff(v * 3, v * 3) - c.mesh.masses[static_cast<std::size_t>(v)] *
                                            (1.0 / (c.config.dt * c.config.dt)),
                c.mesh.masses[static_cast<std::size_t>(v)] / (c.config.dt * c.config.dt));
  }
  return 0;
}
