// _verify/stiffness_calibration_record.cpp
// 刚度标定（docs/plan.md §3 M1.5 欠的活）：把 κ 从"拍的值"换成有物理依据的值。
//
// 标定依据（可解析核对，不引入任何外部经验值）：
//   悬挂布料的顶端竖直边承担整条纱线的重量，其应变是
//       ε_top = (该边下面的总重量) / κ
//   总重量 = N·m·g（N = 竖直段数，m = 内部顶点质量，边法口径）
//   于是    κ = N·m·g / ε_top
//   取目标应变 ε_top = 1%（棉/帆布的弹性伸长量级）反解 κ。
//
// 为什么必须按"应变量级"标定而不是拍一个 κ：精度与代价都由 ε 决定。
//   · 物理上：ε 就是布料的柔软程度（ε 越大越软）
//   · 数值上：PD 的预条件子把横向刚度当成 κ，而真实横向刚度 ∝ 张力 = κ·ε，
//     于是预条件子**高估横向刚度约 1/ε 倍** —— 这是高刚度下外层迭代慢的机理。
//     ε 越小该失配越严重，收敛越慢。标定到物理 ε 就同时把两者摆正。
//
// 本程序输出：给定网格/间距/密度下，若干目标 ε 对应的 κ，以及各自的实测收敛表现。
//
// 说明：诊断与标定工具，不是验收程序。
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

#include "core/mesh/Mesh.h"
#include "core/sim/Integrator.h"
#include "core/sim/Scene.h"

using namespace pd;

namespace {

/// 造悬挂布料（pin 整条上边）。
SimContext makeCloth(int n, Scalar spacing, Scalar density, Scalar stiffness, int iters,
                     Scalar residualTol) {
  SceneConfig cfg;
  cfg.gridNx = cfg.gridNy = n;
  cfg.gridSpacing = spacing;
  cfg.density = density;
  cfg.stiffness = stiffness;
  cfg.maxIterations = iters;
  cfg.relTolerance = 0.0;      // 关掉位移类判据，只看残差与实测应变
  cfg.absTolerance = 0.0;
  cfg.residualTolerance = residualTol;
  SimContext ctx = makeScene(cfg);
  for (int i = 0; i < n; ++i) ctx.mesh.pinned[static_cast<std::size_t>((n - 1) * n + i)] = 1;
  refreshPinPositions(ctx);
  return ctx;
}

Scalar lowestY(const Mesh& m) {
  Scalar r = 1e300;
  for (int v = 0; v < m.vertexCount(); ++v) {
    if (m.isPinned(v)) continue;
    r = std::min(r, m.positions[static_cast<std::size_t>(v)].y);
  }
  return r;
}

}  // namespace

int main(int argc, char** argv) {
  const int n = (argc > 1) ? std::atoi(argv[1]) : 60;
  const Scalar spacing = (argc > 2) ? std::atof(argv[2]) : 0.02;
  const Scalar density = (argc > 3) ? std::atof(argv[3]) : 1.0;
  const Scalar g = 9.81;
  const int frames = 600;

  // 内部顶点质量（边法：4 条边各给 ρℓ²/2）
  const Scalar mInterior = 2.0 * density * spacing * spacing;
  const int N = n - 1;  // 竖直段数

  std::printf("=== 刚度标定记录 ===\n");
  std::printf("网格 %dx%d  间距 %.4g  密度 %.3g  布长 %.4f m\n", n, n, spacing, density,
              N * spacing);
  std::printf("内部顶点质量 m = 2ρℓ² = %.6g kg   竖直段数 N = %d\n", mInterior, N);
  std::printf("竖直链条张力上限 N·m·g = %.6g N\n\n", N * mInterior * g);

  // 实测反解：应变 ≈ A/κ（近线性已验证），先测一点定出 A，再反解目标 κ。
  auto strainAt = [&](Scalar kappa) {
    SimContext c = makeCloth(n, spacing, density, kappa, 40, 1.0e-3);
    for (int s = 0; s < frames * 2; ++s) stepOnce(c);
    return c.mesh.maxRelativeStrain();
  };
  const Scalar kProbe = 1000.0;
  const Scalar epsProbe = strainAt(kProbe);
  const Scalar A = epsProbe * kProbe;  // ε·κ 近似常数
  std::printf("实测标定：κ=%.6g 时 ε_max=%.4g  =>  ε·κ ≈ %.6g\n\n", kProbe,
              static_cast<double>(epsProbe), static_cast<double>(A));

  std::printf("  目标ε_top   反解 κ        实测应变max   实测最低y     稳态残差(m/s²)  迭代/子步\n");
  const double targets[] = {1.0e-2, 3.0e-3, 1.0e-3, 1.0e-4};
  for (double eps : targets) {
    const Scalar kappa = A / eps;
    SimContext c2 = makeCloth(n, spacing, density, kappa, 40, 1.0e-3);
    double iterSum = 0;
    for (int s = 0; s < frames * 2; ++s) iterSum += stepOnce(c2);
    std::printf("  %8.3g   %10.6g   %10.3g   %10.6f   %12.4g   %8.2f\n", eps,
                static_cast<double>(kappa), static_cast<double>(c2.mesh.maxRelativeStrain()),
                static_cast<double>(lowestY(c2.mesh)), static_cast<double>(nonlinearResidual(c2)),
                iterSum / (frames * 2.0));
  }

  std::printf("\n注 1：稳态残差用修正后的 nonlinearResidual（含当前位置重算投影）。\n");
  std::printf("   目标应落到 << g = %.2f m/s²，否则说明外层迭代仍不足。\n", g);
  std::printf("注 2：解析式 ε_top = N·m·g/κ 只计竖直链张力，实测比它大约 %.2g 倍 ——\n",
              static_cast<double>(A / (N * mInterior * g)));
  std::printf("   差值来自为维持宽度而被拉紧的横向边（纱线外张时会拉长水平约束）。\n");
  std::printf("   因此标定以**实测反解**为准，解析式只作初值估计。\n");
  return 0;
}
