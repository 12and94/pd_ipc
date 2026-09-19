// app/bench/main.cpp
// 无窗口基准：跑一个布料场景，打印每阶段耗时、迭代数、应变与能量。
//
// 用法示例：
//   pd_bench --grid 40 40 --steps 600 --stiffness 1e4 --dt 0.008333 --iters 10
//   pd_bench --grid 200 200 --steps 300 --threads 8
#include <chrono>
#include <cstdio>
#include <cstring>
#include <string>

#include "core/assemble/Assembler.h"
#include "core/math/Parallel.h"
#include "core/sim/Integrator.h"
#include "core/sim/Scene.h"

using namespace pd;

namespace {

void printUsage() {
  std::printf(
      "用法: pd_bench [选项]\n"
      "  --grid NX NY        规则网格尺寸（默认 40 40）\n"
      "  --mesh PATH         从 OBJ 载入（覆盖 --grid）\n"
      "  --spacing S         网格间距（默认 0.025）\n"
      "  --steps N           总子步数（默认 600）\n"
      "  --substeps N        每帧子步数（默认 2）\n"
      "  --dt H              子步长（默认 1/120）\n"
      "  --stiffness K       距离约束刚度（默认 1e4）\n"
      "  --density RHO       面密度（默认 1.0）\n"
      "  --damping KD        速度阻尼（默认 0.02）\n"
      "  --iters N           每子步最大 PD 迭代数（默认 10）\n"
      "  --tol T             相对收敛判据（默认 1e-3）\n"
      "  --threads N         线程数（默认自动）\n"
      "  --pin-top           pin 顶边两角（默认行为）\n"
      "  --pin-all           pin 所有顶点（用于验证自由度为 0 的退化情形）\n"
      "  --report N          每 N 步打印一次\n");
}

}  // namespace

int main(int argc, char** argv) {
  SceneConfig cfg;
  int steps = 600;
  int threads = 0;
  bool pinAll = false;

  for (int i = 1; i < argc; ++i) {
    const std::string a = argv[i];
    auto next = [&](Scalar& out) {
      if (i + 1 < argc) out = std::atof(argv[++i]);
    };
    auto nextInt = [&](int& out) {
      if (i + 1 < argc) out = std::atoi(argv[++i]);
    };
    if (a == "--grid" && i + 2 < argc) {
      cfg.gridNx = std::atoi(argv[++i]);
      cfg.gridNy = std::atoi(argv[++i]);
    } else if (a == "--mesh" && i + 1 < argc) {
      cfg.meshPath = argv[++i];
    } else if (a == "--spacing") {
      next(cfg.gridSpacing);
    } else if (a == "--steps") {
      nextInt(steps);
    } else if (a == "--substeps") {
      nextInt(cfg.substepsPerFrame);
    } else if (a == "--dt") {
      next(cfg.dt);
    } else if (a == "--stiffness") {
      next(cfg.stiffness);
    } else if (a == "--density") {
      next(cfg.density);
    } else if (a == "--damping") {
      next(cfg.velocityDamping);
    } else if (a == "--iters") {
      nextInt(cfg.maxIterations);
    } else if (a == "--tol") {
      next(cfg.relTolerance);
    } else if (a == "--threads") {
      nextInt(threads);
    } else if (a == "--pin-all") {
      pinAll = true;
    } else if (a == "--report") {
      nextInt(cfg.reportEvery);
    } else if (a == "--help" || a == "-h") {
      printUsage();
      return 0;
    }
  }

  if (threads > 0) setNumThreads(threads);

  SimContext ctx = makeScene(cfg);
  if (pinAll) {
    for (auto& p : ctx.mesh.pinned) p = 1;
    refreshPinPositions(ctx);
  }

  std::printf("=== PD 布料基准（CPU 方向 1）===\n");
  std::printf("顶点 %d  约束 %d  质量 %.6g  线程 %d\n", ctx.mesh.vertexCount(), ctx.mesh.edgeCount(),
              ctx.mesh.totalMass(), numThreads());
  std::printf("L 自由度 %d  nnz %lld\n", 3 * ctx.mesh.vertexCount(), ctx.solver->stats().nnz);
  std::printf("dt %.6g  子步/帧 %d  刚度 %.6g  PD 迭代上限 %d  容差 %.3g\n\n", cfg.dt,
              cfg.substepsPerFrame, cfg.stiffness, cfg.maxIterations, cfg.relTolerance);

  resetStageTimes(ctx);

  const auto wall0 = std::chrono::steady_clock::now();
  double worstStrain = 0.0;
  int totalIterations = 0;
  for (int s = 0; s < steps; ++s) {
    const int used = stepOnce(ctx);
    totalIterations += used;
    worstStrain = std::max(worstStrain, ctx.mesh.maxRelativeStrain());
    if (cfg.reportEvery > 0 && (s + 1) % cfg.reportEvery == 0) {
      std::printf("  step %6d  迭代 %2d  应变(均值) %.4g  能量 %.6g  速度上限 %.4g\n", s + 1, used,
                  ctx.mesh.meanRelativeStrain(), totalEnergy(ctx), maxSpeed(ctx));
    }
  }
  const double wall = std::chrono::duration<double>(std::chrono::steady_clock::now() - wall0).count();

  const auto& t = ctx.times;
  const double perStepMs = 1000.0 * wall / std::max(1, steps);
  std::printf("--- 结果 ---\n");
  std::printf("子步数 %d   总耗时 %.3f s   单子步 %.3f ms   平均迭代 %.2f\n", steps, wall, perStepMs,
              static_cast<double>(totalIterations) / std::max(1, steps));
  std::printf("最终应变: 最大 %.6g  平均 %.6g\n", ctx.mesh.maxRelativeStrain(),
              ctx.mesh.meanRelativeStrain());
  std::printf("全程最大应变 %.6g\n", worstStrain);
  std::printf("能量: 弹性 %.6g  总能量 %.6g\n", surrogateEnergy(ctx), totalEnergy(ctx));

  std::printf("\n--- 阶段耗时（累计 ms，占比）---\n");
  auto pct = [&](double v) { return wall > 0.0 ? 100.0 * v / wall : 0.0; };
  std::printf("  预测      %10.3f  %5.1f%%\n", t.predict * 1e3, pct(t.predict));
  std::printf("  局部步    %10.3f  %5.1f%%\n", t.localStep * 1e3, pct(t.localStep));
  std::printf("  散射      %10.3f  %5.1f%%\n", t.scatter * 1e3, pct(t.scatter));
  std::printf("  组装+分解 %10.3f  %5.1f%%  （只应在 stamp 变化时发生）\n", t.assemble * 1e3,
              pct(t.assemble));
  std::printf("  全局回代  %10.3f  %5.1f%%\n", t.solve * 1e3, pct(t.solve));
  std::printf("  速度更新  %10.3f  %5.1f%%\n", t.velocity * 1e3, pct(t.velocity));
  std::printf("  其它      %10.3f  %5.1f%%\n", (wall - t.predict - t.localStep - t.scatter -
                                              t.assemble - t.solve - t.velocity) * 1e3,
              pct(wall - t.predict - t.localStep - t.scatter - t.assemble - t.solve - t.velocity));

  const auto& st = ctx.solver->stats();
  std::printf("\n--- 分解复用（方向 1 的核心）---\n");
  std::printf("  符号分解 %d 次  数值分解 %d 次  回代 %d 次\n", st.analyzeCalls, st.factorizeCalls,
              st.solveCalls);
  std::printf("  符号分解耗时 %.3f ms\n", st.lastAnalyzeSeconds * 1e3);
  std::printf("  单次数值分解 %.3f ms（累计 %.3f ms）\n", st.lastFactorizeSeconds * 1e3,
              st.totalFactorizeSeconds * 1e3);
  std::printf("  单次回代 %.3f ms（累计 %.3f ms）\n", st.lastSolveSeconds * 1e3,
              st.totalSolveSeconds * 1e3);
  if (st.factorizeCalls > 0 && st.solveCalls > 0) {
    const double ratio = st.totalSolveSeconds / std::max(1e-12, st.totalFactorizeSeconds);
    std::printf("  收益提示: 本场景回代总耗时是分解的 %.2fx（分解只发生 %d 次）\n", ratio,
                st.factorizeCalls);
  }
  std::printf("  （注：Eigen 的稀疏三角求解是单线程，全局步不会随线程数加速）\n");
  return 0;
}
