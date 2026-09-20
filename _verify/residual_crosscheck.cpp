// _verify/residual_crosscheck.cpp
// 交叉复现：用**外部诊断报告的两个口径**在同一个场景上同时测量，
// 并与项目当前判据对照，定位"10 m/s²"从何而来。
//
// 两个口径（与报告 diagnose.cpp:32-34 逐字一致）：
//   residual_acc = max |g_i + (x_i - x̂_i)(m_i/h²)| / m_i        <- 本子步不动点残差
//   static_acc   = max |g_i/m_i - g_vec|                        <- 静力不平衡（含时间离散偏差）
// 其中 g_i 是**有限差分的距离约束梯度** f = d·κ(1-ℓ/r)（与 κ(Ax - d_c) 代数恒等）。
//
// 对照维度：求解设置。报告用的是旧源码 + abs=1e-6/rel=1e-5；
//           本项目现源码的判据已改为"本步位移量级"，并新增残差判据。
//
// 说明：诊断工具，不是验收程序。
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

#include "core/mesh/Mesh.h"
#include "core/sim/Integrator.h"
#include "core/sim/Scene.h"

using namespace pd;

namespace {

/// 距离约束的有限差分梯度 f = d·κ(1-ℓ/r)，按顶点累加（完全照抄报告的实现）。
std::vector<Vec3> gradFiniteDifference(const Mesh& m) {
  std::vector<Vec3> g(m.positions.size());
  for (const auto& e : m.edges) {
    const Vec3 d = m.positions[static_cast<std::size_t>(e.a)] -
                   m.positions[static_cast<std::size_t>(e.b)];
    const Scalar l = length(d);
    if (l <= 1e-12) continue;
    const Vec3 f = d * (e.stiffness * (1.0 - e.restLength / l));
    g[static_cast<std::size_t>(e.a)] += f;
    g[static_cast<std::size_t>(e.b)] -= f;
  }
  return g;
}

struct Measure {
  Scalar residualAcc = 0.0;  // 报告口径 1
  Scalar staticAcc = 0.0;    // 报告口径 2
  Scalar minY = 0.0;
  Scalar speed = 0.0;
  int used = 0;
  Scalar residualCriterion = -1.0;  // 本项目 nonlinearResidual
};

Measure snapshot(SimContext& s) {
  const Mesh& m = s.mesh;
  const auto g = gradFiniteDifference(m);
  const Scalar invH2 = 1.0 / (s.config.dt * s.config.dt);
  Measure r;
  r.minY = 1e300;
  for (int i = 0; i < m.vertexCount(); ++i) {
    const std::size_t idx = static_cast<std::size_t>(i);
    r.minY = std::min(r.minY, m.positions[idx].y);
    if (m.isPinned(i)) continue;
    const Vec3 resid =
        g[idx] + (m.positions[idx] - s.predicted[idx]) * (m.masses[idx] * invH2);
    r.residualAcc = std::max(r.residualAcc, length(resid) / m.masses[idx]);
    r.staticAcc =
        std::max(r.staticAcc, length(g[idx] / m.masses[idx] - s.config.gravity));
  }
  r.speed = maxSpeed(s);
  return r;
}

/// 完全照抄报告的 cloth()：abs=1e-6、rel=1e-5（其默认值）。
/// residualTol 单独给，用于对照"有无残差判据"。
SimContext makeClothLikeReport(int n, Scalar k, int iters, Scalar absTol, Scalar relTol,
                               Scalar residualTol) {
  SceneConfig c;
  c.gridNx = c.gridNy = n;
  c.gridSpacing = 0.02;
  c.stiffness = k;
  c.maxIterations = iters;
  c.absTolerance = absTol;
  c.relTolerance = relTol;
  c.residualTolerance = residualTol;
  SimContext s = makeScene(c);
  for (int i = 0; i < n; ++i) s.mesh.pinned[static_cast<std::size_t>((n - 1) * n + i)] = 1;
  refreshPinPositions(s);
  return s;
}

void run(const char* label, SimContext& s, int frames) {
  for (int f = 0; f < frames; ++f) stepFrame(s);
  const Measure m = snapshot(s);
  const Scalar crit = nonlinearResidual(s);
  std::printf(
      "  %-34s ymin=%9.6f  residual_acc=%11.5g  static_acc=%10.5g  本项目残差=%10.4g  "
      "提前退出=%d\n",
      label, static_cast<double>(m.minY), static_cast<double>(m.residualAcc),
      static_cast<double>(m.staticAcc), static_cast<double>(crit), s.earlyExitCount);

  // ---- 逐项对照：在同一个顶点上把两个口径的中间量都打出来 ----
  const Mesh& mesh = s.mesh;
  const auto gf = gradFiniteDifference(mesh);
  const Scalar invH2 = 1.0 / (s.config.dt * s.config.dt);
  int worst = -1;
  Scalar worstVal = -1.0;
  for (int i = 0; i < mesh.vertexCount(); ++i) {
    if (mesh.isPinned(i)) continue;
    const std::size_t idx = static_cast<std::size_t>(i);
    const Vec3 resid =
        gf[idx] + (mesh.positions[idx] - s.predicted[idx]) * (mesh.masses[idx] * invH2);
    const Scalar v = length(resid) / mesh.masses[idx];
    if (v > worstVal) {
      worstVal = v;
      worst = i;
    }
  }
  if (worst >= 0) {
    const std::size_t idx = static_cast<std::size_t>(worst);
    const Vec3 dx = mesh.positions[idx] - s.predicted[idx];
    const Vec3 inertia = dx * (mesh.masses[idx] * invH2);
    std::printf(
        "      [逐项] 顶点 %d  m=%.4g  |x-x̂|=%.6e  |g_fd|=%.6e  |inertia|=%.6e  "
        "|g+inertia|/m=%.6e\n",
        worst, static_cast<double>(mesh.masses[idx]), static_cast<double>(length(dx)),
        static_cast<double>(length(gf[idx])), static_cast<double>(length(inertia)),
        static_cast<double>(length(gf[idx] + inertia) / mesh.masses[idx]));

    // ---- 决定性对照：同一顶点上，把该顶点的"本项目口径力"也算出来 ----
    // 本项目口径：f_core = Σ_c κ_c(rel − d_c)（用 ctx.targets 的投影值）
    // 报告口径  ：f_fd   = Σ_c d·κ(1 − ℓ/‖d‖)（显式梯度）
    // 两者代数恒等，若数值不同则必有一方实现有误。
    Vec3 fCore{};
    int incidentEdges = 0;
    for (std::size_t c = 0; c < mesh.edges.size(); ++c) {
      const Edge& e = mesh.edges[c];
      if (e.a != worst && e.b != worst) continue;
      ++incidentEdges;
      const Vec3 rel = mesh.positions[static_cast<std::size_t>(e.a)] -
                       mesh.positions[static_cast<std::size_t>(e.b)];
      const Vec3 kd = s.targets[c] * e.stiffness;
      const Vec3 contrib = rel * e.stiffness - kd;
      if (e.a == worst) fCore += contrib;
      else fCore -= contrib;
    }
    std::printf("      [对照] 该顶点关联边数=%d  |f_core|=%.6e  |f_fd|=%.6e  比值=%.6g\n",
                incidentEdges, static_cast<double>(length(fCore)),
                static_cast<double>(length(gf[idx])),
                static_cast<double>(length(fCore) / std::max(1e-300, length(gf[idx]))));
    // 逐边打印：看 rel 与 d 的模长是否都等于 ℓ，以及 κ(rel−d) 与显式梯度是否一致
    for (std::size_t c = 0; c < mesh.edges.size(); ++c) {
      const Edge& e = mesh.edges[c];
      if (e.a != worst && e.b != worst) continue;
      const Vec3 rel = mesh.positions[static_cast<std::size_t>(e.a)] -
                       mesh.positions[static_cast<std::size_t>(e.b)];
      const Vec3& d = s.targets[c];
      const Vec3 explicitGrad = rel * (e.stiffness * (1.0 - e.restLength / length(rel)));
      const Vec3 coreForm = rel * e.stiffness - d * e.stiffness;
      std::printf(
          "        edge(%d,%d) ℓ=%.6f |rel|=%.9f |d|=%.9f  |κ(rel−d)|=%.6e  "
          "|κ(1−ℓ/|rel|)rel|=%.6e  差=%.3e\n",
          e.a, e.b, static_cast<double>(e.restLength), static_cast<double>(length(rel)),
          static_cast<double>(length(d)), static_cast<double>(length(coreForm)),
          static_cast<double>(length(explicitGrad)),
          static_cast<double>(length(coreForm - explicitGrad)));
    }
  }
}

}  // namespace

int main(int argc, char** argv) {
  const int n = (argc > 1) ? std::atoi(argv[1]) : 60;
  const int frames = (argc > 2) ? std::atoi(argv[2]) : 300;
  const Scalar k = (argc > 3) ? std::atof(argv[3]) : 1.0e5;

  std::printf("=== 交叉复现：10 m/s² 从何而来 ===\n");
  std::printf("场景 grid=%d kappa=%.3g spacing=0.02 h=1/120 2子步/帧 %d帧\n\n", n, k, frames);

  std::printf("--- A. 报告用的设置（旧判据语义：abs=1e-6, rel=1e-5）---\n");
  {
    SimContext s = makeClothLikeReport(n, k, 40, 1.0e-6, 1.0e-5, 0.0);
    run("abs=1e-6 rel=1e-5 无残差判据", s, frames);
  }

  std::printf("\n--- B. 只把判据换成'本步位移量级'（本项目现源码语义）---\n");
  {
    SimContext s = makeClothLikeReport(n, k, 40, 1.0e-12, 1.0e-5, 0.0);
    run("abs=1e-12 rel=1e-5 无残差判据", s, frames);
  }

  std::printf("\n--- C. 再加残差判据（当前默认）---\n");
  {
    SimContext s = makeClothLikeReport(n, k, 40, 1.0e-12, 1.0e-5, 1.0e-3);
    run("abs=1e-12 rel=1e-5 res=1e-3", s, frames);
  }

  std::printf("\n--- D. 强制跑满 40 次迭代（关掉全部提前退出）---\n");
  {
    SimContext s = makeClothLikeReport(n, k, 40, 0.0, 0.0, 0.0);
    run("全部禁用（强制 40 次）", s, frames);
  }

  return 0;
}
