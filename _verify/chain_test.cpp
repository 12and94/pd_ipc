// _verify/chain_test.cpp
// 长链条：N 个顶点沿 x 排成一条线，相邻用距离约束连接（每个内部顶点恰好 2 条约束）。
// 这样把"多约束叠加"这一因素排除，专门检验：
//   1) 自由链（无 pin）：整体自由落体，应变应保持 ≈ 0
//   2) 悬挂链（首端 pin）：应变应 ≈ (累积重量)/(κℓ)，与 N、κ 有明确解析关系
// 解析参考（等长链、每节静止长度 ℓ、每顶点质量 m）：
//   第 k 节（从自由端数起）承受 k 个顶点的重量 ⇒ 应变 ε_k = k·m·g/(κ·ℓ)
//   总伸长 = Σ_{k=1..N-1} ε_k·ℓ = m·g·(N-1)N/(2κ)
#include <cmath>
#include <cstdio>
#include <vector>

#include "core/sim/Integrator.h"
#include "core/sim/Scene.h"

using namespace pd;

namespace {

SimContext makeChain(int n, double k, double mass, double ell, bool pinned, Scalar kd) {
  SceneConfig cfg;
  cfg.dt = 1.0 / 120.0;
  cfg.gravity = Vec3{0.0, -9.81, 0.0};
  cfg.stiffness = k;
  cfg.velocityDamping = kd;
  cfg.maxIterations = 40;
  cfg.relTolerance = 0.0;
  cfg.absTolerance = 1e-16;
  cfg.substepsPerFrame = 1;
  SimContext ctx = makeScene(cfg);

  Mesh& m = ctx.mesh;
  m.positions.clear();
  m.restPositions.clear();
  m.velocities.clear();
  m.masses.clear();
  m.pinned.clear();
  m.pinPositions.clear();
  m.edges.clear();
  m.triangleIndices.clear();

  // 链放在 y 轴上、首端在原点，第 i 个顶点在 y = i*ell
  for (int i = 0; i < n; ++i) {
    const Vec3 p{0.0, i * ell, 0.0};
    m.positions.push_back(p);
    m.restPositions.push_back(p);
    m.velocities.push_back(Vec3{});
    m.masses.push_back(mass);
    m.pinned.push_back(0);
    m.pinPositions.push_back(p);
  }
  for (int i = 0; i + 1 < n; ++i) m.edges.push_back(Edge{i, i + 1, ell, k});
  if (pinned) m.pinned[0] = 1;

  m.buildSparsityPattern();
  ensureBuffers(ctx);
  return ctx;
}

double totalLength(const Mesh& m) {
  double L = 0;
  for (const auto& e : m.edges) {
    L += length(m.positions[static_cast<std::size_t>(e.a)] -
                m.positions[static_cast<std::size_t>(e.b)]);
  }
  return L;
}

}  // namespace

int main() {
  const double ell = 0.05, mass = 0.005, g = 9.81;
  const double restTotal = ell;  // 单节

  std::printf("=== 1) 自由链（无 pin，无重力之外的力）：整体下落，应变应 ≈ 0 ===\n\n");
  std::printf("  %-7s %-9s %-14s %-14s %s\n", "N", "κ", "400步后总长", "理论总长", "判定");
  for (int n : {5, 20, 100}) {
    SimContext c = makeChain(n, 1e4, mass, ell, false, 0.0);
    for (int s = 0; s < 400; ++s) stepOnce(c);
    const double L = totalLength(c.mesh);
    const double want = (n - 1) * restTotal;
    const double err = std::fabs(L - want) / want;
    std::printf("  %-7d %-9.0e %-14.6f %-14.6f %s\n", n, 1e4, L, want,
                err < 1e-3 ? "OK(无形变)" : "NG");
  }

  std::printf("\n=== 2) 悬挂链（首端 pin）：解析 ε_k = k·m·g/(κ·ℓ)，总伸长 = m·g·N(N-1)/(2κ) ===\n\n");
  std::printf("  %-7s %-9s %-14s %-14s %-10s %s\n", "N", "κ", "400步后总长", "理论总长", "相对误差",
              "判定");
  for (int n : {5, 20, 100, 400}) {
    for (double k : {1e4, 1e5}) {
      SimContext c = makeChain(n, k, mass, ell, true, 0.5);  // 加阻尼尽快静止
      for (int s = 0; s < 8000; ++s) stepOnce(c);
      const double L = totalLength(c.mesh);
      const double stretch = mass * g * (n - 1) * n / (2.0 * k);
      const double want = (n - 1) * restTotal + stretch;
      const double rel = std::fabs(L - want) / want;
      std::printf("  %-7d %-9.0e %-14.6f %-14.6f %-10.3e %s\n", n, k, L, want, rel,
                  rel < 0.05 ? "OK" : "NG");
    }
  }

  std::printf("\n  参考：N=100 时理论总伸长 = %.6g（相对 %.3g）\n",
              mass * g * 99 * 100 / (2 * 1e4), mass * g * 99 * 100 / (2 * 1e4) / (99 * ell));
  return 0;
}
