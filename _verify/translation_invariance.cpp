// _verify/translation_invariance.cpp
// 最有针对性的判据：把整个系统平移一段距离，运动应当只是同样的平移。
//
// 若 pin 消元的右端项漏了 κ·q（q 为 pinned 顶点位置），则平移后结果会变差，
// 且 pin 落在原点时错误恰好被掩盖（κ·q = 0）。
//
// 检查内容：
//   1) 单弹簧：pin 在原点 vs pin 平移 1.0 —— 相对位移必须完全相同
//   2) 无重力 + 初始正好在静止长度：任何位置都必须保持不动（严格 0 位移）
//   3) 布料：整体平移后，相对形状（顶点间距离）必须与平移前一致
#include <cmath>
#include <cstdio>
#include <vector>

#include "core/sim/Integrator.h"
#include "core/sim/Scene.h"
#include "core/solver/EigenDirectSolver.h"

using namespace pd;

namespace {

const Scalar kEll = 0.1, kK = 1.0e4, kM = 0.05;

/// 单弹簧：顶点 0 是 pin，放在 pinPos；顶点 1 自由，放在 freePos
SimContext spring(const Vec3& pinPos, const Vec3& freePos, const Vec3& gravity) {
  SimContext ctx;
  ctx.config.dt = 1.0 / 120.0;
  ctx.config.gravity = gravity;
  ctx.config.stiffness = kK;
  ctx.config.velocityDamping = 0.0;
  ctx.config.maxIterations = 40;
  ctx.config.relTolerance = 0.0;
  ctx.config.absTolerance = 1e-18;
  Mesh& m = ctx.mesh;
  m.positions = {pinPos, freePos};
  m.restPositions = m.positions;
  m.velocities = {Vec3{}, Vec3{}};
  m.masses = {kM, kM};
  m.pinned = {1, 0};
  m.pinPositions = m.positions;
  m.edges.push_back(Edge{0, 1, kEll, kK});
  m.buildSparsityPattern();
  ensureBuffers(ctx);
  return ctx;
}

}  // namespace

int main() {
  int fails = 0;
  const Vec3 noGravity{0, 0, 0};

  // ---------- 判据 1：无重力，已处于静止长度 ⇒ 必须严格不动 ----------
  std::printf("=== 1) 无重力 + 初始就在静止长度：必须保持不动（位移严格为 0）===\n\n");
  std::printf("  %-28s %-16s %s\n", "pin 位置", "自由点最大位移", "判定");
  const Vec3 pins[] = {{0, 0, 0}, {1, 0, 0}, {0, 5, 0}, {3, -2, 7}};
  for (const Vec3& p : pins) {
    SimContext c = spring(p, p + Vec3{0, 0, kEll}, noGravity);
    Scalar worst = 0.0;
    for (int s = 0; s < 50; ++s) {
      stepOnce(c);
      worst = std::max(worst, length(c.mesh.positions[1] - (p + Vec3{0, 0, kEll})));
    }
    char buf[128];
    std::snprintf(buf, sizeof(buf), "(%.1f,%.1f,%.1f)", p.x, p.y, p.z);
    const bool ok = worst < 1e-12;
    if (!ok) ++fails;
    std::printf("  %-28s %-16.3e %s\n", buf, worst, ok ? "OK" : "NG");
  }

  // ---------- 判据 2：平移一致性（有重力）----------
  std::printf("\n=== 2) 有重力：pin 在原点 vs 平移后的相对位移必须一致 ===\n\n");
  const Vec3 g{0, -9.81, 0};
  SimContext a = spring(Vec3{0, 0, 0}, Vec3{0, 0.12, 0}, g);
  SimContext b = spring(Vec3{10, 20, 30}, Vec3{10, 20.12, 30}, g);
  for (int s = 0; s < 200; ++s) {
    stepOnce(a);
    stepOnce(b);
  }
  const Vec3 relA = a.mesh.positions[1] - a.mesh.positions[0];
  const Vec3 relB = b.mesh.positions[1] - b.mesh.positions[0];
  const Scalar diff = length(relA - relB);
  std::printf("  pin 在原点 : 相对位置 = (%.6f, %.6f, %.6f)\n", relA.x, relA.y, relA.z);
  std::printf("  pin 平移后 : 相对位置 = (%.6f, %.6f, %.6f)\n", relB.x, relB.y, relB.z);
  std::printf("  两者之差   = %.3e  %s\n", diff, diff < 1e-9 ? "OK（平移不变）" : "NG（平移后行为不同）");
  if (diff >= 1e-9) ++fails;

  // ---------- 判据 3：布料整体平移后的相对形状 ----------
  std::printf("\n=== 3) 布料：平移后相对形状（顶点间距离）必须一致 ===\n\n");
  auto cloth = [&](const Vec3& shift) {
    SceneConfig cfg;
    cfg.gridNx = 8;
    cfg.gridNy = 8;
    cfg.gridSpacing = 0.05;
    cfg.dt = 1.0 / 120.0;
    cfg.stiffness = kK;
    cfg.maxIterations = 20;
    cfg.relTolerance = 1e-8;
    cfg.gravity = g;
    SimContext c = makeScene(cfg);
    for (auto& p : c.mesh.positions) p += shift;
    for (auto& p : c.mesh.pinPositions) p += shift;
    for (auto& p : c.mesh.restPositions) p += shift;
    refreshPinPositions(c);
    for (int s = 0; s < 300; ++s) stepOnce(c);
    return c;
  };
  SimContext c0 = cloth(Vec3{0, 0, 0});
  SimContext c1 = cloth(Vec3{5, 7, -3});

  // 比较所有边长的相对差
  Scalar worstRel = 0.0;
  for (int e = 0; e < c0.mesh.edgeCount(); ++e) {
    const auto& ed = c0.mesh.edges[static_cast<std::size_t>(e)];
    const Scalar l0 = length(c0.mesh.positions[static_cast<std::size_t>(ed.a)] -
                             c0.mesh.positions[static_cast<std::size_t>(ed.b)]);
    const Scalar l1 = length(c1.mesh.positions[static_cast<std::size_t>(ed.a)] -
                             c1.mesh.positions[static_cast<std::size_t>(ed.b)]);
    worstRel = std::max(worstRel, std::fabs(l0 - l1) / std::max(1e-12, l0));
  }
  std::printf("  平移前/后各边长的最大相对差 = %.6e  %s\n", worstRel,
              worstRel < 1e-6 ? "OK（平移不变）" : "NG（平移后形状不同）");
  if (worstRel >= 1e-6) ++fails;

  std::printf("\n  %s（失败 %d 项）\n", fails == 0 ? "全部通过" : "存在失败", fails);
  return fails == 0 ? 0 : 1;
}
