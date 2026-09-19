// _verify/pin_positions.cpp
// 核对被 pin 的顶点到底在什么位置：应当位于顶边整行（z 最大那一行，y=0 平面）
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

#include "core/sim/Integrator.h"
#include "core/sim/Scene.h"

using namespace pd;

int main() {
  const int n = 6;
  SceneConfig cfg;
  cfg.gridNx = n;
  cfg.gridNy = n;
  cfg.gridSpacing = 0.05;
  cfg.stiffness = 1e4;
  SimContext ctx = makeScene(cfg);
  const Mesh& m = ctx.mesh;

  std::printf("=== 网格 %dx%d，间距 %.3f（顶点按 j*nx+i 排列，y 轴朝上）===\n\n", n, n,
              cfg.gridSpacing);
  std::printf("坐标范围: x ∈ [0, %.3f], z ∈ [0, %.3f]\n\n", (n - 1) * cfg.gridSpacing,
              (n - 1) * cfg.gridSpacing);

  std::printf("被 pin 的顶点（Scene.cpp 默认规则）：\n");
  std::printf("  %-6s %-10s %-10s %-10s\n", "索引", "x", "y", "z");
  for (int v = 0; v < m.vertexCount(); ++v) {
    if (!m.isPinned(v)) continue;
    const Vec3& p = m.positions[static_cast<std::size_t>(v)];
    std::printf("  %-6d %-10.4f %-10.4f %-10.4f\n", v, p.x, p.y, p.z);
  }

  std::printf("\n代码里的默认 pin 规则（Scene.cpp）：\n");
  std::printf("    topLeft  = (ny-1)*nx + 0\n");
  std::printf("    topRight = (ny-1)*nx + (nx-1)\n");
  std::printf("  => 期望 pin 在 z = %.3f（最大 z，即 -z 侧观察的远边），x = 0 与 %.3f\n\n",
              (n - 1) * cfg.gridSpacing, (n - 1) * cfg.gridSpacing);

  // 统计每个 z 值上的 pin 数量，看 pin 落在哪一行
  std::printf("按 z 统计 pin 分布（应全部集中在最大 z）：\n");
  for (int j = 0; j < n; ++j) {
    int cnt = 0;
    for (int i = 0; i < n; ++i) {
      const int v = j * n + i;
      if (m.isPinned(v)) ++cnt;
    }
    if (cnt > 0) std::printf("  z=%.3f (第 %d 行): %d 个 pin\n", j * cfg.gridSpacing, j, cnt);
  }

  std::printf("\n按 x 统计（若 pin 落在某一列，说明索引方向被换过）：\n");
  for (int i = 0; i < n; ++i) {
    int cnt = 0;
    for (int j = 0; j < n; ++j) {
      if (m.isPinned(j * n + i)) ++cnt;
    }
    if (cnt > 0) std::printf("  x=%.3f (第 %d 列): %d 个 pin\n", i * cfg.gridSpacing, i, cnt);
  }
  // ---- 跑一段，找出应变最大的几条边 ----
  std::printf("\n=== 跑 300 步后，应变最大的 5 条边 ===\n");
  for (int s = 0; s < 300; ++s) stepOnce(ctx);
  struct Info { int a, b; double len, rest, strain; bool pinnedA, pinnedB; };
  std::vector<Info> infos;
  for (const auto& e : m.edges) {
    const double len = length(m.positions[static_cast<std::size_t>(e.a)] -
                              m.positions[static_cast<std::size_t>(e.b)]);
    infos.push_back({e.a, e.b, len, e.restLength,
                     std::fabs(len - e.restLength) / e.restLength, m.isPinned(e.a) != 0,
                     m.isPinned(e.b) != 0});
  }
  std::sort(infos.begin(), infos.end(),
            [](const Info& x, const Info& y) { return x.strain > y.strain; });
  std::printf("  %-8s %-8s %-12s %-12s %-12s %s\n", "边(a,b)", "当前长", "静止长", "比值", "应变",
              "含pin");
  for (int i = 0; i < 5 && i < (int)infos.size(); ++i) {
    const auto& r = infos[static_cast<std::size_t>(i)];
    std::printf("  (%d,%d)   %-12.5f %-12.5f %-12.5f %-12.5f %s\n", r.a, r.b, r.len, r.rest,
                r.len / r.rest, r.strain,
                (r.pinnedA || r.pinnedB) ? "是" : "否");
  }
  return 0;
}
