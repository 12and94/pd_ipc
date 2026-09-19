// core/sim/Scene.cpp
#include "core/sim/Scene.h"

#include <cmath>

#include "core/assemble/Assembler.h"
#include "core/solver/EigenDirectSolver.h"

namespace pd {

SimContext makeScene(const SceneConfig& config) {
  SimContext ctx;
  ctx.config = config;

  // ---- 网格 ----
  if (!config.meshPath.empty()) {
    if (!Mesh::loadObj(config.meshPath, ctx.mesh, config.stiffness)) {
      // 载入失败时退化为规则网格，保证 bench/viewer 仍可运行。
      ctx.mesh = Mesh::makeGrid(config.gridNx, config.gridNy, config.gridSpacing);
    }
  } else {
    ctx.mesh = Mesh::makeGrid(config.gridNx, config.gridNy, config.gridSpacing);
  }

  Mesh& m = ctx.mesh;
  // 刚度：直接使用配置值，**不做任何按间距的缩放**。
  //
  // 曾经在这里加过 kScale = (spacing0/spacing)²，那是错的：
  //   · 顶点质量已经随间距减小而减小（m ∝ spacing²），刚度无需再缩；
  //   · 那个系数把 κ=1e4 实际变成了 1600（spacing 0.05 时），让所有标定失真。
  // 若将来需要"换分辨率而材料手感不变"，正确做法是把刚度按**质量**参数化
  // （κ = k_mat · m），而不是按长度参数化。
  for (auto& e : m.edges) e.stiffness = config.stiffness;
  const Scalar density = (config.density > 0.0) ? config.density : 1.0;
  if (!m.triangleIndices.empty()) {
    m.computeMassesFromSurfaces(density);
  } else {
    // 规则网格：用静止长度作为面积份额的代理（与网格生成时的口径一致）。
    m.computeMassesFromEdges(density);
  }

  // ---- 默认 pin：规则网格把顶边（j = ny-1）两端的顶点钉住 ----
  // 对 OBJ 载入的网格，调用方可通过 config 之后的接口自行设置 pin。
  if (config.meshPath.empty() && config.gridNx >= 2 && config.gridNy >= 2) {
    const int nx = config.gridNx;
    const int ny = config.gridNy;
    const int topLeft = (ny - 1) * nx + 0;
    const int topRight = (ny - 1) * nx + (nx - 1);
    m.pinned[static_cast<std::size_t>(topLeft)] = 1;
    m.pinned[static_cast<std::size_t>(topRight)] = 1;
  }
  refreshPinPositions(ctx);

  // ---- 求解器 ----
  // 符号分解会推迟到第一次真正组装出 L 之后（见 stepOnce），
  // 这样符号结构必然与真实矩阵一致。
  ensureBuffers(ctx);

  return ctx;
}

void ensureBuffers(SimContext& ctx) {
  const std::size_t n = static_cast<std::size_t>(ctx.mesh.vertexCount());
  const std::size_t dim = 3 * n;
  if (ctx.predicted.size() != n) ctx.predicted.assign(n, Vec3{});
  if (ctx.positionsBeforeStep.size() != n) ctx.positionsBeforeStep.assign(n, Vec3{});
  if (ctx.positionsBeforePrevStep.size() != n) ctx.positionsBeforePrevStep.assign(n, Vec3{});
  if (ctx.targets.size() != static_cast<std::size_t>(ctx.mesh.edgeCount())) {
    ctx.targets.assign(static_cast<std::size_t>(ctx.mesh.edgeCount()), Vec3{});
  }
  if (static_cast<std::size_t>(ctx.b.size()) != dim) ctx.b.resize(dim);
  if (static_cast<std::size_t>(ctx.bBase.size()) != dim) ctx.bBase.resize(dim);
  if (static_cast<std::size_t>(ctx.xSolution.size()) != dim) ctx.xSolution.resize(dim);
  if (!ctx.solver) ctx.solver = std::make_unique<EigenDirectSolver>();
}

void refreshPinPositions(SimContext& ctx) {
  Mesh& m = ctx.mesh;
  for (int v = 0; v < m.vertexCount(); ++v) {
    if (m.isPinned(v)) m.pinPositions[static_cast<std::size_t>(v)] = m.positions[static_cast<std::size_t>(v)];
  }
}

}  // namespace pd
