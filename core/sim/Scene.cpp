// core/sim/Scene.cpp
#include "core/sim/Scene.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>

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

  // ---- 剪切（对角）约束：可选，默认关闭（`SceneConfig::shearStiffness == 0`）----
  //
  // 对角边就是一条普通的距离约束（`Edge`），所以装配、边着色、投影/散射、残差、
  // 以及 Phase 4c 的"三分量并行"都不需要任何改动；两条架构性质也不受影响：
  // ① 左端矩阵仍与位形无关 ⇒ 只分解一次；② 每块仍是标量 × I₃ ⇒ L = Ã ⊗ I₃ 保住。
  //
  // 位置刻意放在**质量计算之后**：这样"带剪切"与"不带剪切"的场景质量完全一致，
  // 对照实验里唯一的变量就是约束集（否则质量一变，动力学跟着变，数字不可比）。
  // 也支持用环境变量 `PD_SHEAR=<刚度>` 临时开启 —— 方便 bench/viewer 做对照，
  // 不必扩命令行解析（同 PD_TRACE_STEP 那类诊断开关的用法）。
  {
    Scalar shear = config.shearStiffness;
    if (const char* env = std::getenv("PD_SHEAR")) {
      const double v = std::atof(env);
      if (v > 0.0) shear = static_cast<Scalar>(v);
    }
    if (shear > 0.0 && config.meshPath.empty()) {
      const int added = ctx.mesh.addShearDiagonals(config.gridNx, config.gridNy, shear);
      std::printf("[scene] 剪切约束：新增对角边 %d 条（刚度 %g）⇒ 总边数 %d\n", added,
                  static_cast<double>(shear), ctx.mesh.edgeCount());
    }
  }

  // ---- 线性（中点）弯曲约束：可选，默认关闭（`SceneConfig::bendStiffness == 0`）----
  //
  // 位置刻意放在**质量计算之后、pin 之前**：
  //   · 在质量之后 —— 弯曲 stencil 不参与质量计算（与剪切同理），所以"带弯曲"与
  //     "不带弯曲"的场景质量逐位相同，对照实验里唯一的变量就是约束集；
  //   · 在 pin 之前 —— stencil 的生成只沿行/列取连续三个顶点，与 pin 标记无关，
  //     放在 pin 前后结果等价；这样排是为了让"网格构造"（质量 → 剪切 → 弯曲）成一段、
  //     pin 只负责标记。
  //
  // 为什么弯曲**没有**对应的投影/散射代码：它的流形是 {x : A_c x = 0}（线性子空间），
  // 投影 p_c ≡ 0 ⇒ 局部步无事可做；右端只多一个与位形无关的常向量（装配期算一次）。
  // 所以这里只需要"生成 stencil"这一步，其余全在装配（Assembler）与残差（Integrator）里。
  // 环境变量 `PD_BEND=<刚度>` 可临时开启（照抄 `PD_SHEAR` 的位置与写法，便于 bench/viewer
  // 做对照而不用扩命令行解析）。
  {
    Scalar bend = config.bendStiffness;
    if (const char* env = std::getenv("PD_BEND")) {
      const double v = std::atof(env);
      if (v > 0.0) bend = static_cast<Scalar>(v);
    }
    if (bend > 0.0 && config.meshPath.empty()) {
      const int added = ctx.mesh.addBendingStencils(config.gridNx, config.gridNy, bend);
      std::printf("[scene] 弯曲约束：新增 stencil %d 条（刚度 %g）\n", added,
                  static_cast<double>(bend));
    }
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
  // 弯曲常向量：**只在真的用了弯曲约束时才分配**。空 vector 的 data() 是 nullptr，
  // 热路径上 `extraRhs ? ... : 0` 直接走"没有额外项"的分支 ⇒ 默认关闭时逐字不变。
  // 注意别写成"每子步都 assign/shrink"：本函数每次 stepOnce 都调用，那会把 O(dim) 的
  // 无谓工作放进热路径（`ensureBuffers` 的既有约定是"幂等且只在尺寸不对时才动"）。
  {
    const std::size_t want = (ctx.mesh.bendCount() > 0) ? dim : std::size_t{0};
    if (ctx.bendRhs.size() != want) ctx.bendRhs.assign(want, Scalar{0});
  }
  if (!ctx.solver) ctx.solver = std::make_unique<EigenDirectSolver>();
}

void refreshPinPositions(SimContext& ctx) {
  Mesh& m = ctx.mesh;
  for (int v = 0; v < m.vertexCount(); ++v) {
    if (m.isPinned(v)) m.pinPositions[static_cast<std::size_t>(v)] = m.positions[static_cast<std::size_t>(v)];
  }
}

}  // namespace pd
