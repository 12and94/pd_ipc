// core/sim/Integrator.cpp
#include "core/sim/Integrator.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>

#include "core/assemble/Assembler.h"
#include "core/energy/DistanceTerm.h"
#include "core/math/Parallel.h"

namespace pd {

namespace {

using Clock = std::chrono::steady_clock;

double secondsSince(const Clock::time_point& t0) {
  return std::chrono::duration<double>(Clock::now() - t0).count();
}

/// 诊断开关：设置环境变量 PD_TRACE_STEP=1 时，逐步打印各阶段的中间量。
/// 用途：当链路出现异常时，不必重新编译就能看到"预测/投影/散射/求解"的数值。
bool traceEnabled() {
  static const bool on = (std::getenv("PD_TRACE_STEP") != nullptr);
  return on;
}

/// 诊断开关：设置 PD_DEBUG_RESIDUAL=1 时，即使残差判据被关闭也计算并打印残差。
bool debugResidual() {
  static const bool on = (std::getenv("PD_DEBUG_RESIDUAL") != nullptr);
  return on;
}

/// 规模阈值：网格小于它时不上并行。
/// 依据：区域本身的开关 + 每迭代若干 barrier 在 18 线程下要几十微秒，
/// 而小网格的每阶段有效工作只有几微秒（实测见 docs/perf.md）。
/// 这是**经验阈值**，Phase 2（图着色消除归约流量）之后要重新标定 ——
/// 那时并行的收益区间会往下移动。用 `if` 子句串行化区域而不是复制代码路径，
/// 是为了保证全项目只有一份流程实现。
constexpr int kMinParallelWork = 512;

#define PD_TRACE(...)                             \
  do {                                            \
    if (traceEnabled()) {                         \
      std::fprintf(stderr, "[step] " __VA_ARGS__); \
      std::fprintf(stderr, "\n");                  \
      std::fflush(stderr);                         \
    }                                             \
  } while (0)

/// 本子步的**非线性残差**（单位 m/s²）：最大顶点上的不平衡力除以该点质量。
///
/// 推导：本子步的不动点条件是（自由顶点，见 docs/plan.md §2.2.1）
///     (m_i/h²)(x_i - x̂_i) - m_i·g - f_i^int = 0
/// 其中弹性力 f_i^int 等于"散射进右端的那一项"的相反数：
///     f_i^int = -Σ_c κ_c d_c(x)   （a 端 +κd_c、b 端 -κd_c）
/// 于是每个顶点的不平衡加速度为
///     r_i = |(m_i/h²)(x_i - x̂_i) + Σ_c κ_c d_c| / m_i
/// 注意这里**不含重力项**：重力已经通过 x̂ = x + hv + h²g 进入，且只进入一次
/// （在右端再补 -Mg 就是把重力算两遍，本项目踩过这个坑）。x 到达不动点时
/// 上式整体为 0，因此它就是"离不动点还有多远"的度量。
///
/// ---- 为什么**必须**把弯曲力算进来（改了残差口径的原因，务必保留）----
/// 残差是**唯一放行判据**（`cfg.residualTolerance > 0` 时位移判据完全不参与判定）。
/// 所以"残差里少算一项力"不是"报告得好看一点"，而是**直接改变了哪些子步被放行**：
///   · 少了弯曲项 ⇒ 一个明显还受弯的子步会被判"已收敛"⇒ **假阴性**；
///   · 仓库的既有立场是"假阴性比假阳性更危险"（见本函数下面那段 targets 过期的说明、
///     以及 docs/perf.md §7.1 的两口径自查），因为假阳性只是浪费迭代，假阴性是把
///     错的位形当成答案交出去，而且**不会有任何症状**。
/// 弯曲力的形式：顶点 v 的力 = Σ_{stencil 含 v} (-k·w_v·(Σ_t w_t x_t))，w = {1,-2,1}。
/// 它**线性且与位形无关地可算**，所以残差里就地重算即可（不需要任何中间量）。
///
/// 为什么不能用"相邻迭代位移"代替它：位移小只说明这一步没怎么动。
/// 高刚度下更新量是被 L 预条件后的梯度，小更新不保证小物理残差
/// （κ=1e5 时 40 次迭代只完成切向运动的 0.46%，而位移看起来"已经很小"）。

/// 顶点 v 的**弯曲力**（A_c 形式，全项目唯一口径）。
///
/// 串行残差（`nonlinearResidual`）与 gather 残差（`vertexResidual`）**共用本函数**：
/// 所以"两条残差口径逐位相同"在这条路径上是**构造性**的，不是靠在两处各写一份、
/// 再靠两张 CSR 的填表顺序一致去论证。这是本次改动刻意选的写法 ——
/// 残差里加项本身就是最容易出"两条路慢慢分叉"的地方。
///
/// 权重 w_v 由关联表带着走（`BendAdjacency::vertexStencilWeight`）：v 是 stencil 的
/// a/b/c 时 w_v 分别是 1/-2/1。三者的符号不同，所以不能只记 stencil 号 ——
/// "按顶点号的大小猜自己是哪一端"在一般情况下就是错的（stencil 的 (a,c) 顺序由生成器
/// 决定，不代表几何位置），那种简化会静默地把中间顶点的 -2 变成 +1。
inline Vec3 bendForceAt(const Mesh& mesh, const BendAdjacency& bendAdj, const Scalar* xSol, int v) {
  Vec3 force{0, 0, 0};
  if (bendAdj.vertexStencils.empty()) return force;  // 没有弯曲约束：零开销、默认关闭时逐字不变
  const std::size_t i = static_cast<std::size_t>(v);
  for (uint32_t k = bendAdj.vertexStart[i]; k < bendAdj.vertexStart[i + 1]; ++k) {
    const BendStencil& s = mesh.bends[bendAdj.vertexStencils[k]];
    const Scalar* xa = xSol + static_cast<std::size_t>(s.a) * 3;
    const Scalar* xb = xSol + static_cast<std::size_t>(s.b) * 3;
    const Scalar* xc = xSol + static_cast<std::size_t>(s.c) * 3;
    // A_s x = x_a - 2 x_b + x_c，逐分量写开（不写 Vector 表达式，避免 vtable 之外的
    // 运算符重载把 `-2*x_b` 算成别的结合顺序；与 surrogateEnergy 里那一处同形）。
    const Vec3 second{xa[0] - 2.0 * xb[0] + xc[0], xa[1] - 2.0 * xb[1] + xc[1],
                      xa[2] - 2.0 * xb[2] + xc[2]};
    // 力 = -∇E = -k·w_v·(A_s x)
    force += second * (-s.stiffness * bendAdj.vertexStencilWeight[k]);
  }
  return force;
}

}  // namespace

Scalar nonlinearResidual(const SimContext& ctx) {
  const Mesh& m = ctx.mesh;
  const int n = m.vertexCount();
  const Scalar invH2 = 1.0 / (ctx.config.dt * ctx.config.dt);
  const BendAdjacency& bendAdj = m.bendAdjacency();

  // 约束力的正确形式：f^int = κ(AᵀA x − Aᵀ d)，其中 d 必须是**当前位置**的投影。
  //
  // **这里踩过一个坑（务必保留说明）**：最初版本直接用了 `ctx.targets`，但那是
  // `stepOnce` 迭代循环**最后一次迭代开始前**算的投影。一旦循环提前退出
  // （位移/残差判据命中），targets 相对当前位置就**过期**了：此时
  //     |κ(rel − d_过期)| ≈ κ·|d_新 − d_旧|  ≫  κ(r − ℓ)
  // 于是残差被算得远小于真实值 —— 实测给出 1.2e-6 m/s² 的假阴性，
  // 而用当前位置重算投影后是 10.28 m/s²。**假阴性比假阳性更危险**。
  // 因此这里就地重算 d_c，不依赖 ctx.targets。
  static std::vector<Scalar> force;
  force.assign(static_cast<std::size_t>(n) * 3, Scalar{0});
  for (std::size_t c = 0; c < m.edges.size(); ++c) {
    const Edge& e = m.edges[c];
    const Vec3 rel = m.positions[static_cast<std::size_t>(e.a)] -
                     m.positions[static_cast<std::size_t>(e.b)];  // A_c x
    const Scalar len = length(rel);
    // d_c = ℓ·unit(x_a − x_b)（当前位置的投影；退化时按 0 处理，与 project() 一致）
    const Vec3 d = (len <= DistanceTerm::kMinLength) ? Vec3{} : rel * (e.restLength / len);
    const Vec3 contrib = rel * e.stiffness - d * e.stiffness;  // κ(A_c x − d_c)
    const std::size_t ia = static_cast<std::size_t>(e.a) * 3;
    const std::size_t ib = static_cast<std::size_t>(e.b) * 3;
    force[ia + 0] += contrib.x;
    force[ia + 1] += contrib.y;
    force[ia + 2] += contrib.z;
    force[ib + 0] -= contrib.x;
    force[ib + 1] -= contrib.y;
    force[ib + 2] -= contrib.z;
  }
  // 弯曲约束的力：**必须算进来**（见本函数顶部"为什么必须把弯曲力算进来"）。
  //
  // 顺手把这条路径与 gather 版做成**逐位相同**：就地按 `bendForceAt(...)` 逐顶点算，
  // 而那个 helper 也正是 gather 版用的**同一份实现** ⇒ 两边不只是"同一个公式"，
  // 而是**同一段代码**（不留"两处各写一份、将来改一处忘另一处"的口子）。
  // 这比"串行版按 stencil 序全局累加"更保守：后者的累加顺序与 gather 版要靠
  // `buildBendAdjacency` 按 stencil 号升序填表来**论证**，而复用同一段代码让这一点
  // 变成构造性的（不需要论证，也不可能不一致）。
  //
  // 代价：每条 stencil 的 A_s x 被它的三个顶点各算一次（共享顶点要重算）。
  // 量级 O(3B) 次乘加，40×40 上 B = 3040，与 3120 条边的距离约束残差同阶 —— 不构成热点。
  // 没有弯曲约束时 vertexStencils 为空，helper 立刻返回零 ⇒ 零开销、行为逐字不变。
  if (!bendAdj.vertexStencils.empty()) {
    // 位置一律取 xSolution（与 gather 口径一致）：它是刚做完回代的结果，
    // 而 m.positions 要等"解包趟"才写回。两者数值相同（解包写的就是 xSolution），
    // 但统一从一处取省掉了任何读序上的疑问。
    const Scalar* xSol = ctx.xSolution.data();
    for (int v = 0; v < n; ++v) {
      const Vec3 fb = bendForceAt(m, bendAdj, xSol, v);
      force[static_cast<std::size_t>(v) * 3 + 0] += fb.x;
      force[static_cast<std::size_t>(v) * 3 + 1] += fb.y;
      force[static_cast<std::size_t>(v) * 3 + 2] += fb.z;
    }
  }

  Scalar worst = 0.0;
  int worstVertex = -1;
  for (int v = 0; v < n; ++v) {
    if (m.isPinned(v)) continue;  // pinned 顶点不是未知量，其行被整行覆盖
    const std::size_t i = static_cast<std::size_t>(v);
    const Scalar mass = m.masses[i];
    if (mass <= 0.0) continue;
    const Vec3 inertia = (m.positions[i] - ctx.predicted[i]) * (mass * invH2);
    const Vec3 total{inertia.x + force[i * 3 + 0], inertia.y + force[i * 3 + 1],
                     inertia.z + force[i * 3 + 2]};
    const Scalar r = length(total) / mass;
    if (r > worst) {
      worst = r;
      worstVertex = v;
    }
  }

  // 调试开关：PD_DEBUG_RESIDUAL=1 时打印最大残差顶点的逐项分解，
  // 并**同时**打印"静力不平衡" |f_int/m - g| —— 后者是另一种口径
  // （只看弹性力与重力是否平衡、不含惯性项），用于和外部诊断对照。
  if (std::getenv("PD_DEBUG_RESIDUAL") != nullptr && worstVertex >= 0) {
    const std::size_t i = static_cast<std::size_t>(worstVertex);
    const Scalar mass = m.masses[i];
    const Vec3 inertia = (m.positions[i] - ctx.predicted[i]) * (mass * invH2);
    const Vec3 f{force[i * 3 + 0], force[i * 3 + 1], force[i * 3 + 2]};
    const Vec3 staticImbalance = f / mass - ctx.config.gravity;
    std::fprintf(stderr,
                 "[residual] 顶点 %d  |x-x̂|=%.3e  |inertia|=%.3e  |f_int|=%.3e  "
                 "|inertia+f_int|/m=%.6e m/s^2  |f_int/m-g|=%.6e m/s^2\n",
                 worstVertex, length(m.positions[i] - ctx.predicted[i]), length(inertia), length(f),
                 worst, length(staticImbalance));
  }
  return worst;
}

// ---------------------------------------------------------------------------
// Phase 3：残差的**逐顶点 gather** 版本（上面 nonlinearResidual 是它的串行对照）
//
// 为什么这里可以用 gather（而散射那边 gather 只是平手）：
//   · 散射写的是 b 的 3 个槽位，gather 要求"每个顶点只写自己的槽位" —— 散射本来
//     就是 scatter（一个顶点被多条边写），gather 化必须把"逐边投影"改成"逐顶点重算
//     投影"，多算一遍投影，换到的只是少几个 barrier。实测平手（docs/perf.md §3.4）。
//   · 残差本来就是 **O(顶点数) 的结果**（每顶点一个 r_i，再取 max），gather 与它的
//     数据流同形：顶点读自己的关联边、只写自己的局部量，天然无冲突、无归约、
//     位级可复现（求和顺序由 CSR 决定，与线程数/分块无关）。
//
// 与 nonlinearResidual 的关系：**同一口径**（同一公式、同一投影、同一退化保护）。
// 求和顺序也一致 —— 关联表（CSR）是**按边号升序**填的（buildAdjacency 用 cursor 顺序遍历
// 边表），所以"顶点 v 的关联边序"与"串行版按全局边序累加进 force[v]"是同一个顺序，
// 两者**逐位相同**（实测 40×40、8 次迭代/子步，无一处偏差）。
// PD_DEBUG_RESIDUAL=1 时 stepOnce 会把两个值都算出来做自查（见那里的 1e-9 告警）。
//
// 端点位置一律从 xSolution 取，而**不是** m.positions：xSolution 是刚做完的回代结果，
// 而 m.positions 要等"解包趟"写回。用 xSolution 就不存在"邻居还没解包"的读序依赖，
// 于是残差能与解包融合在**同一趟**里 —— 这是本次改造省下的主要成本
// （原来残差是 single 里的串行趟，且要它自己的一趟全边遍历）。
// 两者数值上完全一致：解包趟写的就是 m.positions[i] = xSolution[i]（含 pinned 行，
// pin 行的解被 applyPinRhs 覆盖为 q，因此 pinned 端点用 xSolution 读到的也还是 q）。
inline Scalar vertexResidual(const Mesh& mesh, const ConstraintColoring& adj,
                             const BendAdjacency& bendAdj, const Scalar* xSol,
                             const std::vector<Vec3>& predicted, Scalar mass, Scalar invH2, int v) {
  const std::size_t i = static_cast<std::size_t>(v);
  Vec3 force{0, 0, 0};
  const std::size_t begin = adj.vertexStart[i];
  const std::size_t end = adj.vertexStart[i + 1];
  for (std::size_t k = begin; k < end; ++k) {
    const Edge& e = mesh.edges[adj.vertexEdges[k]];
    const Scalar* pa = xSol + static_cast<std::size_t>(e.a) * 3;
    const Scalar* pb = xSol + static_cast<std::size_t>(e.b) * 3;
    const Vec3 rel{pa[0] - pb[0], pa[1] - pb[1], pa[2] - pb[2]};  // A_c x
    const Scalar len = length(rel);
    // d_c 就地按当前位置重算（不能复用 ctx.targets：那是最后一次迭代**开始前**的投影，
    // 提前退出时已过期，会算出远小于真值的假阴性 —— 见 nonlinearResidual 顶部那段）。
    const Vec3 d = (len <= DistanceTerm::kMinLength) ? Vec3{} : rel * (e.restLength / len);
    // 与 nonlinearResidual 逐字同形：两个乘法再相减（不写成 (rel-d)*κ，
    // 那样会多引入一次舍入，让"两个口径"的差从纯求和顺序差异变成额外的乘法顺序差异）。
    const Vec3 contrib = rel * e.stiffness - d * e.stiffness;  // κ(A_c x − d_c)
    if (e.a == v) {
      force += contrib;
    } else {
      force -= contrib;
    }
  }
  // 弯曲力：与串行版**共用同一段代码**（`bendForceAt`），所以两者逐位相同是构造性的，
  // 而不是"靠两张表顺序一致来论证"。漏掉这一项会让残差偏小 ⇒ 假阴性（见顶部说明）。
  // 没有弯曲约束时 vertexStencils 为空，helper 立刻返回零 ⇒ 默认关闭时逐字不变。
  force += bendForceAt(mesh, bendAdj, xSol, v);

  const Vec3 xv{xSol[3 * i + 0], xSol[3 * i + 1], xSol[3 * i + 2]};
  const Vec3 inertia = (xv - predicted[i]) * (mass * invH2);
  return length(inertia + force) / mass;
}

// 关联表（CSR）契约自查：**一次/子步**，放在进并行区域之前。
//
// 为什么不在散射里查（DistanceTerm::checkColoringContract 的邻居位置显然更方便）：
// 实测把这段 O(1) 校验放进 `scatterIntoInRegion` 会让 1 线程的散射阶段
// 100–104 ms → 128–133 ms（**+28 %**，12,000 次调用/300 子步）—— 开销与"几条指令"
// 完全不相称，说明是编译器的代码布局被扰动（详见 DistanceTerm.cpp 末尾那段注记）。
// 而关联表的唯一使用者就是下面的残差 gather，所以由它自己保证前置条件最合适：
// 错了就在**进区域之前**报错，不会出现"跑到一半才发现越界"。
void checkAdjacencyContract(const Mesh& mesh) {
  const ConstraintColoring& coloring = mesh.constraintColoring();
  const std::size_t n = mesh.positions.size();
  if (coloring.vertexStart.size() != n + 1 ||
      coloring.vertexEdges.size() != static_cast<std::size_t>(mesh.edgeCount()) * 2) {
    std::fprintf(stderr,
                 "[residual] 违反前置条件：顶点关联表已过期 —— vertexStart %zu 项（期望 %zu）、"
                 "vertexEdges %zu 项（期望 %zu）。\n"
                 "  残差的 gather 用它当循环边界，规模不对就是越界读（实测表现为卡住不返回）。\n"
                 "  改了拓扑（增删边/改顶点数）之后必须重新调用 Mesh::buildSparsityPattern()。\n",
                 coloring.vertexStart.size(), n + 1, coloring.vertexEdges.size(),
                 static_cast<std::size_t>(mesh.edgeCount()) * 2);
    std::fflush(stderr);
    std::abort();
  }
  // **弯曲 stencil 的关联表也要查**（同一条理由：它是残差 gather 的循环边界）。
  // 这张表有一条额外且很容易踩的失效方式：先建好表再往 `mesh.bends` 里 push
  //（照 Model 里的构造顺序：makeGrid → addShearDiagonals → addBendingStencils 都会重建，
  //  但如果有人绕过这两个入口手工 push，表就停在旧规模）。
  // 症状与"漏建 vertexEdges"完全一样 —— 越界读、卡住不返回，所以宁可在这里早失败。
  const BendAdjacency& bendAdj = mesh.bendAdjacency();
  if (bendAdj.vertexStart.size() != n + 1 ||
      bendAdj.vertexStencils.size() != static_cast<std::size_t>(mesh.bendCount()) * 3 ||
      bendAdj.vertexStencilWeight.size() != bendAdj.vertexStencils.size()) {
    std::fprintf(stderr,
                 "[residual] 违反前置条件：弯曲 stencil 关联表已过期 —— vertexStart %zu 项（期望 %zu）、"
                 "vertexStencils %zu 项（期望 %zu = 3×%d 条 stencil）。\n"
                 "  残差用它算弯曲力，规模不对就是越界读。\n"
                 "  增删弯曲 stencil 之后必须重新调用 Mesh::buildSparsityPattern()"
                 "（正常路径由 Mesh::addBendingStencils 自动完成）。\n",
                 bendAdj.vertexStart.size(), n + 1, bendAdj.vertexStencils.size(),
                 static_cast<std::size_t>(mesh.bendCount()) * 3, mesh.bendCount());
    std::fflush(stderr);
    std::abort();
  }
}

void resetStageTimes(SimContext& ctx) { ctx.times = StageTimes{}; }

int stepOnce(SimContext& ctx) {
  const auto tStep0 = Clock::now();

  Mesh& m = ctx.mesh;
  const SceneConfig& cfg = ctx.config;
  const Scalar h = cfg.dt;
  const Scalar invH = 1.0 / h;
  // 与 nonlinearResidual 里**逐字同形**（1/(h·h) 而不是 invH·invH）：后者是另一种舍入，
  // 会让"串行/并行残差"的对照凭空多出一个 1 ulp 的差，掩盖真正的求和顺序差异。
  const Scalar invH2 = 1.0 / (h * h);
  const int n = m.vertexCount();
  const std::size_t dim = static_cast<std::size_t>(3 * n);

  // 自保护：任何构造路径（makeScene / 测试 / 调试程序）都可能忘记分配缓冲。
  // 这里统一保证尺寸正确，避免越界写（这是实际踩过的坑）。
  ensureBuffers(ctx);
  // 约束着色的懒构建兜底（幂等、O(1) 检查）：**必须在进并行区域之前**做，
  // 否则会在区域里"一边建一边被别的线程读"。正常路径由 buildSparsityPattern 建好。
  m.ensureConstraintColoring();
  // 弯曲 stencil 关联表的懒构建兜底（幂等、O(B)）：同一条理由，**必须在进区域之前**
  //（残差的 gather 用它当循环边界；区域里构建会一边建一边被别的线程读）。
  m.ensureBendAdjacency();
  checkAdjacencyContract(m);  // 残差 gather 的前置条件（一次/子步，不进每迭代热路径）

  // ---------------------------------------------------------------
  // 2) 判定是否需要重新组装 + 数值分解（区域外：它是串行工作，且必须先于迭代）。
  //    正常运行时（拓扑/刚度/h/pin 都没变）这里**一次都不会触发**，
  //    这正是"预分解复用"的收益所在（docs/plan.md §2.3 不变量）。
  // ---------------------------------------------------------------
  {
    const auto t0 = Clock::now();
    const SolverStamp now = computeStamp(m, h, cfg.velocityDamping, ctx.topologyId);
    if (!ctx.stampValid || now != ctx.stamp) {
      assembleLeftHandSide(m, h, cfg.velocityDamping, ctx.L, &ctx.bendRhs);
      // 符号分解的判定与数值分解**不同**：数值分解看数值戳（刚度/mass/h），
      // 符号分解还要看**稀疏结构**。pin 掩码变化会让"两端都自由"的耦合块消失，
      // 结构随之变化 —— 那时必须重新 analyze，否则 factorize 一个结构不同的矩阵后
      // solve 会访问未定义数据（见 IGlobalSolver.h，本项目因此踩过一次访问违例）。
      // 弯曲 stencil 的顶点三元组同理（它引入新的非零块）—— 已进 computeStructureStamp。
      const uint64_t structureNow = computeStructureStamp(m, ctx.topologyId);
      if (!ctx.solver->analyzed() || structureNow != ctx.structureStamp) {
        // 结构变了（或首次）：按**真实矩阵**的结构重新做符号分解，再数值分解。
        ctx.solver->analyze(3 * n, ctx.L);
        ctx.structureStamp = structureNow;
      }
      ctx.solver->factorize(ctx.L);
      ctx.stamp = now;
      ctx.stampValid = true;
      ctx.factorizeCount += 1;
      ctx.times.assemble += secondsSince(t0);
      PD_TRACE("重新组装+分解（第 %d 次）", ctx.factorizeCount);
    } else if (!m.bends.empty()) {
      // **弯曲的常向量必须每子步重算**（O(B)，可忽略），因为 `assembleLeftHandSide`
      // 只在数值 stamp 变化时才跑，而它依赖 `pinPositions` —— 拖拽把手**只改 b、不改 L**
      //（docs/plan.md §2.2.1 的约定，`stampChangesOnlyWhenLeftHandSideValuesChange` 钉住了），
      // 于是"pin 位置变了但 stamp 没变"是常态。沿用旧的常向量会让自由端少掉 pin 的
      // **当前位置**信息（相当于用一步之前的手把位置做消元补偿）。
      // 没有弯曲约束时这一段整个不进（连一次 size 判断都不做）⇒ 默认关闭时行为逐字不变。
      assembleBendingRhs(m, ctx.bendRhs);
    }
  }

  // ---------------------------------------------------------------
  // 3)–5) 一个子步 = **一个并行区域**（docs/parallel-refactor.md §4.1）
  //
  //   改造前：`#pragma omp parallel` 写在 project / scatterInto **函数内部**，
  //   于是每子步要 fork/join 2×K 次（K = 迭代数；40 迭代 → 80 次），
  //   而每次区域开关在 18 线程下要 10–25 µs，比局部步的有效工作还贵。
  //   现在：整个子步只进/出区域一次，各阶段用区域内的 `omp for` 分工，
  //   同步成本从"fork/join"降到"barrier"。
  //
  //   依赖边界（哪能省 barrier、哪必须同步，逐条都给理由）：
  //     · 预测 → 惯性右端：靠右端的 single 隐式 barrier 同步；
  //     · 惯性右端 → 散射(以 bBase 为基值)：同上，single 的 barrier 已保证；
  //     · 局部步 → 散射：散射内部的 single（清缓冲）带隐式 barrier；
  //     · 散射 → pin 覆盖 / 回代：散射的归约趟是 `omp for`（带隐式 barrier）；
  //     · 回代 → 解包+判据：pin 覆盖与回代合在同一个 `omp for` 里，其末尾 barrier 覆盖；
  //     · 迭代末尾：判据 single 的 barrier 让所有线程看到同一个 `ctx.converged`，
  //       因此 `break` 在所有线程上一致（OpenMP 要求如此）。
  //
  //   小网格不上并行：区域本身的开销会超过收益（实测见 docs/perf.md）。
  //   用 `if` 子句把区域"串行化"，而不是再复制一条串行代码路径 ——
  //   这样全项目只有一份流程实现（决策：不做运行时双实现开关）。
  // ---------------------------------------------------------------
  std::vector<Vec3> previous(static_cast<std::size_t>(n));
  const bool useParallel =
      (numThreads() > 1) && (n >= kMinParallelWork) && (m.edgeCount() >= kMinParallelWork);

  // 收敛扫描的"每线程部分最大值"槽位（共享，区域外分配）。
  // 为什么不用 `reduction(max:...)`：MSVC 的经典 OpenMP 2.0 只支持
  // +,*,−,&,^,|,&&,|| 这些运算符，max/min 要 /openmp:llvm 才认（编译器明确报 C7660）。
  // 本阶段**刻意不换运行时** —— 换运行时（vcomp → libomp）会同时改变线程池行为，
  // 那样就分不清"结构改造的收益"与"运行时的差异"了（换运行时留给 Phase 4 单独评估）。
  // max 是精确且可结合的运算，因此"每线程部分值 + 固定顺序合并"与串行扫描**逐位相同**。
  const int mergeSlots = std::max(1, numThreads());
  std::vector<Scalar> partialDiff(static_cast<std::size_t>(mergeSlots), 0.0);
  std::vector<Scalar> partialScale(static_cast<std::size_t>(mergeSlots), 0.0);
  // 残差也一样是"每线程部分最大值 + 固定顺序合并"（max 精确且可结合）。
  // 额外带一个"取到最大值的顶点号"槽位，用于 PD_DEBUG_RESIDUAL 的逐项分解；
  // 相同残差时取**更小**的顶点号（各线程槽位按线程号升序合并 ⇒ 与线程数无关）。
  std::vector<Scalar> partialResidual(static_cast<std::size_t>(mergeSlots), 0.0);
  std::vector<int> partialResidualVertex(static_cast<std::size_t>(mergeSlots), -1);

  // 残差要**并行算**的条件：判据启用了，或开了调试口径。
  // 调试口径（PD_DEBUG_RESIDUAL）需要 worst 顶点的 force[] 逐项分解，串行版顺手就有；
  // 为了不让并行路径为调试多背一份全量 force 缓冲，调试运行里报告串行值、并行值只作对照。
  const bool debugRes = debugResidual();
  const bool wantResidual = (cfg.residualTolerance > 0.0) || debugRes;
  // 两条残差路径的**选择**（不是两个定义 —— 两者数值逐位相同，见 vertexResidual 顶部）：
  //   · 区域内并行时：逐顶点 gather，与解包趟融合（一份工作量摊到 4 个线程上，最快）；
  //   · 区域被 `if` 子句串行化时（单线程 / 小网格）：边序串行版。它每条边只算一次投影，
  //     而 gather 版要让两个端点各算一次 —— 串行时没人替它摊掉这笔多出来的 √与除法，
  //     实测 40×40 会白掉 6.5 %（6.010 → 6.403 ms/子步）。既然两者结果逐位相同，
  //     就没有理由让串行路径去付这笔钱。
  //   · 调试时**两条都算**并互相校对（1e-9 告警），所以此时不受上面这条选择影响。
  const bool gatherResidual = wantResidual && (useParallel || debugRes);
  const bool serialResidual = wantResidual && (!useParallel || debugRes);

  ctx.converged = false;  // 本子步是否达到收敛判据（而非用尽迭代预算）
  ctx.iterationsUsed = 0;

#ifdef _OPENMP
#pragma omp parallel num_threads(numThreads()) if (useParallel)
#endif
  {
    // 诊断开关（默认关闭，不改变任何行为）：PD_AFFINITY=workers|exclusive 时，把"参与全局步
    // 的那几条线程"固定到各自独立的物理核上；workers 模式下其余线程保持**完全自由调度**。
    // 用途仅限标定落点/绑核的影响，结论见 docs/open-issues.md §3。见 core/math/Parallel.h。
    pinSolveThreadsOnce(ctx.solver->parallelComponents());

    // ---------------------------------------------------------------
    // 3) 预测：x̂ = x + h v + h² g
    //    与"保存上一步位置""初始化 previous"融合成一趟（逐顶点，都是只写自己）。
    //    重力只在这里出现一次。
    // ---------------------------------------------------------------
    {
      const auto t0 = Clock::now();
      const Vec3 gh2 = cfg.gravity * (h * h);
      // 不带 nowait：紧接着的 single（惯性右端）需要预测已全部完成。
#ifdef _OPENMP
#pragma omp for
#endif
      for (long long v = 0; v < n; ++v) {
        const std::size_t i = static_cast<std::size_t>(v);
        ctx.positionsBeforePrevStep[i] = ctx.positionsBeforeStep[i];
        const Vec3 p = m.positions[i];
        ctx.positionsBeforeStep[i] = p;
        ctx.predicted[i] = p + m.velocities[i] * h + gh2;
        previous[i] = p;  // 收敛判据的参照：本子步开始时的位置
      }
      if (threadId() == 0) ctx.times.predict += secondsSince(t0);
    }

    // ---------------------------------------------------------------
    // 4a) 右端的惯性部分：bBase = (M/h²) x̂（串行；只有这一项，重力不在这里）
    //     single 的隐式 barrier 同时充当"预测趟已完成"的同步点。
    // ---------------------------------------------------------------
#ifdef _OPENMP
#pragma omp single
#endif
    { assembleInertialRhs(m, ctx.predicted, h, cfg.velocityDamping, ctx.bBase); }

    // ---------------------------------------------------------------
    // 4) PD 迭代：局部步 → 散射 → 覆盖 pin → 全局步 → 判据（循环留在区域内）
    // ---------------------------------------------------------------
    // ---------------------------------------------------------------
    // Chebyshev 加速（Wang 2015 的形式）。**默认关**：不设 PD_CHEB 时，下面的
    // `chebEnabled == false` ⇒ 热路径与改动前**逐字不变**（门禁第 ④ 步因此不受影响）。
    //
    //   PD_CHEB=1          启用；PD_CHEB_RHO=ρ（默认 0.98，"宁小勿大"：估大了会震荡）
    //   PD_CHEB_START=S    前 S 次仍跑朴素 PD（默认 10，论文取法），之后转 Chebyshev
    //
    // 形式：x_k ← ω_k (f(x_k) − x_{k−1}) + x_{k−1}，ω 递推
    //   k ≤ S: ω = 1；k = S+1: ω = 2/(2−ρ²)；k > S+1: ω = 4/(4−ρ²·ω_prev)
    // 注意 ω = 1 时**不做合成**（(a−b)+b 在浮点下 ≠ a）⇒ 预热阶段仍是"逐字朴素 PD"。
    //
    // 口径说明（v1 已知偏差）：非线性残差仍由**未加速解**算出（xSolution 是这一趟的
    // 邻居读源，若就地改写会破坏"邻居读序无关"的前提、引入不确定性）。对慢模态而言
    // 未加速残差是**偏悲观**的估计 ⇒ 残差判据不会过早放行，这是安全的方向。
    // 速度更新用的是 m.positions，也就是**加速后**的位形。
    static const bool chebEnabled = [] {
      const char* s = std::getenv("PD_CHEB");
      return s != nullptr && *s != '\0' && *s != '0';
    }();
    static const Scalar chebRho = [] {
      const char* s = std::getenv("PD_CHEB_RHO");
      return (s != nullptr && *s != '\0') ? static_cast<Scalar>(std::atof(s)) : Scalar{0.98};
    }();
    static const int chebStart = [] {
      const char* s = std::getenv("PD_CHEB_START");
      return (s != nullptr && *s != '\0') ? std::atoi(s) : 10;
    }();
    Scalar chebOmega = Scalar{1};
    for (int k = 0; k < cfg.maxIterations; ++k) {
      // 4b+4c) 融合的"投影 + 散射"（2026-09-22）：按颜色就地算 d_c 并直接累加到右端，
      //    不再物化 targets 中间量（省掉它的写+读，以及两个端点位置的一次重复读取）。
      //    语义与原来"projectInRegion 之后再 scatterIntoInRegion"等价，且两条路共用同一份
      //    算术 helper ⇒ 结果**逐位相同**（见 DistanceTerm.h）。桶归属：并入"散射"桶
      //    （"局部步"不再是独立阶段，这一趟 = 原来的局部步 + 散射）。
      {
        const auto t0 = Clock::now();
        // `bendRhs` 为空（未启用弯曲）时 data() 是 nullptr ⇒ 种子趟走"没有额外项"的分支，
        // 与改动前的 `b[i] = base[i]` 逐字等价（这就是"默认关闭逐字不变"在热路径上的落点）。
        DistanceTerm::projectAndScatterIntoInRegion(m, ctx.bBase.data(), ctx.b.data(), dim,
                                                    ctx.bendRhs.empty() ? nullptr
                                                                        : ctx.bendRhs.data());
        if (threadId() == 0) ctx.times.scatter += secondsSince(t0);
      }

      // 4d) 覆盖 pinned 行 + 4e) 全局步
      //
      // **2026-09-22：这里从 `omp single` 改成了 `omp for`。**
      // 回代占单子步 76 %，一直是 Amdahl 天花板；但本项目**不需要**并行化三角求解本身 ——
      // L 的每一块都是标量 × I₃（L = Ã ⊗ I₃，见 Assembler.cpp 与 DistanceTerm::assembleMatrix），
      // 所以 3n 个方程其实是 3 个**互不相连**的标量系统，可以按 x/y/z 拆成 3 条独立链：
      // 零同步、零竞争、而且每个分量的算术序列与整趟完全一致 ⇒ **结果逐位相同**。
      // 回代是延迟受限的（达成带宽只有流式读写的 9–13 %），3 条链分给 3 条线程正好把
      // "未决访存请求数"提高约 3 倍：实测单次回代快 2.0–2.8×（docs/perf.md §9）。
      // 结构不满足时 parallelComponents() 返回 1，循环只有 1 个迭代，行为与原来的 single 一致。
      {
        const auto t0 = Clock::now();
        const int nComp = ctx.solver->parallelComponents();
#ifdef _OPENMP
#pragma omp for schedule(static, 1)
#endif
        for (int c = 0; c < nComp; ++c) {
          applyPinRhsComponent(m, ctx.b, c, nComp);
          ctx.solver->solveComponent(c, ctx.b, ctx.xSolution);
        }
        // 计时口径随构造一起变了（**引用本桶数字时必须留意**）：
        //   原来 = 唯一执行者自己的执行时间（不含 barrier 等待）；
        //   现在 = 本阶段的**墙钟**（含等最慢的那条线程），与散射/解包桶口径一致。
        // 记账在 barrier 之后由**单线程**完成：求解器的 stats_ 不能由并发的分量线程累加
        // （那是数据竞争），所以走 noteSolveStage()，见 IGlobalSolver.h 的说明。
        if (threadId() == 0) {
          const double stageSeconds = secondsSince(t0);
          ctx.times.solve += stageSeconds;
          ctx.solver->noteSolveStage(stageSeconds);
        }
      }

      // 4f) 解包 + 收敛扫描 + previous 更新 + **非线性残差**（逐顶点融合成一趟）
      //
      //     部分值用**区域内声明的变量**（每线程私有），再写进共享槽位后统一合并 ——
      //     max 是精确运算，合并顺序不影响结果，因此与串行扫描逐位相同。
      //
      //     scale 的定义：取"本子步已发生的位移量级"，
      //     而不是坐标量级 —— 后者会让判据失去平移不变性、并在高刚度下过早放行。
      //
      //     残差为什么能塞进这一趟：它逐顶点 gather（见 vertexResidual），端点位置
      //     取 xSolution —— 与这一趟写的 m.positions[i] 是同一个值，所以"邻居还没解包"
      //     不构成读序依赖，不需要额外的 barrier。这取代了原来 single 里的串行整趟
      //     O(E) 遍历（那是 Phase 3 要消掉的最后一块串行工作）。
      // ω 递推（每迭代一个标量；chebEnabled 为假时整段不执行）
      const bool chebBlend = chebEnabled && (k > chebStart);
      if (chebBlend) {
        if (k == chebStart + 1) chebOmega = Scalar{2} / (Scalar{2} - chebRho * chebRho);
        else chebOmega = Scalar{4} / (Scalar{4} - chebRho * chebRho * chebOmega);
      }      Scalar localDiff = 0.0;
      Scalar localScale = 0.0;
      Scalar localResidual = 0.0;
      int localResidualVertex = -1;
      const auto tScan = Clock::now();
      const Scalar* xSol = ctx.xSolution.data();
      const ConstraintColoring& adj = m.constraintColoring();
      const BendAdjacency& bendAdj = m.bendAdjacency();
#ifdef _OPENMP
#pragma omp for nowait
#endif
      for (long long v = 0; v < n; ++v) {
        const std::size_t i = static_cast<std::size_t>(v);
        const Vec3 xNew{ctx.xSolution[3 * i + 0], ctx.xSolution[3 * i + 1], ctx.xSolution[3 * i + 2]};
        // Chebyshev：把未加速解与"上一步位形"按 ω 合成（ω = 1 时走下面 else 分支，
        // 逐字保持朴素 PD 的算术；previous[i] 正好就是 x_{k−1}，与 Wang 的 prev_X 同义）。
        if (chebBlend) {
          const Vec3 xPrev = previous[i];
          const Vec3 xPut{xPrev.x + (xNew.x - xPrev.x) * chebOmega,
                          xPrev.y + (xNew.y - xPrev.y) * chebOmega,
                          xPrev.z + (xNew.z - xPrev.z) * chebOmega};
          m.positions[i] = xPut;
          localDiff = std::max(localDiff, maxAbsComponent(xPut - xPrev));
          previous[i] = xPut;   // 下一迭代的 x_{k−1} 是**加速后**的位形
        } else {
          m.positions[i] = xNew;
          localDiff = std::max(localDiff, maxAbsComponent(xNew - previous[i]));
          previous[i] = xNew;
        }
        localScale = std::max(localScale, maxAbsComponent(ctx.predicted[i] - ctx.positionsBeforeStep[i]));

        if (gatherResidual) {
          const Scalar mass = m.masses[i];
          // 与串行版同样的跳过条件：pinned 顶点不是未知量（其行被整行覆盖），
          // 零质量顶点不参与（除零）。
          if (!m.isPinned(static_cast<int>(v)) && mass > 0.0) {
            const Scalar r = vertexResidual(m, adj, bendAdj, xSol, ctx.predicted, mass, invH2,
                                            static_cast<int>(v));
            // 严格大于：相同残差时保留更小顶点号（与合并顺序共同保证确定性）。
            if (r > localResidual) {
              localResidual = r;
              localResidualVertex = static_cast<int>(v);
            }
          }
        }
      }
#ifdef _OPENMP
      if (threadId() < mergeSlots) {
        const std::size_t slot = static_cast<std::size_t>(threadId());
        partialDiff[slot] = localDiff;
        partialScale[slot] = localScale;
        partialResidual[slot] = localResidual;
        partialResidualVertex[slot] = localResidualVertex;
      }
      // 显式 barrier：上面的槽位必须全部写完，下面才算得对。
#pragma omp barrier
#pragma omp single
#endif
      {
        Scalar diff = 0.0;
        Scalar scale = 0.0;
        Scalar residualParallel = 0.0;
        int residualVertex = -1;
        for (int t = 0; t < mergeSlots; ++t) {
          diff = std::max(diff, partialDiff[static_cast<std::size_t>(t)]);
          scale = std::max(scale, partialScale[static_cast<std::size_t>(t)]);
          // 合并规则：取更大的残差；相等时取更小的顶点号（槽位按线程号升序）。
          const Scalar r = partialResidual[static_cast<std::size_t>(t)];
          const int rv = partialResidualVertex[static_cast<std::size_t>(t)];
          if (rv >= 0 && (r > residualParallel || (r == residualParallel && residualVertex < 0))) {
            residualParallel = r;
            residualVertex = rv;
          }
        }
        scale = std::max(scale, length(cfg.gravity) * h * h);
        const Scalar rel = (scale > 0.0) ? diff / scale : diff;

        // 注意 absTolerance 用 "> 0" 判定禁用：若只写成 diff <= cfg.absTolerance，
        // 传 0 会因 diff==0 时恒成立而立刻退出（tests/chain 里记录过这个坑）。
        const bool absHit = (cfg.absTolerance > 0.0) && (diff <= cfg.absTolerance);
        const bool relHit = (cfg.relTolerance > 0.0) && (rel <= cfg.relTolerance);

        // 非线性残差判据：直接量"离本子步不动点还有多远"。
        // 位移类判据在高刚度下会失真 —— 更新量是被 L 预条件后的梯度，
        // 小更新不保证小物理残差（κ=1e5 时 40 次迭代只完成切向运动的 0.46%，
        // 而单步位移看起来已经很小）。
        // 注意：调试开关 PD_DEBUG_RESIDUAL 要求**即使判据关闭也要计算残差**，
        // 否则"关掉判据"的对照运行里就看不到残差读数了。
        Scalar residual = gatherResidual ? residualParallel : Scalar{-1};
        bool resHit = false;
        if (serialResidual) {
          // 串行口径（边序累加，每条边只算一次投影）。调试时它还会打印 worst 顶点的逐项分解，
          // 并且此时并行口径也被算了出来 ⇒ 顺手做一次两口径自查。
          const Scalar serial = nonlinearResidual(ctx);
          if (gatherResidual && serial != residual) {
            const Scalar mag = std::max(std::fabs(serial), std::fabs(residual));
            const Scalar deviation = (mag > 0.0) ? std::fabs(serial - residual) / mag : 0.0;
            if (deviation > 1e-9) {
              std::fprintf(stderr,
                           "[residual] 口径不一致：串行 %.9e 并行 %.9e（相对差 %.2e，"
                           "并行的 worst 顶点 %d）\n"
                           "  两者应逐位相同；不一致说明某一边被改坏了。\n",
                           serial, residual, deviation, residualVertex);
            }
          }
          residual = serial;
        }
        if (cfg.residualTolerance > 0.0) {
          const Scalar gScale = std::max<Scalar>(length(cfg.gravity), 1.0e-12);
          resHit = (residual <= cfg.residualTolerance * gScale);
        }
        ctx.lastResidual = residual;
        // 这一趟的时间含 barrier 等待（所有线程都要到齐才算完），因此是"该阶段的墙钟"。
        ctx.times.unpackScan += secondsSince(tScan);

        // ---- 收敛判定 ----
        //
        // **启用残差判据时，残差达标是收敛的必要条件**，位移判据不得单独放行。
        // 理由：位移判据只能说"这一步没怎么动"，而"没怎么动"既可能是收敛，
        // 也可能是停滞（高刚度下更新量是被 L 预条件后的梯度，小更新不保证小残差）。
        // 若让 `absHit || relHit` 独立触发退出，就会出现"残差 4.85 m/s² 而门槛
        // 要求 9.8e-3 却判收敛"的情形 —— 这个组合曾被实测复现。
        // 因此：residualTolerance > 0 时只认 resHit；为 0（未启用）时才退回位移判据。
        const bool converged = (cfg.residualTolerance > 0.0) ? resHit : (absHit || relHit);
        ctx.converged = converged;
        if (converged) ctx.earlyExitCount += 1;  // 累加只能发生在一个线程上
        ctx.iterationsUsed = k + 1;

        PD_TRACE("迭代 %d: |Δx|∞=%.6g 相对=%.6g 残差=%.6g m/s^2", k + 1, diff, rel, residual);
      }
      // single 的隐式 barrier 在所有线程上保证了同一个 `ctx.converged`，
      // 因此这里的 break 是**一致**的（OpenMP 要求各线程以相同方式退出循环）。
      if (ctx.converged) break;
    }

    // ---------------------------------------------------------------
    // 5) 速度更新与阻尼
    //
    //    v_{n+1} = (1 - k_d) · (x_{n+1} - x_n) / h
    //
    //    这就是**全部**的速度更新：位置差除以 h，再乘一个几何衰减因子 (1-k_d)。
    //    被 pin 的顶点速度清零，并把位置数值上再钉一次，保证严格等于 pinPositions
    //    （避免在后续帧累积无意义的动量）。
    //
    //    为什么这里不需要"惯性残差补偿"之类的额外项：
    //    预测步 x̂ = x_n + h·v_n + h²·g 已经把上一步速度带进来了，因此
    //    v_n 的信息只应通过 x̂ 进入本步；再往速度更新里加一项就会把 v_n 计两遍
    //    （曾经按"补偿惯性残差"的思路试过，是错的，已撤掉）。
    //
    //    静止时的固定点：静止意味着 x_{n+1} = x_n = y 且 v = 0，于是 x̂ = y - h²g，
    //    不动点方程 m(y-x̂)/h² + κ(y-ℓ) = 0 退化为纯力平衡
    //        κ(y - ℓ) = m·g   ⇒   y = ℓ - m·g/κ
    //    （重力向量朝 -y，故质量挂在静止长度**下方**，这是能量的极小值）。
    //    注意该结果与 k_d 和 h 都无关 —— 它们是离散格式参数，不影响静力平衡。
    //    验收见 `pd_check` 第 3 项与 `tests/chain/spring_vertical.cpp` 第 3 项。
    // ---------------------------------------------------------------
    {
      const auto t0 = Clock::now();
      const Scalar damp = 1.0 - cfg.velocityDamping;
      // nowait：这是区域里最后一趟循环，区域结束时的隐式 barrier 已足够。
#ifdef _OPENMP
#pragma omp for nowait
#endif
      for (long long v = 0; v < n; ++v) {
        const std::size_t i = static_cast<std::size_t>(v);
        Vec3 vNew = (m.positions[i] - ctx.positionsBeforeStep[i]) * (damp * invH);
        if (m.isPinned(static_cast<int>(v))) {
          vNew = Vec3{};
          m.positions[i] = m.pinPositions[i];  // 数值上再钉一次，保证严格相等
        }
        m.velocities[i] = vNew;
      }
      if (threadId() == 0) ctx.times.velocity += secondsSince(t0);
    }
  }  // 并行区域结束（隐式 barrier：所有线程都已写完）

  ctx.lastStepSeconds = secondsSince(tStep0);
  ctx.times.total += ctx.lastStepSeconds;
  return ctx.iterationsUsed;
}

void stepFrame(SimContext& ctx) {
  const int substeps = std::max(1, ctx.config.substepsPerFrame);
  for (int s = 0; s < substeps; ++s) stepOnce(ctx);
}

Scalar surrogateEnergy(const SimContext& ctx) {
  Scalar e = 0.0;
  for (const auto& edge : ctx.mesh.edges) {
    const Scalar len = length(ctx.mesh.positions[static_cast<std::size_t>(edge.a)] -
                              ctx.mesh.positions[static_cast<std::size_t>(edge.b)]);
    const Scalar d = len - edge.restLength;
    e += 0.5 * edge.stiffness * d * d;
  }
  // 弯曲能量 (k/2)·‖A_s x‖²。**必须一起报**：它是真实存在的一项弹性势能，
  // 只报距离约束那部分会让"能量守恒/单调下降"这类读数在启用弯曲后失去意义。
  // 没有弯曲约束时下面这个循环一次都不进 ⇒ 默认关闭时与改动前逐位相同。
  for (const auto& s : ctx.mesh.bends) {
    // 写法与 `bendForceAt` 里逐字同形：`a - 2*b + c`（**不写成** `(a - b) + (c - b)`，
    // 那样在浮点上不是同一个表达式，位级不同）。
    const Vec3& pa = ctx.mesh.positions[static_cast<std::size_t>(s.a)];
    const Vec3& pb = ctx.mesh.positions[static_cast<std::size_t>(s.b)];
    const Vec3& pc = ctx.mesh.positions[static_cast<std::size_t>(s.c)];
    const Vec3 second{pa.x - 2.0 * pb.x + pc.x, pa.y - 2.0 * pb.y + pc.y,
                      pa.z - 2.0 * pb.z + pc.z};
    e += 0.5 * s.stiffness * lengthSquared(second);
  }
  return e;
}

Scalar totalEnergy(const SimContext& ctx) {
  Scalar kinetic = 0.0;
  Scalar gravityPotential = 0.0;
  const Vec3 g = ctx.config.gravity;
  for (int v = 0; v < ctx.mesh.vertexCount(); ++v) {
    const std::size_t i = static_cast<std::size_t>(v);
    const Scalar mass = ctx.mesh.masses[i];
    kinetic += 0.5 * mass * lengthSquared(ctx.mesh.velocities[i]);
    // 重力势能 U = -m (g · x)
    gravityPotential -= mass * dot(g, ctx.mesh.positions[i]);
  }
  return kinetic + surrogateEnergy(ctx) + gravityPotential;
}

Scalar maxSpeed(const SimContext& ctx) {
  Scalar worst = 0.0;
  for (const auto& v : ctx.mesh.velocities) worst = std::max(worst, length(v));
  return worst;
}

}  // namespace pd
