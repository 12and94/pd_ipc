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
/// 其中距离约束的弹性力 f_i^int 等于"散射进右端的那一项"的相反数：
///     f_i^int = -Σ_c κ_c d_c(x)   （a 端 +κd_c、b 端 -κd_c）
/// 于是每个顶点的不平衡加速度为
///     r_i = |(m_i/h²)(x_i - x̂_i) + Σ_c κ_c d_c| / m_i
/// 注意这里**不含重力项**：重力已经通过 x̂ = x + hv + h²g 进入，且只进入一次
/// （在右端再补 -Mg 就是把重力算两遍，本项目踩过这个坑）。x 到达不动点时
/// 上式整体为 0，因此它就是"离不动点还有多远"的度量。
///
/// 为什么不能用"相邻迭代位移"代替它：位移小只说明这一步没怎么动。
/// 高刚度下更新量是被 L 预条件后的梯度，小更新不保证小物理残差
/// （κ=1e5 时 40 次迭代只完成切向运动的 0.46%，而位移看起来"已经很小"）。
}  // namespace

Scalar nonlinearResidual(const SimContext& ctx) {
  const Mesh& m = ctx.mesh;
  const int n = m.vertexCount();
  const Scalar invH2 = 1.0 / (ctx.config.dt * ctx.config.dt);

  // 约束力的正确形式：f^int = κ(AᵀA x − Aᵀ d) —— **两部分都要**。
  //   推导：真实弹性力 κ(r−ℓ)·unit(x_a−x_b) = κA x − κ ℓ·unit = κA x − κ d_c
  //   按顶点 a、b 展开后，a 端得 +κ(x_a−x_b) − κ d_c，b 端得 −κ(x_a−x_b) + κ d_c。
  //   **只取 κd_c 一侧是错的**（那是散射项，不是力）：本项目实测只取 κd_c 时
  //   残差为 7.07e6 m/s²，而应变 1e-4 的弹簧其力只有 0.2 N 量级 —— 差了 4 个数量级。
  //   原因：x = x̂、v = 0 时 κAᵀA x 与 κAᵀd 各自都是 O(κℓ)，只有两者之差才是 O(κ(r−ℓ))。
  static std::vector<Scalar> force;
  force.assign(static_cast<std::size_t>(n) * 3, Scalar{0});
  for (std::size_t c = 0; c < m.edges.size(); ++c) {
    const Edge& e = m.edges[c];
    const Vec3 rel = m.positions[static_cast<std::size_t>(e.a)] -
                     m.positions[static_cast<std::size_t>(e.b)];  // A_c x
    const Vec3 kd = ctx.targets[c] * e.stiffness;                // κ_c d_c
    const Vec3 contrib = rel * e.stiffness - kd;                 // κ(A_c x − d_c)
    const std::size_t ia = static_cast<std::size_t>(e.a) * 3;
    const std::size_t ib = static_cast<std::size_t>(e.b) * 3;
    force[ia + 0] += contrib.x;
    force[ia + 1] += contrib.y;
    force[ia + 2] += contrib.z;
    force[ib + 0] -= contrib.x;
    force[ib + 1] -= contrib.y;
    force[ib + 2] -= contrib.z;
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

void resetStageTimes(SimContext& ctx) { ctx.times = StageTimes{}; }

int stepOnce(SimContext& ctx) {
  const auto tStep0 = Clock::now();

  Mesh& m = ctx.mesh;
  const SceneConfig& cfg = ctx.config;
  const Scalar h = cfg.dt;
  const Scalar invH = 1.0 / h;
  const int n = m.vertexCount();

  // 自保护：任何构造路径（makeScene / 测试 / 调试程序）都可能忘记分配缓冲。
  // 这里统一保证尺寸正确，避免越界写（这是实际踩过的坑）。
  ensureBuffers(ctx);

  // 保存本步开始时的位置：速度更新必须用它，而不是从上式反推
  // （反推会引入一次额外的舍入，速度就不再是"精确的差分"）。
  ctx.positionsBeforePrevStep = ctx.positionsBeforeStep;
  ctx.positionsBeforeStep = m.positions;

  // ---------------------------------------------------------------
  // 1) 预测：x̂ = x + h v + h² g
  //    （曾用环境变量把这一项反号做过符号自检，已移除。）

  // ---------------------------------------------------------------
  {
    const auto t0 = Clock::now();

    for (int v = 0; v < n; ++v) {
      const std::size_t i = static_cast<std::size_t>(v);
      ctx.predicted[i] = m.positions[i] + m.velocities[i] * h + cfg.gravity * (h * h);
    }
    ctx.times.predict += secondsSince(t0);
  }
  PD_TRACE("预测完成 x̂[1] = %.12g", (n > 1) ? ctx.predicted[1].y : 0.0);

  // ---------------------------------------------------------------
  // 2) 判定是否需要重新组装 + 数值分解。
  //    正常运行时（拓扑/刚度/h/pin 都没变）这里**一次都不会触发**，
  //    这正是"预分解复用"的收益所在（docs/plan.md §2.3 不变量）。
  // ---------------------------------------------------------------
  {
    const auto t0 = Clock::now();
    const SolverStamp now = computeStamp(m, h, cfg.velocityDamping, ctx.topologyId);
    if (!ctx.stampValid || now != ctx.stamp) {
      assembleLeftHandSide(m, h, cfg.velocityDamping, ctx.L);
      // 符号分解的判定与数值分解**不同**：数值分解看数值戳（刚度/mass/h），
      // 符号分解还要看**稀疏结构**。pin 掩码变化会让"两端都自由"的耦合块消失，
      // 结构随之变化 —— 那时必须重新 analyze，否则 factorize 一个结构不同的矩阵后
      // solve 会访问未定义数据（见 IGlobalSolver.h，本项目因此踩过一次访问违例）。
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
    }
  }

  // ---------------------------------------------------------------
  // 3) 右端的惯性部分：b = (M/h²) x̂（保存一份基值，每迭代从它拷贝）
  // ---------------------------------------------------------------
  assembleInertialRhs(m, ctx.predicted, h, cfg.velocityDamping, ctx.bBase);

  // ---------------------------------------------------------------
  // 4) PD 迭代：局部步 → 散射 → 覆盖 pin → 全局步
  // ---------------------------------------------------------------
  std::vector<Vec3> previous = m.positions;
  int usedIterations = 0;

  for (int k = 0; k < cfg.maxIterations; ++k) {
    // 4a) 局部步：逐约束独立闭式投影（可并行，无耦合）
    {
      const auto t0 = Clock::now();
      DistanceTerm::project(m, ctx.targets);
      ctx.times.localStep += secondsSince(t0);
    }

    // 4b) 右端 = 惯性基值 + 约束散射
    ctx.b = ctx.bBase;
    {
      const auto t0 = Clock::now();
      DistanceTerm::scatterInto(m, ctx.targets, ctx.b.data(), static_cast<std::size_t>(ctx.b.size()));
      ctx.times.scatter += secondsSince(t0);
    }

    PD_TRACE("迭代 %d: bBase[4]=%.12g  b[4]=%.12g  L(4,4)=%.12g  L(4,1)=%.12g", k + 1,
             ctx.bBase.size() > 4 ? ctx.bBase[4] : 0.0, ctx.b.size() > 4 ? ctx.b[4] : 0.0,
             (ctx.L.rows() > 4 && ctx.L.cols() > 4) ? ctx.L.coeff(4, 4) : 0.0,
             (ctx.L.rows() > 4 && ctx.L.cols() > 1) ? ctx.L.coeff(4, 1) : 0.0);

    // 4c) 覆盖 pinned 行（pinned 位置严格等于把手位置）
    applyPinRhs(m, ctx.b);

    // 4d) 全局步：一次回代（矩阵不变，分解已复用）
    {
      const auto t0 = Clock::now();
      ctx.solver->solve(ctx.b, ctx.xSolution);
      ctx.times.solve += secondsSince(t0);
    }
    unpackPositions(ctx.xSolution, m.positions);

    // 4e) 收敛判据：以**本子步已发生的位移量级**为参照
    //
    //     scale = max( max|x̂ - xⁿ|∞ , h²|g| )
    //
    //     语义："本次迭代造成的位移，相对于本步已有的位移是否可忽略"。
    //
    //     为什么用位移量级而不是坐标量级（本项目修过的缺陷）：
    //     原先写作 diff / max|x̂|∞（坐标量级），放行的单步位移
    //     = relTolerance × 坐标量级，**与刚度无关**。高刚度下真实位移
    //     （量级 ρgℓ²/κ）远小于该门槛，于是第 1 次迭代就判"收敛"退出：
    //       · 每次子步留下残差，在后续帧累加 -> 布料缓慢漂移；
    //       · 判据不再平移不变 —— 同一场景平移到远处，坐标量级变大，
    //         判据突然变松，收敛质量随摆放位置改变。
    //     二者都是判据自身的问题，与 PD 算法无关。
    //
    //     下限 h²|g| 的作用：完全静止时 x̂ == xⁿ，scale 会退化为 0；
    //     而即使静止，一个子步内重力也必然要引入 ~h²g 量级的位移，
    //     用它做下限可避免"除以接近 0 的参照"而永不收敛。
    {
      usedIterations = k + 1;
      Scalar diff = 0.0;
      Scalar scale = 0.0;
      for (int v = 0; v < n; ++v) {
        const std::size_t i = static_cast<std::size_t>(v);
        diff = std::max(diff, maxAbsComponent(m.positions[i] - previous[i]));
        scale = std::max(scale, maxAbsComponent(ctx.predicted[i] - ctx.positionsBeforeStep[i]));
      }
      scale = std::max(scale, length(cfg.gravity) * h * h);
      previous = m.positions;
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
      Scalar residual = -1.0;
      bool resHit = false;
      const bool wantResidual = (cfg.residualTolerance > 0.0) || debugResidual();
      if (wantResidual) {
        residual = nonlinearResidual(ctx);
        if (cfg.residualTolerance > 0.0) {
          const Scalar gScale = std::max<Scalar>(length(cfg.gravity), 1.0e-12);
          resHit = (residual <= cfg.residualTolerance * gScale);
        }
      }
      ctx.lastResidual = residual;

      PD_TRACE("迭代 %d: |Δx|∞=%.6g 相对=%.6g 残差=%.6g m/s^2", usedIterations, diff, rel, residual);
      if (absHit || relHit || resHit) {
        ctx.earlyExitCount += 1;
        break;
      }
    }
  }
  ctx.iterationsUsed = usedIterations;

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
    for (int v = 0; v < n; ++v) {
      const std::size_t i = static_cast<std::size_t>(v);
      Vec3 vNew = (m.positions[i] - ctx.positionsBeforeStep[i]) * (damp * invH);
      if (m.isPinned(v)) {
        vNew = Vec3{};
        m.positions[i] = m.pinPositions[i];  // 数值上再钉一次，保证严格相等
      }
      m.velocities[i] = vNew;
    }
    ctx.times.velocity += secondsSince(t0);
  }

  ctx.lastStepSeconds = secondsSince(tStep0);
  ctx.times.total += ctx.lastStepSeconds;
  return usedIterations;
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
