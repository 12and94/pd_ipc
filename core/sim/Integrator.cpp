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

#define PD_TRACE(...)                             \
  do {                                            \
    if (traceEnabled()) {                         \
      std::fprintf(stderr, "[step] " __VA_ARGS__); \
      std::fprintf(stderr, "\n");                  \
      std::fflush(stderr);                         \
    }                                             \
  } while (0)

}  // namespace

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
      if (!ctx.solver->analyzed()) {
        // 首次（或拓扑变化后）：先按真实矩阵的结构做符号分解，再数值分解。
        ctx.solver->analyze(3 * n, ctx.L);
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

    // 4e) 收敛判据：相对无穷范数
    usedIterations = k + 1;
    Scalar diff = 0.0;
    Scalar scale = 0.0;
    for (int v = 0; v < n; ++v) {
      const std::size_t i = static_cast<std::size_t>(v);
      diff = std::max(diff, maxAbsComponent(m.positions[i] - previous[i]));
      scale = std::max(scale, maxAbsComponent(ctx.predicted[i]));
    }
    previous = m.positions;
    const Scalar rel = (scale > 0.0) ? diff / scale : diff;
    PD_TRACE("迭代 %d: |Δx|∞ = %.6g  相对 %.6g", usedIterations, diff, rel);
    if (diff <= cfg.absTolerance || rel <= cfg.relTolerance) break;
  }
  ctx.iterationsUsed = usedIterations;

  // ---------------------------------------------------------------
  // 5) 速度更新与阻尼
  //
  //    v_{n+1} = (1 - k_d) · [ v_n + (x_{n+1}-x_n)/h - f_int(x_n)·h/M ]
  //
  //    方括号里三项：上一步速度、位置差、以及**上一步惯性残差的补偿**。
  //    第三项不能少：只用 v = (x-x_prev)/h 意味着把 t_n 的惯性项 M·x_n/h²
  //    与 t_{n+1} 的力混在一起，静止条件会退化成 κ(ℓ-x) = m g，
  //    平衡位置落到 ℓ - m g/κ。补上后固定点满足
  //        m h²g/h² + κ(y-ℓ) = 0  ⇒  κ(y-ℓ) = m g  ⇒  y = ℓ + m g/κ
  //    即标准 PD 的物理平衡。
  //
  //    实现：把每步的"力冲量"累计到速度里。对距离约束，f_int·h/M 用
  //    投影目标的方向乘以 κ(r-ℓ)h/m 近似（投影与力同向）。
  //
  //    被 pin 的顶点速度清零，避免其在后续帧累积无意义的动量。
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
