// _verify/check.cpp
// 唯一的验收程序：只回答"对不对"。不打印推导、不打印中间量。
//
// 检查项（每一项独立，全部对照解析解）：
//   1. 投影长度是否恰好等于静止长度
//   2. 单步解是否等于闭式解 (m ŷ/h² + κℓ)/(m/h² + κ)
//   3. 悬挂弹簧的静止平衡是否等于 ℓ + m g/κ
//   4. pinned 顶点是否严格不动
//   5. 自由落体一步是否为 -h²g
//   6. 布料长时间积分是否有界
//   7. 预分解是否只发生一次（复用是否生效）
//
// 运行：build\Release\pd_check.exe        （退出码 0 = 全部正确，1 = 有错）
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <vector>

#include "core/assemble/Assembler.h"
#include "core/energy/DistanceTerm.h"
#include "core/mesh/Mesh.h"
#include "core/sim/Integrator.h"
#include "core/sim/Scene.h"
#include "core/solver/EigenDirectSolver.h"

using namespace pd;

namespace {

int g_failed = 0;

void item(const char* name, bool ok, double got, double want, double tol) {
  std::printf("  %-34s %s", name, ok ? "OK" : "NG");
  if (!ok) {
    std::printf("     实测 %.10g  期望 %.10g  容差 %.3g", got, want, tol);
    ++g_failed;
  }
  std::printf("\n");
}

const Scalar kH = 1.0 / 120.0;
const Scalar kM = 0.05;
const Scalar kEll = 0.1;
const Scalar kK = 1.0e4;
const Scalar kG = 9.81;

SimContext makeSpring(Scalar y0, Scalar kd) {
  SimContext ctx;
  ctx.config.dt = kH;
  ctx.config.gravity = Vec3{0.0, -kG, 0.0};
  ctx.config.stiffness = kK;
  ctx.config.velocityDamping = kd;
  ctx.config.maxIterations = 40;
  ctx.config.relTolerance = 0.0;
  ctx.config.absTolerance = 1e-18;
  Mesh& m = ctx.mesh;
  m.positions = {Vec3{0.0, 0.0, 0.0}, Vec3{0.0, y0, 0.0}};
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
  // 1. 投影长度
  {
    const Vec3 d = DistanceTerm::projectOne(Vec3{0, 0, 0}, Vec3{0, 0.31, 0}, kEll);
    item("1 投影长度 == 静止长度", std::fabs(length(d) - kEll) < 1e-14, length(d), kEll, 1e-14);
  }

  // 2. 单步解 == 闭式解
  {
    const Scalar y0 = 0.12;
    SimContext ctx = makeSpring(y0, 0.0);
    stepOnce(ctx);
    const Scalar yHat = y0 - kH * kH * kG;
    const Scalar mh2 = kM / (kH * kH);
    // 自由顶点是 1，A_c x = x_a - x_b（a=0 是 pin），d_c.y = -ℓ ⇒ 散射 +κℓ
    const Scalar want = (mh2 * yHat + kK * kEll) / (mh2 + kK);
    const Scalar got = ctx.mesh.positions[1].y;
    item("2 单步解 == 闭式解", std::fabs(got - want) < 1e-12, got, want, 1e-12);
  }

  // 3. 悬挂弹簧平衡 == ℓ - m g/κ
  //
  // 本设置的势能：E(y) = (κ/2)(y-ℓ)² + m g y   （重力向量 g_vec = (0,-g,0)，g 为大小）
  //   dE/dy = κ(y-ℓ) + m g = 0  ⇒  y = ℓ - m g/κ
  // 注意重力项梯度是 +m g（朝 -y），所以平衡在静止长度**下方**。
  // 这一点曾经被写反（写成 ℓ + m g/κ），导致把正确实现误判为 bug，特此注明。
  {
    SimContext ctx = makeSpring(0.12, 0.0);
    for (int s = 0; s < 20000; ++s) stepOnce(ctx);
    const Scalar want = kEll - kM * kG / kK;
    const Scalar got = ctx.mesh.positions[1].y;
    item("3 静止平衡 == ℓ - m g/κ", std::fabs(got - want) < 1e-6 * kEll, got, want, 1e-6 * kEll);
  }

  // 4. pinned 严格不动
  {
    SimContext ctx = makeSpring(0.12, 0.0);
    Scalar worst = 0.0;
    for (int s = 0; s < 200; ++s) {
      stepOnce(ctx);
      worst = std::max(worst, std::fabs(ctx.mesh.positions[0].y));
    }
    item("4 pinned 顶点严格不动", worst == 0.0, worst, 0.0, 0.0);
  }

  // 5. 自由落体
  {
    SimContext ctx;
    ctx.config.dt = kH;
    ctx.config.gravity = Vec3{0.0, -kG, 0.0};
    ctx.config.velocityDamping = 0.0;
    ctx.config.maxIterations = 4;
    Mesh& m = ctx.mesh;
    m.positions = {Vec3{0.0, 0.0, 0.0}};
    m.restPositions = m.positions;
    m.velocities = {Vec3{}};
    m.masses = {1.0};
    m.pinned = {0};
    m.pinPositions = m.positions;
    ensureBuffers(ctx);
    stepOnce(ctx);
    const Scalar want = -kH * kH * kG;
    const Scalar got = ctx.mesh.positions[0].y;
    item("5 自由落体一步 == -h²g", std::fabs(got - want) < 1e-14, got, want, 1e-14);
  }

  // 6. 布料有界
  {
    SceneConfig cfg;
    cfg.gridNx = 12;
    cfg.gridNy = 12;
    cfg.gridSpacing = 0.05;
    cfg.dt = kH;
    cfg.substepsPerFrame = 2;
    cfg.stiffness = kK;
    cfg.maxIterations = 10;
    cfg.relTolerance = 1e-3;
    SimContext cloth = makeScene(cfg);
    for (int s = 0; s < 400; ++s) stepOnce(cloth);
    Scalar worst = 0.0;
    bool finite = true;
    for (const auto& p : cloth.mesh.positions) {
      if (!std::isfinite(p.y)) finite = false;
      worst = std::max(worst, std::fabs(p.y));
    }
    item("6 布料有界", finite && worst < 1.0, worst, 0.0, 1.0);
  }

  // 7. 预分解复用（数值分解次数应为 1）
  {
    SceneConfig cfg;
    cfg.gridNx = 20;
    cfg.gridNy = 20;
    cfg.gridSpacing = 0.02;
    cfg.dt = kH;
    cfg.stiffness = kK;
    cfg.maxIterations = 10;
    SimContext cloth = makeScene(cfg);
    for (int s = 0; s < 200; ++s) stepOnce(cloth);
    const int got = cloth.factorizeCount;
    item("7 数值分解次数 == 1", got == 1, static_cast<double>(got), 1.0, 0.0);
  }

  // 8. 弹性力的方向（能量梯度判据，与方向约定无关）
  //
  // 一维距离约束：顶点 0 在原点、顶点 1 在 +y，r = |y|，U = (κ/2)(r-ℓ)²。
  //
  //   真值：力的方向由"受力点减施力点"给出
  //       F_true = -κ(r-ℓ) · unit(x_free - x_pin)
  //
  //   代码：从行方程反推净弹性项。自由顶点的行方程是
  //       L·y = b        ⇒       (m/h²)·y + κ·y = (m/h²)·ŷ + 散射
  //   把惯性项移到一边，得到"净弹性项 = 散射 − κ·y"：
  //       自由端是 b ⇒ 散射 = -κ·d_c   ⇒  净弹性项 = -κ·d_c - κ·y
  //       自由端是 a ⇒ 散射 = +κ·d_c   ⇒  净弹性项 = +κ·d_c - κ·y
  //
  //   **注意必须走"惯性 = 散射 − κ·y"这一步**：κ·y 是左端算子，不能直接与右端相加
  //   （早先这里写成 κ·x_free ∓ κ·d_c 就是漏了移项，导致误报 NG）。
  {
    bool allOk = true;
    double worstGot = 0.0, worstWant = 0.0;
    const bool verbose = (std::getenv("PD_CHECK_VERBOSE") != nullptr);
    const Scalar ys[] = {0.06, 0.08, 0.0999, 0.1001, 0.12, 0.15};
    for (Scalar y1 : ys) {
      // ---- 情形 1：顶点 0 是 pin、顶点 1 自由（边 (0,1)，自由端是 b）----
      {
        Mesh m;
        m.positions = {Vec3{0, 0, 0}, Vec3{0, y1, 0}};
        m.pinned = {1, 0};
        m.edges.push_back(Edge{0, 1, kEll, kK});
        std::vector<Vec3> t;
        DistanceTerm::project(m, t);
        const Vec3 xFree = m.positions[1];
        const Vec3 xPin = m.positions[0];
        const Scalar r = length(xFree - xPin);
        const Vec3 trueForce = (xFree - xPin) * (-kK * (r - kEll) / r);
        // 净弹性项 = 散射(-κ d_c) − κ·y
        const Vec3 codeForce = (t[0] * (-kK)) - (xFree * kK);
        if (verbose) {
          std::printf("      [v] y=%.4f  free=b  d_c.y=%.4f  scatter.y=%+.1f  κy=%+.1f  net=%+.1f  true=%+.1f\n",
                      y1, t[0].y, (-t[0].y * kK), (xFree.y * kK), codeForce.y, trueForce.y);
        }
        if (std::fabs(codeForce.y - trueForce.y) > 1e-6) allOk = false;
        worstGot = codeForce.y;
        worstWant = trueForce.y;
      }
      // ---- 情形 2：顶点 1 是 pin、顶点 0 自由（边 (0,1)，自由端是 a）----
      {
        Mesh m;
        m.positions = {Vec3{0, y1, 0}, Vec3{0, 0, 0}};
        m.pinned = {0, 1};
        m.edges.push_back(Edge{0, 1, kEll, kK});
        std::vector<Vec3> t;
        DistanceTerm::project(m, t);
        const Vec3 xFree = m.positions[0];
        const Vec3 xPin = m.positions[1];
        const Scalar r = length(xFree - xPin);
        const Vec3 trueForce = (xFree - xPin) * (-kK * (r - kEll) / r);
        // 净弹性项 = 散射(+κ d_c) − κ·y
        const Vec3 codeForce = (t[0] * kK) - (xFree * kK);
        if (verbose) {
          std::printf("      [v] y=%.4f  free=a  d_c.y=%.4f  scatter.y=%+.1f  κy=%+.1f  net=%+.1f  true=%+.1f\n",
                      y1, t[0].y, (t[0].y * kK), (xFree.y * kK), codeForce.y, trueForce.y);
        }
        if (std::fabs(codeForce.y - trueForce.y) > 1e-6) allOk = false;
        worstGot = codeForce.y;
        worstWant = trueForce.y;
      }
    }
    item("8 弹性力符号 == -dU/dy", allOk, worstGot, worstWant, 1e-6);
  }

  std::printf("\n  %s（失败 %d 项）\n", g_failed == 0 ? "全部正确" : "存在错误", g_failed);
  return g_failed == 0 ? 0 : 1;
}
