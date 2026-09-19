// _verify/one_step_trace.cpp
// 把"从 y* = ℓ+mg/κ、v=0 出发走一步"的每个中间量打印出来，
// 与解析式 y1 = (m ŷ/h² + κℓ)/(m/h² + κ) 并排对比。
//
// 判据：解析式在 y* 处给出 y1 = y*（不动点）。若代码给出别的值，差在哪一项一看即知。
#include <cstdio>

#include "core/assemble/Assembler.h"
#include "core/energy/DistanceTerm.h"
#include "core/sim/Integrator.h"
#include "core/sim/Scene.h"
#include "core/solver/EigenDirectSolver.h"

using namespace pd;

int main() {
  const Scalar h = 1.0 / 120.0, m = 0.05, ell = 0.1, k = 1.0e4, g = 9.81;
  const Scalar gAcc = 9.81;                       // |g|
  const Scalar yStar = ell + m * gAcc / k;        // 候选平衡 = 0.10004905

  SimContext ctx;
  ctx.config.dt = h;
  ctx.config.gravity = Vec3{0.0, -gAcc, 0.0};
  ctx.config.stiffness = k;
  ctx.config.velocityDamping = 0.0;
  ctx.config.maxIterations = 1;   // 只做一次迭代，便于逐项核对
  ctx.config.relTolerance = 0.0;
  ctx.config.absTolerance = 0.0;

  Mesh& mesh = ctx.mesh;
  mesh.positions = {Vec3{0.0, 0.0, 0.0}, Vec3{0.0, yStar, 0.0}};
  mesh.restPositions = mesh.positions;
  mesh.velocities = {Vec3{}, Vec3{}};
  mesh.masses = {m, m};
  mesh.pinned = {1, 0};
  mesh.pinPositions = mesh.positions;
  mesh.edges.push_back(Edge{0, 1, ell, k});
  mesh.buildSparsityPattern();
  ensureBuffers(ctx);

  // ---- 解析预期 ----
  const Scalar mh2 = m / (h * h);
  const Scalar yHat = yStar - h * h * gAcc;   // x̂ = x + h·0 + h²g（g 沿 -y）
  const Scalar dcY = -(ell);                   // A_c x = x_a - x_b，a=0(pin), b=1(free) ⇒ d_c.y = -ℓ
  const Scalar scatterAnalytic = -k * dcY;     // b_1 -= κ d_c ⇒ 自由端收到 -κ d_c.y = +κℓ
  const Scalar y1Analytic = (mh2 * yHat + scatterAnalytic) / (mh2 + k);

  std::printf("输入: y* = %.10g, v = 0, h = 1/120, m = %.4g, ℓ = %.4g, κ = %.6g\n\n", yStar, m, ell, k);
  std::printf("解析:\n");
  std::printf("  预测 ŷ            = y* - h²|g| = %.12g\n", yHat);
  std::printf("  d_c.y             = %.12g\n", dcY);
  std::printf("  自由端散射        = -κ·d_c.y = %+.12g\n", scatterAnalytic);
  std::printf("  L_11 = m/h² + κ   = %.12g\n", mh2 + k);
  std::printf("  y1 = (m/h²·ŷ + 散射)/L_11 = %.12g   （应等于 y* = %.12g）\n\n", y1Analytic, yStar);

  // ---- 代码的中间量 ----
  std::vector<Vec3> targets;
  DistanceTerm::project(mesh, targets);
  Eigen::SparseMatrix<Scalar> L;
  assembleLeftHandSide(mesh, h, ctx.config.velocityDamping, L);
  const std::vector<Vec3> pred = {mesh.positions[0] + ctx.config.gravity * (h * h),
                                  mesh.positions[1] + ctx.config.gravity * (h * h)};
  Eigen::VectorXd b;
  assembleInertialRhs(mesh, pred, h, ctx.config.velocityDamping, b);
  const Scalar bInertial = b[4];
  DistanceTerm::scatterInto(mesh, targets, b.data(), static_cast<std::size_t>(b.size()));
  const Scalar scatterCode = b[4] - bInertial;

  std::printf("代码:\n");
  std::printf("  预测 ŷ            = %.12g   %s\n", pred[1].y,
              std::fabs(pred[1].y - yHat) < 1e-15 ? "[与解析一致]" : "[不一致]");
  std::printf("  d_c.y             = %.12g   %s\n", targets[0].y,
              std::fabs(targets[0].y - dcY) < 1e-15 ? "[与解析一致]" : "[不一致]");
  std::printf("  自由端散射        = %+.12g   %s\n", scatterCode,
              std::fabs(scatterCode - scatterAnalytic) < 1e-9 ? "[与解析一致]" : "[不一致]");
  std::printf("  L(4,4)            = %.12g   %s\n", L.coeff(4, 4),
              std::fabs(L.coeff(4, 4) - (mh2 + k)) < 1e-9 ? "[与解析一致]" : "[不一致]");
  std::printf("  b[4]              = %.12g\n", b[4]);
  std::printf("  y1 = b[4]/L(4,4)  = %.12g\n\n", b[4] / L.coeff(4, 4));

  // ---- 真正调 stepOnce，与上面的手工推导对照 ----
  SimContext ctx2;
  ctx2.config = ctx.config;
  {
    Mesh& mm = ctx2.mesh;
    mm.positions = {Vec3{0.0, 0.0, 0.0}, Vec3{0.0, yStar, 0.0}};
    mm.restPositions = mm.positions;
    mm.velocities = {Vec3{}, Vec3{}};
    mm.masses = {m, m};
    mm.pinned = {1, 0};
    mm.pinPositions = mm.positions;
    mm.edges.push_back(Edge{0, 1, ell, k});
    mm.buildSparsityPattern();
    ensureBuffers(ctx2);
  }
  stepOnce(ctx2);
  std::printf("stepOnce（1 次迭代）y1 = %.12g\n", ctx2.mesh.positions[1].y);
  std::printf("诊断结论：逐项对照上表的 [一致/不一致] 标记，差值就落在标为不一致的那一项上\n");
  return 0;
}
