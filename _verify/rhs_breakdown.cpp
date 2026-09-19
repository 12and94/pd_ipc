// 右端残差组成诊断：逐项打印 b 的构成，并与手工复算对照
#include <cstdio>

#include "core/assemble/Assembler.h"
#include "core/energy/DistanceTerm.h"
#include "core/sim/Integrator.h"
#include "core/sim/Scene.h"
#include "core/solver/EigenDirectSolver.h"

using namespace pd;

int main() {
  const Scalar h = 1.0 / 120.0, m = 0.05, ell = 0.1, k = 1.0e4, g = 9.81;
  const Scalar y0 = 0.12;

  SimContext ctx;
  ctx.config.dt = h;
  ctx.config.gravity = Vec3{0.0, -g, 0.0};
  ctx.config.stiffness = k;
  ctx.config.velocityDamping = 0.0;
  ctx.config.maxIterations = 1;
  ctx.config.relTolerance = 0.0;
  ctx.config.absTolerance = 0.0;

  Mesh& mesh = ctx.mesh;
  mesh.positions = {Vec3{0.0, 0.0, 0.0}, Vec3{0.0, y0, 0.0}};
  mesh.restPositions = mesh.positions;
  mesh.velocities = {Vec3{}, Vec3{}};
  mesh.masses = {m, m};
  mesh.pinned = {1, 0};
  mesh.pinPositions = mesh.positions;
  mesh.edges.push_back(Edge{0, 1, ell, k});
  mesh.buildSparsityPattern();
  ensureBuffers(ctx);

  std::printf("配置: y0=%.6f, h=1/120, m=%.4f, ell=%.4f, k=%.6g, g=%.4f\n", y0, m, ell, k, g);
  std::printf("顶点 0 = pinned 在原点, 顶点 1 = 自由, 边 (0,1)\n\n");

  // ---- 手工复现 stepOnce 的第 1~4 步 ----
  const Scalar invDt2 = 1.0 / (h * h);
  const Scalar mOverH2 = m * invDt2;

  std::printf("右端 b 的组成（自由顶点 1 的 y 分量，下标 4）\n");
  std::printf("  [1] 惯性项 (m/h^2) * 预测位置\n");
  std::printf("      m/h^2        = %.6f\n", mOverH2);

  // 预测
  const Vec3 pred1 = mesh.positions[1] + mesh.velocities[1] * h + ctx.config.gravity * (h * h);
  std::printf("      预测位置 yhat = %.10f   (= y0 + h*v + h^2*g)\n", pred1.y);
  std::printf("      => 惯性项     = %.10f\n\n", mOverH2 * pred1.y);

  // 投影
  std::vector<Vec3> targets;
  DistanceTerm::project(mesh, targets);
  std::printf("  [2] 约束散射项 -k * d_c\n");
  std::printf("      边 (a=%d, b=%d), A_c x = x_a - x_b = (%.4f - %.4f) = %.4f\n", mesh.edges[0].a,
              mesh.edges[0].b, mesh.positions[0].y, mesh.positions[1].y,
              mesh.positions[0].y - mesh.positions[1].y);
  std::printf("      d_c.y        = %.10f   (= ell * unit(A_c x))\n", targets[0].y);
  std::printf("      => 散射到顶点1 = -k * d_c.y = %+.10f\n", -k * targets[0].y);
  std::printf("      => 散射到顶点0 = +k * d_c.y = %+.10f  (随后被 pin 行覆盖)\n\n", k * targets[0].y);

  // 实际组装
  std::vector<Vec3> predVec = {mesh.positions[0] + ctx.config.gravity * (h * h), pred1};
  Eigen::VectorXd b;
  assembleInertialRhs(mesh, predVec, h, ctx.config.velocityDamping, b);
  const Scalar bAfterInertial = b[4];
  DistanceTerm::scatterInto(mesh, targets, b.data(), static_cast<std::size_t>(b.size()));
  const Scalar bAfterScatter = b[4];
  applyPinRhs(mesh, b);
  const Scalar bFinal = b[4];

  std::printf("  代码实际组装：\n");
  std::printf("      仅惯性项后   b[4] = %+.10f\n", bAfterInertial);
  std::printf("      加散射后     b[4] = %+.10f   (散射贡献 = %+.10f)\n", bAfterScatter,
              bAfterScatter - bAfterInertial);
  std::printf("      覆盖 pin 后  b[4] = %+.10f\n\n", bFinal);

  Eigen::SparseMatrix<Scalar> L;
  assembleLeftHandSide(mesh, h, ctx.config.velocityDamping, L);
  std::printf("  左端对角 L(4,4) = %+.10f\n", L.coeff(4, 4));
  std::printf("  => y1 = b[4] / L(4,4) = %+.10f\n\n", bFinal / L.coeff(4, 4));

  // 真正跑 stepOnce
  SimContext ctx2;
  ctx2.config = ctx.config;
  {
    Mesh& mm = ctx2.mesh;
    mm.positions = {Vec3{0.0, 0.0, 0.0}, Vec3{0.0, y0, 0.0}};
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
  std::printf("  stepOnce(1 次迭代) = %+.10f\n", ctx2.mesh.positions[1].y);
  return 0;
}
