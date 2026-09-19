// _verify/standard_pd.cpp
// 目的：用标准 PD 公式**独立**组装并求解一次两步系统，与代码组装的结果逐项对照。
// 不做任何补丁式修改，只做对照。
//
// 标准 PD（Bouaziz et al. 2014 / ADMM 形式）：
//   E = 1/2‖x-x̂‖²_{M/h²} + Σ_c (w_c/2)‖A_c x - d_c‖²,  d_c = Π_c(A_c x) = ℓ·unit(A_c x)
//   (M/h² + Σ_c w_c A_cᵀA_c) x = (M/h²) x̂ + Σ_c w_c A_cᵀ d_c
//   A_c x = x_a - x_b ⇒ L_aa += w_c I, L_bb += w_c I, L_ab -= w_c I, L_ba -= w_c I
//                        b_a += w_c d_c, b_b -= w_c d_c
#include <cstdio>

#include <Eigen/Dense>

#include "core/assemble/Assembler.h"
#include "core/energy/DistanceTerm.h"
#include "core/sim/Integrator.h"
#include "core/sim/Scene.h"
#include "core/solver/EigenDirectSolver.h"

using namespace pd;

int main() {
  const Scalar h = 1.0 / 120.0, m = 0.05, ell = 0.1, k = 1.0e4, g = 9.81;
  const Scalar y1 = 0.12;  // 自由顶点初始高度（顶点 0 pinned 在原点）

  SimContext ctx;
  ctx.config.dt = h;
  ctx.config.gravity = Vec3{0.0, -g, 0.0};
  ctx.config.stiffness = k;
  Mesh& mesh = ctx.mesh;
  mesh.positions = {Vec3{0.0, 0.0, 0.0}, Vec3{0.0, y1, 0.0}};
  mesh.restPositions = mesh.positions;
  mesh.velocities = {Vec3{}, Vec3{}};
  mesh.masses = {m, m};
  mesh.pinned = {1, 0};
  mesh.pinPositions = mesh.positions;
  mesh.edges.push_back(Edge{0, 1, ell, k});
  mesh.buildSparsityPattern();
  ensureBuffers(ctx);

  const Scalar mh2 = m / (h * h);
  const Scalar yHat = y1 - h * h * g;  // 预测（初速 0）

  // ---- 标准 PD：独立组装并求解（只取 y 分量，3 自由度可解耦） ----
  // 顶点 0 = pinned（未知量被固定为 0），顶点 1 = 自由。
  const Scalar Acx = 0.0 - y1;                 // A_c x = x_a - x_b，a=0, b=1
  const Scalar unitY = Acx / std::fabs(Acx);   // = -1
  const Scalar dcY = unitY * ell;              // d_c 的 y 分量 = -0.1

  // 代码里的 pinned 行：L_00 = 1, b_0 = 0；自由行：L_11 = m/h² + κ, b_1 = (m/h²)ŷ + (-κ d_c)
  const Scalar L11 = mh2 + k;
  const Scalar b1 = mh2 * yHat - k * dcY;
  const Scalar ySolve = b1 / L11;
  std::printf("标准 PD 独立求解（一步）:\n");
  std::printf("  A_c x = %.4g  ⇒ unit = %.4g  ⇒ d_c.y = %.4g\n", Acx, unitY, dcY);
  std::printf("  L11 = m/h² + κ = %.6g\n", L11);
  std::printf("  b1  = (m/h²)ŷ - κ d_c.y = %.6g * %.6g - %.6g * %.6g = %.6g\n", mh2, yHat, k, dcY, b1);
  std::printf("  y   = %.10g\n", ySolve);

  // ---- 代码组装 ----
  std::vector<Vec3> targets;
  DistanceTerm::project(mesh, targets);
  Eigen::SparseMatrix<Scalar> L;
  assembleLeftHandSide(mesh, h, 0.0, L);
  std::vector<Vec3> pred = mesh.positions;
  pred[1].y = yHat;
  Eigen::VectorXd b;
  assembleInertialRhs(mesh, pred, h, 0.0, b);
  const Scalar bInertial = b[4];
  DistanceTerm::scatterInto(mesh, targets, b.data(), static_cast<std::size_t>(b.size()));
  applyPinRhs(mesh, b);

  std::printf("\n代码组装:\n");
  std::printf("  d_c.y = %.4g   L(4,4) = %.6g   L(4,1) = %.6g\n", targets[0].y, L.coeff(4, 4),
              L.coeff(4, 1));
  std::printf("  惯性 b[4] = %.6g   散射 = %.6g   最终 b[4] = %.6g\n", bInertial, b[4] - bInertial,
              b[4]);

  // ---- 用稠密解核验代码的 L 与 b ----
  Eigen::MatrixXd Ld = Eigen::MatrixXd(L);
  Eigen::VectorXd bd = b;
  const Eigen::VectorXd x = Ld.fullPivLu().solve(bd);
  std::printf("\n代码 L·x = b 的稠密解: y = %.10g\n", x[4]);
  std::printf("标准 PD 独立解:        y = %.10g\n", ySolve);
  std::printf("两者之差 = %.3g  %s\n", std::fabs(x[4] - ySolve),
              std::fabs(x[4] - ySolve) < 1e-12 ? "【一致】" : "【不一致 → 组装偏离标准 PD】");

  // ---- 多次迭代后的不动点（代码） ----
  SimContext ctx2;
  ctx2.config = ctx.config;
  ctx2.config.maxIterations = 60;
  {
    Mesh& mm = ctx2.mesh;
    mm.positions = {Vec3{0.0, 0.0, 0.0}, Vec3{0.0, y1, 0.0}};
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
  std::printf("\n代码 stepOnce（60 次迭代）y = %.10g\n", ctx2.mesh.positions[1].y);
  std::printf("（注意：y 变了之后 d_c 也变，故它与上面 一步固定 d_c 的解不同是正常的）\n");
  return 0;
}
