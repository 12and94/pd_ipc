// _verify/bend_audit.cpp
// 线性（中点）弯曲约束的**独立诊断程序**（`pd_bendaudit`）。
//
// 为什么单独写它：弯曲约束没有局部步、没有散射，它的全部影响就是"左端块 + 右端常向量"，
// 所以一旦出问题，"物理不像"这种描述完全没有分辨力 —— 必须逐条把代数对上去：
//   A. 三条链 + 只开弯曲 + 只 pin 两端：平衡点必须是"三点共线"（A_c x = 0，能量的全局极小）。
//      —— 这一条同时验证「左端块」「右端常向量」「pin 行覆盖」三者的符号是否自洽：
//         任一处符号错，解都会落到"中间点被推离中点"的位置上。
//   B. 线性解必须让**线性**残差到机器精度（L x == b）：验证右端装配的每一项都进去了。
//   C. 弯曲残差公式（-k·w_v·A_c x）与**能量的有限差分梯度**一致：验证残差口径本身没错。
//   D. 网格上"平直静止位形 ⇒ 弯曲力恒为 0"（A_c x = 0 的平凡情形）。
//   E. 网格上"给中间点一个横向扰动 ⇒ 弯曲力方向指向直线"（回复力方向）。
//   F. `L = Ã ⊗ I₃`：跨 x/y/z 分量的块必须仍为 0（否则三分量并行会静默失效）。
//
// 注意：本程序不改生产代码，只读 core 的公开接口（与 `_verify/solve_audit.cpp` 同一风格）。
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

#include "core/assemble/Assembler.h"
#include "core/energy/DistanceTerm.h"
#include "core/mesh/Mesh.h"
#include "core/solver/IGlobalSolver.h"
#include "core/solver/EigenDirectSolver.h"

using namespace pd;

namespace {

int g_failed = 0;

void item(const char* name, bool ok, double got, double want, double tol) {
  std::printf("  %-46s %s   got=%.6e want=%.6e tol=%.1e\n", name, ok ? "OK" : "NG", got, want, tol);
  if (!ok) ++g_failed;
}

/// 手工搭一条 n 顶点的直链（只有弯曲 stencil，没有距离约束），pin 两端。
/// 顶点沿 x 等距排列，中间顶点沿 y 扰动 eps。
Mesh makeChain(int n, Scalar spacing, Scalar eps, Scalar kBend) {
  Mesh m;
  m.positions.resize(static_cast<std::size_t>(n));
  m.restPositions.resize(static_cast<std::size_t>(n));
  m.velocities.assign(static_cast<std::size_t>(n), Vec3{});
  m.masses.assign(static_cast<std::size_t>(n), 1.0);
  m.pinned.assign(static_cast<std::size_t>(n), 0);
  m.pinPositions.resize(static_cast<std::size_t>(n));
  for (int i = 0; i < n; ++i) {
    const Vec3 p{static_cast<Scalar>(i) * spacing, (i == n / 2) ? eps : Scalar{0}, Scalar{0}};
    m.positions[static_cast<std::size_t>(i)] = p;
    m.restPositions[static_cast<std::size_t>(i)] = p;
  }
  for (int i = 1; i + 1 < n; ++i) {
    BendStencil s;
    s.a = i - 1;
    s.b = i;
    s.c = i + 1;
    s.stiffness = kBend;
    m.bends.push_back(s);
  }
  m.pinned[0] = 1;
  m.pinned[static_cast<std::size_t>(n - 1)] = 1;
  m.pinPositions = m.positions;
  m.buildSparsityPattern();
  return m;
}

}  // namespace

int main() {
  std::printf("=== 线性（中点）弯曲约束审计 ===\n");

  // ---- 0：求解器的自检（对角正定矩阵，解必须是 b/d）----
  // 必须先排除"求解器/调用方式本身有问题"这种可能，否则后面的诊断全是噪声。
  {
    Eigen::SparseMatrix<Scalar> D(3, 3);
    D.insert(0, 0) = 2.0;
    D.insert(1, 1) = 4.0;
    D.insert(2, 2) = 8.0;
    D.makeCompressed();
    Eigen::VectorXd bb(3);
    bb << 2.0, 8.0, 24.0;
    EigenDirectSolver s;
    s.analyze(3, D);
    s.factorize(D);
    Eigen::VectorXd xx;
    s.solve(bb, xx);
    std::printf("\n[0] 求解器自检：diag(2,4,8)·x = (2,8,24) ⇒ x 应为 (1,2,3)，实得 (%.6g, %.6g, %.6g)\n",
                xx.size() > 2 ? static_cast<double>(xx[0]) : 0.0,
                xx.size() > 2 ? static_cast<double>(xx[1]) : 0.0,
                xx.size() > 2 ? static_cast<double>(xx[2]) : 0.0);
    const bool ok = xx.size() == 3 && std::fabs(xx[0] - 1.0) < 1e-12 && std::fabs(xx[1] - 2.0) < 1e-12 &&
                    std::fabs(xx[2] - 3.0) < 1e-12;
    item("0  求解器对角自检", ok, ok ? 0.0 : 1.0, 0.0, 0.0);
  }

  const Scalar kBend = 1000.0;
  const Scalar h = 1.0 / 120.0;

  // ---- A：三条链、只开弯曲、只 pin 两端。隐式欧拉一步的平衡点是解析可算的 ----
  {
    std::printf("\n[A] 三条链 + 只 pin 两端：与隐式欧拉解析解对照（$y=\\hat y\\,(m/h^2)/(m/h^2+4k)$）\n");
    const int n = 3;
    const Scalar spacing = 0.1;
    Mesh m = makeChain(n, spacing, 0.03, kBend);  // 中间顶点初始抬高 0.03

    Eigen::SparseMatrix<Scalar> L;
    std::vector<Scalar> bendRhs;
    assembleLeftHandSide(m, h, 0.0, L, &bendRhs);

    // 这一小节的矩阵只有 9×9，直接打出来人工核对最省事
    // （矩阵是手工推导的，看数字比看代码快 —— 本轮就是靠它定位到"pinned 列没清掉"的）。
    std::printf("   [数值核对] 组装出的 L（dim=9）：\n");
    for (int r = 0; r < 9; ++r) {
      std::printf("     row %d:", r);
      for (int c = 0; c < 9; ++c) {
        const Scalar v = L.coeff(r, c);
        if (v != 0.0) std::printf(" (%d)=%.10g", c, static_cast<double>(v));
      }
      std::printf("\n");
    }
    // 期望：v0/v2 是 pinned ⇒ 对角恰好 1、该行其余为 0（自由行里也不能留 pinned 列！）；
    //       v1 自由 ⇒ 对角 = m/h² + k·w_b² = 14400 + 4000 = 18400。
    std::printf("     期望：pinned 行 = (对角 1, 其余 0)；自由行 v1 对角 = m/h² + 4k = %.6g\n",
                static_cast<double>(m.masses[1] / (h * h) + 4.0 * kBend));
    bool pinRowsAreUnit = true;
    for (int v = 0; v < n; ++v) {
      if (!m.isPinned(v)) continue;
      for (int c = 0; c < 9; ++c) {
        for (int d = 0; d < 3; ++d) {
          // 只有**对角**（列号 == 行号）是 1，其余全为 0（自由行里也不能留 pinned 列）。
          const Scalar want = (c == 3 * v + d) ? 1.0 : 0.0;
          const Scalar got = L.coeff(3 * v + d, c);
          if (std::fabs(got) != want) {
            pinRowsAreUnit = false;
            std::printf("     [违例] L(%d,%d) = %.17g（期望 %.17g）\n", 3 * v + d, c,
                        static_cast<double>(got), static_cast<double>(want));
          }
        }
      }
    }
    item("A0 pinned 行必须是精确单位行（自由行里不得留 pinned 列）", pinRowsAreUnit, 0.0, 0.0, 0.0);

    bool decoupled = true;
    for (int r = 0; r < 9 && decoupled; ++r) {
      for (int c = 0; c < 9 && decoupled; ++c) {
        if (r % 3 != c % 3 && std::fabs(L.coeff(r, c)) > 0.0) decoupled = false;
      }
    }
    item("A1 跨 x/y/z 分量的块必须为 0（L = Ã⊗I₃）", decoupled, 0.0, 0.0, 0.0);

    // 实际求解：pin 行由 applyPinRhs 覆盖，自由行由 linear solve 得到
    std::vector<Vec3> predicted = m.positions;
    Eigen::VectorXd b;
    assembleInertialRhs(m, predicted, h, 0.0, b);
    for (int i = 0; i < 9; ++i) b[i] += bendRhs[static_cast<std::size_t>(i)];
    applyPinRhs(m, b);

    EigenDirectSolver solver;
    solver.analyze(9, L);
    solver.factorize(L);
    Eigen::VectorXd x;
    solver.solve(b, x);

    // ---- 这**不是**"一步就回到 y=0"：隐式欧拉一步的解是
    //        (m/h² + 4k)·y_b = (m/h²)·ŷ_b  ⇒  y_b = ŷ_b·(m/h²)/(m/h² + 4k)
    //      因为惯性项也在被极小化的目标里（速度不能一步归零）。
    //      这里 m/h² = 14400、4k = 4000 ⇒ y_b = 0.03·14400/18400 = 0.023478…。
    //      这个解析值就是正确性判据：它把弯曲力的大小、符号、以及惯性项的参与一起验了
    //     （本次实现过程中曾把 pinned 列留在自由行里，症状正是这个值完全不对 +
    //      线性残差 8.9e-2；详见 Assembler.cpp (2b) 的注记）。
    {
      const double want = 0.03 * (14400.0 / (14400.0 + 4.0 * static_cast<double>(kBend)));
      const double got = static_cast<double>(x[4]);
      std::printf("     解：中间顶点 y = %.9g（解析 = ŷ·(m/h²)/(m/h²+4k) = %.9g）\n", got, want);
      item("A2 一步隐式欧拉解 = ŷ·(m/h²)/(m/h²+4k)", std::fabs(got - want) < 1e-12, got, want, 1e-12);
    }

    std::printf("   [数值核对] 解的最大分量 = %.9g；线性残差 |Lx-b|∞ = %.3e\n",
                x.cwiseAbs().maxCoeff(), (L * x - b).cwiseAbs().maxCoeff());

    // B：线性残差必须到机器精度（这正是"右端装配对了"的判据）
    const double lr = (L * x - b).cwiseAbs().maxCoeff() / std::max(1.0, b.cwiseAbs().maxCoeff());
    item("B  线性残差 |Lx-b|∞/|b|∞", lr < 1e-12, lr, 0.0, 1e-12);

    // C：弯曲力与能量有限差分梯度一致（用非平直位形，保证场不为零）
    Mesh m2 = makeChain(n, spacing, 0.03, kBend);
    auto bendEnergy = [&](const std::vector<Vec3>& pos) {
      Scalar e = 0.0;
      for (const auto& s : m2.bends) {
        const Vec3 second = pos[static_cast<std::size_t>(s.a)] - pos[static_cast<std::size_t>(s.b)] * 2.0 +
                            pos[static_cast<std::size_t>(s.c)];
        e += 0.5 * s.stiffness * lengthSquared(second);
      }
      return e;
    };
    // 解析力：v 的力 = -k·w_v·A_c x
    auto analyticForce = [&](const std::vector<Vec3>& pos, int v) {
      Vec3 f{0, 0, 0};
      for (const auto& s : m2.bends) {
        int w = 0;
        if (s.a == v) w = 1;
        else if (s.b == v) w = -2;
        else if (s.c == v) w = 1;
        else continue;
        const Vec3 second = pos[static_cast<std::size_t>(s.a)] - pos[static_cast<std::size_t>(s.b)] * 2.0 +
                            pos[static_cast<std::size_t>(s.c)];
        f += second * (-s.stiffness * static_cast<Scalar>(w));
      }
      return f;
    };
    const Scalar eps = 1e-7;
    double worstFd = 0.0;
    for (int comp = 0; comp < 3; ++comp) {
      std::vector<Vec3> pp = m2.positions, pm = m2.positions;
      pp[1][comp] += eps;
      pm[1][comp] -= eps;
      const double fd = (bendEnergy(pp) - bendEnergy(pm)) / (2.0 * eps);
      const Vec3 f = analyticForce(m2.positions, 1);
      const double ana = f[comp];
      worstFd = std::max(worstFd, std::fabs(fd + ana));  // 力 = -dE/dx ⇒ dE/dx + f = 0
    }
    item("C  弯曲力 == -dE/dx（有限差分）", worstFd < 1e-6 * kBend * spacing,
         worstFd, 0.0, 1e-6 * kBend * spacing);
  }

  // ---- D/E：网格上的弯曲力（平直 ⇒ 0；横向扰动 ⇒ 指向直线）；F：跨分量必须仍为 0 ----
  {
    std::printf("\n[D/E/F] 20x20 规则网格上的弯曲力方向与结构\n");
    const int nx = 20, ny = 20;
    const Scalar spacing = 0.05;
    Mesh m = Mesh::makeGrid(nx, ny, spacing);
    const int added = m.addBendingStencils(nx, ny, kBend);
    // 行方向 (nx-2) 条 × ny 行 + 列方向 (ny-2) 条 × nx 列。20×20 ⇒ 2*18*20 = 720。
    const int expect = (nx - 2) * ny + (ny - 2) * nx;
    item("D0 stencil 计数 = (nx-2)*ny + (ny-2)*nx", added == expect, added, expect, 0.0);

    // 平直位形：A_c x = 0 ⇒ 力恒为 0（网格本身就在平面上）
    auto maxBendForce = [&](const Mesh& mm) {
      double worst = 0.0;
      for (int v = 0; v < mm.vertexCount(); ++v) {
        Vec3 f{0, 0, 0};
        const BendAdjacency& adj = mm.bendAdjacency();
        for (uint32_t k = adj.vertexStart[static_cast<std::size_t>(v)];
             k < adj.vertexStart[static_cast<std::size_t>(v) + 1]; ++k) {
          const BendStencil& s = mm.bends[adj.vertexStencils[k]];
          const Vec3& pa = mm.positions[static_cast<std::size_t>(s.a)];
          const Vec3& pb = mm.positions[static_cast<std::size_t>(s.b)];
          const Vec3& pc = mm.positions[static_cast<std::size_t>(s.c)];
          const Vec3 second = pa - pb * 2.0 + pc;
          f += second * (-s.stiffness * adj.vertexStencilWeight[k]);
        }
        worst = std::max(worst, static_cast<double>(length(f)));
      }
      return worst;
    };
    // 注：20×20 用 spacing=0.05，浮点上不是精确值 ⇒ 二阶差分只到 ~1e-13，不能断言"严格 0"。
    // 那条"平直位形严格为 0"的断言在 `tests/primitives` 里用二进制精确坐标（0.25 的整数倍）做，
    // 这里只做量级检查。
    item("D  平直静止位形下弯曲力 == 0（< 1e-10）", maxBendForce(m) < 1e-10, maxBendForce(m), 0.0,
         1e-10);

    // 给网格内部某个顶点一个 +y 扰动，检查它的回复力指向 -y（把点拉回平面）
    const int vp = 10 * nx + 10;
    m.positions[static_cast<std::size_t>(vp)].y += 0.01;
    {
      Vec3 f{0, 0, 0};
      const BendAdjacency& adj = m.bendAdjacency();
      for (uint32_t k = adj.vertexStart[static_cast<std::size_t>(vp)];
           k < adj.vertexStart[static_cast<std::size_t>(vp) + 1]; ++k) {
        const BendStencil& s = m.bends[adj.vertexStencils[k]];
        const Vec3& pa = m.positions[static_cast<std::size_t>(s.a)];
        const Vec3& pb = m.positions[static_cast<std::size_t>(s.b)];
        const Vec3& pc = m.positions[static_cast<std::size_t>(s.c)];
        const Vec3 second = pa - pb * 2.0 + pc;
        f += second * (-s.stiffness * adj.vertexStencilWeight[k]);
      }
      std::printf("     顶点 v%d 抬高 0.01 后的弯曲力 = (%.6g, %.6g, %.6g)\n", vp,
                  static_cast<double>(f.x), static_cast<double>(f.y), static_cast<double>(f.z));
      item("E  横向位移的回复力必须指向 -y", f.y < 0.0, f.y, -1.0, 0.0);
    }

    // F：装配出的 L 的跨分量块必须仍为 0
    Eigen::SparseMatrix<Scalar> L;
    std::vector<Scalar> bendRhs;
    assembleLeftHandSide(m, h, 0.0, L, &bendRhs);
    double worstCross = 0.0;
    for (int k = 0; k < L.outerSize(); ++k) {
      for (Eigen::SparseMatrix<Scalar>::InnerIterator it(L, k); it; ++it) {
        if (it.row() % 3 != it.col() % 3) worstCross = std::max(worstCross, std::fabs(it.value()));
      }
    }
    item("F  跨 x/y/z 分量的填充必须为 0", worstCross == 0.0, worstCross, 0.0, 0.0);

    std::printf("     L nnz = %lld；bendRhs 非零项 = ", (long long)L.nonZeros());
    int nz = 0;
    for (Scalar v : bendRhs) {
      if (v != 0.0) ++nz;
    }
    std::printf("%d / %zu\n", nz, bendRhs.size());
  }

  // ---- G：三条链 + 距离约束 + 两端 pin，一个子步后线性解必须满足 Lx = b ----
  {
    std::printf("\n[G] 距离约束 + 弯曲 混合：线性解必须满足 Lx = b（右端每项都进去了）\n");
    const int n = 3;
    Mesh m = makeChain(n, 0.1, 0.03, kBend);
    for (int i = 0; i + 1 < n; ++i) {
      Edge e;
      e.a = i;
      e.b = i + 1;
      e.restLength = 0.1;
      e.stiffness = 1.0e4;
      m.edges.push_back(e);
    }
    m.buildSparsityPattern();
    Eigen::SparseMatrix<Scalar> L;
    std::vector<Scalar> bendRhs;
    assembleLeftHandSide(m, h, 0.0, L, &bendRhs);
    std::vector<Vec3> predicted = m.positions;
    Eigen::VectorXd b;
    assembleInertialRhs(m, predicted, h, 0.0, b);
    // 距离约束的散射：b = bBase + Σκd_c + bendRhs（与生产同序）
    std::vector<Vec3> targets;
    DistanceTerm::project(m, targets);
    std::vector<Scalar> bv(b.data(), b.data() + b.size());
    DistanceTerm::scatterInto(m, targets, bv.data(), bv.size());
    for (std::size_t i = 0; i < bv.size(); ++i) b[i] = bv[i] + bendRhs[i];
    applyPinRhs(m, b);
    EigenDirectSolver solver;
    solver.analyze(static_cast<int>(b.size()), L);
    solver.factorize(L);
    Eigen::VectorXd x;
    solver.solve(b, x);
    const double lr = (L * x - b).cwiseAbs().maxCoeff() / std::max(1.0, b.cwiseAbs().maxCoeff());
    item("G  线性残差 |Lx-b|∞/|b|∞", lr < 1e-12, lr, 0.0, 1e-12);
  }

  // ---- H：尺度律：等效抗弯刚度对"网格间距 s"与"刚度 k"的幂次 ----
  //
  // 为什么需要它：文档里曾有 8 处写"等效抗弯刚度 ∝ s⁴/k"。那个写法**对 k 的关系是反的**
  // （k 越大当然越硬），幂次也可疑。这里用**确定性**的办法把它测准 —— 不跑迭代、不涉及时间：
  //
  //   ① 给一块**物理尺寸固定**的布（L×L）一个**给定的光滑形状** φ(x,y) = δ·(x/L)²；
  //      它的二阶导 φ_xx = 2δ/L² 是常数 ⇒ 二阶差分对它是**精确**的，没有离散化误差。
  //   ② 直接算弯曲能 E = (k/2)·Σ_c ‖A_c φ‖²（行、列两向 stencil 都算）。
  //   ③ 由连续板能量 E_cont = (D/2)·∫(φ_xx² + φ_yy²)dA 反解"等效抗弯刚度"
  //        D_eff = 2E / ∫(φ_xx² + φ_yy²)dA = 2E / (4δ²/L²)
  //   ④ 只改分辨率（N，s = L/(N−1)）或只改 k，看 D_eff 怎么变。
  //
  // 解析预期（同一套离散，可精确推导）：E = (k·s⁴/2)·(N−2)·N·φ_xx²
  //   ⇒ D_eff = k·s²·N(N−2)/(N−1)²  ≈  k·s²
  //   ⇒ **D_eff ∝ k·s²**（二维板，每单位面积）：k 越大越硬（线性）；**网格越细越软**，
  //      要保持同样的手感必须 **k ∝ 1/s²** —— 与文档原写的 `s⁴/k` 方向相反、幂次也不同。
  {
    std::printf("\n[H] 尺度律：等效抗弯刚度 D_eff 对 s / k 的幂次（给定形状直接算能量）\n");
    const double L = 1.0;
    const double delta = 0.1;
    const double phiXX = 2.0 * delta / (L * L);
    const double integral2D = phiXX * phiXX * L * L;   // ∫(φ_xx² + φ_yy²)dA（φ_yy = 0）

    // 二维板：给定 φ(x,y)=δ(x/L)²，返回 D_eff = 2E/∫(…)
    auto measureSheet = [&](int N, Scalar k) -> double {
      const double s = L / (N - 1);
      Mesh m = Mesh::makeGrid(N, N, static_cast<Scalar>(s));
      m.addBendingStencils(N, N, k);   // 行、列都加（Standard）
      for (int j = 0; j < N; ++j) {
        for (int i = 0; i < N; ++i) {
          const double x = i * s;
          const double phi = delta * (x / L) * (x / L);
          m.positions[static_cast<std::size_t>(j * N + i)] = Vec3{x, static_cast<Scalar>(phi), 0.0};
        }
      }
      double energy = 0.0;
      for (const BendStencil& st : m.bends) {
        const Vec3 second = m.positions[static_cast<std::size_t>(st.a)] -
                            m.positions[static_cast<std::size_t>(st.b)] * 2.0 +
                            m.positions[static_cast<std::size_t>(st.c)];
        const double sq = static_cast<double>(second.x) * second.x +
                          static_cast<double>(second.y) * second.y +
                          static_cast<double>(second.z) * second.z;
        energy += 0.5 * static_cast<double>(st.stiffness) * sq;
      }
      return 2.0 * energy / integral2D;
    };

    // 一维纤维（单行、手工搭链）：返回 EI_eff = 2E/∫φ''²dx
    // ⚠️ 这里**不能**用 `Mesh::makeGrid(N, 1, s)`：`makeGrid` 在 `ny < 2` 时直接返回**空网格**
    //    （Mesh.cpp:14），接着按 N 写 positions 就是越界写（实测 0xC0000005）。
    auto measureFiber = [&](int N, Scalar k) -> double {
      const double s = L / (N - 1);
      Mesh m;
      m.positions.resize(static_cast<std::size_t>(N));
      m.restPositions.resize(static_cast<std::size_t>(N));
      m.velocities.assign(static_cast<std::size_t>(N), Vec3{});
      m.masses.assign(static_cast<std::size_t>(N), 1.0);
      m.pinned.assign(static_cast<std::size_t>(N), 0);
      m.pinPositions.resize(static_cast<std::size_t>(N));
      for (int i = 0; i < N; ++i) {
        const double x = i * s;
        const double phi = delta * (x / L) * (x / L);
        m.positions[static_cast<std::size_t>(i)] = Vec3{x, static_cast<Scalar>(phi), 0.0};
        m.restPositions[static_cast<std::size_t>(i)] = Vec3{x, 0.0, 0.0};
      }
      for (int i = 1; i + 1 < N; ++i) {
        BendStencil st;
        st.a = i - 1;
        st.b = i;
        st.c = i + 1;
        st.stiffness = k;
        m.bends.push_back(st);
      }
      m.buildSparsityPattern();
      double energy = 0.0;
      for (const BendStencil& st : m.bends) {
        const Vec3 second = m.positions[static_cast<std::size_t>(st.a)] -
                            m.positions[static_cast<std::size_t>(st.b)] * 2.0 +
                            m.positions[static_cast<std::size_t>(st.c)];
        const double sq = static_cast<double>(second.x) * second.x +
                          static_cast<double>(second.y) * second.y +
                          static_cast<double>(second.z) * second.z;
        energy += 0.5 * static_cast<double>(st.stiffness) * sq;
      }
      return 2.0 * energy / (phiXX * phiXX * L);
    };

    const Scalar k0 = 1000.0;
    std::printf("  %-6s %-11s %-15s %-20s %-10s\n", "N", "s", "D_eff(二维板)", "D_eff/(k·s²)", "拟合幂次");
    double prev = 0.0, prevS = 1.0;
    const int ns[] = {8, 16, 32, 64};
    for (int N : ns) {
      const double s = L / (N - 1);
      const double d = measureSheet(N, k0);
      const double ratio = d / (static_cast<double>(k0) * s * s);
      char expo[32] = "-";
      if (prev > 0.0) {
        const double p = std::log(d / prev) / std::log(s / prevS);
        std::snprintf(expo, sizeof(expo), "%.4f", p);
      }
      std::printf("  %-6d %-11.6f %-15.6e %-20.6f %-10s\n", N, s, d, ratio, expo);
      prev = d;
      prevS = s;
    }
    // 断言 1：与解析式 D_eff = k·s²·N(N−2)/(N−1)² 逐位吻合（抛物线的二阶差分是精确的）
    for (int N : ns) {
      const double s = L / (N - 1);
      const double want = static_cast<double>(k0) * s * s * N * (N - 2) / ((N - 1.0) * (N - 1.0));
      char name[80];
      std::snprintf(name, sizeof(name), "H  二维板 N=%d 与解析式一致", N);
      const double got = measureSheet(N, k0);
      item(name, std::fabs(got - want) <= 1e-12 * std::fabs(want), got, want, 1e-12);
    }
    // 断言 2：对 k 线性（k 加倍 ⇒ D_eff 加倍）
    {
      const double a = measureSheet(16, k0), b2 = measureSheet(16, 2.0 * k0);
      item("H  对 k 线性（2k ⇒ 2×D_eff）", std::fabs(b2 / a - 2.0) < 1e-12, b2 / a, 2.0, 1e-12);
    }
    // 断言 3：拟合幂次的**渐近值**（二维 2 / 一维 3）。容差放到 0.05：
    //   有限 N 下拟合值带一个 O(1/N) 偏置 —— 例如一维的精确式是 EI = k·s³·(N−2)/(N−1)，
    //   在 N=16→64 上拟合出 2.963 而不是 3.000（解析可算，见断言 4），这不是误差。
    {
      const double s1 = L / 15.0, s2 = L / 63.0;
      const double p2 = std::log(measureSheet(64, k0) / measureSheet(16, k0)) / std::log(s2 / s1);
      const double p1 = std::log(measureFiber(64, k0) / measureFiber(16, k0)) / std::log(s2 / s1);
      item("H  二维板拟合幂次 ≈ 2（渐近）", std::fabs(p2 - 2.0) < 0.05, p2, 2.0, 0.05);
      item("H  一维纤维拟合幂次 ≈ 3（渐近）", std::fabs(p1 - 3.0) < 0.05, p1, 3.0, 0.05);
    }
    // 断言 4：与**精确离散公式**逐位吻合（比"拟合幂次"强得多）：
    //   二维板 D_eff = k·s²·N(N−2)/(N−1)²；一维纤维 EI_eff = k·s³·(N−2)/(N−1)。
    for (int N : ns) {
      const double s = L / (N - 1);
      const double wantF = static_cast<double>(k0) * s * s * s * (N - 2) / (N - 1.0);
      char name[80];
      std::snprintf(name, sizeof(name), "H  一维纤维 N=%d 与解析式一致", N);
      const double gotF = measureFiber(N, k0);
      item(name, std::fabs(gotF - wantF) <= 1e-12 * std::fabs(wantF), gotF, wantF, 1e-12);
    }
    std::printf("  ⇒ 结论：等效抗弯刚度 ∝ **k·s²**（二维板，每单位面积）—— 对 k **线性**（不是倒数！）；\n"
                "     一维纤维 ∝ **k·s³**。保持同样手感需要 **k ∝ 1/s²**（二维板）：\n"
                "     网格越细，同一个 k 越软 ⇒ 必须**调大** k。文档原先写的「∝ s⁴/k」方向与幂次都不对。\n");
  }

  std::printf("\n  %s（失败 %d 项）\n", g_failed == 0 ? "全部正确" : "存在错误", g_failed);
  return g_failed == 0 ? 0 : 1;
}
