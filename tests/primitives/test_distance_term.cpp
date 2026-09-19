// tests/primitives/test_distance_term.cpp
// 距离约束的"定义级"不变量与装配系数断言。
//
// 为什么单独有这个文件（docs/plan.md 决策 D15 / 铁律 §10.2）：
//   投影返回向量的长度必须**恰好**等于静止长度 —— 这是投影的定义。
//   先把这类原语级断言钉死，再看残差、能量、垂度这些派生量；
//   否则一个量级错误会被一堆互相印证的派生指标掩盖。
#include "core/assemble/Assembler.h"
#include "core/energy/DistanceTerm.h"
#include "core/mesh/Mesh.h"
#include "tests/primitives/test_harness.h"

#include <cmath>

using namespace pd;

// ---------------------------------------------------------------------------
// 1) 投影的定义级不变量
// ---------------------------------------------------------------------------
TEST(projectionLengthIsExactlyRestLength) {
  const Scalar rest = 0.1;
  // 覆盖：静止、拉伸 10 倍、压缩到 0.1 倍、斜向、极短边。
  const std::vector<std::pair<Vec3, Vec3>> cases = {
      {{0, 0, 0}, {rest, 0, 0}},
      {{0, 0, 0}, {rest * 10.0, 0, 0}},
      {{0, 0, 0}, {rest * 0.1, 0, 0}},
      {{0, 0, 0}, {0.3, -0.4, 0.5}},
      {{1.0, 2.0, 3.0}, {1.0 + rest, 2.0, 3.0}},
      {{0, 0, 0}, {rest * 1e-8, 0, 0}},  // 未触发退化保护，仍应给出长度 rest
  };
  for (const auto& c : cases) {
    const Vec3 d = DistanceTerm::projectOne(c.first, c.second, rest);
    CHECK_NEAR(length(d), rest, 1e-15);
  }
}

TEST(projectionIsCollinearAndForward) {
  const Scalar rest = 0.25;
  const Vec3 a{0.1, -0.2, 0.3};
  const Vec3 b{-0.4, 0.5, -0.6};
  const Vec3 d = DistanceTerm::projectOne(a, b, rest);
  const Vec3 raw = a - b;
  const Vec3 rawUnit = normalized(raw);
  const Vec3 projUnit = normalized(d);
  // 方向必须与 (x_a - x_b) 一致（同向，不是反向）。
  CHECK_NEAR(dot(rawUnit, projUnit), 1.0, 1e-14);
}

TEST(projectionIsIdentityAtRestConfiguration) {
  const Scalar rest = 0.37;
  // 两点距离恰好等于静止长度时，投影结果必须与原向量**逐分量严格相等**。
  const Vec3 a{1.0, 0.0, 0.0};
  const Vec3 b{0.63, 0.0, 0.0};  // a - b = (0.37, 0, 0)，长度即 rest
  const Vec3 d = DistanceTerm::projectOne(a, b, rest);
  const Vec3 raw = a - b;
  CHECK(d.x == raw.x && d.y == raw.y && d.z == raw.z);
}

TEST(projectionHandlesDegenerateEdge) {
  const Vec3 p{1.0, 2.0, 3.0};
  const Vec3 d = DistanceTerm::projectOne(p, p, 0.1);
  // 退化边方向未定义：约定返回零向量，且长度断言不适用（由调用方跳过）。
  CHECK(d.x == 0.0 && d.y == 0.0 && d.z == 0.0);
}

// ---------------------------------------------------------------------------
// 2) 梯度与有限差分一致（用真实能量定义，与代码里的组装独立）
// ---------------------------------------------------------------------------
TEST(gradientMatchesFiniteDifference) {
  const Scalar k = 1234.0;
  const Scalar rest = 0.2;
  const Vec3 pa{0.1, 0.3, -0.2};
  const Vec3 pb{-0.4, 0.25, 0.35};

  auto energy = [&](const Vec3& a, const Vec3& b) {
    const Scalar r = length(a - b);
    const Scalar d = r - rest;
    return 0.5 * k * d * d;
  };

  const Scalar eps = 1e-7;
  for (int comp = 0; comp < 3; ++comp) {
    Vec3 ap = pa, am = pa;
    ap[comp] += eps;
    am[comp] -= eps;
    const Scalar fd = (energy(ap, pb) - energy(am, pb)) / (2.0 * eps);

    // 解析梯度：∂E/∂x_a = k(r-ℓ) * u，其中 u = (x_a-x_b)/r。
    const Vec3 diff = pa - pb;
    const Scalar r = length(diff);
    const Vec3 grad = diff * (k * (r - rest) / r);
    CHECK_NEAR(grad[comp], fd, 1e-6 * std::fabs(fd) + 1e-9);
  }
}

// ---------------------------------------------------------------------------
// 3) 装配系数：右端必须是 2κ d，矩阵块必须是 ±κ I
//    系数 2 来自 prox 的定义（见 docs/plan.md §2.2.1）；漏掉它会让等效刚度减半，
//    而且是一个"静默"错误 —— 所以这里逐元素断言。
// ---------------------------------------------------------------------------
TEST(scatterUsesStandardPdSign) {
  Mesh m;
  m.positions = {Vec3{0, 0, 0}, Vec3{0, 0.3, 0}};
  m.restPositions = m.positions;
  m.velocities = {Vec3{}, Vec3{}};
  m.masses = {0.1, 0.1};
  m.pinned = {0, 0};
  m.pinPositions = m.positions;
  const Scalar k = 500.0, rest = 0.1;
  m.edges.push_back(Edge{0, 1, rest, k});

  std::vector<Vec3> targets;
  DistanceTerm::project(m, targets);
  CHECK_NEAR(length(targets[0]), rest, 1e-15);

  // 标准 PD 的投影：d_c = ℓ · unit(x_a - x_b)。
  // 本案 a=0 在原点、b=1 在 +y ⇒ x_a - x_b = -0.3ŷ ⇒ d_c.y = -ℓ。
  CHECK_NEAR(targets[0].y, -rest, 1e-15);

  std::vector<Scalar> b(6, 0.0);
  DistanceTerm::scatterInto(m, targets, b.data(), b.size());

  // 标准 PD 的散射：b_a += κ d_c，b_b -= κ d_c。
  //   b_a(=顶点0).y = +κ·d_c.y = -κℓ
  //   b_b(=顶点1).y = -κ·d_c.y = +κℓ
  CHECK_NEAR(b[1], -k * rest, 1e-12);  // 顶点 0
  CHECK_NEAR(b[4], +k * rest, 1e-12);  // 顶点 1
  // 其它分量为 0。
  CHECK(b[0] == 0.0 && b[2] == 0.0 && b[3] == 0.0 && b[5] == 0.0);
}

TEST(assembledMatrixHasConstantBlocks) {
  Mesh m;
  m.positions = {Vec3{0, 0, 0}, Vec3{0.123, 0.3, -0.05}};  // 任意位形
  m.restPositions = {Vec3{0, 0, 0}, Vec3{0, 0.1, 0}};
  m.velocities = {Vec3{}, Vec3{}};
  m.masses = {0.1, 0.2};
  m.pinned = {0, 0};
  m.pinPositions = m.positions;
  const Scalar k = 777.0;
  m.edges.push_back(Edge{0, 1, 0.1, k});
  m.buildSparsityPattern();

  const Scalar dt = 1.0 / 100.0;
  Eigen::SparseMatrix<Scalar> L;
  assembleLeftHandSide(m, dt, 0.0, L);

  // 标准 PD：L = M/h² + Σ κ A_cᵀA_c
  //   对角块 +κ I，耦合块 -κ I，惯性项用未缩放的 M/h²（阻尼为 0 时 M_γ = M）。
  CHECK_NEAR(L.coeff(0, 0), m.masses[0] / (dt * dt) + k, 1e-9);
  CHECK_NEAR(L.coeff(4, 4), m.masses[1] / (dt * dt) + k, 1e-9);
  // 非对角块：-κ I
  CHECK_NEAR(L.coeff(4, 1), -k, 1e-9);
  CHECK_NEAR(L.coeff(1, 4), -k, 1e-9);
  // 矩阵必须对称
  CHECK_NEAR(L.coeff(1, 4) - L.coeff(4, 1), 0.0, 1e-12);
}

TEST(pinnedRowIsExactIdentityRow) {
  Mesh m;
  m.positions = {Vec3{1, 2, 3}, Vec3{1, 2, 3.3}};
  m.restPositions = m.positions;
  m.velocities = {Vec3{}, Vec3{}};
  m.masses = {0.1, 0.1};
  m.pinned = {1, 0};
  m.pinPositions = {Vec3{1, 2, 3}, Vec3{1, 2, 3.3}};
  m.edges.push_back(Edge{0, 1, 0.3, 400.0});
  m.buildSparsityPattern();

  const Scalar dt = 1.0 / 120.0;
  Eigen::SparseMatrix<Scalar> L;
  assembleLeftHandSide(m, dt, 0.0, L);

  // pinned 行：对角 1、其余为 0。
  CHECK_NEAR(L.coeff(0, 0), 1.0, 1e-15);
  CHECK_NEAR(L.coeff(0, 1), 0.0, 1e-15);
  CHECK_NEAR(L.coeff(0, 4), 0.0, 1e-15);
  CHECK_NEAR(L.coeff(1, 0), 0.0, 1e-15);

  // 自由顶点 1 的对角**必须包含**这条约束的 +κ：
  // 顶点 0 被 pin，但弹簧仍挂在顶点 1 上，这份刚度不能丢（丢了整个系统会明显偏软）。
  // 期望值 = m/h² + κ = 0.1·120² + 400 = 1440 + 400 = 1840。
  CHECK_NEAR(L.coeff(4, 4), m.masses[1] / (dt * dt) + 400.0, 1e-9);
  // 双向耦合项因为有一端被 pin 而消失（那一行被整行覆盖）。
  CHECK_NEAR(L.coeff(4, 1), 0.0, 1e-15);
  CHECK_NEAR(L.coeff(1, 4), 0.0, 1e-15);
}

// ---------------------------------------------------------------------------
// 4) 规则网格的约束集自检：4-邻域边数 = 2*nx*ny - nx - ny，且无重复边
// ---------------------------------------------------------------------------
TEST(gridEdgeCountAndNoDuplicates) {
  for (const auto& dims : {std::pair<int, int>{2, 2}, {5, 4}, {13, 7}, {40, 40}}) {
    Mesh m = Mesh::makeGrid(dims.first, dims.second, 0.05);
    CHECK(m.vertexCount() == dims.first * dims.second);
    const int expected = 2 * dims.first * dims.second - dims.first - dims.second;
    CHECK_MSG(m.edgeCount() == expected,
              "边数应为 " + std::to_string(expected) + "，实际 " + std::to_string(m.edgeCount()));
    CHECK(m.countDuplicateEdges() == 0);
  }
}

// ---------------------------------------------------------------------------
// 5) 位置打包/解包往返一致
// ---------------------------------------------------------------------------
TEST(packUnpackRoundTrip) {
  std::vector<Vec3> in = {{1.5, -2.25, 3.125}, {0, 0, 0}, {-1e-9, 1e9, 0.5}};
  Eigen::VectorXd packed;
  packPositions(in, packed);
  std::vector<Vec3> out;
  unpackPositions(packed, out);
  CHECK(out.size() == in.size());
  for (std::size_t i = 0; i < in.size(); ++i) {
    CHECK(out[i].x == in[i].x && out[i].y == in[i].y && out[i].z == in[i].z);
  }
}

TEST_MAIN("primitives/distance_term")
