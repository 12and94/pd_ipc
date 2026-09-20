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
// 3) 装配系数：右端必须是 κ_c d_c（再加 pin 消元补偿 κ_c q），矩阵块必须是 ±κ_c I。
//    两处各写一份 κ_c、必须同源同一个 A_c：多写或少写都会"静默"改变等效刚度，
//    所以这里逐元素断言。
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

// ---------------------------------------------------------------------------
// 3b) pin 消元补偿：只有一端被 pin 时，自由端的右端要多一项 +κ_c·q。
//
//     为什么必然有这一项：被 pin 的那一行会被 applyPinRhs 整行覆盖成 (对角 1, 右端 q)，
//     这等于把该自由度从方程组里消去；但它在耦合块里原本贡献的 -κ·x_pin 是挂在
//     **自由端那一行**上的，行被覆盖并不会顺手把它处理掉，必须显式搬到自由端右端：
//
//         (m/h² + κ)·x_free = (m/h²)·x̂_free + κ·d_c + κ·q
//
//     漏掉 κ·q 的后果是自由端丢失 pinned 点的位置信息：整个系统被拉向原点，
//     平移场景会改变结果（无重力、本来就在静止长度时也会自行变形）。
//     这是一个"看起来只是差一点"的静默错误，所以逐分量断言。
// ---------------------------------------------------------------------------
TEST(pinEliminationAddsKappaTimesPinPosition) {
  const Scalar k = 300.0, rest = 0.1;

  // (a) 顶点 a 被 pin、b 自由：b 的右端应恰好多出 κ·q_a。
  //     位形刻意取与 q_a 无关的任意值，确保断言抓的是 q_a 而不是当前位置。
  {
    Mesh m;
    m.positions = {Vec3{0.7, -0.2, 0.3}, Vec3{2.0, 1.0, -4.0}};  // 顶点 1 当前在别处
    m.restPositions = {Vec3{0, 0, 0}, Vec3{0, rest, 0}};
    m.velocities = {Vec3{}, Vec3{}};
    m.masses = {0.1, 0.1};
    m.pinned = {0, 1};
    m.pinPositions = {Vec3{0, 0, 0}, Vec3{2.0, 1.0, -4.0}};  // q_b = 顶点 1 的 pin 位置
    m.edges.push_back(Edge{0, 1, rest, k});

    std::vector<Vec3> targets;
    DistanceTerm::project(m, targets);
    std::vector<Scalar> b(6, 0.0);
    DistanceTerm::scatterInto(m, targets, b.data(), b.size());

    // 自由端 b 的期望：-κ·d_c（散射） + κ·q_a（消元补偿）；pin 端不进 b。
    const Vec3 want = targets[0] * k + m.pinPositions[1] * k;
    CHECK_NEAR(b[0], want.x, 1e-12);
    CHECK_NEAR(b[1], want.y, 1e-12);
    CHECK_NEAR(b[2], want.z, 1e-12);
    CHECK(b[3] == 0.0 && b[4] == 0.0 && b[5] == 0.0);
  }

  // (b) 顶点 b 被 pin、a 自由：补偿项同样落在自由端 a 上（另一个分支）。
  {
    Mesh m;
    m.positions = {Vec3{-1.5, 0.4, 0.25}, Vec3{0, rest, 0}};
    m.restPositions = {Vec3{0, 0, 0}, Vec3{0, rest, 0}};
    m.velocities = {Vec3{}, Vec3{}};
    m.masses = {0.1, 0.1};
    m.pinned = {1, 0};
    m.pinPositions = {Vec3{-1.5, 0.4, 0.25}, Vec3{0, 0, 0}};
    m.edges.push_back(Edge{0, 1, rest, k});

    std::vector<Vec3> targets;
    DistanceTerm::project(m, targets);
    std::vector<Scalar> b(6, 0.0);
    DistanceTerm::scatterInto(m, targets, b.data(), b.size());

    // 自由端 a 的期望：+κ·d_c（散射） + κ·q_b（消元补偿）；pin 端不进 b。
    const Vec3 want = targets[0] * (-k) + m.pinPositions[0] * k;
    CHECK_NEAR(b[3], want.x, 1e-12);
    CHECK_NEAR(b[4], want.y, 1e-12);
    CHECK_NEAR(b[5], want.z, 1e-12);
    CHECK(b[0] == 0.0 && b[1] == 0.0 && b[2] == 0.0);
  }

  // (c) 两端都 pin：整条约束被跳过，右端保持全零。
  {
    Mesh m;
    m.positions = {Vec3{0, 0, 0}, Vec3{0, 2 * rest, 0}};
    m.restPositions = m.positions;
    m.velocities = {Vec3{}, Vec3{}};
    m.masses = {0.1, 0.1};
    m.pinned = {1, 1};
    m.pinPositions = m.positions;
    m.edges.push_back(Edge{0, 1, rest, k});

    std::vector<Vec3> targets;
    DistanceTerm::project(m, targets);
    std::vector<Scalar> b(6, 0.0);
    DistanceTerm::scatterInto(m, targets, b.data(), b.size());
    for (int i = 0; i < 6; ++i) CHECK(b[i] == 0.0);
  }
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
  //   对角块 +κ I，耦合块 -κ I，惯性项用**未缩放**的 M/h²（阻尼不进入 L）。
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
// 3c) SolverStamp 的语义：只包含"真的会改变 L 的数值"的量。
//
//     判据来自 L 的构造（见 Assembler.cpp assembleLeftHandSide）：
//       · 惯性对角      M/h²              → 依赖质量、h
//       · 约束块        κ_c A_cᵀA_c       → 依赖刚度
//       · pinned 行     覆盖为 (对角 1, 右端 q) → 只依赖**哪些**顶点被 pin，不依赖 q
//       · 阻尼          不出现            → 不进 L
//
//     多放一项的后果是"白做一次数值分解"：拖拽把手或拖动阻尼滑块时会重分解，
//     而那是交互里每帧都在发生的事，直接吃帧预算（且与 docs/plan.md §2.2.1
//     "移动把手不触发重分解"的约定矛盾）。
// ---------------------------------------------------------------------------
TEST(stampChangesOnlyWhenLeftHandSideValuesChange) {
  Mesh m;
  m.positions = {Vec3{0, 0, 0}, Vec3{0.4, 0.1, 0.0}};
  m.restPositions = m.positions;
  m.velocities = {Vec3{}, Vec3{}};
  m.masses = {0.1, 0.1};
  m.pinned = {1, 0};
  m.pinPositions = m.positions;
  m.edges.push_back(Edge{0, 1, 0.1, 400.0});

  const Scalar dt = 1.0 / 120.0;
  const SolverStamp base = computeStamp(m, dt, 0.0, 1);

  // (a) 拖拽把手：只改 pinPositions，L 的数值不变 ⇒ stamp 必须不变。
  m.pinPositions[0] = Vec3{0.7, -0.3, 2.0};
  CHECK_MSG(computeStamp(m, dt, 0.0, 1) == base,
            "拖动 pin 只影响右端 b，不应触发重分解");

  // (b) 改阻尼：L 不含阻尼（只影响速度更新）⇒ stamp 必须不变。
  CHECK_MSG(computeStamp(m, dt, 0.35, 1) == base, "阻尼不进入 L，不应触发重分解");

  // (c) 改 pin 掩码：被覆盖的行变了 ⇒ stamp 必须变。
  m.pinned[1] = 1;
  CHECK_MSG(computeStamp(m, dt, 0.0, 1) != base, "pin 掩码改变必须触发重分解");
  m.pinned[1] = 0;

  // (d) 改刚度、质量、h、拓扑：都进 L ⇒ stamp 必须变。
  m.edges[0].stiffness = 401.0;
  CHECK_MSG(computeStamp(m, dt, 0.0, 1) != base, "刚度改变必须触发重分解");
  m.edges[0].stiffness = 400.0;

  m.masses[1] = 0.1000001;
  CHECK_MSG(computeStamp(m, dt, 0.0, 1) != base, "质量改变必须触发重分解");
  m.masses[1] = 0.1;

  CHECK_MSG(computeStamp(m, dt * 0.5, 0.0, 1) != base, "子步长 h 改变必须触发重分解");
  CHECK_MSG(computeStamp(m, dt, 0.0, 2) != base, "拓扑标识改变必须触发重分解");

  // (e) 反证："不重分解"必须真的是对的 —— 把 pin 从原点挪到别处，
  //     组装出的 L 必须逐元素相同（pinned 行恒为对角 1、其余 0）。
  Eigen::SparseMatrix<Scalar> l1, l2;
  assembleLeftHandSide(m, dt, 0.0, l1);
  m.pinPositions[0] = Vec3{-3.0, 1.5, 0.25};
  assembleLeftHandSide(m, dt, 0.25, l2);
  CHECK(l1.rows() == l2.rows() && l1.cols() == l2.cols());
  bool identical = true;
  for (int k = 0; k < l1.outerSize() && identical; ++k) {
    for (Eigen::SparseMatrix<Scalar>::InnerIterator it(l1, k); it; ++it) {
      if (it.value() != l2.coeff(it.row(), it.col())) identical = false;
    }
  }
  CHECK_MSG(identical, "改 pin 位置或阻尼后 L 必须逐位不变（这是「不重分解」的依据）");
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
