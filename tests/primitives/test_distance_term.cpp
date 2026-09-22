// tests/primitives/test_distance_term.cpp
// 距离约束的"定义级"不变量与装配系数断言。
//
// 为什么单独有这个文件（docs/plan.md 决策 D15 / 铁律 §10.2）：
//   投影返回向量的长度必须**恰好**等于静止长度 —— 这是投影的定义。
//   先把这类原语级断言钉死，再看残差、能量、垂度这些派生量；
//   否则一个量级错误会被一堆互相印证的派生指标掩盖。
#include "core/assemble/Assembler.h"
#include "core/energy/DistanceTerm.h"
#include "core/math/Parallel.h"
#include "core/mesh/Mesh.h"
#include "core/solver/EigenDirectSolver.h"
#include "tests/primitives/test_harness.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <memory>
#include <string>

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

// ---------------------------------------------------------------------------
// 6) 并行散射分支的对照 —— 此前**完全无覆盖**，这是它的第一个测试
//
// 为什么必须单独钉它：`scatterInto` 在 `ne < 256` 或单线程时走**串行分支**，
// 而本文件此前所有散射测试（以及 _verify/direction_audit）用的都是**单条边**的网格，
// 于是 **OpenMP 分支从未被任何测试执行过** —— 只有 bench/viewer 会走到它。
// 这正是"并行散射每线程持有一份全维缓冲、每次调用清零并全量累加 P*dim"
// 这个 O(P*dim) 代价长期没被发现的原因（见 README §5 第 7 项）。
//
// 判据：并行分支与串行分支、以及一份**独立写出的**参考实现，必须给出同一个右端。
// ---------------------------------------------------------------------------
namespace {

/// 20x20 规则网格：400 顶点、760 条边（> 256，保证触发并行分支）。
/// 位形做成确定性的非平凡扰动：既避开退化边，也避开对称（对称会让"符号错"
/// 这类缺陷在不同位置互相抵消，掩盖问题）。
Mesh makePerturbedGrid(int n, Scalar spacing, Scalar stiffness) {
  Mesh m = Mesh::makeGrid(n, n, spacing);
  for (auto& e : m.edges) e.stiffness = stiffness;
  for (int v = 0; v < m.vertexCount(); ++v) {
    const Scalar i = static_cast<Scalar>(v % n);
    const Scalar j = static_cast<Scalar>(v / n);
    Vec3& p = m.positions[static_cast<std::size_t>(v)];
    p.x += spacing * 0.20 * std::sin(1.7 * i + 0.5 * j);
    p.y += spacing * 0.30 * std::cos(0.9 * i - 1.3 * j);
    p.z += spacing * 0.25 * std::sin(0.4 * i + 2.1 * j);
  }
  // 钉两个对角顶点：边表里同时存在"两端自由""a 端 pin""b 端 pin"三种情形，
  // 于是 pin 消元补偿的两个分支都会被走到。
  m.pinned[0] = 1;
  m.pinned[static_cast<std::size_t>(m.vertexCount() - 1)] = 1;
  m.pinPositions = m.positions;  // 契约：有 pinned 就必须填好 pinPositions
  return m;
}

/// 独立参考实现：按**边序号顺序**把贡献加进 b。刻意与 core 的实现分开写，
/// 避免"自己验证自己"。语义取自 docs/plan.md §2.2.1：
///   b_a += κ d_c、b_b -= κ d_c；单端 pinned 时给自由端补 +κ q。
std::vector<Scalar> referenceScatter(const Mesh& m, const std::vector<Vec3>& targets) {
  std::vector<Scalar> b(static_cast<std::size_t>(3 * m.vertexCount()), 0.0);
  for (std::size_t c = 0; c < m.edges.size(); ++c) {
    const Edge& e = m.edges[c];
    const bool aPinned = m.isPinned(e.a);
    const bool bPinned = m.isPinned(e.b);
    if (aPinned && bPinned) continue;
    const Vec3 contrib = targets[c] * e.stiffness;
    const std::size_t ia = static_cast<std::size_t>(e.a) * 3;
    const std::size_t ib = static_cast<std::size_t>(e.b) * 3;
    if (!aPinned) {
      b[ia + 0] += contrib.x;
      b[ia + 1] += contrib.y;
      b[ia + 2] += contrib.z;
    }
    if (!bPinned) {
      b[ib + 0] -= contrib.x;
      b[ib + 1] -= contrib.y;
      b[ib + 2] -= contrib.z;
    }
    if (aPinned && !bPinned) {
      const Vec3& q = m.pinPositions[static_cast<std::size_t>(e.a)];
      b[ib + 0] += e.stiffness * q.x;
      b[ib + 1] += e.stiffness * q.y;
      b[ib + 2] += e.stiffness * q.z;
    } else if (bPinned && !aPinned) {
      const Vec3& q = m.pinPositions[static_cast<std::size_t>(e.b)];
      b[ia + 0] += e.stiffness * q.x;
      b[ia + 1] += e.stiffness * q.y;
      b[ia + 2] += e.stiffness * q.z;
    }
  }
  return b;
}

double maxAbsValue(const std::vector<Scalar>& v) {
  double w = 0.0;
  for (Scalar x : v) w = std::max(w, std::fabs(static_cast<double>(x)));
  return w;
}

double maxAbsDiff(const std::vector<Scalar>& a, const std::vector<Scalar>& b) {
  double w = 0.0;
  for (std::size_t i = 0; i < a.size(); ++i) {
    w = std::max(w, std::fabs(static_cast<double>(a[i]) - static_cast<double>(b[i])));
  }
  return w;
}

}  // namespace

TEST(parallelScatterMatchesSerialReference) {
  const int savedThreads = numThreads();  // 必须还原：本套件其余测试共用同一进程
  Mesh m = makePerturbedGrid(20, 0.02, 2.0e3);

  // 前置条件：网格必须真的能触发并行分支（ne >= 256），否则这个测试等于没测。
  CHECK_MSG(m.edgeCount() >= 256,
            "边数 " + std::to_string(m.edgeCount()) + " < 256，不会走并行分支");
  CHECK(m.countDuplicateEdges() == 0);

  std::vector<Vec3> targets;
  DistanceTerm::project(m, targets);
  const std::vector<Scalar> ref = referenceScatter(m, targets);
  const std::size_t dim = static_cast<std::size_t>(3 * m.vertexCount());

  std::vector<Scalar> bSerial(dim, 0.0);
  setNumThreads(1);
  DistanceTerm::scatterInto(m, targets, bSerial.data(), bSerial.size());

  std::vector<Scalar> bParallel(dim, 0.0);
  setNumThreads(18);
  DistanceTerm::scatterInto(m, targets, bParallel.data(), bParallel.size());
  setNumThreads(savedThreads);

  const double scale = std::max(1.0, maxAbsValue(ref));
  const double tol = 1e-12 * scale;
  const double dSerialRef = maxAbsDiff(bSerial, ref);
  const double dParallelRef = maxAbsDiff(bParallel, ref);
  const double dParallelSerial = maxAbsDiff(bParallel, bSerial);

  char buf[360];
  std::snprintf(buf, sizeof(buf),
                "20x20 网格(%d 边, 含 2 个 pin) 最大差: 串行vs参考 %.3g, 并行vs参考 %.3g, "
                "并行vs串行 %.3g   (参考容差 %.3g, 量级 %.3g)",
                m.edgeCount(), dSerialRef, dParallelRef, dParallelSerial, tol, scale);
  CHECK_MSG(dSerialRef <= tol, buf);   // 参考实现按边序号序累加；顺序不同，只承诺 ≤1e-12
  CHECK_MSG(dParallelRef <= tol, buf);
  // 并行与 1 线程现在**必须位级相同**：两者都按颜色序累加（Phase 2 之后全项目只有这一条顺序）。
  CHECK_MSG(dParallelSerial == 0.0, "并行与 1 线程必须位级相同（同为颜色序累加）");

  std::printf("    [注] 并行 vs 1 线程位级相同: %s；vs 边序参考最大差 %.3g（顺序不同，只承诺 ≤1e-12）\n",
              (dParallelSerial == 0.0) ? "是" : "否", dParallelRef);
}

TEST(scatterIsThreadCountIndependent) {
  const int savedThreads = numThreads();
  Mesh m = makePerturbedGrid(20, 0.02, 2.0e3);
  std::vector<Vec3> targets;
  DistanceTerm::project(m, targets);
  const std::size_t dim = static_cast<std::size_t>(3 * m.vertexCount());

  std::vector<Scalar> base(dim, 0.0);
  setNumThreads(1);
  DistanceTerm::scatterInto(m, targets, base.data(), base.size());

  const double scale = std::max(1.0, maxAbsValue(base));
  double worst = 0.0;
  int runs = 0;
  int bitwiseSame = 0;

  for (const int threads : {2, 4, 8, 18}) {
    setNumThreads(threads);
    // 同一线程数重复 4 次：`schedule(dynamic, 64)` 的分块分配会随运行变化，
    // 因此"重复运行也位级相同"是一条有分辨力的断言（Phase 2 之前它不成立）。
    for (int rep = 0; rep < 4; ++rep) {
      std::vector<Scalar> b(dim, 0.0);
      DistanceTerm::scatterInto(m, targets, b.data(), b.size());
      const double d = maxAbsDiff(b, base);
      worst = std::max(worst, d);
      ++runs;
      if (d == 0.0) ++bitwiseSame;
      char buf[288];
      std::snprintf(buf, sizeof(buf),
                    "%d 线程第 %d 次: 与 1 线程结果的最大差 %.3g —— Phase 2 之后"
                    "累加顺序固定为颜色序，必须**位级相同**",
                    threads, rep + 1, d);
      CHECK_MSG(d == 0.0, buf);
    }
  }
  setNumThreads(savedThreads);

  std::printf("    [注] 跨线程数/重复运行共 %d 次，最大差 %.3g（量级 %.3g）；位级相同 %d 次\n",
              runs, worst, scale, bitwiseSame);
  CHECK_MSG(bitwiseSame == runs, "全部运行都必须与 1 线程结果位级相同");
}

// ---------------------------------------------------------------------------
// 7) 约束着色（Phase 2）：同色边两两不共享顶点 —— 散射因此不需要任何归约
//
// 这是"按颜色分组执行散射"的全部依据：同色内并发写 b 无冲突、无需私有缓冲、
// 无需原子加；而"同一顶点每色最多被写一次"又让求和顺序固定为颜色序 ⇒ 位级可复现。
// ---------------------------------------------------------------------------
TEST(constraintColoringIsProperAndComplete) {
  for (const auto& dims : {std::pair<int, int>{2, 2}, {5, 5}, {20, 20}, {40, 40}}) {
    Mesh m = Mesh::makeGrid(dims.first, dims.second, 0.05);
    const ConstraintColoring& col = m.constraintColoring();
    const int ne = m.edgeCount();
    const std::string tag = std::to_string(dims.first) + "x" + std::to_string(dims.second);

    CHECK_MSG(col.built, tag + "：着色应当已构建（buildSparsityPattern 会一并建好）");
    CHECK_MSG(static_cast<int>(col.colorOf.size()) == ne, tag + "：colorOf 必须覆盖每条边");
    CHECK_MSG(col.colorCount >= 1, tag + "：至少要有一色");

    bool complete = true;
    for (uint32_t c : col.colorOf) complete = complete && (c < static_cast<uint32_t>(col.colorCount));
    CHECK_MSG(complete, tag + "：每条边都必须被着色");

    // 分组（CSR）必须与 colorOf 自洽
    bool grouped = (static_cast<int>(col.order.size()) == ne) &&
                   (static_cast<int>(col.start.size()) == col.colorCount + 1) && col.begin(0) == 0 &&
                   static_cast<int>(col.end(col.colorCount - 1)) == ne;
    for (int c = 0; c < col.colorCount && grouped; ++c) {
      for (uint32_t k = col.begin(c); k < col.end(c); ++k) {
        grouped = grouped && (col.colorOf[col.order[k]] == static_cast<uint32_t>(c));
      }
    }
    CHECK_MSG(grouped, tag + "：按颜色分组的 CSR 必须与 colorOf 自洽");

    // 恰当性：同一颜色里不允许两条边共享顶点
    std::vector<int> owner(static_cast<std::size_t>(m.vertexCount()), -1);
    int conflicts = 0;
    for (int c = 0; c < col.colorCount; ++c) {
      std::fill(owner.begin(), owner.end(), -1);
      for (uint32_t k = col.begin(c); k < col.end(c); ++k) {
        const Edge& e = m.edges[col.order[k]];
        const std::size_t a = static_cast<std::size_t>(e.a);
        const std::size_t b = static_cast<std::size_t>(e.b);
        if (owner[a] >= 0 || owner[b] >= 0) ++conflicts;
        owner[a] = 1;
        owner[b] = 1;
      }
    }
    CHECK_MSG(conflicts == 0,
              tag + "：同色边不允许共享顶点（实测冲突 " + std::to_string(conflicts) + " 处）");

    // 下界：任何恰当边着色至少要用 Δ（最大顶点度）种颜色。
    // 只断言**下界**：MIS 迭代法的颜色数上界我援引不出可靠证明，写下来就是假门槛
    //（实测值用 [注] 打出来，见 docs/perf.md）。
    std::vector<int> deg(static_cast<std::size_t>(m.vertexCount()), 0);
    for (const Edge& e : m.edges) {
      ++deg[static_cast<std::size_t>(e.a)];
      ++deg[static_cast<std::size_t>(e.b)];
    }
    int maxDeg = 0;
    for (int d : deg) maxDeg = std::max(maxDeg, d);
    CHECK_MSG(col.colorCount >= maxDeg, tag + "：颜色数 " + std::to_string(col.colorCount) +
                                            " 小于最大顶点度 " + std::to_string(maxDeg));
    std::printf("    [注] %s：边 %d、颜色 %d（Δ=%d）\n", tag.c_str(), ne, col.colorCount, maxDeg);
  }
}

TEST(constraintColoringIsDeterministicAndTopologyOnly) {
  const int savedThreads = numThreads();
  Mesh m = makePerturbedGrid(20, 0.02, 2.0e3);  // 760 条边 ⇒ 走并行标记路径
  const ConstraintColoring reference = m.constraintColoring();

  // (a) 线程数不影响：1 线程与 18 线程各自重算，结果必须逐位相同
  ConstraintColoring c1;
  setNumThreads(1);
  buildConstraintColoring(m, c1);
  ConstraintColoring c18;
  setNumThreads(18);
  buildConstraintColoring(m, c18);
  setNumThreads(savedThreads);

  CHECK_MSG(c1.colorCount == c18.colorCount, "两次着色的颜色数必须相同");
  CHECK_MSG(c1.colorOf == c18.colorOf, "1 线程与 18 线程的 colorOf 必须逐位相同");
  CHECK_MSG(c1.order == c18.order, "1 线程与 18 线程的颜色分组必须逐位相同");
  CHECK_MSG(reference.colorOf == c1.colorOf && reference.order == c1.order,
            "buildSparsityPattern 建的着色与直接调用 buildConstraintColoring 必须一致");

  // (b) 位形 / 刚度 / pin 掩码都不影响着色（它是纯拓扑量，与 SolverStamp 关心
  //     的"是不是要重分解"无关，所以改 pin 不该触发重着色）
  for (auto& p : m.positions) p = p + Vec3{1.5, -2.5, 3.5};
  for (auto& e : m.edges) e.stiffness *= 3.0;
  for (std::size_t v = 0; v < m.pinned.size(); v += 7) m.pinned[v] = 1;
  ConstraintColoring after;
  buildConstraintColoring(m, after);
  CHECK_MSG(after.colorOf == reference.colorOf, "位形/刚度/pin 变了，着色不应改变");
  CHECK_MSG(after.order == reference.order, "位形/刚度/pin 变了，分组不应改变");
  CHECK_MSG(after.colorCount == reference.colorCount, "位形/刚度/pin 变了，颜色数不应改变");
}

// ---------------------------------------------------------------------------
// 8) 全局步按 x/y/z 三分量拆分（L = Ã ⊗ I₃，2026-09-22）
//
// 为什么必须钉住这三条：这条改造把全局步从"一次 3n 三角求解"变成"3 条互不相干的链"，
// 收益全来自"分给 3 条线程"（实测单次回代快 2.0–2.8×），而它之所以是**零风险**改造，
// 全靠下面两个性质 —— 任一条破了，它就从"不改变任何物理结果"变成"悄悄改了物理"：
//   ① L 的每一块都是标量 × I₃（惯性 m/h²·I、约束 ±κ·I、pin 行 I）
//      ⇒ 图是 3 个互不相连的连通分量 ⇒ 消元填充不可能跨分量 ⇒ 因子按分量块对角，
//      于是 solveComponent(c) 只碰自己那 1/3 的自由度（无竞争、无同步）；
//   ② 各分量的算术序列与整趟完全一致（列顺序、i>j 筛选、累加顺序都不变）
//      ⇒ 结果**逐位相同**（不是"1e-12 内相同"）。
// 另外 pin 覆盖也必须能按分量拆 —— 否则并行路径会漏写 2/3 的 pin 行。
// ---------------------------------------------------------------------------
namespace {

/// 组装并分解一个 20x20 扰动网格（含 pin），返回可直接求解的求解器。
std::unique_ptr<EigenDirectSolver> makeFactorizedSolver(const Mesh& m, Scalar dt,
                                                        Eigen::SparseMatrix<Scalar>& lOut) {
  assembleLeftHandSide(m, dt, 0.0, lOut);
  auto solver = std::make_unique<EigenDirectSolver>();
  solver->analyze(3 * m.vertexCount(), lOut);
  solver->factorize(lOut);
  return solver;
}

/// 一个非平凡右端：b = (M/h²)·x̂（与生产同款）再覆盖 pin 行。
Eigen::VectorXd makeComponentTestRhs(const Mesh& m, Scalar dt) {
  std::vector<Vec3> predicted = m.positions;
  for (std::size_t i = 0; i < predicted.size(); ++i) {
    predicted[i] = predicted[i] + Vec3{1e-3 * static_cast<Scalar>(i % 13),
                                       -2e-3 * static_cast<Scalar>(i % 7),
                                       3e-3 * static_cast<Scalar>(i % 5)};
  }
  Eigen::VectorXd b;
  assembleInertialRhs(m, predicted, dt, 0.0, b);
  applyPinRhs(m, b);
  return b;
}

bool sameBits(const Eigen::VectorXd& a, const Eigen::VectorXd& b) {
  if (a.size() != b.size() || a.size() == 0) return false;
  return std::memcmp(a.data(), b.data(),
                     sizeof(Scalar) * static_cast<std::size_t>(a.size())) == 0;
}

}  // namespace

TEST(globalSolveSplitsIntoThreeComponentChains) {
  Mesh m = makePerturbedGrid(20, 0.02, 2.0e3);
  Eigen::SparseMatrix<Scalar> L;
  auto solver = makeFactorizedSolver(m, 1.0 / 120.0, L);
  CHECK_MSG(solver->parallelComponents() == 3,
            "L = Ã ⊗ I₃（图是 3 个连通分量）时应当能拆成 3 个分量，实际 " +
                std::to_string(solver->parallelComponents()));
}

TEST(componentSolveMatchesFullSolveBitwise) {
  Mesh m = makePerturbedGrid(20, 0.02, 2.0e3);
  Eigen::SparseMatrix<Scalar> L;
  auto solver = makeFactorizedSolver(m, 1.0 / 120.0, L);
  const int n = 3 * m.vertexCount();
  const Eigen::VectorXd b = makeComponentTestRhs(m, 1.0 / 120.0);
  CHECK_MSG(b.cwiseAbs().maxCoeff() > 0.0, "右端不应全零，否则这个测试等于没测");

  Eigen::VectorXd xFull(n), xSplit(n);
  solver->solve(b, xFull);
  for (int c = 0; c < solver->parallelComponents(); ++c) {
    solver->solveComponent(c, b, xSplit);
  }
  CHECK_MSG(sameBits(xFull, xSplit),
            "按分量拆分求解必须与整趟**逐位相同**（每个分量的算术序列没变）");
}

TEST(componentPinOverwriteMatchesFullPinOverwrite) {
  Mesh m = makePerturbedGrid(20, 0.02, 2.0e3);
  const Eigen::VectorXd base = makeComponentTestRhs(m, 1.0 / 120.0);
  Eigen::VectorXd whole = base;
  applyPinRhs(m, whole);
  Eigen::VectorXd split = base;
  for (int c = 0; c < 3; ++c) applyPinRhsComponent(m, split, c, 3);
  CHECK_MSG(sameBits(whole, split), "三分量分别覆盖 pin 行必须与整份覆盖逐位相同");
  Eigen::VectorXd single = base;
  applyPinRhsComponent(m, single, 0, 1);
  CHECK_MSG(sameBits(whole, single), "components == 1 时 applyPinRhsComponent 必须等价于 applyPinRhs");
}

TEST(fusedProjectScatterMatchesSplitPath) {
  // Phase 4d 的回归锚点：融合路径（`projectAndScatterIntoInRegion`）与两趟路径
  // （`project` → `scatterInto`）必须**逐位相同**。
  //
  // 为什么这条断言重要：融合之所以是"零数值风险"，靠的是两条路**共用同一份算术 helper**
  // （`projectEdgeValue` / `scatterEdgeValue`）。这条断言把那个结构性保证钉死 ——
  // 将来谁把 helper 拆开、或改成不等价的写法（例如把颜色序换成边序），它会立刻变红。
  Mesh m = makePerturbedGrid(20, 0.02, 2.0e3);
  const std::size_t dim = static_cast<std::size_t>(3 * m.vertexCount());

  std::vector<Vec3> targets;
  DistanceTerm::project(m, targets);  // 两趟路径用的投影

  // 一个非平凡的基值（右端 = base + 散射；base 不能全零，否则测试没有分辨力）
  std::vector<Scalar> base(dim, 0.0);
  for (std::size_t i = 0; i < dim; ++i) {
    base[i] = 0.5 + 0.25 * static_cast<Scalar>(i % 7);
  }
  std::vector<Scalar> bSplit = base;
  std::vector<Scalar> bFused = base;

  DistanceTerm::scatterInto(m, targets, bSplit.data(), dim);  // 两趟路径（内部自己开区域）

  // 融合路径是"区域内"函数，必须在并行区域里调用（与 stepOnce 的用法一致）。
#ifdef _OPENMP
#pragma omp parallel num_threads(4)
#endif
  {
    DistanceTerm::projectAndScatterIntoInRegion(m, base.data(), bFused.data(), dim);
  }

  CHECK_MSG(maxAbsValue(bSplit) > 0.0, "右端不应全零，否则这个测试等于没测");
  CHECK_MSG(std::memcmp(bSplit.data(), bFused.data(), sizeof(Scalar) * dim) == 0,
            "融合路径与两趟路径必须逐位相同（两条路共用同一份算术 helper）");
}

TEST(shearDiagonalsAreAddedOnceAndKeepEdgeUniqueness) {
  // 剪切约束（对角边）的构造断言。它本身是普通 `Edge`，所以真正要钉的是**网格构造**：
  //   · 默认（不调用本函数）时边表与历史完全一致 —— 所有基线不受影响；
  //   · 每个四边形恰好多一条对角线，且**不产生重复边**（4-邻域边表里没有对角线，
  //     但这条断言是防止将来有人改成交叉对角线或两种对角线都加时忘了去重）。
  const int nx = 20;
  const int ny = 20;
  Mesh base = Mesh::makeGrid(nx, ny, 0.05);
  CHECK(base.edgeCount() == 2 * nx * ny - nx - ny);  // 760：默认约束集不变

  Mesh sheared = Mesh::makeGrid(nx, ny, 0.05);
  const int added = sheared.addShearDiagonals(nx, ny, 1.0e3);
  CHECK_MSG(added == (nx - 1) * (ny - 1), "每个四边形应恰好多一条对角边");
  CHECK_MSG(sheared.edgeCount() == base.edgeCount() + (nx - 1) * (ny - 1),
            "总边数 = 4-邻域边 + 对角边");
  CHECK_MSG(sheared.countDuplicateEdges() == 0, "对角边不应与已有边重复");
}

// ---------------------------------------------------------------------------
// 9) 线性（中点）弯曲约束：构造 / 平直零力 / 三点解析 / 位形无关
//
// 约束形式（`core/mesh/Mesh.h` 的 BendStencil）：沿行/列取连续三个顶点 (a,b,c)，
//   A_s x = x_a - 2 x_b + x_c，E_s = (k/2)‖A_s x‖²
// 它没有局部步、没有散射（投影恒为 0），全部影响只在这三处：
//   ① 网格生成（函数体自动重建稀疏结构/着色/关联表）；② 左端块 k·w_p·w_q·I₃ + 右端常向量；
//   ③ 残差（残差是唯一放行判据，口径必须含弯曲力 —— 见 Integrator.cpp）。
// 本节钉的是其中"能解析算出来"的部分。
// ---------------------------------------------------------------------------
namespace {

/// 手工搭一条三点直链（只有弯曲 stencil，没有距离约束），pin 两端。
/// 顶点沿 +x 等距排列；中间顶点沿 y 抬高 `eps`。
Mesh makeThreePointChain(Scalar spacing, Scalar eps, Scalar kBend) {
  const int n = 3;
  Mesh m;
  m.positions.resize(static_cast<std::size_t>(n));
  m.velocities.assign(static_cast<std::size_t>(n), Vec3{});
  m.masses.assign(static_cast<std::size_t>(n), 1.0);
  m.pinned.assign(static_cast<std::size_t>(n), 0);
  m.pinPositions.resize(static_cast<std::size_t>(n));
  for (int i = 0; i < n; ++i) {
    const Vec3 p{static_cast<Scalar>(i) * spacing, (i == 1) ? eps : Scalar{0}, Scalar{0}};
    m.positions[static_cast<std::size_t>(i)] = p;
    m.pinPositions[static_cast<std::size_t>(i)] = p;
  }
  m.restPositions = m.positions;
  BendStencil s;
  s.a = 0;
  s.b = 1;
  s.c = 2;
  s.stiffness = kBend;
  m.bends.push_back(s);
  m.pinned[0] = 1;
  m.pinned[2] = 1;
  m.buildSparsityPattern();
  return m;
}

/// 独立重算某个顶点的弯曲力（不复用 core 的实现，避免"自己验证自己"）：
///   f_v = Σ_{stencil 含 v} -k·w_v·(x_a - 2x_b + x_c)
Vec3 referenceBendForce(const Mesh& m, int v) {
  Vec3 f{0, 0, 0};
  for (const BendStencil& s : m.bends) {
    const int vertex[3] = {s.a, s.b, s.c};
    const Scalar weight[3] = {1.0, -2.0, 1.0};
    for (int t = 0; t < 3; ++t) {
      if (vertex[t] != v) continue;
      const Vec3 second = m.positions[static_cast<std::size_t>(s.a)] -
                          m.positions[static_cast<std::size_t>(s.b)] * 2.0 +
                          m.positions[static_cast<std::size_t>(s.c)];
      f += second * (-s.stiffness * weight[t]);
    }
  }
  return f;
}

}  // namespace

TEST(bendingStencilCountAndVertexIndicesAreExact) {
  // ① stencil 数量与顶点索引正确。
  //
  // 计数是可解析的：行方向每行 (nx-2) 条 × ny 行，列方向每列 (ny-2) 条 × nx 列，
  // 合计 (nx-2)*ny + (ny-2)*nx。**端点不生成**（中心差分需要两侧邻居），
  // 所以边界一圈没有弯曲约束 —— 这是标准做法（自由端/自然边界条件），不是漏了。
  const int nx = 20;
  const int ny = 14;
  const Scalar spacing = 0.05;
  const Scalar k = 777.0;

  Mesh base = Mesh::makeGrid(nx, ny, spacing);
  CHECK_MSG(base.bendCount() == 0, "默认（不调用 addBendingStencils）时不应有任何弯曲 stencil");

  Mesh m = Mesh::makeGrid(nx, ny, spacing);
  const int added = m.addBendingStencils(nx, ny, k);
  const int expected = (nx - 2) * ny + (ny - 2) * nx;
  CHECK_MSG(added == expected, "弯曲 stencil 数应为 (nx-2)*ny + (ny-2)*nx = " +
                                   std::to_string(expected) + "，实际 " + std::to_string(added));
  CHECK(m.bendCount() == expected);

  // 逐条核对顶点三元组：前 (nx-2)*ny 条是行方向、其余是列方向的。
  // 索引约定 v = j*nx + i（与 makeGrid 一致）。
  bool indicesOk = true;
  for (int j = 0; j < ny && indicesOk; ++j) {
    for (int i = 1; i + 1 < nx && indicesOk; ++i) {
      const int idx = j * (nx - 2) + (i - 1);
      const BendStencil& s = m.bends[static_cast<std::size_t>(idx)];
      indicesOk = (s.a == j * nx + i - 1) && (s.b == j * nx + i) && (s.c == j * nx + i + 1);
    }
  }
  CHECK_MSG(indicesOk, "行方向 stencil 的顶点三元组必须是同行的 (i-1, i, i+1)");

  const int rowCount = (nx - 2) * ny;
  indicesOk = true;
  for (int i = 0; i < nx && indicesOk; ++i) {
    for (int j = 1; j + 1 < ny && indicesOk; ++j) {
      const int idx = rowCount + i * (ny - 2) + (j - 1);
      const BendStencil& s = m.bends[static_cast<std::size_t>(idx)];
      indicesOk = (s.a == (j - 1) * nx + i) && (s.b == j * nx + i) && (s.c == (j + 1) * nx + i);
    }
  }
  CHECK_MSG(indicesOk, "列方向 stencil 的顶点三元组必须是同列的 (j-1, j, j+1)");

  bool stiffnessOk = true;
  for (const BendStencil& s : m.bends) {
    if (s.stiffness != k) stiffnessOk = false;
  }
  CHECK_MSG(stiffnessOk, "每条 stencil 的刚度都必须等于传入值");

  // 端点不含弯曲：中间顶点（b）的行列索引都必须落在 1..n-2 里。
  bool noBoundary = true;
  for (int r = 0; r < rowCount; ++r) {
    const int i = m.bends[static_cast<std::size_t>(r)].b % nx;
    if (i < 1 || i > nx - 2) noBoundary = false;
  }
  for (int r = rowCount; r < m.bendCount(); ++r) {
    const int j = m.bends[static_cast<std::size_t>(r)].b / nx;
    if (j < 1 || j > ny - 2) noBoundary = false;
  }
  CHECK_MSG(noBoundary, "端点（i=0/nx-1 或 j=0/ny-1）不应生成弯曲 stencil（自由端边界）");

  // 函数会自动重建拓扑级数据（照抄 addShearDiagonals 的做法）：
  // 关联表规模必须与 stencil 数一致，否则残差的 gather 会越界读。
  const BendAdjacency& adj = m.bendAdjacency();
  CHECK_MSG(adj.vertexStart.size() == static_cast<std::size_t>(m.vertexCount()) + 1,
            "vertexStart 必须覆盖全部顶点（buildSparsityPattern 应被自动调用）");
  CHECK_MSG(adj.vertexStencils.size() == static_cast<std::size_t>(m.bendCount()) * 3,
            "vertexStencils 必须是 3×stencil 数（每个 stencil 记三份）");
  CHECK(adj.vertexStencilWeight.size() == adj.vertexStencils.size());
  // 每个顶点的关联项按 stencil 号升序（这是"串行残差 == gather 残差逐位相同"的依据）
  bool sorted = true;
  for (int v = 0; v < m.vertexCount() && sorted; ++v) {
    for (uint32_t idx = adj.vertexStart[static_cast<std::size_t>(v)];
         idx + 1 < adj.vertexStart[static_cast<std::size_t>(v) + 1]; ++idx) {
      if (adj.vertexStencils[idx] >= adj.vertexStencils[idx + 1]) sorted = false;
    }
  }
  CHECK_MSG(sorted, "顶点关联 stencil 的 CSR 必须按 stencil 号升序（残差两口径逐位相同的前提）");
}

TEST(bendingForceIsZeroOnFlatConfiguration) {
  // ② 平直静止位形下弯曲力为零（A_c x = 0 ⇒ ∇E = 0）。
  // 这是"弯曲项不会给平直布料凭空加力"的定义级判据 —— 加错符号/系数会立刻在平坦
  // 位形上留下非零力，而那一项在视觉上几乎看不出来（布料会"自己抖"而已）。
  //
  // **用二进制可精确表示的坐标**（0.25 的整数倍）：这样
  // `x_a - 2.0*x_b + x_c` 是精确的 0，可以断言"**严格**等于 0"而不是"小于某阈值"。
  // 这条区别是必要的：`makeGrid(nx, ny, 0.05)` 的坐标是 i*0.05，浮点上并不精确
  // （0.05 不可精确表示），二阶差分只会到 ~1e-13 量级 —— 那种装置**验不出**系数错
  // （把 2 写成 1.9999 也能过），所以这里刻意不用它。
  Mesh m = Mesh::makeGrid(6, 5, 0.25);
  const int added = m.addBendingStencils(6, 5, 5.0e3);
  CHECK(added == (6 - 2) * 5 + (5 - 2) * 6);  // 20 + 18 = 38

  double worst = 0.0;
  for (int v = 0; v < m.vertexCount(); ++v) {
    worst = std::max(worst, static_cast<double>(length(referenceBendForce(m, v))));
  }
  CHECK_MSG(worst == 0.0, "平直静止位形下弯曲力必须严格为 0（A_c x = 0，坐标为二进制精确值）");

  // 反证（让这条断言有分辨力）：把中间一个顶点的系数算错一点点，力就必须非零。
  // 这里用"人为地把一个顶点抬高 1e-9"来构造非零二阶差分 —— 若上面的断言是
  // "力恒为 0"这种恒真式，这一条会失败。
  Mesh m2 = Mesh::makeGrid(6, 5, 0.25);
  m2.addBendingStencils(6, 5, 5.0e3);
  m2.positions[static_cast<std::size_t>(2 * 6 + 3)].y += 1e-9;
  double worst2 = 0.0;
  for (int v = 0; v < m2.vertexCount(); ++v) {
    worst2 = std::max(worst2, static_cast<double>(length(referenceBendForce(m2, v))));
  }
  CHECK_MSG(worst2 > 0.0, "抬 1e-9 之后弯曲力必须非零（否则上面那条断言没有分辨力）");
}

TEST(bendingRestoringForceMatchesAnalyticExpression) {
  // ③ 三点链的小规模解析：回复力方向与大小等于 -k·w_v·(A_c x)（逐分量），
  //    并验证"上抬中间点 ⇒ 力朝 -y"（把点拉回直线）。
  const Scalar spacing = 0.1;
  const Scalar k = 2.0e3;
  const Scalar eps = 0.01;
  Mesh m = makeThreePointChain(spacing, eps, k);

  // A_c x = x_a - 2 x_b + x_c = (0,0,0) - 2(0.1, eps, 0) + (0.2,0,0) = (0, -2eps, 0)
  // 中间顶点 b 的权重是 -2 ⇒ 力 = -k·(-2)·(0, -2eps, 0) = (0, -4k·eps, 0)。
  const Vec3 second{0.0, -2.0 * eps, 0.0};
  const Vec3 wantMid = second * (-k * (-2.0));
  // 端点 a 与 c 各拿 -k·(+1)·second = (0, +2k·eps, 0)（两侧对称，把两端朝中间拉）。
  const Vec3 wantEnd = second * (-k * 1.0);

  for (int v = 0; v < m.vertexCount(); ++v) {
    const Vec3 f = referenceBendForce(m, v);
    const Vec3 want = (v == 1) ? wantMid : wantEnd;
    CHECK_NEAR(f.x, want.x, 1e-12);
    CHECK_NEAR(f.y, want.y, 1e-12);
    CHECK_NEAR(f.z, want.z, 1e-12);
  }
  CHECK_MSG(wantMid.y < 0.0, "上抬中间点后它的弯曲力必须朝 -y（回复力，不是推得更远）");
  CHECK_NEAR(wantMid.y, -4.0 * k * eps, 1e-12);
  CHECK_NEAR(wantEnd.y, +2.0 * k * eps, 1e-12);

  // 左端矩阵：自由行 v1 的对角 = m/h² + 4k（w_b² = 4）；
  // **pinned 行必须恰好是单位行** —— 这一条同时钉住"pinned 列必须从自由行里拿掉"：
  // 漏了它矩阵会变成不定的（行列式变负），而求解器**不报错**、给出全错的解。
  const Scalar h = 1.0 / 120.0;
  Eigen::SparseMatrix<Scalar> L;
  std::vector<Scalar> bendRhs;
  assembleLeftHandSide(m, h, 0.0, L, &bendRhs);
  CHECK_NEAR(L.coeff(0, 0), 1.0, 1e-15);
  CHECK_NEAR(L.coeff(0, 3), 0.0, 0.0);
  CHECK_NEAR(L.coeff(0, 6), 0.0, 0.0);
  CHECK_NEAR(L.coeff(3, 0), 0.0, 0.0);  // 自由行里也不能留 pinned 列
  CHECK_NEAR(L.coeff(3, 3), m.masses[1] / (h * h) + 4.0 * k, 1e-9);
  CHECK_NEAR(L.coeff(4, 4), m.masses[1] / (h * h) + 4.0 * k, 1e-9);

  // 右端常向量：两个端点都被 pin、w=1 ⇒ C = 1·q_a + 1·q_c = 0.2（x 分量）⇒
  // 中间顶点的右端补偿 = -k·w_b·C = -2000·(-2)·0.2 = **+800**（x 分量）。
  // 注意符号：对 **a 端**（w_a=+1）被 pin 的情形补偿是负的，对 **中间顶点**（w_b=-2）
  // 就是正的 —— 所以"补偿项的方向"必须按 w_p 定，不能一律写成 -k·C。
  CHECK_NEAR(bendRhs[3], -k * (-2.0) * 0.2, 1e-9);
  CHECK_NEAR(bendRhs[3], 800.0, 1e-9);  // k = 2000、间距 0.1 时的具体值
  CHECK_NEAR(bendRhs[4], 0.0, 0.0);
  CHECK_NEAR(bendRhs[5], 0.0, 0.0);
  // pin 行与端点行不应收到常向量
  CHECK(bendRhs[0] == 0.0 && bendRhs[1] == 0.0 && bendRhs[2] == 0.0);
  CHECK(bendRhs[6] == 0.0 && bendRhs[7] == 0.0 && bendRhs[8] == 0.0);
}

TEST(bendingLeftHandSideIsConfigurationIndependent) {
  // ④ 位形无关：改位形后 L **逐位不变**。这是"只分解一次"的前提
  //    （docs/plan.md §2.3 不变量 1），也是选"中点形式"而不是二面角形式的全部理由。
  const int nx = 12;
  const int ny = 12;
  const Scalar k = 1.0e3;
  const Scalar h = 1.0 / 120.0;

  Mesh m = Mesh::makeGrid(nx, ny, 0.05);
  m.addBendingStencils(nx, ny, k);
  m.pinned[5 * nx + 5] = 1;
  m.pinPositions = m.positions;

  Eigen::SparseMatrix<Scalar> l1;
  std::vector<Scalar> rhs1;
  assembleLeftHandSide(m, h, 0.0, l1, &rhs1);

  // 大幅扰动位形（含横向起伏）：L 必须逐位不变。
  for (int v = 0; v < m.vertexCount(); ++v) {
    const Scalar i = static_cast<Scalar>(v % nx);
    const Scalar j = static_cast<Scalar>(v / nx);
    m.positions[static_cast<std::size_t>(v)] +=
        Vec3{0.03 * std::sin(i), 0.05 * std::cos(j), 0.02 * std::sin(i + j)};
  }
  Eigen::SparseMatrix<Scalar> l2;
  std::vector<Scalar> rhs2;
  assembleLeftHandSide(m, h, 0.0, l2, &rhs2);

  CHECK(l1.rows() == l2.rows() && l1.cols() == l2.cols());
  CHECK(l1.nonZeros() == l2.nonZeros());
  bool identical = true;
  for (int kk = 0; kk < l1.outerSize() && identical; ++kk) {
    for (Eigen::SparseMatrix<Scalar>::InnerIterator it(l1, kk); it; ++it) {
      if (it.value() != l2.coeff(it.row(), it.col())) identical = false;
    }
  }
  CHECK_MSG(identical, "改位形后弯曲的左端块必须逐位不变（这是「只分解一次」的前提）");

  // 常向量只依赖 **pin 位置**：位形变了它必须不变。
  Mesh m3 = Mesh::makeGrid(nx, ny, 0.05);
  m3.addBendingStencils(nx, ny, k);
  m3.pinned[5 * nx + 5] = 1;
  m3.pinPositions = m3.positions;
  std::vector<Scalar> rA;
  assembleBendingRhs(m3, rA);
  CHECK_MSG(rA == rhs1, "同一场景下两次装配的常向量必须逐位相同");
  for (int v = 0; v < m3.vertexCount(); ++v) {
    m3.positions[static_cast<std::size_t>(v)] += Vec3{0.1, -0.2, 0.3};
  }
  std::vector<Scalar> rB;
  assembleBendingRhs(m3, rB);
  CHECK_MSG(rA == rB, "弯曲常向量只依赖 pin 位置 ⇒ 改位形必须逐位不变");

  // 量级与方向的可分辨性：非平坦位形下弯曲力随刚度单调变强、方向始终朝 -y。
  Mesh soft = makeThreePointChain(0.1, 0.02, 10.0);
  Mesh hard = makeThreePointChain(0.1, 0.02, 10.0e3);
  const double softY = static_cast<double>(referenceBendForce(soft, 1).y);
  const double hardY = static_cast<double>(referenceBendForce(hard, 1).y);
  CHECK_MSG(hardY < softY && softY < 0.0,
            "弯曲力必须随刚度单调变强，且方向始终朝 -y（回复力）");
}

// ---------------------------------------------------------------------------
// 10) 弯曲 stencil 的采样变体（降本开关）：单向 / 隔行 / 棋盘
//
// 动机：弯曲的代价几乎全在**消元填充**上（因子 nnz 2.77×），而填充由"每个未知量被多少
// 约束耦合"决定 ⇒ 减少 stencil 数量是直接打在填充上的降本手段（docs/perf.md §13）。
// 本节钉的是"变体没有偷偷改变结构性质"：
//   · 默认选项必须与不传选项**逐条相同**（默认路径逐位不变，历史基线才算数）；
//   · 各变体生成的 stencil 仍然**只有连续三个顶点的形状**（不会退化成别的约束）；
//   · 隔行采样是完整采样的**严格子集**（不是重排、不是新增）；
//   · 棋盘里每个中心顶点**恰好出现一次**，方向严格按 `(i+j)` 的奇偶；
//   · 最重要的一条：**每个 3×3 块仍是 标量 × I₃** ⇒ `L = Ã ⊗ I₃` 与 Phase 4c 的
//     三分量并行回代对这些变体同样成立（这不是"文档承诺"，这里逐块验）。
// ---------------------------------------------------------------------------
namespace {

/// 方向判定：行 stencil 的三个顶点同行且依次相差 1；列 stencil 依次相差 nx。
/// 返回 0 = 行、1 = 列、-1 = 两者都不是（这种形状不该出现）。
int bendDirectionKind(const BendStencil& s, int nx) {
  const bool row = (s.b - s.a == 1) && (s.c - s.b == 1);
  const bool col = (s.b - s.a == nx) && (s.c - s.b == nx);
  if (row && !col) return 0;
  if (col && !row) return 1;
  return -1;
}

/// 稀疏结构里"每个块是 标量 × I₃"的逐块自检（即 `L = Ã ⊗ I₃`）。
bool allBlocksAreScalarTimesIdentity(const Mesh& m, const Eigen::SparseMatrix<Scalar>& L) {
  for (const BlockEntry& e : m.blockPattern()) {
    const int r = 3 * e.row;
    const int c = 3 * e.col;
    const Scalar d0 = L.coeff(r + 0, c + 0);
    const Scalar d1 = L.coeff(r + 1, c + 1);
    const Scalar d2 = L.coeff(r + 2, c + 2);
    if (d0 != d1 || d1 != d2) return false;
    for (int u = 0; u < 3; ++u) {
      for (int v = 0; v < 3; ++v) {
        if (u != v && L.coeff(r + u, c + v) != 0.0) return false;
      }
    }
  }
  return true;
}

}  // namespace

TEST(bendingSamplingVariantsAreExactSubsetsAndKeepTensorStructure) {
  const int nx = 20;
  const int ny = 14;
  const Scalar spacing = 0.05;
  const Scalar k = 777.0;

  // ---- ① 默认选项 == 不传选项：逐条、按序相同（默认路径逐位不变）----
  Mesh plain = Mesh::makeGrid(nx, ny, spacing);
  plain.addBendingStencils(nx, ny, k);
  Mesh dflt = Mesh::makeGrid(nx, ny, spacing);
  dflt.addBendingStencils(nx, ny, k, BendGenOptions{});
  CHECK_MSG(plain.bendCount() == dflt.bendCount(),
            "不传选项与传默认 BendGenOptions 的 stencil 数必须相同");
  bool identicalOrder = (plain.bendCount() == dflt.bendCount());
  for (int i = 0; i < plain.bendCount() && identicalOrder; ++i) {
    const BendStencil& p = plain.bends[static_cast<std::size_t>(i)];
    const BendStencil& q = dflt.bends[static_cast<std::size_t>(i)];
    if (p.a != q.a || p.b != q.b || p.c != q.c || p.stiffness != q.stiffness) {
      identicalOrder = false;
    }
  }
  CHECK_MSG(identicalOrder,
            "默认选项必须与不传选项**按序逐条相同**（生成次序是残差逐位可复现的前提）");

  // ---- ② 单向：只有该方向的 stencil，数量可解析，端点仍不生成 ----
  const int rowOnlyExpected = (nx - 2) * ny;
  const int colOnlyExpected = (ny - 2) * nx;
  for (int variant = 0; variant < 2; ++variant) {
    const bool wantRow = (variant == 0);
    const BendSampling sampling = wantRow ? BendSampling::RowsOnly : BendSampling::ColsOnly;
    Mesh m = Mesh::makeGrid(nx, ny, spacing);
    const int added = m.addBendingStencils(nx, ny, k, BendGenOptions{sampling, 1});
    const int expect = wantRow ? rowOnlyExpected : colOnlyExpected;
    CHECK_MSG(added == expect, std::string(wantRow ? "只取行" : "只取列") +
                                   "方向的 stencil 数应为 " + std::to_string(expect) + "，实际 " +
                                   std::to_string(added));
    bool dirOk = true;
    bool boundaryOk = true;
    for (const BendStencil& s : m.bends) {
      const int kind = bendDirectionKind(s, nx);
      if (kind != (wantRow ? 0 : 1)) dirOk = false;
      const int i = s.b % nx;
      const int j = s.b / nx;
      // 端点不生成：约束的是**沿 stencil 方向**的那个索引 —— 行 stencil 要求 i 是内点
      //（中心差分需要左右邻居），但它落在 j=0 / ny-1 那一行是**合法**的：那条曲率只沿 i。
      // 所以"自由端边界"是逐方向成立的，不是"必须两个索引都是内点"。
      if (wantRow ? (i < 1 || i > nx - 2) : (j < 1 || j > ny - 2)) boundaryOk = false;
      if (s.stiffness != k) dirOk = false;
    }
    CHECK_MSG(dirOk, "单向变体必须只生成该方向的 stencil（形状 = 连续三个顶点，刚度原样传递）");
    CHECK_MSG(boundaryOk, "单向变体仍不得在端点生成 stencil（自由端边界不因降本而改变）");
    // 拓扑级数据必须被自动重建（漏了会越界读）
    CHECK(m.bendAdjacency().vertexStencils.size() ==
          static_cast<std::size_t>(m.bendCount()) * 3);
  }

  // ---- ③ 隔行采样是完整采样的严格子集（步长 2、步长 3 都是）----
  auto triples = [](const Mesh& m) {
    std::vector<std::array<int, 3>> out;
    for (const BendStencil& s : m.bends) out.push_back({s.a, s.b, s.c});
    std::sort(out.begin(), out.end());
    return out;
  };
  const std::vector<std::array<int, 3>> full = triples(plain);
  for (int stride = 2; stride <= 3; ++stride) {
    Mesh sub = Mesh::makeGrid(nx, ny, spacing);
    const int added = sub.addBendingStencils(nx, ny, k, BendGenOptions{BendSampling::Standard, stride});
    const int perRow = (nx - 2 + stride - 1) / stride;  // ceil((nx-2)/stride)
    const int perCol = (ny - 2 + stride - 1) / stride;
    const int expect = perRow * ny + perCol * nx;
    CHECK_MSG(added == expect, "步长 " + std::to_string(stride) + " 的 stencil 数应为 " +
                                   std::to_string(expect) + "，实际 " + std::to_string(added));
    const std::vector<std::array<int, 3>> got = triples(sub);
    bool isSubset = std::includes(full.begin(), full.end(), got.begin(), got.end());
    CHECK_MSG(isSubset, "步长 " + std::to_string(stride) +
                            " 的结果必须是完整采样的**严格子集**（只是少取，不是重排/新增）");
    std::vector<std::array<int, 3>> uniq = got;
    uniq.erase(std::unique(uniq.begin(), uniq.end()), uniq.end());
    CHECK_MSG(uniq.size() == got.size(), "隔行采样不得产生重复 stencil（重复会双倍加刚度）");
  }

  // ---- ④ 棋盘：每个中心顶点恰好一次，方向严格按 (i+j) 的奇偶 ----
  {
    Mesh m = Mesh::makeGrid(nx, ny, spacing);
    const int added = m.addBendingStencils(nx, ny, k, BendGenOptions{BendSampling::Checkerboard, 1});
    CHECK(added == m.bendCount());
    CHECK(added > 0);

    std::vector<int> centerHits(static_cast<std::size_t>(m.vertexCount()), 0);
    bool parityOk = true;
    for (const BendStencil& s : m.bends) {
      const int i = s.b % nx;
      const int j = s.b / nx;
      const int kind = bendDirectionKind(s, nx);
      centerHits[static_cast<std::size_t>(s.b)] += 1;
      const bool even = ((i + j) % 2) == 0;
      if (even ? (kind != 0) : (kind != 1)) parityOk = false;
    }
    CHECK_MSG(parityOk, "棋盘变体里 (i+j) 为偶的中心必须取行方向、为奇的取列方向");
    int centers = 0;
    bool atMostOnce = true;
    for (int v = 0; v < m.vertexCount(); ++v) {
      if (centerHits[static_cast<std::size_t>(v)] > 1) atMostOnce = false;
      if (centerHits[static_cast<std::size_t>(v)] == 1) ++centers;
    }
    CHECK_MSG(atMostOnce, "棋盘变体里每个中心顶点最多被取一次（不得同一顶点拿两个方向）");
    CHECK_MSG(centers == added, "棋盘变体的每个 stencil 都有各自不同的中心顶点");
    // 数量按"定义"独立数一遍：偶数和对的合法中心取行 + 奇数和对的合法中心取列。
    // 注意它**不等于**单向的 (nx-2)*ny —— 行方向只覆盖 j 全部而 i 内点、列方向只覆盖 i 全部
    // 而 j 内点，两者可取的中心集合不同，差的就是四个角附近的那些格子。
    int checkerExpected = 0;
    for (int j = 0; j < ny; ++j) {
      for (int i = 1; i + 1 < nx; ++i) {
        if (((i + j) & 1) == 0) ++checkerExpected;
      }
    }
    for (int i = 0; i < nx; ++i) {
      for (int j = 1; j + 1 < ny; ++j) {
        if (((i + j) & 1) != 0) ++checkerExpected;
      }
    }
    CHECK_MSG(added == checkerExpected, "棋盘变体的 stencil 数应为 " +
                                            std::to_string(checkerExpected) + "，实际 " +
                                            std::to_string(added));
    CHECK_MSG(added >= plain.bendCount() / 2 - 2 && added <= plain.bendCount() / 2 + 2,
              "棋盘变体的规模应约为完整采样的一半（各方向各覆盖一半中心）");
  }

  // ---- ⑤ 结构性质：四种变体都必须仍是 标量 × I₃（`L = Ã ⊗ I₃`）----
  const Scalar h = 1.0 / 120.0;
  const BendGenOptions variants[] = {
      BendGenOptions{BendSampling::Standard, 1},     BendGenOptions{BendSampling::RowsOnly, 1},
      BendGenOptions{BendSampling::ColsOnly, 1},     BendGenOptions{BendSampling::Checkerboard, 1},
      BendGenOptions{BendSampling::Standard, 2},     BendGenOptions{BendSampling::Checkerboard, 2},
      BendGenOptions{BendSampling::RowsOnly, 3},
  };
  for (const BendGenOptions& opt : variants) {
    Mesh m = Mesh::makeGrid(nx, ny, spacing);
    m.addBendingStencils(nx, ny, k, opt);
    m.pinned[5 * nx + 5] = 1;
    m.pinPositions = m.positions;
    Eigen::SparseMatrix<Scalar> L;
    std::vector<Scalar> rhs;
    assembleLeftHandSide(m, h, 0.0, L, &rhs);
    const std::string tag = std::string(bendSamplingName(opt.sampling)) + "/" +
                            std::to_string(opt.stride);
    CHECK_MSG(L.rows() == 3 * m.vertexCount(),
              "变体 " + tag + " 的左端矩阵规模必须仍是 3N");
    CHECK_MSG(allBlocksAreScalarTimesIdentity(m, L),
              "变体 " + tag + " 的每个 3×3 块必须仍是 标量 × I₃（⊗ 结构与三分量并行的前提）");
  }
}

TEST(stridedBendingSamplingCannotSeeCoarseBendingSoItIsRejected) {
  // 为什么"隔行采样"这条路被否掉 —— 这里把原因钉成**可构造的反例**，而不是只写在文档里。
  //
  // 道理：一条 stencil 在中心顶点上要求 `x_i = (x_{i-1} + x_{i+1})/2`（局部共线）。
  // 只在**奇数**中心采样时，偶数顶点那一层子格子完全不受约束：把奇数顶点取成两侧偶数
  // 顶点的**中点**，则所有被采样约束的 `A_c x` **恰好为 0**，而偶数子格子可以任意弯 ——
  // 弯曲自由度整体搬到了没被采样的子格子上。所以隔行采样不是"弯曲刚度小一半"，
  // 而是**对粗尺度弯曲几乎没有阻抗**（实测：k 从 1e3 加到 1e5，应变只变 0.4 %；
  // 而完整采样变 76 %。见 docs/perf.md §13）。
  //
  // 本测试用的数字都是二进制精确值（0.5 的整数倍），所以断言是"严格等于 0"而不是阈值。
  const int nx = 9;
  const int ny = 2;
  const Scalar k = 1000.0;
  const Scalar c = 1.0;  // 偶数子格子的二次剖面 y = c·(i/2)²

  Mesh m = Mesh::makeGrid(nx, ny, 1.0);
  m.addBendingStencils(nx, ny, k, BendGenOptions{BendSampling::Standard, 2});
  // stride = 2 ⇒ 行中心只取 i = 1,3,5,7（ny = 2 ⇒ 列方向没有内点，自然 0 条）
  CHECK_MSG(m.bendCount() == 8, "9x2 网格、步长 2 只应生成两行 × 4 个行中心 = 8 条");
  for (const BendStencil& s : m.bends) {
    CHECK((s.b % nx) % 2 == 1);  // 中心必须都是奇数 i
  }

  // 偶数顶点：粗尺度**强弯曲**的二次剖面（0, c, 4c, 9c, 16c）；
  // 奇数顶点：取两侧偶数顶点的中点 ⇒ 恰好满足全部被采样的约束。
  const Scalar evenProfile[5] = {0.0, 1.0, 4.0, 9.0, 16.0};
  for (int j = 0; j < ny; ++j) {
    for (int i = 0; i < nx; ++i) {
      Scalar y = 0.0;
      if (i % 2 == 0) {
        y = evenProfile[i / 2] * c;
      } else {
        const int kk = (i - 1) / 2;
        y = 0.5 * (evenProfile[kk] + evenProfile[kk + 1]) * c;  // 0.5, 2.5, 6.5, 12.5
      }
      m.positions[static_cast<std::size_t>(j * nx + i)] = Vec3{static_cast<Scalar>(i), y, 0.0};
    }
  }

  // ① 被采样的约束严格全为 0（"看不见"这个位形）。
  Scalar worstSampled = 0.0;
  for (const BendStencil& s : m.bends) {
    const Vec3 second = m.positions[static_cast<std::size_t>(s.a)] -
                        m.positions[static_cast<std::size_t>(s.b)] * 2.0 +
                        m.positions[static_cast<std::size_t>(s.c)];
    worstSampled = std::max(worstSampled, static_cast<Scalar>(length(second)));
  }
  CHECK_MSG(worstSampled == 0.0, "隔行采样下必须能构造出「被采样约束严格全为 0」的位形");

  // ② 但它**明显不是直线**（否则这个反例没有意义）：
  // 中点的 y = 4c，而两端点连线在中点的值 = (0 + 16c)/2 = 8c ⇒ 偏离 4c（非常大）。
  const Scalar midY = m.positions[static_cast<std::size_t>(4)].y;
  const Scalar chordY = 0.5 * (m.positions[0].y + m.positions[static_cast<std::size_t>(nx - 1)].y);
  CHECK_NEAR(midY, 4.0 * c, 1e-15);
  CHECK_NEAR(chordY, 8.0 * c, 1e-15);
  CHECK_MSG(std::abs(midY - chordY) > 1.0, "该位形是强弯曲的（中点偏离弦 4c），不是直线");

  // ③ 同一个位形在**完整采样**下必须被"看见"：偶数中心 i = 2 处的二阶差分为
  // y_1 - 2y_2 + y_3 = 0.5c - 2c + 2.5c = c ≠ 0 ⇒ 完整采样确实能阻抗它。
  // 这正是"隔行采样不是完整采样的近似，而是另一个（退化的）约束集"的证据。
  const Scalar secondAtEvenCenter = m.positions[1].y - 2.0 * m.positions[2].y + m.positions[3].y;
  CHECK_NEAR(secondAtEvenCenter, c, 1e-15);
  CHECK_MSG(secondAtEvenCenter != 0.0,
            "同一个位形在完整采样下必须非零（完整采样若也看不见它，那这条断言就没有分辨力）");

  // ④ 顺带钉住"完整采样下中心必须覆盖全部内点"：正是这条保证了上面的位形无处可藏。
  // nx = 9 ⇒ 行中心 i = 1..7 共 7 个/行，其中偶数中心 3 个（2/4/6），2 行共 6 个。
  Mesh full = Mesh::makeGrid(nx, ny, 1.0);
  full.addBendingStencils(nx, ny, k, BendGenOptions{BendSampling::Standard, 1});
  int evenCenters = 0;
  for (const BendStencil& s : full.bends) {
    if ((s.b % nx) % 2 == 0) ++evenCenters;
  }
  CHECK_MSG(evenCenters == 6,
            "完整采样必须覆盖全部内点中心（含偶数中心 6 个）—— 那类位形因此无处可藏");
}

TEST_MAIN("primitives/distance_term")