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

TEST_MAIN("primitives/distance_term")