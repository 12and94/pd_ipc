// core/mesh/Mesh.cpp
#include "core/mesh/Mesh.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <unordered_set>

namespace pd {

Mesh Mesh::makeGrid(int nx, int ny, Scalar spacing) {
  Mesh m;
  if (nx < 2 || ny < 2) return m;

  m.positions.resize(static_cast<std::size_t>(nx) * ny);
  m.restPositions.resize(m.positions.size());
  m.velocities.assign(m.positions.size(), Vec3{});
  m.pinned.assign(m.positions.size(), 0);
  m.pinPositions.resize(m.positions.size());

  // 网格铺在 XZ 平面上（y 轴朝上，与渲染的相机约定一致）。
  for (int j = 0; j < ny; ++j) {
    for (int i = 0; i < nx; ++i) {
      const int v = j * nx + i;
      const Vec3 p{static_cast<Scalar>(i) * spacing, 0.0, static_cast<Scalar>(j) * spacing};
      m.positions[static_cast<std::size_t>(v)] = p;
      m.restPositions[static_cast<std::size_t>(v)] = p;
    }
  }

  // 4-邻域距离约束：水平边 + 竖直边。对角边刻意不加（约束越少越好）。
  auto addEdge = [&](int a, int b) {
    Edge e;
    e.a = a;
    e.b = b;
    e.restLength = length(m.restPositions[static_cast<std::size_t>(b)] -
                          m.restPositions[static_cast<std::size_t>(a)]);
    e.stiffness = 0.0;  // 由调用方或场景统一设置
    m.edges.push_back(e);
  };
  for (int j = 0; j < ny; ++j) {
    for (int i = 0; i < nx; ++i) {
      const int v = j * nx + i;
      if (i + 1 < nx) addEdge(v, v + 1);
      if (j + 1 < ny) addEdge(v, v + nx);
    }
  }

  // 质量按"每个顶点的面积份额"给出；对规则网格等价于 density * spacing^2。
  m.computeMassesFromEdges(1.0);
  m.buildSparsityPattern();
  return m;
}

bool Mesh::loadObj(const std::string& path, Mesh& out, Scalar stiffness) {
  std::ifstream in(path);
  if (!in) return false;

  std::vector<Vec3> verts;
  std::vector<std::array<int, 3>> tris;
  std::string line;
  while (std::getline(in, line)) {
    std::istringstream ss(line);
    std::string tag;
    ss >> tag;
    if (tag == "v") {
      Scalar x = 0, y = 0, z = 0;
      ss >> x >> y >> z;
      verts.push_back(Vec3{x, y, z});
    } else if (tag == "f") {
      // 支持 "f a/b/c" 与 "f a//c" 两种写法，索引从 1 开始（负索引不支持）。
      std::array<int, 3> idx{};
      int count = 0;
      std::string token;
      while (count < 3 && ss >> token) {
        const std::size_t slash = token.find('/');
        const std::string first = (slash == std::string::npos) ? token : token.substr(0, slash);
        try {
          idx[static_cast<std::size_t>(count)] = std::stoi(first) - 1;
        } catch (...) {
          return false;
        }
        ++count;
      }
      if (count == 3) tris.push_back(idx);
    }
  }
  if (verts.empty() || tris.empty()) return false;

  out = Mesh{};
  out.positions = verts;
  out.restPositions = verts;
  out.velocities.assign(verts.size(), Vec3{});
  out.pinned.assign(verts.size(), 0);
  out.pinPositions = verts;
  out.triangleIndices = tris;

  // 由三角形的三条边去重得到距离约束集合。
  std::unordered_set<uint64_t> seen;
  auto key = [](int a, int b) {
    const uint32_t lo = static_cast<uint32_t>(std::min(a, b));
    const uint32_t hi = static_cast<uint32_t>(std::max(a, b));
    return (static_cast<uint64_t>(lo) << 32) | hi;
  };
  for (const auto& t : tris) {
    for (int k = 0; k < 3; ++k) {
      const int a = t[static_cast<std::size_t>(k)];
      const int b = t[static_cast<std::size_t>((k + 1) % 3)];
      if (a == b) continue;
      if (!seen.insert(key(a, b)).second) continue;
      Edge e;
      e.a = a;
      e.b = b;
      e.restLength = length(out.restPositions[static_cast<std::size_t>(b)] -
                            out.restPositions[static_cast<std::size_t>(a)]);
      e.stiffness = stiffness;
      out.edges.push_back(e);
    }
  }

  out.computeMassesFromSurfaces(1.0);
  out.buildSparsityPattern();
  return true;
}

void Mesh::computeMassesFromSurfaces(Scalar density) {
  const int n = vertexCount();
  masses.assign(static_cast<std::size_t>(n), 0.0);

  // 需要三角形面表；若调用方没有提供，则退化为按边计质量。
  if (triangleIndices.empty()) {
    computeMassesFromEdges(density);
    return;
  }
  for (const auto& t : triangleIndices) {
    const Vec3& p0 = restPositions[static_cast<std::size_t>(t[0])];
    const Vec3& p1 = restPositions[static_cast<std::size_t>(t[1])];
    const Vec3& p2 = restPositions[static_cast<std::size_t>(t[2])];
    const Scalar area = 0.5 * length(cross(p1 - p0, p2 - p0));
    for (int k = 0; k < 3; ++k) masses[static_cast<std::size_t>(t[static_cast<std::size_t>(k)])] +=
        density * area / 3.0;
  }
  // 退化保护：质量必须严格为正，否则全局矩阵失去正定性保证。
  Scalar minPositive = 1e300;
  for (Scalar m : masses) {
    if (m > 1e-12) minPositive = std::min(minPositive, m);
  }
  if (minPositive > 1e299) minPositive = 1.0;
  for (Scalar& m : masses) {
    if (m < 1e-12) m = minPositive;
  }
}

void Mesh::computeMassesFromEdges(Scalar density) {
  // 用静止长度作为"面积份额"的代理：规则网格下与三角形法一致（每个顶点 4 条半长边）。
  const int n = vertexCount();
  masses.assign(static_cast<std::size_t>(n), 0.0);
  for (const auto& e : edges) {
    const Scalar share = density * e.restLength * e.restLength * 0.5;
    masses[static_cast<std::size_t>(e.a)] += share;
    masses[static_cast<std::size_t>(e.b)] += share;
  }
  Scalar minPositive = 1e300;
  for (Scalar m : masses) {
    if (m > 1e-12) minPositive = std::min(minPositive, m);
  }
  if (minPositive > 1e299) minPositive = 1.0;
  for (Scalar& m : masses) {
    if (m < 1e-12) m = minPositive;
  }
}

void Mesh::buildSparsityPattern() {
  // 结构 = 对角块 + 每条约束引入的两个非对角块。
  // 该结构只依赖拓扑，与位形、刚度无关，因此只需构建一次。
  blocks_.clear();
  const int n = vertexCount();
  blocks_.reserve(static_cast<std::size_t>(n) + 2 * edges.size());
  for (int i = 0; i < n; ++i) blocks_.push_back(BlockEntry{i, i});
  for (const auto& e : edges) {
    blocks_.push_back(BlockEntry{e.a, e.b});
    blocks_.push_back(BlockEntry{e.b, e.a});
  }
  std::sort(blocks_.begin(), blocks_.end(), [](const BlockEntry& p, const BlockEntry& q) {
    return p.row != q.row ? p.row < q.row : p.col < q.col;
  });
  blocks_.erase(std::unique(blocks_.begin(), blocks_.end(),
                            [](const BlockEntry& p, const BlockEntry& q) {
                              return p.row == q.row && p.col == q.col;
                            }),
                blocks_.end());
}

Scalar Mesh::totalMass() const {
  Scalar s = 0.0;
  for (Scalar m : masses) s += m;
  return s;
}

Scalar Mesh::maxRelativeStrain() const {
  Scalar worst = 0.0;
  for (const auto& e : edges) {
    const Scalar len = length(positions[static_cast<std::size_t>(e.a)] -
                              positions[static_cast<std::size_t>(e.b)]);
    if (e.restLength <= 0.0) continue;
    worst = std::max(worst, std::fabs(len - e.restLength) / e.restLength);
  }
  return worst;
}

Scalar Mesh::meanRelativeStrain() const {
  if (edges.empty()) return 0.0;
  Scalar sum = 0.0;
  for (const auto& e : edges) {
    const Scalar len = length(positions[static_cast<std::size_t>(e.a)] -
                              positions[static_cast<std::size_t>(e.b)]);
    if (e.restLength <= 0.0) continue;
    sum += std::fabs(len - e.restLength) / e.restLength;
  }
  return sum / static_cast<Scalar>(edges.size());
}

int Mesh::countDuplicateEdges() const {
  std::unordered_set<uint64_t> seen;
  int dup = 0;
  for (const auto& e : edges) {
    const uint32_t lo = static_cast<uint32_t>(std::min(e.a, e.b));
    const uint32_t hi = static_cast<uint32_t>(std::max(e.a, e.b));
    const uint64_t k = (static_cast<uint64_t>(lo) << 32) | hi;
    if (!seen.insert(k).second) ++dup;
  }
  return dup;
}

}  // namespace pd
