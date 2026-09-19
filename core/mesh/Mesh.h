// core/mesh/Mesh.h
// 布料网格：顶点、距离约束（弹簧边）、质量、pinned 标记，以及全局矩阵的稀疏结构。
//
// 首期约束集刻意只保留"距离约束 + pin"（见 docs/plan.md 决策 D1）：
// 约束越少，链路里每个数字都能手算核对。
#pragma once

#include <array>
#include <cstdint>
#include <string>
#include <vector>

#include "core/math/Vec3.h"
#include "core/solver/SparsePattern.h"

namespace pd {

/// 一条距离约束（弹簧）：连接顶点 a 与 b，静止长度 restLength，刚度 stiffness。
/// 首期每个顶点对最多一条边（构造时去重），因此不需要多点约束的权重系数。
struct Edge {
  int a = -1;
  int b = -1;
  Scalar restLength = 0.0;
  Scalar stiffness = 0.0;  ///< 即文档中的 κ_c（能量系数，不是半能量系数）
};

class Mesh {
 public:
  // ---- 顶点数据 ----
  std::vector<Vec3> positions;      ///< 当前位置 x（每步求解后的结果）
  std::vector<Vec3> velocities;     ///< 当前速度 v
  std::vector<Vec3> restPositions;  ///< 静止位置（用于计算质量与可选的可视化）
  std::vector<Scalar> masses;       ///< 顶点质量 m_v
  std::vector<uint8_t> pinned;      ///< 1 = 被钉住（不参与求解，位置由 pinPositions 决定）
  std::vector<Vec3> pinPositions;   ///< pinned 顶点的目标位置（拖拽把手时会改变）

  // ---- 约束 ----
  std::vector<Edge> edges;

  /// 三角形面表（仅从 OBJ 载入时有值）。距离约束由它去重生成，质量也可由它计算。
  std::vector<std::array<int, 3>> triangleIndices;

  int vertexCount() const { return static_cast<int>(positions.size()); }
  int edgeCount() const { return static_cast<int>(edges.size()); }
  bool isPinned(int v) const { return pinned[static_cast<std::size_t>(v)] != 0; }

  // ---- 构造 ----
  /// 生成规则网格布料，只建立 4-邻域距离约束。
  /// 顶点数 = nx*ny，边数应为 2*nx*ny - nx - ny（测试里会断言这个计数）。
  static Mesh makeGrid(int nx, int ny, Scalar spacing);

  /// 从 OBJ 读取三角网格，再按"共享边的两个顶点"建立距离约束。
  /// 返回 false 表示文件不可读或没有可用几何。
  static bool loadObj(const std::string& path, Mesh& out, Scalar stiffness);

  /// 计算顶点质量：m_v = density * Σ_{三角形含 v} 面积 / 3。
  /// 注意质量必须来自**静止几何**，与当前位形无关（否则全局矩阵会随位形变化）。
  void computeMassesFromSurfaces(Scalar density);

  /// 建立顶点质量：`m_v = density * 每个顶点的面积份额`。
  /// 对只给出边的模型（例如两顶点弹簧），调用方可直接设置 masses。
  void computeMassesFromEdges(Scalar density);

  /// 从边表构建 3N x 3N 稀疏结构（唯一非零块列表，行主序）。
  /// 该结构**只依赖拓扑**，因此只需构建一次即可长期复用。
  void buildSparsityPattern();

  const std::vector<BlockEntry>& blockPattern() const { return blocks_; }

  // ---- 统计与校验 ----
  Scalar totalMass() const;
  /// 累计距离约束违反量：Σ_c |‖x_a - x_b‖ - l_c| / l_c（用于验收与诊断）。
  Scalar maxRelativeStrain() const;
  Scalar meanRelativeStrain() const;
  /// 边去重计数校验：检查是否存在重复顶点对（返回重复条数）。
  int countDuplicateEdges() const;

 private:
  std::vector<BlockEntry> blocks_;  ///< 稀疏结构：所有非零 3x3 块的位置
};

}  // namespace pd
