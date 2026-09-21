// core/mesh/ConstraintColoring.h
// 约束（边）的图着色：**同色边两两不共享顶点**。
//
// 为什么要它（docs/parallel-refactor.md Phase 2）：
//   散射的冲突源是"同一顶点被多条边累加"。用"每线程私有全维缓冲 + 全量归约"能避开冲突，
//   但代价是 O(P·dim) 的内存流量（实测 40×40 每次 1.38 MB，是有效工作的 ~9 倍 ≈ P/2）。
//   按颜色分组执行则**完全不需要归约**，而且同一顶点的累加顺序被固定成颜色序 0,1,…,C-1
//   ⇒ 结果与线程数、分块、调度无关（位级可复现）。
//
// 着色对象是**约束**，不是顶点：冲突图的节点是边，两节点相邻 ⟺ 共享一个顶点
// （即网格图的**线图**）。用并行的极大独立集迭代（Luby 风格）着色，
// 权重取"边索引的确定性哈希"——因此同一网格在任何时刻、任何线程数下都得到**同一套颜色**。
//
// 这是**拓扑级数据**：与位形、刚度、阻尼、pin 掩码都无关，与 `Mesh::blockPattern()`
// 同生命周期（由 `Mesh::buildSparsityPattern()` 一并构建）。
#pragma once

#include <cstdint>
#include <vector>

namespace pd {

class Mesh;  // 前向声明：着色结果本身不需要 Mesh 的完整定义（避免头文件循环依赖）

/// 约束着色结果。
///
/// 三份数据：
///   · colorOf[e]  —— 边 e 的颜色（0..colorCount-1）
///   · order[]     —— 按颜色分组的边索引（CSR：第 c 色的边是 order[start[c] .. start[c+1])）
///   · start[]     —— 大小 colorCount+1 的颜色区间
/// 分组存放是为了让"按颜色扫描"成为顺序访问（缓存友好）。
struct ConstraintColoring {
  std::vector<uint32_t> colorOf;
  std::vector<uint32_t> start;  ///< 颜色区间（CSR），大小 colorCount+1
  std::vector<uint32_t> order;  ///< 按颜色分组的边索引
  int colorCount = 0;
  /// 是否已经构建过。**空网格（0 条边）也是"已构建"**，用它而不是 colorCount 来判断
  /// "是否需要构建"，才能让懒构建的检查恒为 O(1) 且不重复构建。
  bool built = false;

  /// 顶点 → 关联边（CSR，拓扑级）。着色用它做冲突判定；
  /// **散射的 gather 路径也用它**（每个顶点读自己的关联边、只写自己的槽位）。
  /// 放这里是因为两者生命周期完全相同（都只依赖拓扑），共享一份缓存与失效规则。
  std::vector<uint32_t> vertexStart;   ///< 大小 vertexCount+1
  std::vector<uint32_t> vertexEdges;   ///< 每条关联边出现两次（两端各一次）

  /// 有没有颜色类可用（0 条边时为 false，散射的逐色循环会自然跳过）。
  bool empty() const { return colorCount <= 0; }

  /// 第 c 色的边在 order 中的区间。
  uint32_t begin(int c) const { return start[static_cast<std::size_t>(c)]; }
  uint32_t end(int c) const { return start[static_cast<std::size_t>(c) + 1]; }
};

/// 并行地给约束着色（Jones–Plassmann 风格：优先级 + "最小可用颜色"）。
///
/// 算法（每轮并行一趟，串行应用一趟）：
///   1. 给每条边一个**确定性优先级**（索引的哈希）与全序比较 (优先级, 索引)；
///   2. 每轮里，"所有邻居都已着色"或"在剩余集合中优先级高于所有未着色邻居"的边
///      取**最小可用颜色**（避开已着色邻居用掉的颜色）；
///   3. 未着色的边进入下一轮，直到全部着色。
///
/// 为什么不是"纯极大独立集迭代"（每轮一种新颜色）：那样颜色数会明显多于必要值 ——
/// 实测 40×40 网格用了 **15** 色（理论最少 4），而每次散射的颜色轮数 = 颜色数，
/// 每色一个 barrier，颜色数直接变成热路径成本。加入"取最小可用颜色"后颜色数 ≤ Δ+1。
///
/// 终止性与正确性：
///   · 每轮至少着色一条边（剩余集合里优先级最大者必是局部极大）⇒ 必然终止；
///   · 同一轮被着色的边两两不相邻（相邻两条不可能同时满足上面两个条件之一）
///     ⇒ 候选颜色只需避开"上一轮及更早"的邻居颜色，最终一定是**恰当着色**；
///   · 每个顶点的颜色数不超过它的度数 ⇒ 总颜色数 ≤ Δ+1。
///
/// 确定性：优先级是边索引的纯函数、每轮按索引序迭代、并行部分只做"逐元素写"，
/// 所以结果与线程数完全无关（同网格两次着色逐位相同）。
///
/// 复杂度：O(轮数 · E · Δ)，一次性成本，不进帧预算。
void buildConstraintColoring(const Mesh& mesh, ConstraintColoring& out);

}  // namespace pd
