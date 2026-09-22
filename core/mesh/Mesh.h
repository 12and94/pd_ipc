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
#include "core/mesh/ConstraintColoring.h"
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

/// **线性（中点）弯曲约束**：沿网格一行/一列取**连续三个顶点** (a, b, c)（b 在中间），
///
///     A_c x = x_a - 2 x_b + x_c        （二阶中心差分，即离散化的二阶导数）
///     E_c(x) = (k/2)·‖A_c x‖²          （k = stiffness）
///
/// 等价说法：在**顶点 b 与两邻居的中点之间**挂一根静止长度为 0 的弹簧
/// （中点 = (x_a+x_c)/2，所以 x_b - 中点 = -A_c x/2，能量 = k/2·‖A_c x‖² 差一个系数约定）。
/// 它是离散化的欧拉–伯努利弯曲能：`A_c x` 就是曲率的离散代理。
///
/// ---- 选这个形式（而不是二面角弯曲）的唯一理由：它让下面四条性质全部保住 ----
///
///   1. **A_c 是常数矩阵**（三个权重 {1,-2,1} 与位形无关）⇒ 左端矩阵与位形无关
///      ⇒ 仍然"只分解一次"（项目的核心前提，见 docs/plan.md §2.3）。
///   2. 每个 3×3 块是 `k·w_p·w_q·I₃`（w 取 {1,-2,1}）⇒ **仍是标量 × I₃**
///      ⇒ `L = Ã ⊗ I₃` 结构保住 ⇒ Phase 4c 的"按 x/y/z 三分量并行回代"继续有效
///      （docs/perf.md §9；`pd_solvecomp` 的跨分量填充必须仍为 0）。
///   3. 约束的流形是 `{x : A_c x = 0}`（一个**线性子空间**）⇒ 局部步的投影 `p_c ≡ 0`
///      ⇒ **局部步无事可做、散射不需要加任何项、边着色一个字都不用改**。
///   4. 右端只多一个**与位形无关的常向量**，而且只有 stencil 里含 pinned 顶点时才非零
///      （见 Assembler.cpp 的推导）⇒ 它是"装配期算一次"的量，不进每迭代热路径。
///
/// 对比二面角弯曲（docs/plan.md §2.2.5）：那个形式的 ∇θ 依赖位形 ⇒ 每帧局部步都要重算、
/// 投影不是 0、右端每迭代都变。本形式用"曲率是二阶差分"这个更朴素的离散化换掉了那一切。
///
/// `b` 是中间顶点（受弯的那个），不是端点顺序意义上的"中间值"。
struct BendStencil {
  int a = -1;
  int b = -1;
  int c = -1;
  Scalar stiffness = 0.0;  ///< 即文档中的 k（能量系数，不是半能量系数）
};

/// **弯曲 stencil 的采样形状**（`addBendingStencils` 的代价–质量折中开关）。
///
/// 背景：弯曲的代价几乎全在**消元填充**上（40×40 / k=1000：L nnz +78 %，而因子 nnz +177 %，
/// 单子步 +73 %，见 docs/perf.md §12）。填充由"每个未知量被多少约束耦合"决定，而弯曲
/// stencil 的耦合半径是 2（连着三个顶点），所以**减少 stencil 数量是直接打在填充上的手段**，
/// 而且它完全不碰架构性质：stencil 变少 ⇒ ⊗ 结构、只分解一次、三分量并行全都自动成立。
///
/// 代价是**弯曲响应的各向异性/离散化误差**：stencil 少了，等效抗弯刚度也变软，
/// 要按采样密度把刚度放大回去（密度 1/s ⇒ 刚度 ×s 量级），这一步必须实测标定。
///
/// 注意：无论取哪种采样，**端点依旧不生成**（自由端边界，见下）。
enum class BendSampling {
  /// 行、列都取，中心顶点全取 —— 今天的行为，也是默认值。
  Standard = 0,
  /// 只取**行方向**（沿 i 的连续三个顶点）：stencil 数 ≈ 1/2，只有横向纤维抗弯。
  RowsOnly = 1,
  /// 只取**列方向**（沿 j 的连续三个顶点）：stencil 数 ≈ 1/2，只有纵向纤维抗弯。
  ColsOnly = 2,
  /// **棋盘交替**：每个中心顶点只在"行/列"里取一个方向，方向由 `(i+j)` 的奇偶决定
  /// ⇒ stencil 数 ≈ 1/2，但两个方向都还在（各覆盖约一半的顶点），各向异性比单向小。
  Checkerboard = 3,
};

/// `BendSampling` 的可读名字（用于日志/审计打印）。
const char* bendSamplingName(BendSampling s);

/// 弯曲 stencil 的生成选项。默认值 = `Standard` + 步长 1，与不传本结构时**逐位相同**。
struct BendGenOptions {
  BendSampling sampling = BendSampling::Standard;

  /// 中心顶点的采样步长（`>= 1`）：`2` = 隔一个取一个，stencil 数 ≈ 1/2；`3` ≈ 1/3。
  /// 行方向按 `i` 计数、列方向按 `j` 计数，都从第一个可用中心（`i == 1` / `j == 1`）起算，
  /// 因此 `stride = 2` 的结果**是** `stride = 1` 结果的**严格子集**（测试里有断言）。
  int stride = 1;
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

  /// 线性（中点）弯曲约束（见上面的 `BendStencil`）。**默认空 = 关闭**。
  std::vector<BendStencil> bends;

  /// 三角形面表（仅从 OBJ 载入时有值）。距离约束由它去重生成，质量也可由它计算。
  std::vector<std::array<int, 3>> triangleIndices;

  int vertexCount() const { return static_cast<int>(positions.size()); }
  int edgeCount() const { return static_cast<int>(edges.size()); }
  int bendCount() const { return static_cast<int>(bends.size()); }
  bool isPinned(int v) const { return pinned[static_cast<std::size_t>(v)] != 0; }

  // ---- 构造 ----
  /// 生成规则网格布料，只建立 4-邻域距离约束。
  /// 顶点数 = nx*ny，边数应为 2*nx*ny - nx - ny（测试里会断言这个计数）。
  static Mesh makeGrid(int nx, int ny, Scalar spacing);

  /// **剪切约束（可选）**：给规则网格的每个四边形加一条**对角边**，刚度单独给。
  ///
  /// 为什么这样实现剪切：对角边本身就是一条普通距离约束（`Edge`），所以装配、边着色、
  /// 投影/散射、残差、以及 Phase 4c 的"三分量并行"全都**不需要任何改动**；
  /// 而它恰好提供抗剪切 —— 四边形此前可以自由塌成平行四边形，加了对角线就不能了。
  ///
  /// 约定：索引 `v = j*nx + i`；每个四边形取"左下 → 右上"那条对角线；静止长度取
  /// **静止位形**下的实际距离（规则网格上即 spacing·√2）。4-邻域边表里不存在对角线，
  /// 因此不会产生重复边。
  ///
  /// **调用后本函数会自动重建稀疏结构与约束着色/关联表**（拓扑级数据），
  /// 所以调用方不需要额外动作、也不会因为漏掉重建而踩越界读（2026-09-22 踩过一次）。
  /// 返回新增的边数。
  int addShearDiagonals(int nx, int ny, Scalar stiffness);

  /// **线性（中点）弯曲约束（可选）**：沿网格的每一行与每一列生成"连续三个顶点"的
  /// 弯曲 stencil（见 `BendStencil`），刚度单独给。
  ///
  /// 生成规则（索引约定 `v = j*nx + i`，与 `makeGrid` 一致）：
  ///   · **行方向**：固定 j，对 i = 1 .. nx-2 生成 (j*nx+i-1, j*nx+i, j*nx+i+1) —— 共 (nx-2) 个 × ny 行；
  ///   · **列方向**：固定 i，对 j = 1 .. ny-2 生成 ((j-1)*nx+i, j*nx+i, (j+1)*nx+i) —— 共 (ny-2) 个 × nx 列；
  ///   · 合计 `(nx-2)*ny + (ny-2)*nx` 条。
  ///
  /// **端点不生成**（这是本次实现的明确取舍，不是漏了）：网格边界上的顶点只有一侧邻居，
  /// 中心差分需要两侧 → 边界取不到二阶差分。标准做法是让边界成为**自由端**（零弯矩、
  /// 自然边界条件）。所以边界的一圈不会自己挺起来；要"夹住边界"得另加约束
  ///（例如把边界一圈 pin 住，或将来把单侧差分加进来）。
  ///
  /// **调用后本函数会自动重建稀疏结构与约束着色/关联表**（拓扑级数据），
  /// 照抄 `addShearDiagonals` 的做法：漏掉重建会踩越界读，那条注释记录了实际踩过的坑。
  /// `opt` 只影响**取哪些 stencil**（见 `BendGenOptions`）；默认值即"行列全取、步长 1"。
  /// 返回新增的 stencil 数。
  int addBendingStencils(int nx, int ny, Scalar stiffness,
                         const BendGenOptions& opt = BendGenOptions{});

  /// 从 OBJ 读取三角网格，再按"共享边的两个顶点"建立距离约束。
  /// 返回 false 表示文件不可读或没有可用几何。
  static bool loadObj(const std::string& path, Mesh& out, Scalar stiffness);

  /// 计算顶点质量：m_v = density * Σ_{三角形含 v} 面积 / 3。
  /// 注意质量必须来自**静止几何**，与当前位形无关（否则全局矩阵会随位形变化）。
  void computeMassesFromSurfaces(Scalar density);

  /// 建立顶点质量：`m_v = density * 每个顶点的面积份额`。
  /// 对只给出边的模型（例如两顶点弹簧），调用方可直接设置 masses。
  void computeMassesFromEdges(Scalar density);

  /// 从边表构建 3N x 3N 稀疏结构（唯一非零块列表，行主序），
  /// **并顺带建好约束着色**（见下面的说明）。
  /// 两者都只依赖拓扑，因此生命周期相同、一起失效：**改拓扑后必须重新调用本函数**。
  void buildSparsityPattern();

  const std::vector<BlockEntry>& blockPattern() const { return blocks_; }

  /// 约束着色（同色边互不相邻），供散射按颜色分组执行 —— 见 `ConstraintColoring.h`。
  /// 由 `buildSparsityPattern()` 一并构建：**它也是拓扑级数据**（与位形、刚度、阻尼、
  /// pin 掩码都无关），所以任何构造路径只要调了 `buildSparsityPattern()` 就已经有了。
  const ConstraintColoring& constraintColoring() const { return coloring_; }

  /// **顶点 → 关联弯曲 stencil** 的 CSR（含该顶点在自己 stencil 里的权重 w ∈ {1,-2,1}）。
  /// 同样由 `buildSparsityPattern()` 一并构建（拓扑级数据，与着色共享同一套失效规则）。
  /// 唯一使用者是**残差**：弯曲力必须算进残差（残差是唯一放行判据，见 Integrator.cpp）。
  const BendAdjacency& bendAdjacency() const { return bendAdjacency_; }

  /// 确保着色已构建（幂等，O(1) 检查）。给"手搭 Mesh 但没调 buildSparsityPattern()"
  /// 的调用方兜底（测试、`_verify/` 程序里很常见）：**必须在进入并行区域之前调用**，
  /// 否则会在区域里一边构建一边被别的线程读。散射的两条入口都已经在处理。
  void ensureConstraintColoring() const;

  /// 确保"顶点 → 弯曲 stencil"关联表已构建（幂等，O(B)）。
  ///
  /// 为什么与 `ensureConstraintColoring()` 分开、不能合成一个：两者的失效形态不同 ——
  /// 着色/边关联表由**边表**驱动，弯曲关联表由 **`bends`** 驱动。一个只加边、另一个只加
  /// stencil，合并后就没法判断该不该重建（"看边数没变所以跳过"会把新加的 stencil 漏掉）。
  ///
  /// 触发条件用 `vertexStart.empty()`：`buildSparsityPattern()` 建出来的表至少有
  /// "vertexCount+1 个零"，所以"空 vector"唯一地表示"从没建过"。
  /// **必须在进入并行区域之前调用**（同上）。
  void ensureBendAdjacency();

  // ---- 统计与校验 ----
  Scalar totalMass() const;
  /// 累计距离约束违反量：Σ_c |‖x_a - x_b‖ - l_c| / l_c（用于验收与诊断）。
  Scalar maxRelativeStrain() const;
  Scalar meanRelativeStrain() const;
  /// 边去重计数校验：检查是否存在重复顶点对（返回重复条数）。
  int countDuplicateEdges() const;

 private:
  std::vector<BlockEntry> blocks_;         ///< 稀疏结构：所有非零 3x3 块的位置
  mutable ConstraintColoring coloring_;    ///< 约束着色（拓扑级；可懒构建，故 mutable）
  BendAdjacency bendAdjacency_;            ///< 顶点→stencil 关联表（拓扑级，同上）
};

}  // namespace pd
