// core/sim/Scene.h
// 场景配置与被仿真的整体状态（网格 + 求解器 + 计时统计）。
#pragma once

#include <cstdint>
#include <memory>
#include <string>

#include "core/mesh/Mesh.h"
#include "core/solver/IGlobalSolver.h"

namespace pd {

struct SceneConfig {
  // ---- 时间 ----
  Scalar dt = 1.0 / 120.0;  ///< 子步长 h
  int substepsPerFrame = 2; ///< 每帧子步数（帧预算自适应可调档）
  Vec3 gravity{0.0, -9.81, 0.0};

  // ---- 材料 ----
  /// 距离约束刚度 κ（能量系数）。**默认值 2300 来自刚度标定，不是拍的。**
  ///
  /// 标定依据（见 docs/pd-convergence.md §4 与 _verify/stiffness_calibration_record.cpp）：
  /// 悬挂布料的顶端应变 ε 与 κ 满足实测关系 `ε·κ ≈ ρg·(N·ℓ)·(...)`，在本项目
  /// 默认网格（60×60、间距 0.02、密度 1）下实测 **ε·κ ≈ 23.05**。取 ε = 1%
  /// （棉/帆布的弹性伸长量级）即得 κ ≈ 2305。
  ///
  /// **为什么不能凭手感调大**：κ 同时决定精度与代价。
  ///   · 物理上 ε = 1% 才是布料；ε = 0.01% 是"比钢还刚"（旧默认 1e5 即此量级）
  ///   · 数值上 PD 的预条件子把横向刚度当作 κ，而真实横向刚度 ∝ 张力 = κ·ε，
  ///     故预条件子**高估横向刚度约 1/ε 倍** —— ε 越小失配越大、外层迭代越慢。
  ///     标定到物理 ε 同时把两者摆正。
  /// 实测代价（60×60，40 次迭代上限）：κ=1e5 时稳态残差 > g（未平衡）；
  /// κ=2305 时残差约 0.1g，且放宽门槛后 1 次迭代即可（58 FPS）。
  Scalar stiffness = 2.3e3;
  Scalar density = 1.0;      ///< 面密度（决定顶点质量）
  Scalar velocityDamping = 0.02;  ///< 速度阻尼系数 k_d

  // ---- PD 迭代 ----
  int maxIterations = 10;          ///< 每子步最大 PD 迭代数

  /// 相对收敛判据：`|x^{k+1}-x^k|∞ <= relTolerance · scale`。
  ///
  /// **scale 的定义见 Integrator.cpp：它取"本子步已发生的位移量级"**
  ///   scale = max( max|x̂ - xⁿ|∞ , h²|g| )
  /// 而不是坐标量级。判据是"本次迭代相对本步已有位移是否可忽略"。
  ///
  /// 为什么不用坐标量级（这是本项目修过的一个缺陷）：
  /// 若写成 `diff <= relTolerance · max|x̂|∞`（坐标量级），则
  ///   · 放行的单步位移 = relTolerance × 坐标量级，**与刚度无关**；
  ///   · 高刚度下真实位移远小于这个门槛，第 1 次迭代就判"收敛"退出，
  ///     留下残差并在后续帧累加，表现为布料缓慢漂移；
  ///   · 更糟的是**判据不再平移不变**：把同一场景平移到远处，坐标量级变大会
  ///     让判据突然变松，收敛质量随摆放位置变化。
  Scalar relTolerance = 1.0e-3;

  /// 绝对收敛判据（单位 m）：单次迭代位移小于它即认为收敛。
  /// **≤0 表示禁用**（不设单独标志位，就用 "> 0" 判断，避免"传 0 反而恒成立"的坑）。
  ///
  /// 注意 1e-12 在 3600 顶点规模下实际不可达，所以真正起作用的是 relTolerance。
  /// 曾经把它放宽到 1e-6 想让高刚度问题"好转"，那是错的：
  /// 近静止时 relTolerance·scale ≈ 1e-5 × 6.8e-4 ≈ 6.8e-9 m，比 1e-6 紧约 147 倍，
  /// 于是绝对判据反而成了**主导**的提前退出条件，把问题掩盖掉而不是解决。
  /// 收敛不足是外层迭代能力问题，不是门槛松紧问题（见 tests/chain/convergence_criterion.cpp）。
  Scalar absTolerance = 1.0e-12;

  /// **非线性残差判据（推荐）**：不平衡力相对于重力量级的比值上限。≤0 表示禁用。
  ///
  /// 残差定义为"当前位形离本子步不动点还有多远"，单位 m/s²：
  ///     r_i = |(m_i/h²)(x_i - x̂_i) - m_i·g - f_i^int| / m_i
  /// 其中 f^int 是距离约束的弹性力。不动点处 r = 0（这是**位移判据做不到**的：
  /// 位移小只说明"这一步没怎么动"，不说明"离正确解近"）。
  ///
  /// 为什么必须要有它：高刚度下 PD 外层迭代收敛极慢，收缩因子
  /// q = κ/(κ + m/h²)。工作点 κ=1e5、m/h²=11.52 时 q ≈ 0.99988，
  /// 40 次迭代只能完成切向运动的 0.46%（见 tests/chain/convergence_criterion.cpp）。
  /// 此时"相邻迭代位移小"完全无法反映求解质量，必须直接量不平衡力。
  ///
  /// 默认 1e-3 的含义：不平衡力不超过重力的千分之一。物理上可忽略，
  /// 且比"位移判据"更贴近"这个子步解好了没有"。
  ///
  /// **注意这里有两个不同的"默认"（2026-09-21 复核时点明的口径）**：
  ///   · 库默认（本字段）：**1e-3**，面向"要精度"的无窗口运行与测试；
  ///   · 查看器默认：**3e-1**（`app/viewer/main.cpp` 显式覆盖，实时档）。
  /// 两者不是松紧微调的关系：60×60 布料的稳态残差大约只降到 0.3·g 量级，
  /// 1e-3·g 在 40 次迭代预算内**达不到**，会跑满迭代（物理 ~47 ms/帧、约 20 FPS）。
  /// 因此文档里"60 FPS / 迭代 1 次"的性能数字对应的是**查看器的 0.3 档**，
  /// 不是本默认值；引用这些数字时必须同时写明门槛（见 docs/code-map.md §3）。
  Scalar residualTolerance = 1.0e-3;

  // ---- 规模 ----
  int gridNx = 40;
  int gridNy = 40;
  Scalar gridSpacing = 0.025;
  std::string meshPath;  ///< 非空则从 OBJ 载入，忽略 gridNx/gridNy

  /// **剪切（对角）约束的刚度**：> 0 时给规则网格的每个四边形加一条对角边（见
  /// `Mesh::addShearDiagonals`），提供抗剪切。**默认 0 = 不启用** —— 此时网格、边表、
  /// 装配矩阵与之前完全一致（所有历史基线与断言都不受影响）。
  /// 注：它复用"距离约束"的整套实现（装配/着色/投影/散射/三分量并行），所以不影响
  /// "与位形无关、只分解一次"与 `L = Ã ⊗ I₃` 这两条架构性质。
  Scalar shearStiffness = 0.0;

  /// **线性（中点）弯曲约束的刚度 k**：> 0 时沿每一行/列给连续三个顶点生成弯曲 stencil
  /// （见 `Mesh::addBendingStencils` 与 `BendStencil`），提供抗弯刚度 —— 布料会自己
  /// 挺起来/保留褶皱，而不是像面条一样只靠距离约束。
  ///
  /// **默认 0 = 不启用**：此时 `mesh.bends` 为空、装配与残差都不多算任何一项，
  /// 行为与改动前**逐字相同**（已用 `build/_baseline/pd_bench_pre_fuse.exe` 逐字比对）。
  ///
  /// 两条架构性质不受影响（这是选"中点形式"而不是二面角形式的理由）：
  /// ① A_c 是常数矩阵 ⇒ L 与位形无关 ⇒ 仍然只分解一次；
  /// ② 每块是 `k·w_p·w_q·I₃` ⇒ 仍是标量 × I₃ ⇒ `L = Ã ⊗ I₃` 与三分量并行保住。
  ///
  /// ⚠️ **量纲提醒（局限，不是 bug）**：`‖A_c x‖²` 里有两个位移差分，所以这条约束的
  /// 等效抗弯刚度 **∝ s⁴ / k**（s = 网格间距）—— 弯曲刚度是**网格相关**的，
  /// 换间距/分辨率必须重新标定，而且它会把弯曲刚度与间距耦合成"同一个 k 在不同网格上
  /// 手感不同"。要尺度无关得改成 `k·s⁴`（或按质量参数化，与 κ 的处理原则一致）。
  Scalar bendStiffness = 0.0;

  /// **弯曲 stencil 的采样形状**（见 `BendSampling`）：默认 `Standard` = 行列全取、
  /// 步长 1，即"照旧"。它只决定**取哪些** stencil，不改变任何架构性质（stencil 变少 ⇒
  /// ⊗ 结构、只分解一次、三分量并行都自动成立），代价是各向异性与等效刚度变软：
  /// 采样密度降到 1/s，就要把 `bendStiffness` 大致乘 s 补回来（实测标定见 docs/perf.md §13）。
  BendSampling bendSampling = BendSampling::Standard;

  /// 弯曲 stencil 中心顶点的采样步长（>= 1；2 = 隔一个取一个）。默认 1。
  int bendStride = 1;

  // ---- 诊断 ----
  bool verbose = false;
  int reportEvery = 0;  ///< >0 时每 N 步打印一次状态
};

/// 单个阶段耗时统计（单位：秒）。
struct StageTimes {
  double predict = 0.0;
  double localStep = 0.0;
  double scatter = 0.0;
  double assemble = 0.0;   ///< 重新数值分解（只发生在 stamp 变化时）
  double solve = 0.0;      ///< 前代 + 回代（每迭代一次，热路径）
  double unpackScan = 0.0; ///< 解包 + 收敛扫描（逐顶点融合的一趟；含 barrier 等待）
  double velocity = 0.0;
  double total = 0.0;
};

struct SimContext {
  SceneConfig config;
  Mesh mesh;
  std::unique_ptr<IGlobalSolver> solver;

  uint64_t topologyId = 1;      ///< 拓扑变更时自增
  SolverStamp stamp;            ///< 上次分解时的判定戳
  bool stampValid = false;      ///< 是否已经分解过至少一次
  uint64_t structureStamp = 0;  ///< 上次符号分解时的**稀疏结构**戳（见 Assembler.h）

  int factorizeCount = 0;       ///< 实际发生的数值分解次数（测试要断言它）
  int iterationsUsed = 0;       ///< 最近一个子步用掉的 PD 迭代数
  int earlyExitCount = 0;       ///< 累计"未跑满 maxIterations 就判收敛退出"的子步数
                                ///< （诊断用：判断收敛判据是否过于宽松的关键指标）
  bool converged = false;       ///< 最近一个子步是否**达到收敛判据**（而非用尽迭代预算）
                                ///< 区分"解好了"与"迭代用完了"：后者在报告里必须写清，
                                ///< 不能当成收敛。见 docs/pd-convergence.md。
  Scalar lastResidual = -1.0;   ///< 最近一次迭代的**非线性残差**（m/s²，最大不平衡力/质量）
                                ///< 负值表示本子步未启用残差判据（residualTolerance <= 0）

  std::vector<Vec3> predicted;  ///< 预测位置 x̂
  std::vector<Vec3> positionsBeforeStep;      ///< 本子步开始时的位置（速度更新用）
  std::vector<Vec3> positionsBeforePrevStep;  ///< 上一步开始时的位置（阻尼项需要）
  std::vector<Vec3> targets;    ///< 局部步的投影目标 d_c
  Eigen::VectorXd b;            ///< 右端
  Eigen::VectorXd bBase;        ///< 惯性右端 (M/h²)x̂，每迭代从它拷贝
  Eigen::VectorXd xSolution;    ///< 全局步解（按 [x,y,z] 连续存放）
  Eigen::SparseMatrix<Scalar> L;///< 全局矩阵（数值只在 stamp 变化时重装）

  /// **线性（中点）弯曲的常向量右端**（尺寸 3N，按自由度索引）。
  /// 与位形无关（只有 stencil 里含 pinned 顶点时才非零），每子步由
  /// `assembleBendingRhs()` 重算一次（O(B)，可忽略），进区域前完成 ⇒ 无竞争。
  ///
  /// 为什么是**每子步**重算而不是"随 L 一起只算一次"：它依赖 `pinPositions`，而拖拽把手
  /// **只改 b、不改 L**（docs/plan.md §2.2.1 的约定），所以"pin 位置变了但 stamp 没变"
  /// 是常态 —— 那时若沿用旧的常向量，自由端的右端就少了 pin 的当前位置信息。
  /// 没有弯曲约束时它恒为空 vector，`data()` 为 nullptr，热路径上等于不存在。
  std::vector<Scalar> bendRhs;

  StageTimes times;
  double lastStepSeconds = 0.0;
};

/// 按当前网格规模准备所有临时缓冲（幂等）。
/// 任何构造 SimContext 的路径（makeScene / 单元测试 / 调试程序）都应调用它，
/// 否则 stepOnce 会在未分配的缓冲上越界 —— 这是实际踩过的坑，故做成显式且自保护的接口。
void ensureBuffers(SimContext& ctx);

/// 按配置构建场景：生成或载入网格、设置质量/pin/刚度、创建求解器。
SimContext makeScene(const SceneConfig& config);

/// 把网格中 pinned 标记过的顶点写入 pinPositions（构造场景时调用一次）。
void refreshPinPositions(SimContext& ctx);

}  // namespace pd
