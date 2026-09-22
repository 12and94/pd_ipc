// core/solver/IGlobalSolver.h
// 全局步求解器的抽象。三条实现路线（CPU 直接法 / GPU Chebyshev / GPU 局部+CPU 全局）
// 的差异全部集中在这里（见 docs/plan.md §1.1 的架构约定）。
//
// 三个操作粒度必须分开（docs/plan.md §2.4(2)）：
//   analyze   —— 符号分解与排序，只依赖稀疏结构（拓扑）=> 通常只做一次
//   factorize —— 数值分解，依赖矩阵数值（拓扑/刚度/h/pin）=> 只在 stamp 变化时做
//   solve     —— 前代 + 回代 => 每次全局步都做，是热路径
#pragma once

#include <vector>

#include <Eigen/Sparse>

#include "core/math/Vec3.h"
#include "core/solver/SparsePattern.h"

namespace pd {

/// 触发重新数值分解的条件集合。
/// L 的数值只依赖这些量（与位形无关），因此可以用"逐字段精确比较"判断是否需要重分解：
/// 相比浮点容差比较，精确比较不会因为"几乎没变"而误判，语义也更明确。
///
/// **判据原则：只放"真的会改变 L 的数值"的量。**
/// 多放一项的代价是每次变动都白做一次数值分解（拖动滑块时会卡帧）；
/// 少放一项的代价是矩阵过期、解错误且不报错。因此每个字段都要能说清它改的是
/// L 的哪一部分（见 Assembler.cpp `computeStamp` 的逐项说明）。
struct SolverStamp {
  uint64_t topologyId = 0;      ///< 拓扑标识（顶点数、边数、结构版本）→ 决定稀疏结构
  uint64_t stiffnessBits = 0;   ///< 刚度集合的位模式哈希（进入约束块 κ_c A_cᵀA_c）
  uint64_t massBits = 0;        ///< 质量集合的位模式哈希（进入惯性对角 M/h²）
  uint64_t dtBits = 0;          ///< 子步长 h 的位模式（h 进入 M/h²）
  uint64_t pinBits = 0;         ///< pin **掩码**哈希（决定哪些行被覆盖、哪些约束块被跳过）
  uint64_t dampingBits = 0;     ///< 阻尼系数的位模式。**仅为诊断记录，不参与 == 比较**：
                                ///< 阻尼只影响速度更新（右端/时间推进），不影响 L。

  bool operator==(const SolverStamp& o) const {
    // 注意：dampingBits 刻意不参与比较（它不影响 L）。
    return topologyId == o.topologyId && stiffnessBits == o.stiffnessBits &&
           massBits == o.massBits && dtBits == o.dtBits && pinBits == o.pinBits;
  }
  bool operator!=(const SolverStamp& o) const { return !(*this == o); }
};

class IGlobalSolver {
 public:
  virtual ~IGlobalSolver() = default;

  /// 符号分解：给定自由度数与**真实矩阵的结构**（而非理想化的块模式）。
  ///
  /// 为什么必须用真实矩阵：理想块模式（对角块 + 每条约束两个耦合块）会把
  /// "结构上恒为零"的位置也算进去（例如某顶点只与 pinned 顶点相连时，它的
  /// 非对角块会被跳过）。Eigen 的 SimplicialLDLT 在符号分解阶段就把这些位置
  /// 计入消去树，随后 factorize 一个结构不同的矩阵会让它在 solve 阶段访问
  /// 未定义的数据（实测表现为访问违例）。
  /// 用真实矩阵的结构做符号分解可确保两者严格一致。
  virtual void analyze(int n, const Eigen::SparseMatrix<Scalar>& patternMatrix) = 0;

  /// 数值分解：给定已组装好的 L。
  virtual void factorize(const Eigen::SparseMatrix<Scalar>& L) = 0;

  /// 是否已经完成符号分解（决定是否需要先调 analyze）。
  virtual bool analyzed() const = 0;

  /// 求解 L x = b。支持多右端（Eigen 的 solve(B)）以便将来批处理。
  virtual void solve(const Eigen::VectorXd& b, Eigen::VectorXd& x) = 0;

  /// **可拆分的独立分量个数**（1 = 只能整趟求解，调用方退回单线程）。
  ///
  /// 为什么会有多个（2026-09-22 的新发现，实测见 docs/perf.md §9）：
  /// 本项目 L 的每一块都是**标量 × I₃** —— 惯性项 (m_v/h²) 出现在 (3v+d, 3v+d)、
  /// 距离约束块 ±κ 出现在 (3a+d, 3b+d)（d 不变）、pin 行 (3v+d, 3v+d) = 1
  /// （见 Assembler.cpp 与 DistanceTerm::assembleMatrix）。也就是 L = Ã ⊗ I₃：
  /// 3n 个方程其实是 3 个**互不相连**的标量系统（x/y/z 完全解耦），而且消元产生的
  /// 填充不可能跨连通分量 ⇒ 因子是 3 份结构同构的标量因子。
  ///   · 于是"一次 3n 三角求解"可以拆成 3 条**互不相干**的链：零同步、零竞争；
  ///   · 每个分量的算术序列与整趟完全一致（列顺序、i>j 筛选、累加顺序都不变）
  ///     ⇒ **结果逐位相同**（`pd_solvecomp` 与测试里的断言都钉住这条性质）；
  ///   · 回代是延迟受限的（达成带宽只有流式读写的 9–13 %，见 docs/solver-feasibility.md §2），
  ///     瓶颈是"未决访存请求数"不足 —— 3 条独立链分给 3 条线程正好把它提高约 3 倍。
  ///     实测单次回代快 **2.0–2.8×**（40×40 / 100×100 / 200×200，见 docs/perf.md §9）。
  ///
  /// 结构不满足时（自由度不是 3 的倍数、或填充跨分量 —— 将来加入各向异性/非 ⊗ 的
  /// 能量项就会这样）返回 1，调用方自动退回整趟求解。
  virtual int parallelComponents() const { return 1; }

  /// 求解第 c 个分量（c ∈ [0, parallelComponents())）：只读写属于该分量的自由度。
  ///
  /// 契约：
  ///   · 不同 c 之间**内存完全不相交** ⇒ 可从多条线程并发调用（同步由调用方负责）；
  ///   · 调用前 x 必须已经是长度 n 的向量（拆分路径里不允许 resize —— 多线程下不安全）；
  ///   · x 的其它分量保持不动；
  ///   · parallelComponents() == 1 时等价于 solve()（此时允许 x 为空、自动 resize）。
  virtual void solveComponent(int c, const Eigen::VectorXd& b, Eigen::VectorXd& x) {
    (void)c;
    solve(b, x);
  }

  /// 记一次"全局步阶段"的耗时（秒）。
  ///
  /// **为什么需要它**：拆成多条线程之后，`solveComponent()` 会被并发调用，
  /// 如果让它自己 `stats_.totalSolveSeconds += ...`，那就是多线程非原子读-改-写
  /// （数据竞争，C++ 层面是 UB —— 2026-09-22 本轮真的这么写错过了，现已改掉）。
  /// 所以约定：**调用方在 `parallelComponents() > 1` 时，于该阶段的 barrier 之后
  /// 由单线程调用一次本函数**，把整个阶段的墙钟交给求解器记账。
  /// · `parallelComponents() == 1` 时本函数应为无操作（那条路上 `solve()` 自己已经记过，
  ///   不做这个判断就会重复计数）；
  /// · 默认实现是无操作，其它后端可以只依赖 `solve()` 自己的记账。
  /// 记账后 `Stats::totalSolveSeconds` 的语义是**阶段墙钟的累计**（拆分时即"每个全局步的
  /// 求解阶段用时之和"，与 `StageTimes::solve` 同口径），不再是"各分量工作量之和"。
  virtual void noteSolveStage(double seconds) { (void)seconds; }

  /// 诊断：符号分解 / 数值分解 / 回代的累计次数与最近一次耗时（秒）。
  struct Stats {
    int analyzeCalls = 0;
    int factorizeCalls = 0;
    int solveCalls = 0;
    double lastAnalyzeSeconds = 0.0;
    double lastFactorizeSeconds = 0.0;
    double lastSolveSeconds = 0.0;
    double totalFactorizeSeconds = 0.0;
    double totalSolveSeconds = 0.0;
    long long nnz = 0;         ///< analyze 时传入的结构矩阵的 nnz（方便核对"结构是多少"）
    long long factorNnz = 0;   ///< 因子的 nnz（诊断：填充有多大 —— 代价全在这里）
  };
  virtual const Stats& stats() const = 0;
};

}  // namespace pd
