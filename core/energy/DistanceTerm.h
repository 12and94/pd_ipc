// core/energy/DistanceTerm.h
// 首期唯一的能量项：距离约束（弹簧）。
//
// 能量形式（docs/plan.md §2.2.1）：
//   A_c x = x_a - x_b
//   E_c(x) = (κ_c / 2) * ( ‖x_a - x_b‖ - l_c )²
//
// 局部步（近端算子）有闭式解：
//   d_c = Π_c(A_c x) = l_c * (x_a - x_b) / ‖x_a - x_b‖
//
// 全局步的矩阵块是**常数**（与位形无关）：
//   L_aa += κ_c I,  L_bb += κ_c I,  L_ab -= κ_c I,  L_ba -= κ_c I
// 右端的散射量是 κ_c d_c：b_a += κ_c d_c，b_b -= κ_c d_c。
// （矩阵块与散射必须同源同一个 A_c；两处各写一份 κ_c，没有额外的系数 2。）
#pragma once

#include <vector>

#include <Eigen/Sparse>

#include "core/math/Vec3.h"
#include "core/mesh/Mesh.h"

namespace pd {

class DistanceTerm {
 public:
  /// 退化保护阈值：边长小于该值时跳过投影（避免除零）。
  static constexpr Scalar kMinLength = 1e-12;

  /// 局部步：对每条约束求投影目标 d_c。
  /// 输入 positions，输出 targets（与 mesh.edges 一一对应）。
  static void project(const Mesh& mesh, std::vector<Vec3>& targets);

  /// 单条约束的投影（供原语测试直接调用）。
  static Vec3 projectOne(const Vec3& pa, const Vec3& pb, Scalar restLength);

  /// 把 κ_c d_c 散射进右端（含 pin 消元补偿，见 .cpp 内的推导）。
  /// b 指向长度 3*vertexCount 的连续缓冲（按 [x,y,z] 存放）。
  /// 用裸指针而不是 Eigen 类型，是为了让 GPU 后端也能共用同一份语义。
  /// 内部使用"线程私有缓冲 + 固定顺序归约"，不做浮点原子加（决策 D6）。
  ///
  /// 前置条件：若 mesh 里有 pinned 顶点，`mesh.pinPositions` 必须已按顶点数填充
  /// （正常流程由 makeScene / refreshPinPositions 保证）。违反时本函数会立即
  /// 报错退出，而不是靠越界读把问题变成难以定位的崩溃。
  static void scatterInto(const Mesh& mesh, const std::vector<Vec3>& targets, Scalar* b,
                          std::size_t dim);

  /// 组装全局矩阵的数值部分（不含 M/h² 与 pinned 行——那些由组装器统一处理）。
  /// 调用方需保证 triplets 已清空；本函数只做追加。
  static void assembleMatrix(const Mesh& mesh, const std::vector<Vec3>& /*unused*/,
                             std::vector<Eigen::Triplet<Scalar>>& triplets);

  // ---------------------------------------------------------------------------
  // 区域内版本（"一个子步一个区域"方案的接口，见 docs/parallel-refactor.md §4.1）
  //
  // **契约：必须在已有的 `#pragma omp parallel` 区域内调用。** 它们只做 `omp for`，
  // 不再自己开区域 —— 这正是把每子步 2K 次 fork/join 降成 1 次的关键。
  // 在区域外调用属于误用（OpenMP 的孤儿 worksharing 构造行为未定义），
  // 因此这里不设运行时断言（串行化的区域 omp_in_parallel() 也返回 false，
  // 设断言会把合法用法误杀）；正确用法由 `Integrator::stepOnce` 与下面的
  // 独立入口各自保证：独立入口自己开区域，区域内版本由 stepOnce 在区域内调用。
  //
  // 为什么保留独立入口（project / scatterInto）：测试、`_verify/` 程序与诊断
  // 工具都直接调它们，且小网格上不该为一次投影付 fork/join 的代价。
  // ---------------------------------------------------------------------------

  /// 区域内版本：把投影写进 targets（逐边写，无跨线程依赖）。
  static void projectInRegion(const Mesh& mesh, std::vector<Vec3>& targets);

  /// 区域内版本：`b[i] = base[i] + Σ_c κ_c (A_c^T d_c)[i]`（含 pin 消元补偿）。
  ///
  /// 实现方式（Phase 2 之后）：按**约束着色**逐色执行 —— 同色边两两不共享顶点，
  /// 因此并发写 b 无冲突，**不需要私有缓冲、不需要原子加、也不需要归约**。
  /// 求和顺序是全项目唯一的那一条：**颜色序 0,1,…,C-1**（同一顶点每色最多被写一次），
  /// 所以结果与线程数、分块、调度完全无关，可位级复现。
  ///
  /// 代价：每次调用有 C 个 barrier（C = 颜色数）；换掉的是旧方案 O(P·dim) 的归约流量。
  /// 详见 `docs/parallel-refactor.md` Phase 2 与 `docs/perf.md`。
  ///
  /// base 与 b 必须都能访问 3*vertexCount 个元素；`base == b` 表示"累加到现有值"
  /// （独立入口 `scatterInto` 复用这条实现时用）。
  static void scatterIntoInRegion(const Mesh& mesh, const std::vector<Vec3>& targets,
                                  const Scalar* base, Scalar* b, std::size_t dim);
};

}  // namespace pd
