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
};

}  // namespace pd
