// _verify/force_vs_strain.cpp
// 逐条约束对比：散射项 vs 真实弹力 κ(r-ℓ)
//   真实弹力大小 = κ·|r - ℓ|
//   散射项在顶点上的力 = |κ·d_c| = κ·ℓ      （d_c 长度恒为 ℓ）
// 若两者比值恒为 ℓ/|r-ℓ|，说明散射项没有携带应变信息。
#include <cmath>
#include <cstdio>
#include <vector>

#include "core/energy/DistanceTerm.h"
#include "core/mesh/Mesh.h"

using namespace pd;

int main() {
  const Scalar ell = 0.1, k = 1.0e4;
  std::printf("=== 单约束：散射项 vs 真实弹力 ===\n");
  std::printf("κ=%.0f, ℓ=%.3f\n\n", k, ell);
  std::printf("  %-10s %-12s %-14s %-14s %-12s %s\n", "当前长度r", "应变(r-ℓ)/ℓ", "真实弹力κ(r-ℓ)",
              "散射项|κd_c|", "比值", "说明");

  for (Scalar r : {0.09, 0.099, 0.1, 0.105, 0.11, 0.15, 0.2}) {
    Mesh m;
    m.positions = {Vec3{0, 0, 0}, Vec3{0, r, 0}};
    m.pinned = {0, 0};
    m.edges.push_back(Edge{0, 1, ell, k});
    std::vector<Vec3> t;
    DistanceTerm::project(m, t);

    const Scalar trueForce = k * (r - ell);      // 真实弹力（带符号）
    const Scalar scatter = k * length(t[0]);     // 散射项大小 = κℓ
    const Scalar ratio = (std::fabs(trueForce) > 1e-12) ? scatter / std::fabs(trueForce) : 0.0;
    std::printf("  %-10.4f %-12.4f %-14.4f %-14.4f %-12.4g %s\n", r, (r - ell) / ell, trueForce,
                scatter, ratio, ratio > 1.001 ? "散射项偏大" : (ratio < 0.999 ? "散射项偏小" : "相等"));
  }

  std::printf("\n结论要点：\n");
  std::printf("  · 真实弹力 ∝ (r-ℓ)，在 r=ℓ 处为 0；\n");
  std::printf("  · 散射项 |κ·d_c| = κ·ℓ 恒为常数（%.1f），与 r 无关。\n", k * ell);
  std::printf("  => 散射项只在 |r-ℓ| 恰好等于 ℓ×某值时与真实弹力同量级，其余情况都偏差很大。\n");
  return 0;
}
