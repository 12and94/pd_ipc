// _verify/direction_audit.cpp
// 投影方向与弹性力的三项审计。全部按正确公式重写（早期版本把左端算子当力用，结论是错的）。
//
// 判据：
//   A. 全局一致性：方向只由边的端点顺序 (a,b) 决定，没有任何特例；pin 配置不影响它
//   B. 物理合理性：行方程反推的"净弹性项" == 能量梯度的负值 -dU/dx
//   C. 标准 PD 一致性：与独立写出的 A_c x = x_a - x_b、d_c = ℓ·unit(A_c x)、
//      b_a += κ d_c、b_b -= κ d_c 逐位对照
//
// 净弹性项的求法（关键，早期版本在这里出错）：
//     自由顶点行：  (m/h²)·y + κ·y = (m/h²)·ŷ + 散射
//     移项：        m(y-ŷ)/h² = 散射 − κ·y
//     ⇒ 净弹性项 = 散射 − κ·y        （κ·y 是左端算子，必须先移项再比较）
#include <cmath>
#include <cstdio>
#include <vector>

#include "core/energy/DistanceTerm.h"
#include "core/mesh/Mesh.h"

using namespace pd;

namespace {

int gFail = 0;

void item(const char* name, bool ok, const char* detail = "") {
  std::printf("  %-56s %s %s\n", name, ok ? "OK" : "NG", detail);
  if (!ok) ++gFail;
}

const Scalar kEll = 0.1;
const Scalar kK = 1.0e4;

/// 独立标准 PD：返回顶点对 (a,b) 上的散射贡献（a 收 +κ d_c，b 收 -κ d_c）
struct StandardPd {
  Vec3 dcOnA;
  Vec3 dcOnB;
};

StandardPd standardScatter(const Vec3& xa, const Vec3& xb, Scalar ell, Scalar k) {
  const Vec3 Acx = xa - xb;
  const Scalar r = length(Acx);
  StandardPd out;
  if (r < 1e-15) return out;
  const Vec3 dc = Acx * (ell / r);
  out.dcOnA = dc * k;
  out.dcOnB = dc * (-k);
  return out;
}

/// 真值：距离约束的弹性力（能量梯度负值，方向取"受力点减施力点"）
Vec3 trueForce(const Vec3& xFree, const Vec3& xOther, Scalar ell, Scalar k) {
  const Vec3 d = xFree - xOther;
  const Scalar r = length(d);
  if (r < 1e-15) return Vec3{0, 0, 0};
  return d * (-k * (r - ell) / r);
}

}  // namespace

int main() {
  // ---------- A. 全局一致性 ----------
  std::printf("A. 全局一致性（方向只由端点顺序决定，无特例）\n");
  {
    Mesh m1;  // a 在原点、b 在 +y
    m1.positions = {Vec3{0, 0, 0}, Vec3{0, 0.15, 0}};
    m1.pinned = {0, 0};
    m1.edges.push_back(Edge{0, 1, kEll, kK});
    std::vector<Vec3> t1;
    DistanceTerm::project(m1, t1);

    Mesh m2;  // 同一几何，端点顺序相反
    m2.positions = {Vec3{0, 0.15, 0}, Vec3{0, 0, 0}};
    m2.pinned = {0, 0};
    m2.edges.push_back(Edge{0, 1, kEll, kK});
    std::vector<Vec3> t2;
    DistanceTerm::project(m2, t2);

    item("d_c 沿 x_a - x_b（顺序 a=原点,b=上方 ⇒ 沿 -y）",
         std::fabs(t1[0].y + kEll) < 1e-14);
    item("交换端点顺序 ⇒ d_c 整体反号", std::fabs(t2[0].y - kEll) < 1e-14);
    item("长度恒为 ℓ（与顺序无关）",
         std::fabs(length(t1[0]) - kEll) < 1e-15 && std::fabs(length(t2[0]) - kEll) < 1e-15);

    // pin 配置不应影响方向
    std::vector<Vec3> tFree, tPinA, tPinB;
    m1.pinned = {0, 0};
    DistanceTerm::project(m1, tFree);
    m1.pinned = {1, 0};
    DistanceTerm::project(m1, tPinA);
    m1.pinned = {0, 1};
    DistanceTerm::project(m1, tPinB);
    item("三种 pin 配置下 d_c 完全相同（无按 pin 翻转的特例）",
         tPinA[0].y == tFree[0].y && tPinB[0].y == tFree[0].y);
  }

  // ---------- B. 物理合理性（净弹性项 vs 能量梯度）----------
  std::printf("\nB. 物理合理性（净弹性项 == -dU/dx）\n");
  {
    bool allOk = true;
    double worstGot = 0, worstWant = 0;
    const Scalar ys[] = {0.06, 0.08, 0.0999, 0.1001, 0.12, 0.15};
    for (Scalar y : ys) {
      // 自由端是 b
      {
        Mesh m;
        m.positions = {Vec3{0, 0, 0}, Vec3{0, y, 0}};
        m.pinned = {1, 0};
        m.edges.push_back(Edge{0, 1, kEll, kK});
        std::vector<Vec3> t;
        DistanceTerm::project(m, t);
        std::vector<Scalar> b(6, 0.0);
        DistanceTerm::scatterInto(m, t, b.data(), b.size());
        const Vec3 xFree = m.positions[1];
        const Vec3 scatter{b[3], b[4], b[5]};
        const Vec3 net = scatter - xFree * kK;              // 净弹性项 = 散射 − κy
        const Vec3 want = trueForce(xFree, m.positions[0], kEll, kK);
        if (std::fabs(net.y - want.y) > 1e-6) { allOk = false; worstGot = net.y; worstWant = want.y; }
      }
      // 自由端是 a
      {
        Mesh m;
        m.positions = {Vec3{0, y, 0}, Vec3{0, 0, 0}};
        m.pinned = {0, 1};
        m.edges.push_back(Edge{0, 1, kEll, kK});
        std::vector<Vec3> t;
        DistanceTerm::project(m, t);
        std::vector<Scalar> b(6, 0.0);
        DistanceTerm::scatterInto(m, t, b.data(), b.size());
        const Vec3 xFree = m.positions[0];
        const Vec3 scatter{b[0], b[1], b[2]};
        const Vec3 net = scatter - xFree * kK;
        const Vec3 want = trueForce(xFree, m.positions[1], kEll, kK);
        if (std::fabs(net.y - want.y) > 1e-6) { allOk = false; worstGot = net.y; worstWant = want.y; }
      }
    }
    char buf[128];
    std::snprintf(buf, sizeof(buf), "(最后不符: got=%+.4g want=%+.4g)", worstGot, worstWant);
    item("拉伸/压缩两侧、两种端点顺序，净弹性项均等于 -dU/dx", allOk, allOk ? "" : buf);
  }

  // ---------- C. 标准 PD 一致性 ----------
  std::printf("\nC. 标准 PD 一致性（与独立实现逐位对照）\n");
  {
    bool same = true;
    for (int trial = 0; trial < 8; ++trial) {
      const Vec3 xa{0.01 * trial, 0.02 + 0.03 * trial, -0.015 * trial};
      const Vec3 xb{-0.005 * trial, 0.0, 0.01 * trial};
      Mesh m;
      m.positions = {xa, xb};
      m.pinned = {0, 0};
      m.edges.push_back(Edge{0, 1, kEll, kK});
      std::vector<Vec3> t;
      DistanceTerm::project(m, t);
      std::vector<Scalar> b(6, 0.0);
      DistanceTerm::scatterInto(m, t, b.data(), b.size());

      const StandardPd ref = standardScatter(xa, xb, kEll, kK);
      const Vec3 gotA{b[0], b[1], b[2]};
      const Vec3 gotB{b[3], b[4], b[5]};
      if (length(gotA - ref.dcOnA) > 1e-9 || length(gotB - ref.dcOnB) > 1e-9) same = false;
    }
    item("8 组任意位形下散射与独立标准 PD 实现一致", same);
  }

  std::printf("\n  %s（失败 %d 项）\n", gFail == 0 ? "三项判据全部通过" : "存在不通过项", gFail);
  return gFail == 0 ? 0 : 1;
}
