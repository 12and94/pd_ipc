// core/math/Vec3.h
// 三维向量与基础几何运算。
// 设计说明：CPU 版本统一使用 double（见 docs/plan.md 决策 D8），
// 目的是让"参考实现"具备可对照机器精度的能力；float32 只会在 GPU 后端出现。
#pragma once

#include <cmath>

namespace pd {

using Scalar = double;

struct Vec3 {
  Scalar x = 0.0, y = 0.0, z = 0.0;

  Vec3() = default;
  Vec3(Scalar x_, Scalar y_, Scalar z_) : x(x_), y(y_), z(z_) {}

  Scalar& operator[](int i) { return i == 0 ? x : (i == 1 ? y : z); }
  Scalar operator[](int i) const { return i == 0 ? x : (i == 1 ? y : z); }
};

// ---- 基本运算 ----
inline Vec3 operator+(const Vec3& a, const Vec3& b) { return {a.x + b.x, a.y + b.y, a.z + b.z}; }
inline Vec3 operator-(const Vec3& a, const Vec3& b) { return {a.x - b.x, a.y - b.y, a.z - b.z}; }
inline Vec3 operator-(const Vec3& a) { return {-a.x, -a.y, -a.z}; }
inline Vec3 operator*(const Vec3& a, Scalar s) { return {a.x * s, a.y * s, a.z * s}; }
inline Vec3 operator*(Scalar s, const Vec3& a) { return a * s; }
inline Vec3 operator/(const Vec3& a, Scalar s) { return {a.x / s, a.y / s, a.z / s}; }

inline Vec3& operator+=(Vec3& a, const Vec3& b) { a.x += b.x; a.y += b.y; a.z += b.z; return a; }
inline Vec3& operator-=(Vec3& a, const Vec3& b) { a.x -= b.x; a.y -= b.y; a.z -= b.z; return a; }
inline Vec3& operator*=(Vec3& a, Scalar s) { a.x *= s; a.y *= s; a.z *= s; return a; }

inline Scalar dot(const Vec3& a, const Vec3& b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
inline Vec3 cross(const Vec3& a, const Vec3& b) {
  return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}
inline Scalar lengthSquared(const Vec3& a) { return dot(a, a); }
inline Scalar length(const Vec3& a) { return std::sqrt(dot(a, a)); }

/// 归一化。长度过小时返回零向量（调用方负责处理退化情形）。
inline Vec3 normalized(const Vec3& a, Scalar eps = 1e-300) {
  const Scalar len = length(a);
  if (len <= eps) return Vec3{0, 0, 0};
  return a / len;
}

/// 逐分量取最大值，用于无穷范数。
inline Scalar maxAbsComponent(const Vec3& a) {
  return std::fmax(std::fabs(a.x), std::fmax(std::fabs(a.y), std::fabs(a.z)));
}

}  // namespace pd
