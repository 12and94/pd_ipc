// tests/primitives/test_harness.h
// 极简测试框架：不引入第三方依赖（见 docs/plan.md 决策 D15）。
//
// 用法：
//   TEST(名字) { CHECK(条件); CHECK_NEAR(实际值, 期望值, 容差); }
//   int main 由 TEST_MAIN() 生成。
//
// 设计意图：断言必须"打印出实际值"，这样失败时不需要再去加打印。
#pragma once

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <functional>
#include <string>
#include <vector>

namespace pdtest {

struct TestCase {
  std::string name;
  std::function<void()> fn;
};

inline std::vector<TestCase>& registry() {
  static std::vector<TestCase> r;
  return r;
}

inline int& failureCount() {
  static int f = 0;
  return f;
}

inline int& checkCount() {
  static int c = 0;
  return c;
}

struct Registrar {
  Registrar(const std::string& name, std::function<void()> fn) {
    registry().push_back(TestCase{name, std::move(fn)});
  }
};

inline void reportFailure(const char* file, int line, const std::string& what) {
  ++failureCount();
  std::printf("    [FAIL] %s:%d  %s\n", file, line, what.c_str());
}

inline int runAll(const char* suiteName) {
  std::printf("=== %s: %zu 个测试 ===\n", suiteName, registry().size());
  int failedTests = 0;
  for (auto& t : registry()) {
    const int before = failureCount();
    std::printf("  [run ] %s\n", t.name.c_str());
    t.fn();
    if (failureCount() > before) {
      ++failedTests;
      std::printf("  [BAD ] %s\n", t.name.c_str());
    } else {
      std::printf("  [ok  ] %s\n", t.name.c_str());
    }
  }
  std::printf("--- %s: %zu 个断言, %d 个失败测试, %d 个失败断言 ---\n", suiteName, checkCount(),
              failedTests, failureCount());
  return failedTests == 0 ? 0 : 1;
}

}  // namespace pdtest

#define TEST(name)                                                        \
  static void test_##name();                                              \
  static ::pdtest::Registrar registrar_##name(#name, test_##name);        \
  static void test_##name()

#define CHECK(cond)                                                                    \
  do {                                                                                 \
    ++::pdtest::checkCount();                                                          \
    if (!(cond)) ::pdtest::reportFailure(__FILE__, __LINE__, "CHECK 失败: " #cond);      \
  } while (0)

#define CHECK_NEAR(actual, expected, tol)                                                        \
  do {                                                                                           \
    ++::pdtest::checkCount();                                                                    \
    const double a_ = static_cast<double>(actual);                                               \
    const double e_ = static_cast<double>(expected);                                              \
    const double t_ = static_cast<double>(tol);                                                   \
    if (!(std::fabs(a_ - e_) <= t_)) {                                                            \
      char buf_[512];                                                                             \
      std::snprintf(buf_, sizeof(buf_), "CHECK_NEAR 失败: %s = %.17g, 期望 %.17g, 差 %.3g (容差 %.3g)", \
                    #actual, a_, e_, std::fabs(a_ - e_), t_);                                     \
      ::pdtest::reportFailure(__FILE__, __LINE__, buf_);                                          \
    }                                                                                             \
  } while (0)

#define CHECK_LE(actual, bound)                                                     \
  do {                                                                              \
    ++::pdtest::checkCount();                                                       \
    const double a_ = static_cast<double>(actual);                                  \
    const double b_ = static_cast<double>(bound);                                   \
    if (!(a_ <= b_)) {                                                              \
      char buf_[512];                                                               \
      std::snprintf(buf_, sizeof(buf_), "CHECK_LE 失败: %s = %.17g > %.17g", #actual, a_, b_); \
      ::pdtest::reportFailure(__FILE__, __LINE__, buf_);                            \
    }                                                                               \
  } while (0)

#define CHECK_MSG(cond, msg)                                                      \
  do {                                                                            \
    ++::pdtest::checkCount();                                                     \
    if (!(cond)) ::pdtest::reportFailure(__FILE__, __LINE__, std::string("CHECK 失败: ") + (msg)); \
  } while (0)

#define TEST_MAIN(suiteName)                            \
  int main() { return ::pdtest::runAll(suiteName); }
