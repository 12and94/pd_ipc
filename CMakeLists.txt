cmake_minimum_required(VERSION 3.20)
project(pd_ipc LANGUAGES CXX)

# ---------------------------------------------------------------------------
# 全局编译设置
#   统一使用 MSVC（见 docs/environment.md：PATH 上的 MinGW g++ 8.1 不支持 -std=c++20）
# ---------------------------------------------------------------------------
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif()

# ---------------------------------------------------------------------------
# 依赖：全部使用本机已安装的库，不引入包管理器（docs/environment.md §4）
# ---------------------------------------------------------------------------
set(EIGEN_ROOT "D:/tools/Eigen/eigen-3.4.0" CACHE PATH "Eigen 3 的头文件根目录")
if(NOT EXISTS "${EIGEN_ROOT}/Eigen/Dense")
  message(FATAL_ERROR "找不到 Eigen：${EIGEN_ROOT}/Eigen/Dense 不存在。请用 -DEIGEN_ROOT=<路径> 指定。")
endif()

set(GLFW_ROOT "D:/Vulkan/glfw-3.3.8.bin.WIN64" CACHE PATH "GLFW 预编译包根目录")

find_package(OpenMP QUIET)

# ---------------------------------------------------------------------------
# 核心库：不依赖窗口系统，保证 bench 与 CI 永远可用
# ---------------------------------------------------------------------------
add_library(pd_core STATIC
  core/math/Parallel.cpp
  core/mesh/Mesh.cpp
  core/energy/DistanceTerm.cpp
  core/assemble/Assembler.cpp
  core/solver/EigenDirectSolver.cpp
  core/sim/Scene.cpp
  core/sim/Integrator.cpp
)
target_include_directories(pd_core PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
# Eigen 用 SYSTEM 包含，避免它的内部警告污染构建输出
target_include_directories(pd_core SYSTEM PUBLIC "${EIGEN_ROOT}")

if(MSVC)
  # /bigobj：模板实例较多时避免 obj 段数超限
  target_compile_options(pd_core PUBLIC /W4 /permissive- /bigobj /utf-8)
  # Eigen 与 MSVC 的已知组合：定义这个宏可避免部分对齐相关的告警
  target_compile_definitions(pd_core PUBLIC EIGEN_MPL2_ONLY)
endif()

if(OpenMP_CXX_FOUND)
  target_link_libraries(pd_core PUBLIC OpenMP::OpenMP_CXX)
  message(STATUS "OpenMP 已启用：并行局部步与归约可用")
else()
  message(WARNING "未找到 OpenMP：将退化为单线程（功能正确，性能下降）")
endif()

# ---------------------------------------------------------------------------
# 测试
# ---------------------------------------------------------------------------
enable_testing()

add_executable(test_primitives tests/primitives/test_distance_term.cpp)
target_link_libraries(test_primitives PRIVATE pd_core)
add_test(NAME primitives COMMAND test_primitives)

add_executable(test_spring_vertical tests/chain/spring_vertical.cpp)
target_link_libraries(test_spring_vertical PRIVATE pd_core)
add_test(NAME spring_vertical COMMAND test_spring_vertical)

# ---------------------------------------------------------------------------
add_executable(pd_std_pd _verify/standard_pd.cpp)
target_link_libraries(pd_std_pd PRIVATE pd_core)

add_executable(pd_steady _verify/steady.cpp)
target_link_libraries(pd_steady PRIVATE pd_core)

add_executable(pd_one_step _verify/one_step_trace.cpp)
target_link_libraries(pd_one_step PRIVATE pd_core)


# ---------------------------------------------------------------------------
# 验收程序：只输出"正确/错误"（见 _verify/check.cpp）
# ---------------------------------------------------------------------------
add_executable(pd_trans _verify/translation_invariance.cpp)
target_link_libraries(pd_trans PRIVATE pd_core)
add_executable(pd_pinpos _verify/pin_positions.cpp)
target_link_libraries(pd_pinpos PRIVATE pd_core)
add_executable(pd_chain _verify/chain_test.cpp)
target_link_libraries(pd_chain PRIVATE pd_core)
add_executable(pd_springfree _verify/spring_free.cpp)
target_link_libraries(pd_springfree PRIVATE pd_core)
add_executable(pd_fvs _verify/force_vs_strain.cpp)
target_link_libraries(pd_fvs PRIVATE pd_core)
add_executable(pd_forceaudit _verify/force_audit.cpp)
target_link_libraries(pd_forceaudit PRIVATE pd_core)
add_executable(pd_initaudit _verify/init_audit.cpp)
target_link_libraries(pd_initaudit PRIVATE pd_core)
add_executable(pd_hang _verify/hang_test.cpp)
target_link_libraries(pd_hang PRIVATE pd_core)
add_executable(pd_kappafree _verify/kappa_free.cpp)
target_link_libraries(pd_kappafree PRIVATE pd_core)
add_executable(pd_kappa _verify/kappa_effect.cpp)
target_link_libraries(pd_kappa PRIVATE pd_core)
add_executable(pd_calib _verify/stiffness_calibration.cpp)
target_link_libraries(pd_calib PRIVATE pd_core)
add_executable(pd_scale _verify/scale_sweep.cpp)
target_link_libraries(pd_scale PRIVATE pd_core)
add_executable(pd_diraudit _verify/direction_audit.cpp)
target_link_libraries(pd_diraudit PRIVATE pd_core)
add_executable(pd_dampstab _verify/damping_stability.cpp)
target_link_libraries(pd_dampstab PRIVATE pd_core)
add_executable(pd_rhs _verify/rhs_breakdown.cpp)
target_link_libraries(pd_rhs PRIVATE pd_core)
add_executable(pd_iterscan _verify/iteration_scan.cpp)
target_link_libraries(pd_iterscan PRIVATE pd_core)
add_executable(pd_check _verify/check.cpp)
target_link_libraries(pd_check PRIVATE pd_core)

# ---------------------------------------------------------------------------
# 无窗口基准（CI 与性能基线用）
# ---------------------------------------------------------------------------
add_executable(pd_bench app/bench/main.cpp)
target_link_libraries(pd_bench PRIVATE pd_core)

# ---------------------------------------------------------------------------
# 可交互查看器（GLFW + OpenGL 固定管线；不引入 ImGui，HUD 用线段数字自绘）
# ---------------------------------------------------------------------------
if(EXISTS "${GLFW_ROOT}/include/GLFW/glfw3.h")
  set(GLFW_LIB "${GLFW_ROOT}/lib-vc2022/glfw3.lib")
  if(NOT EXISTS "${GLFW_LIB}")
    set(GLFW_LIB "${GLFW_ROOT}/lib-vc2022/glfw3_mt.lib")
  endif()
  if(EXISTS "${GLFW_LIB}")
    add_executable(pd_viewer app/viewer/main.cpp)
    target_link_libraries(pd_viewer PRIVATE pd_core "${GLFW_LIB}" opengl32)
    target_include_directories(pd_viewer PRIVATE "${GLFW_ROOT}/include")
    # GLFW 静态库需要链接系统库
    target_link_libraries(pd_viewer PRIVATE user32 gdi32 shell32)
    # 静态链接 GLFW 时需要此宏
    target_compile_definitions(pd_viewer PRIVATE GLFW_STATIC)
    message(STATUS "查看器已启用，GLFW: ${GLFW_LIB}")
  else()
    message(WARNING "找到 GLFW 头文件但找不到库，跳过 pd_viewer")
  endif()
else()
  message(STATUS "未找到 GLFW，跳过 pd_viewer（bench 与测试不受影响）")
endif()
