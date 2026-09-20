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
  Scalar stiffness = 1.0e4;  ///< 距离约束刚度 κ（默认值来自刚度标定，见 docs/perf.md）
  Scalar density = 1.0;      ///< 面密度（决定顶点质量）
  Scalar velocityDamping = 0.02;  ///< 速度阻尼系数 k_d

  // ---- PD 迭代 ----
  int maxIterations = 10;          ///< 每子步最大 PD 迭代数
  Scalar relTolerance = 1.0e-3;    ///< 相对收敛判据
  Scalar absTolerance = 1.0e-12;   ///< 绝对收敛判据

  // ---- 规模 ----
  int gridNx = 40;
  int gridNy = 40;
  Scalar gridSpacing = 0.025;
  std::string meshPath;  ///< 非空则从 OBJ 载入，忽略 gridNx/gridNy

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

  std::vector<Vec3> predicted;  ///< 预测位置 x̂
  std::vector<Vec3> positionsBeforeStep;      ///< 本子步开始时的位置（速度更新用）
  std::vector<Vec3> positionsBeforePrevStep;  ///< 上一步开始时的位置（阻尼项需要）
  std::vector<Vec3> targets;    ///< 局部步的投影目标 d_c
  Eigen::VectorXd b;            ///< 右端
  Eigen::VectorXd bBase;        ///< 惯性右端 (M/h²)x̂，每迭代从它拷贝
  Eigen::VectorXd xSolution;    ///< 全局步解（按 [x,y,z] 连续存放）
  Eigen::SparseMatrix<Scalar> L;///< 全局矩阵（数值只在 stamp 变化时重装）

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
