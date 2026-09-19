# 代码地图与验证方法

> 本文档只讲两件事：**每部分代码在哪**、**怎么运行并看结果**。
> 计划与设计见 `docs/plan.md`，进度与现象记录见 `README.md`。

---

## 1. 先看这一条命令

```powershell
cd D:\dsh_workspace\pd_ipc
.\tools\build.ps1
.\build\Release\pd_check.exe
```

输出就是"对/错"，逐项列出：

```
  1 投影长度 == 静止长度     OK
  2 单步解 == 闭式解           OK
  3 静止平衡 == ℓ + m g/κ     NG     实测 0.09995095  期望 0.10004905  容差 1e-07
  4 pinned 顶点严格不动        OK
  5 自由落体一步 == -h²g      OK
  6 布料有界                     OK
  7 数值分解次数 == 1          OK

  存在错误（失败 1 项）
```

退出码 0 = 全对，1 = 有错（可用于 CI）。
**当前状态：6 项 OK，第 3 项 NG** —— 这就是唯一未解决的问题：
悬挂弹簧的静止平衡落在 $\ell-mg/\kappa$ 而不是 $\ell+mg/\kappa$。

---

## 2. 代码地图（按数据流顺序）

每一部分都标了**它负责什么**和**改这里会影响哪一项验收**。

### 2.1 数学与并行原语

| 文件 | 内容 |
|---|---|
| `core/math/Vec3.h` | 三维向量与点乘/叉乘/长度。全项目统一 `double` |
| `core/math/Parallel.h` | `parallelFor`；`ThreadLocalBuffers`（线程私有缓冲 + 固定顺序归约，不用浮点原子加，保证结果可复现） |
| `core/math/Parallel.cpp` | 线程数（默认为物理核数减 2，因本机 20 逻辑核含 E-core） |

### 2.2 网格与场景

| 文件 | 内容 | 相关验收 |
|---|---|---|
| `core/mesh/Mesh.h/.cpp` | 顶点/速度/质量/pin 标记、边表（距离约束）、稀疏结构 `buildSparsityPattern()`、规则网格生成、OBJ 载入 | 6、7 |
| `core/sim/Scene.h` | `SceneConfig`（dt / 重力 / 刚度 / 密度 / 阻尼 / 迭代数 / 容差）、`SimContext`（网格 + 求解器 + 所有临时缓冲 + 阶段耗时）、`StageTimes` | — |
| `core/sim/Scene.cpp` | `makeScene()` 按配置建场景；`ensureBuffers()` 统一分配临时缓冲；`refreshPinPositions()` | 6、7 |

### 2.3 约束内核（**符号问题就在这里**）

| 文件 | 内容 | 相关验收 |
|---|---|---|
| `core/energy/DistanceTerm.h` | 接口：`project` / `projectOne` / `scatterInto` / `assembleMatrix` | 1、2、3 |
| `core/energy/DistanceTerm.cpp` | **三处与方向有关的约定都在这个文件里**：<br>① `project()` 第 ~47 行：`raw = freeIsB ? (xb - xa) : (xa - xb)`（投影方向）<br>② `scatterInto()` 第 ~78/114/152 行：`contribution = targets * (-e.stiffness)`（右端散射符号）<br>③ `assembleMatrix()` 第 ~192 行：对角 `+κ`、耦合 `-κ`（矩阵块符号） | 1、2、3 |

> **排查建议**：这三处的符号必须成对。看代码时请把它们放在一起看，
> 单看任何一处都是"对的"。第 3 项 NG 的根因就在这一组的配对上。

### 2.4 全局矩阵组装与求解器

| 文件 | 内容 | 相关验收 |
|---|---|---|
| `core/assemble/Assembler.h/.cpp` | `computeStamp()`（判定何时需要重新分解：拓扑/刚度/质量/h/pin/阻尼）、`assembleLeftHandSide()`（$M/h^2+\sum\kappa A^\top A$ + pinned 行覆盖）、`assembleInertialRhs()`（右端 $(M/h^2)\hat x$，**重力不在这里**）、`applyPinRhs()`（pinned 行右端置为把手位置）、`pack/unpackPositions` | 3、7 |
| `core/solver/IGlobalSolver.h` | 求解器抽象：`analyze` / `factorize` / `solve` 三段 + `SolverStamp`。三条路线的差异都收在这里 | 7 |
| `core/solver/EigenDirectSolver.h/.cpp` | `Eigen::SimplicialLDLT` 实现；符号分解用**真实矩阵结构**（不能用理想块模式，否则 solve 阶段会访问越界） | 7 |
| `core/solver/SparsePattern.h` | `BlockEntry`（3×3 非零块的位置） | — |

### 2.5 时间积分（**阻尼与速度更新在这里**）

| 文件 | 内容 | 相关验收 |
|---|---|---|
| `core/sim/Integrator.h` | `stepOnce` / `stepFrame` / `surrogateEnergy` / `totalEnergy` / `maxSpeed` | 2、3、5、6 |
| `core/sim/Integrator.cpp` | `stepOnce()` 的五个阶段，逐个编号：<br>**1) 预测** $\hat x=x+hv+h^2g$（第 ~60 行）——**重力只在这里出现一次**<br>**2) 判定+数值分解**（第 ~85 行）——只在 stamp 变化时执行<br>**3) 惯性右端**（第 ~105 行）<br>**4) PD 迭代**：局部步 → 散射 → 覆盖 pin → 全局回代 → 收敛判据（第 ~115–158 行）<br>**5) 速度更新与阻尼**（第 ~160–190 行） | 2、3、5、6 |

> **排查建议**：第 3 项 NG 与第 5 阶段的速度更新形式直接相关。
> 目前形式是 $v_{n+1}=(1-k_d)(x_{n+1}-x_n)/h$。
> 把 $\hat y=y-h^2g$ 与行方程联立可得 $\kappa(y-\ell)=-mg$（错误）；
> 要让它是 $+mg$，第 5 阶段或右端二者之一需要改。

### 2.6 应用

| 文件 | 内容 |
|---|---|
| `app/bench/main.cpp` | 无窗口基准：分阶段耗时、迭代数、应变、**分解复用统计**（符号/数值分解与回代次数） |
| `app/viewer/main.cpp` | GLFW + OpenGL 固定管线查看器；HUD 用线段自绘数字（不依赖 ImGui/字体）；键盘可调刚度/迭代/子步/线程 |

### 2.7 验证程序（`_verify/`）

| 文件 | 用途 |
|---|---|
| `_verify/check.cpp` | **唯一验收程序**：7 项检查，只输出 OK/NG |
| `_verify/one_step_trace.cpp` | 逐项打印一步的中间量（预测/投影/散射/对角/解），与解析式并排对照。**排查符号问题的主要工具** |
| `_verify/standard_pd.cpp` | 用标准 PD 公式**独立组装并求解**，与代码组装逐项比对 |
| `_verify/steady.cpp` | 多步稳态：从候选平衡出发是否保持不动 |
| `_verify/fixed_point.cpp` | 枚举"矩阵块符号 × 散射符号"，看哪个让不动点满足物理力平衡 |
| `_verify/variants.cpp` | 枚举"投影方向 × 矩阵块形式"，同时看平衡与布料稳定性 |

这些程序原先只在需要时手工构建，现在 `pd_check` 已纳入常规构建；
其余几个如需运行，在 `CMakeLists.txt` 里给它们各加两行即可（照抄 `pd_check` 的写法）。

---

## 3. 运行与验证命令

```powershell
# 构建（Release）
.\tools\build.ps1
.\tools\build.ps1 -Clean      # 清理后重建

# ---- 验收（只看对错）----
.\build\Release\pd_check.exe          # 8 项，退出码 0/1
#   1 投影长度 == 静止长度        5 自由落体一步 == -h²g
#   2 单步解 == 闭式解            6 布料有界
#   3 静止平衡 == ℓ - m g/κ       7 数值分解次数 == 1
#   4 pinned 顶点严格不动         8 弹性力方向 == -dU/dx

# ---- 测试套件（打印每条断言）----
.\build\Release\test_primitives.exe         # 45 断言
.\build\Release\test_spring_vertical.exe    # 124 断言

# ---- 专项验证 ----
.\build\Release\pd_diraudit.exe       # 投影方向：全局一致性 / 物理合理性 / 标准 PD 一致性
.\build\Release\pd_dampstab.exe       # 小阻尼 + 长时间稳定性（单弹簧 + 布料）
.\build\Release\pd_scale.exe          # 网格规模扫描：稳定性、耗时、分解复用
.\build\Release\pd_calib.exe          # 刚度标定扫描
.\build\Release\pd_kappafree.exe      # 无 pin/无重力下 κ 是否控制应变
.\build\Release\pd_std_pd.exe         # 代码组装 vs 标准 PD 独立求解
.\build\Release\pd_rhs.exe            # 右端 b 的逐项组成
.\build\Release\pd_one_step.exe       # 一步的每个中间量 vs 解析式
.\build\Release\pd_iterscan.exe       # 迭代次数扫描（看是否收敛/震荡）

# ---- 性能 ----
.\build\Release\pd_bench.exe --grid 40 40 --steps 600 --iters 10

# ---- 实时渲染（GLFW + OpenGL，可截图）----
.\build\Release\pd_viewer.exe                                   # 交互窗口
.\build\Release\pd_viewer.exe --grid 60 --frames 300 --shot out.png   # 有界运行 + 截图
```

查看器参数：`--grid N` 网格边长、`--stiffness K` 刚度、`--frames N` 跑够帧数自动退出、
`--shot FILE.png` 截图（自写 PNG 编码，无图像库依赖）、`--pin-corners` 只钉两角（默认钉整条上边）。
窗口内交互：左键拖拽旋转 / 滚轮缩放 / `SPACE` 暂停 / `S` 单步 / `G` 重力 / `R` 重置 /
`[` `]` 刚度 / `-` `=` 迭代数 / `,` `.` 子步 / `;` `'` 线程数 / `ESC` 退出。

**实测（RTX 4070 Ti 主机，CPU 18 线程）**：60×60 布料（3600 顶点 / 7080 约束），
vsync 下稳定 **60 FPS**，物理 **1.6–2.1 ms/帧**（PD 迭代 1–2 次即达收敛判据）。

调试环境变量：

```powershell
$env:PD_TRACE_STEP = "1"      # stepOnce 逐步打印：预测 / 迭代残差 / 求解结果
$env:PD_CHECK_VERBOSE = "1"   # pd_check 第 8 项打印逐点明细
```

---

## 4. 当前状态一句话

投影、矩阵装配、pin 处理、自由落体、布料稳定性、预分解复用**都正确**；
**唯一 NG 是第 3 项**：悬挂弹簧的静止平衡位置（$\ell-mg/\kappa$ vs 正确的 $\ell+mg/\kappa$）。
相关代码集中在 `core/energy/DistanceTerm.cpp` 的三处方向约定
与 `core/sim/Integrator.cpp` 第 5 阶段的速度更新。
