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
  3 静止平衡 == ℓ - m g/κ     OK
  4 pinned 顶点严格不动        OK
  5 自由落体一步 == -h²g      OK
  6 布料有界                     OK
  7 数值分解次数 == 1          OK
  8 弹性力符号 == -dU/dy        OK

  全部正确（失败 0 项）
```

退出码 0 = 全对，1 = 有错（可用于 CI）。

**关于第 3 项为什么是 $\ell-mg/\kappa$（曾经的判断错误，务必注意）**

本设置的势能是 $E(y)=\frac{\kappa}{2}(y-\ell)^2+m\,g\,y$（重力向量 $g_{\rm vec}=(0,-g,0)$，坐标 $y$ 向上），
$\mathrm{d}E/\mathrm{d}y=\kappa(y-\ell)+mg=0$ 给出 $y=\ell-mg/\kappa$。
**这是一个受压构型**：自由点在 $+\hat y$ 侧、pin 在原点，重力把质量往下拉，
弹性力 $\kappa(y-\ell)$ 在 $y<\ell$ 时朝 $+y$，两者反向抵消。

曾经把它误判为"代码 bug"（以为平衡应在 $\ell+mg/\kappa$），并据此改了实现，反而把正确的散射符号改坏。
要点：$y=\ell+mg/\kappa$ 是**能量的极大值**，不是平衡点；
"悬挂在固定点下方"对应坐标 $y=-\ell-mg/\kappa$、长度 $\ell+mg/\kappa$ —— 长度与坐标不能混用。

**真正需要修的缺陷（本轮已修）**：pin 消元的右端补偿项 $\kappa q$ 缺失。
pinned 行被覆盖为（对角 1，右端 $q$）等价于消去该自由度，必须把原耦合项 $-\kappa x_{\rm pin}$
的贡献移到自由端右端，否则自由点会被拉向原点。判据见 `pd_trans`（平移一致性）。

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

> **`pinPositions` 契约（踩过坑）**：一旦 `mesh.pinned` 里有 1，`pinPositions` 就必须按顶点数填好，
> 因为散射阶段要读 `pinPositions[pin]` 做消元补偿。正常流程由 `makeScene` / `refreshPinPositions`
> 保证；**手工搭 Mesh 的代码（测试、`_verify/` 程序）极易只写 `pinned` 而漏掉它**。
> 漏掉后越界读会表现为随机崩溃（0xC0000005，`pd_diraudit` 就这样挂过一次）；
> 现在 `scatterInto` 会在入口直接查这个前置条件并报错退出，而不是让它变成难定位的崩溃。

### 2.3 约束内核（方向约定与 pin 消元补偿都在这里）

| 文件 | 内容 | 相关验收 |
|---|---|---|
| `core/energy/DistanceTerm.h` | 接口：`project` / `projectOne` / `scatterInto` / `assembleMatrix` | 1、2、3 |
| `core/energy/DistanceTerm.cpp` | **三处与方向有关的约定都在这个文件里**：<br>① `project()`：`raw = x_a - x_b`，`d_c = ℓ·raw/‖raw‖`——方向只由端点顺序决定，**没有任何 pin 分支**（早期那个"强制指向自由端"的翻转已删除，见 `README.md` §4.1）<br>② `scatterInto()`（串行 / OpenMP / 非 OpenMP 三个分支）：`contribution = targets[c] * e.stiffness`，`b_a += contribution`、`b_b -= contribution`；某一端被 pin 时改为给自由端 `b_free += κ·q`（pin 消元补偿，见 `README.md` §4.2）<br>③ `assembleMatrix()`：对角 `+κ`、耦合 `-κ`，只跳过**两端都 pin** 的约束<br>三个分支（串行 / OpenMP / 非 OpenMP）必须保持一致——改一处就要改三处。 | 1、2、3 |

> **排查建议**：这三处的符号必须成对，看代码时请放在一起看，单看任何一处都"像是对的"。
> 另外散射里还有一段 **pin 消元补偿**（b_free += κ·q），它是独立的一项，容易漏 —— 见 §4 的说明。

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

> **说明**：目前形式是 $v_{n+1}=(1-k_d)(x_{n+1}-x_n)/h$。
> 重力只通过第 1 阶段的 $\hat x=x+hv+h^2g$ 进入，右端**不需要**再补外力项
> （曾经补过 $-Mg$，那是把重力算了两遍，已撤掉）。

### 2.6 应用

| 文件 | 内容 |
|---|---|
| `app/bench/main.cpp` | 无窗口基准：分阶段耗时、迭代数、应变、**分解复用统计**（符号/数值分解与回代次数） |
| `app/viewer/main.cpp` | GLFW + OpenGL 固定管线查看器；HUD 用线段自绘数字（不依赖 ImGui/字体）；键盘可调刚度/迭代/子步/线程 |

### 2.7 验证程序（`_verify/`）

| 文件 | 用途 |
|---|---|
| `_verify/check.cpp` | **端到端验收程序**：8 项检查，只输出 OK/NG |
| `_verify/translation_invariance.cpp` | 平移不变性：无重力时零位移、pin 平移后相对形状不变 |
| `_verify/chain_test.cpp` | 自由链/悬挂链与解析式对照 |
| `_verify/direction_audit.cpp` | 投影方向三项判据（全局一致 / 物理合理 / 与标准 PD 逐位一致） |

> **以下程序是排查期留下的，前提假设部分已过期，不是受支持的验收集**（见 `README.md` §4）：
> `one_step_trace`、`standard_pd`、`steady`、`rhs_breakdown`、`force_audit`、`kappa_effect`、
> `hang_test`、`init_audit`、`iteration_scan`、`scale_sweep`、`stiffness_calibration` 等。
> 它们记录了当时的推理路径，可作追溯用，但结论请以 `pd_check` 与两个测试套件为准。

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
#   4 pinned 顶点严格不动         8 弹性力符号 == -dU/dy

# ---- 测试套件（打印每条断言）----
.\build\Release\test_primitives.exe         # 68 断言（12 个测试）
.\build\Release\test_spring_vertical.exe    # 139 断言（8 个测试）

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
`--shot FILE.png` 截图（自写 PNG 编码，无图像库依赖）、`--pin-single` 只钉顶边中点、`--pin-corners` 只钉两角（默认钉整条上边）、`--iters N` `--tol T` `--damping K`。
窗口内交互：左键拖拽旋转 / 滚轮缩放 / `SPACE` 暂停 / `S` 单步 / `G` 重力 / `R` 重置 /
`[` `]` 刚度 / `-` `=` 迭代数 / `,` `.` 子步 / `;` `'` 线程数 / `ESC` 退出。

**实测**（i7-12700F，18 线程）：60×60 布料（3600 顶点 / 7080 约束），`--frames 600` 平均
**59.4 FPS**（vsync 封顶 60）、稳态物理 **1.52 ms/帧**、迭代 1 次、全程数值分解 1 次。
开局需要 30–40 次迭代收敛（物理 10–40 ms/帧），约 1 s 后进入稳态。

调试环境变量：

```powershell
$env:PD_TRACE_STEP = "1"      # stepOnce 逐步打印：预测 / 迭代残差 / 求解结果
$env:PD_CHECK_VERBOSE = "1"   # pd_check 第 8 项打印逐点明细
```

---

## 4. 当前状态

**验收 8/8 全过**，另有三组独立验证：

| 验证 | 内容 | 结果 |
|---|---|---|
| `pd_check` | 8 项验收 | 全过 |
| `pd_trans` | 平移一致性：pin 在任意位置下静止位移严格为 0；平移系统后相对形状不变 | 全过 |
| `pd_chain` | 长链条对照解析解：自由链长度精确不变；悬挂链伸长与 $mg\,N(N-1)/(2\kappa)$ 吻合 | 全过 |
| `pd_diraudit` | 方向三项：全局一致 / 物理合理 / 与标准 PD 逐位一致 | 全过 |
| `test_primitives` + `test_spring_vertical` | 68 + 139 断言 | 全绿 |

**实测性能**（i7-12700F，18 线程）：60×60 布料（3600 顶点 / 7080 约束）
查看器 `--frames 600` 平均 59.4 FPS（vsync 封顶）、稳态物理 1.52 ms/帧，数值分解全程只发生 1 次。

**已知不足**（不是缺陷，是尚未实现）：

1. 材质模型只有 4-邻域距离约束，**没有对角连接、没有弯曲项** —— 四边形可自由剪切成菱形，
   抗剪切/抗弯不足，折叠与自穿插都无法阻止；
2. **没有碰撞处理**（IPC 暂缓，见 `docs/design-discussion.md`）；
3. **刚度标定未最终确定**：$\kappa$ 与顶点质量的量级需按目标网格配好，
   否则布料在自重下会明显伸长。正确做法是把刚度按质量参数化（$\kappa=k_{\rm mat}\cdot m$），
   而不是按长度（曾经按 $\mathrm{spacing}^2$ 缩放过，是错的，已撤掉）。

**文档注意**：`README.md` 与 `docs/plan.md` 里关于"重力只出现在预测里"的结论是对的，
早期文档曾写过"右端要补 $-Mg$"，那**是错的**（等于把重力算两遍），已全部更正。
若在别处（含 git 历史或 `_verify/` 旧程序的注释）再看到这个说法，以本节为准。
