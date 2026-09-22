# 代码地图与验证方法

> 本文档只讲两件事：**每部分代码在哪**、**怎么运行并看结果**。
> 计划与设计见 `docs/plan.md`，进度与现象记录见 `README.md`，
> **收敛问题的成因与刚度标定见 `docs/pd-convergence.md`**，
> **并行改造方案见 `docs/parallel-refactor.md`（Phase 1 已实施）、性能记录见 `docs/perf.md`**。

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
| `core/mesh/Mesh.h/.cpp` | 顶点/速度/质量/pin 标记、边表（距离约束）、稀疏结构 `buildSparsityPattern()`、**约束着色**（同一入口一并构建）、规则网格生成、OBJ 载入 | 6、7 |
| `core/mesh/ConstraintColoring.h/.cpp` | 约束（边）的图着色：同色边两两不共享顶点 —— 散射据此按颜色分组执行，**无需归约/原子加**。并行 Jones–Plassmann（优先级哈希 + 最小可用颜色），颜色数 ≤ Δ+1，同色内按顶点序排列（避免假共享）。**拓扑级数据**，随 `buildSparsityPattern()` 构建 | 散射的正确性与性能 |
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
| `core/energy/DistanceTerm.h` | 接口：`project` / `projectOne` / `scatterInto` / `assembleMatrix`，以及**区域内版本** `projectInRegion` / `scatterIntoInRegion`（"一个子步一个区域"方案用，契约是"必须在已有并行区域内调用"，见 `docs/parallel-refactor.md` §4.1） | 1、2、3 |
| `core/energy/DistanceTerm.cpp` | **三处与方向有关的约定都在这个文件里**（`projectEdge` / `scatterEdge` / `assembleMatrix` 三个内核函数，串行与并行两条路径共用同一份语义，不再各写一遍）：<br>① `projectEdge()`：`raw = x_a - x_b`，`d_c = ℓ·raw/‖raw‖`——方向只由端点顺序决定，**没有任何 pin 分支**（早期那个"强制指向自由端"的翻转已删除，见 `README.md` §4.1）<br>② `scatterEdge()`：`contribution = targets[c] * e.stiffness`，`out_a += contribution`、`out_b -= contribution`；某一端被 pin 时改为给自由端 `out_free += κ·q`（pin 消元补偿，见 `README.md` §4.2）<br>③ `assembleMatrix()`：对角 `+κ`、耦合 `-κ`，只跳过**两端都 pin** 的约束<br>求和规则只有一条：**按约束着色逐色累加**（同色边互不相邻 ⇒ 同一个顶点在每一色里最多被写一次 ⇒ 只读改写、无竞争、无归约），顺序与线程数无关，见 `core/mesh/ConstraintColoring.h` | 1、2、3 |

> **排查建议**：这三处的符号必须成对，看代码时请放在一起看，单看任何一处都"像是对的"。
> 另外散射里还有一段 **pin 消元补偿**（b_free += κ·q），它是独立的一项，容易漏 —— 见 §4 的说明。

> **并行分支在哪被验证**：`scatterInto` 在"边数 < 256 或单线程"时走串行回退，所以此前所有
> 散射测试（都是单条边网格）**从未执行过 OpenMP 分支**。2026-09-21 补了两个测试
> （`tests/primitives/test_distance_term.cpp` 里的 `parallelScatterMatchesSerialReference`
> 与 `scatterIsThreadCountIndependent`），用 20×20/760 边的扰动网格把并行分支纳入验收；
> 实测并行与串行**位级不同**（最大差 7.11e-15），详见 `README.md` §3.2 与 `docs/parallel-refactor.md` §2.4。
>
> **哪条路径会并行**：`stepOnce` 用一个 `if` 子句决定是否真的 fork
> （`numThreads() > 1 && 顶点数 ≥ 512 && 边数 ≥ 512`，阈值是经验值，Phase 2 后要重标）；
> 小网格走"串行化的区域"—— 全项目只有一份流程实现，没有运行时双实现开关。

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
| `core/sim/Integrator.cpp` | `stepOnce()` 的流程（**Phase 1 之后：一个子步只进入一个并行区域**，见 `docs/parallel-refactor.md` §4.1），按阶段编号：<br>**1) 判定+数值分解**（区域外，串行）——只在 stamp 变化时执行<br>**2) 预测** $\hat x=x+hv+h^2g$（区域内的 `omp for`；**融合**了"保存上一步位置"与"初始化 previous"）——**重力只在这里出现一次**<br>**3) 惯性右端** $b_{\text{base}}=(M/h^2)\hat x$（`omp single`）<br>**4) PD 迭代**：局部步(`projectInRegion`) → 散射(`scatterIntoInRegion`，融合了基值) → 覆盖 pin+全局回代（同一个 `single`） → 解包+收敛扫描+previous+**非线性残差**（融合成一趟，`omp for nowait` + 每线程部分最大值；**Phase 3 之后残差是逐顶点 gather，不再是串行趟**，见 `docs/perf.md` §7） → 判据（`omp single`，只做合并与判定）<br>**5) 速度更新与阻尼**（`omp for nowait`） | 2、3、5、6 |

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
| `_verify/solve_audit.cpp` | **全局步求解器可行性审计**（`pd_solveaudit`）：延迟/带宽性质、因子层集与关键路径、换排序/分解的对照，以及"层调度并行回代"的原型（含正确性自证）。结论见 `docs/solver-feasibility.md`；**只测量，不改生产代码** |

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
.\build\Release\test_primitives.exe         # 125 断言（16 个测试；含并行散射分支与约束着色的对照）
.\build\Release\test_spring_vertical.exe    # 155 断言（11 个测试）
.\build\Release\test_convergence_criterion.exe  # 15 断言（4 个测试，收敛判据与外层迭代质量）

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
# 收敛门槛：--residual-tol T 设的是**真正放行的那一个**（T×|g|，默认 1e-3）；
# --tol 只改位移判据，而残差判据启用时位移判据不参与判定（调它不改变迭代数）。
# 例：门槛与代价的因果（同工况 40×40 / 300 子步 / 上限 40，2026-09-21 实测）——
.\build\Release\pd_bench.exe --grid 40 40 --steps 300 --iters 40 --residual-tol 100
#   → 平均迭代 3.75、单子步 1.4 ms：门槛极松，几乎立刻"达到判据"
.\build\Release\pd_bench.exe --grid 40 40 --steps 300 --iters 40 --residual-tol 0.3
#   → 平均迭代 40.00、单子步 12.2 ms、**收敛 0/300**、末次残差 8.64 m/s²：
#     该工况连 0.3·g 都达不到（与"bench 的 κ=1e4 在 40×40 下不收敛"的注记一致）
.\build\Release\pd_bench.exe --grid 40 40 --steps 300 --iters 40 --residual-tol 0
#   → 换判据而非收紧：_residual-tol 0 是**禁用**残差判据、回退到位移判据 --tol
#     （平均迭代 39.32、收敛 41/300）
# 注：bench 与 viewer 的工况不同（间距 0.025 vs 0.02、pin 两角 vs 整条上边），
#   所以 bench 不能直接复现 viewer 的帧率数字，只能用来观察门槛的因果。
# 输出里有"收敛：N/M 子步达到判据（放行判据 …）；末次残差 …"：
# 平均迭代数相同，**"达到判据"与"跑满上限"性质完全不同**，报告时必须区分。

# ---- 实时渲染（GLFW + OpenGL，可截图）----
.\build\Release\pd_viewer.exe                                   # 交互窗口
.\build\Release\pd_viewer.exe --grid 60 --frames 300 --shot out.png   # 有界运行 + 截图
```

查看器参数：`--grid N` 网格边长、`--stiffness K` 刚度（默认 2300）、`--frames N` 跑够帧数自动退出、
`--shot FILE.png` 截图（自写 PNG 编码，无图像库依赖）、`--pin-single` 只钉顶边中点、`--pin-corners` 只钉两角（默认钉整条上边）、`--iters N` `--damping K` `--dt H`、
**`--residual-tol T`（真正的放行门槛，默认 0.3；`<=0` 禁用并回退到位移判据 `--tol T`）**、
`--no-early-exit`（关掉全部提前退出，优先级最高）。
窗口内交互：左键拖拽旋转 / 滚轮缩放 / `SPACE` 暂停 / `S` 单步 / `G` 重力 / `R` 重置 /
`[` `]` 刚度 / `-` `=` 迭代数 / `,` `.` 子步 / `;` `'` 线程数 / `ESC` 退出。

**实测**（i7-12700F，18 线程）：60×60 布料（3600 顶点 / 7080 约束），`--frames 600` 平均
**59.4 FPS**（vsync 封顶 60）、稳态物理 **1.52 ms/帧**、迭代 1 次、全程数值分解 1 次。

**门槛必须一起写**：上面这组数字是查看器默认档（`--residual-tol 0.3`，即门槛 0.3·g）的稳态。
查看器**覆盖**了库默认值 —— `SceneConfig::residualTolerance` 默认是 **1e-3**，
在 40 次迭代预算内压不到，会跑满迭代（物理约 47 ms/帧、约 20 FPS）。见 `core/sim/Scene.h`。

**瞬态（2026-09-21 复跑更正）**：原先这里写"开局需要 30–40 次迭代收敛，约 1 s 后进入稳态"，
实测并不成立：

| 阶段 | 实测（`--grid 60 --frames 600`，默认 0.3 档） |
|---|---|
| 第 1 秒 | 迭代 26、判定"收敛"、残差 2.941 —— 这一小段与旧描述相符 |
| 约 1–10 s | 残差升到 **15.31** 后缓慢回落（8.98 → 9.07 → 8.17 → 5.08 → 4.18 → 3.31），**始终高于 0.3·g = 2.943**；40 次迭代**全部用尽**（HUD 标"用尽迭代"），物理 45–52 ms/帧、约 20 FPS |
| 约 600 帧（≈10 s 仿真时间） | 残差跨过门槛 → 迭代 1 次、物理 1.85 ms/帧、60 FPS、应变 0.0100 |

两次独立运行（600 / 700 帧）的残差序列**逐位相同**，说明这是确定的物理行为而非抖动，
切换点由子步数唯一决定。**稳态数字本身复现无误**（迭代 1、1.85 ms/帧、60 FPS、应变 1%、分解 1 次）。

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
| `test_primitives` + `test_spring_vertical` + `test_convergence_criterion` | 125 + 155 + 15 断言 | 全绿 |

**实测性能**（i7-12700F，18 线程）：60×60 布料（3600 顶点 / 7080 约束）
查看器 `--frames 400` 实测 60 FPS（vsync 封顶）、稳态物理 1.9 ms/帧、迭代 1 次、
应变 1%、数值分解全程只发生 1 次。（数字为 2026-09-20 刚度标定后的复测；
标定前记录为"59.4 FPS / 1.52 ms/帧"，那组数字来自会过早退出迭代的旧判据，
可比性有限 —— 见 `docs/pd-convergence.md`。）

> **这组数字的门槛是 `--residual-tol 0.3`，不是库默认 1e-3**，且指的是**稳态**：
> 2026-09-21 复跑时 `--frames 600` 的前约 10 s 仍在瞬态（40 次迭代用尽、物理 45–52 ms/帧、
> 约 20 FPS），之后才落到上表状态。所以引用时必须写清"门槛 + 是稳态还是启动段"（详见 §3）。
> 同一轮复跑还把**并行净减速**测出来了（散射阶段 18 线程比 1 线程慢 10×），
> 属 M2 未处理事项，现状与数据见 `README.md` §5 第 7 项。

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
