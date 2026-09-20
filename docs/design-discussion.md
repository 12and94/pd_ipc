# PD + IPC 布料仿真器：三条实现路线的可行性分析

> 讨论稿（未开始实现）。目标：判断三个方向是否合理、能否共用一套工程、以及 IPC 与 PD 结合是否可行。
> GPU 部分按 Vulkan compute（compute shader + VkBuffer/SSBO）考虑，不涉及 CUDA/DX12。
>
> **状态（已确认决策，见 `docs/plan.md`）**
> - 精度目标：**实时交互**（60Hz 帧预算，float32 可用，接受近似）。
> - 网格规模：**任意大小**（需自适应：小网格走 CPU、大网格走 GPU）。
> - 碰撞范围（初期）：**仅布料–刚体**，跳过自碰撞。
> - 摩擦：**先无摩擦**。
> - 方向 2 的全局步：**Chebyshev 半迭代**（不是裸 Jacobi）。
> - IPC 相关实现**暂缓**；本文 §5、§6 保留为后续参考，不进入当前计划。
> - 起始目标：**纯 CPU 方向 1**。
> - **铁律**：本项目禁止参考 `dp_ipc` 工程（本项目前身）——不读其源码、不搬运、不引用其结论、不与其对照、不依赖其路径。**其他同类实现不设限，可读可引用**。正确性只来自独立推导 + 原语不变量 + 解析解对照。见 `docs/plan.md` §10.1。

---

## 0. 结论先行

| 问题 | 结论 |
|---|---|
| 三个方向都合理吗？ | 是。但它们是**同一算法家族的三档后端**，不是三个独立项目；共享率估计 70–85%。 |
| 值得都做吗？ | 值得，但不要并列开工。建议顺序：①CPU 多线程（参考实现 + 正确性基线）→ ③GPU 局部/CPU 全局（性能性价比最高、最容易加碰撞）→ ②全 GPU（目标形态）。 |
| IPC + PD 合适吗？ | 合适，且是当前热门方向。关键认识：**IPC 是接触模型，PD 是求解器**，二者正交，不是竞争关系。 |
| 有现成工作吗？ | 有：PD-with-dry-frictional-contact (2020)、Penetration-free PD on GPU (2022)、Subspace-Preconditioned GPU PD with Contact (SA'23)、GIPC/StiffGIPC (2024–2025) 等。见 §6。 |
| 最大风险 | PD 的局部步是 **prox 算子**而非 IPC 的原始能量；把非线性 barrier 塞进局部步会破坏 PD 的固定矩阵优势与收敛性。需要"冻结刚度 / filter / 外层线搜索"三件套。见 §5。 |
| 工程上最省力的路径 | 不要把 PD 和 IPC 写成两个系统。做成 `IncrementalPotential + SolverBackend`，让 Newton+PCG 和 PD 成为同一个能量泛函下的两种后端。见 §4.3。 |

---

## 1. 统一数学框架（三条路线的共同底座）

### 1.1 PD 的局部–全局迭代

隐式欧拉下最大化消散，记预测位置 $\hat x = x^t + h v^t + h^2 g$。PD 引入辅助变量 $d_c$（约束 $C_c$ 的投影目标），用**固定权重**的代理能量：

$$\tilde g(x,d) = \frac{1}{2}x^\top M x + \underbrace{\sum_c \frac{w_c}{2}\|A_c x - d_c\|_F^2}_{\text{代理弹性}} \;-\; x^\top M \hat x$$

- **局部步**：$d_c \leftarrow \arg\min_{d}\ \frac{w_c}{2}\|A_c x - d\|^2 + \Psi_c(d)$，逐约束/逐单元完全并行、无耦合。对弹簧/三角形拉伸、ARAP 四面体、弯曲等都存在闭式解。
- **全局步**：解**常数**系统
  $$L x = M\hat x + \sum_c w_c A_c^\top d_c,\qquad L = M + \sum_c w_c A_c^\top A_c$$
  $L$ 是加权的 cotan/sqrt-Laplacian 组合 + 质量对角阵，**与当前位形无关**。

这一点是三条路线共享的关键性质：$L$ 的结构和数值只依赖拓扑、权重 $w_c$、质量，**不依赖 $x$**。于是：

- 路线 1 可以**一次符号+数值分解，跨帧复用**；
- 路线 2 的 SpMV/迭代求解器只需缓存 CSR 与对角预条件子；
- 路线 3 的两端只需传 $b$ 和 $x$，不传矩阵。

### 1.2 接触的两种挂法（决定了 IPC 怎么进来）

**挂法 A：接触进入全局矩阵**（IPC/Newton 的原生方式）

$$E(x) = \tfrac12\|x-\hat x\|_M^2 + h^2\Psi(x) + h^2\kappa B(x) + h^2 D(x), \qquad H = M + h^2\nabla^2\Psi + h^2\kappa\nabla^2 B + \dots$$

问题：$B$ 的 Hessian 随 $x$ 变化（近接触时量级 ~$\kappa/d^2$，$\kappa$ 可达 $10^9$），矩阵每步都变 → CPU 需要反复符号/数值分解，GPU 上 Jacobi 条件数崩坏。这是"IPC 直接 + 纯迭代 GPU 求解"困难的根本原因。

**挂法 B：接触作为局部步的 prox / filter**（PD 的原生方式）

在局部步里对每个顶点求
$$y_v \leftarrow \operatorname*{arg\,min}_{y\in \mathcal F_v(x)}\ \tfrac{\alpha_v}{2}\|y - z_v\|^2 \quad(\mathcal F_v:\ \text{由接触构造的可行集})$$
全局矩阵里**不出现** barrier Hessian，$L$ 依旧常数。这正是所有"PD + 接触"论文的做法（PD dry friction 2020、penetration-free PD 2022、subspace-preconditioned GPU PD 2023）。

> **结论：路线 2 与 IPC 结合得最自然。** PD 的"矩阵固定"恰好是 GPU Jacobi 能跑得动的前提，而把接触做成 filter/projection 恰好避免了 barrier Hessian 污染矩阵。这不是巧合，而是二者互补。

### 1.3 IPC 侧需要什么

IPC 的增量势能（其官方 toolkit 文档给出的形式）：

$$x^{t+1} = \arg\min_x\ \tfrac12 (x-\hat x)^\top M (x-\hat x) + h^2\Psi(x) + h^2\kappa B(x) + h^2 D(v(x))$$

其中 barrier $B$ 作用于**非相邻**图元对的**无符号距离**，$D$ 是平滑的摩擦耗散项，配合 CCD 线搜索保证无穿透。需要的积木：碰撞候选（broad/narrow phase）、距离与最近点（PT/EE/PE/PP）、CCD、barrier 的梯度/Hessian、摩擦的切向算子。这些都在 [ipc-toolkit](https://ipctk.xyz/tutorials/simulation.html) 里有成熟实现，可作为**参考实现与正确性 oracle**，但它是 CPU/Eigen 且以 robust 优先，不适合直接当性能内核。

---

## 2. 三条路线逐条评估

### 路线 1：CPU 多线程——预测/局部投影并行，全局用稀疏 Cholesky 预分解后串行求解

**合理性：高。** 这是学术 PD 实现的标准形态，也是全项目的地基。

要点与坑：

1. **"预分解"要分清三种含义**（这是最常见的误解来源）：
   - *符号分解*（elimination tree / fill-reducing ordering）：只要拓扑不变就永久复用。
   - *数值分解*：只在 $w_c$、质量、约束集变化时才需要重做。PD 里 $w_c$ 固定 → 通常**只在换材料/换约束集时做一次**。
   - *每帧求解*：只是三角回代（forward/backward substitution）。
   建议实现成 `Solver::analyze() / factorize() / solve()` 三段，并把"何时 invalidate"显式建模（材料、$\kappa$ 是否进入矩阵、约束集是否变化）。

2. **回代不是"串行"**：稀疏 Cholesky 的求解器本质是稀疏三角求解，有依赖但要乘**多右端项**，可以按超节点/消去树并行。建议：
   - 用 CHOLMOD/Eigen SimplicialLDLT（原型）→ 后期换 SuiteSparseQR 级别的？不，用 CHOLMOD + 多线程 BLAS（OpenBLAS/MKL）；
   - 把**全局步的多次求解批处理**：PD 每个时间步通常只需 1 次全局求解（若做 inner iteration 预热则多右端一起解）。

3. **真正的瓶颈排序**（30k–100k 顶点级布料，单帧典型量级）：
   - 全局回代：几 ms（远小于局部步，除非网格极大）
   - **局部步投影**：与约束数成正比，容易并行，OpenMP/TBB 的 `parallel_for` 即可
   - **碰撞检测（一旦加 IPC）**：broad phase + narrow phase + **CCD** — 这才是 CPU 端的头号开销，必须并行 BVH 遍历与批量 CCD
   - 组装/归约（stencil 散射）：注意原子竞争，建议按颜色/分块局部缓冲再归约

4. **该路线独有的价值**：它是唯一能给出"参考解"的实现——能量单调性、投影残差、与 Newton+PCG 的收敛对比都要靠它。**不要跳过**。

5. 可选增强：并行局部步 + 串行全局步本身可用流水线重叠（局部步算 $d$ 的同时做上一次的全局求解），但收益有限；更大的收益来自**多重网格**或**子空间加速**。

### 路线 2：GPU 全程——compute shader 做预测/局部投影，全局步在 GPU 上 Jacobi 迭代

**合理性：高（是目标形态），但"Jacobi"是其中最弱的一环。**

要点与坑：

1. **Jacobi 的收敛性**：对 $L x=b$（$L$ 是加权 Laplacian + 质量），Jacobi 的谱半径 $\approx 1 - O(h^2)$，细网格/强各向异性下迭代数会到 **几十到几百**，且每步都是全内存带宽的 SpMV。学术共识与经验：
   - **Chebyshev 半迭代**：几乎不增加成本（只需要 $\lambda_{\max},\lambda_{\min}$ 估计），收敛显著优于 Jacobi，是 PD GPU 实现的事实标准（PD 原文即推荐 Chebyshev + Jacobi 预条件）。
   - **分块（$3\times3$）Jacobi / 逐顶点逆对角**：比标量 Jacobi 好。
   - **图着色 Gauss–Seidel / 红黑**：迭代数更少，但并行度下降、需要多趟 barrier 同步；在 Vulkan 里三趟 pass + 屏障，仍可实现。
   - **多重网格预条件**：迭代数可降到常数级，但实现复杂度高（需要网格层次构造、延拓/限制算子），且 cloth 的约束图与几何网格不完全一致。列为后期项。

2. **Vulkan compute 的硬约束（这条路线真正的工程难点）**：
   - **浮点原子加**：core Vulkan 的 `atomicAdd` 只保证整数；shader float atomics 需要 `VK_EXT_shader_atomic_float`（桌面 NVIDIA/AMD 基本支持，Intel/MoltenVK 不稳）。**不要依赖 float atomic** 做组装/归约。替代方案：
     (a) **散射到每单元/每约束的私有缓冲**，再做一次并行归约树；
     (b) 序列化为 **定点整数** 后用 `atomicAdd(int64)`（排序/求和单调，可复现）；
     (c) warp/subgroup 内 `subgroupAdd` 归约 + 每 warp 一个槽位，最后第二趟归约。
   - **精度**：`shaderFloat64` 在消费级 GPU 上要么不支持要么极慢（1/16~1/32）。**必须用 float32 设计**，并认真做"精度策略"：以**相对残差**为收敛判据、对 $M\hat x$ 与 $b$ 做补偿求和（Kahan）、避免 $\kappa$ 巨大数与位置小量直接相减（改用**间隙量** $d-\hat d$ 作为主变量）。
   - **描述符/内存管理**：为每个语义单独的 storage buffer（位置、速度、单元索引、约束目标 $d$、$b$、残差、碰撞对列表、BVH 节点），用 descriptor array 或"大 buffer + 偏移"以规避描述符数量上限；管线用 specialization constants 编译"拉伸/弯曲/接触"的变体而不是运行时分支。
   - **同步**：一次 substep 内多个 pass（预测 → 组装 → Chebyshev 迭代 → 局部投影 → 收敛检查）之间用 `vkCmdPipelineBarrier`，全部记录进同一 command buffer；用 2–3 个 in-flight frame 重叠 CPU 侧逻辑。**收敛检查需要回读**：把残差写成 counter/flag，用一小块 host-visible buffer 只拷 4–8 字节，避免全量回读。

3. **每步迭代次数与收敛判定**：PD 论文/TOD 通常固定 inner iteration 数（如 5–20）而非解到收敛；GPU 版本同样如此。但**加 IPC 后必须有真正的收敛判据**（见 §5.5），否则无穿透保证会退化。

4. **和路线 3 的性能分水岭**：路线 2 只有在"顶点数 ≳ 50k–100k、且帧预算 < 5ms"时才明显胜出；小网格下 GPU 本可闲置，CPU 直接求解反而更稳更快。

### 路线 3：GPU 做预测/局部投影，传回 CPU 求解全局

**合理性：非常高，是被低估的实用路线；在"要加 IPC"的前提下它甚至可能是最优起点。**

1. **同步/Copy 代价可接受**：传递量只有 $b$（$3N$ 浮点，$N=100\text{k}$ 时 1.2 MB）与 $x$。PCIe 4.0 x16 双向各 ~25 GB/s → 单次传输几十 μs。真正的风险是**强制同步点**（fence 等待 + staging buffer 往返）。处理办法：
   - **双缓冲 / 跨帧流水**：CPU 求解第 $k$ 帧的同时，GPU 做第 $k$ 帧的碰撞检测/候选筛选与第 $k+1$ 帧的投影预计算；
   - 把 CPU 全局求解放进独立线程，用事件同步而不是忙等。

2. **为什么它在加 IPC 后仍可能胜出**：
   - IPC 的 barrier 组装、CCD、摩擦切向算子、线搜索的**控制流复杂**，CPU 上调试成本低一个数量级；
   - 一旦 barrier Hessian 进入矩阵（挂法 A），矩阵每步都变，CPU 的 CHOLMOD 虽然要重分解，但仍比 GPU 迭代稳健得多；
   - 接触区域通常只占网格很小一部分，**局部稀疏更新**（只更新受影响的块 + 增量分解 / 低秩更新）在 CPU 上更容易做。

3. **建议的定位**：路线 3 作为**"IPC 首个可跑版本"**，随后把全局步逐步迁移到 GPU（即演进为路线 2）。这也意味着路线 3 不是浪费——它和路线 2 共享除全局求解器以外的全部代码。

### 2.4 三路线对比表

| 维度 | ① CPU 多线程 + 直接法 | ② 全 GPU + Jacobi/Chebyshev | ③ GPU 局部 + CPU 全局 |
|---|---|---|---|
| 主要目的 | 正确性基线、可调参考 | 最终实时形态 | 过渡形态 / 碰撞调试友好 |
| 全局求解 | 预分解 Cholesky，回代 | Chebyshev/Jacobi，多次 SpMV | CPU Cholesky（可重分解） |
| 迭代数敏感度 | 不敏感（精确解） | 敏感（几十~几百） | 不敏感 |
| IPC barrier 进矩阵 | 尚可（重分解代价） | 差（条件数崩坏） | 好 |
| 数据传输 | 无 | 无 | 每子步 ~MB 级，可流水 |
| 主要风险 | 局部步/CCD 并行效率 | float32 精度、原子、收敛慢 | 同步 stall |
| 实现难度 | 低 | 高 | 中 |
| 建议顺序 | 1st | 3rd | 2nd |

---

## 3. 三者能共用多少？——共享清单

**可以共用（估计占代码量 70–85%）**

| 模块 | 说明 | 三条路线 |
|---|---|---|
| 网格与拓扑（顶点/边/面/四面体、邻接、拉普拉斯权重） | 一次做好 | 1/2/3 |
| 约束集与投影算子（弹簧、三角拉伸、弯曲、ARAP） | 局部步的核 | 1/2/3 |
| `project_*` 内核 | CPU 上 OpenMP、GPU 上 compute shader，**数学与接口完全一致** | 1/2/3 |
| 组装：$b = M\hat x + \sum w_c A_c^\top d_c$ 与 $L$ 的 CSR 结构 | $L$ 结构只依赖拓扑 | 1/2/3 |
| 稀疏求解器抽象（`ISparseSolver`） | 具体后端不同 | 1/2/3（实现不同） |
| 碰撞候选/距离/CCD/barrier/摩擦 | IPC 部分 | 1/2/3 |
| 时间积分、子步、收敛判据、统计埋点 | 完全共用 | 1/2/3 |
| 场景/OBJ 加载、渲染、回归测试 | 完全共用 | 1/2/3 |

**必须分叉的部分（15–30%）**

- 全局求解器后端：CHOLMOD 回代 / GPU Chebyshev / CPU Cholesky-on-demand。
- 内存域与数据布局：CPU SoA vs GPU SSBO + staging。
- 数值精度策略：CPU 可 float64；GPU 必须 float32 并额外做残差/补偿处理。
- IPC 接触的求值位置：CPU 上直接算 vs GPU 上算（挂法 A 与 B 的差异）。

**结论：应当做成一个项目、一套核心库、可切换后端。** 三个 `main` 或者三个场景配置，共享 90% 代码。

---

## 4. 建议的工程架构

### 4.1 目录草案

```
pd_ipc/
  core/          mesh, topology, constraints, projection kernels(数学), energy model
  assembler/     CSR 结构, b 与 L 的组装(MPI 无关的 CPU 版 + GPU 版)
  solver/        ISparseSolver, IGlobalSolver
    cpu_chol/    analyze/factorize/solve(symbolic 复用)
    gpu_cheb/    Vulkan compute: SpMV, Chebyshev, 残差归约
  contact/       broadphase, narrow, distance, ccd, barrier, friction   (IPC)
  ipc/           能量/梯度/Hessian 组装, PSD 投影, 线搜索, κ 自适应
  vk/            Vulkan 基础设施: device/queue/allocator/pipeline/descriptor/pass graph
  app/           场景、渲染、UI、帧循环
  tests/         (a) 单元测试  (b) 物理回归(能量单调/穿透/收敛)  (c) 与 ipc-toolkit 对比
```

### 4.2 关键接口（示意）

```cpp
struct IEnergyField {                 // 弹性 / 接触 都实现它
  virtual double energy(const Vec& x) const = 0;                 // 全局步线搜索用
  virtual Vec    gradient(const Vec& x) const = 0;               // Newton 用
  virtual void   hessian(const Vec& x, SpMat& H) const = 0;      // Newton 用
  virtual void   localProjection(const Vec& x, Vec& d) const = 0;// PD 局部步 (prox)
};

struct IGlobalSolver {                // 路线 1/2/3 的差异集中在这里
  virtual void analyze(const SpMat& L) = 0;
  virtual void factorize(const SpMat& L) = 0;
  virtual void solve(const Vec& b, Vec& x) = 0;
  virtual bool matrixChanged() const = 0;   // 由约束权重/κ 是否进矩阵决定
};
```

`IEnergyField` 上同时具备 `gradient/hessian` 与 `localProjection`，就能让 **Newton+PCG** 与 **PD** 共享同一份能量与接触代码——这正是"IPC 与 PD 同一项目"的落点。

### 4.3 强烈建议的一条架构决策

**不要让 PD 取代 IPC，而是让 PD 成为 IPC 能量泛函下的一种求解后端。**

即优化器层：

```
IncrementalPotential E(x) = 惯性 + h²Ψ(x) + h²κB(x) + h²D(x)
   ├─ SolverBackend::NewtonPCG      (参考 & 接触正确性基线)
   ├─ SolverBackend::ProjectiveDynamics (路线 1/2/3 的局部步 + 固定矩阵)
   └─ SolverBackend::GaussNewton(可选, GIPC 思路)
```

好处：
- 接触/CCD/摩擦/候选筛选**只写一份**；
- 可以直接做对照实验（同一场景，PD vs Newton 的迭代数/能量/穿透量），这正是评估"PD+IPC 是否可靠"的唯一途径；
- 避免"PD 版和 IPC 版各有一套碰撞检测"的经典灾难。

---

## 5. IPC × PD 结合：可行，但有 5 个必须先想清楚的技术点

### 5.1 概念澄清（避免方向性误解）

- **IPC** = 一种接触模型（barrier 势 + 摩擦势 + CCD 保证无穿透），以及围绕它的一套**优化式时间积分**。它本身不规定用哪种优化器；论文里的 projected Newton 只是其中一种（arXiv 综述也指出 barrier 兼容 preconditioned projected Newton 等多种优化方法）。
- **PD** = 一种求解器（局部–全局，固定系统矩阵）。
- 所以"IPC 和 PD 结合"= **用 PD 求解带 IPC barrier 的增量势能**。可行、已有工作、且互补（见 §1.2）。

### 5.2 难点 1：局部步是 prox，不是能量本身

PD 局部步最小化的是 $\frac{w_c}{2}\|A_cx-d\|^2 + \Psi_c(d)$；对 barrier 项就是求近端算子 $\mathrm{prox}_{\kappa B}$。这有三个后果：

- 单约束、单顶点的 barrier prox **没有闭式解**，需要小规模牛顿（1D 或 3D），成本高于弹性投影；
- 顶点可能同时触及多个接触（PT/EE 混合），**投影之间互斥/耦合**，不能"每个接触独立投影再平均"；
- 必须"取可行集内最优"，而不是"沿法向推一下"。

### 5.3 难点 2：filter / 可行集构造（这也是主流解法）

把"多接触投影"退化为**约束可行集** $\mathcal F_v$：对顶点 $v$，收集所有潜在接触，构造一个（凸的）局部可行域，然后

$$y_v = \Pi_{\mathcal F_v}(z_v),\quad \mathcal F_v \approx \{\,y : n_i^\top (y - p_i) \ge \hat d_i \ \forall i\,\}\ \text{或球近似} \cup \text{半空间近似}$$

- 每个接触给一个**半空间/球**近似（由距离与法向给出），交集是凸的 → 投影有闭式或 1D 牛顿；
- 与 CCD 结合给出"安全余量" $\hat d_i > 0$，天然得到无穿透的局部步（这是 penetration-free PD 2022 的思路）；
- 代价：近似 → 与真实 barrier 的动力学有偏差；$\hat d$ 取太大则动力学失真，取太小则局部牛顿迭代爆炸。**这个折中是这类方法的核心研究点，也是我们最该做实验的地方。**

### 5.4 难点 3：收敛性与保证的降级

- Newton+IPC 有明确的下降性/单调性论证（能量单调下降、CCD 线搜索保证可行）。
- PD 的局部–全局只对**代理能量** $\tilde g$ 单调，对真实增量势 $E$ **没有**同样保证；局部步取 prox 后，全局解可能落到能量更差的点。
- 因此 IPC 的三项保证（无穿透 / 无翻转 / 有能量界）在 PD 下**必须重新论证或降级**：
  - 无穿透：靠 $\mathcal F_v$ 的余量 + 外层 CCD 线搜索（沿全局步方向做**真实能量** $E$ 的过滤线搜索）来保；
  - 无翻转：靠参考映射（reference map）或对单元做 PSD/翻转检测投影；
  - 能量单调：只能实测（$E(x_{k+1}) \le E(x_k)$ 作为回归测试项，不成立就调 $\hat d$ / 加外层线搜索）。
- **务实建议：保留一层外层过滤线搜索**（对真实 $E$ 做，含 CCD 上界），即使这意味着 PD 不再"纯"。

### 5.5 难点 4：barrier 刚度 $\kappa$ 与刚性问题

- $\kappa$ 自适应（IPC 原论文做法）在 PD 下要重新设计：$\kappa$ 不进 $L$，只影响 $\hat d$ 与局部 prox 的"硬度"。可以把 $\kappa$ 的物理含义映射为 filter 余量 $\hat d \leftarrow \hat d(\kappa, d_{\text{hat}})$；
- $d_{\text{hat}}$ 很小、相对速度大时，prox 越来越接近硬约束 → 局部牛顿迭代数上升。需要在局部步做**热启动**（上一迭代的 $d$ 作为初值，投影解在迭代间变化很小）；
- 与摩擦 $D(v(x))$ 的耦合：PD dry frictional contact (2020) 给出了局部步解摩擦滑动的做法（基于摩擦锥投影）；IPC 的平滑摩擦势需要 $\epsilon_v$ 正则化，也可做成局部 prox。**摩擦是必须早做决定的部分**，因为它决定局部步的接口形状。

### 5.6 难点 5：性能现实

- 若坚持"接触完全用局部步、矩阵恒定"，则全局步很便宜（GPU 友好），但**局部步会变成新的瓶颈**（多次小牛顿 + 邻域查询 + 热启动数据依赖）；
- 若允许"接触进矩阵"（挂法 A），CPU 直接法要频繁重分解（或用 block 更新），GPU 迭代法条件数崩坏 → 与路线 2 冲突。
- 所以：**路线 1/3 更适合先做 IPC；路线 2 上的 IPC 需要先把 §5.3 的 filter 方案做扎实。**

---

## 6. 相关工作与可利用的现成件

| 工作 | 对本项目的用处 |
|---|---|
| [IPC (Li et al. 2020)](https://ipc-sim.github.io/) | barrier + CCD + 摩擦的原始定义；能量/梯度/Hessian 与 κ 自适应的规范 |
| [IPC Toolkit](https://ipctk.xyz/tutorials/simulation.html) | 碰撞候选、距离、CCD、barrier/friction 势、`CollisionMesh`（含 codim 顶点、位移映射）、Hessian assembler 复用、PSD 投影（CLAMP/ABS）。**建议用作正确性 oracle 与算法参考，而非性能内核** |
| [Projective Dynamics with Dry Frictional Contact (Ly et al. 2020)](https://dl.acm.org/doi/abs/10.1145/3386569.3392396) | PD 局部步做接触+摩擦投影的范式 |
| [Penetration-free Projective Dynamics on the GPU (SIGGRAPH 2022)](https://dl.acm.org/doi/abs/10.1145/3528223.3530069) | GPU 上 PD + 无穿透的局部步构造（filter 思路的直接参考） |
| [Subspace-Preconditioned GPU PD with Contact for Cloth (SA 2023)](https://dl.acm.org/doi/fullHtml/10.1145/3610548.3618157) | 与路线 2 最接近：GPU PD + 接触 + 子空间预条件（用于提升 GPU 迭代收敛） |
| [GIPC (2024)](https://dl.acm.org/doi/pdf/10.1145/3643028) / [StiffGIPC (2025)](https://dlnext.acm.org/doi/full/10.1145/3735126) / [Advancing GPU IPC (2024)](https://ar5iv.labs.arxiv.org/html/2411.06224) | IPC barrier 的 Gauss–Newton / GPU 化最新进展；对"barrier 能否 GPU 化"给出肯定答案与具体技巧 |
| [Robust and Efficient Penetration-Free Elastodynamics without Barriers (2025)](https://arxiv.org/pdf/2512.12151v1) | 无 barrier 的替代路线；提到 barrier 兼容 projected Newton/预条件等多种优化方法 |
| [ADMM ⊇ Projective Dynamics (2017)](https://doi.org/10.1109/TVCG.2017.2730875) | PD = ADMM 特例；理解"约束/接触如何并入局部步"的理论基础 |
| [DiffCloth (TOG)](https://people.csail.mit.edu/liyifei/uploads/diffcloth-highres-tog.pdf) | 含干摩擦接触的 PD 布料实现细节（可借鉴的工程写法） |

> 备注：ACM DL / arXiv 的 PDF 在本环境下多数无法直接抓取（403 或 content-type 限制），上表结论结合了可抓取的 IPC Toolkit 文档与检索摘要。若需要精确公式（如 subspace PD 的全局步构造、GIPC 的 GN 细节），建议由你提供 PDF 或在联网环境下补读。

---

## 7. 分阶段路线图（建议）

**阶段 0：骨架与参考实现（路线 1）**
- 网格/拓扑/约束/投影算子 + `ISparseSolver`（CHOLMOD 或 Eigen SimplicialLDLT）
- 无碰撞：拉伸/弯曲布料，验证能量单调、PD 与 Newton+PCG 解一致性
- 输出基线性能表（顶点数 × 迭代数 × 各阶段耗时）

**阶段 1：Vulkan 基础设施 + GPU 局部步（路线 3 的前半）**
- device/allocator/pipeline/pass 抽象；float32；禁止依赖 float atomic 的归约方案
- 把 `localProjection` 与组装搬到 compute shader，与 CPU 版逐位对照（容差内一致）
- 先只替换局部步，全局步仍走 CPU Cholesky → 得到路线 3 的完整形态

**阶段 2：IPC 接触（先 CPU，后 GPU）**
- 候选筛选 + 距离 + CCD + barrier 能量/梯度/Hessian，`ipc-toolkit` 做 oracle
- 在 **Newton+PCG 后端**上先跑通（保正确性），再切到 PD 后端
- 落地 §5.3 filter 局部步 + §5.4 外层过滤线搜索；回归测试：穿透量 = 0、能量单调、收敛率

**阶段 3：路线 2（GPU 全局步）**
- Chebyshev 半迭代 + 逐顶点块 Jacobi 预条件；固定迭代数 vs 残差判据对比
- 若迭代数过大 → 子空间预条件（SA'23）或几何多重网格
- 与路线 3 做同场景 benchmark，确定"何时该留在 CPU 求解"

**阶段 4：打磨**
- 摩擦、自碰撞、多分辨率/子空间加速、渲染与交互、回归基准固化

---

## 8. 待你确认的开放问题

1. **精度目标**：实时交互（<16ms，允许 float32 与近似）还是离线/高质量（可 float64 + 更多迭代）？这直接决定路线 2 的可行性边界。
2. **网格规模**：目标顶点数区间（10k / 100k / 1M）？决定 GPU 全局步是否必要。
3. **GPU 硬件**：是否支持 `VK_EXT_shader_atomic_float`、`shaderFloat64`？（需要的话我可以在你机器上枚举设备能力，只读探测，不动代码。）
4. **接触范围**：只做布料–刚体/布料–布料，还是要自碰撞？自碰撞会让 CCD/候选筛选的并行与 filter 构造难度上一个台阶。
5. **摩擦**：现在就要，还是先做无摩擦（无摩擦可大幅简化局部步）？
6. **是否引入 ipc-toolkit 作为依赖**（省大量 CCD/距离代码，但引入 Eigen 依赖与 CPU-only 内核），还是自研并只把它当 oracle？
7. **PD 的技术选型**：约束用弹簧/三角拉伸（经典 PD）还是 ARAP/子空间？是否接受"外层对真实能量做线搜索"（这会偏离纯 PD，但换来 IPC 的保证）。
