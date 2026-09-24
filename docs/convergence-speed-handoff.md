# 运行效率与公开数据的对照：接手说明

> 本文只记录**可以核实的材料、命令、参数、数值与操作过程**。不含成因判断，也不含任务安排。
> 全部为中文（本仓库约定）。

---

## 1. 当前待查的现象（陈述）

本项目的 PD（Projective Dynamics）迭代收敛速度，与公开论文报告的数字**不一致**。已有的一次性对照记录如下。

### 1.1 公开侧的数字（来源：Wang 2015, *Chebyshev Semi-iterative Approach for Accelerating Projective and Position-based Dynamics*）

| 来源 | 场景 | 方法 | 公开说法 |
|---|---|---|---|
| Fig. 9(b) | 桌布 10K 顶点 / 20K 单元 | 未加速 PD + 直接法 | ~8–10 次迭代到 1e-2 |
| Fig. 9(d) | 渔网 20K 顶点 | 同上 | ~10 次到 1e-2 |
| Table 1 | 同上例子 | 同上 | 每帧 10 次迭代；ρ = 0.9992–0.9999 |
| Fig. 3 | 刚度扫描（未标网格规模） | 直接法 | κ=1e4 时到 1e-3 约 1.3K 次；κ≥1e6 平台在 ~4e-2 |

读图方式：论文 PDF → `pdftoppm` 转图 → 人工读数（放大截图存于 `%TEMP%\pd_fig\`）。
论文只定义 `e(k) = ∇E(q^(k))`，**未写明纵轴的归一化方式**；其发布源码的 Direct 路径不调用 `Get_Error`；
Fig. 3 图注区分 "linear system 的误差" 与 "dynamical system 的误差"；正文称某例子"每步切 8 个子步"。

### 1.2 我们侧的数字（2 范数相对残差 `r_k/r_0`，`r_0` 取该子步起始）

> **口径声明（2026-09-24 复核补充）**：下表全部数字均在 **`PD_SHEAR=3000000`** 下测得；
> 不设该变量时同一命令给出的是**另一套数**（见 §3 的警告）。下表"锚点"列指出该行的归一化分母
> 用的是**我们锚点 r₀**（当前位形处）还是**论文锚点 e₀ = ∇E(预测子位形)**——
> 两者不是同一个量，混用会读错。

| 场景 | 相位 | 锚点 | 到 1e-1 | 到 1e-2 | 每步收缩率 |
|---|---|---|---|---|---|
| 100×100（10 000 顶点）κ=2305、h=1/120、顶边 pin | 瞬态（第 60 子步） | r₀ | 2 次 | 7 次 | 0.86 |
| 同上 | 稳态（第 600 子步） | r₀ | 未达到 | 未达到（48 次内） | 0.9977 |
| 60×60 | 瞬态 / 稳态 | r₀ | 4 次 / 未达到 | 17–18 次 / 未达到 | 0.83 / 0.9975 |
| 按论文参数复刻：101×101、s=0.01、m=1、κ=3e6、对角 3e6、h=1/30、两角 pin | **第 20 子步** | r₀ | 未达到 | 未达到（48 次内，r/r₀=0.166） | 0.938 |
| 同上 | **第 200 子步** | r₀ | 未达到 | 未达到（48 次内，r/r₀=0.151） | 0.942 |
| 同上 | **第 200 子步** | **e₀（论文锚点）** | **46 次** | 未达到（48 次内） | — |

加速后（`PD_CHEB=1`，ρ=0.98，S=10）：

| 场景 | 朴素 | 加速后 | 残差改善 |
|---|---|---|---|
| 60×60、300 子步、每子步 20 / 40 / 80 / 160 次迭代 | 12.33 / 8.59 / 6.14 / 2.65 | 9.38 / 7.41 / 3.52 / 0.684 | 1.31× / 1.16× / 1.74× / 3.87× |
| 硬度扫描（80 次迭代）κ=230.5 / 2305 / 23050 | 0.237 / 0.919 / 8.685 | 0.117 / 0.708 / 7.460 | 2.03× / 1.30× / 1.16× |
| 60×60 顶边 pin 稳态、子步内 `r_40/r_0` | 0.798 | 0.355 | 2.25× |

#### 1.3 已记录的对照与实验（文档位置）

`docs/pd-convergence.md`：

| 小节 | 内容 |
|---|---|
| §8.1–§8.3 | 上表的原始记录（公开侧 / 未加速 / 加速后） |
| §8.6 | 外部只读复核给出的两条数字修正与实验清单 |
| §8.7 | 两个已执行的对照实验（复刻算子的谱；"预测子位形处 ∇E" 的锚点测量）及其数值 |
| §8.8 | 两个已执行的对照实验（8–10 次窗口的几何平均；轨迹误差口径 K→err 表）及其数值 |
| §4.4b–§4.4d | 读他们源码得到的事实、复刻与 h 扫描、两种残差定义的数值对照 |
| §4.5 | 迭代算子的谱（`ρ(T)`、`λ_min(B)`、`λ_max(B)`、`ρ_cheb`）与实测收缩率的对照 |
| §3.3 | 刚度标定与预条件失配的**既有**记录（仓库原有内容） |

`docs/chebyshev.md`：§7.4 v1 实测；§7.5 真实区间法的两份实现与其数值结果。

---

## 2. 残差的计算方式（两侧定义）

### 2.1 我们侧

- 生产实现：`core/sim/Integrator.cpp` 的 `nonlinearResidual(ctx)`（返回每顶点 `|∇E|/m` 的 **max**）。
  定义式：`∇E = M/h²·(x − x̂) + Σ_c κ_c (A_c x − d_c)`，其中 `d_c = ℓ_c · unit(A_c x)` 为**当前位形**的投影，
  pinned 顶点跳过（不是未知量），零质量顶点跳过。
- 诊断工具 `_verify/residual_trace.cpp`（`pd_restrace`）**复算**每顶点残差，输出 **max 与 2 范数**两种，
  并把复算的 max 与生产 `nonlinearResidual` **逐位互校**（不一致则退出码 4）。
- 归一化口径：`r_k/r_0`，`r_0` = **该子步起始**（"解全局步之前"）的残差。
- 工具另打印两种 r₀ 与"以预测子位形处 ∇E 为分母"的门槛表（`E1` 行）。

### 2.2 他们侧（源码事实）

- `Chebyshev/lib/PROJECTIVE_MESH.h` 的 `Get_Error(TYPE t)`：
  `Error_i = (M_i + fixed_i)/t² · old_X_i + Σ_c κ_c · d_c − Σ_j A_ij X_j − A_ii X_i`，即 `b − A·X`。
- `old_X` 的赋值位置：`Chebyshev/lib/DYNAMIC_MESH.h` 的 `Begin_Constraints()`（`memcpy(old_X, X, ...)`），
  它由 `Update(TYPE t, int iterations)` 在**预测子 `Update(t)` 之后**调用。
- `Direct_Constraints(TYPE* next_X, TYPE t)` 的右端：`bx(i) = c·old_X[i]`，再对 `all_VL[index] != -1` 的邻居累加
  `d·(spring_k·all_VL[index]/|d|)`（`c = (M_i + fixed_i)/t²`）；随后对 x/y/z 三个标量系统各解一次。
- 数值对照记录：`docs/pd-convergence.md` §4.4d（60×60 稳态下 `b−A·x` 与我们的 `∇E` 的 max / 2 范数及其相对差）。

---

## 3. 复现我们的数字：工具与命令

构建：`.\tools\build.ps1`（Release，输出到 `build\Release\`）。
验收：`.\tools\check.ps1`（exit 0 = 通过；其第 ④ 步比对默认约束集的物理输出与 `build\_baseline\pd_bench_pre_fuse.exe`）。

> **命令里的可执行文件与参数由门禁 `check.ps1` 第 ⑥ 步自动反查**（2026-09-24 新增）：文档里出现的
> `pd_*.exe` 必须在 CMake 目标表里存在、且每个 `--flag` 必须在该程序接受的表内 ——
> 所以"改了程序名/参数名却忘了改文档"会直接红。**但它查不出"少写一句环境变量"**（例如漏 `$env:PD_SHEAR`），
> 那类只能靠上面的口径声明。

| 程序 | 作用 |
|---|---|
| `pd_restrace` | 子步内"残差–迭代"曲线；`--grid --spacing --density --stiffness --dt --damping --pin top\|corners\|soft-corners --steps --iters --settle-iters --traj [--traj-steps --traj-ref]` |
| `pd_itspectrum` | 迭代算子 `T = I − A⁻¹(M/h²+N)`、`B = I − T` 的谱与 `ρ_cheb`；`--grid --spacing --density --stiffness --dt --pin --steps --iters` |
| `pd_bench` | 端到端跑固定子步数并打印末次残差/应变/能量；`--grid Nx Ny --steps --iters --stiffness` |
| `pd_viewer` | 交互查看器；`--preset preview\|realtime\|accurate`、`--cheb` 类开关见 `docs/code-map.md` |
| `pd_sharefactor` / `pd_solvecomp` / `pd_bendaudit` / `pd_check` | 因子共享、分量拆分、弯曲、几何等自检 |

环境变量（默认关）：`PD_CHEB=1`（Chebyshev）、`PD_CHEB_RHO`（默认 0.98）、`PD_CHEB_START`（默认 10）、
`PD_CHEB_MODE=interval`（实验性，启用时打印警告）、`PD_CHEB_M/D`、`PD_SHEAR`（对角剪切刚度）、`PD_BEND`（弯曲刚度）、
`PD_AFFINITY=exclusive`（长跑必须绑定）、`PD_DEBUG_RESIDUAL=1`。

**复刻论文参数的完整命令**（**必须每条都在 `PD_SHEAR` 生效的环境下跑** —— 见下面那段警告）：

```powershell
$env:PD_SHEAR='3000000'    # 对角（剪切）刚度 3e6：§1.2 与 §8.7/§8.8 的**全部**数字都在这个设置下测得
# 1) 子步内残差–迭代曲线（论文 Fig.9 口径的可比读数）
.\build\Release\pd_restrace.exe  --grid 101 --spacing 0.01 --density 5128 --stiffness 3e6 --dt 0.0333333 --damping 0 --pin corners --steps 20 --iters 48
# 2) 同一算子的谱
.\build\Release\pd_itspectrum.exe --grid 101 --spacing 0.01 --density 5128 --stiffness 3e6 --dt 0.0333333 --pin corners --steps 20 --iters 200
# 3) 轨迹误差口径（K 次/子步 vs 参考解 512 次/子步）
.\build\Release\pd_restrace.exe  --grid 101 --spacing 0.01 --density 5128 --stiffness 3e6 --dt 0.0333333 --damping 0 --pin corners --traj --traj-steps 20 --traj-ref 512
# 4) 第 200 子步（"46 次到 1e-1"那一行）
.\build\Release\pd_restrace.exe  --grid 101 --spacing 0.01 --density 5128 --stiffness 3e6 --dt 0.0333333 --damping 0 --pin corners --steps 200 --iters 48
Remove-Item Env:\PD_SHEAR
```

> ⚠ **漏设 `PD_SHEAR` 会得到另一套完全不同的数**（2026-09-24 复核时实测）：不设时第 20 子步
> 得 e₀=4160.49 / r₀=6845.36 / |v|∞=3.667，设了得 **e₀=8685.79 / r₀=5874.58 / |v|∞=3.584**
> —— 即本文 §1.2 与 §8.7/§8.8 记的那一套。**剪切对角边是这批对照的唯一变量，参数量之前务必先确认它开着。**

质量换算（实测）：带对角边时本引擎内部顶点质量 ≈ **1.95 · ρ · s²**，
故 `--density 5128` 时 m ≈ 1（与论文每顶点 m=1 对齐）。
> 注：该系数**只在 `PD_SHEAR` 开启时有意义**（它依赖对角边的拓扑），且 `Scene.cpp` 刻意把剪切
> 加在**质量计算之后**，所以开关剪切**不改变质量**、只改变约束集。

工具口径提醒（复核时实测，容易误读）：
- `pd_restrace` 的 **`--iters` 只决定曲线画多长**；前面那几百个沉降子步跑的是 `--settle-iters`（默认 **40**）。
  所以第 1/2/3/4 条的"到子步末的位形"都由 40 次/子步决定，与 `--iters` 无关。
- 打印的门槛表有两张：**我们的锚点 r₀** 与 **论文锚点 e₀**。§8.2 里"46 次到 1e-1"那一条出自**后者**。

---

## 4. 材料位置

| 材料 | 位置 |
|---|---|
| 我们的文档 | `docs/pd-convergence.md`（§1–§8）、`docs/chebyshev.md`、`docs/perf.md`、`docs/open-issues.md`、`docs/code-map.md`、`_perf/HANDOFF.md` |
| 我们的诊断工具源码 | `_verify/residual_trace.cpp`、`_verify/iteration_spectrum.cpp`、`_verify/share_factor.cpp` |
| 论文文本（pdftotext 输出） | `%TEMP%\pd_fig\cheb2.txt` |
| 论文图像（放大截图） | `%TEMP%\pd_fig\`（`fig9_zoom.png`、`fig9_cd.png`、`fig3_zoom.png` 等） |
| 论文源码副本（未进本仓库） | `%TEMP%\chebcode\x\Chebyshev_sim\`（下载自论文页面的 `Wang-2015-CSI.zip`） |
| 论文源码的**可构建副本** | `%TEMP%\chebsim\`（见第 5 节的操作记录） |
| 历史报告 | `build\_perf\REPRO_REPORT.md`（含 §9–§15 的各轮记录） |

论文源码的许可为"仅教育/研究用途"，本仓库**未收录**其代码；`%TEMP%` 下的副本不进版本控制。

---

## 5. 他们源码的构建与运行记录（可复制的操作）

### 5.1 环境

- 工程文件：`Chebyshev/CPU_Projective_Square/new.sln` / `new.vcxproj`（`PlatformToolset=v120`，Win32）。
- 本机 MSVC：`D:\tools\Microsoft Visual Studio\2022\Community\`（v143）。
- 依赖：工程只声明 `glew32.lib`；zip 内含 `glut32.dll ×6`、`glew32.dll`、`glew32.lib`、`glew32s.lib`、内置 Eigen。
- **两次构建障碍与当时采用的做法**：
  1. 缺 `GL/glew.h` / `gl/glut.h` → 使用本机既有头：`D:\tools\anaconda3\envs\zsfzz\Library\include`（内含 `GL/glut.h`），
     并新建 shim `%TEMP%\chebsim\glincl\GL\glew.h`（内容为 `#include <GL/glut.h>`）。
  2. 内置 2015 版 Eigen 报 `scalar_product_traits<double,float>` 等错误 → 以 `CLOTHING<double>` 实例化。

### 5.2 无窗口运行方式（不启动他们的 OpenGL 驱动）

自建 `%TEMP%\chebsim\CPU_Projective_Square\my_main.cpp`：

```cpp
#include "CLOTHING.h"
#include <cstdio>
#include <cstdlib>
int main(int argc, char** argv) {
  const double h = 1.0/30.0;
  int steps = (argc > 1) ? atoi(argv[1]) : 3;
  int iters = (argc > 2) ? atoi(argv[2]) : 40;
  CLOTHING<double> cloth;
  cloth.gravity = -9.8;          // 该成员在发布源码中未见赋值语句
  cloth.Initialize(h);
  printf("VERTS %d\n", cloth.number);
  for (int s = 0; s < steps; ++s) { printf("STEP %d\n", s); cloth.Update(h, iters); }
  return 0;
}
```

编译脚本 `%TEMP%\chebsim\CPU_Projective_Square\build.bat`：

```bat
call "D:\tools\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars32.bat" >nul
cd /d "%~dp0"
cl /nologo /EHsc /O2 /I "..\glincl" /I "D:\tools\anaconda3\envs\zsfzz\Library\include" /I ".." /I "..\eigen" /I "..\lib" my_main.cpp /Fe:cheb_test.exe
```

### 5.3 对 `%TEMP%\chebsim\lib\PROJECTIVE_MESH.h` 做过的两处改动（仅该副本）

1. 在 `Update(TYPE t, int iterations)` 内，把被注释的 `Direct_Constraints(next_X, t);` 打开、把
   `Jacobi_Constraints(next_X, t);` 注释掉（即改为 Fig. 9(b) 的 Direct 分支）。
2. 在 `Direct_Constraints` 内 `Eigen::VectorXd z = solver.solve(bz);` 之后插入逐迭代日志：
   调用 `Get_Error(t)`，随后对 `Error[0 .. number*3-1]` 求 2 范数并 `printf("CHEBLOG %d %.10g\n", ++n, v)`。

---

## 6. 读他们源码得到的事实清单

| 项 | 值 / 事实 | 获取方式 |
|---|---|---|
| 顶点数 / 三角形数 | 10201 / 20000 | 运行 `CLOTHING<double>` 后打印 `number` / `t_number` |
| 边表长度 `e_number` | 30200 | 同上 |
| 质量 `M[i]` | 1（`Compute_Mass` / `Set_Mass` 全树无调用点） | 打印 |
| 固定方式 | `fixed[i] = 100000`，进 `(M+fixed)/h²`；`Update(t)` 中 `if (fixed[i]) continue;` | 打印 + 源码 |
| `spring_k` / `bending_k` | 3e6 / 1e4 | 打印 + `CLOTHING.h` |
| `rho` | 0.9999 | 打印 + `CLOTHING.h` |
| 网格 | `Make_A_Plane(101,101,-0.5,0.5)`，间距 0.01，平面 1×1 | `CLOTHING.h` |
| 邻接表 | 把弹簧邻居与弯曲邻居混排；`all_VV[0..5] = 1 101 102 0 2 3`（含对角 102） | 打印 |
| `all_VL[]` | `0.01 0.01 -1 0.01 0.01 -1` ⇒ 对角项为 `-1`，在右端被 `if (all_VL[index]==-1) continue;` 跳过 | 打印 |
| `all_VC[0..3]` | 9.60409e7 / 1.51109e7 / 9.1209e6 / 1.51209e7 | 打印 |
| `time_step` | 1/30（`OPENGL_DRIVER.h`：`float time_step=1/30.0;`） | 源码 |
| `gravity` | 全树（`*.h/*.cpp/*.cu`）未见赋值语句；显式设为 -9.8 后运行 | 源码检索 + 运行 |
| 发布版迭代 | `Update(t, iterations)` 内 `Jacobi_Constraints` 生效、`Direct_Constraints` 被注释 | 源码 |
| 发布版 Chbyshev 递推 | `l<=10` 时 `omega=1`；`l==11` 时 `2/(2−ρ²)`；其后 `4/(4−ρ²·ω)`；另有 `under_relax=1` | 源码 |
| 方布场景是否含接触 | 驱动 `OPENGL_DRIVER.h` 内无桌子/接触对象（仅有 `glutWireSphere` 画标记）；论文只在连衣裙例子提"每子步末处理碰撞与摩擦" | 源码检索 |

---

## 7. 当前状态（现象记录）

**他们代码在无窗口 harness 中**：可编译、可运行（10201 顶点、`Update(h, iters)` 正常返回）。

运行的中间量（`h=1/30`、`gravity=-9.8`、两角 pin、`CLOTHING<double>`）：

| 状态 | 应变 max | 应变 mean | `|V|max` | `|Error|max` |
|---|---|---|---|---|
| 初始（0 步） | 0 | 0 | 0 | 0 |
| 1 步、1 次迭代 | 0.1814 | 3.91e-4 | 0.302913 | **0** |
| 1 步、10 次迭代 | 0.04527 | 2.84e-4 | 0.330758 | **0** |
| 1 步、40 次迭代 | 0.02683 | 1.51e-4 | 0.354145 | **0** |

其中 `|Error|max` 恒为精确 0（`%.10g` 下逐次输出 `0`；对 `Error` 数组统计"非零元素个数" = 0）。

他们侧尚未取得可用于逐迭代对照的中间数据序列。

---

## 8. 仓库既有约定（与本次无关的常规门槛，供查阅）

- 文档全中文；`*.ps1` 必须带 UTF-8 BOM，其余跟踪文本 UTF-8 无 BOM + LF（`tools/check.ps1` 第 ③ 步）。
- 不读取邻居仓库 `dp_ipc`。
- 默认约束集的物理输出与存档基线逐字一致（`tools/check.ps1` 第 ④ 步）；改动默认路径需重新标定基线。
- 提交信息用中文；长提交信息写入临时文件后用 `git commit -F`。
