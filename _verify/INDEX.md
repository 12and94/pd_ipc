# `_verify/` 程序索引：谁在验收、谁只是证据

> 为什么需要这份索引（2026-09-24 建）：本目录此前只有"源码本身"这一种索引维度，
> 于是**分不清哪个程序是验收集、哪个是排查期的证据、哪个连构建都没进**。
> `README.md` §4 末尾那句"其中一些的前提假设已经过期，仅供追溯"是被动的提示，
> 需要读者自己逐个读源码才能判断。实测代价：`_verify/residual_audit.cpp` 被
> `docs/pd-convergence.md` §3.1 作为"残差公式正确性的判定依据"点名引用，
> 但它**没有进 CMake 的构建表** —— 照文档复现会卡在"找不到 exe"。
>
> 本文件只做**标注**，不搬动、不删除任何程序：那 29 个文件是排查证据链，
> 按本仓库"历史结论保留并标注"的纪律应当留在原地。

## 0. 三态定义（按**构建状态**分，这是唯一能自动核对的口径）

| 标记 | 含义 | 判据 |
|---|---|---|
| ✅ **验收集** | `tools/check.ps1` 会跑它，非 0 退出即门禁失败 | 见 §1 |
| 🔷 **可构建 · 证据** | 进 CMake、能编能跑，但**不在门禁里**；用途是历史测量与对照 | `CMakeLists.txt` 有 `add_executable`，但 `check.ps1` 不跑 |
| ⬜ **未构建 · 遗留** | 有源码，**没有** `add_executable`；门禁第 ⑥ 步会对被文档引用的那几个报红 | 见 §3 |

抽查方式：`tools/check.ps1` 第 ⑥ 步会把所有文档里出现的 `xxx.exe` 与 CMake 目标表对账。

## 1. ✅ 验收集（9 个，全部进 CMake + 全部门禁跑）

| 程序 | 源码 | 作用 |
|---|---|---|
| `pd_check` | `check.cpp` | 端到端 8 项验收（闭式解 / 平衡 / pin / 自由落体 / 有界 / 分解 1 次 / 力符号） |
| `test_primitives` | `tests/primitives/` | 距离约束定义级不变量、并行散射、着色、三分量拆分（233 断言） |
| `test_spring_vertical` | `tests/chain/` | 两顶点弹簧逐步对照独立递推参考解（155 断言） |
| `test_convergence_criterion` | `tests/chain/` | 收敛判据的性质断言（15 断言） |
| `pd_trans` | `translation_invariance.cpp` | 平移不变性 |
| `pd_chain` | `chain_test.cpp` | 自由链 / 悬挂链与解析解对照 |
| `pd_diraudit` | `direction_audit.cpp` | 投影方向三项判据 |
| `pd_bendaudit` | `bend_audit.cpp` | 弯曲约束的代数审计（22 条） |
| `pd_solvecomp` | `solve_components.cpp` | 三分量拆分的结构自检 + 逐位比对 |

另外 `check.ps1` 还会跑三个**一致性自检**（不单独算验收集，但同样会红）：
`pd_sharefactor`（因子切片同构）、`pd_restrace`（复算残差与生产逐位一致）、`pd_bench --preset`（三档预设自报）。

## 2. 🔷 可构建 · 证据（17 个：进 CMake，不在门禁里）

按用途分组。**引用它们的数字时必须带工况与机时**（跨会话不可比，见 `docs/perf.md` 的约定）。
"文档"列写 `—` 的，表示**在本索引建立时没有任何文档提到该程序名**（用途只能读源码）；
写不出文档位置的，其"作用"列按源码注释归纳，**未经独立核对**。

### 2.1 收敛 / 谱 / 残差（文献对照与成因定位的主力）

| 程序 | 源码 | 作用 | 主要文档 |
|---|---|---|---|
| `pd_restrace` | `residual_trace.cpp` | 子步内"残差–迭代"曲线（论文 Fig.9 口径）、E1 论文锚点、`--traj` 轨迹误差口径 | `pd-convergence.md` §8.2/§8.7/§8.8、`convergence-speed-handoff.md` §3 |
| `pd_itspectrum` | `iteration_spectrum.cpp` | 迭代算子 `T`/`B` 的谱、`ρ_cheb`、两种残差定义的对照 | `pd-convergence.md` §4.4d/§4.5/§8.7 |
| `pd_std_pd` | `standard_pd.cpp` | "标准 PD"的参考写法对照 | 追溯 |
| `pd_one_step` | `one_step_trace.cpp` | 单步逐项追踪 | 追溯 |
| `pd_rhs` | `rhs_breakdown.cpp` | 右端项分解 | 追溯 |
| `pd_iterscan` | `iteration_scan.cpp` | 按迭代数的扫描 | `pd-convergence.md` §4.4 |
| `pd_steady` | `steady.cpp` | 稳态判定与读数 | `pd-convergence.md` §3.2 |
| `pd_dampstab` | `damping_stability.cpp` | **阻尼稳定性**（文档零提及，用途待确认） | — |

### 2.2 刚度标定与 κ 影响

| 程序 | 源码 | 作用 | 主要文档 |
|---|---|---|---|
| `pd_calib` | `stiffness_calibration.cpp` | 刚度标定（ε·κ 恒定性的来源） | `pd-convergence.md` §4.1 |
| `pd_kappa` | `kappa_effect.cpp` | κ 的影响扫描 | `pd-convergence.md` §4.4 |
| `pd_kappafree` | `kappa_free.cpp` | 自由端的 κ 效应 | 追溯 |

### 2.3 力 / 能量 / 质量

| 程序 | 源码 | 作用 | 主要文档 |
|---|---|---|---|
| `pd_forceaudit` | `force_audit.cpp` | 逐顶点力的分解 | `pd-convergence.md` §3.1 |
| `pd_fvs` | `force_vs_strain.cpp` | 力 vs 应变 | `pd-convergence.md` §4.1 |
| `pd_initaudit` | `init_audit.cpp` | 初始量（质量 / 静止长度 / 应变）审计 | 追溯 |

### 2.4 几何 / pin / 链

| 程序 | 源码 | 作用 | 主要文档 |
|---|---|---|---|
| `pd_pinpos` | `pin_positions.cpp` | pin 位置的处理 | 追溯 |
| `pd_springfree` | `spring_free.cpp` | 自由弹簧 | 追溯 |
| `pd_hang` | `hang_test.cpp` | 悬挂链 | `pd-convergence.md` §3.3 |

### 2.5 性能（求解器与内存）

| 程序 | 源码 | 作用 | 主要文档 |
|---|---|---|---|
| `pd_scale` | `scale_sweep.cpp` | 规模扫描 | 追溯 |

### 2.6 新增（2026-09-24，Phase 4c/4d 的对照）

| 程序 | 源码 | 作用 | 主要文档 |
|---|---|---|---|
| `pd_sharefactor` | `share_factor.cpp` | 三分量因子"只存一份 vs 各存一份" | `perf.md` §16 |
| `pd_solveaudit` | `solve_audit.cpp` | 全局步求解器的延迟 / 带宽性质、换排序对照 | `perf.md` §17 |

> ⚠ **上表"文档提及 0 次"的 6 个程序**（`pd_dampstab`、`pd_kappa`… 中凡标注"追溯"与"—"的）
> 在本文件建立时**没有任何文档提到它们的文件名**。它们不是孤儿（进 CMake、能跑），
> 但**用途只能靠读源码**。要处置它们时，二选一：
> ① 在对应文档里补一句"这个证据由 `xxx.exe` 产生"；② 在本文档的"作用"列补齐后即可。

## 3. ⬜ 未构建 · 遗留（3 个：有源码，**没有** `add_executable`）

| 程序 | 源码 | 被文档引用处 | 处置建议 |
|---|---|---|---|
| — | `residual_audit.cpp` | `pd-convergence.md` §3.1（**作为判定依据点名**） | 建议**补进 CMake**：它锁住"残差公式在三个解析位形上正确"这条已被引用的判据 |
| — | `residual_crosscheck.cpp` | `pd-convergence.md` §2.2 / §7 | 同上，或改为"手工编译"的明确标注 |
| — | `stiffness_calibration_record.cpp` | `pd-convergence.md` §4.1 | 同上 |

**现状**：`docs/pd-convergence.md` §7 给的是**手工编译命令**
（`cl /nologo /std:c++20 ... _verify\residual_audit.cpp ...`），
即"这三个不进构建表"是有意的做法，不是遗漏。本索引的作用是把这件事**显式化**，
使引用它们的文档不必依赖读者自己去发现编译方式不同。

## 4. 维护约定

- 新增 `_verify/*.cpp` 时：**要么**在 `CMakeLists.txt` 里加 `add_executable`，
  **要么**在本文 §3 记一行（连同引用它的文档位置）。
- 新增文档时若引用某个程序：确保程序的 `.exe` 名与 `CMakeLists.txt` 的目标名一致 ——
  `tools/check.ps1` 第 ⑥ 步会校验这一点。
- 本文件不记录"哪个程序被删了"：本目录按约定不删除历史证据。
