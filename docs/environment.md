# 开发环境记录（只读探查结果）

> 探查方式：PowerShell 只读查询 + 自建 Vulkan C++ 探针（`tools/vk_probe.cpp`）。
> 未修改任何系统配置、未安装任何软件。WMI/CIM 查询被沙箱策略拒绝，故部分信息改用注册表与 .NET API 获取。

---

## 1. 主机

| 项 | 值 |
|---|---|
| OS | Windows 10 Pro，build 10.0.26200 |
| CPU | 12th Gen Intel Core i7-12700F |
| 逻辑处理器 | **20**（i7-12700F 为 8 P-core（HT）+ 4 E-core；**线程池不能假设同构**） |
| 内存 | 未取到（CIM 拒绝访问）；GPU 报告 host-visible heap 63.91 GB，据此推断物理内存 ≥ 32 GB |
| GPU | **NVIDIA GeForce RTX 4070 Ti**，驱动 591.86，独显 |
| 显示驱动 Vulkan | device apiVersion **1.4.325**；loader 1.4.321 |

**对架构的影响**：20 线程里有 4 个 E-core。局部步是纯吞吐型负载，E-core 参与仍然有收益；但**线程数不应硬编码为 `ProcessorCount`**，并且要避免假设"20 个等价线程的线性加速"（基准表的加速比会低于 20x，属正常）。

---

## 2. GPU / Vulkan 能力（探针实测，`tools/vk_probe.cpp`）

### 2.1 特性

| 能力 | 结果 | 对项目的意义 |
|---|---|---|
| `shaderFloat64` | **YES** | RTX 4070 Ti（AD104）fp64 为 1/64 吞吐，**不作为热路径**，但可用于少量诊断内核或残差累加 |
| `shaderInt64` | YES | 定点归约路径可用 |
| `shaderFloat16` / `shaderInt8` | YES | 存储压缩可选（非必需） |
| `storageBuffer16BitAccess` | YES | 同上 |
| `shaderBufferInt64Atomics` / `shaderSharedInt64Atomics` | YES | **归约可用 64 位整数原子**（定点化求和） |
| `bufferDeviceAddress` | YES | 可用裸指针式传参，减少描述符数量 |
| `timelineSemaphore` / `hostQueryReset` | YES | 异步 CPU↔GPU 流水与计时 |
| **`VK_EXT_shader_atomic_float`** | **YES** | float 原子加可用（但**不依赖**它做正确性路径，见下） |
| **`VK_EXT_shader_atomic_float2`** | **NO** | **不可用**：没有 fp64 原子加、没有 float 原子 min/max |
| `VK_KHR_cooperative_matrix` / `VK_NV_cooperative_matrix` | YES | 与本项目无关（无矩阵乘负载） |
| `VK_KHR_synchronization2` / `maintenance4/5/6` / `dynamic_rendering` | YES | 可用现代同步 API，命令缓冲记录更简洁 |

### 2.2 Subgroup

| 项 | 值 |
|---|---|
| `subgroupSize` | **32**（NVIDIA 固定 32） |
| compute 阶段支持 | YES |
| 支持的运算 | BASIC / VOTE / ARITHMETIC / BALLOT / SHUFFLE / SHUFFLE_RELATIVE / CLUSTERED / QUAD / PARTITIONED_NV 全部 YES（mask `0x7ff`） |
| `quadOperationsInAllStages` | YES |

**结论**：`subgroupAdd`/`subgroupBallot` 可用 → **归约走 subgroup 树 + 每 workgroup 一个槽位**是最优路径；`VK_EXT_shader_atomic_float` 存在，可作为"快速路径"开关，但**正确性基线不依赖它**（因为 E-core 混构 CPU 侧同样不用原子加，两侧保持同一套语义，便于逐位对照）。

### 2.3 计算限制

| 限制 | 值 | 备注 |
|---|---|---|
| `maxComputeSharedMemorySize` | 48 KB | **不是 64 KB 以上**：SpMV 的行分块与共享内存暂存要按 48 KB 设计 |
| `maxComputeWorkGroupInvocations` | 1024 | |
| `maxComputeWorkGroupSize` | 1024×1024×64 | |
| `maxStorageBufferRange` | 4095 MB | 单 buffer 上限，250k 顶点仅需 ~10 MB，无压力 |
| `maxDescriptorSetStorageBuffers` | 1048576 | 描述符几乎无上限（配合 descriptor indexing） |
| `maxPerStageDescriptorStorageBuffers` | 1048576 | |
| `maxPushConstantsSize` | 256 B | 迭代参数（步长、Chebyshev 系数、元素数）走 push constant |
| `maxBoundDescriptorSets` | 32 | |

### 2.4 显存

| heap | 大小 | 属性 |
|---|---|---|
| heap 0 | **11.71 GB** | device local |
| heap 1 | **63.91 GB** | host visible（含 staging 与系统内存映射区） |

**结论**：显存充裕，无需做显存分页/流式；staging 双缓冲代价可忽略。

---

## 3. 工具链

| 工具 | 版本 | 位置 / 备注 |
|---|---|---|
| MSVC | 19.40.33811 (VS 2022 Community) | `D:\tools\Microsoft Visual Studio\2022\Community`（**非** `C:\Program Files`） |
| CMake | 3.30.0 | 在 PATH；另有 `D:\tools\cmake-3.20.2-windows-x86_64`（不用） |
| Ninja | **未安装** | → CMake 用 `Visual Studio 17 2022` 生成器，或自行放置 ninja.exe |
| vswhere | **未在标准位置** | 位于 `C:\Program Files (x86)\Microsoft Visual Studio\Installer\`（**不存在**）→ 构建脚本不要依赖它 |
| Windows SDK | 10.0.22621.0 | 头 `D:\Windows Kits\10\Include\10.0.22621.0`，库 `D:\Windows Kits\10\Lib\10.0.22621.0` |
| Vulkan SDK | **1.3.290.0** | `C:\VulkanSDK\1.3.290.0`；`VULKAN_SDK`/`VK_SDK_PATH` 均已设置 |
| glslc | shaderc v2023.8 (v2024.1-9-g3ac03b8) | `C:\VulkanSDK\1.3.290.0\Bin\glslc.exe` |
| glslangValidator | 11:14.3.0 | 同 Bin |
| Git | 2.45.0 | 可用 |
| Python | 3.12.7 | `D:\tools\anaconda3` |
| CUDA | 12.6 | 已安装（本项目不使用，仅备查） |
| MinGW g++ | **8.1.0** | 在 PATH，但**不支持 `-std=c++20`**（只到 c++2a）→ **不要用它**，统一 MSVC |

### 3.1 验证过的编译命令（可直接复用）

本机 PATH 缺少 `INCLUDE`/`LIB`（非开发者命令提示符），手工编译需要显式设置：

```powershell
$vs='D:\tools\Microsoft Visual Studio\2022\Community'; $kits='D:\Windows Kits\10'
$sdkVer=(Get-ChildItem "$kits\Include" -Directory | Sort-Object Name -Descending | Select-Object -First 1).Name
$env:INCLUDE="$vs\VC\Tools\MSVC\14.40.33807\include;$kits\Include\$sdkVer\ucrt;$kits\Include\$sdkVer\shared;$kits\Include\$sdkVer\um"
$env:LIB="$vs\VC\Tools\MSVC\14.40.33807\lib\x64;$kits\Lib\$sdkVer\ucrt\x64;$kits\Lib\$sdkVer\um\x64"
cl /nologo /EHsc /std:c++17 /O2 /I"$env:VULKAN_SDK\Include" tools\vk_probe.cpp `
   /Fe:tools\vk_probe.exe /link /LIBPATH:"$env:VULKAN_SDK\Lib" vulkan-1.lib
```

**踩坑记录（对新工程的构建脚本有直接价值）**：
1. `link` 找不到 `kernel32.lib` → 必须把 Windows SDK 的 `um\x64`、`ucrt\x64` 加进 `LIB`；
2. `cl /Fe:tools\x.exe` 在 `tools\` 目录不存在时**不报错也不产出 exe**，且 `/Fe` 的目录部分行为不可靠 → 构建脚本应显式 `mkdir` 并分步 `cl /c` + `link /OUT:`；
3. Vulkan 头里 `subgroupSize` 在 `VkPhysicalDeviceProperties`，`supportedStages/supportedOperations` 在 **`VkPhysicalDeviceSubgroupProperties`**（不在 `VkPhysicalDeviceVulkan12Features`）；16-bit storage 特性在 `VkPhysicalDevice16BitStorageFeatures`（不在 Vulkan12 结构里）；float16/int8 在 `VkPhysicalDeviceShaderFloat16Int8Features`。探针首版因此编译失败，已修正。

---

## 4. 依赖可用性

| 依赖 | 状态 | 位置 / 备注 |
|---|---|---|
| **Eigen 3.4.0** | **已有** | `D:\tools\Eigen\eigen-3.4.0`（含 `Eigen/Dense`、`Eigen/Sparse`）→ CMake 用 `Eigen3::Eigen` 或直接 `include_directories` |
| **OpenMP** | 可用 | MSVC `/openmp`（2.0 语义）或 `/openmp:llvm`；E-core 混构下**动态调度**更稳 |
| **GLFW 3.3.8** | 已有 | `D:\Vulkan\glfw-3.3.8.bin.WIN64`，含 `lib-vc2022`（匹配 MSVC）、`include` |
| **Vulkan** | 已有 | SDK 1.3.290.0；`find_package(Vulkan)` 可用 |
| glslc / glslangValidator | 已有 | shader 编译无缺口 |
| **VMA** (`vk_mem_alloc.h`) | **未找到** | SDK 内无；需从 GitHub 取单头文件（或先不用，直接 `vkAllocateMemory` 手写子分配） |
| **CHOLMOD / SuiteSparse** | **未安装** | 未找到 `cholmod.h`/`lib`；计划 D4/Q4 的"直接法预分解"**先用 Eigen `SimplicialLDLT`**，CHOLMOD 作为可选增量（需自行构建 SuiteSparse） |
| ImGui | 未找到 | M3 交互 UI 需取单仓（或先用极简自绘 HUD） |
| Catch2 / doctest | 未找到 | 取单头 doctest 即可（无构建负担） |
| Embree | **已有** | `D:\tools\embree`（`include`/`lib`/`bin`）—— 可选：CPU 侧 BVH/射线查询加速（布料–刚体距离查询），非必需 |
| Boost 1.86 | 已有 | 本项目不需要 |
| GLM | 已有 | `D:\tools\g-truc-glm-0af55cc`；可用但**推荐自写最小 Vec3/Mat3**（避免与 Eigen 的 ABI/类型混用） |
| vcpkg | 未安装 | 不引入包管理器，依赖全部走本地路径 |

---

## 5. 工程约定与禁止事项

### 5.1 【铁律】禁止参考 `dp_ipc` 工程
本项目**只隔离 `D:\dsh_workspace\dp_ipc` 这一个工程**（本项目的**前身**：同作者、同假设、同 PD 内核的早期形态）。对其他项目（`learn_cloth`、`about_houdini` 等）**不设限**——源码、文档、方案、参数都可以参考、引用与对照。针对 `dp_ipc` 的禁止项包括但不限于：

- 不读其源码（求解器、碰撞、约束、网格、渲染等任何模块）；
- 不采用其方案、推导、默认参数或"已验证结论"；
- 不引用其文档、调试记录或测试用例；
- 不与其对照，也不在其工作目录下构建或运行任何东西；
- 不在本项目的代码、注释、文档、提交信息中出现对该工程路径的引用或依赖。

**理由**：本项目要求从标准 PD 公式与文献出发、独立推导并逐项验证。以 `dp_ipc` 为参照，等于让一个已知含未修缺陷的早期版本（缺陷 A/B/C 见 `README.md` §4）替本轮推导背书，把"验证"退化为"复现"，并失去独立的第二意见。本项目的正确性来源只能是**数学推导 + 原语不变量断言 + 与解析解的对照**（见 `docs/plan.md` §2.2、M1）。

**执行方式**：规定写进 `docs/contributing.md` §1，并由 CI 扫描**仅针对 `dp_ipc` 的引用**（命中即失败）。旧规则写的是"扫描本工作区其他工程的路径引用"，范围过大且会误伤本项目自身路径，已收窄。构建配置仍一律使用本项目内的相对路径，保证本仓库可独立构建。

### 5.2 其他项目约定
- 依赖全部使用本机已安装的库（见 §4），不引入包管理器，不联网获取（除单头测试/UI 库一次获取外，均在本项目内记录来源）。
- 所有"已验证"的结论必须能指向一条可重跑的测试断言；不接受口头结论。

---

## 6. 对本计划的确认与修正

| 计划中的假设 | 探查结论 | 动作 |
|---|---|---|
| 方向 2 用 Chebyshev + 不依赖 float 原子 | 确认可行；`VK_EXT_shader_atomic_float` **存在**但 `_float2` 不存在 | 归约主线仍用 subgroup 树 + 定点 int64；float 原子仅作可选快速路径 |
| 48 KB 共享内存够用 | 确认 48 KB | SpMV 行分块按 48 KB 设计；不要把整个预条件子放进共享内存 |
| 需要 `shaderFloat64` 做 GPU 精确累加 | 硬件支持但吞吐 1/64 | 热路径 baseline 仍 float32；fp64 仅用于 GPU 侧诊断/残差校验内核 |
| CPU 侧 20 线程线性加速预期 | 20 逻辑核含 4 E-core | 基准表按 1/2/4/8/16 线程呈现，不承诺线性；线程数运行时可配 |
| Eigen 需获取 | 本地已有 3.4.0 | CMake 直接指向 `D:\tools\Eigen\eigen-3.4.0` |
| CHOLMOD 可用 | **未安装** | M2 的"预分解复用"先用 Eigen `SimplicialLDLT`；CHOLMOD 降级为可选项（Q4 需重新确认是否值得自行构建） |
| 渲染可用 GLFW | 本地有 GLFW 3.3.8 (vc2022) | M3 可直接用 GLFW；ImGui 需自行获取（渲染后端不受 §5.1 影响 —— 那里只锁 `dp_ipc`） |
| 编译器 | MSVC 19.40 + CMake 3.30 | 统一 MSVC；**禁用 PATH 上的 MinGW g++ 8.1** |

---

## 7. 复现方式

```powershell
# 重新生成能力报告（只读，不修改系统状态）
# 当前用 §3.1 的命令手工编译；M0 会提供 tools/build_probe.ps1 封装
cl /nologo /EHsc /std:c++17 /O2 /I"$env:VULKAN_SDK\Include" tools\vk_probe.cpp `
   /Fe:tools\vk_probe.exe /link /LIBPATH:"$env:VULKAN_SDK\Lib" vulkan-1.lib
.\tools\vk_probe.exe
```

探查脚本与探针源码均随工程保留：`tools/vk_probe.cpp`。若更换机器，重跑该探针即可刷新本文档第 2 节。
