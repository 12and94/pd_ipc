# _perf/refresh_report.ps1
#
# 把性能看板的数据刷新到当前代码：
#   ① 线程扫描日志 build/_perf/logs/<工况>_t<线程>_r<轮次>.txt（4 工况 × 5 线程数 × 3 轮）
#      —— 参数与历史日志逐字保持一致（子步数 / iters / 残差容差），否则曲线不可比。
#      **循环顺序是"轮次在外、工况与线程数在内"（逐轮交替）**：同一轮里 5 个线程档背靠背测完，
#      会话内的机器漂移才会对各档公平。旧写法是"一个线程档连跑 2 次再换下一档"，
#      18 线程永远最后测 —— 漂移会系统性地偏袒/坑某一档，看板的"最优线程"标签因此不可信。
#      （2026-09-24 改；当时的症状：100×100 档标出"最优=8 线程"，而同会话交替复测是 4T 更快。）
#   ② Phase 4c 的 A/B 日志 build/_baseline/phase4c_ab.log —— 与 phase3_ab.log **同格式**，
#      所以报告里的解析器（make_report.js 的 parsePhase3）可以直接复用。
#
# 用法： .\_perf\refresh_report.ps1            # 默认每例交替 5 次取最小值
#        .\_perf\refresh_report.ps1 -Reps 3
#        .\_perf\refresh_report.ps1 -SweepOnly  # 只重采 ① 扫描 + ④ 四套约束集（②③ 不动）
#
# 纪律（都是踩坑换来的，见 docs/perf.md §0 与 §6.2）：
#   · 交替测量、重复取最小值（单次运行波动可达 2×）；
#   · 每张表的条件（运行时 / 线程数 / 工况）必须写清楚；
#   · 本脚本自己核对两个二进制链接的是 vcomp 而不是 libomp。
#
# ⚠️ 本文件在 build/ 下（已被 gitignore），但**必须带 UTF-8 BOM** —— PowerShell 5.1 读
#    无 BOM 的中文 .ps1 会按 GBK 解码而报语法错。改完请复查前 3 字节是 EF BB BF。
#   `-SweepOnly`：只重采 ① 线程扫描与 ④ 四套约束集，**不动** ②③ 的 A/B 日志。
#   改扫描口径后用这一档：既保证 ① 与 ④ 仍是同一轮会话（看板 §二 与 §十 要对得上），
#   又不把 §八/§九/§十一 的 A/B 面板拖进新会话（那些面板靠"同轮交替取最小"自证，重采只会让
#   文档里记的百分比与看板脱钩 —— 分工是"文档记录当时那一轮、看板显示最近一轮"）。
param([int]$Reps = 5, [int]$SweepReps = 3, [switch]$SkipSweep, [switch]$SweepOnly)

$ErrorActionPreference = 'Stop'
$root = Split-Path -Parent $PSScriptRoot   # 本脚本在 <仓库根>/_perf/ 下（2026-09-24 从 _perf/ 迁入版本控制）
Set-Location $root

# 全部用**绝对路径**：PowerShell 的当前位置与 .NET 的 [Environment]::CurrentDirectory
# 是两回事，[IO.File]::WriteAllText 用相对路径会写到别处去。
$NEW = Join-Path $root 'build\Release\pd_bench.exe'
$OLD = Join-Path $root 'build\_baseline\pd_bench_pre_comp.exe'
$LOGDIR = Join-Path $root 'build\_perf\logs'
$ABLOG = Join-Path $root 'build\_baseline\phase4c_ab.log'
$UTF8 = New-Object Text.UTF8Encoding($false)   # 显式无 BOM：看板的 JS 按 UTF-8 读

if (-not (Test-Path $OLD)) { throw "缺少基线二进制 $OLD（本次改造的对照；见 build/_baseline/）" }
if (-not (Test-Path $LOGDIR)) { New-Item -ItemType Directory -Path $LOGDIR | Out-Null }

# ---- 运行时核对：必须都是 vcomp（踩过 CMake 缓存残留 /openmp:llvm 的坑）----
foreach ($exe in @($NEW, $OLD)) {
  $bytes = [IO.File]::ReadAllBytes($exe)
  $ascii = [Text.Encoding]::ASCII.GetString($bytes)
  $n = ([regex]::Matches($ascii, 'libomp')).Count
  if ($n -ne 0) { throw "$exe 里出现了 $n 次 libomp —— 不是 vcomp，先查 CMake 缓存！" }
  "运行时核对：$exe ⇒ vcomp（libomp 命中 0 次）"
}

# ---- 工况参数：与 build/_perf/logs 里的历史日志逐字一致 ----
$workloads = @(
  @{ id = '40x40_it40';   args = @('--grid', '40', '40', '--steps', '300', '--iters', '40', '--residual-tol', '0.3') },
  @{ id = '60x60_it2';    args = @('--grid', '60', '60', '--steps', '300', '--iters', '2', '--residual-tol', '0.3') },
  @{ id = '100x100_it10'; args = @('--grid', '100', '100', '--steps', '20', '--iters', '10', '--residual-tol', '0.001') },
  @{ id = '200x200_it10'; args = @('--grid', '200', '200', '--steps', '8', '--iters', '10', '--residual-tol', '0.001') }
)

function RunBench($exe, $a) {
  return ((& $exe @a 2>&1) -join "`n")
}

# 从一次运行输出里取：单子步 / 平均迭代 / 散射 / 解包+判据 / 全局回代
function Metrics($txt) {
  $head = [regex]::Match($txt, '单子步\s+([\d.]+)\s+ms\s+平均迭代\s+([\d.]+)')
  $st = @{}
  foreach ($line in ($txt -split "`n")) {
    $m = [regex]::Match($line, '^\s*(\S+)\s+([\d.]+)\s+[\d.]+%')
    if ($m.Success) { $st[$m.Groups[1].Value] = [double]$m.Groups[2].Value }
  }
  return [pscustomobject]@{
    total   = [double]$head.Groups[1].Value
    iter    = [double]$head.Groups[2].Value
    scatter = [double]$st['散射']
    unpack  = [double]$st['解包+判据']
    solve   = [double]$st['全局回代']
  }
}

# ---- ① 线程扫描日志 ----
# 逐轮交替：**轮次在最外层**，工况与线程数在内层。这样"取最小值"是跨轮取，
# 而每一轮里 5 个线程档的测量时刻彼此接近 ⇒ 漂移不再系统性地偏袒某个档位。
if ($SkipSweep) {
  '=== ① 线程扫描：已按 -SkipSweep 跳过 ==='
} else {
  # 先删掉本目录里旧的扫描日志（只匹配 <工况>_t<线程>_r<轮次>.txt 这一种命名）：
  # 换了轮次数之后，上一轮的残留文件会混进报告侧的"取最小值"，而它来自另一个会话。
  $stale = Get-ChildItem -Path $LOGDIR -Filter '*_t*_r*.txt' -File -ErrorAction SilentlyContinue
  if ($stale) { $stale | Remove-Item -Force; "  清理旧扫描日志 $($stale.Count) 个" }
  "=== ① 线程扫描（$($workloads.Count) 工况 × 5 线程 × $SweepReps 轮，逐轮交替）==="
  foreach ($r in 1..$SweepReps) {
    foreach ($w in $workloads) {
      foreach ($t in 1, 2, 4, 8, 18) {
        $a = $w.args + @('--threads', "$t")
        $txt = RunBench $NEW $a
        [IO.File]::WriteAllText((Join-Path $LOGDIR "$($w.id)_t${t}_r${r}.txt"), $txt, $UTF8)
      }
    }
    "  第 $r 轮完成（$($workloads.Count) 工况 × 5 档）"
  }
}

# ---- ② Phase 4c 的 A/B 日志（同一会话内交替、各取最小值）----
if ($SweepOnly) {
  '=== ② Phase 4c A/B：已按 -SweepOnly 跳过（看板 §六 保留它自己那一轮）==='
} else {
"=== ② Phase 4c A/B（每例交替各 $Reps 次取最小）==="
$cases = @(
  @{ name = '40x40_it40_t1';   args = @('--grid', '40', '40', '--steps', '300', '--iters', '40', '--residual-tol', '0.3', '--threads', '1') },
  @{ name = '40x40_it40_t4';   args = @('--grid', '40', '40', '--steps', '300', '--iters', '40', '--residual-tol', '0.3', '--threads', '4') },
  @{ name = '40x40_it40_t18';  args = @('--grid', '40', '40', '--steps', '300', '--iters', '40', '--residual-tol', '0.3', '--threads', '18') },
  @{ name = '60x60_it2_t4';    args = @('--grid', '60', '60', '--steps', '300', '--iters', '2', '--residual-tol', '0.3', '--threads', '4') },
  @{ name = '100x100_it10_t4'; args = @('--grid', '100', '100', '--steps', '20', '--iters', '10', '--residual-tol', '0.001', '--threads', '4') },
  @{ name = '200x200_it10_t4'; args = @('--grid', '200', '200', '--steps', '8', '--iters', '10', '--residual-tol', '0.001', '--threads', '4') }
)

$rows = @()
foreach ($c in $cases) {
  $bv = 1e30; $bs = 1e30; $bu = 1e30; $bp = 1e30; $bi = -1
  $nv = 1e30; $ns = 1e30; $nu = 1e30; $np = 1e30; $ni = -1
  for ($r = 0; $r -lt $Reps; $r++) {
    $o = Metrics (RunBench $OLD $c.args)
    if ($o.total -lt $bv) { $bv = $o.total; $bs = $o.scatter; $bu = $o.unpack; $bp = $o.solve; $bi = $o.iter }
    $n = Metrics (RunBench $NEW $c.args)
    if ($n.total -lt $nv) { $nv = $n.total; $ns = $n.scatter; $nu = $n.unpack; $np = $n.solve; $ni = $n.iter }
  }
  $delta = ($nv / $bv - 1) * 100
  $rows += ('{0,-22}{1,9:F3}{2,9:F3}{3,8:F1}%{4,8:F2}{5,8:F2}{6,10:F1}{7,10:F1}{8,10:F1}{9,10:F1}{10,10:F1}{11,10:F1}' -f `
             $c.name, $bv, $nv, $delta, $bi, $ni, $bs, $ns, $bu, $nu, $bp, $np)
  "  $($c.name)：$([Math]::Round($bv,3)) → $([Math]::Round($nv,3)) ms   $([Math]::Round($delta,1))%"
}

$out = @()
$out += 'Phase 4c A/B：全局步按 x/y/z 三分量并行（L = Ã ⊗ I₃） vs 整趟串行回代'
$out += "基线 pd_bench_pre_comp.exe / 新版 pd_bench.exe，每例交替各 $Reps 次，取最小值（ms/子步）"
$out += ''
$out += '用例                        基线        新版       变化   基线迭代   新版迭代     基线散射     新版散射     基线解包     新版解包     基线回代     新版回代'
$out += $rows
[IO.File]::WriteAllText($ABLOG, ($out -join "`r`n") + "`r`n", $UTF8)
"=== 写入 $ABLOG ==="
}

# ---- ③ 可选物理（剪切 / 弯曲）与弯曲采样变体的日志 ----
# 看板的 §八/§九/§十一 直接读这三个日志，这里顺手刷新，免得"看板里的数字是上个月的"。
# 三个脚本各自内部就是"同轮交替、取最小值"，所以直接用它们的默认 Reps。
if ($SweepOnly) {
  '=== ③ 附加约束与采样变体的 A/B 日志：已按 -SweepOnly 跳过 ==='
} else {
"=== ③ 附加约束与采样变体的 A/B 日志（每例交替各 $Reps 次取最小）==="
foreach ($s in @(
  @{ script = 'ab_shear.ps1';         log = 'shear_ab.log' },
  @{ script = 'ab_bend.ps1';          log = 'bend_ab.log' },
  @{ script = 'ab_bend_variants.ps1'; log = 'bend_variants_ab.log' })) {
  $p = Join-Path $PSScriptRoot $s.script
  if (Test-Path $p) {
    & $p -Reps $Reps | Out-Null
    "  $($s.script) ⇒ build\_baseline\$($s.log)"
  } else {
    "  跳过 $($s.script)（不存在）"
  }
}
}

# ---- ④ 四套约束集的阶段构成（看板 §十）----
# 默认 / 只剪切 / 只弯曲 / 两者都开，各跑 3 次取单子步最小值；日志进 build/_perf/logs_cs/
# （原样的 pd_bench 输出，看板侧直接用 parseBench 解析）。
"=== ④ 四套约束集的阶段构成（3 工况 × 4 套 × 3 次）==="
$cs = Join-Path $PSScriptRoot 'constraint_sets.ps1'
if (Test-Path $cs) {
  & $cs -Reps 3 | Select-Object -Last 1
} else {
  '  跳过 constraint_sets.ps1（不存在）'
}
