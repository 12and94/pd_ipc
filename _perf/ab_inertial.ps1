param(
  [int]$Reps = 5,
  [int]$Grid = 40,
  [int]$Iters = 40,
  [int]$Steps = 300
)

# _perf/ab_inertial.ps1 -- "惯性右端并行化"（omp single → omp for）的成对 A/B。
#
# ⚠️ **结论：否（2026-09-23）**，实测无收益，源码已回退。见 docs/perf.md §14。
#    ⇒ **现在跑本脚本，两端是同一份代码，差值应 ≈ 0**；要重做实验需先把
#      `core/sim/Integrator.cpp` 里那一步改回 `assembleInertialRhsInRegion()`（§14 有说明）。
#    决定性证据是"1 线程"这一档：单线程时 for 与串行循环等价，却已经测出 −1.7 %，
#    说明那点差别来自代码生成/布局，不是并行；去掉这个地板后 4T/18T 只剩 −1.2 %/−0.6 %。
#
# 背景：一个子步里，惯性右端 bBase = (M/h²)x̂ 原本是 `omp single`（一条线程干、其余空转等）。
# 本实验把它改成 `omp for`（按顶点切给所有线程）。正确性上没有风险：
# 逐顶点只写自己、无归并、算术逐字未变 ⇒ 与串行版**逐位相同**（已单独验证）。
# 同步点数也不变（原来靠 single 的隐式 barrier，现在靠 for 的隐式 barrier）。
# 所以唯一变量就是"这一步由几条线程干" —— 收益/代价能干净归因。
#
# 纪律（docs/perf.md §0）：**同一轮会话内交替、重复 N 次取最小值**。
# 线程数 1 是**"机制失效"对照**：单线程时 `omp for` 与串行等价，理论上差值应为 0 ——
# 若这里就测出差别，说明这点差别与并行无关（本实验正是靠它识破了 1.7 % 的布局抖动）。
#
# 注意：本文件在 build/ 下（被 gitignore），且**必须带 UTF-8 BOM** ——
# PowerShell 5.1 读无 BOM 的中文 .ps1 会按 GBK 解码而报语法错。

$ErrorActionPreference = "Stop"
$root = Split-Path -Parent $PSScriptRoot   # 本脚本在 <仓库根>/_perf/ 下（2026-09-24 从 _perf/ 迁入版本控制）
Set-Location $root
$exeOld = Join-Path $root "build\_baseline\pd_bench_pre_inertpar.exe"   # 串行（single）
$exeNew = Join-Path $root "build\Release\pd_bench.exe"                  # 并行（for）
$log = Join-Path $root "build\_baseline\inertial_ab.log"

function Run-Bench([string]$exe, [int]$threads) {
  $psi = New-Object System.Diagnostics.ProcessStartInfo
  $psi.FileName = $exe
  $psi.Arguments = "--grid $Grid $Grid --steps $Steps --iters $Iters --residual-tol 0.3 --threads $threads"
  $psi.RedirectStandardOutput = $true
  $psi.UseShellExecute = $false
  $p = [System.Diagnostics.Process]::Start($psi)
  $t = $p.StandardOutput.ReadToEnd()
  $p.WaitForExit()
  return $t
}

function Parse-Stats([string]$text) {
  $res = [ordered]@{}
  $lines = $text -split "`n"
  # ⚠️ 必须**显式抓"单子步"后面那个数**。旧写法是"在子步数行里取第一个 x.y"，
  # 抓到的是 `总耗时 0.802 s` 里的 **0.802（秒）**，于是"单子步(ms)"整列标着秒数
  # （差 1000/300 = 3.33 倍）。比值不受影响、绝对值全错 —— 2026-09-23 发现并修。
  $stepLine = ($lines | Select-String "子步数" | Select-Object -First 1).ToString()
  $mPer = [regex]::Match($stepLine, "单子步\s+([\d.]+)")
  if ($mPer.Success) { $res.perStep = [double]$mPer.Groups[1].Value }
  $mWall = [regex]::Match($stepLine, "总耗时\s+([\d.]+)")
  if ($mWall.Success) { $res.wallSeconds = [double]$mWall.Groups[1].Value }
  foreach ($row in @(
    @("predictMs", "预测"),
    @("projScatterMs", "投影\+散射"),
    @("assembleMs", "组装\+分解"),
    @("solveMs", "全局回代"),
    @("unpackMs", "解包\+判据"),
    @("velocityMs", "速度更新"),
    @("otherMs", "其它"))) {
    $m = $lines | Select-String ("^\s+" + $row[1]) | Select-Object -First 1
    if ($m) {
      $nums = [regex]::Matches($m.ToString(), "[\d]+\.\d+")
      if ($nums.Count -ge 1) { $res[$row[0]] = [double]$nums[0].Value }
    }
  }
  $res.lastResidual = (($lines | Select-String "末次残差" | Select-Object -First 1) -replace '.*末次残差\s+([\d.eE+-]+).*','$1')
  $res.strainMax = (($lines | Select-String "最终应变:" | Select-Object -First 1) -replace '.*最大\s+([\d.eE+-]+).*','$1')
  return $res
}

$cases = @(1, 4, 18)
$runs = @()
foreach ($threads in $cases) {
  foreach ($rep in 1..$Reps) {
    Write-Host "threads=$threads rep $rep/$Reps ..." -ForegroundColor Cyan
    $o = Parse-Stats (Run-Bench $exeOld $threads)
    $n = Parse-Stats (Run-Bench $exeNew $threads)
    $runs += [pscustomobject]@{ threads = $threads; rep = $rep; mode = "single"; stats = $o }
    $runs += [pscustomobject]@{ threads = $threads; rep = $rep; mode = "for"; stats = $n }
  }
}

function Min-Stat([int]$threads, [string]$mode, [string]$key) {
  return ($runs | Where-Object { $_.threads -eq $threads -and $_.mode -eq $mode } |
          ForEach-Object { $_.stats[$key] } | Measure-Object -Minimum).Minimum
}

$lines = New-Object System.Collections.Generic.List[string]
$lines.Add("=== 惯性右端并行化（omp single → omp for）成对 A/B ===")
$lines.Add("工况：${Grid}x${Grid} / iters ${Iters} / $Steps 子步 / 同轮交替 $Reps 次取最小值")
$lines.Add("基线 exe：build\_baseline\pd_bench_pre_inertpar.exe（single）")
$lines.Add("实验 exe：build\Release\pd_bench.exe（for）")
$lines.Add("时间：$(Get-Date -Format 'yyyy-MM-dd HH:mm:ss')")
foreach ($threads in $cases) {
  $lines.Add("")
  $lines.Add("--- $threads 线程 ---")
  $row = "{0,-10}" -f "各轮"
  foreach ($rep in 1..$Reps) { $row += "{0,10}" -f "rep$rep" }
  $lines.Add($row)
  foreach ($mode in @("single", "for")) {
    $row = "{0,-10}" -f $mode
    foreach ($rep in 1..$Reps) {
      $v = ($runs | Where-Object { $_.threads -eq $threads -and $_.rep -eq $rep -and $_.mode -eq $mode }).stats.perStep
      $row += ("{0,10:F3}" -f $v)
    }
    $lines.Add($row)
  }
  $so = Min-Stat $threads "single" "perStep"
  $sn = Min-Stat $threads "for" "perStep"
  $pct = if ($so -ne 0) { 100.0 * ($sn - $so) / $so } else { 0.0 }
  $lines.Add(("取最小：single {0:F3} ms   for {1:F3} ms   变化 {2:+0.0;-0.0}%" -f $so, $sn, $pct))
  $keys = @("predictMs", "projScatterMs", "assembleMs", "solveMs", "unpackMs", "velocityMs", "otherMs")
  $names = @("预测", "投影+散射", "组装+分解", "全局回代", "解包+判据", "速度更新", "其它")
  for ($i = 0; $i -lt $keys.Count; ++$i) {
    $a = Min-Stat $threads "single" $keys[$i]
    $b = Min-Stat $threads "for" $keys[$i]
    $p2 = if ($a -ne 0) { 100.0 * ($b - $a) / $a } else { 0.0 }
    $lines.Add(("    {0,-12} single {1,10:F3}   for {2,10:F3}   {3:+0.0;-0.0}%" -f $names[$i], $a, $b, $p2))
  }
  $lines.Add(("    物理校验：末次残差 {0} vs {1}；应变 max {2} vs {3}" -f `
    (Min-Stat $threads "single" "lastResidual"), (Min-Stat $threads "for" "lastResidual"), `
    (Min-Stat $threads "single" "strainMax"), (Min-Stat $threads "for" "strainMax")))
}
[IO.File]::WriteAllText($log, ($lines -join "`r`n") + "`r`n", (New-Object System.Text.UTF8Encoding $false))
$lines | ForEach-Object { $_ }
Write-Host "日志：$log" -ForegroundColor Green
