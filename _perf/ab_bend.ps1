param(
  [int]$Reps = 5,
  [int]$Grid = 40,
  [int]$Iters = 40,
  [int]$Steps = 300,
  [int]$Threads = 4,
  [string]$Bend = "1000"
)

# _perf/ab_bend.ps1 -- 弯曲约束（PD_BEND）开/关的成对 A/B。
#
# 纪律（docs/perf.md §0）：**同一轮会话内交替、重复多次取最小值**。
# 本机读数在会话内会漂移（同一个二进制同工况测到过 155 ms → 123 ms），
# 单次对照能给出完全相反的结论。
#
# 注意：本文件在 build/ 下（被 gitignore），且**必须带 UTF-8 BOM** ——
# PowerShell 5.1 读无 BOM 的中文 .ps1 会按 GBK 解码而报语法错。

$ErrorActionPreference = "Stop"
$root = Split-Path -Parent $PSScriptRoot   # 本脚本在 <仓库根>/_perf/ 下（2026-09-24 从 _perf/ 迁入版本控制）
Set-Location $root
$exe = Join-Path $root "build\Release\pd_bench.exe"
$outDir = Join-Path $root "build\_baseline"
$log = Join-Path $outDir "bend_ab.log"

function Run-Bench([string]$bendValue) {
  $psi = New-Object System.Diagnostics.ProcessStartInfo
  $psi.FileName = $exe
  $psi.Arguments = "--grid $Grid $Grid --steps $Steps --iters $Iters --residual-tol 0.3 --threads $Threads"
  $psi.RedirectStandardOutput = $true
  $psi.UseShellExecute = $false
  $psi.EnvironmentVariables["PD_BEND"] = $bendValue
  $p = [System.Diagnostics.Process]::Start($psi)
  $t = $p.StandardOutput.ReadToEnd()
  $p.WaitForExit()
  return $t
}

function Parse-Stats([string]$text) {
  $res = [ordered]@{}
  # ⚠️ 必须**显式抓"单子步"后面那个数**（旧写法抓到的是 `总耗时` 的秒数，
  # 于是"单子步(ms)"整列差 3.33 倍；比值不受影响）—— 2026-09-23 发现并修。
  $stepLine = ($text -split "`n" | Select-String "子步数" | Select-Object -First 1).ToString()
  $mPer = [regex]::Match($stepLine, "单子步\s+([\d.]+)")
  if ($mPer.Success) { $res.perStep = [double]$mPer.Groups[1].Value }
  # 阶段桶与 nnz 逐项取出来（表格式输出，用正则抓数字最稳）
  $lines = $text -split "`n"
  foreach ($row in @(
    @("predict", "预测"),
    @("projScatter", "投影\+散射"),
    @("assemble", "组装\+分解"),
    @("solve", "全局回代"),
    @("unpack", "解包\+判据"),
    @("velocity", "速度更新"),
    @("other", "其它"))) {
    $m = $lines | Select-String ("^\s+" + $row[1]) | Select-Object -First 1
    if ($m) {
      $nums = [regex]::Matches($m.ToString(), "[\d]+\.\d+")
      if ($nums.Count -ge 1) { $res[$row[0] + "Ms"] = [double]$nums[0].Value }
    }
  }
  $m = $lines | Select-String "L 自由度" | Select-Object -First 1
  if ($m) {
    $mn = [regex]::Match($m.ToString(), "nnz\s+(\d+)")
    if ($mn.Success) { $res.lNnz = [int]$mn.Groups[1].Value }
    $mf = [regex]::Match($m.ToString(), "因子 nnz\s+(\d+)")
    if ($mf.Success) { $res.factorNnz = [int]$mf.Groups[1].Value }
  }
  $m = $lines | Select-String "符号分解.*数值分解.*回代" | Select-Object -First 1
  if ($m) {
    $nums = [regex]::Matches($m.ToString(), "\d+")
    if ($nums.Count -ge 3) { $res.solveCalls = [int]$nums[2].Value }
  }
  $res.lastResidual = ($lines | Select-String "末次残差" | Select-Object -First 1) -replace '.*末次残差\s+([\d.e+-]+).*','$1'
  $res.strainMax = ($lines | Select-String "最终应变:" | Select-Object -First 1) -replace '.*最大\s+([\d.e+-]+).*','$1'
  $res.energy = ($lines | Select-String "^能量:" | Select-Object -First 1) -replace '.*弹性\s+([\d.e+-]+).*','$1'
  return $res
}

$runs = @()
foreach ($rep in 1..$Reps) {
  Write-Host "rep $rep/$Reps ..." -ForegroundColor Cyan
  $offText = Run-Bench ""
  $onText = Run-Bench $Bend
  $runs += [pscustomobject]@{ rep = $rep; mode = "off"; stats = (Parse-Stats $offText); text = $offText }
  $runs += [pscustomobject]@{ rep = $rep; mode = "on"; stats = (Parse-Stats $onText); text = $onText }
}

$lines = New-Object System.Collections.Generic.List[string]
$lines.Add("=== 弯曲约束（PD_BEND=$Bend）开/关成对 A/B ===")
$lines.Add("工况：${Grid}x${Grid} / iters ${Iters} / $Steps 子步 / $Threads 线程 / 同轮交替 $Reps 次取最小值")
$lines.Add("时间：$(Get-Date -Format 'yyyy-MM-dd HH:mm:ss')")
$lines.Add("")
$lines.Add("各轮单子步（ms）：")
foreach ($rep in 1..$Reps) {
  $o = ($runs | Where-Object { $_.rep -eq $rep -and $_.mode -eq "off" }).stats.perStep
  $n = ($runs | Where-Object { $_.rep -eq $rep -and $_.mode -eq "on" }).stats.perStep
  $lines.Add(("  rep {0}: off {1:F3}  on {2:F3}  变化 {3:+0.0;-0.0}%" -f $rep, $o, $n, (100.0 * ($n - $o) / $o)))
}
$lines.Add("")
$keys = @("perStep", "predictMs", "projScatterMs", "assembleMs", "solveMs", "unpackMs", "velocityMs", "otherMs")
$names = @("单子步(ms)", "预测", "投影+散射", "组装+分解", "全局回代", "解包+判据", "速度更新", "其它")
$lines.Add("取最小值对照：")
$lines.Add(("{0,-14} {1,12} {2,12} {3,12}" -f "项", "关闭", "开启", "变化"))
for ($i = 0; $i -lt $keys.Count; ++$i) {
  $k = $keys[$i]
  $off = ($runs | Where-Object { $_.mode -eq "off" } | ForEach-Object { $_.stats[$k] } | Measure-Object -Minimum).Minimum
  $on = ($runs | Where-Object { $_.mode -eq "on" } | ForEach-Object { $_.stats[$k] } | Measure-Object -Minimum).Minimum
  $pct = if ($off -ne 0) { 100.0 * ($on - $off) / $off } else { 0.0 }
  $lines.Add(("{0,-14} {1,12:F3} {2,12:F3} {3,12:+0.0;-0.0}%" -f $names[$i], $off, $on, $pct))
}
$lines.Add("")
foreach ($mode in @("off", "on")) {
  $r = $runs | Where-Object { $_.mode -eq $mode } | Select-Object -First 1
  $lines.Add(("{0,-4}: L nnz = {1}  因子 nnz = {2}  回代次数 = {3}  末次残差 = {4}  最终应变 max = {5}  弹性能 = {6}" -f `
    $mode, $r.stats.lNnz, $r.stats.factorNnz, $r.stats.solveCalls, $r.stats.lastResidual, $r.stats.strainMax, $r.stats.energy))
}
[IO.File]::WriteAllText($log, ($lines -join "`r`n") + "`r`n", (New-Object System.Text.UTF8Encoding $false))
$lines | ForEach-Object { $_ }
Write-Host "日志：$log" -ForegroundColor Green
