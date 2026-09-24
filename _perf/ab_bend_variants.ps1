param(
  [int]$Reps = 5,
  [int]$Grid = 40,
  [int]$Iters = 40,
  [int]$Steps = 300,
  [int]$Threads = 4
)

# _perf/ab_bend_variants.ps1 -- 弯曲 stencil 的**采样变体**降本对照。
#
# 问题：弯曲的代价几乎全在消元填充上（40x40/k=1000：L nnz +78 %，因子 nnz +177 %，
# 单子步 +73 %，见 docs/perf.md §12）。填充由"每个未知量被多少约束耦合"决定，而弯曲
# stencil 的耦合半径是 2 ⇒ 少取 stencil 是直接打在填充上的手段。
#
# 本脚本对下面这些配置做**同一轮会话内交替、重复 N 次取最小**的对照（纪律见 docs/perf.md §0）：
#   · off            —— 不启用弯曲（基准）
#   · std/s1         —— 今天的行为（行列全取、步长 1）
#   · rows / cols    —— 单向（只取行 / 只取列）
#   · checker        —— 棋盘（每个中心只在行/列里取一个方向，(i+j) 奇偶决定）
#   · s2 / s3        —— 隔行采样（步长 2 / 3）
# 每一档同时跑"原刚度"与"按采样密度补偿后的刚度"（密度 1/s ⇒ 刚度 x s），
# 因为降本后的变体**等效抗弯刚度会变软**，只比时间不比响应是没有意义的对照。
#
# 注意：本文件在 build/ 下（被 gitignore），且**必须带 UTF-8 BOM** ——
# PowerShell 5.1 读无 BOM 的中文 .ps1 会按 GBK 解码而报语法错。

$ErrorActionPreference = "Stop"
$root = Split-Path -Parent $PSScriptRoot   # 本脚本在 <仓库根>/_perf/ 下（2026-09-24 从 _perf/ 迁入版本控制）
Set-Location $root
$exe = Join-Path $root "build\Release\pd_bench.exe"
$outDir = Join-Path $root "build\_baseline"
$log = Join-Path $outDir "bend_variants_ab.log"

# label, sampling, stride, stiffness（sampling = "" 表示不启用弯曲）
$configs = @(
  @{ label = "off";        sampling = "";        stride = 1; k = 0 },
  @{ label = "std/s1";     sampling = "standard"; stride = 1; k = 1000 },
  @{ label = "std/s2";     sampling = "standard"; stride = 2; k = 1000 },
  @{ label = "std/s3";     sampling = "standard"; stride = 3; k = 1000 },
  @{ label = "std/s2 x2k"; sampling = "standard"; stride = 2; k = 2000 },
  @{ label = "std/s3 x3k"; sampling = "standard"; stride = 3; k = 3000 },
  @{ label = "rows/s1";    sampling = "rows";     stride = 1; k = 1000 },
  @{ label = "rows/s1 x2k";sampling = "rows";     stride = 1; k = 2000 },
  @{ label = "cols/s1";    sampling = "cols";     stride = 1; k = 1000 },
  @{ label = "cols/s1 x2k";sampling = "cols";     stride = 1; k = 2000 },
  @{ label = "checker/s1"; sampling = "checker";  stride = 1; k = 1000 },
  @{ label = "checker/s1 x2k"; sampling = "checker"; stride = 1; k = 2000 },
  @{ label = "checker/s2"; sampling = "checker";  stride = 2; k = 1000 },
  @{ label = "checker/s2 x2k"; sampling = "checker"; stride = 2; k = 2000 }
)

function Run-Bench([hashtable]$cfg) {
  $psi = New-Object System.Diagnostics.ProcessStartInfo
  $psi.FileName = $exe
  $psi.Arguments = "--grid $Grid $Grid --steps $Steps --iters $Iters --residual-tol 0.3 --threads $Threads"
  $psi.RedirectStandardOutput = $true
  $psi.UseShellExecute = $false
  if ($cfg.k -gt 0) {
    $psi.EnvironmentVariables["PD_BEND"] = "$($cfg.k)"
    $psi.EnvironmentVariables["PD_BEND_SAMPLING"] = $cfg.sampling
    $psi.EnvironmentVariables["PD_BEND_STRIDE"] = "$($cfg.stride)"
  }
  $p = [System.Diagnostics.Process]::Start($psi)
  $t = $p.StandardOutput.ReadToEnd()
  $p.WaitForExit()
  return $t
}

function Parse-Stats([string]$text) {
  $res = [ordered]@{}
  $lines = $text -split "`n"
  # ⚠️ 必须**显式抓"单子步"后面那个数**（旧写法抓到的是 `总耗时` 的秒数，
  # 于是"单子步(ms)"整列差 3.33 倍；比值不受影响）—— 2026-09-23 发现并修。
  $stepLine = ($lines | Select-String "子步数" | Select-Object -First 1).ToString()
  $mPer = [regex]::Match($stepLine, "单子步\s+([\d.]+)")
  if ($mPer.Success) { $res.perStep = [double]$mPer.Groups[1].Value }
  foreach ($row in @(
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
  $m = $lines | Select-String "弯曲约束" | Select-Object -First 1
  if ($m) {
    $ms = [regex]::Match($m.ToString(), "stencil\s+(\d+)\s+条")
    if ($ms.Success) { $res.stencils = [int]$ms.Groups[1].Value } else { $res.stencils = 0 }
  } else {
    $res.stencils = 0
  }
  $m = $lines | Select-String "L 自由度" | Select-Object -First 1
  if ($m) {
    $mn = [regex]::Match($m.ToString(), "nnz\s+(\d+)")
    if ($mn.Success) { $res.lNnz = [int]$mn.Groups[1].Value }
    $mf = [regex]::Match($m.ToString(), "因子 nnz\s+(\d+)")
    if ($mf.Success) { $res.factorNnz = [int]$mf.Groups[1].Value }
  }
  $m = $lines | Select-String "末次残差" | Select-Object -First 1
  if ($m) { $res.lastResidual = [regex]::Match($m.ToString(), "([\d.eE+-]+)\s+m/s").Groups[1].Value }
  $m = $lines | Select-String "最终应变" | Select-Object -First 1
  if ($m) {
    $mm = [regex]::Match($m.ToString(), "最大\s+([\d.eE+-]+)\s+平均\s+([\d.eE+-]+)")
    if ($mm.Success) {
      $res.strainMax = [double]$mm.Groups[1].Value
      $res.strainAvg = [double]$mm.Groups[2].Value
    }
  }
  $m = $lines | Select-String "全程最大应变" | Select-Object -First 1
  if ($m) { $res.peakStrain = [double]([regex]::Match($m.ToString(), "([\d.eE+-]+)").Groups[1].Value) }
  $m = $lines | Select-String "^能量:" | Select-Object -First 1
  if ($m) {
    $me = [regex]::Match($m.ToString(), "弹性\s+([\d.eE+-]+)")
    if ($me.Success) { $res.elastic = [double]$me.Groups[1].Value }
  }
  return $res
}

$runs = @()
foreach ($rep in 1..$Reps) {
  Write-Host "rep $rep/$Reps ..." -ForegroundColor Cyan
  foreach ($cfg in $configs) {
    $text = Run-Bench $cfg
    $runs += [pscustomobject]@{ rep = $rep; label = $cfg.label; k = $cfg.k; stats = (Parse-Stats $text); text = $text }
  }
}

function Min-Stat([string]$label, [string]$key) {
  return ($runs | Where-Object { $_.label -eq $label } | ForEach-Object { $_.stats[$key] } |
          Measure-Object -Minimum).Minimum
}

$lines = New-Object System.Collections.Generic.List[string]
$lines.Add("=== 弯曲 stencil 采样变体的降本对照（PD_BEND_SAMPLING / PD_BEND_STRIDE）===")
$lines.Add("工况：${Grid}x${Grid} / iters ${Iters} / $Steps 子步 / $Threads 线程 / 同轮交替 $Reps 次取最小值")
$lines.Add("时间：$(Get-Date -Format 'yyyy-MM-dd HH:mm:ss')")
$lines.Add("")
$lines.Add("各轮单子步（ms）：")
$hdr = "{0,-16}" -f "配置"
foreach ($rep in 1..$Reps) { $hdr += "{0,10}" -f "rep$rep" }
$lines.Add($hdr)
foreach ($cfg in $configs) {
  $row = "{0,-16}" -f $cfg.label
  foreach ($rep in 1..$Reps) {
    $v = ($runs | Where-Object { $_.rep -eq $rep -and $_.label -eq $cfg.label }).stats.perStep
    $row += ("{0,10:F3}" -f $v)
  }
  $lines.Add($row)
}
$lines.Add("")
$lines.Add("取最小值对照（off = 不启用弯曲；单子步/各阶段为 ms；应变与能量为质量指标）")
$lines.Add("")
$lines.Add(("{0,-16} {1,-9} {2,7} {3,10} {4,10} {5,9} {6,9} {7,8} {8,9} {9,9} {10,9} {11,10}" -f `
  "配置", "刚度", "stencil", "单子步", "vs off", "回代", "投影散射", "组装", "解包", "应变max", "应变avg", "弹性能"))
$offPerStep = Min-Stat "off" "perStep"
foreach ($cfg in $configs) {
  $k = "off"; if ($cfg.k -gt 0) { $k = "$($cfg.k)" }
  $ps = Min-Stat $cfg.label "perStep"
  $pct = "n/a"
  if ($offPerStep -ne 0) { $pct = ("{0:+0.0;-0.0}%" -f (100.0 * ($ps - $offPerStep) / $offPerStep)) }
  $lines.Add(("{0,-16} {1,-9} {2,7} {3,10:F3} {4,10} {5,9:F3} {6,9:F3} {7,8:F2} {8,9:F2} {9,9:F5} {10,9:F5} {11,10:F4}" -f `
    $cfg.label, $k, (Min-Stat $cfg.label "stencils"), $ps, $pct, `
    (Min-Stat $cfg.label "solveMs"), (Min-Stat $cfg.label "projScatterMs"), `
    (Min-Stat $cfg.label "assembleMs"), (Min-Stat $cfg.label "unpackMs"), `
    (Min-Stat $cfg.label "strainMax"), (Min-Stat $cfg.label "strainAvg"), `
    (Min-Stat $cfg.label "elastic")))
}
$lines.Add("")
$lines.Add("规模与峰值：")
$lines.Add(("{0,-16} {1,10} {2,11} {3,12} {4,12}" -f "配置", "L nnz", "因子 nnz", "全程峰值应变", "末次残差"))
foreach ($cfg in $configs) {
  $lines.Add(("{0,-16} {1,10} {2,11} {3,12:F5} {4,12}" -f `
    $cfg.label, (Min-Stat $cfg.label "lNnz"), (Min-Stat $cfg.label "factorNnz"), `
    (Min-Stat $cfg.label "peakStrain"), (Min-Stat $cfg.label "lastResidual")))
}
[IO.File]::WriteAllText($log, ($lines -join "`r`n") + "`r`n", (New-Object System.Text.UTF8Encoding $false))
$lines | ForEach-Object { $_ }
Write-Host "日志：$log" -ForegroundColor Green
