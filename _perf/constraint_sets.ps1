param(
  [int]$Reps = 2
)

# _perf/constraint_sets.ps1 -- 四套约束集下的"阶段构成"对照。
#
# 看板的第二节（分函数层级）画的是**默认约束集**（只有距离约束 + pin）的阶段构成。
# 本脚本把同一件事对另外三套配置也采一遍，于是可以回答：
#   · 只加剪切，时间挪到哪一块了？
#   · 只加弯曲，时间挪到哪一块了？
#   · 两个都加，是简单相加还是有交叉效应（同一批填充/同一批 barrier）？
#
# ⚠️ 四套配置的**物理不同**（约束集不同），所以这里只比"时间与构成"，
#    绝不比应变/能量。每个 (工况, 配置) 各跑 $Reps 次、报告取最小值。
#
# 输出：build/_perf/logs_cs/<工况id>__<配置id>_r<轮次>.txt（原样的 pd_bench 输出，
#       所以看板侧的 parseBench 可以原样复用，不需要新解析器）。
#
# 用法： .\_perf\constraint_sets.ps1 -Reps 2
#
# 注意：本文件在 build/ 下（被 gitignore），且**必须带 UTF-8 BOM** ——
# PowerShell 5.1 读无 BOM 的中文 .ps1 会按 GBK 解码而报语法错。

$ErrorActionPreference = 'Stop'
$root = Split-Path -Parent $PSScriptRoot   # 本脚本在 <仓库根>/_perf/ 下（2026-09-24 从 _perf/ 迁入版本控制）
Set-Location $root
$exe = Join-Path $root 'build\Release\pd_bench.exe'
$outDir = Join-Path $root 'build\_perf\logs_cs'
$UTF8 = New-Object Text.UTF8Encoding($false)
if (-not (Test-Path $outDir)) { New-Item -ItemType Directory -Path $outDir | Out-Null }

# 刚度取值与看板 §八/§九 的 A/B 保持一致，这样"阶段构成"与"开/关代价"两处能互相印证。
$SHEAR = '10000'
$BEND = '1000'

$configs = @(
  @{ id = 'none';       shear = '';      bend = '';      label = '默认（距离 + pin）' },
  @{ id = 'shear';      shear = $SHEAR;  bend = '';      label = '只加剪切' },
  @{ id = 'bend';       shear = '';      bend = $BEND;   label = '只加弯曲' },
  @{ id = 'both';       shear = $SHEAR;  bend = $BEND;   label = '剪切 + 弯曲' }
)

# 工况参数与 refresh_report.ps1 里的线程扫描**逐字一致**（这样 §二 与本节可比）。
$workloads = @(
  @{ id = '40x40_it40';   args = @('--grid', '40', '40', '--steps', '300', '--iters', '40', '--residual-tol', '0.3') },
  @{ id = '60x60_it2';    args = @('--grid', '60', '60', '--steps', '300', '--iters', '2', '--residual-tol', '0.3') },
  @{ id = '100x100_it10'; args = @('--grid', '100', '100', '--steps', '20', '--iters', '10', '--residual-tol', '0.001') }
)
$THREADS = 4   # 默认档

function Run-Cfg($args2, [string]$shear, [string]$bend) {
  $psi = New-Object System.Diagnostics.ProcessStartInfo
  $psi.FileName = $exe
  $psi.Arguments = ($args2 + @('--threads', "$THREADS") -join ' ')
  $psi.RedirectStandardOutput = $true
  $psi.UseShellExecute = $false
  if ($shear -ne '') { $psi.EnvironmentVariables['PD_SHEAR'] = $shear }
  if ($bend -ne '') { $psi.EnvironmentVariables['PD_BEND'] = $bend }
  $p = [System.Diagnostics.Process]::Start($psi)
  $t = $p.StandardOutput.ReadToEnd()
  $p.WaitForExit()
  return $t
}

# 只从输出里取"单子步 + 各阶段累计"做控制台摘要（看板侧走 make_report.js 的 parseBench）
function Metrics([string]$txt) {
  $head = [regex]::Match($txt, '子步数\s+(\d+)\s+总耗时\s+([\d.]+)\s+s\s+单子步\s+([\d.]+)\s+ms\s+平均迭代\s+([\d.]+)')
  $st = @{}
  foreach ($line in ($txt -split "`n")) {
    $m = [regex]::Match($line, '^\s*(\S+)\s+([\d.]+)\s+[\d.]+%')
    if ($m.Success) { $st[$m.Groups[1].Value] = [double]$m.Groups[2].Value }
  }
  return [pscustomobject]@{
    steps = [int]$head.Groups[1].Value
    perStep = [double]$head.Groups[3].Value
    iter = [double]$head.Groups[4].Value
    scatter = [double]$st['投影+散射']
    solve = [double]$st['全局回代']
    unpack = [double]$st['解包+判据']
  }
}

$summary = @()
foreach ($w in $workloads) {
  foreach ($c in $configs) {
    $best = $null
    for ($r = 1; $r -le $Reps; $r++) {
      $txt = Run-Cfg $w.args $c.shear $c.bend
      [IO.File]::WriteAllText((Join-Path $outDir "$($w.id)__$($c.id)_r$r.txt"), $txt, $UTF8)
      $m = Metrics $txt
      if ($null -eq $best -or $m.perStep -lt $best.perStep) { $best = $m }
    }
    $summary += [pscustomobject]@{
      workload = $w.id; cfg = $c.id; label = $c.label
      perStep = $best.perStep; iter = $best.iter
      scatter = $best.scatter / $best.steps; solve = $best.solve / $best.steps
      unpack = $best.unpack / $best.steps
    }
    "  $($w.id) / $($c.label)：$([Math]::Round($best.perStep,3)) ms/子步（迭代 $([Math]::Round($best.iter,2))）"
  }
}

""
"=== 汇总（ms/子步，已按子步数折算；4 线程） ==="
foreach ($w in $workloads) {
  $rows = $summary | Where-Object { $_.workload -eq $w.id }
  $base = ($rows | Where-Object { $_.cfg -eq 'none' }).perStep
  "--- $($w.id) ---"
  foreach ($r in $rows) {
    $d = 100.0 * ($r.perStep - $base) / $base
    # ⚠️ `-f` 的**对齐**必须是整数，不能把带符号的格式串塞进对齐位（`{2,+5:F1}` 会报
    #    "Input string was not in a correct format"）。带符号用格式串自身：
    #    `{2,7:+#,0.0;-#,0.0}`（对齐 7、格式含正/负两个 section）。
    "  {0,-10} 单子步 {1,7:F3}  ({2,7:+#,0.0;-#,0.0} %)   投影+散射 {3,7:F3}   全局回代 {4,7:F3}   解包+判据 {5,7:F3}" -f `
      $r.label, $r.perStep, $d, $r.scatter, $r.solve, $r.unpack
  }
}
"日志目录：$outDir"
