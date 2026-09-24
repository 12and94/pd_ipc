# _perf/ab_sharefactor.ps1 —— 「因子只存一份」上线前后的成对 A/B（同一轮会话内交替、取最小值）
#
# 对照：
#   基线 = build\_baseline\pd_bench_pre_sharefactor.exe（改动前：CSC 因子里三块各存一份，各扫自己的分量）
#   新版 = build\Release\pd_bench.exe               （改动后：再存一份紧凑单份，三条线程读同一份）
#
# 纪律（docs/perf.md §0，都是踩坑换来的）：
#   · 同一轮会话内**交替**跑、重复取**最小值**（单次运行波动可达 2×）；
#   · 只比较**物理完全相同**的运行 —— 本脚本会把两侧的"收敛统计 / 最终应变 / 全程峰值 / 能量"
#     抽出来逐字比对，不同就标出来（那说明两侧不是同一次物理，数字不可比）；
#   · 换过 OpenMP 开关后要核对二进制实际链接的运行时（这里每次跑都核对 vcomp/libomp）。
#
# ⚠️ 本文件必须带 UTF-8 BOM（PowerShell 5.1 读无 BOM 的中文 .ps1 会按 GBK 解码而语法错）。
# ⚠️ 日志写 build\_baseline\sharefactor_ab.log；**看板不读它**（同 solve_ab.log：一次性对照，
#    没有画进看板的板块；要画的话按 _perf/make_report.js 的 parseTwoWay/parsePhase3 加一节）。
#
# 用法： .\_perf\ab_sharefactor.ps1 -Reps 5
param([int]$Reps = 5)

$ErrorActionPreference = 'Stop'
$root = Split-Path -Parent $PSScriptRoot   # 本脚本在 <仓库根>/_perf/ 下（2026-09-24 起）
$OLD = Join-Path $root 'build\_baseline\pd_bench_pre_sharefactor.exe'
$NEW = Join-Path $root 'build\Release\pd_bench.exe'
$LOG = Join-Path $root 'build\_baseline\sharefactor_ab.log'
$UTF8 = New-Object Text.UTF8Encoding($false)

if (-not (Test-Path $OLD)) { throw "缺少基线二进制 $OLD（改动前归档；见 docs/perf.md §16）" }
if (-not (Test-Path $NEW)) { throw "缺少 build\Release\pd_bench.exe（先跑 .\tools\build.ps1）" }

# 运行时核对（换过 OpenMP 开关后必做）：exe 里不应出现 libomp
foreach ($p in @($OLD, $NEW)) {
  $ascii = [Text.Encoding]::GetEncoding('iso-8859-1').GetString([IO.File]::ReadAllBytes($p))
  $rt = @()
  if ($ascii -match 'vcomp') { $rt += 'vcomp' }
  if ($ascii -match 'libomp') { $rt += 'libomp' }
  "运行时 {0}: {1}" -f (Split-Path -Leaf $p), (($rt -join '+') -replace '^$', '未识别')
}

function Metrics([string]$txt) {
  $head = [regex]::Match($txt, '单子步\s+([\d.]+)\s+ms\s+平均迭代\s+([\d.]+)')
  $solve = 0.0
  foreach ($line in ($txt -split "`n")) {
    $m = [regex]::Match($line, '^\s*全局回代\s+([\d.]+)\s+[\d.]+%')
    if ($m.Success) { $solve = [double]$m.Groups[1].Value }
  }
  # 物理（逐字可比的两侧）：收敛统计 / 最终应变 / 全程峰值 / 能量
  $phys = (($txt -split "`r?`n") | Where-Object { $_ -match '收敛：|最终应变|全程最大应变|能量:' }) -join "`n"
  return [pscustomobject]@{
    total = if ($head.Success) { [double]$head.Groups[1].Value } else { $null }
    iter  = if ($head.Success) { [double]$head.Groups[2].Value } else { $null }
    solve = $solve
    phys  = $phys
  }
}

function RunBench([string]$exe, [string[]]$a) {
  return ((& $exe @a 2>&1) -join "`n")
}

$cases = @(
  @{ name = '40x40_it40_t4';   args = @('--grid', '40', '40', '--steps', '300', '--iters', '40', '--residual-tol', '0.3', '--threads', '4') },
  @{ name = '60x60_it2_t4';    args = @('--grid', '60', '60', '--steps', '300', '--iters', '2', '--residual-tol', '0.3', '--threads', '4') },
  @{ name = '100x100_it10_t4'; args = @('--grid', '100', '100', '--steps', '300', '--iters', '10', '--residual-tol', '0.3', '--threads', '4') },
  @{ name = '200x200_it10_t4'; args = @('--grid', '200', '200', '--steps', '100', '--iters', '10', '--residual-tol', '0.3', '--threads', '4') }
)

$rows = @()
$bad = 0
foreach ($c in $cases) {
  $bv = 1e30; $bs = 0.0; $bphys = ''
  $nv = 1e30; $ns = 0.0; $nphys = ''
  for ($r = 0; $r -lt $Reps; $r++) {
    $o = Metrics (RunBench $OLD $c.args)
    if ($null -ne $o.total -and $o.total -lt $bv) { $bv = $o.total; $bs = $o.solve; $bphys = $o.phys }
    $n = Metrics (RunBench $NEW $c.args)
    if ($null -ne $n.total -and $n.total -lt $nv) { $nv = $n.total; $ns = $n.solve; $nphys = $n.phys }
  }
  $delta = ($nv / $bv - 1) * 100
  $dSolve = if ($bs -gt 0) { ($ns / $bs - 1) * 100 } else { 0 }
  $same = ($bphys -ceq $nphys)
  if (-not $same) { $bad++ }
  $rows += ('{0,-18}{1,9:F3}{2,9:F3}{3,9:F1}%{4,11:F1}{5,11:F1}{6,10:F1}%{7,8}' -f `
             $c.name, $bv, $nv, $delta, $bs, $ns, $dSolve, $(if ($same) { 'same' } else { 'DIFF' }))
  "  {0}：{1:F3} → {2:F3} ms/子步  {3:F1}%（回代 {4:F1} → {5:F1} ms，{6:F1}%）物理{7}" -f `
    $c.name, $bv, $nv, $delta, $bs, $ns, $dSolve, $(if ($same) { '相同' } else { '**不同，数字不可比**' })
}

$out = @()
$out += '「因子只存一份」上线前后 A/B：紧凑单份（三条线程读同一份） vs 三块各存一份'
$out += "基线 pd_bench_pre_sharefactor.exe / 新版 pd_bench.exe，每例交替各 $Reps 次，取最小值（ms/子步）"
$out += ''
$out += '用例                   基线       新版      变化   基线回代    新版回代   回代变化  物理'
$out += $rows
if ($bad -gt 0) { $out += ""; $out += "⚠️ 有 $bad 个用例两侧物理不同 ⇒ 那几行不可比（见 docs/perf.md §0 纪律 2）" }
[IO.File]::WriteAllText($LOG, ($out -join "`r`n") + "`r`n", $UTF8)
"日志：$LOG"
if ($bad -gt 0) { exit 2 }
