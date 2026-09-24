# Phase 3 A/B 测量：并行残差（gather，与解包融合）vs 基线（串行残差）
#
# 纪律（见 docs/perf.md §0）：
#   · 只在**同一会话内交替**跑 A/B（A B A B A B），取各自最小值 —— 跨会话的读数不可比；
#   · 两个二进制的 OpenMP 运行时必须相同（本脚本先自查 vcomp/libomp），否则读数无意义；
#   · 物理必须同一口径：打印"平均迭代"与残差读数，用来确认两次运行的物理一致。
#
# 基线 = build/_baseline/pd_bench_phase2.exe（着色散射 + single 里串行残差）
# 新版 = build/Release/pd_bench.exe            （着色散射 + 逐顶点 gather 并行残差）
#
# 三个 PowerShell 5.1 的坑（都踩过，改这个文件时别改回去）：
#   1) 变量名**大小写不敏感**：累加器若叫 $a/$b，就会把上面的 $A/$B（exe 路径）覆盖成
#      空数组，调用时 $exe 变成空串，报 "The expression after '&' ... (:String)"。
#      所以这里叫 $exeBase/$exeNew 与 $rowsBase/$rowsNew。
#   2) 5.1 读**没有 BOM** 的 .ps1 会按 ANSI（本机 GBK）解码 ⇒ 中文注释/字串全部乱码、
#      还会把引号吃掉导致语法错。本文件因此**必须带 UTF-8 BOM**（build/ 已被 gitignore，
#      不影响仓库"UTF-8 无 BOM"的约定）。注意：用某些编辑器/工具改完会掉 BOM，改完请复查：
#        $b=[IO.File]::ReadAllBytes('_perf\ab_residual.ps1'); '{0:X2}{1:X2}{2:X2}' -f $b[0],$b[1],$b[2]   # 应为 EFBBBF
#   3) `-f` 的实参**不能**用反引号换行接在**方法调用的括号里**（实测：`$list.Add('{0}{1}' -f `⏎ 'a','b')`
#      只把第一个实参交给 -f，报 FormatError "Index ... less than the size of the argument list"）。
#      所以一律先 `$s = $fmt -f ...`（语句级 + 格式串放变量里），再 Add($s)。

param([int]$Reps = 3)

$ErrorActionPreference = 'Stop'
$root = Split-Path -Parent $PSScriptRoot   # 本脚本在 <仓库根>/_perf/ 下（2026-09-24 从 _perf/ 迁入版本控制）   # 仓库根
$exeBase = Join-Path $root 'build\_baseline\pd_bench_phase2.exe'
$exeNew = Join-Path $root 'build\Release\pd_bench.exe'
foreach ($path in @($exeBase, $exeNew)) {
  if (-not (Test-Path $path)) { throw "缺二进制：$path" }
}

# 运行时自查：换过 /openmp 开关之后必须核对（docs/perf.md §6.2）
foreach ($path in @($exeBase, $exeNew)) {
  $s = [System.Text.Encoding]::GetEncoding('iso-8859-1').GetString([System.IO.File]::ReadAllBytes($path))
  $rt = @()
  if ($s -match 'vcomp') { $rt += 'vcomp' }
  if ($s -match 'libomp') { $rt += 'libomp' }
  Write-Output ("运行时 {0}: {1}" -f (Split-Path -Leaf $path), ($rt -join '+'))
}

$cases = @(
  @{ name = '40x40_it40_t1';    args = @('--grid','40','40','--steps','300','--iters','40','--residual-tol','0.3','--threads','1') },
  @{ name = '40x40_it40_t4';    args = @('--grid','40','40','--steps','300','--iters','40','--residual-tol','0.3','--threads','4') },
  @{ name = '40x40_it40_t18';   args = @('--grid','40','40','--steps','300','--iters','40','--residual-tol','0.3','--threads','18') },
  @{ name = '60x60_it2_t4';     args = @('--grid','60','60','--steps','2000','--iters','2','--residual-tol','0.3','--threads','4') },
  @{ name = '100x100_it10_t4';  args = @('--grid','100','100','--steps','300','--iters','10','--residual-tol','0.3','--threads','4') },
  @{ name = '200x200_it10_t4';  args = @('--grid','200','200','--steps','100','--iters','10','--residual-tol','0.3','--threads','4') },
  @{ name = '200x200_it10_t18'; args = @('--grid','200','200','--steps','100','--iters','10','--residual-tol','0.3','--threads','18') }
)

function Invoke-Bench([string]$exe, [string[]]$extra) {
  # 具名绑定 + `@extra` 数组展开；直接 & 调用（本机 PATH 上没有别的 shell 假设）。
  $out = & $exe @extra 2>&1
  $text = ($out | Out-String)
  $ms = [double][regex]::Match($text, '单子步\s+([0-9.]+)\s*ms').Groups[1].Value
  $it = [regex]::Match($text, '平均迭代\s+([0-9.]+)').Groups[1].Value
  $un = [double][regex]::Match($text, '解包\+判据\s+([0-9.]+)').Groups[1].Value
  $sc = [double][regex]::Match($text, '散射\s+([0-9.]+)').Groups[1].Value
  $so = [double][regex]::Match($text, '全局回代\s+([0-9.]+)').Groups[1].Value
  return [pscustomobject]@{ ms = $ms; it = $it; unpack = $un; scatter = $sc; solve = $so }
}

$lines = New-Object System.Collections.Generic.List[string]
$lines.Add('Phase 3 A/B：并行残差（逐顶点 gather + 与解包融合） vs 串行残差')
$lines.Add("基线 pd_bench_phase2.exe / 新版 pd_bench.exe，每例交替各 $Reps 次，取最小值（ms/子步）")
$lines.Add('')
$headerFmt = '{0,-18} {1,9} {2,9} {3,8} {4,7} {5,7} {6,9} {7,9} {8,9} {9,9} {10,8} {11,8}'
$header = $headerFmt -f `
  '用例','基线','新版','变化','基线迭代','新版迭代','基线散射','新版散射','基线解包','新版解包','基线回代','新版回代'
$lines.Add($header)

$rowFmt = '{0,-18} {1,9:F3} {2,9:F3} {3,7:F1}% {4,7} {5,7} {6,9:F1} {7,9:F1} {8,9:F1} {9,9:F1} {10,8:F1} {11,8:F1}'

foreach ($case in $cases) {
  $rowsBase = @(); $rowsNew = @()
  for ($i = 0; $i -lt $Reps; $i++) {
    $rowsBase += (Invoke-Bench -exe $exeBase -extra $case.args)
    $rowsNew += (Invoke-Bench -exe $exeNew -extra $case.args)
  }
  $msBase = ($rowsBase | Measure-Object -Property ms -Minimum).Minimum
  $msNew = ($rowsNew | Measure-Object -Property ms -Minimum).Minimum
  $unBase = ($rowsBase | Measure-Object -Property unpack -Minimum).Minimum
  $unNew = ($rowsNew | Measure-Object -Property unpack -Minimum).Minimum
  $soBase = ($rowsBase | Measure-Object -Property solve -Minimum).Minimum
  $soNew = ($rowsNew | Measure-Object -Property solve -Minimum).Minimum
  $scBase = ($rowsBase | Measure-Object -Property scatter -Minimum).Minimum
  $scNew = ($rowsNew | Measure-Object -Property scatter -Minimum).Minimum
  $itBase = ($rowsBase | Where-Object { $_.ms -eq $msBase } | Select-Object -First 1).it
  $itNew = ($rowsNew | Where-Object { $_.ms -eq $msNew } | Select-Object -First 1).it
  $delta = 100.0 * ($msNew - $msBase) / $msBase
  $row = $rowFmt -f `
    $case.name, $msBase, $msNew, $delta, $itBase, $itNew, $scBase, $scNew, $unBase, $unNew, $soBase, $soNew
  Write-Output $row
  $lines.Add($row)
}
$lines.Add('')
$lines.Add('变化 = (新版 - 基线)/基线：负数表示新版更快。解包/回代是各自最小 ms 那次运行的阶段累计值。')
$lines.Add('迭代数应逐例相同 —— 残差两个口径实测逐位相同，因此放行判据不会分叉。')

$log = Join-Path $root 'build\_baseline\phase3_ab.log'
[System.IO.File]::WriteAllText($log, ($lines -join "`r`n"), (New-Object System.Text.UTF8Encoding($false)))
Write-Output "日志：$log"
