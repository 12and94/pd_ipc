# 全局步回代改造的 A/B/C 对照：基线（Eigen solve） / 变体 A（另存 CSR 按行 gather） / 变体 B（scatter 形式）
#
# 纪律（见 docs/perf.md §0）：同一会话内**交替**跑、取最小值；每个二进制都核对 OpenMP 运行时。
# 注意本文件必须带 UTF-8 BOM（PowerShell 5.1 读无 BOM 的中文 .ps1 会按 GBK 解码而语法错）；
# 变量名大小写不敏感，所以路径变量别叫 $a/$b（会被累加器覆盖）。

param([int]$Reps = 3)

$ErrorActionPreference = 'Stop'
$root = Split-Path -Parent $PSScriptRoot   # 本脚本在 <仓库根>/_perf/ 下（2026-09-24 从 _perf/ 迁入版本控制）
$exeEigen = Join-Path $root 'build\_baseline\pd_bench_phase3.exe'   # 基线：Eigen 的 solve
$exeCsr   = Join-Path $root 'build\_baseline\pd_bench_a_csr.exe'    # 变体 A：CSR 前代
$exeFinal = Join-Path $root 'build\Release\pd_bench.exe'            # 变体 B：scatter 前代（当前）
foreach ($path in @($exeEigen, $exeCsr, $exeFinal)) {
  if (-not (Test-Path $path)) { throw "缺二进制：$path" }
}

foreach ($path in @($exeEigen, $exeCsr, $exeFinal)) {
  $s = [System.Text.Encoding]::GetEncoding('iso-8859-1').GetString([System.IO.File]::ReadAllBytes($path))
  $rt = @()
  if ($s -match 'vcomp') { $rt += 'vcomp' }
  if ($s -match 'libomp') { $rt += 'libomp' }
  Write-Output ("运行时 {0}: {1}" -f (Split-Path -Leaf $path), ($rt -join '+'))
}

$cases = @(
  @{ name = '40x40_it40_t4';   args = @('--grid','40','40','--steps','300','--iters','40','--residual-tol','0.3','--threads','4') },
  @{ name = '60x60_it2_t4';    args = @('--grid','60','60','--steps','1000','--iters','2','--residual-tol','0.3','--threads','4') },
  @{ name = '100x100_it10_t4'; args = @('--grid','100','100','--steps','200','--iters','10','--residual-tol','0.3','--threads','4') },
  @{ name = '200x200_it10_t4'; args = @('--grid','200','200','--steps','80','--iters','10','--residual-tol','0.3','--threads','4') }
)

function Invoke-Bench([string]$exe, [string[]]$extra) {
  $out = & $exe @extra 2>&1
  $text = ($out | Out-String)
  $ms = [double][regex]::Match($text, '单子步\s+([0-9.]+)\s*ms').Groups[1].Value
  $sv = [double][regex]::Match($text, '全局回代\s+([0-9.]+)').Groups[1].Value
  $it = [regex]::Match($text, '平均迭代\s+([0-9.]+)').Groups[1].Value
  return [pscustomobject]@{ ms = $ms; solve = $sv; it = $it }
}

$lines = New-Object System.Collections.Generic.List[string]
$lines.Add("全局步回代：基线(Eigen solve) / 变体A(CSR 前代) / 变体B(scatter 前代)　交替各 $Reps 次取最小")
$lines.Add('')
$fmt = '{0,-18} {1,9} {2,9} {3,9} {4,8} {5,8} {6,9} {7,9} {8,9} {9,8}'
$header = $fmt -f '用例','基线ms','A ms','B ms','A/B ms','B/B ms','基线回代','A 回代','B 回代','迭代'
$lines.Add($header)
Write-Output $header

foreach ($case in $cases) {
  $rowsE = @(); $rowsA = @(); $rowsB = @()
  for ($i = 0; $i -lt $Reps; $i++) {
    $rowsE += (Invoke-Bench -exe $exeEigen -extra $case.args)
    $rowsA += (Invoke-Bench -exe $exeCsr   -extra $case.args)
    $rowsB += (Invoke-Bench -exe $exeFinal -extra $case.args)
  }
  $msE = ($rowsE | Measure-Object -Property ms -Minimum).Minimum
  $msA = ($rowsA | Measure-Object -Property ms -Minimum).Minimum
  $msB = ($rowsB | Measure-Object -Property ms -Minimum).Minimum
  $svE = ($rowsE | Measure-Object -Property solve -Minimum).Minimum
  $svA = ($rowsA | Measure-Object -Property solve -Minimum).Minimum
  $svB = ($rowsB | Measure-Object -Property solve -Minimum).Minimum
  $it = ($rowsB | Where-Object { $_.ms -eq $msB } | Select-Object -First 1).it
  $row = $fmt -f $case.name, ('{0:F3}' -f $msE), $msA, $msB,
      ('{0:F1}%' -f (100 * ($msA - $msE) / $msE)), ('{0:F1}%' -f (100 * ($msB - $msE) / $msE)),
      ('{0:F1}' -f $svE), $svA, $svB, $it
  Write-Output $row
  $lines.Add($row)
}
$lines.Add('')
$lines.Add('A/B ms = (变体 - 基线)/基线；负数为更快。回代列是各自最小 ms 那次运行的阶段累计值。')

$log = Join-Path $root 'build\_baseline\solve_ab.log'
[System.IO.File]::WriteAllText($log, ($lines -join "`r`n"), (New-Object System.Text.UTF8Encoding($false)))
Write-Output "日志：$log"
