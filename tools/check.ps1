# tools/check.ps1 —— 一条命令跑完全部验收（M4 回归固化，2026-09-24 建）
#
# 为什么需要它：验收集本来是"照着 HANDOFF §2 的清单逐条手跑"，于是"今天到底跑全了没有"
# 只能靠人记；而文档里的数字（403 断言 / 22 条 / 8 项 / 233 …）都是**手抄**的，漂了没人知道。
# 同一类债在 2026-09-23 与 09-24 连还了两轮（断言计数三种说法、README 查看器数字自相矛盾、
# 脚本 BOM 丢失、看板口径与文档冲突…）。本脚本把六件事一次做完：
#   ① 跑完验收集：任一程序非 0 退出即失败；
#   ② 数字权威化：把每个程序**自报的计数**抽出来，再与文档里"当前计数"的口径句逐处比对；
#   ③ 编码门禁：所有跟踪的 *.ps1 必须带 UTF-8 BOM；所有跟踪文本必须 UTF-8 无 BOM + LF；
#   ④ 结构性回归：默认约束集的物理输出必须与存档基线 pd_bench_pre_fuse.exe **逐字相同**；
#      看板 HTML 存在时再跑一遍 _perf/check_report.js（DOM 桩真跑渲染 + 关键数字反查）；
#   ⑤ 迭代预算三档预设（docs/perf.md §18）的自报值 vs 文档表；
#   ⑥ 文档命令一致性（2026-09-24 新增）：文档里出现的 `xxx.exe` 必须在 CMake 构建表内，
#      且每个 `--flag` 必须在**该程序**接受的表内 —— 抓"程序没构建""参数名写错"这一类漂移；
#      另外反查两条"结构性前提"（pd_solvecomp / pd_sharefactor）：⊗ 结构不被破坏、
#      三分量因子切片同构且逐位相同、单份求解 == 各份求解。
#
# 用法： .\tools\check.ps1                 # 全部（约 10 秒，不含 A/B 测量）
#        .\tools\check.ps1 -NoPhysics      # 跳过物理基线比对（存档二进制不在时）
#        .\tools\check.ps1 -NoBuildCheck   # 跳过编码门禁（只想快速跑程序时）
# 退出码：0 = 全部通过；1 = 有失败项（逐条打印）。
#
# 纪律：本脚本只**读**仓库与二进制，不改任何文件；发现数字对不上时是"报错"，不是"帮你改文档"。
# ⚠️ 本文件必须带 UTF-8 BOM（PowerShell 5.1 读无 BOM 的中文 .ps1 会按 GBK 解码而语法错）。
param([switch]$NoPhysics, [switch]$NoBuildCheck)

$ErrorActionPreference = 'Stop'
$root = Split-Path -Parent $PSScriptRoot
Set-Location $root

$script:Fail = @()
$script:Warn = @()

function Step([string]$t) { Write-Host "`n=== $t ===" -ForegroundColor Cyan }
function Good([string]$m) { Write-Host "  [ok]   $m" -ForegroundColor Green }
function Bad([string]$m) { $script:Fail += $m; Write-Host "  [FAIL] $m" -ForegroundColor Red }
function Note([string]$m) { Write-Host "         $m" -ForegroundColor DarkGray }

function RunExe([string]$name, [string[]]$exeArgs) {
  $exe = Join-Path $root "build\Release\$name.exe"
  if (-not (Test-Path $exe)) {
    Bad "缺可执行文件 build\Release\$name.exe（先跑 .\tools\build.ps1）"
    return $null
  }
  # ⚠️ 必须在本地把 ErrorActionPreference 降下来：PS 5.1 会把**原生命令写到 stderr 的任何一行**
  # 包装成 NativeCommandError，而本脚本开头是 'Stop' ⇒ 一行无害的告警就会让整个门禁中断。
  # （第一次就是这么被 pd_bendaudit 拦住 的：求解器对合成矩阵打了一行"退回旧路径"的说明。）
  $eap = $ErrorActionPreference
  $ErrorActionPreference = 'Continue'
  try {
    $out = (& $exe @exeArgs 2>&1 | Out-String)
    $code = $LASTEXITCODE
  } finally {
    $ErrorActionPreference = $eap
  }
  return [pscustomobject]@{ Name = $name; Code = $code; Out = $out }
}

# 从程序输出里抓一个数：抓不到就是**自报格式变了**（那也是失败 —— 否则文档数字会悄悄失去来源）
function Grab([string]$text, [string]$pattern, [string]$what) {
  $m = [regex]::Match($text, $pattern)
  if (-not $m.Success) { Bad "没抓到自报数字：$what（模式：$pattern）"; return $null }
  return [int]$m.Groups[1].Value
}

# ---------------------------------------------------------------- ① 验收集
Step '① 验收集：全部程序必须 exit 0'

$runs = [ordered]@{}
$runs['pd_check']                   = RunExe 'pd_check'                   @()
$runs['test_primitives']            = RunExe 'test_primitives'            @()
$runs['test_spring_vertical']       = RunExe 'test_spring_vertical'       @()
$runs['test_convergence_criterion'] = RunExe 'test_convergence_criterion' @()
$runs['pd_trans']                   = RunExe 'pd_trans'                   @()
$runs['pd_chain']                   = RunExe 'pd_chain'                   @()
$runs['pd_diraudit']                = RunExe 'pd_diraudit'                @()
$runs['pd_bendaudit']               = RunExe 'pd_bendaudit'               @()
$runs['pd_solvecomp']               = RunExe 'pd_solvecomp'               @('--grid', '40', '40', '--reps', '50', '--threads', '4')

foreach ($k in $runs.Keys) {
  if ($null -eq $runs[$k]) { continue }
  if ($runs[$k].Code -eq 0) { Good "$k  exit=0" } else { Bad "$k  exit=$($runs[$k].Code)" }
}

# ---------------------------------------------------------------- ② 数字权威化
Step '② 数字权威化：程序自报的计数 == 文档里"当前计数"的口径句'

$counts = [ordered]@{}
if ($runs['test_primitives']) {
  $counts.primitives = Grab $runs['test_primitives'].Out '(\d+) 个断言, \d+ 个失败测试' 'test_primitives 断言数'
}
if ($runs['test_spring_vertical']) {
  $counts.spring = Grab $runs['test_spring_vertical'].Out '(\d+) 个断言, \d+ 个失败测试' 'test_spring_vertical 断言数'
}
if ($runs['test_convergence_criterion']) {
  $counts.conv = Grab $runs['test_convergence_criterion'].Out '(\d+) 个断言, \d+ 个失败测试' 'test_convergence_criterion 断言数'
}
if ($counts.primitives -and $counts.spring -and $counts.conv) {
  $counts.total = $counts.primitives + $counts.spring + $counts.conv
}
if ($runs['pd_check']) {
  $counts.check = Grab $runs['pd_check'].Out '全部正确（(\d+) 项，失败 \d+ 项）' 'pd_check 项数'
}
if ($runs['pd_bendaudit']) {
  $counts.bendaudit = Grab $runs['pd_bendaudit'].Out '断言 (\d+) 条，失败 \d+ 项' 'pd_bendaudit 断言总数'
}
if ($counts.Count) {
  ($counts.Keys | ForEach-Object { "$_=$($counts[$_])" }) -join '  ' | ForEach-Object { Note "自报：$_" }
}

# 文档里"当前状态"的口径句。**只列现状句** —— 历史段落（例如"2026-09-21 复跑三套 295 断言"）
# 故意不列：它们本来就该保留当时的数字。若某句确实被重写，请同步更新这张表（不要删条目）。
$docRules = @()
if ($counts.total)      { $docRules += @{ file='README.md';          pat='三个测试套件 \*\*(\d+)\*\* 个断言';             want=$counts.total;      desc='三套测试合计' } }
if ($counts.primitives) { $docRules += @{ file='README.md';          pat='`test_primitives`（27 个测试 / (\d+) 个断言';    want=$counts.primitives; desc='test_primitives' } }
if ($counts.spring)     { $docRules += @{ file='README.md';          pat='`test_spring_vertical`（11 个测试 / (\d+) 个断言'; want=$counts.spring;   desc='test_spring_vertical' } }
if ($counts.check)      { $docRules += @{ file='README.md';          pat='`pd_check` (\d+) 项验收全通过';                  want=$counts.check;      desc='pd_check' } }
if ($counts.total)      { $docRules += @{ file='docs\code-map.md';   pat='# 合计 (\d+) 断言';                              want=$counts.total;      desc='三套测试合计' } }
if ($counts.primitives) { $docRules += @{ file='docs\code-map.md';   pat='test_primitives\.exe[ \t]+#[ \t]+(\d+) 断言';    want=$counts.primitives; desc='test_primitives' } }
if ($counts.total)      { $docRules += @{ file='docs\plan.md';       pat='三套测试 \*\*(\d+) 断言\*\*全绿';                 want=$counts.total;      desc='三套测试合计' } }
if ($counts.total)      { $docRules += @{ file='_perf\HANDOFF.md';   pat='三套测试 \*\*(\d+) 断言\*\*';                     want=$counts.total;      desc='三套测试合计' } }
if ($counts.primitives) { $docRules += @{ file='_perf\HANDOFF.md';   pat='test_primitives\.exe[ \t]+#[ \t]+(\d+) 断言';     want=$counts.primitives; desc='test_primitives' } }
if ($counts.spring)     { $docRules += @{ file='_perf\HANDOFF.md';   pat='test_spring_vertical\.exe[ \t]+#[ \t]+(\d+) 断言'; want=$counts.spring;    desc='test_spring_vertical' } }
if ($counts.conv)       { $docRules += @{ file='_perf\HANDOFF.md';   pat='test_convergence_criterion\.exe # (\d+) 断言';   want=$counts.conv;       desc='test_convergence_criterion' } }
if ($counts.bendaudit)  { $docRules += @{ file='_perf\HANDOFF.md';   pat='弯曲约束的代数审计（断言 (\d+) 条';                want=$counts.bendaudit;  desc='pd_bendaudit' } }

foreach ($r in $docRules) {
  $path = Join-Path $root $r.file
  if (-not (Test-Path $path)) { Bad "$($r.file) 不存在"; continue }
  $text = [IO.File]::ReadAllText($path, [Text.Encoding]::UTF8)
  $m = [regex]::Match($text, $r.pat)
  if (-not $m.Success) {
    Bad "$($r.file)：$($r.desc) 的口径句不见了（若确实重写了措辞，请同步更新 tools/check.ps1 的规则表）"
  } elseif ([int]$m.Groups[1].Value -ne $r.want) {
    Bad "$($r.file)：$($r.desc) 写的是 $($m.Groups[1].Value)，程序自报 $($r.want)"
  } else {
    Good "$($r.file)：$($r.desc) = $($m.Groups[1].Value)（与程序自报一致）"
  }
}

# "233 + 155 + 15 = 403" 这种四数并列式：两处写法不同（code-map 给首尾两数加粗、plan 整句加粗），
# 两种形式都认；两个文件各必须命中一处，四个数都要与程序自报一致。
if ($counts.total) {
  $forms = @(
    '\*\*(\d+)\*\* \+ (\d+) \+ (\d+) = \*\*(\d+) 断言\*\*',
    '\*\*(\d+) \+ (\d+) \+ (\d+) = (\d+) 断言\*\*'
  )
  foreach ($f in 'docs\code-map.md', 'docs\plan.md') {
    $text = [IO.File]::ReadAllText((Join-Path $root $f), [Text.Encoding]::UTF8)
    $m = $null
    foreach ($p in $forms) { $m = [regex]::Match($text, $p); if ($m.Success) { break } }
    if (-not $m -or -not $m.Success) {
      Bad "$f：'233 + 155 + 15 = 403 断言' 那条口径句不见了"
    } else {
      $got = @([int]$m.Groups[1].Value, [int]$m.Groups[2].Value, [int]$m.Groups[3].Value, [int]$m.Groups[4].Value)
      $want = @($counts.primitives, $counts.spring, $counts.conv, $counts.total)
      if (($got -join ',') -ne ($want -join ',')) { Bad "$f：写成 $($got -join ' + ')，程序自报 $($want -join ' + ')" }
      else { Good "$f：$($got -join ' + ')，与程序自报一致" }
    }
  }
}

# pd_solvecomp 没有 item() 计数，改用**结构性判据**反查（比数列断言更强）
if ($runs['pd_solvecomp']) {
  $o = $runs['pd_solvecomp'].Out
  if ($o -match '跨分量的 i≠j 项（填充跨分量）：0') { Good 'pd_solvecomp：跨分量填充 == 0' }
  else { Bad 'pd_solvecomp：跨分量填充不为 0（⊗ 结构被破坏 ⇒ 三分量并行会静默退化）' }
  if ($o -match '最大逐元素差 0\.000e\+00') { Good 'pd_solvecomp：三分量拆 vs 整趟 逐位相同' }
  else { Bad 'pd_solvecomp：三分量拆与整趟不再逐位相同' }
}

# 因子"三份完全一致"这条前提：只有 pd_sharefactor 在查（它是"只存一份"与任何"共用因子"优化的前提）。
# 这里只跑 40×40 / reps 1 / 不测冷档 —— 只取断言，不取时间（时间在 docs/perf.md §16 里另记）。
$sf = RunExe 'pd_sharefactor' @('--grid', '40', '40', '--reps', '1', '--flush-mb', '0')
if ($sf) {
  if ($sf.Code -ne 0) { Bad "pd_sharefactor exit=$($sf.Code)（前提自检或逐位比对不通过）" }
  elseif ($sf.Out -match '不一致的列 = 0；对角 1/D 不一致 = 0') { Good 'pd_sharefactor：三分量因子切片同构且数值逐位相同' }
  else { Bad 'pd_sharefactor：三分量切片的"同构 + 逐位相同"前提不成立' }
  if ($sf.Out -match '最大逐元素差 0\.000e\+00') { Good 'pd_sharefactor：单份求解 == 各份求解（逐位）' }
  else { Bad 'pd_sharefactor：单份求解与各份求解不再逐位相同' }
}
# 残差-迭代曲线工具（_verify/residual_trace.cpp，pd_restrace）：它给出的 2 范数曲线依赖
# 「复算每顶点残差」（生产的 nonlinearResidual 只返回 max 范数），而复算必须与生产**逐位一致**，
# 否则曲线不可信。这里用最小工况把这条自校验钉住。
$rt = RunExe 'pd_restrace' @('--grid', '20', '--steps', '5', '--iters', '4')
if ($rt) {
  if ($rt.Code -ne 0) { Bad "pd_restrace exit=$($rt.Code)（复算残差与生产不一致 ⇒ 曲线不可信）" }
  elseif ($rt.Out -match '一致 ✓') { Good 'pd_restrace：复算残差与生产 nonlinearResidual 逐位一致（2 范数曲线可信）' }
  else { Bad 'pd_restrace：没看到「一致 ✓」自校验行' }
}

# ---------------------------------------------------------------------------
# ⑥ 迭代预算三档预设（docs/perf.md §18）：唯一权威是 core/sim/Scene.h 的 iterationPresetTable()。
# 做法与其他"数字权威化"一致：让**程序自己报**每个档位的生效值，再要求文档那张表写着同样的数。
# 这样"档位被改"和"文档没跟着改"都会红。顺手验两条容易悄悄坏的性质：
#   ① --preset 真的生效（而不是被忽略）；② 不写 --preset 时**维持库默认**（10 次 / 1e-3）——
#   后者是全部历史基线（含第 ④ 步的物理逐字比对）的前提，谁改了默认值这里就会响。
$presetWant = @{ 'preview' = @(20, '0.3'); 'realtime' = @(40, '0.3'); 'accurate' = @(80, '0.001') }
$presetGot = @{}
foreach ($name in @('preview', 'realtime', 'accurate')) {
  $r = RunExe 'pd_bench' @('--preset', $name, '--grid', '20', '20', '--steps', '1')
  if (-not $r) { continue }
  if ($r.Code -ne 0) { Bad "pd_bench --preset $name exit=$($r.Code)"; continue }
  $mi = [regex]::Match($r.Out, 'PD 迭代上限 (\d+)')
  $mt = [regex]::Match($r.Out, '残差容差 ([\d.eE+-]+)')
  if (-not $mi.Success -or -not $mt.Success) { Bad "没抓到 --preset $name 的自报值（表头格式变了？）"; continue }
  $gotIters = [int]$mi.Groups[1].Value
  $gotTol = [double]$mt.Groups[1].Value
  $want = $presetWant[$name]
  if ($gotIters -ne $want[0] -or [math]::Abs($gotTol - [double]$want[1]) -gt 1e-12) {
    Bad "档位 $name 自报 $gotIters 次 / $gotTol，与门禁基准 $($want[0]) 次 / $($want[1]) 不一致"
  } else {
    $presetGot[$name] = @($gotIters, $want[1])
  }
}
if ($presetGot.Count -eq 3) {
  Good '三档自报与基准一致（preview 20/0.3、realtime 40/0.3、accurate 80/1e-3）'
} else {
  Bad "三档没有全部自报成功（只拿到 $($presetGot.Count) 个）"
}
$rd = RunExe 'pd_bench' @('--grid', '20', '20', '--steps', '1')
if ($rd -and $rd.Code -eq 0) {
  if ($rd.Out -match 'PD 迭代上限 10' -and $rd.Out -match '残差容差 0\.001') {
    Good '不写 --preset 时维持库默认（10 次 / 1e-3）⇒ 历史基线口径不变'
  } else {
    Bad '不写 --preset 的默认值变了 —— 会破坏全部历史基线的时间口径'
  }
}
# 文档那张表必须写着程序自报的数（容忍 1e-3 的几种写法）
$perfPath = Join-Path $root 'docs\perf.md'
if (Test-Path $perfPath) {
  $perfLines = [IO.File]::ReadAllLines($perfPath, [Text.Encoding]::UTF8)
  foreach ($name in @('preview', 'realtime', 'accurate')) {
    if (-not $presetGot.ContainsKey($name)) { continue }
    $iters = [string]$presetGot[$name][0]
    $tolStr = $presetGot[$name][1]
    $tolAlt = if ($tolStr -eq '0.001') { @('1e-3', '1E-3', '0.001') } else { @($tolStr) }
    $found = $false
    foreach ($line in $perfLines) {
      if ($line -notmatch ('`' + $name + '`')) { continue }
      if ($line -notmatch ('\|\s*' + $iters + '\s*\|')) { continue }
      foreach ($alt in $tolAlt) { if ($line.Contains($alt)) { $found = $true } }
    }
    if ($found) { Good "perf.md §18 档位表：$name = $iters 次 / $tolStr" }
    else { Bad "perf.md §18 档位表里 $name 那一行与程序自报（$iters 次 / $tolStr）不一致" }
  }
}

# ---------------------------------------------------------------- ③ 编码门禁
if (-not $NoBuildCheck) {
  Step '③ 编码门禁：*.ps1 必须带 BOM；跟踪文本必须 UTF-8 无 BOM + LF'
  $tracked = & git -C $root ls-files
  $textExt = @('.md', '.cpp', '.h', '.js', '.txt', '.cmake', '.json', '.ps1', '.gitignore', '.gitattributes')
  $checked = 0; $badBom = 0; $badEol = 0
  foreach ($f in $tracked) {
    $ext = [IO.Path]::GetExtension($f)
    $isText = ($textExt -contains $ext) -or ($f -eq 'CMakeLists.txt') -or ($f -eq '.gitignore') -or ($f -eq '.gitattributes')
    if (-not $isText) { continue }
    $bytes = [IO.File]::ReadAllBytes((Join-Path $root $f))
    $hasBom = ($bytes.Length -ge 3 -and $bytes[0] -eq 0xEF -and $bytes[1] -eq 0xBB -and $bytes[2] -eq 0xBF)
    $checked++
    if ($ext -eq '.ps1') {
      if (-not $hasBom) { Bad "编码：$f 是 .ps1 但没有 UTF-8 BOM（PS 5.1 会按 GBK 解码 ⇒ 语法错）"; $badBom++ }
    } else {
      if ($hasBom) { Bad "编码：$f 有 BOM（仓库文本要求 UTF-8 无 BOM）"; $badBom++ }
      $crlf = 0
      for ($i = 0; $i -lt $bytes.Length - 1; $i++) { if ($bytes[$i] -eq 0x0D -and $bytes[$i + 1] -eq 0x0A) { $crlf++; break } }
      if ($crlf -gt 0) { Bad "编码：$f 含 CRLF 行尾（仓库要求 LF）"; $badEol++ }
    }
  }
  if ($badBom -eq 0 -and $badEol -eq 0) { Good "编码：$checked 个跟踪文本文件全部合规（*.ps1 带 BOM；其余无 BOM + LF）" }
}

# ---------------------------------------------------------------- ④ 物理基线逐字比对
if (-not $NoPhysics) {
  Step '④ 默认约束集物理输出 vs 存档基线（逐字相同）'
  $old = Join-Path $root 'build\_baseline\pd_bench_pre_fuse.exe'
  $new = Join-Path $root 'build\Release\pd_bench.exe'
  if (-not (Test-Path $old)) {
    $script:Warn += "没有存档基线 $old ⇒ 跳过逐字比对（它不在版本控制里，见 _perf/HANDOFF.md）"
    Note '跳过：存档二进制不存在'
  } elseif (-not (Test-Path $new)) {
    Bad '缺 build\Release\pd_bench.exe'
  } else {
    $a = @('--grid', '40', '40', '--steps', '300', '--iters', '40', '--residual-tol', '0.3', '--threads', '4')
    $pick = {
      param($s)
      (($s -split "`r?`n") | Where-Object { $_ -match '收敛：|最终应变|全程最大应变|能量:' }) -join "`n"
    }
    $o1 = & $pick ((& $old @a 2>&1) -join "`n")
    $o2 = & $pick ((& $new @a 2>&1) -join "`n")
    if ($o1 -ceq $o2) {
      Good '默认约束集物理输出逐字相同（末次残差 / 应变 / 峰值 / 能量）'
      ($o2 -split "`n") | ForEach-Object { Note $_.Trim() }
    } else {
      Bad '默认约束集的物理输出与存档基线**不一致** ⇒ 要么物理真的变了（要重新标定基线），要么默认值被改了'
      Compare-Object ($o1 -split "`n") ($o2 -split "`n") | ForEach-Object { Note ("{0} {1}" -f $_.SideIndicator, $_.InputObject) }
    }
  }
}

# ---------------------------------------------------------------- ⑤ 看板自检
Step '⑤ 看板自检（HTML 存在才跑）'
$html = Join-Path $root 'build\_perf\perf-report.html'
if (Test-Path $html) {
  $r = & node '_perf\check_report.js' 2>&1 | Out-String
  if ($LASTEXITCODE -eq 0) { Good '看板自检通过（几何 + 关键数字反查 + 采集时间行）' }
  else { Bad "看板自检失败（exit=$LASTEXITCODE）"; $r -split "`n" | Select-Object -Last 12 | ForEach-Object { Note $_.Trim() } }
} else {
  Note '没有看板 HTML，跳过（node _perf\make_report.js 生成）'
}

# ---------------------------------------------------------------- ⑥ 文档命令一致性
# 为什么需要它（2026-09-24）：上面五步都在比"文档里写的**数字**"，但没有一步在比
# "文档里写的**命令**"。于是一整类漂移可以活很久而不红：
#   · 文档点名了一个诊断程序，而它根本没进 CMake 的构建表（照文档复现会卡住）——
#     实测抓到 `_verify/residual_audit.cpp` 等 3 个；
#   · 文档给的参数在当前程序里不存在/写错（实测抓到 `--stiffness` 漏写导致整行数字对不上）。
# 做法：把文档里出现的 `xxx.exe` 与它的 `--flag` 抽出来，比对
#   ① 该程序必须在 CMakeLists.txt 的 add_executable 表里；② 每个 --flag 必须在该程序接受的表内。
# **本步查不出的**（写在文档里而不是藏起来）：漏写 `$env:PD_SHEAR` 这类"命令对但环境缺一句"——
# 那类只能靠口径声明，见 docs/convergence-speed-handoff.md §3 的警告。
Step '⑥ 文档命令一致性：程序必须已构建、参数必须被接受'

$allowedFlags = @{
  'pd_bench'       = @('--grid','--steps','--iters','--stiffness','--residual-tol','--tol','--preset','--threads')
  'pd_viewer'      = @('--grid','--iters','--stiffness','--residual-tol','--tol','--preset','--shot',
                       '--pin-corners','--pin-single','--frames','--damping','--no-early-exit','--dt','--size','--pos')
  'pd_restrace'    = @('--grid','--iters','--steps','--settle-iters','--stiffness','--spacing','--density',
                       '--soft-mass','--dt','--damping','--traj','--traj-steps','--traj-ref','--pin')
  'pd_itspectrum'  = @('--grid','--spacing','--density','--stiffness','--dt','--damping','--pin','--steps','--iters')
  'pd_sharefactor' = @('--grid','--stiffness','--reps','--threads','--flush-mb')
  'pd_solvecomp'   = @('--grid','--reps','--threads')
}

# CMake 建出来的可执行文件表（含别名：依赖缺失时会被跳过，但那属于构建期警告，不在这里报）
$cmakeText = [IO.File]::ReadAllText((Join-Path $root 'CMakeLists.txt'), [Text.Encoding]::UTF8)
$builtExes = @{}
foreach ($m in [regex]::Matches($cmakeText, 'add_executable\(\s*(\w+)')) { $builtExes[$m.Groups[1].Value] = $true }

$docSources = @('README.md') +
  (Get-ChildItem (Join-Path $root 'docs')   -Filter *.md -ErrorAction SilentlyContinue | ForEach-Object { "docs\$($_.Name)" }) +
  (Get-ChildItem (Join-Path $root '_perf')  -Filter *.md -ErrorAction SilentlyContinue | ForEach-Object { "_perf\$($_.Name)" })

# ── **来源豁免表**：这些 .exe 被文档引用，但**本来就不该**出现在 CMakeLists.txt 里。
# 第 ⑥ 步只对"应该由本项目构建"的程序报红；其余按来源分类豁免，避免把"正确的历史记录"判成漂移。
# 判据不是"漏了就豁免" —— 每条都写了它是**哪来的**。新增文档引用外部工具时在此加一行。
$externallyBuilt = @{
  # 外部工具链（不是本项目的目标）
  'glslc'                  = 'Vulkan SDK 的着色器编译器（docs/environment.md）'
  'ninja'                  = '构建工具本身（docs/environment.md）'
  'vk_probe'               = 'tools/ 下手工编译的 Vulkan 探测程序（docs/environment.md）'
  'x'                      = '论文源码里的工程名，不是本项目目标（handoff §4）'
  # 论文发布代码的本地 harness（%TEMP%\chebsim\，见 handoff §5，不进版本控制）
  'cheb_test'              = '论文源码副本的 harness（handoff §5.2）'
  # 存档基线二进制：	ools/check.ps1 第 ④ 步正是拿它做逐字比对
  'pd_bench_pre_fuse'      = '存档基线（build\_baseline\，第 ④ 步的比对对象），不重新构建'
  # 历史性能对照二进制：一次性测量用，见 docs/perf.md 的记录
  'pd_bench_a_csr'         = 'perf.md 的历史 A/B 二进制'
  'pd_bench_ctrl_single'   = 'perf.md 的历史 A/B 二进制'
  'pd_bench_phase2'        = 'perf.md 的历史 A/B 二进制'
  'pd_bench_pre_comp'      = 'perf.md 的历史 A/B 二进制'
  'pd_bench_pre_sharefactor' = 'perf.md 的历史 A/B 二进制'
}
$refExes = @{}      # exe -> 首次出现的文件
$badFlags = @{}
foreach ($rel in $docSources) {
  $p = Join-Path $root $rel
  if (-not (Test-Path $p)) { continue }
  $text = [IO.File]::ReadAllText($p, [Text.Encoding]::UTF8)
  foreach ($m in [regex]::Matches($text, '([A-Za-z_][A-Za-z0-9_]*)\.exe')) {
    $name = $m.Groups[1].Value
    if (-not $refExes.ContainsKey($name)) { $refExes[$name] = $rel }
  }
  # 抽参数：只取"程序名之后、行尾之前"的 --word，且只对该程序有表时才检查
  foreach ($cm in [regex]::Matches($text, '([A-Za-z_][A-Za-z0-9_]*)\.exe([^\r\n]*)')) {
    $name = $cm.Groups[1].Value
    if (-not $allowedFlags.ContainsKey($name)) { continue }
    if ($externallyBuilt.ContainsKey($name)) { continue }   # 豁免项不查参数
    foreach ($fm in [regex]::Matches($cm.Groups[2].Value, '(?<![\w-])(--[A-Za-z][A-Za-z0-9-]*)')) {
      $flag = $fm.Groups[1].Value
      if ($allowedFlags[$name] -notcontains $flag) {
        $badFlags["$rel :: $name $flag"] = $true
      }
    }
  }
}

$missingExe = @($refExes.Keys | Where-Object { -not $builtExes.ContainsKey($_) -and -not $externallyBuilt.ContainsKey($_) } | Sort-Object)
if ($missingExe.Count -eq 0) {
  Good "文档引用的 $($refExes.Count) 个程序全部在 CMake 构建表内（另有 $($externallyBuilt.Count) 个按来源豁免）"
} else {
  foreach ($n in $missingExe) {
    Bad "文档（$($refExes[$n])）引用了 $n.exe，但 CMakeLists.txt 没有它 ⇒ 照文档复现会卡住"
  }
}
if ($badFlags.Count -eq 0) {
  Good '文档里出现的 --flag 全部在该程序接受的表内'
} else {
  foreach ($k in ($badFlags.Keys | Sort-Object)) {
    Bad "文档参数对不上：$k（该程序不接受它，或文档写错了名字）"
  }
}

if ($script:Warn.Count) { foreach ($w in $script:Warn) { Write-Host "警告：$w" -ForegroundColor Yellow } }
if ($script:Fail.Count) {
  Write-Host "结论：**未通过**（$($script:Fail.Count) 项失败）" -ForegroundColor Red
  $script:Fail | ForEach-Object { Write-Host "  - $_" -ForegroundColor Red }
  exit 1
}
Write-Host '结论：全部通过' -ForegroundColor Green
exit 0
