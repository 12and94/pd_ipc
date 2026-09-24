param(
  [switch]$DebugTools,
  [switch]$Clean,
  [string]$Config = "Release",
  [string]$BuildDir = "build"
)

# tools/build.ps1 -- 统一构建入口，把本机环境的坑固定下来。
#
#   .\tools\build.ps1                  # 配置 + Release 构建
#   .\tools\build.ps1 -DebugTools      # 额外构建 _debug/ 下的调试程序
#   .\tools\build.ps1 -Clean           # 先删除 build 目录
#
# 坑 1：源码是 UTF-8（含中文注释），必须让编译器按 UTF-8 读取（/utf-8）。
#       否则某些汉字的字节序列会吃掉行尾的 "*/"，报出莫名其妙的语法错误。
#       CMake 里已加 /utf-8；这里同时设置 CL 环境变量，覆盖手工编译的场景。
# 坑 2：PATH 上的 MinGW g++ 8.1 不支持 -std=c++20，必须用 MSVC。
# 坑 3：本机没有 Ninja、没有 vswhere，只能用 "Visual Studio 17 2022" 生成器。

$ErrorActionPreference = "Stop"
$root = Split-Path -Parent $PSScriptRoot
Set-Location $root

$env:CL = "/utf-8"
$env:VSLANG = "1033"   # 让 MSVC 输出英文诊断，便于检索

if ($Clean -and (Test-Path $BuildDir)) {
  Write-Host "clean $BuildDir ..." -ForegroundColor Yellow
  Remove-Item -Recurse -Force $BuildDir
}

$configureArgs = @("-S", ".", "-B", $BuildDir, "-G", "Visual Studio 17 2022", "-A", "x64")
if ($DebugTools) { $configureArgs += "-DPD_BUILD_DEBUG_TOOLS=ON" }

Write-Host "configure ..." -ForegroundColor Cyan
& cmake @configureArgs
if ($LASTEXITCODE -ne 0) { throw "cmake configure failed" }

Write-Host "build ($Config) ..." -ForegroundColor Cyan
& cmake --build $BuildDir --config $Config
if ($LASTEXITCODE -ne 0) { throw "build failed" }

Write-Host "done. executables in $BuildDir\$Config\" -ForegroundColor Green
Get-ChildItem "$BuildDir\$Config" -Filter *.exe | ForEach-Object { "  $($_.Name)" }
