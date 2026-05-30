param(
  [ValidateSet(
    "Setting_I.r",
    "Setting_II_df=3.r",
    "Setting_II_df=5.r",
    "Setting_III.r",
    "Setting_IV.r",
    "Setting_SUPP.R"
  )]
  [string]$Setting = "Setting_I.r",

  [int]$Rep = 1,
  [int]$BootstrapRuns = 0,
  [int]$MaxAttempts = 5,
  [int]$BootstrapMaxAttempts = 5,
  [double]$BiasCutoff = 15,
  [int]$Seed = 1,
  [string]$OutputDir = "",
  [string]$LogFile = "",
  [switch]$RunBootstrap
)

$ErrorActionPreference = "Stop"

$scriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
$repoRoot = Resolve-Path (Join-Path $scriptDir "..")

if ($OutputDir -eq "") {
  $OutputDir = Join-Path $repoRoot "simulation_output"
}

$requiredPackages = @(
  "dplyr",
  "nlme",
  "Deriv",
  "stringr",
  "LaplacesDemon",
  "purrr",
  "MASS",
  "mvtnorm",
  "Matrix",
  "here",
  "R.utils",
  "xtable"
)

$packageCsv = ($requiredPackages -join ",")
$checkExpr = @"
pkgs <- strsplit('$packageCsv', ',')[[1]]
missing <- pkgs[!(pkgs %in% rownames(installed.packages()))]
if (length(missing)) {
  cat(paste(missing, collapse='\n'))
  quit(status=2)
}
"@

$missingOutput = & Rscript -e $checkExpr 2>&1
if ($LASTEXITCODE -eq 2) {
  Write-Error ("Missing R packages:`n" + ($missingOutput -join "`n") + "`nInstall them before running the simulation.")
}
if ($LASTEXITCODE -ne 0) {
  Write-Error ($missingOutput -join "`n")
}

New-Item -ItemType Directory -Force -Path $OutputDir | Out-Null

if ($LogFile -eq "") {
  $timestamp = Get-Date -Format "yyyyMMdd_HHmmss"
  $logName = "{0}_{1}.log" -f ([IO.Path]::GetFileNameWithoutExtension($Setting)), $timestamp
  $LogFile = Join-Path $OutputDir $logName
}

$env:SIM_REP = [string]$Rep
$env:SIM_SEED = [string]$Seed
$env:SIM_K_RUNS = [string]$BootstrapRuns
$env:SIM_MAX_ATTEMPTS = [string]$MaxAttempts
$env:SIM_BOOTSTRAP_MAX_ATTEMPTS = [string]$BootstrapMaxAttempts
$env:SIM_RBIAS_REMOVE_CUTOFF = [string]$BiasCutoff
$env:SIM_OUTPUT_DIR = (Resolve-Path $OutputDir).Path
$env:SIM_SKIP_BOOTSTRAP = if ($RunBootstrap) { "false" } else { "true" }

Write-Host "Repository: $repoRoot"
Write-Host "Setting: $Setting"
Write-Host "SIM_REP=$env:SIM_REP"
Write-Host "SIM_SEED=$env:SIM_SEED"
Write-Host "SIM_K_RUNS=$env:SIM_K_RUNS"
Write-Host "SIM_SKIP_BOOTSTRAP=$env:SIM_SKIP_BOOTSTRAP"
Write-Host "SIM_MAX_ATTEMPTS=$env:SIM_MAX_ATTEMPTS"
Write-Host "SIM_BOOTSTRAP_MAX_ATTEMPTS=$env:SIM_BOOTSTRAP_MAX_ATTEMPTS"
Write-Host "SIM_RBIAS_REMOVE_CUTOFF=$env:SIM_RBIAS_REMOVE_CUTOFF"
Write-Host "SIM_OUTPUT_DIR=$env:SIM_OUTPUT_DIR"
Write-Host "LogFile=$LogFile"

Push-Location $repoRoot
try {
  $elapsed = Measure-Command {
    $oldErrorActionPreference = $ErrorActionPreference
    $oldNativeErrorPreference = $null
    if (Get-Variable -Name PSNativeCommandUseErrorActionPreference -Scope Global -ErrorAction SilentlyContinue) {
      $oldNativeErrorPreference = $global:PSNativeCommandUseErrorActionPreference
      $global:PSNativeCommandUseErrorActionPreference = $false
    }
    $ErrorActionPreference = "Continue"
    try {
      & Rscript (Join-Path "simulation" $Setting) 2>&1 | Tee-Object -FilePath $LogFile
      $rExitCode = $LASTEXITCODE
    }
    finally {
      $ErrorActionPreference = $oldErrorActionPreference
      if ($null -ne $oldNativeErrorPreference) {
        $global:PSNativeCommandUseErrorActionPreference = $oldNativeErrorPreference
      }
    }
    if ($rExitCode -ne 0) {
      throw "Rscript failed with exit code $rExitCode"
    }
  }
  Write-Host ("ElapsedSeconds={0:n1}" -f $elapsed.TotalSeconds)
}
finally {
  Pop-Location
}
