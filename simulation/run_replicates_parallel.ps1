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

  [int]$Replicates = 10,
  [int]$Jobs = 4,
  [int]$BootstrapRuns = 0,
  [int]$MaxAttempts = 5,
  [int]$BootstrapMaxAttempts = 5,
  [double]$BiasCutoff = 1000,
  [int]$BaseSeed = 1,
  [string]$OutputRoot = "",
  [switch]$RunBootstrap
)

$ErrorActionPreference = "Stop"

$scriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
$repoRoot = Resolve-Path (Join-Path $scriptDir "..")
$runner = Join-Path $scriptDir "run_simulation.ps1"

if ($OutputRoot -eq "") {
  $settingName = [IO.Path]::GetFileNameWithoutExtension($Setting)
  $timestamp = Get-Date -Format "yyyyMMdd_HHmmss"
  $OutputRoot = Join-Path $repoRoot ("simulation_output_parallel\{0}_{1}" -f $settingName, $timestamp)
}

New-Item -ItemType Directory -Force -Path $OutputRoot | Out-Null

$replicateJobs = @()
for ($i = 1; $i -le $Replicates; $i++) {
  while (($replicateJobs | Where-Object { $_.State -eq "Running" }).Count -ge $Jobs) {
    Start-Sleep -Seconds 5
    $completed = $replicateJobs | Where-Object { $_.State -ne "Running" -and -not $_.PSBeginTimeReported }
    foreach ($job in $completed) {
      $job | Add-Member -NotePropertyName PSBeginTimeReported -NotePropertyValue $true -Force
      Receive-Job -Job $job
      if ($job.State -ne "Completed") {
        throw "A replicate job failed. See job output above."
      }
    }
  }

  $repOutput = Join-Path $OutputRoot ("rep_{0:D3}" -f $i)
  $seed = $BaseSeed + $i - 1

  $argList = @(
    "-NoProfile",
    "-ExecutionPolicy", "Bypass",
    "-File", $runner,
    "-Setting", $Setting,
    "-Rep", "1",
    "-BootstrapRuns", [string]$BootstrapRuns,
    "-MaxAttempts", [string]$MaxAttempts,
    "-BootstrapMaxAttempts", [string]$BootstrapMaxAttempts,
    "-BiasCutoff", [string]$BiasCutoff,
    "-Seed", [string]$seed,
    "-OutputDir", $repOutput
  )
  if ($RunBootstrap) {
    $argList += "-RunBootstrap"
  }

  Write-Host "Starting replicate $i seed=$seed output=$repOutput"
  $replicateJobs += Start-Job -ScriptBlock {
    param($argsForPowerShell)
    & powershell @argsForPowerShell
    if ($LASTEXITCODE -ne 0) {
      throw "Replicate process failed with exit code $LASTEXITCODE"
    }
  } -ArgumentList (, $argList)
}

while (($replicateJobs | Where-Object { $_.State -eq "Running" }).Count -gt 0) {
  Start-Sleep -Seconds 5
  $completed = $replicateJobs | Where-Object { $_.State -ne "Running" -and -not $_.PSBeginTimeReported }
  foreach ($job in $completed) {
    $job | Add-Member -NotePropertyName PSBeginTimeReported -NotePropertyValue $true -Force
    Receive-Job -Job $job
    if ($job.State -ne "Completed") {
      throw "A replicate job failed. See job output above."
    }
  }
}

foreach ($job in $replicateJobs | Where-Object { -not $_.PSBeginTimeReported }) {
  Receive-Job -Job $job
  if ($job.State -ne "Completed") {
    throw "A replicate job failed. See job output above."
  }
}

if ($replicateJobs.Count -gt 0) {
  Remove-Job -Job $replicateJobs
}
Write-Host "Parallel replicates complete. OutputRoot=$OutputRoot"
