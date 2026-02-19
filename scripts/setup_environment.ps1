param(
  [string]$VenvDir = ".venv"
)

$ErrorActionPreference = "Stop"

if (-not (Get-Command python -ErrorAction SilentlyContinue)) {
  Write-Host "Python is not available on PATH." -ForegroundColor Red
  exit 1
}

python -m venv $VenvDir

$activateScript = Join-Path $VenvDir "Scripts/Activate.ps1"
if (-not (Test-Path $activateScript)) {
  Write-Host "Could not find venv activation script at $activateScript" -ForegroundColor Red
  exit 1
}

. $activateScript
python -m pip install --upgrade pip
pip install -r requirements.txt

Write-Host "Python environment ready at $VenvDir." -ForegroundColor Green
Write-Host "Activate with: $activateScript"
Write-Host "For simulation build dependencies (MPI/GSL), use WSL2 Ubuntu and run scripts/setup_environment.sh there." -ForegroundColor Yellow
