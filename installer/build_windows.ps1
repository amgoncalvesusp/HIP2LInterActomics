<# Build the native Windows bundle and, when available, the Inno Setup installer. #>

[CmdletBinding()]
param(
    [string]$Python = "python",
    [switch]$InstallBuildDependencies,
    [string]$OutputRoot = "build\release-windows"
)

Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

if ($env:OS -ne "Windows_NT") {
    throw "Este build precisa ser executado no Windows."
}

$RepoRoot = Resolve-Path (Join-Path $PSScriptRoot "..")
$OutputRoot = [System.IO.Path]::GetFullPath((Join-Path $RepoRoot $OutputRoot))
$WorkDir = Join-Path $OutputRoot "pyinstaller-work"
$BundleRoot = Join-Path $OutputRoot "bundle"
$InstallerRoot = Join-Path $OutputRoot "installer"
$OriginalPath = $env:PATH

New-Item -ItemType Directory -Force -Path $OutputRoot, $WorkDir, $BundleRoot, $InstallerRoot | Out-Null
Push-Location $RepoRoot
try {
    if ($InstallBuildDependencies) {
        & $Python -m pip install --upgrade pip
        & $Python -m pip install -r requirements.txt "pyinstaller>=6.10,<7"
        if ($LASTEXITCODE -ne 0) { throw "Falha ao instalar dependencias de build." }
    }

    & $Python -c "import PyQt6, jinja2, matplotlib, numpy, reportlab, scipy, sklearn, PyInstaller"
    if ($LASTEXITCODE -ne 0) {
        throw "Dependencias ausentes. Rode novamente com -InstallBuildDependencies."
    }

    # Conda stores OpenSSL, libffi and expat in Library\bin. Make that directory
    # visible while PyInstaller resolves transitive dependencies of Python's
    # extension modules; the resolved DLLs are then copied into the bundle.
    $PythonPrefix = (& $Python -c "import sys; print(sys.prefix)").Trim()
    $CondaLibraryBin = Join-Path $PythonPrefix "Library\bin"
    if (Test-Path -LiteralPath $CondaLibraryBin) {
        $env:PATH = "$CondaLibraryBin;$PythonPrefix;$(Join-Path $PythonPrefix 'DLLs');$env:PATH"
    }

    & $Python -m PyInstaller --noconfirm --clean --workpath $WorkDir --distpath $BundleRoot HIP2LInterActomics.spec
    if ($LASTEXITCODE -ne 0) { throw "PyInstaller falhou." }

    $Exe = Join-Path $BundleRoot "HIP2LInterActomics\HIP2LInterActomics.exe"
    if (-not (Test-Path -LiteralPath $Exe)) { throw "Executavel nao gerado: $Exe" }
    $EnvironmentFile = Join-Path (Split-Path -Parent $Exe) "environment.yml"
    if (-not (Test-Path -LiteralPath $EnvironmentFile)) {
        throw "environment.yml nao foi incluido no bundle: $EnvironmentFile"
    }
    foreach ($CliFile in @("hipplinteractomics_terminal.py", "hipplinteractomics_multiple_run.py")) {
        $CliPath = Join-Path (Split-Path -Parent $Exe) $CliFile
        if (-not (Test-Path -LiteralPath $CliPath)) {
            throw "$CliFile nao foi incluido no bundle: $CliPath"
        }
    }

    $IsccCandidates = @(
        "$env:LOCALAPPDATA\Programs\Inno Setup 6\ISCC.exe",
        "$env:ProgramFiles(x86)\Inno Setup 6\ISCC.exe",
        "$env:ProgramFiles\Inno Setup 6\ISCC.exe"
    )
    $Iscc = $IsccCandidates | Where-Object { $_ -and (Test-Path -LiteralPath $_) } | Select-Object -First 1
    if ($Iscc) {
        & $Iscc "/DBundleDir=$(Split-Path -Parent $Exe)" "/DOutputDir=$InstallerRoot" "installer\HIP2LInterActomics.iss"
        if ($LASTEXITCODE -ne 0) { throw "Inno Setup falhou." }
    } else {
        Write-Warning "Inno Setup 6 nao encontrado; o bundle .exe foi gerado, mas o Setup.exe nao."
    }

    Write-Host "Executavel: $Exe"
    Write-Host "Ambiente terminal-only: $EnvironmentFile"
    if ($Iscc) { Write-Host "Instalador: $(Join-Path $InstallerRoot 'HIP2LInterActomics-Setup.exe')" }
}
finally {
    $env:PATH = $OriginalPath
    Pop-Location
}
