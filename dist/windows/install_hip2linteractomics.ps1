<#
Instala HIP2LInterActomics em uma maquina Windows limpa.

O instalador cria dois ambientes:
  - luna-gui: executa a interface grafica.
  - luna-env: executa LUNA, RDKit, Open Babel e PyMOL.

Snap nao e um formato nativo do Windows. Para Windows, use este instalador
PowerShell ou rode a versao Linux/Snap dentro de WSL2.
#>

[CmdletBinding()]
param(
    [string]$GuiEnv = "luna-gui",
    [string]$LunaEnv = "luna-env",
    [string]$CondaRoot = "$env:USERPROFILE\.hip2linteractomics\miniforge3",
    [switch]$GuiOnly,
    [switch]$SkipLuna,
    [switch]$NoShortcut
)

Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

$ScriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
$RepoRoot = Resolve-Path (Join-Path $ScriptDir "..\..")

function Write-Step {
    param([string]$Message)
    Write-Host ""
    Write-Host "[HIP2L] $Message" -ForegroundColor Cyan
}

function Resolve-Conda {
    $candidates = @()
    foreach ($varName in @("HIP2LINTERACTOMICS_GUI_CONDA", "LUNA_GUI_CONDA", "CONDA_EXE", "MAMBA_EXE")) {
        $value = [Environment]::GetEnvironmentVariable($varName)
        if ($value) { $candidates += $value }
    }

    $cmd = Get-Command conda -ErrorAction SilentlyContinue
    if ($cmd) { $candidates += $cmd.Source }

    $candidates += @(
        (Join-Path $CondaRoot "Scripts\conda.exe"),
        "$env:USERPROFILE\miniforge3\Scripts\conda.exe",
        "$env:USERPROFILE\miniconda3\Scripts\conda.exe",
        "$env:USERPROFILE\anaconda3\Scripts\conda.exe",
        "$env:LOCALAPPDATA\miniforge3\Scripts\conda.exe",
        "$env:LOCALAPPDATA\miniconda3\Scripts\conda.exe",
        "$env:LOCALAPPDATA\anaconda3\Scripts\conda.exe",
        "$env:ProgramData\miniforge3\Scripts\conda.exe",
        "$env:ProgramData\miniconda3\Scripts\conda.exe",
        "$env:ProgramData\anaconda3\Scripts\conda.exe"
    )

    foreach ($candidate in $candidates) {
        if ($candidate -and (Test-Path -LiteralPath $candidate)) {
            return (Resolve-Path -LiteralPath $candidate).Path
        }
    }
    return $null
}

function Install-Miniforge {
    New-Item -ItemType Directory -Force -Path (Split-Path -Parent $CondaRoot) | Out-Null
    $downloadDir = Join-Path $env:TEMP "hip2linteractomics-installer"
    New-Item -ItemType Directory -Force -Path $downloadDir | Out-Null
    $installer = Join-Path $downloadDir "Miniforge3-Windows-x86_64.exe"
    $url = "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Windows-x86_64.exe"

    Write-Step "Conda nao encontrado. Baixando Miniforge: $url"
    Invoke-WebRequest -Uri $url -OutFile $installer

    Write-Step "Instalando Miniforge em $CondaRoot"
    $args = @(
        "/S",
        "/InstallationType=JustMe",
        "/RegisterPython=0",
        "/AddToPath=0",
        "/D=$CondaRoot"
    )
    $process = Start-Process -FilePath $installer -ArgumentList $args -Wait -PassThru
    if ($process.ExitCode -ne 0) {
        throw "Falha ao instalar Miniforge. Exit code: $($process.ExitCode)"
    }
}

function Invoke-Conda {
    param([Parameter(ValueFromRemainingArguments = $true)][string[]]$Args)
    & $script:CondaExe @Args
    if ($LASTEXITCODE -ne 0) {
        throw "Conda falhou: $($Args -join ' ')"
    }
}

function Test-CondaEnv {
    param([string]$Name)
    & $script:CondaExe run -n $Name python -V *> $null
    return $LASTEXITCODE -eq 0
}

function Conda-CreateOrInstall {
    param(
        [string]$Name,
        [string[]]$Packages
    )
    if (Test-CondaEnv $Name) {
        Invoke-Conda install -n $Name --override-channels -c conda-forge -y @Packages
    } else {
        Invoke-Conda create -n $Name --override-channels -c conda-forge -y @Packages
    }
}

function Install-GuiEnv {
    Write-Step "Criando/atualizando ambiente $GuiEnv"
    Conda-CreateOrInstall -Name $GuiEnv -Packages @("python=3.11", "pip")
    Invoke-Conda run -n $GuiEnv python -m pip install --upgrade pip
    Invoke-Conda run -n $GuiEnv python -m pip install -r (Join-Path $RepoRoot "requirements.txt")
}

function Install-LunaEnv {
    Write-Step "Criando/atualizando ambiente $LunaEnv"
    Conda-CreateOrInstall -Name $LunaEnv -Packages @(
        "python=3.9", "pip", "rdkit", "openbabel", "pymol-open-source",
        "biopython=1.79", "numpy", "pandas", "scipy", "scikit-learn",
        "matplotlib", "seaborn", "networkx"
    )
    Invoke-Conda run -n $LunaEnv python -s -m pip install pdbecif "mmh3<4" xopen colorlog
    Invoke-Conda run -n $LunaEnv python -s -m pip install --no-build-isolation -U luna
    Invoke-Conda run -n $LunaEnv python -s (Join-Path $RepoRoot "luna_gui\core\_luna_patch.py")
}

function New-Launcher {
    $desktop = [Environment]::GetFolderPath("Desktop")
    if (-not $desktop) { return }
    $launcher = Join-Path $desktop "HIP2LInterActomics.cmd"
    $runBat = Join-Path $RepoRoot "dist\windows\run_gui.bat"
    $guiPython = (& $script:CondaExe run -n $GuiEnv python -c "import sys; print(sys.executable)" | Select-Object -Last 1)
    $content = @(
        "@echo off",
        "set HIP2LINTERACTOMICS_GUI_CONDA=$script:CondaExe",
        "set HIP2LINTERACTOMICS_GUI_PYTHON=$guiPython",
        "call `"$runBat`""
    )
    Set-Content -Path $launcher -Value $content -Encoding ASCII
    Write-Step "Atalho criado em $launcher"
}

if (-not (Test-Path -LiteralPath (Join-Path $RepoRoot "requirements.txt"))) {
    throw "Rode este instalador a partir de uma copia completa do repositorio."
}

$script:CondaExe = Resolve-Conda
if (-not $script:CondaExe) {
    Install-Miniforge
    $script:CondaExe = Join-Path $CondaRoot "Scripts\conda.exe"
}

if (-not (Test-Path -LiteralPath $script:CondaExe)) {
    throw "Conda nao encontrado apos instalacao: $script:CondaExe"
}

Write-Step "Usando conda: $script:CondaExe"
Invoke-Conda config --set channel_priority flexible

Install-GuiEnv
if (-not $GuiOnly -and -not $SkipLuna) {
    Install-LunaEnv
} else {
    Write-Step "Instalacao do luna-env pulada por opcao do usuario."
}

if (-not $NoShortcut) {
    New-Launcher
}

Write-Step "Instalacao concluida."
Write-Host "Abra com:"
Write-Host "  $RepoRoot\dist\windows\run_gui.bat"
Write-Host ""
Write-Host "Paralelismo no Windows nativo: nproc fica limitado a 1 por estabilidade do LUNA."
Write-Host "Para usar varios nucleos, execute a versao Linux em WSL2 ou em uma maquina Linux."
