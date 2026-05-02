@echo off
setlocal
REM ============================================================
REM  HIP2LInterActomics_GUI launcher - Windows
REM ============================================================
REM  Sets a clean PATH so the luna-gui env's DLLs are picked up
REM  instead of any other conda/Qt install on the system (avoids
REM  conflicts with Schrodinger PyMOL3, Miniconda3, etc.).
REM
REM  Prerequisites (one-time):
REM    conda create -n luna-gui python=3.11 -y
REM    conda activate luna-gui
REM    pip install -r ..\..\requirements.txt
REM ============================================================

REM Resolve repo root (two levels above this script)
for %%I in ("%~dp0..\..") do set "REPO_ROOT=%%~fI"

REM Preferred explicit Python path for the GUI env

REM Allow user override
set "PYTHON_EXE="
if not "%HIP2LINTERACTOMICS_GUI_PYTHON%"=="" set "PYTHON_EXE=%HIP2LINTERACTOMICS_GUI_PYTHON%"
if not defined PYTHON_EXE if not "%LUNA_GUI_PYTHON%"=="" set "PYTHON_EXE=%LUNA_GUI_PYTHON%"

REM If already running from an activated luna-gui env, reuse it.
if not defined PYTHON_EXE (
    if /I "%CONDA_DEFAULT_ENV%"=="luna-gui" (
        if exist "%CONDA_PREFIX%\python.exe" set "PYTHON_EXE=%CONDA_PREFIX%\python.exe"
    )
)

if not defined PYTHON_EXE (
    for %%P in (
        "%USERPROFILE%\.conda\envs\luna-gui\python.exe"
        "%USERPROFILE%\miniconda3\envs\luna-gui\python.exe"
        "%USERPROFILE%\anaconda3\envs\luna-gui\python.exe"
        "%USERPROFILE%\miniforge3\envs\luna-gui\python.exe"
        "%LOCALAPPDATA%\Schrodinger\PyMOL3\envs\luna-gui\python.exe"
        "%LOCALAPPDATA%\miniconda3\envs\luna-gui\python.exe"
        "%LOCALAPPDATA%\anaconda3\envs\luna-gui\python.exe"
        "%LOCALAPPDATA%\miniforge3\envs\luna-gui\python.exe"
        "C:\ProgramData\miniconda3\envs\luna-gui\python.exe"
        "C:\ProgramData\anaconda3\envs\luna-gui\python.exe"
        "C:\ProgramData\miniforge3\envs\luna-gui\python.exe"
    ) do (
        if exist "%%~P" (
            set "PYTHON_EXE=%%~P"
            goto :python_found
        )
    )
)

if not defined PYTHON_EXE (
    echo [ERRO] python do env luna-gui nao encontrado.
    echo Defina HIP2LINTERACTOMICS_GUI_PYTHON ou ative o env luna-gui antes de chamar este launcher.
    pause
    exit /b 1
)

if not exist "%PYTHON_EXE%" (
    echo [ERRO] python informado para o env luna-gui nao existe:
    echo %PYTHON_EXE%
    echo Defina HIP2LINTERACTOMICS_GUI_PYTHON para o python.exe correto do env luna-gui.
    pause
    exit /b 1
)

:python_found
for %%I in ("%PYTHON_EXE%") do set "ENV_DIR=%%~dpI"
if "%ENV_DIR:~-1%"=="\" set "ENV_DIR=%ENV_DIR:~0,-1%"

set "PATH=%ENV_DIR%;%ENV_DIR%\Library\mingw-w64\bin;%ENV_DIR%\Library\usr\bin;%ENV_DIR%\Library\bin;%ENV_DIR%\Scripts;%ENV_DIR%\bin;C:\Windows\System32;C:\Windows;C:\Windows\System32\Wbem"

"%PYTHON_EXE%" "%REPO_ROOT%\run.py"
pause
