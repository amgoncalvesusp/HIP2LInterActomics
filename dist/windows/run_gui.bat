@echo off
setlocal
REM ============================================================
REM  LUNA GUI launcher - Windows
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
set "PYTHON_EXE=C:\Users\adria\AppData\Local\Schrodinger\PyMOL3\envs\luna-gui\python.exe"

REM Allow user override
if not "%LUNA_GUI_PYTHON%"=="" set "PYTHON_EXE=%LUNA_GUI_PYTHON%"

if not exist "%PYTHON_EXE%" (
    for %%P in (
        "C:\Users\adria\AppData\Local\Schrodinger\PyMOL3\envs\luna-gui\python.exe"
        "%USERPROFILE%\miniconda3\envs\luna-gui\python.exe"
        "%USERPROFILE%\anaconda3\envs\luna-gui\python.exe"
        "C:\ProgramData\miniconda3\envs\luna-gui\python.exe"
        "C:\ProgramData\anaconda3\envs\luna-gui\python.exe"
    ) do (
        if exist "%%~P" (
            set "PYTHON_EXE=%%~P"
            goto :python_found
        )
    )
    echo [ERRO] python do env luna-gui nao encontrado.
    echo Esperado em: %PYTHON_EXE%
    echo Defina LUNA_GUI_PYTHON ou ajuste este launcher para o caminho correto.
    pause
    exit /b 1
)

:python_found
for %%I in ("%PYTHON_EXE%") do set "ENV_DIR=%%~dpI"
if "%ENV_DIR:~-1%"=="\" set "ENV_DIR=%ENV_DIR:~0,-1%"

set "PATH=%ENV_DIR%;%ENV_DIR%\Library\mingw-w64\bin;%ENV_DIR%\Library\usr\bin;%ENV_DIR%\Library\bin;%ENV_DIR%\Scripts;%ENV_DIR%\bin;C:\Windows\System32;C:\Windows;C:\Windows\System32\Wbem"

"%PYTHON_EXE%" "%REPO_ROOT%\run.py"
pause
