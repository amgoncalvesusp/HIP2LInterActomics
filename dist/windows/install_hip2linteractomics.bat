@echo off
setlocal
REM Wrapper for the PowerShell installer.
powershell.exe -NoProfile -ExecutionPolicy Bypass -File "%~dp0install_hip2linteractomics.ps1" %*
if errorlevel 1 (
    echo.
    echo [ERRO] A instalacao falhou.
    pause
    exit /b 1
)
echo.
echo Instalacao concluida.
pause
