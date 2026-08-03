@echo off
start "" /wait "%~dp0..\HIP2LInterActomics.exe" --terminal %*
exit /b %errorlevel%
