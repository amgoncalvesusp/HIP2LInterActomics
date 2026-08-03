@echo off
start "" /wait "%~dp0..\HIP2LInterActomics.exe" --multiple-run %*
exit /b %errorlevel%
