HIP2LInterActomics - Distribuicao Windows
========================================

Conteudo
--------

  install_hip2linteractomics.ps1  Instalador para maquina Windows nova
  install_hip2linteractomics.bat  Wrapper para executar o PowerShell installer
  run_gui.bat                     Launcher da GUI

Instalacao Recomendada
----------------------

Abra PowerShell ou Anaconda Prompt na raiz do repositorio e rode:

  .\dist\windows\install_hip2linteractomics.ps1

Se a politica de execucao bloquear scripts, use:

  powershell.exe -NoProfile -ExecutionPolicy Bypass -File .\dist\windows\install_hip2linteractomics.ps1

O instalador:

  - detecta Conda/Miniforge/Miniconda;
  - instala Miniforge em %USERPROFILE%\.hip2linteractomics\miniforge3 se Conda nao existir;
  - cria ou atualiza luna-gui;
  - cria ou atualiza luna-env;
  - instala LUNA, RDKit, Open Babel, PyMOL e dependencias analiticas;
  - aplica o patch de compatibilidade usado pela GUI;
  - cria um atalho HIP2LInterActomics.cmd na area de trabalho.

Uso
---

  .\dist\windows\run_gui.bat

Snap No Windows
---------------

O formato .snap nao e um instalador nativo do Windows. Para Windows, use o
instalador PowerShell acima. Se voce precisa do pacote .snap ou de paralelismo
real com nproc > 1, use Linux ou WSL2.

Paralelismo No Windows
----------------------

No Windows nativo, o LUNA usa multiprocessing em modo spawn. Esse modo causa
falhas conhecidas em fluxos internos do LUNA quando nproc > 1. Por estabilidade,
a GUI e o script terminal limitam nproc para 1 no Windows nativo.

Para usar varios nucleos:

  - rode a versao Linux em uma maquina Linux;
  - ou use WSL2/Ubuntu e execute o instalador Linux.
