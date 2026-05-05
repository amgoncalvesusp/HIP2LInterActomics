HIP²LInterActomics — Distribuição Windows
================================

Conteúdo:
  run_gui.bat  — launcher que isola o PATH e chama o Python correto do env luna-gui

Instalação (uma única vez):
  1) Instale o Miniconda: https://docs.conda.io/en/latest/miniconda.html
  2) Abra um PowerShell / Prompt e rode:
       conda create -n luna-gui python=3.11 -y
       conda activate luna-gui
       pip install -r ..\..\requirements.txt
  3) O launcher tenta detectar o python.exe do env luna-gui automaticamente.
     Se precisar forçar outro caminho, defina a variável:
       set LUNA_GUI_PYTHON=C:\caminho\para\python.exe

Uso:
  Duplo-clique em run_gui.bat  — ou rode pelo terminal:
       .\run_gui.bat

A primeira vez que abrir, vá à aba "1. Setup" e clique em
"Instalar LUNA" para criar o ambiente luna-env (onde o LUNA de
fato roda; é separado do luna-gui).

Limitações conhecidas no Windows:
  - nproc (paralelização) fica travado em 1 devido a um bug do
    LUNA no modo spawn do multiprocessing do Windows. Para rodar
    milhares de ligantes em paralelo, use Linux ou WSL2.
