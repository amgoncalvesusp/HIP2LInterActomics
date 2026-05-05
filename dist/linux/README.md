# HIP²LInterActomics — Distribuição Linux

## Conteúdo
- `run_gui.sh` — launcher principal no Linux; ativa o env `luna-gui`, valida dependências da GUI e inicia a aplicação

## Instalação (uma única vez)

```bash
# 1) Instale Miniconda se ainda não tiver
# https://docs.conda.io/en/latest/miniconda.html

# 2) Crie o env da GUI
conda create -n luna-gui python=3.11 -y
conda activate luna-gui
pip install -r ../../requirements.txt

# 3) Torne o launcher executável
chmod +x run_gui.sh
```

## Uso

```bash
./run_gui.sh
```

Na primeira execução, vá à aba **1. Setup** e clique em **Instalar LUNA**
para criar o ambiente `luna-env` (onde o LUNA de fato roda — é separado
do `luna-gui`).

O launcher Linux agora também verifica se o env da GUI contém:
- `PyQt6`
- `matplotlib`
- `numpy`
- `scipy`

Se faltar alguma dependência, ele falha cedo e informa o comando correto.

## Dependências do sistema (apt / dnf)

Para o PyQt6 e o PyMOL funcionarem corretamente, instale os pacotes
gráficos padrão:

```bash
# Debian / Ubuntu
sudo apt install libxcb-cursor0 libxkbcommon-x11-0 libgl1 libegl1

# Fedora / RHEL
sudo dnf install xcb-util-cursor libxkbcommon-x11 mesa-libGL mesa-libEGL
```

## Vantagens do Linux sobre Windows

- **Paralelização funcional**: use `nproc` = número de núcleos disponíveis
  na aba Executar. O bug de pickling do `spawn` no Windows não afeta o
  Linux (usa `fork`).
- **Melhor performance** para bibliotecas grandes (milhares de ligantes).
- **Sem conflitos de DLL** entre múltiplas instalações de conda.
- **Melhor fluxo de visualização**: a GUI exporta gráficos em `PNG`, `SVG` e `PDF`, e também exporta clusters em `CSV`.

## Variável de ambiente opcional

Se você usa um nome diferente para o env:

```bash
HIP2LINTERACTOMICS_GUI_ENV=meu-env ./run_gui.sh
```

Se preferir apontar diretamente para um Python específico:

```bash
HIP2LINTERACTOMICS_GUI_PYTHON=/caminho/para/python ./run_gui.sh
```
