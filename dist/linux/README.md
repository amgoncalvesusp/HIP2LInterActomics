# HIP2LInterActomics - Distribuicao Linux

## Conteudo

- `install_hip2linteractomics.sh`: instalador para maquina Linux nova.
- `run_gui.sh`: launcher principal da GUI.
- `build_snap.sh`: gera o pacote `.snap` usando `snapcraft`.

## Instalacao Recomendada

Em uma copia completa do repositorio:

```bash
chmod +x dist/linux/install_hip2linteractomics.sh
./dist/linux/install_hip2linteractomics.sh
```

O instalador:

- detecta Conda/Miniforge/Miniconda;
- instala Miniforge em `~/.hip2linteractomics/miniforge3` se Conda nao existir;
- cria ou atualiza `luna-gui`;
- cria ou atualiza `luna-env`;
- instala LUNA, RDKit, Open Babel, PyMOL e dependencias analiticas;
- aplica o patch de compatibilidade usado pela GUI;
- cria um atalho desktop opcional.

Para instalar apenas a GUI e deixar o `luna-env` para a aba `1. Inicio`:

```bash
./dist/linux/install_hip2linteractomics.sh --gui-only
```

## Uso

```bash
./dist/linux/run_gui.sh
```

Variaveis opcionais:

```bash
HIP2LINTERACTOMICS_GUI_ENV=meu-env ./dist/linux/run_gui.sh
HIP2LINTERACTOMICS_GUI_PYTHON=/caminho/para/python ./dist/linux/run_gui.sh
HIP2LINTERACTOMICS_GUI_CONDA=/caminho/para/conda ./dist/linux/run_gui.sh
```

## Dependencias De Sistema

```bash
# Debian / Ubuntu
sudo apt install libxcb-cursor0 libxkbcommon-x11-0 libgl1 libegl1

# Fedora / RHEL
sudo dnf install xcb-util-cursor libxkbcommon-x11 mesa-libGL mesa-libEGL
```

## Snap

O formato `.snap` e exclusivo do Linux. Para gerar:

```bash
sudo snap install snapcraft --classic
./dist/linux/build_snap.sh
```

Instalacao local do arquivo gerado:

```bash
sudo snap install ./hip2linteractomics_0.1_amd64.snap --classic --dangerous
hip2linteractomics
```

O Snap empacota a GUI. O runtime pesado do LUNA continua no ambiente Conda
`luna-env`, criado pela aba `1. Inicio` ou pelo instalador Linux.

## Paralelismo

No Linux, `nproc` pode usar varios nucleos. Este e o ambiente recomendado para
bibliotecas grandes de ligantes e para analises longas.
