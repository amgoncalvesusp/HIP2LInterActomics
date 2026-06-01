# HIP²LInterActomics

`HIP²LInterActomics` é uma interface gráfica e também um fluxo executável por terminal para preparar, calcular, organizar e interpretar interações intermoleculares usando o LUNA. O aplicativo foi pensado para análise racional de triagem virtual, poses de docking, complexos hidratados e trajetórias de dinâmica molecular representadas como frames ou conformações.

O programa transforma arquivos estruturais em tabelas, mapas de calor, fingerprints de interação, sessões PyMOL, análises de importância de features e relatórios exportáveis. A ideia central é permitir que o usuário saia de um conjunto de complexos proteína-ligante ou proteína-proteína e chegue a uma leitura química interpretável: quais contatos aparecem, onde aparecem, quais resíduos participam, quais fingerprints são mais informativos e quais ligantes compartilham padrões de interação.

## O Que O Aplicativo Faz

- Prepara complexos de docking separando proteína, ligante e moléculas de água.
- Preserva águas estruturais quando a análise hidratada está ativada.
- Mantém numeração de resíduos e cadeias sempre que possível durante conversões.
- Executa o LUNA para calcular interações intermoleculares.
- Gera fingerprints de interação `HIFP`, `EIFP` e `FIFP`.
- Calcula matriz de similaridade de Tanimoto entre ligantes.
- Agrupa ligantes por perfis de interação.
- Filtra sessões PyMOL por binding modes ou tipos de interação.
- Gera mapas de calor ligante x resíduo e mapas completos com múltiplos tipos de interação por célula.
- Calcula importância de fingerprints por modelos supervisionados de classificação ou regressão.
- Atribui classes e níveis de shell aos fingerprints, tratando colisões e features não confiáveis.
- Gera gráficos e relatórios HTML/PDF para inspeção e documentação.
- Oferece um script terminal para rodar o mesmo fluxo por arquivo de configuração.

## Fluxo De Trabalho

1. Configure o ambiente Conda e instale o `luna-env` pela aba `1. Início`.
2. Informe proteína, ligantes, água e modo de trajetória/poses na aba `2. Projeto`.
3. Escolha fingerprints, similaridade, rótulos, filtros e configurações LUNA na aba `3. Análises`.
4. Execute o cálculo na aba `4. Executar`.
5. Explore tabelas, mapas de calor, clusters, análises FP e sessões PyMOL na aba `5. Resultados`.
6. Reabra projetos já calculados usando o `workdir` salvo.

## Estrutura Do Repositório

```text
luna_gui/
  luna_gui/                         Código-fonte da interface e dos runtimes
  tests/                            Testes automatizados principais
  dist/linux/run_gui.sh             Launcher Linux
  dist/windows/run_gui.bat          Launcher Windows
  hipplinteractomics_terminal.py    Execução completa por arquivo de configuração
  requirements.txt                  Dependências da GUI
  run.py                            Launcher simples: python run.py
```

Pastas de resultados, workdirs de teste e arquivos temporários gerados durante análises podem ser grandes e não devem ser versionados.

## Instalação Rápida

O projeto usa dois ambientes separados:

- `luna-gui`: ambiente da interface gráfica.
- `luna-env`: ambiente interno onde o LUNA é instalado e executado.

Essa separação evita conflitos entre PyQt6, matplotlib, PyMOL, OpenBabel e as dependências do LUNA.

### 1. Instale Miniconda Ou Anaconda

Instale Miniconda ou Anaconda e confirme que o comando `conda` funciona:

```bash
conda --version
```

### 2. Crie O Ambiente Da GUI

No Windows, abra PowerShell ou Anaconda Prompt:

```powershell
cd D:\projeto_LUNA_GUIDE\luna_gui
conda create -n luna-gui python=3.11 -y
conda activate luna-gui
pip install -r requirements.txt
python run.py
```

No Linux:

```bash
cd /caminho/para/luna_gui
conda create -n luna-gui python=3.11 -y
conda activate luna-gui
pip install -r requirements.txt
python run.py
```

### 3. Instale O Ambiente LUNA Pela GUI

Ao abrir o aplicativo pela primeira vez:

1. Vá para `1. Início`.
2. Clique em `Verificar novamente`.
3. Clique em `Instalar LUNA`.
4. Aguarde a criação do ambiente `luna-env`.

Depois disso, a GUI passa a usar o `luna-env` para rodar os cálculos do LUNA.

## Uso Pelo Launcher

### Windows

```powershell
cd D:\projeto_LUNA_GUIDE\luna_gui\dist\windows
.\run_gui.bat
```

Se o launcher não encontrar o Python da GUI automaticamente:

```powershell
set LUNA_GUI_PYTHON=C:\caminho\para\python.exe
.\run_gui.bat
```

### Linux

```bash
cd /caminho/para/luna_gui/dist/linux
chmod +x run_gui.sh
./run_gui.sh
```

Se o ambiente da GUI tiver outro nome:

```bash
HIP2LINTERACTOMICS_GUI_ENV=meu-env ./run_gui.sh
```

Se quiser apontar diretamente para um Python:

```bash
HIP2LINTERACTOMICS_GUI_PYTHON=/caminho/para/python ./run_gui.sh
```

## Dependências Gráficas No Linux

Para PyQt6, PyMOL e OpenGL funcionarem corretamente no Linux:

```bash
# Debian / Ubuntu
sudo apt install libxcb-cursor0 libxkbcommon-x11-0 libgl1 libegl1

# Fedora / RHEL
sudo dnf install xcb-util-cursor libxkbcommon-x11 mesa-libGL mesa-libEGL
```

Em WSL2, use um servidor gráfico compatível ou o WSLg atualizado.

## Execução Por Terminal

O arquivo `hipplinteractomics_terminal.py` permite rodar o fluxo completo sem abrir a interface gráfica. Ele recebe um arquivo JSON ou um dicionário Python literal com os mesmos campos salvos pela GUI.

Exemplo:

```bash
conda activate luna-gui
python hipplinteractomics_terminal.py config_projeto.json
```

O arquivo de configuração pode conter, por exemplo:

```json
{
  "protein_file": "/caminho/para/receptor.pdb",
  "ligand_file": "/caminho/para/ligantes",
  "selected_ligands": "ALL",
  "trajectory_analysis": false,
  "workdir": "/caminho/para/workdir",
  "out_ifp": true,
  "ifp_type": "ALL",
  "ifp_levels": 6,
  "ifp_radius": 2.0,
  "ifp_length": 4096,
  "ifp_bit": true,
  "sim_matrix": true,
  "out_pse": true,
  "include_waters": true,
  "add_h": true,
  "ph": 7.4,
  "nproc": 1,
  "overwrite": true,
  "fp_labels_csv": "",
  "fp_labels_id_column": "",
  "fp_labels_column": "",
  "fp_label_task": "regression",
  "fp_use_otsu_threshold": true
}
```

Os outputs são gerados no mesmo formato usado pela GUI, então um projeto calculado por terminal pode ser aberto depois em `5. Resultados`.

## Análises FP

A aba `5. Resultados > Análises FP` integra os fingerprints brutos do LUNA com informações de shells, colisões, classes de interação, níveis e modelos de importância.

O fluxo geral é:

1. Carregar matriz `entries x FP`.
2. Carregar detalhes de shells e ocorrências de fingerprints.
3. Classificar a natureza dominante de cada feature.
4. Resolver colisões por prevalência, z-score e Otsu quando necessário.
5. Atribuir nível de shell predominante.
6. Desconsiderar features não confiáveis por classe ou nível.
7. Treinar modelos separados por nível assinado.
8. Calcular importância das features em cada modelo.
9. Transformar z-scores dos coeficientes de importância em p-values pela equação de Keiser and Hert [1].
10. Gerar tabelas, gráficos e mapas de calor usando labels no formato `FP_nível`.

Quando um projeto antigo é carregado e a matriz de fingerprints ainda não possui níveis assinados, a GUI tenta reconstruir os labels e regravar os CSVs compatíveis com a análise atual.

## Principais Resultados Gerados

- Tabelas de interações por ligante, frame ou pose.
- Mapas de calor por tipo de interação.
- Mapa de calor completo ligantes x resíduos.
- Estatísticas globais de interações.
- Gráficos por átomo do ligante em análises de trajetória/poses.
- Matrizes de fingerprints.
- Matrizes de similaridade.
- Clusters hierárquicos.
- Dashboards de importância de fingerprints.
- Sessões PyMOL filtradas.
- Sessões PyMOL para shells de fingerprints.
- Relatórios HTML e PDF.

## Solução De Problemas

### Conda Não É Encontrado

Verifique se `conda` está no PATH:

```bash
conda info
```

Se necessário, configure a variável usada pela GUI:

```bash
# Linux
export HIP2LINTERACTOMICS_GUI_CONDA=/caminho/para/conda

# Windows PowerShell
$env:HIP2LINTERACTOMICS_GUI_CONDA="C:\caminho\para\conda.exe"
```

### PyMOL Ou OpenGL Falha No Linux

Instale as dependências gráficas do sistema e confirme que a sessão X11/Wayland está funcionando. Em ambientes remotos, WSL ou containers, erros nativos de PyMOL/OpenGL geralmente indicam bibliotecas gráficas ausentes ou incompatíveis.

### Paralelização No Windows

No Windows, alguns fluxos do LUNA têm limitações com `multiprocessing` em modo `spawn`. Para bibliotecas muito grandes, prefira Linux ou WSL2 e aumente `nproc` na aba `4. Executar`.

### Projetos Antigos

Projetos calculados antes das atribuições de nível em fingerprints podem ser carregados normalmente. A GUI tenta atualizar as matrizes e recalcular os dashboards FP no formato atual.

## Testes

Para rodar os testes principais:

```bash
conda activate luna-gui
python -m pytest tests
```

Para uma verificação rápida dos módulos ligados a LUNA/PyMOL:

```bash
python -B -m pytest tests/test_pymol_launcher.py tests/test_luna_api_runner.py
```

## Nome Do Projeto

O nome público do software é `HIP²LInterActomics`. O nome técnico do repositório pode continuar como `LUNA_GUI` ou `HIP2LInterActomics_GUI`, e o pacote Python interno permanece `luna_gui` para manter compatibilidade com imports e projetos salvos.

## Referência

[1] Keiser, M. J.; Hert, J. In Chemogenomics: Methods and Applications; Jacoby, E., Ed.; Methods in Molecular Biology; Humana Press: Totowa, NJ, 2009; pp 195-205.
