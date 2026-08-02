# HIP²LInterActomics

Fluxos gráfico e headless para preparação molecular, execução do LUNA e análise de interações intermoleculares.

[Português](#português) · [Español](#español) · [English](#english)

## Downloads

- [Windows 10/11 - instalador `.exe`](../../releases/latest/download/HIP2LInterActomics-Setup.exe)
- [Linux x86_64 - aplicativo `.AppImage`](../../releases/latest/download/HIP2LInterActomics-x86_64.AppImage)

Os instaladores incluem a interface e suas dependências. Não é necessário instalar Python ou PyQt6. O instalador Windows cria automaticamente atalhos no menu Iniciar e na área de trabalho; a AppImage registra o aplicativo e cria o atalho na área de trabalho no primeiro lançamento. Na primeira análise, o próprio aplicativo orienta a instalação do ambiente científico LUNA.

---

## Português

### Visão geral

HIP²LInterActomics oferece uma interface PyQt6 e dois comandos de terminal independentes da interface:

- <code>hipplinteractomics-terminal</code>: executa um projeto LUNA a partir de JSON ou de argumentos diretos.
- <code>hipplinteractomics-multiple-run</code>: gera o produto cartesiano de configurações e executa cada projeto em série.
- <code>hip2linteractomics</code>: abre a interface gráfica.

O modo de terminal é realmente headless: não importa PyQt, não cria <code>QApplication</code> e não inicia loop gráfico. O backend Matplotlib é fixado em <code>Agg</code> e o Qt, quando alcançado indiretamente por alguma dependência, recebe <code>QT_QPA_PLATFORM=offscreen</code>.

A aplicação produz tabelas de interações, fingerprints HIFP/EIFP/FIFP, matrizes de Tanimoto, agrupamentos, mapas de calor, análise de importância de features, sessões PyMOL e relatórios HTML/PDF.

As melhorias científicas propostas, ainda não implementadas, estão organizadas em [Scientific Methodology Roadmap](docs/SCIENTIFIC_METHODOLOGY.md).

### Requisitos

- Windows 10/11 de 64 bits ou Linux x86_64.
- Python 3.9 ou superior; Python 3.11 é recomendado.
- Miniconda, Anaconda ou Miniforge para criar o ambiente científico <code>luna-env</code>.
- Acesso à internet na primeira instalação do LUNA.

São usados dois ambientes separados:

- ambiente do aplicativo: PyQt6, Matplotlib, NumPy, SciPy e scikit-learn;
- <code>luna-env</code>: LUNA, RDKit, Open Babel, PyMOL, Biopython e dependências químicas.

Os comandos headless podem ser usados em servidores sem X11, Wayland ou desktop. PyMOL interativo continua exigindo uma sessão gráfica; a geração normal do pipeline não exige a GUI do aplicativo.

### Opção Terminal-Only

O arquivo <code>environment.yml</code> prepara o aplicativo e a pilha científica no mesmo <code>luna-env</code>, sem a GUI PyQt6 do aplicativo. LUNA é instalado logo após a ativação, em uma segunda fase, porque seu metadata exige <code>pdbecif</code> já presente. PyMOL open source permanece no ambiente para gerar shells e sessões <code>.pse</code> diretamente com <code>pymol -cq</code>. O pacote Conda do PyMOL pode trazer Qt/PyQt5 como dependência transitiva própria; isso não inicializa a interface HIP²LInterActomics. <code>matplotlib-base</code> mantém os demais relatórios em <code>Agg</code>.

Conteúdo completo:

~~~yaml
name: luna-env

channels:
  - conda-forge
  - nodefaults

dependencies:
  # Runtime compatible with the current LUNA release.
  - python=3.9
  - pip>=23

  # Molecular stack; PyMOL is retained for headless PSE export.
  - biopython=1.79
  - rdkit
  - openbabel
  - pymol-open-source
  - networkx

  # Numerical analysis and static, headless reporting.
  - numpy<2
  - pandas
  - scipy
  - scikit-learn
  - matplotlib-base
  - seaborn-base

  # Phase 1: install pip prerequisites before LUNA metadata is evaluated.
  # Phase 2 is documented below: pip install --no-build-isolation -U luna.
  - pip:
      - pdbecif
      - mmh3<4
      - xopen
      - colorlog
~~~

Windows, no PowerShell ou Anaconda Prompt:

~~~powershell
conda env create -f environment.yml
conda activate luna-env
python -m pip install --no-build-isolation -U luna
python luna_gui\core\_luna_patch.py
python -m pip install --no-deps -e .
hipplinteractomics-terminal --help
hipplinteractomics-multiple-run --help
hipplinteractomics-terminal --config projeto.json
~~~

Linux:

~~~bash
conda env create -f environment.yml
conda activate luna-env
python -m pip install --no-build-isolation -U luna
python luna_gui/core/_luna_patch.py
python -m pip install --no-deps -e .
hipplinteractomics-terminal --help
hipplinteractomics-multiple-run --help
hipplinteractomics-terminal --config projeto.json
~~~

Execute os comandos a partir da raiz clonada do repositório; o comando <code>pip install --no-deps -e .</code> registra as CLIs sem reinstalar a pilha resolvida pelo Conda. Para atualizar um ambiente existente após alterar o YAML:

~~~bash
conda env update -n luna-env -f environment.yml --prune
~~~

### Instalação do pacote e comandos globais

O arquivo <code>pyproject.toml</code> registra automaticamente os três comandos. A instalação cria wrappers em <code>Scripts\</code> no Windows e <code>bin/</code> no Linux.

#### Windows — ambiente virtual

~~~powershell
py -3.11 -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip install ".[gui]"
hipplinteractomics-terminal --help
hipplinteractomics-multiple-run --help
hip2linteractomics
~~~

#### Linux — ambiente virtual

~~~bash
python3.11 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install ".[gui]"
hipplinteractomics-terminal --help
hipplinteractomics-multiple-run --help
hip2linteractomics
~~~

#### Instalação isolada e disponível no PATH com pipx

Windows:

~~~powershell
py -m pip install --user pipx
py -m pipx ensurepath
pipx install ".[gui]"
~~~

Linux:

~~~bash
python3 -m pip install --user pipx
python3 -m pipx ensurepath
pipx install ".[gui]"
~~~

Abra um novo terminal após <code>ensurepath</code>. O <code>pipx</code> expõe todos os comandos do pacote sem misturar dependências com o Python do sistema.

Para desenvolvimento editável:

~~~bash
python -m pip install -e .
~~~

Para construir e instalar a distribuição oficial:

~~~bash
python -m pip install build
python -m build
python -m pip install dist/hip2linteractomics-1.1.0-py3-none-any.whl
~~~

### Preparação do LUNA

A forma mais simples é abrir <code>hip2linteractomics</code>, acessar a aba inicial e escolher <strong>Instalar LUNA</strong>. O comando headless também aceita o Python de um ambiente já preparado:

~~~bash
hipplinteractomics-terminal --config projeto.json \
  --python-exe /caminho/para/luna-env/bin/python
~~~

No PowerShell:

~~~powershell
hipplinteractomics-terminal --config projeto.json --python-exe C:\caminho\para\luna-env\python.exe
~~~

Quando <code>--python-exe</code> não é fornecido, o programa procura o Conda e resolve o ambiente indicado por <code>--env-name</code>, cujo padrão é <code>luna-env</code>.

### Pré-processamento de complexos

A interface em <strong>Projeto > Preparar arquivos de complexos</strong> e o terminal usam o mesmo processador em lote. Uma pasta homogênea de MOL2 gera <code>proteinas_pdb/</code> e <code>ligantes_mol2/</code>; uma pasta de PDB/ENT gera <code>proteinas_pdb/</code> e <code>ligantes_sdf/</code>. Águas permanecem com a proteína.

~~~bash
hipplinteractomics-terminal \
  --prepare-complexes /dados/complexos \
  --prepare-output /dados/preparados \
  --python-exe /caminho/para/luna-env/bin/python
~~~

Para MOL2, <code>--last-protein-atom N</code> substitui a detecção automática. Para PDB/ENT, o conversor usa primeiro o Open Babel do <code>luna-env</code> e depois RDKit como fallback, com ambiente isolado e pools numéricos limitados a um thread para evitar picos de memória. O resumo e eventuais avisos ficam em <code>preprocess.log</code>.

### Comando headless individual

Gere um JSON completo:

~~~bash
hipplinteractomics-terminal --write-template projeto.json
~~~

Valide e mostre o comando LUNA sem executar:

~~~bash
hipplinteractomics-terminal --config projeto.json --dry-run
~~~

Execute:

~~~bash
hipplinteractomics-terminal --config projeto.json
~~~

Regenere somente artefatos de um projeto concluído:

~~~bash
hipplinteractomics-terminal --config projeto.json --results-only
~~~

Argumentos diretos podem substituir qualquer valor do JSON:

~~~bash
hipplinteractomics-terminal \
  --protein-file /dados/receptor.pdb \
  --ligand-file /dados/ligantes \
  --workdir /resultados/ensaio-01 \
  --ifp-type EIFP \
  --ifp-levels 2 \
  --ifp-radius 10 \
  --ifp-length 2048 \
  --ifp-format bin \
  --sim-matrix \
  --nproc 4
~~~

Um JSON mínimo tem esta forma:

~~~json
{
  "protein_file": "/dados/receptor.pdb",
  "ligand_file": "/dados/ligantes",
  "selected_ligands": "ALL",
  "workdir": "/resultados/base",
  "out_ifp": true,
  "ifp_type": "EIFP",
  "ifp_levels": 2,
  "ifp_radius": 10,
  "ifp_length": 2048,
  "ifp_bit": true,
  "sim_matrix": true,
  "nproc": 1,
  "overwrite": false,
  "python_exe": "/opt/conda/envs/luna-env/bin/python",
  "terminal_results": true
}
~~~

Também é aceita a forma aninhada com objetos <code>project</code> e <code>terminal</code>. Argumentos da CLI têm precedência sobre o arquivo. Use <code>hipplinteractomics-terminal --help</code> para a lista completa de entradas, saídas, filtros, parâmetros de fingerprint, runtime e exportação.

### Execução múltipla em série

O orquestrador executa duas etapas estritas:

1. cria todos os JSONs do produto cartesiano;
2. chama <code>hipplinteractomics-terminal</code> uma vez para cada JSON, aguardando o término antes de iniciar o seguinte.

Exemplo portável para PowerShell e Bash:

~~~bash
hipplinteractomics-multiple-run \
  --base-config projeto-base.json \
  --bits 1024 2048 \
  --formats bin cnt \
  --levels-growth 2:10 3:5 6:2 \
  --config-dir generated_configs
~~~

A sintaxe literal solicitada também é aceita no PowerShell:

~~~powershell
hipplinteractomics-multiple-run --base-config projeto-base.json --bits '[1024, 2048]' --formats '["bin", "cnt"]' --levels-growth '[(2,10), (3,5), (6,2)]' --config-dir generated_configs
~~~

As combinações acima geram 12 execuções. Cada configuração recebe:

- <code>ifp_length</code> a partir de <code>--bits</code>;
- <code>ifp_bit=true</code> para <code>bin</code> e <code>false</code> para <code>cnt</code>;
- <code>ifp_levels</code> e <code>ifp_radius</code> a partir de cada par <code>(nível, growth_ratio)</code>;
- um <code>workdir</code> exclusivo, evitando colisão entre resultados.

Os JSONs são gravados atomicamente com temporários únicos antes do primeiro subprocesso. A saída combinada de stdout/stderr é mostrada no terminal e salva com buffer de 64 KiB em <code>generated_configs/logs/</code>, sem <code>flush()</code> por linha. O estado fica em <code>generated_configs/pipeline_summary.json</code>: ele é criado antes da execução e atualizado atomicamente após cada run. Por padrão o lote para no primeiro exit code diferente de zero; <code>--continue-on-error</code> processa os itens restantes.

Ao repetir o mesmo comando, o orquestrador compara o SHA-256 de cada configuração com o resumo e pula somente runs marcados como <code>completed</code> cujo conteúdo continua idêntico. Runs <code>failed</code>, <code>interrupted</code> ou modificados são executados novamente. Assim, uma queda após a 30ª de 50 simulações retoma na 31ª sem uma opção adicional:

~~~bash
# Primeira execução ou retomada: use exatamente o mesmo comando.
hipplinteractomics-multiple-run \
  --base-config projeto-base.json \
  --bits 1024 2048 \
  --formats bin cnt \
  --levels-growth 2:10 3:5 6:2 \
  --config-dir generated_configs
~~~

<code>SIGTERM</code> e <code>SIGINT</code> são encaminhados ao subprocesso ativo; o run é marcado como <code>interrupted</code>, o resumo é persistido e o processo retorna <code>128 + sinal</code>. Em Slurm, <code>nproc</code> deriva de <code>SLURM_CPUS_PER_TASK</code> ou <code>SLURM_JOB_CPUS_PER_NODE</code>. Em PBS/Torque, usa <code>PBS_NUM_PPN</code> ou os slots locais de <code>PBS_NODEFILE</code>. Fora de um scheduler, usa o pedido limitado por <code>os.cpu_count()</code>; no Windows nativo, a proteção do LUNA mantém <code>nproc=1</code>. A lista estática de ligantes é lida uma vez pelo orquestrador e materializada nos JSONs filhos.

Para chamar um arquivo ou executável específico:

~~~bash
hipplinteractomics-multiple-run \
  --base-config projeto-base.json \
  --bits 2048 \
  --formats bin \
  --levels-growth 2:10 \
  --terminal-executable ./hipplinteractomics_terminal.py
~~~

### Execução direta pelo código-fonte

A instalação do pacote é recomendada. Durante desenvolvimento, os equivalentes são:

~~~bash
python hipplinteractomics_terminal.py --config projeto.json
python hipplinteractomics_multiple_run.py \
  --base-config projeto-base.json \
  --bits 1024 2048 \
  --formats bin cnt \
  --levels-growth 2:10 3:5
python run.py
~~~

### Pacotes desktop

Windows:

~~~powershell
.\installer\build_windows.ps1 -InstallBuildDependencies
~~~

Linux:

~~~bash
bash installer/build_linux.sh
~~~

O PyInstaller não faz compilação cruzada. O Windows produz <code>HIP2LInterActomics-Setup.exe</code> com Inno Setup 6 e cria o atalho da área de trabalho durante a instalação. O Linux produz <code>HIP2LInterActomics-x86_64.AppImage</code> e um bundle <code>.tar.gz</code>; no primeiro lançamento, a AppImage registra o menu de aplicativos e cria o atalho usando seu caminho atual. Em uma tag de versão, o workflow <code>.github/workflows/build-installers.yml</code> publica automaticamente o instalador Windows e a AppImage na release do GitHub.

#### Inclusão física do YAML nas distribuições

A distribuição Python usa:

~~~text
# MANIFEST.in
include environment.yml
~~~

~~~toml
[tool.setuptools.data-files]
"share/hip2linteractomics" = ["environment.yml"]
~~~

O bundle nativo usa a entrada <code>("environment.yml", ".")</code> em <code>HIP2LInterActomics.spec</code>. No Windows, a regra recursiva de <code>[Files]</code> do Inno Setup copia o YAML para <code>{app}\environment.yml</code>. No Linux, <code>installer/build_linux.sh</code> valida o arquivo no bundle antes de criar o <code>.tar.gz</code>. Os dois scripts de build falham explicitamente se o asset estiver ausente.


### Testes

~~~bash
python -m pytest -q
python -m compileall -q luna_gui tests tools \
  hipplinteractomics_terminal.py hipplinteractomics_multiple_run.py
~~~

### Estrutura relevante

~~~text
hipplinteractomics_terminal.py       CLI headless individual
hipplinteractomics_multiple_run.py   matriz JSON e execução serial
pyproject.toml                       metadados e comandos instaláveis
luna_gui/                            pacote e interface PyQt6
tests/                               testes automatizados
installer/                           builds PyInstaller/Inno Setup
.github/workflows/                   CI multiplataforma
~~~

Dados experimentais, ambientes Conda, workdirs, logs e diretórios de build não devem ser versionados.

### Solução de problemas

- <strong>Comando não encontrado:</strong> ative o ambiente virtual ou execute novamente <code>pipx ensurepath</code> e abra outro terminal.
- <strong>Conda não encontrado:</strong> use <code>--conda-exe</code>, <code>--python-exe</code> ou defina <code>HIP2LINTERACTOMICS_GUI_CONDA</code>.
- <strong>LUNA não instalado:</strong> crie o <code>luna-env</code> pela interface ou informe um Python que consiga executar <code>import luna</code>.
- <strong>Falha multinúcleo no Windows:</strong> use <code>--nproc 1</code>; para paralelismo real, use Linux ou WSL2.
- <strong>Lote interrompido:</strong> consulte o arquivo <code>.log</code> da combinação e <code>pipeline_summary.json</code>.
- <strong>PyMOL/OpenGL:</strong> a visualização interativa exige sessão gráfica, mesmo quando o cálculo é headless.

### Licença e referência

Licença MIT; consulte [LICENSE](LICENSE).

Keiser, M. J.; Hert, J. In <em>Chemogenomics: Methods and Applications</em>; Jacoby, E., Ed.; Methods in Molecular Biology; Humana Press: Totowa, NJ, 2009; pp. 195–205.

---


## Español

### Descripción general

HIP²LInterActomics ofrece una aplicación PyQt6 y dos comandos de terminal independientes de la interfaz:

- <code>hipplinteractomics-terminal</code> ejecuta un proyecto LUNA desde JSON y/o argumentos directos.
- <code>hipplinteractomics-multiple-run</code> genera una matriz cartesiana de configuraciones y la procesa en serie.
- <code>hip2linteractomics</code> abre la interfaz gráfica.

El flujo de terminal es completamente headless: no importa PyQt, no crea <code>QApplication</code> ni inicia un loop gráfico. Matplotlib utiliza <code>Agg</code> y cualquier acceso indirecto a Qt recibe <code>QT_QPA_PLATFORM=offscreen</code>.

La aplicación puede producir tablas de interacciones, fingerprints HIFP/EIFP/FIFP, matrices de Tanimoto, clusters, mapas de calor, análisis supervisados de importancia, sesiones PyMOL e informes HTML/PDF.

Las mejoras científicas propuestas, todavía no implementadas, están organizadas en [Scientific Methodology Roadmap](docs/SCIENTIFIC_METHODOLOGY.md).

### Requisitos

- Windows 10/11 de 64 bits o Linux x86_64.
- Python 3.9 o superior; se recomienda Python 3.11.
- Miniconda, Anaconda o Miniforge para el entorno científico <code>luna-env</code>.
- Acceso a internet durante la primera instalación de LUNA.

El entorno de la aplicación contiene PyQt6 y la pila de análisis/visualización. El entorno separado <code>luna-env</code> contiene LUNA, RDKit, Open Babel, PyMOL, Biopython y las dependencias químicas.

Los comandos headless funcionan en servidores sin X11, Wayland o escritorio. La visualización interactiva con PyMOL sí requiere una sesión gráfica.

### Opción Terminal-Only

El archivo <code>environment.yml</code> prepara la aplicación y la pila científica en el mismo <code>luna-env</code>, sin la GUI PyQt6 de la aplicación. LUNA se instala después de activar el entorno, en una segunda fase, porque su metadata requiere que <code>pdbecif</code> ya exista. PyMOL open source permanece disponible para generar shells y sesiones <code>.pse</code> con <code>pymol -cq</code>. El paquete Conda de PyMOL puede incorporar Qt/PyQt5 como dependencia transitiva propia, pero no inicializa la interfaz HIP²LInterActomics. <code>matplotlib-base</code> mantiene los demás informes en <code>Agg</code>.

~~~yaml
name: luna-env

channels:
  - conda-forge
  - nodefaults

dependencies:
  # Runtime compatible with the current LUNA release.
  - python=3.9
  - pip>=23

  # Molecular stack; PyMOL is retained for headless PSE export.
  - biopython=1.79
  - rdkit
  - openbabel
  - pymol-open-source
  - networkx

  # Numerical analysis and static, headless reporting.
  - numpy<2
  - pandas
  - scipy
  - scikit-learn
  - matplotlib-base
  - seaborn-base

  # Phase 1: install pip prerequisites before LUNA metadata is evaluated.
  # Phase 2 is documented below: pip install --no-build-isolation -U luna.
  - pip:
      - pdbecif
      - mmh3<4
      - xopen
      - colorlog
~~~

Windows, desde PowerShell o Anaconda Prompt:

~~~powershell
conda env create -f environment.yml
conda activate luna-env
python -m pip install --no-build-isolation -U luna
python luna_gui\core\_luna_patch.py
python -m pip install --no-deps -e .
hipplinteractomics-terminal --help
hipplinteractomics-multiple-run --help
hipplinteractomics-terminal --config proyecto.json
~~~

Linux:

~~~bash
conda env create -f environment.yml
conda activate luna-env
python -m pip install --no-build-isolation -U luna
python luna_gui/core/_luna_patch.py
python -m pip install --no-deps -e .
hipplinteractomics-terminal --help
hipplinteractomics-multiple-run --help
hipplinteractomics-terminal --config proyecto.json
~~~

Ejecute estos pasos desde la raíz clonada; <code>pip install --no-deps -e .</code> registra las CLIs sin reinstalar la pila resuelta por Conda. Para actualizar el entorno:

~~~bash
conda env update -n luna-env -f environment.yml --prune
~~~


### Instalación del paquete y comandos globales

<code>pyproject.toml</code> registra los tres comandos. La instalación crea wrappers nativos en <code>Scripts\</code> en Windows y <code>bin/</code> en Linux.

#### Windows — entorno virtual

~~~powershell
py -3.11 -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip install ".[gui]"
hipplinteractomics-terminal --help
hipplinteractomics-multiple-run --help
hip2linteractomics
~~~

#### Linux — entorno virtual

~~~bash
python3.11 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install ".[gui]"
hipplinteractomics-terminal --help
hipplinteractomics-multiple-run --help
hip2linteractomics
~~~

#### Instalación aislada en PATH con pipx

Windows:

~~~powershell
py -m pip install --user pipx
py -m pipx ensurepath
pipx install ".[gui]"
~~~

Linux:

~~~bash
python3 -m pip install --user pipx
python3 -m pipx ensurepath
pipx install ".[gui]"
~~~

Abra otra terminal después de <code>ensurepath</code>. Para desarrollo editable utilice <code>python -m pip install -e .</code>.

Para construir e instalar el wheel oficial:

~~~bash
python -m pip install build
python -m build
python -m pip install dist/hip2linteractomics-1.1.0-py3-none-any.whl
~~~

### Preparación de LUNA

La opción más sencilla es ejecutar <code>hip2linteractomics</code> y seleccionar <strong>Instalar LUNA</strong> en la primera pestaña. También puede indicar un entorno ya preparado:

~~~bash
hipplinteractomics-terminal --config proyecto.json \
  --python-exe /ruta/a/luna-env/bin/python
~~~

PowerShell:

~~~powershell
hipplinteractomics-terminal --config proyecto.json --python-exe C:\ruta\a\luna-env\python.exe
~~~

Sin <code>--python-exe</code>, el comando resuelve Conda y el entorno indicado por <code>--env-name</code>, cuyo valor predeterminado es <code>luna-env</code>.

### Preprocesamiento de complejos

La interfaz en <strong>Proyecto > Preparar archivos de complejos</strong> y el terminal utilizan el mismo procesador por lotes. Una carpeta homogénea de MOL2 genera <code>proteinas_pdb/</code> y <code>ligantes_mol2/</code>; una carpeta de PDB/ENT genera <code>proteinas_pdb/</code> y <code>ligantes_sdf/</code>. Las aguas permanecen con la proteína.

~~~bash
hipplinteractomics-terminal \
  --prepare-complexes /datos/complejos \
  --prepare-output /datos/preparados \
  --python-exe /ruta/a/luna-env/bin/python
~~~

Para MOL2, <code>--last-protein-atom N</code> sustituye la detección automática. Para PDB/ENT, el conversor prueba primero Open Babel desde <code>luna-env</code> y después RDKit, con un entorno aislado y un solo thread numérico para evitar picos de memoria. El resumen queda en <code>preprocess.log</code>.

### Ejecución headless individual

~~~bash
hipplinteractomics-terminal --write-template proyecto.json
hipplinteractomics-terminal --config proyecto.json --dry-run
hipplinteractomics-terminal --config proyecto.json
hipplinteractomics-terminal --config proyecto.json --results-only
~~~

Los argumentos directos sustituyen los valores del JSON:

~~~bash
hipplinteractomics-terminal \
  --protein-file /datos/receptor.pdb \
  --ligand-file /datos/ligandos \
  --workdir /resultados/ejecucion-01 \
  --ifp-type EIFP \
  --ifp-levels 2 \
  --ifp-radius 10 \
  --ifp-length 2048 \
  --ifp-format bin \
  --sim-matrix \
  --nproc 4
~~~

JSON mínimo:

~~~json
{
  "protein_file": "/datos/receptor.pdb",
  "ligand_file": "/datos/ligandos",
  "selected_ligands": "ALL",
  "workdir": "/resultados/base",
  "out_ifp": true,
  "ifp_type": "EIFP",
  "ifp_levels": 2,
  "ifp_radius": 10,
  "ifp_length": 2048,
  "ifp_bit": true,
  "sim_matrix": true,
  "nproc": 1,
  "overwrite": false,
  "python_exe": "/opt/conda/envs/luna-env/bin/python",
  "terminal_results": true
}
~~~

También se aceptan objetos anidados <code>project</code> y <code>terminal</code>. Ejecute <code>hipplinteractomics-terminal --help</code> para consultar todas las entradas, salidas, opciones de fingerprint, filtros, runtime y exportación.

### Ejecución múltiple en serie

El orquestador siempre realiza dos etapas ordenadas: primero escribe todos los JSON del producto cartesiano y después inicia un subprocess síncrono por configuración.

~~~bash
hipplinteractomics-multiple-run \
  --base-config proyecto-base.json \
  --bits 1024 2048 \
  --formats bin cnt \
  --levels-growth 2:10 3:5 6:2 \
  --config-dir generated_configs
~~~

También se admite la sintaxis literal solicitada:

~~~powershell
hipplinteractomics-multiple-run --base-config proyecto-base.json --bits '[1024, 2048]' --formats '["bin", "cnt"]' --levels-growth '[(2,10), (3,5), (6,2)]' --config-dir generated_configs
~~~

Este ejemplo genera 12 ejecuciones. Cada una recibe un directorio exclusivo y asigna bits a <code>ifp_length</code>, formato a <code>ifp_bit</code>, nivel a <code>ifp_levels</code> y growth ratio a <code>ifp_radius</code>.

Todos los JSON se escriben atómicamente con nombres temporales únicos antes de iniciar la ejecución. stdout/stderr se muestran y se guardan con un búfer de 64 KiB en <code>generated_configs/logs/</code>, sin <code>flush()</code> por línea. <code>generated_configs/pipeline_summary.json</code> se crea antes del primer run y se actualiza atómicamente después de cada resultado. El lote se detiene con el primer exit code distinto de cero, salvo que se utilice <code>--continue-on-error</code>.

Al repetir exactamente el mismo comando, el orquestador compara el SHA-256 de cada configuración y omite solo los runs <code>completed</code> cuyo contenido no cambió. Los estados <code>failed</code> e <code>interrupted</code>, así como cualquier JSON modificado, se ejecutan nuevamente. Por ejemplo, una caída después de la simulación 30 de 50 reanuda automáticamente en la 31:

~~~bash
hipplinteractomics-multiple-run \
  --base-config proyecto-base.json \
  --bits 1024 2048 \
  --formats bin cnt \
  --levels-growth 2:10 3:5 6:2 \
  --config-dir generated_configs
~~~

<code>SIGTERM</code> y <code>SIGINT</code> se transmiten al subprocess activo; el run queda como <code>interrupted</code>, el resumen se persiste y el proceso devuelve <code>128 + señal</code>. En Slurm, <code>nproc</code> se obtiene de <code>SLURM_CPUS_PER_TASK</code> o <code>SLURM_JOB_CPUS_PER_NODE</code>. En PBS/Torque utiliza <code>PBS_NUM_PPN</code> o los slots locales de <code>PBS_NODEFILE</code>. Sin scheduler, la solicitud queda limitada por <code>os.cpu_count()</code>; en Windows nativo, la protección de LUNA conserva <code>nproc=1</code>. La lista estática de ligandos se lee una sola vez y se materializa en los JSON hijos.

Utilice <code>--terminal-executable PATH</code> para seleccionar un script o ejecutable específico.

### Ejecución desde el código fuente

~~~bash
python hipplinteractomics_terminal.py --config proyecto.json
python hipplinteractomics_multiple_run.py \
  --base-config proyecto-base.json \
  --bits 1024 2048 \
  --formats bin cnt \
  --levels-growth 2:10 3:5
python run.py
~~~

### Builds de escritorio

Windows:

~~~powershell
.\installer\build_windows.ps1 -InstallBuildDependencies
~~~

Linux:

~~~bash
bash installer/build_linux.sh
~~~

PyInstaller no realiza compilación cruzada. Windows produce <code>HIP2LInterActomics-Setup.exe</code> con Inno Setup 6 y crea el acceso directo del escritorio durante la instalación. Linux produce <code>HIP2LInterActomics-x86_64.AppImage</code> y un bundle <code>.tar.gz</code>; en el primer inicio, la AppImage se registra en el menú y crea el acceso directo con su ruta actual. Para una etiqueta de versión, GitHub Actions publica automáticamente el instalador de Windows y la AppImage en la release.

#### Inclusión física del YAML

<code>MANIFEST.in</code> contiene <code>include environment.yml</code> y <code>pyproject.toml</code> instala el asset en <code>share/hip2linteractomics</code>. PyInstaller usa <code>("environment.yml", ".")</code>. La regla recursiva de Inno Setup lo copia a <code>{app}\environment.yml</code> y el script Linux valida el archivo antes de crear el <code>.tar.gz</code>. Ambos builds fallan si falta el asset.


### Pruebas

~~~bash
python -m pytest -q
python -m compileall -q luna_gui tests tools \
  hipplinteractomics_terminal.py hipplinteractomics_multiple_run.py
~~~

### Estructura relevante

~~~text
hipplinteractomics_terminal.py       CLI headless individual
hipplinteractomics_multiple_run.py   matriz JSON y ejecución serial
pyproject.toml                       metadatos y comandos instalables
luna_gui/                            paquete e interfaz PyQt6
tests/                               pruebas automatizadas
installer/                           builds PyInstaller/Inno Setup
.github/workflows/                   CI multiplataforma
~~~

No deben versionarse datos experimentales, entornos Conda, workdirs, logs ni intermediarios de build.

### Solución de problemas

- <strong>Comando no encontrado:</strong> active el entorno virtual o ejecute <code>pipx ensurepath</code> y abra otra terminal.
- <strong>Conda no encontrado:</strong> use <code>--conda-exe</code>, <code>--python-exe</code> o <code>HIP2LINTERACTOMICS_GUI_CONDA</code>.
- <strong>LUNA no instalado:</strong> cree <code>luna-env</code> desde la GUI o seleccione un Python capaz de ejecutar <code>import luna</code>.
- <strong>Error multiproceso en Windows:</strong> use <code>--nproc 1</code>; utilice Linux o WSL2 para paralelismo real.
- <strong>Lote interrumpido:</strong> revise el log de la combinación y <code>pipeline_summary.json</code>.
- <strong>PyMOL/OpenGL:</strong> la visualización interactiva necesita una sesión gráfica aunque el cálculo sea headless.

### Licencia y referencia

Licencia MIT; consulte [LICENSE](LICENSE).

Keiser, M. J.; Hert, J. In <em>Chemogenomics: Methods and Applications</em>; Jacoby, E., Ed.; Methods in Molecular Biology; Humana Press: Totowa, NJ, 2009; pp. 195–205.



## English

### Overview

HIP²LInterActomics provides a PyQt6 desktop application and two UI-independent terminal commands:

- <code>hipplinteractomics-terminal</code> runs one LUNA project from JSON and/or direct arguments.
- <code>hipplinteractomics-multiple-run</code> generates a Cartesian configuration matrix and processes it serially.
- <code>hip2linteractomics</code> opens the desktop interface.

The terminal path is fully headless. It does not import PyQt, create a <code>QApplication</code>, or start an event loop. Matplotlib uses <code>Agg</code>, and indirect Qt access is constrained with <code>QT_QPA_PLATFORM=offscreen</code>.

The application can produce interaction tables, HIFP/EIFP/FIFP fingerprints, Tanimoto matrices, clusters, heatmaps, supervised feature-importance analyses, PyMOL sessions, and HTML/PDF reports.

Proposed scientific improvements that are not yet implemented are organized in the [Scientific Methodology Roadmap](docs/SCIENTIFIC_METHODOLOGY.md).

### Requirements

- 64-bit Windows 10/11 or x86_64 Linux.
- Python 3.9 or newer; Python 3.11 is recommended.
- Miniconda, Anaconda, or Miniforge for the scientific <code>luna-env</code>.
- Internet access during the initial LUNA installation.

The application environment contains PyQt6 and the analysis/visualization stack. The separate <code>luna-env</code> contains LUNA, RDKit, Open Babel, PyMOL, Biopython, and chemistry dependencies.

Headless commands work on servers without X11, Wayland, or a desktop. Interactive PyMOL visualization still requires a graphical session.

### Terminal-Only option

<code>environment.yml</code> prepares the application and scientific stack in one <code>luna-env</code>, without the application's PyQt6 GUI. LUNA is installed immediately after activation in a second phase because its metadata requires <code>pdbecif</code> to exist first. Open-source PyMOL remains available to generate shells and <code>.pse</code> sessions through <code>pymol -cq</code>. PyMOL's Conda package may bring Qt/PyQt5 as its own transitive dependency, but it does not initialize the HIP²LInterActomics interface. <code>matplotlib-base</code> keeps all other reports on <code>Agg</code>.

~~~yaml
name: luna-env

channels:
  - conda-forge
  - nodefaults

dependencies:
  # Runtime compatible with the current LUNA release.
  - python=3.9
  - pip>=23

  # Molecular stack; PyMOL is retained for headless PSE export.
  - biopython=1.79
  - rdkit
  - openbabel
  - pymol-open-source
  - networkx

  # Numerical analysis and static, headless reporting.
  - numpy<2
  - pandas
  - scipy
  - scikit-learn
  - matplotlib-base
  - seaborn-base

  # Phase 1: install pip prerequisites before LUNA metadata is evaluated.
  # Phase 2 is documented below: pip install --no-build-isolation -U luna.
  - pip:
      - pdbecif
      - mmh3<4
      - xopen
      - colorlog
~~~

Windows PowerShell or Anaconda Prompt:

~~~powershell
conda env create -f environment.yml
conda activate luna-env
python -m pip install --no-build-isolation -U luna
python luna_gui\core\_luna_patch.py
python -m pip install --no-deps -e .
hipplinteractomics-terminal --help
hipplinteractomics-multiple-run --help
hipplinteractomics-terminal --config project.json
~~~

Linux:

~~~bash
conda env create -f environment.yml
conda activate luna-env
python -m pip install --no-build-isolation -U luna
python luna_gui/core/_luna_patch.py
python -m pip install --no-deps -e .
hipplinteractomics-terminal --help
hipplinteractomics-multiple-run --help
hipplinteractomics-terminal --config project.json
~~~

Run these steps from the cloned repository root; <code>pip install --no-deps -e .</code> registers the CLIs without reinstalling the stack resolved by Conda. Update an existing environment with:

~~~bash
conda env update -n luna-env -f environment.yml --prune
~~~


### Package installation and global commands

<code>pyproject.toml</code> registers all three commands. Installation creates native wrappers under <code>Scripts\</code> on Windows and <code>bin/</code> on Linux.

#### Windows virtual environment

~~~powershell
py -3.11 -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip install ".[gui]"
hipplinteractomics-terminal --help
hipplinteractomics-multiple-run --help
hip2linteractomics
~~~

#### Linux virtual environment

~~~bash
python3.11 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install ".[gui]"
hipplinteractomics-terminal --help
hipplinteractomics-multiple-run --help
hip2linteractomics
~~~

#### Isolated PATH installation with pipx

Windows:

~~~powershell
py -m pip install --user pipx
py -m pipx ensurepath
pipx install ".[gui]"
~~~

Linux:

~~~bash
python3 -m pip install --user pipx
python3 -m pipx ensurepath
pipx install ".[gui]"
~~~

Open a new terminal after <code>ensurepath</code>. For editable development use <code>python -m pip install -e .</code>.

Build and install the official wheel with:

~~~bash
python -m pip install build
python -m build
python -m pip install dist/hip2linteractomics-1.1.0-py3-none-any.whl
~~~

### Preparing LUNA

The simplest setup is to run <code>hip2linteractomics</code> and select <strong>Install LUNA</strong> on the first tab. A prepared environment can be selected directly:

~~~bash
hipplinteractomics-terminal --config project.json \
  --python-exe /path/to/luna-env/bin/python
~~~

PowerShell:

~~~powershell
hipplinteractomics-terminal --config project.json --python-exe C:\path\to\luna-env\python.exe
~~~

Without <code>--python-exe</code>, the command resolves Conda and the environment selected by <code>--env-name</code>, which defaults to <code>luna-env</code>.

### Complex preprocessing

The <strong>Project > Prepare complex files</strong> dialog and the terminal use the same batch processor. A homogeneous MOL2 folder produces <code>proteinas_pdb/</code> and <code>ligantes_mol2/</code>; a PDB/ENT folder produces <code>proteinas_pdb/</code> and <code>ligantes_sdf/</code>. Waters remain with the protein.

~~~bash
hipplinteractomics-terminal \
  --prepare-complexes /data/complexes \
  --prepare-output /data/prepared \
  --python-exe /path/to/luna-env/bin/python
~~~

For MOL2, <code>--last-protein-atom N</code> overrides automatic detection. For PDB/ENT, conversion tries Open Babel from <code>luna-env</code> first and RDKit second, with an isolated environment and one numerical thread to avoid memory spikes. The summary is written to <code>preprocess.log</code>.

### Single headless run

~~~bash
hipplinteractomics-terminal --write-template project.json
hipplinteractomics-terminal --config project.json --dry-run
hipplinteractomics-terminal --config project.json
hipplinteractomics-terminal --config project.json --results-only
~~~

Direct options override values loaded from JSON:

~~~bash
hipplinteractomics-terminal \
  --protein-file /data/receptor.pdb \
  --ligand-file /data/ligands \
  --workdir /results/run-01 \
  --ifp-type EIFP \
  --ifp-levels 2 \
  --ifp-radius 10 \
  --ifp-length 2048 \
  --ifp-format bin \
  --sim-matrix \
  --nproc 4
~~~

Minimal JSON:

~~~json
{
  "protein_file": "/data/receptor.pdb",
  "ligand_file": "/data/ligands",
  "selected_ligands": "ALL",
  "workdir": "/results/base",
  "out_ifp": true,
  "ifp_type": "EIFP",
  "ifp_levels": 2,
  "ifp_radius": 10,
  "ifp_length": 2048,
  "ifp_bit": true,
  "sim_matrix": true,
  "nproc": 1,
  "overwrite": false,
  "python_exe": "/opt/conda/envs/luna-env/bin/python",
  "terminal_results": true
}
~~~

Nested <code>project</code> and <code>terminal</code> objects are also accepted. Run <code>hipplinteractomics-terminal --help</code> for every input, output, fingerprint, filter, runtime, and export option.

### Serial batch execution

The orchestrator always completes two ordered stages: it first writes every Cartesian-product JSON and then launches one synchronous subprocess per configuration.

~~~bash
hipplinteractomics-multiple-run \
  --base-config base-project.json \
  --bits 1024 2048 \
  --formats bin cnt \
  --levels-growth 2:10 3:5 6:2 \
  --config-dir generated_configs
~~~

The requested literal syntax is also supported:

~~~powershell
hipplinteractomics-multiple-run --base-config base-project.json --bits '[1024, 2048]' --formats '["bin", "cnt"]' --levels-growth '[(2,10), (3,5), (6,2)]' --config-dir generated_configs
~~~

This example generates 12 runs. Each run receives a unique work directory and maps bits to <code>ifp_length</code>, format to <code>ifp_bit</code>, level to <code>ifp_levels</code>, and growth ratio to <code>ifp_radius</code>.

All JSON files are written atomically through unique temporary names before execution starts. Combined stdout/stderr is streamed and saved with a 64 KiB buffer under <code>generated_configs/logs/</code>, without line-by-line <code>flush()</code>. <code>generated_configs/pipeline_summary.json</code> is created before the first run and atomically updated after every result. Processing stops on the first nonzero exit code unless <code>--continue-on-error</code> is supplied.

Repeating the same command compares each configuration's SHA-256 digest with the summary and skips only unchanged runs marked <code>completed</code>. Runs marked <code>failed</code> or <code>interrupted</code>, and configurations whose contents changed, execute again. If a 50-run batch stops after run 30, the same command resumes with run 31:

~~~bash
hipplinteractomics-multiple-run \
  --base-config base-project.json \
  --bits 1024 2048 \
  --formats bin cnt \
  --levels-growth 2:10 3:5 6:2 \
  --config-dir generated_configs
~~~

<code>SIGTERM</code> and <code>SIGINT</code> are relayed to the active child; the run is recorded as <code>interrupted</code>, the summary is persisted, and the process exits with <code>128 + signal</code>. Under Slurm, <code>nproc</code> is derived from <code>SLURM_CPUS_PER_TASK</code> or <code>SLURM_JOB_CPUS_PER_NODE</code>. PBS/Torque uses <code>PBS_NUM_PPN</code> or the local slots in <code>PBS_NODEFILE</code>. Without a scheduler, the requested value is capped by <code>os.cpu_count()</code>; native Windows retains LUNA's safety limit of <code>nproc=1</code>. The static ligand list is read once by the orchestrator and materialized into every child JSON.

Use a specific source script or packaged executable with <code>--terminal-executable PATH</code>.

### Running from the source tree

~~~bash
python hipplinteractomics_terminal.py --config project.json
python hipplinteractomics_multiple_run.py \
  --base-config base-project.json \
  --bits 1024 2048 \
  --formats bin cnt \
  --levels-growth 2:10 3:5
python run.py
~~~

### Desktop builds

Windows:

~~~powershell
.\installer\build_windows.ps1 -InstallBuildDependencies
~~~

Linux:

~~~bash
bash installer/build_linux.sh
~~~

PyInstaller does not cross-compile. Windows produces <code>HIP2LInterActomics-Setup.exe</code> with Inno Setup 6 and creates the desktop shortcut during installation. Linux produces <code>HIP2LInterActomics-x86_64.AppImage</code> and a native <code>.tar.gz</code> bundle; on first launch, the AppImage registers itself in the application menu and creates a shortcut using its current path. For a version tag, GitHub Actions automatically publishes the Windows installer and AppImage to the GitHub release.

#### Physical YAML inclusion

<code>MANIFEST.in</code> contains <code>include environment.yml</code>, while <code>pyproject.toml</code> installs the asset under <code>share/hip2linteractomics</code>. PyInstaller uses <code>("environment.yml", ".")</code>. Inno Setup's recursive rule copies it to <code>{app}\environment.yml</code>, and the Linux build validates the file before creating the <code>.tar.gz</code>. Both builds fail when the asset is absent.


### Tests

~~~bash
python -m pytest -q
python -m compileall -q luna_gui tests tools \
  hipplinteractomics_terminal.py hipplinteractomics_multiple_run.py
~~~

### Relevant repository structure

~~~text
hipplinteractomics_terminal.py       single-run headless CLI
hipplinteractomics_multiple_run.py   JSON matrix and serial runner
pyproject.toml                       package metadata and installed commands
luna_gui/                            PyQt6 package and interface
tests/                               automated tests
installer/                           PyInstaller/Inno Setup builds
.github/workflows/                   cross-platform CI
~~~

Do not commit experimental data, Conda environments, work directories, logs, or build intermediates.

### Troubleshooting

- <strong>Command not found:</strong> activate the virtual environment or run <code>pipx ensurepath</code> and open a new terminal.
- <strong>Conda not found:</strong> use <code>--conda-exe</code>, <code>--python-exe</code>, or <code>HIP2LINTERACTOMICS_GUI_CONDA</code>.
- <strong>LUNA not installed:</strong> create <code>luna-env</code> through the GUI or select a Python that can execute <code>import luna</code>.
- <strong>Windows multiprocessing failure:</strong> use <code>--nproc 1</code>; use Linux or WSL2 for actual parallelism.
- <strong>Interrupted batch:</strong> inspect the combination log and <code>pipeline_summary.json</code>.
- <strong>PyMOL/OpenGL:</strong> interactive visualization requires a graphical session even when calculations are headless.

### License and reference

MIT License; see [LICENSE](LICENSE).

Keiser, M. J.; Hert, J. In <em>Chemogenomics: Methods and Applications</em>; Jacoby, E., Ed.; Methods in Molecular Biology; Humana Press: Totowa, NJ, 2009; pp. 195–205.

---
