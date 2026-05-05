# Manual de Instalação e Uso — HIP²LInterActomics

> Versão da interface: 1.0  
> Sistema operacional: Windows 10/11 ou Linux (Ubuntu 20.04+, Fedora 36+)  
> Nível necessário: **iniciante** — nenhum conhecimento de linha de comando é exigido (salvo a etapa de instalação inicial).

---

## Índice

1. [O que é o LUNA e o que esta GUI faz](#1-o-que-é-o-luna-e-o-que-esta-gui-faz)
2. [O que você NÃO precisa fazer](#2-o-que-você-não-precisa-fazer)
3. [Requisitos](#3-requisitos)
4. [Instalação — passo a passo](#4-instalação--passo-a-passo)
   - 4.1 Windows
   - 4.2 Linux
5. [Primeira execução — aba Setup](#5-primeira-execução--aba-setup)
6. [Fluxo básico de uso](#6-fluxo-básico-de-uso)
7. [Abas em detalhe](#7-abas-em-detalhe)
   - 7.1 Setup
   - 7.2 Projeto
   - 7.3 Análises
   - 7.4 Executar
   - 7.5 Resultados
   - 7.6 Histórico
8. [Preparar arquivos de docking (wizard)](#8-preparar-arquivos-de-docking-wizard)
9. [Opções avançadas (interações e filtros)](#9-opções-avançadas-interações-e-filtros)
10. [Diferenças entre Windows e Linux](#10-diferenças-entre-windows-e-linux)
11. [Perguntas frequentes](#11-perguntas-frequentes)
12. [Solução de problemas](#12-solução-de-problemas)

---

## 1. O que é o LUNA e o que esta GUI faz

**LUNA** (_Ligand-protein ANAlysis toolkit_) é uma biblioteca Python para calcular interações
proteína-ligante em larga escala. A partir de uma proteína (arquivo PDB) e de um conjunto de
ligantes (arquivo MOL2 ou SDF), o LUNA:

- Detecta e classifica as interações moleculares (ligações de hidrogênio, π-stacking,
  interações iônicas, cátion-π, hidrofóbicas, etc.)
- Gera **Interaction Fingerprints (IFP)** — vetores numéricos que codificam o padrão
  de interação de cada ligante com a proteína
- Calcula **matrizes de similaridade** de Tanimoto entre fingerprints
- Exporta **sessões PyMOL (.pse)** para visualização 3-D

Esta **GUI** (interface gráfica) abstrai completamente a linha de comando do LUNA.
Você configura tudo por menus e botões; a GUI traduz para os comandos corretos internamente.

---

## 2. O que você NÃO precisa fazer

| Tarefa | Status |
|---|---|
| Baixar ou instalar o LUNA manualmente | ✅ A GUI instala automaticamente |
| Criar o ambiente conda `luna-env` | ✅ A GUI cria e configura ao clicar em "Instalar LUNA" |
| Editar arquivos de configuração | ✅ Tudo é feito pelos campos da GUI |
| Entender a API Python do LUNA | ✅ A GUI cuida disso |

O **único passo manual** é instalar o Miniconda e criar o ambiente da GUI em si
(descrito abaixo — são três comandos que você executa uma única vez).

---

## 3. Requisitos

### Hardware
- Processador: qualquer CPU moderna (≥ 2 núcleos recomendado)
- RAM: mínimo 8 GB (16 GB recomendado para bibliotecas grandes)
- Espaço em disco: ~4 GB (para o ambiente conda com LUNA + dependências)
- Conexão com a internet na **primeira execução** (para baixar os pacotes)

### Software
- **Windows 10/11** (64 bits) ou **Linux** (Ubuntu 20.04+, Fedora 36+ ou equivalente)
- **Miniconda** (gerenciador de pacotes Python/conda) — veja a seção 4

### Arquivos de entrada que você fornece
- Arquivo da proteína em formato **PDB** (`.pdb`)
- Arquivo de ligantes em formato **MOL2** (`.mol2`) ou **SDF** (`.sdf`)
  - Pode ser um arquivo multi-molécula (todos os ligantes em um só arquivo)
  - Ou uma pasta com um arquivo `.mol2` por ligante

---

## 4. Instalação — passo a passo

### 4.1 Windows

#### Passo 1 — Instalar o Miniconda

1. Acesse: `https://docs.conda.io/en/latest/miniconda.html`
2. Baixe o instalador **Miniconda3 — Windows — 64-bit**
3. Execute o `.exe` e siga o instalador (opções padrão são suficientes)
4. Ao final, o Windows terá um atalho chamado **"Anaconda Prompt"** no Menu Iniciar

#### Passo 2 — Criar o ambiente da GUI

Abra o **Anaconda Prompt** (Menu Iniciar → procure por "Anaconda Prompt") e
execute os três comandos abaixo, um por vez:

```
conda create -n luna-gui python=3.11 -y
conda activate luna-gui
pip install -r requirements.txt
```

Aguarde cada comando concluir antes de digitar o próximo.
O arquivo `requirements.txt` instala também `numpy` e `scipy`, necessários
para as exportações gráficas e para o clustering hierárquico na aba Resultados.

#### Passo 3 — Launcher Windows

O arquivo `windows\run_gui.bat` agora tenta localizar automaticamente o
`python.exe` do ambiente `luna-gui`.

Se for necessário forçar manualmente outro executável, abra um Prompt e rode:
```
set HIP2LINTERACTOMICS_GUI_PYTHON=C:\Users\SeuNome\miniconda3\envs\luna-gui\python.exe
windows\run_gui.bat
```

#### Passo 4 — Iniciar a GUI

Dê um duplo-clique em `windows\run_gui.bat`.

Uma janela preta abrirá brevemente e em seguida a GUI será exibida.

---

### 4.2 Linux

#### Passo 1 — Instalar o Miniconda

```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
# Siga as instruções; ao final, feche e reabra o terminal
```

#### Passo 2 — Instalar dependências gráficas do sistema

```bash
# Ubuntu / Debian
sudo apt install libxcb-cursor0 libxkbcommon-x11-0 libgl1 libegl1

# Fedora / RHEL
sudo dnf install xcb-util-cursor libxkbcommon-x11 mesa-libGL mesa-libEGL
```

#### Passo 3 — Criar o ambiente da GUI

```bash
conda create -n luna-gui python=3.11 -y
conda activate luna-gui
pip install -r requirements.txt
```

O arquivo `requirements.txt` instala também `numpy` e `scipy`, necessários
para as exportações gráficas e para o clustering hierárquico na aba Resultados.

#### Passo 4 — Tornar o launcher executável

```bash
chmod +x linux/run_gui.sh
```

#### Passo 5 — Iniciar a GUI

```bash
./linux/run_gui.sh
```

---

## 5. Primeira execução — aba Setup

Quando a GUI abrir pela primeira vez, você estará na aba **"1. Setup"**.

Esta aba tem dois propósitos:
1. Verificar se o Miniconda/conda está disponível no seu sistema
2. Instalar o LUNA e todas as suas dependências no ambiente `luna-env`

### O que a GUI instala automaticamente

Ao clicar em **"Instalar LUNA"**, a GUI executa (sem intervenção do usuário):

1. Cria o ambiente conda `luna-env` com Python 3.9
2. Instala via conda-forge: RDKit, OpenBabel, PyMOL (open-source), Biopython, NumPy, pandas, scipy, matplotlib, seaborn, networkx
3. Instala via pip: pdbecif, mmh3, xopen, colorlog
4. Instala o LUNA em si via pip
5. Aplica um patch automático para compatibilidade com Windows (correção de overflow de inteiros)

**Este processo leva entre 10 e 30 minutos** dependendo da velocidade da internet e do computador.
Todo o progresso é exibido no painel de log da aba Setup em tempo real.

### Atenção

- **Não feche a GUI** durante a instalação
- Se a instalação for interrompida, você pode clicar em "Instalar LUNA" novamente — a GUI
  verifica o que já foi instalado
- Após a instalação, o botão ficará verde e mostrará o caminho do Python do `luna-env`

---

## 6. Fluxo básico de uso

```
[Setup] → [Projeto] → [Análises] → [Executar] → [Resultados]
```

1. **Setup** (uma única vez) — instalar o LUNA
2. **Projeto** — escolher proteína, ligantes e pasta de saída
3. **Análises** — escolher quais cálculos realizar (IFP, matriz de similaridade, sessões PyMOL)
4. **Executar** — iniciar o LUNA e acompanhar o progresso no log
5. **Resultados** — visualizar fingerprints, heatmaps, estatísticas, clusters, abrir sessões PyMOL e exportar gráficos

---

## 7. Abas em detalhe

### 7.1 Aba Setup

| Elemento | Função |
|---|---|
| **Detectar conda** | Verifica se o Miniconda está instalado e mostra o caminho |
| **Instalar LUNA** | Instala/atualiza o LUNA no `luna-env` automaticamente |
| **Verificar LUNA** | Confirma se o LUNA está operacional sem reinstalar |
| Painel de log | Exibe todo o progresso da instalação em tempo real |

---

### 7.2 Aba Projeto

Esta aba define as entradas para a análise.

| Campo | O que inserir |
|---|---|
| **Proteína (PDB)** | Clique em "Procurar..." e selecione o arquivo `.pdb` da proteína já preparada (com hidrogênios, sem ligante de referência se não for desejado) |
| **Ligantes (MOL2/SDF)** | Clique em "Arquivo..." para um único arquivo multi-molécula, ou "Pasta MOL2..." para consolidar todos os `.mol2` de uma pasta |
| **Diretório de trabalho** | Pasta onde os resultados serão salvos. Pode ser uma pasta nova — a GUI a criará |
| **Incluir águas (HOH)** | Marque se quiser incluir moléculas de água nas interações (análise "hidratada"). Ativa o modo avançado automaticamente |
| **Lista de ligantes** | Após carregar o arquivo de ligantes, todos os nomes aparecem aqui. Use "Selecionar tudo" para analisar a biblioteca completa |

**Filtro de ligantes:** o campo de texto no topo da lista filtra por nome em tempo real.
Útil quando a biblioteca tem centenas de compostos e você quer selecionar apenas um subconjunto.

**Botão "Preparar arquivos de docking...":** se seus arquivos MOL2 vêm diretamente
de um programa de docking (como GOLD) com proteína e ligante no mesmo arquivo,
use este wizard antes de prosseguir. Veja a [Seção 8](#8-preparar-arquivos-de-docking-wizard).

---

### 7.3 Aba Análises

Define o que o LUNA calculará para cada ligante.

#### Fingerprints de interação (IFP)

Gera vetores numéricos descrevendo o padrão de interação de cada ligante.

| Parâmetro | Descrição | Valor padrão |
|---|---|---|
| **Tipo** | EIFP (baseado em ECFP), HIFP (hierárquico), FIFP (baseado em fragmentos) | EIFP |
| **Levels** | Número de camadas de expansão radial | 2 |
| **Radius step** | Incremento do raio a cada nível (Å) | 5.73 |
| **Length** | Tamanho do vetor fingerprint (bits) | 4096 |
| **Bit fingerprint** | Marque para fingerprint binário (0/1). Desmarcado = fingerprint de contagem | Desmarcado |

#### Matriz de similaridade (Tanimoto)

Gera uma matriz N×N com os coeficientes de Tanimoto entre todos os pares de ligantes.
Útil para identificar agrupamentos (clusters) de moléculas com perfil de interação similar.

#### Clusters

A partir da matriz de similaridade, a GUI também executa **clustering hierárquico**
dos ligantes e mostra:

- dendrograma;
- matriz reordenada por cluster;
- tabela `ligante → cluster`;
- exportação de clusters em `CSV`.

#### Exportar sessões PyMOL (.pse)

Gera um arquivo `.pse` por ligante, pronto para abrir no PyMOL e visualizar as interações em 3D.

**Filtrar por tipo de interação:** se quiser PSEs contendo apenas um tipo específico
(ex: só "Cation-pi"), marque a caixa abaixo e selecione os tipos na lista.

#### Filtrar por binding modes (.cfg)

Para filtrar resultados com base em padrões de interação pré-definidos.
Use o botão "Editor visual" para criar o arquivo de configuração.

#### Opções avançadas

Veja a [Seção 9](#9-opções-avançadas-interações-e-filtros).

---

### 7.4 Aba Executar

| Campo | Descrição |
|---|---|
| **Núcleos (--nproc)** | Quantos núcleos do processador usar em paralelo. No Windows, mantenha em **1** (limitação técnica). No Linux, use quantos núcleos tiver disponíveis |
| **Sobrescrever projeto** | Se a pasta de trabalho já contiver resultados anteriores, marque para substituí-los |
| **▶ Executar LUNA** | Inicia a análise. A GUI mostra o comando exato no campo abaixo |
| **■ Cancelar** | Interrompe a análise em andamento |
| Painel de log | Exibe saída do LUNA em tempo real. Útil para acompanhar o progresso e identificar erros |

**Tempo estimado:**
- 1 ligante: segundos a poucos minutos
- 100 ligantes: 10–30 minutos (Linux, 4 núcleos)
- 1.000 ligantes: 1–3 horas (Linux, 8 núcleos)
- No Windows (nproc=1), multiplicar por 4–8×

Quando concluído com sucesso, a GUI muda automaticamente para a aba **Resultados**.

---

### 7.5 Aba Resultados

Após uma análise, esta aba é populada automaticamente. Ela contém quatro sub-abas:

#### Fingerprints

Tabela com o conteúdo do arquivo `ifp.csv` gerado pelo LUNA:
- Coluna `ligand_id`: identificador do ligante
- Coluna `on_bits`: posições do vetor fingerprint com valor ≠ 0

Ajuste "Linhas a exibir" para navegar por bibliotecas grandes.

#### Matriz de similaridade

Heatmap colorido (escala de 0 a 1) com os valores de Tanimoto entre todos os ligantes.
- Cores quentes (amarelo/verde) = alta similaridade de interação
- Cores frias (azul/roxo) = baixa similaridade

Disponível apenas se "Matriz de similaridade" foi marcada na aba Análises.

#### Heatmap por tipo

Para cada **tipo de interação** (Hydrogen bond, Cation-pi, Hydrophobic, etc.),
exibe um heatmap com **ligantes × resíduos**, mostrando onde cada ligante
interage com quais aminoácidos da proteína.

1. Clique em **"Calcular heatmap"** (requer o `luna-env`)
2. Use o menu suspenso para alternar entre tipos de interação
3. Colunas de resíduos sem nenhuma interação são ocultadas automaticamente

Este é o gráfico mais biologicamente informativo — revela o perfil de interação
de cada ligante com cada sítio da proteína.

#### Exportações

Na aba Resultados, a GUI permite:

- exportar o gráfico visível em `PNG`, `SVG` ou `PDF`;
- exportar os clusters calculados em `CSV`;
- embutir os principais gráficos e clusters em um relatório HTML.

#### Sessões PyMOL

Lista todos os arquivos `.pse` gerados. Dê duplo-clique em um arquivo
para abri-lo diretamente no PyMOL.

Disponível apenas se "Exportar sessões PyMOL" foi marcada na aba Análises.

**Botão "Exportar relatório HTML":** gera um relatório HTML autocontido
(todos os gráficos embutidos) com os principais resultados da análise,
incluindo clusters quando disponíveis.

---

### 7.6 Aba Histórico

Lista todos os projetos executados anteriormente (workdirs).
- Dê um duplo-clique ou clique em **"Carregar projeto selecionado"** para reabrir
  qualquer projeto anterior — a GUI recarrega a configuração e os resultados automaticamente
- Clique em **"Remover da lista"** para limpar entradas antigas (não apaga os arquivos)

---

## 8. Preparar arquivos de docking (wizard)

Programas de docking molecular como **GOLD** e **AutoDock Vina** frequentemente
geram arquivos MOL2 onde **proteína e ligante estão combinados no mesmo arquivo**.
O LUNA precisa desses dois componentes separados.

Este wizard automatiza exatamente essa separação.

### Quando usar

Use este wizard **antes** de carregar os arquivos na aba Projeto, quando:
- Seus arquivos `.mol2` foram gerados pelo GOLD (saída típica: `protein_ligand.mol2`)
- O arquivo contém átomos da proteína (numerados de 1 a N) seguidos dos átomos do ligante

### Como usar

1. Na aba **Projeto**, clique em **"Preparar arquivos de docking..."**
2. Em **"Pasta de origem"**, selecione a pasta com os arquivos MOL2 do docking
3. Em **"Último átomo da proteína"**, insira o número do último átomo que pertence à proteína
   - Para descobrir esse número: abra um dos arquivos MOL2 em um editor de texto
     e localize onde os átomos da proteína terminam e os do ligante começam
   - Exemplo: se os átomos da proteína vão de 1 a 4068, insira `4068`
4. Clique em **"Executar preparação"**
5. Serão criadas duas subpastas:
   - `proteinas_pdb/` — um arquivo `.pdb` por pose de docking
   - `ligantes_mol2/` — um arquivo `_ligand.mol2` por pose
6. Se a opção **"Usar as pastas geradas como entradas do projeto"** estiver marcada
   (padrão), a GUI automaticamente consolida os ligantes e os carrega na aba Projeto

### O que o wizard faz internamente

- Remove **átomos LP** (lone pairs) que causam erros no OpenBabel
- **Renumera** os átomos do ligante a partir de 1 (necessário para o LUNA)
- Recalcula o cabeçalho do MOL2 (número de átomos e ligações)
- Gera os arquivos PDB da proteína com formatação padrão

---

## 9. Opções avançadas (interações e filtros)

As opções avançadas ficam na seção inferior da aba **Análises**, dentro da caixa
**"Opções avançadas — DefaultInteractionConfig + InteractionCalculator"**.

Ao marcar esta caixa (ou alterar qualquer campo dela), a GUI **automaticamente**
muda o motor de execução de CLI para a API Python do LUNA — isso permite
controle fino, porém o tempo de execução pode ser ligeiramente maior.

### DefaultInteractionConfig — thresholds de distância

| Campo | Significado |
|---|---|
| **HB max D–A dist (Å)** | Distância máxima doador–aceitador para classificar como ligação de hidrogênio. Padrão LUNA: ~3.9 Å. O notebook de referência usa 4.0 Å |
| **HB min D–H–A ângulo (°)** | Ângulo mínimo D–H–A para ligação de hidrogênio |
| **Hidrofóbica max dist (Å)** | Distância máxima C–C para contato hidrofóbico |
| **π–π max dist (Å)** | Distância máxima centro–centro para π-stacking |

Deixe em `(padrão)` para usar os valores originais do LUNA.

### InteractionCalculator — flags de comportamento

| Flag | Descrição |
|---|---|
| **add_proximal** | Inclui contatos baseados apenas em proximidade (sem classificação química) |
| **add_atom_atom** | Inclui interações atômicas genéricas (sem característica de grupo funcional) |
| **add_dependent_inter** _(ligado por padrão)_ | Inclui interações dependentes de outras, como pontes de hidrogênio mediadas por água e pontes de sal |
| **add_h2o_pairs_with_no_target** _(ligado por padrão)_ | Inclui pares água–átomo onde a água não tem alvo definido |
| **ignore_self_inter** _(ligado por padrão)_ | Ignora interações intramoleculares (recomendado para análise proteína-ligante) |

---

## 10. Diferenças entre Windows e Linux

| Característica | Windows | Linux |
|---|---|---|
| **Paralelização (nproc)** | Fixo em **1** — bug técnico do multiprocessing no modo `spawn` | Funciona com quantos núcleos quiser |
| **Tempo de análise** | ~4–8× mais lento que Linux para mesma carga | Referência |
| **Compatibilidade** | Total para análises de até centenas de ligantes | Total, incluindo bibliotecas de milhares |
| **Launcher** | `windows\run_gui.bat` (duplo-clique) | `linux/run_gui.sh` (terminal) |
| **DLL conflicts** | PATH sanitizado pelo launcher | Sem problema |
| **Validação do env da GUI** | Manual | Automática no launcher |

**Recomendação:** para análises de produção com > 500 ligantes, use Linux ou WSL2
(Windows Subsystem for Linux).

---

## 11. Perguntas frequentes

**Q: Preciso comprar o LUNA?**
Não. O LUNA é um software de código aberto, gratuito, disponível no PyPI.
A GUI o instala automaticamente.

**Q: Preciso ter o PyMOL instalado separadamente?**
Não — a GUI instala o `pymol-open-source` dentro do `luna-env` automaticamente
(para gerar os arquivos `.pse`). Para **abrir** os arquivos `.pse` e visualizá-los,
você precisa de qualquer versão do PyMOL instalada no seu sistema
(inclusive a versão gratuita `pymol-open-source`).

**Q: Meu arquivo de ligantes tem 5.000 moléculas. A GUI aguenta?**
Sim. O LUNA foi projetado para análises em larga escala. No Linux com nproc > 1,
milhares de ligantes são factíveis em poucas horas. No Windows, mantenha expectativas
realistas para bibliotecas muito grandes (horas a dias com nproc=1).

**Q: Posso rodar a GUI em um servidor remoto (SSH)?**
Sim, com encaminhamento de X11 (`ssh -X`) ou usando um servidor VNC/X2Go.
A GUI usa Qt que precisa de um servidor gráfico.

**Q: Meu arquivo PDB tem múltiplos modelos. O que o LUNA usa?**
O LUNA usa o primeiro modelo do PDB. Prepare o arquivo manualmente se necessário
(remova `MODEL`/`ENDMDL` extras).

**Q: O que é o arquivo `entries.txt`?**
É gerado automaticamente pela GUI com os nomes dos ligantes selecionados.
O LUNA usa esse arquivo para saber quais moléculas do MOL2 analisar.
Você não precisa criá-lo manualmente.

**Q: Posso reutilizar uma análise anterior e adicionar novos ligantes?**
Sim — use o recurso **Fork de projeto existente** na aba Análises.
Aponte para a pasta do projeto anterior e a GUI reutilizará os resultados já calculados.

---

## 12. Solução de problemas

### A GUI não abre / erro de DLL (Windows)

Verifique se o `python.exe` do env `luna-gui` realmente importa o PyQt6:
```
conda activate luna-gui
python -c "import sys; print(sys.executable); import PyQt6.QtWidgets; print('PyQt6 OK')"
```

Se necessário, force o launcher a usar esse executável:
```
set HIP2LINTERACTOMICS_GUI_PYTHON=C:\Users\SeuNome\miniconda3\envs\luna-gui\python.exe
windows\run_gui.bat
```

### "conda não encontrado" na aba Setup

O Miniconda não está instalado ou não está no PATH.
Reinstale o Miniconda e certifique-se de marcar **"Add to PATH"** durante a instalação
(ou abra sempre pelo **Anaconda Prompt**).

### Instalação do LUNA trava na etapa "conda create"

Verifique sua conexão com a internet. Canais conda-forge precisam de acesso à internet.

### LUNA gera 0 interações para meu ligante

Causas comuns:
1. O ligante está posicionado longe da proteína no PDB/MOL2 — verifique visualmente no PyMOL
2. O arquivo MOL2 do ligante tem átomos LP não removidos — use o wizard de preparação
3. O nome do ligante no `entries.txt` não bate com o nome no MOL2 — verifique na
   lista de ligantes da aba Projeto

### "UnicodeEncodeError" no log (Windows)

Isso foi corrigido automaticamente. Se ainda aparecer, reabra a GUI — ela
força UTF-8 no processo filho (LUNA) via variável de ambiente.

### "cannot pickle '_thread.lock'" no log (Windows)

Reduza `nproc` para **1** na aba Executar. Este é um bug do multiprocessing
do LUNA no Windows que não afeta a análise — apenas impede paralelização.

### IFP gerado mas resultados parecem errados

Um IFP com poucas interações pode ser esperado para ligantes pequenos ou
mal posicionados. Verifique visualmente abrindo a sessão PyMOL do ligante.

---

## Estrutura de arquivos gerados

Após uma análise, a pasta de trabalho conterá:

```
<workdir>/
├── .luna_gui.json              ← configuração do projeto (reabrir pelo Histórico)
├── entries.txt                 ← lista de ligantes analisados
├── _luna_api_runner.py         ← (apenas modo avançado) script API gerado
├── _luna_api_params.json       ← (apenas modo avançado) parâmetros em JSON
└── results/
    ├── fingerprints/
    │   └── ifp.csv             ← fingerprints de todos os ligantes
    ├── pse/
    │   ├── ligante1.pse        ← sessão PyMOL por ligante
    │   └── ligante2.pse
    └── <ligante>.pkl.gz        ← resultados completos por ligante (para o Heatmap)
sim_matrix.csv                  ← matriz de similaridade Tanimoto (se habilitado)
```

---

*Para dúvidas sobre o LUNA em si, consulte a documentação oficial:*  
*https://luna-toolkit.readthedocs.io*
