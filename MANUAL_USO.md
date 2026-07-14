# Manual de Uso do HIP2LInterActomics

Este manual explica como instalar, abrir, configurar, executar e interpretar os resultados do HIP2LInterActomics. Ele foi escrito para uso direto por pesquisadores que querem analisar interações proteína-ligante, poses de docking, complexos hidratados ou trajetórias de dinâmica molecular com o LUNA sem precisar controlar a API manualmente.

## Sumário

1. [Visão Geral](#visão-geral)
2. [Instalação](#instalação)
3. [Abrindo o Aplicativo](#abrindo-o-aplicativo)
4. [Fluxo Recomendado](#fluxo-recomendado)
5. [Aba 1. Início](#aba-1-início)
6. [Aba 2. Projeto](#aba-2-projeto)
7. [Aba 3. Análises](#aba-3-análises)
8. [Aba 4. Executar](#aba-4-executar)
9. [Aba 5. Resultados](#aba-5-resultados)
10. [Aba 6. Histórico](#aba-6-histórico)
11. [Execução por Terminal](#execução-por-terminal)
12. [Arquivos Gerados](#arquivos-gerados)
13. [Boas Práticas](#boas-práticas)
14. [Solução de Problemas](#solução-de-problemas)

## Visão Geral

O HIP2LInterActomics é uma interface gráfica e um fluxo terminal para executar análises do LUNA em complexos moleculares. O aplicativo organiza a preparação dos dados, a execução do LUNA e a interpretação dos resultados em tabelas, gráficos, mapas de calor, fingerprints de interação, sessões PyMOL e relatórios.

O programa trabalha com dois ambientes Conda separados:

- `luna-gui`: ambiente da interface gráfica, com PyQt6, matplotlib, numpy, scipy e scikit-learn.
- `luna-env`: ambiente usado para executar LUNA, RDKit, Open Babel, PyMOL e dependências químicas.

Essa separação reduz conflitos entre bibliotecas gráficas e bibliotecas nativas usadas pelo LUNA.

## Instalação

### Linux

Use o instalador completo em uma cópia do repositório:

```bash
chmod +x dist/linux/install_hip2linteractomics.sh
./dist/linux/install_hip2linteractomics.sh
```

O instalador detecta Conda, instala Miniforge se necessário, cria `luna-gui`, cria `luna-env`, instala LUNA/RDKit/Open Babel/PyMOL e aplica o patch de compatibilidade do projeto.

Dependências gráficas recomendadas:

```bash
# Debian / Ubuntu
sudo apt install libxcb-cursor0 libxkbcommon-x11-0 libgl1 libegl1

# Fedora / RHEL
sudo dnf install xcb-util-cursor libxkbcommon-x11 mesa-libGL mesa-libEGL
```

### Windows

Use o instalador PowerShell:

```powershell
powershell.exe -NoProfile -ExecutionPolicy Bypass -File .\dist\windows\install_hip2linteractomics.ps1
```

O instalador cria os ambientes `luna-gui` e `luna-env`, instala as dependências e cria um atalho `HIP2LInterActomics.cmd` na área de trabalho.

### Snap Linux

O formato `.snap` é exclusivo do Linux. Para gerar o pacote:

```bash
sudo snap install snapcraft --classic
./dist/linux/build_snap.sh
```

Depois, instale localmente:

```bash
sudo snap install ./hip2linteractomics_1.0.0_amd64.snap --classic --dangerous
hip2linteractomics
```

No Windows, `.snap` não é o formato apropriado. Use o instalador PowerShell ou rode a versão Linux em WSL2.

## Abrindo o Aplicativo

### Linux

```bash
./dist/linux/run_gui.sh
```

Variáveis úteis:

```bash
HIP2LINTERACTOMICS_GUI_ENV=luna-gui ./dist/linux/run_gui.sh
HIP2LINTERACTOMICS_GUI_CONDA=/caminho/para/conda ./dist/linux/run_gui.sh
HIP2LINTERACTOMICS_GUI_PYTHON=/caminho/para/python ./dist/linux/run_gui.sh
```

### Windows

```powershell
.\dist\windows\run_gui.bat
```

Se precisar forçar o Python da GUI:

```powershell
$env:HIP2LINTERACTOMICS_GUI_PYTHON="C:\caminho\para\luna-gui\python.exe"
.\dist\windows\run_gui.bat
```

## Fluxo Recomendado

1. Abra a aba `1. Início` e verifique Conda/LUNA.
2. Configure entradas na aba `2. Projeto`.
3. Escolha fingerprints, filtros e análises na aba `3. Análises`.
4. Execute o LUNA na aba `4. Executar`.
5. Interprete resultados, gráficos, FP e PyMOL na aba `5. Resultados`.
6. Use `6. Histórico` para reabrir projetos.

## Aba 1. Início

Essa aba verifica o ambiente e instala o LUNA quando necessário.

Campos e botões principais:

- `Verificar novamente`: refaz a detecção de Conda, `luna-env` e LUNA.
- `Baixar Miniconda`: abre a página oficial de instalação do Conda.
- `Instalar LUNA`: cria ou atualiza o ambiente `luna-env`.
- Painel de log: mostra cada etapa da verificação e instalação.

Quando a instalação termina, a GUI passa a ter um Python de cálculo separado, normalmente dentro de `luna-env`.

## Aba 2. Projeto

Essa aba define os arquivos de entrada e a pasta de trabalho.

### Entradas

- Proteína/receptor: normalmente um arquivo `.pdb`.
- Ligantes: arquivo `.mol2`, `.sdf`, `.sd`, `.mol`, `.pdb`, `.ent` ou pasta com ligantes.
- Workdir: pasta onde todos os outputs serão salvos.
- Lista de ligantes: permite selecionar todos ou apenas parte da biblioteca.

### Análise de Trajetória ou Poses

Marque essa opção quando cada entrada representa um frame de dinâmica molecular ou uma pose de docking. Nesse modo, os gráficos respeitam a ordem de frames/poses e não reordenam moléculas por rótulos externos.

Quando esse modo está ativo, a aba `5. Resultados > Estatísticas` também pode mostrar o gráfico por átomos do ligante, com percentuais de interação por tipo.

### Pré-processamento

Use o pré-processamento quando os arquivos precisam ser separados ou normalizados antes da análise.

Para arquivos PDB, o aplicativo tenta respeitar:

- `ATOM`: proteína.
- `HETATM`: ligantes, águas e outros componentes.
- `HOH`, `WAT` ou `WTM`: moléculas de água.

Quando a distinção `ATOM/HETATM` não é suficiente, a GUI usa o nome do resíduo. Resíduos que não são aminoácidos essenciais/especiais e não são água são tratados como candidatos a ligante.

Problemas recuperáveis são registrados no `.log`, e o fluxo continua quando isso não compromete o cálculo das interações.

## Aba 3. Análises

Essa aba define o que será calculado pelo LUNA.

### Fingerprints de Interação

Parâmetros principais:

- Tipo: `HIFP`, `EIFP`, `FIFP` ou todos.
- Níveis: número de níveis/shells considerados.
- Raio: passo radial usado na geração dos shells.
- Tamanho: número de posições do vetor de fingerprint.
- Binário ou contagem: define se o valor final é presença/ausência ou número de ocorrências.

Os seeds de importância dos fingerprints são salvos como `seed_ifp_[H/E/F]_importance.txt`, quando usados na etapa de importância.

### Matriz de Similaridade

Gera matrizes de similaridade entre ligantes a partir dos fingerprints. Elas alimentam o mapa de similaridade e os clusters.

### Rótulos

Permite carregar arquivo CSV/TSV com rótulos externos.

- Regressão: gráficos de moléculas podem ser ordenados em ordem decrescente do valor do rótulo.
- Classificação: gráficos de moléculas podem ser agrupados por classe.
- Trajetória/poses: a ordem dos frames/poses tem prioridade sobre rótulos.

### Filtrar por Binding Modes

Permite criar filtros de interação por meio do editor visual. Use essa opção para selecionar moléculas que cumprem padrões específicos de interação com resíduos ou tipos de contato.

### Filtrar PSE por Tipo de Interação

Controla quais tipos de interação aparecem nas sessões PyMOL geradas. O checkbox geral seleciona ou desseleciona todos os tipos.

### Arquivo Completo de Interações

Permite apontar para um `.cfg` do LUNA com configuração completa de interações. Use quando precisar alterar parâmetros finos além dos controles rápidos da GUI.

## Aba 4. Executar

Essa aba executa o projeto.

Campos principais:

- `Núcleos (--nproc)`: número de núcleos usados.
- `Sobrescrever projeto existente`: permite reusar a mesma pasta e substituir saídas anteriores.
- `Executar LUNA`: inicia a análise.
- `Cancelar`: interrompe a execução.
- Log: mostra mensagens do processo em tempo real.

No Linux, `nproc` pode usar múltiplos núcleos. No Windows nativo, a GUI limita `nproc=1` porque o LUNA usa `multiprocessing` em modo `spawn`, que pode falhar com `nproc > 1`. Para paralelismo real no Windows, use WSL2/Linux.

## Aba 5. Resultados

Essa aba carrega resultados gerados pelo LUNA, projetos antigos e workdirs calculados por terminal.

### Estatísticas

Mostra estatísticas globais de interações e gráficos de distribuição por tipo/resíduo.

Quando a análise de trajetória/poses está ativa, a GUI também mostra o gráfico por átomos do ligante, com:

- eixo X: átomos pesados da molécula pequena;
- eixo Y: porcentagem;
- barras: tipos de interação realizados por cada átomo.

Quando possível, a estrutura 2D do ligante com IDs de átomos é exibida para ajudar a interpretar o gráfico.

### Mapa de Calor por Tipo

Mostra mapas de calor separados por tipo de interação, relacionando ligantes/poses com resíduos.

### Mapa de Calor Completo Ligantes x Resíduos

Mostra uma visão integrada de múltiplos tipos de interação em uma matriz ligante x resíduo.

### Triagem Virtual Racional

Essa aba substitui a antiga aba de sessões PyMOL. Ela permite aplicar filtros dinâmicos de binding modes aos resultados e selecionar ligantes com padrões de interação desejados.

### Fingerprints

Mostra as matrizes de fingerprints geradas pelo LUNA. Quando o projeto usa o fluxo atual, os labels seguem o formato:

```text
[fingerprint]_[nível]
```

Exemplo:

```text
1994_2
```

### Matriz de Similaridade

Mostra a similaridade entre ligantes calculada a partir dos fingerprints.

### Clusters

Mostra agrupamentos hierárquicos de ligantes com perfis de interação semelhantes.

### Análises FP

Essa aba calcula e mostra a importância dos fingerprints.

O fluxo geral é:

1. Carrega a matriz `entries x FP`.
2. Carrega detalhes dos shells e ocorrências.
3. Classifica a natureza/tipo da feature.
4. Verifica colisões por classe.
5. Atribui a classe predominante quando confiável.
6. Aplica z-score ou Otsu para definir limiares.
7. Marca features instáveis como `unreliable feature by class`.
8. Busca níveis/shells das features confiáveis por classe.
9. Resolve colisões por nível.
10. Marca features instáveis como `unreliable feature by level`.
11. Salva labels no formato `FP_nível`.
12. Separa features por nível.
13. Treina um modelo por nível.
14. Calcula importância das features em cada modelo.
15. Transforma z-scores dos coeficientes de importância em p-values pela equação de Keiser and Hert [1].
16. Seleciona features importantes por p-value.
17. Gera tabelas, gráficos e mapas com labels `FP_nível`.

Termos importantes:

- Features confiáveis por classe: bits cuja classe/tipo de interação pode ser atribuída de forma confiável.
- Features confiáveis por nível: bits cuja origem por nível/shell pode ser atribuída de forma confiável.
- Elegíveis para importância: features confiáveis por classe e por nível, usadas nos modelos.
- Features importantes: elegíveis que passam pelo critério estatístico configurado.

### Sessão FP em PyMOL

Gera sessões PyMOL para visualizar os shells que originaram um fingerprint em uma molécula específica.

Seletores pesquisáveis:

- Tipo de fingerprint.
- Feature/fingerprint.
- Molécula.

O PSE gerado inclui a indicação do número do shell para facilitar a inspeção visual.

## Aba 6. Histórico

Mostra workdirs usados anteriormente.

Você pode:

- reabrir projetos já calculados;
- carregar resultados antigos;
- remover entradas do histórico sem apagar os arquivos do disco.

## Execução por Terminal

O script `hipplinteractomics_terminal.py` executa o mesmo fluxo da GUI a partir de um arquivo JSON ou de um dicionário Python literal.

Exemplo:

```bash
conda activate luna-gui
python hipplinteractomics_terminal.py meu_projeto.json
```

Exemplo mínimo de configuração:

```json
{
  "protein_file": "/caminho/receptor.pdb",
  "ligand_file": "/caminho/ligantes",
  "selected_ligands": "ALL",
  "workdir": "/caminho/workdir",
  "trajectory_analysis": false,
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
  "overwrite": true
}
```

Projetos calculados por terminal podem ser abertos depois na aba `5. Resultados`.

## Arquivos Gerados

Estrutura típica:

```text
workdir/
  .luna_gui.json
  entries.txt
  _luna_api_runner.py
  _luna_api_params.json
  hipplinteractomics_terminal.log
  sim_matrix_E.csv
  sim_matrix_H.csv
  sim_matrix_F.csv
  results/
    analysis_summary.json
    residue_matrix.json
    fingerprints/
      ifp.csv
      ifp_E.csv
      ifp_H.csv
      ifp_F.csv
      fp_detail_E.json
      fp_detail_H.json
      fp_detail_F.json
      seed_ifp_E_importance.txt
      seed_ifp_H_importance.txt
      seed_ifp_F_importance.txt
    interactions/
    pse/
    fp_sessions/
```

Nem todos os arquivos aparecem em todos os projetos. A presença depende das opções marcadas nas abas `2. Projeto` e `3. Análises`.

## Boas Práticas

- Use Linux ou WSL2 para bibliotecas grandes.
- No Windows nativo, mantenha `nproc=1`.
- Guarde o `workdir` completo para reabrir resultados.
- Prefira arquivos PDB/MOL2/SDF bem preparados e com nomes consistentes.
- Para análises hidratadas, revise se águas relevantes estão nomeadas como `HOH`, `WAT` ou `WTM`.
- Para projetos com rótulos externos, confira se a coluna de ID bate com os IDs dos ligantes.
- Para comparar projetos, mantenha os mesmos parâmetros de fingerprint.

## Solução de Problemas

### Conda não encontrado

Confira:

```bash
conda info
```

Ou aponte manualmente:

```bash
export HIP2LINTERACTOMICS_GUI_CONDA=/caminho/para/conda
```

No Windows:

```powershell
$env:HIP2LINTERACTOMICS_GUI_CONDA="C:\caminho\para\conda.exe"
```

### A GUI não abre

Verifique se o ambiente `luna-gui` consegue importar PyQt6:

```bash
conda activate luna-gui
python -c "import PyQt6.QtWidgets; print('PyQt6 OK')"
```

### PyMOL/OpenGL falha no Linux

Instale bibliotecas gráficas:

```bash
sudo apt install libxcb-cursor0 libxkbcommon-x11-0 libgl1 libegl1
```

Em servidor remoto, confirme se existe sessão gráfica, X11, Wayland, VNC, X2Go ou WSLg.

### PSE de fingerprint falha

Isso geralmente indica problema nativo do PyMOL/OpenGL ou dependências incompatíveis no `luna-env`. Verifique:

```bash
conda activate luna-env
python -c "import pymol; print('PyMOL OK')"
python -c "import luna; print(luna.__file__)"
```

### Nenhuma interação foi detectada

Possíveis causas:

- ligante está longe do sítio de ligação;
- receptor ou ligante sem hidrogênios adequados;
- arquivo de entrada mal separado;
- entries não batem com os nomes dos ligantes;
- filtro de binding modes está muito restritivo.

### Projetos antigos sem níveis nos fingerprints

Ao carregar resultados antigos, a GUI tenta reconstruir os níveis assinados e regravar matrizes compatíveis. Se isso não for possível, a análise FP pode ignorar features sem nível confiável.

## Referência

[1] Keiser, M. J.; Hert, J. In Chemogenomics: Methods and Applications; Jacoby, E., Ed.; Methods in Molecular Biology; Humana Press: Totowa, NJ, 2009; pp 195-205.
