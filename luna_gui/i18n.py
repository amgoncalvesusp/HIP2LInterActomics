"""Small runtime translation layer for the LUNA GUI.

The GUI was originally written with direct Portuguese strings.  This module
keeps those strings as the source language and translates widgets and
matplotlib text at runtime, so the scientific/data values remain untouched.
"""
from __future__ import annotations

import re
from typing import Any, Callable


LANGUAGES = {
    "pt": "Português",
    "en": "English",
    "es": "Español",
}

_current_language = "pt"
_qt_hooks_installed = False
_mpl_hooks_installed = False


TRANSLATIONS: dict[str, dict[str, str]] = {
    # Menus, tabs and app-level labels.
    "&Arquivo": {"en": "&File", "es": "&Archivo"},
    "Aj&uda": {"en": "&Help", "es": "Ay&uda"},
    "&Idioma": {"en": "&Language", "es": "&Idioma"},
    "Abrir projeto...": {"en": "Open project...", "es": "Abrir proyecto..."},
    "Salvar projeto": {"en": "Save project", "es": "Guardar proyecto"},
    "Sair": {"en": "Exit", "es": "Salir"},
    "Sobre": {"en": "About", "es": "Acerca de"},
    "Sobre o LUNA GUI": {"en": "About LUNA GUI", "es": "Acerca de LUNA GUI"},
    "Documentação LUNA": {"en": "LUNA documentation", "es": "Documentación de LUNA"},
    "Pronto": {"en": "Ready", "es": "Listo"},
    "Defina um workdir antes de salvar.": {"en": "Set a workdir before saving.", "es": "Define un workdir antes de guardar."},
    "Projeto salvo em": {"en": "Project saved in", "es": "Proyecto guardado en"},
    "Sim": {"en": "Yes", "es": "Sí"},
    "Não": {"en": "No", "es": "No"},
    "Padrão LUNA": {"en": "LUNA default", "es": "Predeterminado LUNA"},
    "Idioma alterado para Português": {"en": "Language changed to Portuguese", "es": "Idioma cambiado a portugués"},
    "Idioma alterado para English": {"en": "Language changed to English", "es": "Idioma cambiado a inglés"},
    "Idioma alterado para Español": {"en": "Language changed to Spanish", "es": "Idioma cambiado a español"},
    "1. Setup": {"en": "1. Setup", "es": "1. Configuración"},
    "2. Projeto": {"en": "2. Project", "es": "2. Proyecto"},
    "3. Análises": {"en": "3. Analyses", "es": "3. Análisis"},
    "4. Executar": {"en": "4. Run", "es": "4. Ejecutar"},
    "5. Resultados": {"en": "5. Results", "es": "5. Resultados"},
    "6. Histórico": {"en": "6. History", "es": "6. Historial"},

    # Common controls.
    "Procurar...": {"en": "Browse...", "es": "Examinar..."},
    "Arquivo...": {"en": "File...", "es": "Archivo..."},
    "Pasta...": {"en": "Folder...", "es": "Carpeta..."},
    "Pasta MOL2/SDF": {"en": "MOL2/SDF folder", "es": "Carpeta MOL2/SDF"},
    "Selecionar tudo": {"en": "Select all", "es": "Seleccionar todo"},
    "Nenhum": {"en": "None", "es": "Ninguno"},
    "Atualizar": {"en": "Refresh", "es": "Actualizar"},
    "Cancelar": {"en": "Cancel", "es": "Cancelar"},
    "Fechar": {"en": "Close", "es": "Cerrar"},
    "OK": {"en": "OK", "es": "OK"},
    "Erro": {"en": "Error", "es": "Error"},
    "Sem dados": {"en": "No data", "es": "Sin datos"},
    "Processando...": {"en": "Processing...", "es": "Procesando..."},
    "Processando... (pode levar alguns minutos)": {"en": "Processing... (this may take a few minutes)", "es": "Procesando... (puede tardar unos minutos)"},

    # Setup tab.
    "Verificando ambiente...": {"en": "Checking environment...", "es": "Verificando entorno..."},
    "O que é o LUNA GUI": {"en": "What is LUNA GUI", "es": "Qué es LUNA GUI"},
    (
        "O LUNA GUI é uma interface gráfica para preparar, executar e interpretar análises "
        "baseadas no LUNA. Ele ajuda a avaliar resultados de virtual screening, comparar "
        "poses ou trajetórias de dinâmica molecular, construir modelos de machine learning "
        "com fingerprints de interação e extrair pharmacophoric features de um alvo "
        "específico. Nesta aba você prepara o ambiente: verifica o Conda, instala ou atualiza "
        "o ambiente separado 'luna-env' e garante que as dependências analíticas, como "
        "scikit-learn, estejam disponíveis."
    ): {
        "en": (
            "LUNA GUI is a graphical interface for preparing, running and interpreting analyses "
            "based on LUNA. It helps evaluate virtual screening results, compare poses or "
            "molecular dynamics trajectories, build machine-learning models from interaction "
            "fingerprints and extract pharmacophoric features from a specific target. In this "
            "tab you prepare the environment: check Conda, install or update the isolated "
            "'luna-env' environment and ensure analytical dependencies, such as scikit-learn, "
            "are available."
        ),
        "es": (
            "LUNA GUI es una interfaz gráfica para preparar, ejecutar e interpretar análisis "
            "basados en LUNA. Ayuda a evaluar resultados de virtual screening, comparar poses "
            "o trayectorias de dinámica molecular, construir modelos de machine learning con "
            "fingerprints de interacción y extraer pharmacophoric features de un blanco "
            "específico. En esta pestaña preparas el entorno: verifica Conda, instala o "
            "actualiza el entorno separado 'luna-env' y garantiza que dependencias analíticas, "
            "como scikit-learn, estén disponibles."
        ),
    },
    "Verificar novamente": {"en": "Check again", "es": "Verificar de nuevo"},
    "Baixar Miniconda": {"en": "Download Miniconda", "es": "Descargar Miniconda"},
    "Instalar LUNA (cria luna-env)": {"en": "Install LUNA (creates luna-env)", "es": "Instalar LUNA (crea luna-env)"},
    "Refaz a detecção do Conda, do ambiente 'luna-env' e da instalação do LUNA.": {"en": "Runs Conda, 'luna-env' and LUNA detection again.", "es": "Repite la detección de Conda, 'luna-env' y la instalación de LUNA."},
    "Abre a página oficial do Miniconda para instalar o Conda no sistema.": {"en": "Opens the official Miniconda page to install Conda on this system.", "es": "Abre la página oficial de Miniconda para instalar Conda en el sistema."},
    "Mostra o passo a passo da verificação e da instalação do ambiente.": {"en": "Shows the environment check and installation log.", "es": "Muestra el paso a paso de la verificación e instalación del entorno."},

    # Project tab.
    (
        "Nesta aba você define o modo do projeto, escolhe o diretório de trabalho, pode "
        "pré-processar arquivos de docking e então informa proteína e ligantes para a análise."
    ): {
        "en": "In this tab you define the project mode, choose the working directory, optionally pre-process docking files, then provide the protein and ligands for the analysis.",
        "es": "En esta pestaña defines el modo del proyecto, eliges el directorio de trabajo, puedes preprocesar archivos de docking y luego indicas la proteína y los ligandos para el análisis.",
    },
    "Entradas": {"en": "Inputs", "es": "Entradas"},
    "Modo do projeto:": {"en": "Project mode:", "es": "Modo del proyecto:"},
    "Fork de projeto existente (desmarcado = projeto novo)": {"en": "Fork an existing project (unchecked = new project)", "es": "Fork de un proyecto existente (desmarcado = proyecto nuevo)"},
    "Projeto fonte:": {"en": "Source project:", "es": "Proyecto fuente:"},
    "Diretório de trabalho:": {"en": "Working directory:", "es": "Directorio de trabajo:"},
    "Pré-processamento:": {"en": "Pre-processing:", "es": "Preprocesamiento:"},
    "Preparar arquivos de docking...": {"en": "Prepare docking files...", "es": "Preparar archivos de docking..."},
    "Proteína (PDB):": {"en": "Protein (PDB):", "es": "Proteína (PDB):"},
    "Ligantes (MOL2/SDF):": {"en": "Ligands (MOL2/SDF):", "es": "Ligandos (MOL2/SDF):"},
    "Incluir águas (HOH) — análise hidratada": {"en": "Include waters (HOH) - hydrated analysis", "es": "Incluir aguas (HOH) - análisis hidratado"},
    "Águas detectadas nos inputs: 0": {"en": "Waters detected in inputs: 0", "es": "Aguas detectadas en los inputs: 0"},
    "Análise de trajetória de dinâmica molecular (entradas = frames)": {"en": "Molecular dynamics trajectory analysis (entries = frames)", "es": "Análisis de trayectoria de dinámica molecular (entradas = frames)"},
    "Ligantes detectados": {"en": "Detected ligands", "es": "Ligandos detectados"},
    (
        "Depois de carregar os ligantes, marque apenas os que devem entrar na análise. "
        "O filtro ajuda a localizar nomes específicos sem apagar a lista."
    ): {
        "en": "After loading ligands, check only those that should enter the analysis. The filter helps locate specific names without clearing the list.",
        "es": "Después de cargar los ligandos, marca solo los que deben entrar en el análisis. El filtro ayuda a localizar nombres específicos sin borrar la lista.",
    },
    "Filtrar por nome (texto livre)...": {"en": "Filter by name (free text)...", "es": "Filtrar por nombre (texto libre)..."},
    "0 ligantes": {"en": "0 ligands", "es": "0 ligandos"},
    "Pasta com proteínas PDB": {"en": "Folder with PDB proteins", "es": "Carpeta con proteínas PDB"},
    "Selecionar pasta com proteínas PDB": {"en": "Select folder with PDB proteins", "es": "Seleccionar carpeta con proteínas PDB"},
    "Selecionar proteína": {"en": "Select protein", "es": "Seleccionar proteína"},
    "Selecionar ligantes": {"en": "Select ligands", "es": "Seleccionar ligandos"},
    "Pasta com arquivos MOL2/SDF": {"en": "Folder with MOL2/SDF files", "es": "Carpeta con archivos MOL2/SDF"},
    "Consolidado": {"en": "Consolidated", "es": "Consolidado"},
    "Próximo passo": {"en": "Next step", "es": "Próximo paso"},

    # Analyses tab.
    (
        "Aqui você escolhe quais saídas o LUNA vai gerar. Se estiver em dúvida, deixe os "
        "valores padrão e ative apenas o que pretende inspecionar depois: fingerprints, "
        "similaridade, sessões PyMOL ou filtros específicos."
    ): {
        "en": "Here you choose which outputs LUNA will generate. If unsure, keep the defaults and enable only what you plan to inspect later: fingerprints, similarity, PyMOL sessions or specific filters.",
        "es": "Aquí eliges qué salidas generará LUNA. Si tienes dudas, deja los valores por defecto y activa solo lo que quieras inspeccionar después: fingerprints, similaridad, sesiones PyMOL o filtros específicos.",
    },
    "Preparação de hidrogênios": {"en": "Hydrogen preparation", "es": "Preparación de hidrógenos"},
    (
        "Define se o LUNA deve protonar as estruturas antes da análise. "
        "Quando ativado, o pH abaixo é usado na adição de hidrogênios."
    ): {
        "en": "Defines whether LUNA should protonate structures before analysis. When enabled, the pH below is used when adding hydrogens.",
        "es": "Define si LUNA debe protonar las estructuras antes del análisis. Cuando está activado, el pH indicado se usa al agregar hidrógenos.",
    },
    "Adicionar hidrogênios antes da análise": {"en": "Add hydrogens before analysis", "es": "Agregar hidrógenos antes del análisis"},
    "pH usado na adição de hidrogênios pelo LUNA.": {"en": "pH used by LUNA when adding hydrogens.", "es": "pH usado por LUNA al agregar hidrógenos."},
    "Fingerprints de interação (IFP)": {"en": "Interaction fingerprints (IFP)", "es": "Fingerprints de interacción (IFP)"},
    (
        "Fingerprints transformam as interações em uma tabela comparável entre ligantes. "
        "É a saída mais útil para análise posterior e costuma ser a opção principal."
    ): {
        "en": "Fingerprints transform interactions into a ligand-comparable table. This is usually the most useful output for downstream analysis and often the main option.",
        "es": "Los fingerprints transforman las interacciones en una tabla comparable entre ligandos. Es la salida más útil para análisis posteriores y suele ser la opción principal.",
    },
    "Todos (H + E + F)": {"en": "All (H + E + F)", "es": "Todos (H + E + F)"},
    "Tipo:": {"en": "Type:", "es": "Tipo:"},
    "Saída IFP:": {"en": "IFP output:", "es": "Salida IFP:"},
    "Rótulos para importância de fingerprints": {"en": "Labels for fingerprint importance", "es": "Rótulos para importancia de fingerprints"},
    (
        "Se marcado, a aba FP análises usa este CSV para treinar a importância das features. "
        "Informe o arquivo, a coluna do ID do ligante e a coluna que contém os rótulos/classes."
    ): {
        "en": "When checked, the FP analyses tab uses this CSV to train feature importance. Provide the file, the ligand ID column and the column containing labels/classes.",
        "es": "Cuando está marcado, la pestaña de análisis FP usa este CSV para entrenar la importancia de las features. Indica el archivo, la columna del ID del ligando y la columna con rótulos/clases.",
    },
    "Arquivo CSV:": {"en": "CSV file:", "es": "Archivo CSV:"},
    "Coluna do ligand_id:": {"en": "ligand_id column:", "es": "Columna ligand_id:"},
    "Coluna de rótulo:": {"en": "Label column:", "es": "Columna de rótulo:"},
    "Tarefa:": {"en": "Task:", "es": "Tarea:"},
    "Regressão": {"en": "Regression", "es": "Regresión"},
    "Classificação": {"en": "Classification", "es": "Clasificación"},
    "Aplicar Otsu tambem a interacoes/residuos": {"en": "Apply Otsu also to interactions/residues", "es": "Aplicar Otsu también a interacciones/residuos"},
    "Matriz de similaridade (Tanimoto)": {"en": "Similarity matrix (Tanimoto)", "es": "Matriz de similaridad (Tanimoto)"},
    (
        "A matriz de similaridade ajuda a ver quais ligantes se comportam de forma parecida. "
        "Ela também alimenta os gráficos e clusters da aba de resultados."
    ): {
        "en": "The similarity matrix helps identify ligands with similar interaction behavior. It also feeds charts and clusters in the results tab.",
        "es": "La matriz de similaridad ayuda a ver qué ligandos se comportan de forma parecida. También alimenta los gráficos y clusters de la pestaña de resultados.",
    },
    "Saída:": {"en": "Output:", "es": "Salida:"},
    "Exportar sessões PyMOL (.pse)": {"en": "Export PyMOL sessions (.pse)", "es": "Exportar sesiones PyMOL (.pse)"},
    "As sessões PyMOL permitem abrir a proteína e os ligantes já com as interações destacadas visualmente.": {
        "en": "PyMOL sessions open the protein and ligands with interactions already highlighted visually.",
        "es": "Las sesiones PyMOL permiten abrir la proteína y los ligandos con las interacciones ya resaltadas visualmente.",
    },
    "Pasta de saída:": {"en": "Output folder:", "es": "Carpeta de salida:"},
    "Filtrar por binding modes (.cfg)": {"en": "Filter by binding modes (.cfg)", "es": "Filtrar por binding modes (.cfg)"},
    "Use este filtro quando você já tem regras de binding modes e quer limitar quais poses ou interações entram na análise.": {
        "en": "Use this filter when you already have binding-mode rules and want to limit which poses or interactions enter the analysis.",
        "es": "Usa este filtro cuando ya tienes reglas de binding modes y quieres limitar qué poses o interacciones entran en el análisis.",
    },
    "Editor visual": {"en": "Visual editor", "es": "Editor visual"},
    "Filtrar PSE por tipo de interação (opcional)": {"en": "Filter PSE by interaction type (optional)", "es": "Filtrar PSE por tipo de interacción (opcional)"},
    "Se marcado, LUNA gera sessões PyMOL contendo apenas os tipos selecionados.": {
        "en": "When checked, LUNA generates PyMOL sessions containing only the selected types.",
        "es": "Cuando está marcado, LUNA genera sesiones PyMOL que contienen solo los tipos seleccionados.",
    },
    "Opções avançadas — DefaultInteractionConfig + InteractionCalculator": {"en": "Advanced options - DefaultInteractionConfig + InteractionCalculator", "es": "Opciones avanzadas - DefaultInteractionConfig + InteractionCalculator"},
    (
        "Use este quadro apenas para limitar globalmente as distâncias máximas do cálculo "
        "ou ajustar flags do InteractionCalculator. Se quiser mudar parâmetros específicos "
        "de cada interação, use o quadro do arquivo completo .cfg logo abaixo."
    ): {
        "en": "Use this section only to globally limit maximum calculation distances or adjust InteractionCalculator flags. To change interaction-specific parameters, use the full .cfg file section below.",
        "es": "Usa este cuadro solo para limitar globalmente las distancias máximas del cálculo o ajustar flags de InteractionCalculator. Para cambiar parámetros específicos de cada interacción, usa el cuadro del archivo .cfg completo abajo.",
    },
    "Distância máxima global do cálculo de interações": {"en": "Global maximum distance for interaction calculation", "es": "Distancia máxima global del cálculo de interacciones"},
    (
        "Se definido, qualquer distância máxima do LUNA acima deste valor será limitada "
        "para este projeto. Distâncias padrão menores do que o valor informado serão mantidas."
    ): {
        "en": "If set, any LUNA maximum distance above this value is capped for this project. Default distances below the informed value are kept.",
        "es": "Si se define, cualquier distancia máxima de LUNA por encima de este valor será limitada para este proyecto. Distancias por defecto menores que el valor informado se mantienen.",
    },
    "Distância máxima global (A):": {"en": "Global maximum distance (A):", "es": "Distancia máxima global (A):"},
    "(sem limite global)": {"en": "(no global limit)", "es": "(sin límite global)"},
    "<b>Flags do InteractionCalculator</b>": {"en": "<b>InteractionCalculator flags</b>", "es": "<b>Flags de InteractionCalculator</b>"},
    "Arquivo completo de interações (.cfg)": {"en": "Full interaction file (.cfg)", "es": "Archivo completo de interacciones (.cfg)"},
    (
        "Use este quadro quando quiser controlar, em um único arquivo .cfg, "
        "quais interações o LUNA deve calcular e todos os parâmetros de cada uma. "
        "Os ajustes do quadro avançado acima continuam sendo aplicados por cima deste arquivo."
    ): {
        "en": "Use this section when you want one .cfg file to control which interactions LUNA should calculate and all parameters for each one. The advanced settings above are still applied on top of this file.",
        "es": "Usa este cuadro cuando quieras controlar, en un único archivo .cfg, qué interacciones debe calcular LUNA y todos los parámetros de cada una. Los ajustes avanzados anteriores siguen aplicándose sobre este archivo.",
    },
    "(opcional: config_interacoes.cfg)": {"en": "(optional: interaction_config.cfg)", "es": "(opcional: config_interacciones.cfg)"},
    "<a href='open_default_interaction_config'>Abrir exemplo padrão utilizado pelo LUNA</a>": {"en": "<a href='open_default_interaction_config'>Open the default example used by LUNA</a>", "es": "<a href='open_default_interaction_config'>Abrir el ejemplo por defecto usado por LUNA</a>"},
    "0.05 Flexivel": {"en": "0.05 Flexible", "es": "0.05 Flexible"},
    "0.02 Medio": {"en": "0.02 Medium", "es": "0.02 Medio"},
    "0.01 Conservador": {"en": "0.01 Conservative", "es": "0.01 Conservador"},

    # Run tab.
    (
        "Esta aba executa o projeto com as opções já escolhidas. O número de núcleos controla "
        "o paralelismo, e o log abaixo mostra exatamente o que o LUNA está fazendo."
    ): {
        "en": "This tab runs the project with the selected options. The number of cores controls parallelism, and the log below shows exactly what LUNA is doing.",
        "es": "Esta pestaña ejecuta el proyecto con las opciones ya elegidas. El número de núcleos controla el paralelismo y el log muestra exactamente qué está haciendo LUNA.",
    },
    "Opções de execução": {"en": "Run options", "es": "Opciones de ejecución"},
    "Núcleos (--nproc):": {"en": "Cores (--nproc):", "es": "Núcleos (--nproc):"},
    "Sobrescrever projeto existente (--overwrite)": {"en": "Overwrite existing project (--overwrite)", "es": "Sobrescribir proyecto existente (--overwrite)"},
    "▶ Executar LUNA": {"en": "▶ Run LUNA", "es": "▶ Ejecutar LUNA"},
    "■ Cancelar": {"en": "■ Cancel", "es": "■ Cancelar"},
    "Aguardando execução": {"en": "Waiting to run", "es": "Esperando ejecución"},
    "Mostra o progresso informado pelo LUNA durante a etapa atual.": {"en": "Shows the progress reported by LUNA during the current stage.", "es": "Muestra el progreso informado por LUNA durante la etapa actual."},
    "LUNA não pronto": {"en": "LUNA not ready", "es": "LUNA no está listo"},
    "Configuração inválida": {"en": "Invalid configuration", "es": "Configuración inválida"},
    "Alerta de hidrogênios (Add_H)": {"en": "Hydrogen warning (Add_H)", "es": "Alerta de hidrógenos (Add_H)"},
    "Concluído": {"en": "Completed", "es": "Concluido"},
    "0% - iniciando": {"en": "0% - starting", "es": "0% - iniciando"},
    "100% - concluído": {"en": "100% - completed", "es": "100% - concluido"},
    "Falhou": {"en": "Failed", "es": "Falló"},

    # Results tab and sub-tabs.
    (
        "Use esta aba para revisar o que o LUNA gerou. Você pode carregar um workdir já existente, "
        "visualizar tabelas e gráficos, exportar figuras e agrupar ligantes por similaridade."
    ): {
        "en": "Use this tab to review what LUNA generated. You can load an existing workdir, inspect tables and charts, export figures and cluster ligands by similarity.",
        "es": "Usa esta pestaña para revisar lo que LUNA generó. Puedes cargar un workdir existente, visualizar tablas y gráficos, exportar figuras y agrupar ligandos por similaridad.",
    },
    "Workdir:": {"en": "Workdir:", "es": "Workdir:"},
    "(usa o workdir do projeto atual)": {"en": "(uses the current project workdir)", "es": "(usa el workdir del proyecto actual)"},
    "Carregar resultados": {"en": "Load results", "es": "Cargar resultados"},
    "Exportar gráfico atual...": {"en": "Export current chart...", "es": "Exportar gráfico actual..."},
    "Exportar relatório HTML": {"en": "Export HTML report", "es": "Exportar reporte HTML"},
    "Gerar relatório PDF": {"en": "Generate PDF report", "es": "Generar reporte PDF"},
    "Fingerprints": {"en": "Fingerprints", "es": "Fingerprints"},
    "Pré-visualização do arquivo `ifp.csv`. Cada linha representa um ligante e cada coluna descreve partes do fingerprint gerado.": {
        "en": "Preview of the `ifp.csv` file. Each row represents a ligand and each column describes parts of the generated fingerprint.",
        "es": "Vista previa del archivo `ifp.csv`. Cada fila representa un ligando y cada columna describe partes del fingerprint generado.",
    },
    "Linhas a exibir:": {"en": "Rows to show:", "es": "Filas a mostrar:"},
    "Matriz de similaridade": {"en": "Similarity matrix", "es": "Matriz de similaridad"},
    "Mapa de calor da similaridade entre ligantes. Valores mais altos indicam fingerprints mais parecidos.": {
        "en": "Heatmap of ligand similarity. Higher values indicate more similar fingerprints.",
        "es": "Mapa de calor de la similaridad entre ligandos. Valores más altos indican fingerprints más parecidos.",
    },
    "Estatísticas": {"en": "Statistics", "es": "Estadísticas"},
    "Resume quantas interações de cada tipo aparecem no conjunto analisado. É uma visão global do perfil químico observado.": {
        "en": "Summarizes how many interactions of each type appear in the analyzed set. It is a global view of the observed chemical profile.",
        "es": "Resume cuántas interacciones de cada tipo aparecen en el conjunto analizado. Es una visión global del perfil químico observado.",
    },
    "Calcular estatísticas (usa luna-env)": {"en": "Calculate statistics (uses luna-env)", "es": "Calcular estadísticas (usa luna-env)"},
    "Heatmap por tipo": {"en": "Heatmap by type", "es": "Heatmap por tipo"},
    "Mostra, para um tipo de interação escolhido, quais resíduos da proteína aparecem associados aos ligantes.": {
        "en": "Shows, for the selected interaction type, which protein residues are associated with ligands.",
        "es": "Muestra, para un tipo de interacción elegido, qué residuos de la proteína aparecen asociados a los ligandos.",
    },
    "Calcular heatmap (usa luna-env)": {"en": "Calculate heatmap (uses luna-env)", "es": "Calcular heatmap (usa luna-env)"},
    "Método:": {"en": "Method:", "es": "Método:"},
    "Clusters:": {"en": "Clusters:", "es": "Clusters:"},
    "Atualizar clusters": {"en": "Update clusters", "es": "Actualizar clusters"},
    "Agrupa ligantes com base na matriz de similaridade. O dendrograma mostra a relação hierárquica, e a tabela lista o cluster atribuído a cada ligante.": {
        "en": "Clusters ligands based on the similarity matrix. The dendrogram shows the hierarchical relationship, and the table lists the cluster assigned to each ligand.",
        "es": "Agrupa ligandos con base en la matriz de similaridad. El dendrograma muestra la relación jerárquica y la tabla lista el cluster atribuido a cada ligando.",
    },
    "Exportar clusters CSV": {"en": "Export clusters CSV", "es": "Exportar clusters CSV"},
    "Abrir no PyMOL": {"en": "Open in PyMOL", "es": "Abrir en PyMOL"},
    "Lista as sessões PyMOL exportadas pelo LUNA. Abra um arquivo para inspecionar visualmente interações e poses.": {
        "en": "Lists the PyMOL sessions exported by LUNA. Open a file to visually inspect interactions and poses.",
        "es": "Lista las sesiones PyMOL exportadas por LUNA. Abre un archivo para inspeccionar visualmente interacciones y poses.",
    },
    "Sessões PyMOL": {"en": "PyMOL sessions", "es": "Sesiones PyMOL"},
    "Estatisticas": {"en": "Statistics", "es": "Estadísticas"},
    "Heatmap completo ligantes x residuos": {"en": "Complete ligand x residue heatmap", "es": "Heatmap completo ligandos x residuos"},
    "Sessoes em PyMOL": {"en": "PyMOL sessions", "es": "Sesiones en PyMOL"},
    "FP analises": {"en": "FP analyses", "es": "Análisis FP"},
    "FP sessao em PyMOL": {"en": "FP PyMOL session", "es": "Sesión FP en PyMOL"},
    "Filtragens dinâmicas por binding modes": {"en": "Dynamic filtering by binding modes", "es": "Filtrados dinámicos por binding modes"},
    (
        "Gera novas sessões PyMOL em uma subpasta separada usando regras .cfg de binding modes. "
        "As sessões originais não são alteradas; você pode criar, apagar e recriar filtros quantas vezes precisar."
    ): {
        "en": "Generates new PyMOL sessions in a separate subfolder using .cfg binding-mode rules. Original sessions are not changed; you can create, delete and recreate filters as often as needed.",
        "es": "Genera nuevas sesiones PyMOL en una subcarpeta separada usando reglas .cfg de binding modes. Las sesiones originales no se alteran; puedes crear, borrar y recrear filtros cuantas veces necesites.",
    },
    "Regras .cfg:": {"en": ".cfg rules:", "es": "Reglas .cfg:"},
    "Nome do filtro:": {"en": "Filter name:", "es": "Nombre del filtro:"},
    "Gerar filtragem": {"en": "Generate filter", "es": "Generar filtrado"},
    "Apagar filtragem selecionada": {"en": "Delete selected filter", "es": "Borrar filtrado seleccionado"},
    "Arquivos .pse da filtragem selecionada:": {"en": ".pse files from the selected filter:", "es": "Archivos .pse del filtrado seleccionado:"},
    "Selecione uma filtragem para listar os arquivos.": {"en": "Select a filter to list its files.", "es": "Selecciona un filtrado para listar los archivos."},
    "Visao:": {"en": "View:", "es": "Vista:"},
    "Visão:": {"en": "View:", "es": "Vista:"},
    "Totais do projeto": {"en": "Project totals", "es": "Totales del proyecto"},
    "Todos os frames": {"en": "All frames", "es": "Todos los frames"},
    "Todos os ligantes": {"en": "All ligands", "es": "Todos los ligandos"},
    "Atualizar heatmaps": {"en": "Update heatmaps", "es": "Actualizar heatmaps"},
    (
        "Mostra, em cada par ligante x residuo, todas as interacoes presentes em faixas coloridas. "
        "As cores indicam a classe da interacao, nao a contagem."
    ): {
        "en": "Shows all interactions present in each ligand x residue pair as colored bands. Colors indicate interaction class, not count.",
        "es": "Muestra, en cada par ligando x residuo, todas las interacciones presentes como franjas de color. Los colores indican la clase de interacción, no el conteo.",
    },
    "Carregar analises de FP": {"en": "Load FP analyses", "es": "Cargar análisis FP"},
    (
        "Classifica os fingerprints com base nos shells reais do LUNA, "
        "atribui classes confiaveis via limiar por z-score e mostra as features mais relevantes."
    ): {
        "en": "Classifies fingerprints based on real LUNA shells, assigns reliable classes through a z-score threshold and shows the most relevant features.",
        "es": "Clasifica los fingerprints con base en los shells reales de LUNA, asigna clases confiables mediante umbral de z-score y muestra las features más relevantes.",
    },
    "Algoritmo:": {"en": "Algorithm:", "es": "Algoritmo:"},
    "Corte p-value:": {"en": "p-value cutoff:", "es": "Corte p-value:"},
    "Classes importantes": {"en": "Important classes", "es": "Clases importantes"},
    "Frequencia por classe": {"en": "Frequency by class", "es": "Frecuencia por clase"},
    "Cobertura e importancia": {"en": "Coverage and importance", "es": "Cobertura e importancia"},
    "Heatmap importância": {"en": "Importance heatmap", "es": "Heatmap de importancia"},
    "Frequencia interacoes": {"en": "Interaction frequency", "es": "Frecuencia de interacciones"},
    "Interações prevalentes": {"en": "Prevalent interactions", "es": "Interacciones prevalentes"},
    "Heatmap interações": {"en": "Interaction heatmap", "es": "Heatmap de interacciones"},
    "Fingerprint:": {"en": "Fingerprint:", "es": "Fingerprint:"},
    "Molecula:": {"en": "Molecule:", "es": "Molécula:"},
    "Gerar sessao": {"en": "Generate session", "es": "Generar sesión"},
    (
        "Gera uma sessao PyMOL para o fingerprint escolhido, recuperando os shells "
        "que originaram aquele bit para um ligante especifico."
    ): {
        "en": "Generates a PyMOL session for the selected fingerprint, recovering the shells that originated that bit for a specific ligand.",
        "es": "Genera una sesión PyMOL para el fingerprint elegido, recuperando los shells que originaron ese bit para un ligando específico.",
    },
    "Abrir sessao": {"en": "Open session", "es": "Abrir sesión"},

    # Tables.
    "Feature": {"en": "Feature", "es": "Feature"},
    "Moleculas": {"en": "Molecules", "es": "Moléculas"},
    "Cobertura (%)": {"en": "Coverage (%)", "es": "Cobertura (%)"},
    "Classe prevalente (%)": {"en": "Prevalent class (%)", "es": "Clase prevalente (%)"},
    "Z-score classe": {"en": "Class z-score", "es": "Z-score clase"},
    "Classe atribuida": {"en": "Assigned class", "es": "Clase asignada"},
    "Importancia": {"en": "Importance", "es": "Importancia"},
    "Z-score Importance": {"en": "Importance z-score", "es": "Z-score importancia"},
    "Colisoes": {"en": "Collisions", "es": "Colisiones"},
    "Perfil da base": {"en": "Dataset profile", "es": "Perfil de la base"},
    "Ligante": {"en": "Ligand", "es": "Ligando"},
    "Ordem": {"en": "Order", "es": "Orden"},

    # Binding mode editor and preparation dialog.
    "Editor de Binding Modes": {"en": "Binding Modes Editor", "es": "Editor de Binding Modes"},
    "Tipo de interação": {"en": "Interaction type", "es": "Tipo de interacción"},
    "accept_only (lista)": {"en": "accept_only (list)", "es": "accept_only (lista)"},
    "+ Adicionar tipo": {"en": "+ Add type", "es": "+ Agregar tipo"},
    "− Remover linha": {"en": "- Remove row", "es": "- Remover fila"},
    "Carregar .cfg": {"en": "Load .cfg", "es": "Cargar .cfg"},
    "Salvar .cfg": {"en": "Save .cfg", "es": "Guardar .cfg"},
    "Preparar arquivos de docking": {"en": "Prepare docking files", "es": "Preparar archivos de docking"},
    "Pasta de origem:": {"en": "Source folder:", "es": "Carpeta de origen:"},
    "Pasta de saída:": {"en": "Output folder:", "es": "Carpeta de salida:"},
    "Último átomo da proteína": {"en": "Last protein atom", "es": "Último átomo de la proteína"},
    "Detectar automaticamente": {"en": "Detect automatically", "es": "Detectar automáticamente"},
    "Informar manualmente": {"en": "Enter manually", "es": "Informar manualmente"},
    "Detectar agora": {"en": "Detect now", "es": "Detectar ahora"},
    "(selecione a pasta de origem primeiro)": {"en": "(select the source folder first)", "es": "(selecciona primero la carpeta de origen)"},
    "Último átomo da proteína:": {"en": "Last protein atom:", "es": "Último átomo de la proteína:"},
    "Usar as pastas geradas como entradas do projeto ao fechar": {"en": "Use generated folders as project inputs when closing", "es": "Usar las carpetas generadas como entradas del proyecto al cerrar"},
    "Executar preparação": {"en": "Run preparation", "es": "Ejecutar preparación"},
    "Aguardando execução.": {"en": "Waiting to run.", "es": "Esperando ejecución."},

    # Plot text.
    "Similaridade de Tanimoto": {"en": "Tanimoto similarity", "es": "Similaridad de Tanimoto"},
    "Sem dados de interação": {"en": "No interaction data", "es": "Sin datos de interacción"},
    "Sem interações para esta visão.": {"en": "No interactions for this view.", "es": "Sin interacciones para esta vista."},
    "Aminoácidos": {"en": "Amino acids", "es": "Aminoácidos"},
    "Residuos": {"en": "Residues", "es": "Residuos"},
    "Resíduos": {"en": "Residues", "es": "Residuos"},
    "Ligantes": {"en": "Ligands", "es": "Ligandos"},
    "Frames": {"en": "Frames", "es": "Frames"},
    "% de ligantes": {"en": "% of ligands", "es": "% de ligandos"},
    "% de frames": {"en": "% of frames", "es": "% de frames"},
    "Total (todas as entradas)": {"en": "Total (all entries)", "es": "Total (todas las entradas)"},
    "Total neste ligante": {"en": "Total for this ligand", "es": "Total en este ligando"},
    "Contagem por tipo de interação": {"en": "Count by interaction type", "es": "Conteo por tipo de interacción"},
    "Interações por tipo": {"en": "Interactions by type", "es": "Interacciones por tipo"},
    "Tipos de interacao por par ligante x residuo": {"en": "Interaction types by ligand x residue pair", "es": "Tipos de interacción por par ligando x residuo"},
    "Clique para ocultar/mostrar": {"en": "Click to hide/show", "es": "Click para ocultar/mostrar"},
    "Sem matriz de residuos disponivel": {"en": "No residue matrix available", "es": "Sin matriz de residuos disponible"},
    "Sem matriz completa disponivel": {"en": "No complete matrix available", "es": "Sin matriz completa disponible"},
    "Sem dados de interacao para o heatmap completo": {"en": "No interaction data for the complete heatmap", "es": "Sin datos de interacción para el heatmap completo"},
    "Sem features importantes para resumir.": {"en": "No important features to summarize.", "es": "Sin features importantes para resumir."},
    "Sem features importantes para plotar.": {"en": "No important features to plot.", "es": "Sin features importantes para graficar."},
    "Sem features importantes para gerar o heatmap.": {"en": "No important features to generate the heatmap.", "es": "Sin features importantes para generar el heatmap."},
    "Sem interacoes prevalentes confiaveis nas features importantes.": {"en": "No reliable prevalent interactions among important features.", "es": "Sin interacciones prevalentes confiables en las features importantes."},
    "Todas as interacoes estao ocultas pela legenda.": {"en": "All interactions are hidden by the legend.", "es": "Todas las interacciones están ocultas por la leyenda."},
    "Sem interacoes prevalentes para gerar o heatmap.": {"en": "No prevalent interactions to generate the heatmap.", "es": "Sin interacciones prevalentes para generar el heatmap."},
    "% Important features": {"en": "% Important features", "es": "% Features importantes"},
    "Distribuicao das classes entre as features mais importantes": {"en": "Class distribution among the most important features", "es": "Distribución de clases entre las features más importantes"},
    "Assignment frequency of each class (%)": {"en": "Assignment frequency of each class (%)", "es": "Frecuencia de asignación de cada clase (%)"},
    "Feature id": {"en": "Feature id", "es": "ID de feature"},
    "Frequencia de atribuicao de classes nas features importantes": {"en": "Class assignment frequency in important features", "es": "Frecuencia de asignación de clases en features importantes"},
    "% Fingerprints containing the feature": {"en": "% Fingerprints containing the feature", "es": "% Fingerprints que contienen la feature"},
    "Cobertura das features importantes e importancia do modelo": {"en": "Important feature coverage and model importance", "es": "Cobertura de features importantes e importancia del modelo"},
    "Mapa de presenca das features importantes por classe": {"en": "Presence map of important features by class", "es": "Mapa de presencia de features importantes por clase"},
    "Assignment frequency of each interaction (%)": {"en": "Assignment frequency of each interaction (%)", "es": "Frecuencia de asignación de cada interacción (%)"},
    "Frequencia de atribuicao da interacao prevalente nas features importantes": {"en": "Prevalent interaction assignment frequency in important features", "es": "Frecuencia de asignación de la interacción prevalente en features importantes"},
    "Numero de ligantes": {"en": "Number of ligands", "es": "Número de ligandos"},
    "Interacao e residuo prevalentes nas features importantes": {"en": "Prevalent interaction and residue in important features", "es": "Interacción y residuo prevalentes en features importantes"},
    "Interacoes prevalentes das features importantes por ligante": {"en": "Prevalent interactions of important features by ligand", "es": "Interacciones prevalentes de features importantes por ligando"},
    "Clustering hierárquico": {"en": "Hierarchical clustering", "es": "Clustering jerárquico"},
    "Distância": {"en": "Distance", "es": "Distancia"},
    "Matriz reordenada por cluster": {"en": "Matrix reordered by cluster", "es": "Matriz reordenada por cluster"},

    # Report text.
    "Relatório LUNA GUI": {"en": "LUNA GUI report", "es": "Reporte LUNA GUI"},
    "Como interpretar as análises": {"en": "How to interpret the analyses", "es": "Cómo interpretar los análisis"},
    "Distribuição de interações": {"en": "Interaction distribution", "es": "Distribución de interacciones"},
    "Resumo das FP análises": {"en": "FP analyses summary", "es": "Resumen de análisis FP"},
    "Atribuição de clusters": {"en": "Cluster assignment", "es": "Asignación de clusters"},
    "Interações mais frequentes": {"en": "Most frequent interactions", "es": "Interacciones más frecuentes"},
    "Resíduos mais frequentes": {"en": "Most frequent residues", "es": "Residuos más frecuentes"},
    "Sem interações contabilizadas.": {"en": "No counted interactions.", "es": "Sin interacciones contabilizadas."},
    "Sem resíduos contabilizados.": {"en": "No counted residues.", "es": "Sin residuos contabilizados."},
    "Proteína": {"en": "Protein", "es": "Proteína"},
    "Incluir águas": {"en": "Include waters", "es": "Incluir aguas"},
    "Adicionar H": {"en": "Add H", "es": "Agregar H"},
    "Fingerprint de contagem": {"en": "Count fingerprint", "es": "Fingerprint de conteo"},
    "Sessões PyMOL": {"en": "PyMOL sessions", "es": "Sesiones PyMOL"},
    "Filtro binding modes": {"en": "Binding modes filter", "es": "Filtro binding modes"},
    "Config. interações": {"en": "Interaction config", "es": "Config. interacciones"},
    "Distância máxima global": {"en": "Global maximum distance", "es": "Distancia máxima global"},
    "Rótulos FP": {"en": "FP labels", "es": "Rótulos FP"},
    "Tarefa FP": {"en": "FP task", "es": "Tarea FP"},
    "Núcleos": {"en": "Cores", "es": "Núcleos"},
    (
        "Este relatório resume os resultados carregados na aba 5.Resultados, os parâmetros usados no projeto e uma leitura guiada dos gráficos. Use-o como material de triagem: padrões fortes indicam hipóteses para inspeção estrutural, não uma conclusão automática de afinidade."
    ): {
        "en": "This report summarizes the results loaded in tab 5.Results, the project parameters and a guided reading of the charts. Use it as screening material: strong patterns indicate hypotheses for structural inspection, not an automatic affinity conclusion.",
        "es": "Este reporte resume los resultados cargados en la pestaña 5.Resultados, los parámetros usados en el proyecto y una lectura guiada de los gráficos. Úsalo como material de triaje: patrones fuertes indican hipótesis para inspección estructural, no una conclusión automática de afinidad.",
    },
    (
        "Estatísticas de interação: contam quantas vezes cada tipo de contato foi observado. Barras altas indicam mecanismos recorrentes, como ligações de hidrogênio, contatos hidrofóbicos ou interações iônicas. Compare tipos dominantes com resíduos frequentes para separar padrões químicos reais de ruído de pose."
    ): {
        "en": "Interaction statistics: count how many times each contact type was observed. Tall bars indicate recurring mechanisms, such as hydrogen bonds, hydrophobic contacts or ionic interactions. Compare dominant types with frequent residues to separate real chemical patterns from pose noise.",
        "es": "Estadísticas de interacción: cuentan cuántas veces se observó cada tipo de contacto. Barras altas indican mecanismos recurrentes, como puentes de hidrógeno, contactos hidrofóbicos o interacciones iónicas. Compara tipos dominantes con residuos frecuentes para separar patrones químicos reales del ruido de pose.",
    },
    (
        "Heatmap por tipo: cruza ligantes e resíduos para uma interação escolhida. Células mais intensas indicam mais ocorrências daquele contato. Colunas densas sugerem resíduos-chave; linhas densas sugerem ligantes com muitos contatos daquele tipo."
    ): {
        "en": "Heatmap by type: crosses ligands and residues for a selected interaction. More intense cells indicate more occurrences of that contact. Dense columns suggest key residues; dense rows suggest ligands with many contacts of that type.",
        "es": "Heatmap por tipo: cruza ligandos y residuos para una interacción elegida. Celdas más intensas indican más ocurrencias de ese contacto. Columnas densas sugieren residuos clave; filas densas sugieren ligandos con muchos contactos de ese tipo.",
    },
    (
        "Matriz de similaridade: compara ligantes pelos fingerprints de interação. Valores próximos de 1 indicam perfis de interação semelhantes; valores baixos indicam modos de interação distintos, mesmo quando as moléculas parecem estruturalmente parecidas."
    ): {
        "en": "Similarity matrix: compares ligands by interaction fingerprints. Values close to 1 indicate similar interaction profiles; low values indicate distinct interaction modes, even when molecules look structurally similar.",
        "es": "Matriz de similaridad: compara ligandos por fingerprints de interacción. Valores cercanos a 1 indican perfiles de interacción semejantes; valores bajos indican modos de interacción distintos, incluso cuando las moléculas parecen estructuralmente parecidas.",
    },
    (
        "Cada barra representa a contagem de uma classe de interação. Use barras dominantes para identificar forças químicas recorrentes e barras raras para procurar contatos específicos que podem diferenciar poucos ligantes."
    ): {
        "en": "Each bar represents the count of an interaction class. Use dominant bars to identify recurring chemical forces and rare bars to look for specific contacts that may differentiate a few ligands.",
        "es": "Cada barra representa el conteo de una clase de interacción. Usa barras dominantes para identificar fuerzas químicas recurrentes y barras raras para buscar contactos específicos que pueden diferenciar pocos ligandos.",
    },
    (
        "Cada célula compara dois ligantes pelo fingerprint de interação. Tons mais intensos indicam maior similaridade; blocos ao longo da diagonal sugerem famílias com modos de interação semelhantes."
    ): {
        "en": "Each cell compares two ligands by interaction fingerprint. More intense tones indicate higher similarity; blocks along the diagonal suggest families with similar interaction modes.",
        "es": "Cada celda compara dos ligandos por fingerprint de interacción. Tonos más intensos indican mayor similaridad; bloques a lo largo de la diagonal sugieren familias con modos de interacción semejantes.",
    },

    # Standard interaction labels used in plots.  Data keys remain unchanged.
    "Hydrogen bond": {"pt": "Ligação de hidrogênio", "en": "Hydrogen bond", "es": "Puente de hidrógeno"},
    "Weak hydrogen bond": {"pt": "Ligação de hidrogênio fraca", "en": "Weak hydrogen bond", "es": "Puente de hidrógeno débil"},
    "Hydrophobic": {"pt": "Hidrofóbica", "en": "Hydrophobic", "es": "Hidrofóbica"},
    "Halogen bond": {"pt": "Ligação de halogênio", "en": "Halogen bond", "es": "Enlace de halógeno"},
    "Chalcogen bond": {"pt": "Ligação de calcogênio", "en": "Chalcogen bond", "es": "Enlace de calcógeno"},
    "Ionic": {"pt": "Iônica", "en": "Ionic", "es": "Iónica"},
    "Salt bridge": {"pt": "Ponte salina", "en": "Salt bridge", "es": "Puente salino"},
    "Cation-pi": {"pt": "Cátion-pi", "en": "Cation-pi", "es": "Catión-pi"},
    "Pi-stacking": {"pt": "Empilhamento pi", "en": "Pi-stacking", "es": "Apilamiento pi"},
    "Face-to-face": {"pt": "Face-a-face", "en": "Face-to-face", "es": "Cara a cara"},
    "Edge-to-face": {"pt": "Borda-a-face", "en": "Edge-to-face", "es": "Borde a cara"},
    "Water-bridged hydrogen bond": {"pt": "Ponte de hidrogênio mediada por água", "en": "Water-bridged hydrogen bond", "es": "Puente de hidrógeno mediado por agua"},
    "Disulfide bond": {"pt": "Ponte dissulfeto", "en": "Disulfide bond", "es": "Puente disulfuro"},
    "Metal coordination": {"pt": "Coordenação metálica", "en": "Metal coordination", "es": "Coordinación metálica"},
    "Van der Waals": {"pt": "Van der Waals", "en": "Van der Waals", "es": "Van der Waals"},
    "Proximal": {"pt": "Proximal", "en": "Proximal", "es": "Proximal"},
}

_REVERSE_TRANSLATIONS: dict[str, str] = {}
for source, values in TRANSLATIONS.items():
    _REVERSE_TRANSLATIONS[source] = source
    for translated in values.values():
        _REVERSE_TRANSLATIONS[translated] = source


def language() -> str:
    return _current_language


def set_language(lang: str) -> str:
    global _current_language
    lang = str(lang or "pt").lower()
    if lang not in LANGUAGES:
        lang = "pt"
    _current_language = lang
    return lang


def source_text(value: Any) -> str:
    text = "" if value is None else str(value)
    return _REVERSE_TRANSLATIONS.get(text, text)


def t(value: Any, lang: str | None = None) -> str:
    text = "" if value is None else str(value)
    source = source_text(text)
    target = lang or _current_language
    translated = TRANSLATIONS.get(source, {}).get(target)
    if translated is not None:
        return translated
    if target == "pt":
        return source
    return _translate_patterns(source, target)


def _translate_patterns(text: str, lang: str) -> str:
    patterns: list[tuple[re.Pattern[str], Callable[[re.Match[str]], str]]] = [
        (
            re.compile(r"^(\d+) de (\d+) ligantes selecionados$"),
            lambda m: (
                f"{m.group(1)} of {m.group(2)} selected ligands"
                if lang == "en"
                else f"{m.group(1)} de {m.group(2)} ligandos seleccionados"
            ),
        ),
        (
            re.compile(r"^Águas detectadas nos inputs: (\d+)$"),
            lambda m: (
                f"Waters detected in inputs: {m.group(1)}"
                if lang == "en"
                else f"Aguas detectadas en los inputs: {m.group(1)}"
            ),
        ),
        (
            re.compile(r"^(\d+) entradas processadas$"),
            lambda m: (
                f"{m.group(1)} processed entries"
                if lang == "en"
                else f"{m.group(1)} entradas procesadas"
            ),
        ),
        (
            re.compile(r"^Gerado em (.+)\.$"),
            lambda m: (
                f"Generated on {m.group(1)}."
                if lang == "en"
                else f"Generado el {m.group(1)}."
            ),
        ),
        (
            re.compile(r"^Entradas processadas: (.+)\.$"),
            lambda m: (
                f"Processed entries: {m.group(1)}."
                if lang == "en"
                else f"Entradas procesadas: {m.group(1)}."
            ),
        ),
        (
            re.compile(r"^(\d+) entradas · (\d+) tipos$"),
            lambda m: (
                f"{m.group(1)} entries · {m.group(2)} types"
                if lang == "en"
                else f"{m.group(1)} entradas · {m.group(2)} tipos"
            ),
        ),
        (
            re.compile(r"^Frame: (.+)$"),
            lambda m: f"Frame: {m.group(1)}",
        ),
        (
            re.compile(r"^Ligante: (.+)$"),
            lambda m: (
                f"Ligand: {m.group(1)}" if lang == "en" else f"Ligando: {m.group(1)}"
            ),
        ),
        (
            re.compile(r"^Interações por tipo — (.+)$"),
            lambda m: (
                f"Interactions by type - {m.group(1)}"
                if lang == "en"
                else f"Interacciones por tipo - {m.group(1)}"
            ),
        ),
        (
            re.compile(r"^Resíduos × ligantes — (.+)$"),
            lambda m: (
                f"Residues x ligands - {m.group(1)}"
                if lang == "en"
                else f"Residuos x ligandos - {m.group(1)}"
            ),
        ),
        (
            re.compile(r"^Sem ocorrencias de '(.+)'$"),
            lambda m: (
                f"No occurrences of '{m.group(1)}'"
                if lang == "en"
                else f"Sin ocurrencias de '{m.group(1)}'"
            ),
        ),
        (
            re.compile(r"^(\d+) ligantes em ordem hierárquica$"),
            lambda m: (
                f"{m.group(1)} ligands in hierarchical order"
                if lang == "en"
                else f"{m.group(1)} ligandos en orden jerárquico"
            ),
        ),
        (
            re.compile(r"^(\d+) arquivos \.pse nesta filtragem\.$"),
            lambda m: (
                f"{m.group(1)} .pse files in this filter."
                if lang == "en"
                else f"{m.group(1)} archivos .pse en este filtrado."
            ),
        ),
        (
            re.compile(r"^Sessao salva em (.+)$"),
            lambda m: (
                f"Session saved to {m.group(1)}"
                if lang == "en"
                else f"Sesión guardada en {m.group(1)}"
            ),
        ),
        (
            re.compile(r"^Falhou \(exit code (\d+)\)$"),
            lambda m: (
                f"Failed (exit code {m.group(1)})"
                if lang == "en"
                else f"Falló (exit code {m.group(1)})"
            ),
        ),
    ]
    for pattern, repl in patterns:
        match = pattern.match(text)
        if match:
            return repl(match)
    return text


def install_translation_hooks() -> None:
    install_qt_hooks()
    install_matplotlib_hooks()


def install_qt_hooks() -> None:
    global _qt_hooks_installed
    if _qt_hooks_installed:
        return
    try:
        from PyQt6.QtCore import QObject
        from PyQt6.QtGui import QAction
        from PyQt6.QtWidgets import (
            QAbstractButton,
            QComboBox,
            QGroupBox,
            QLabel,
            QLineEdit,
            QMenu,
            QProgressBar,
            QTabWidget,
            QTableWidget,
            QWidget,
        )
    except Exception:
        return

    def wrap_simple(cls, method_name: str, source_attr: str) -> None:
        original = getattr(cls, method_name)
        if getattr(original, "_luna_i18n_wrapped", False):
            return

        def wrapped(self, value, *args, **kwargs):
            if isinstance(value, str):
                source = source_text(value)
                try:
                    setattr(self, source_attr, source)
                except Exception:
                    pass
                value = t(source)
            return original(self, value, *args, **kwargs)

        wrapped._luna_i18n_wrapped = True  # type: ignore[attr-defined]
        setattr(cls, method_name, wrapped)

    wrap_simple(QLabel, "setText", "_luna_i18n_text")
    wrap_simple(QAbstractButton, "setText", "_luna_i18n_text")
    wrap_simple(QGroupBox, "setTitle", "_luna_i18n_title")
    wrap_simple(QLineEdit, "setPlaceholderText", "_luna_i18n_placeholder")
    wrap_simple(QWidget, "setToolTip", "_luna_i18n_tooltip")
    wrap_simple(QAction, "setText", "_luna_i18n_text")
    wrap_simple(QAction, "setToolTip", "_luna_i18n_tooltip")
    wrap_simple(QMenu, "setTitle", "_luna_i18n_title")
    wrap_simple(QProgressBar, "setFormat", "_luna_i18n_format")

    combo_add_item = QComboBox.addItem

    def add_item(self, *args):
        args = list(args)
        text_index = 0 if args and isinstance(args[0], str) else 1 if len(args) > 1 and isinstance(args[1], str) else None
        user_data = args[text_index + 1] if text_index is not None and len(args) > text_index + 1 else None
        if text_index is not None and user_data is not None:
            source = source_text(args[text_index])
            args[text_index] = t(source)
        return combo_add_item(self, *args)

    QComboBox.addItem = add_item

    combo_set_item_text = QComboBox.setItemText

    def set_item_text(self, index: int, text: str):
        if self.itemData(index) is not None:
            text = t(source_text(text))
        return combo_set_item_text(self, index, text)

    QComboBox.setItemText = set_item_text

    tab_add_tab = QTabWidget.addTab

    def add_tab(self, *args):
        args = list(args)
        text_index = 1 if len(args) > 1 and isinstance(args[1], str) else 2 if len(args) > 2 and isinstance(args[2], str) else None
        if text_index is not None:
            args[text_index] = t(source_text(args[text_index]))
        return tab_add_tab(self, *args)

    QTabWidget.addTab = add_tab

    tab_set_text = QTabWidget.setTabText

    def set_tab_text(self, index: int, text: str):
        return tab_set_text(self, index, t(source_text(text)))

    QTabWidget.setTabText = set_tab_text

    table_set_headers = QTableWidget.setHorizontalHeaderLabels

    def set_horizontal_header_labels(self, labels):
        return table_set_headers(self, [t(source_text(label)) for label in labels])

    QTableWidget.setHorizontalHeaderLabels = set_horizontal_header_labels

    _qt_hooks_installed = True


def install_matplotlib_hooks() -> None:
    global _mpl_hooks_installed
    if _mpl_hooks_installed:
        return
    try:
        from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
        from matplotlib.text import Text
    except Exception:
        return

    original_set_text = Text.set_text
    if not getattr(original_set_text, "_luna_i18n_wrapped", False):

        def set_text(self, value):
            if isinstance(value, str):
                source = source_text(value)
                setattr(self, "_luna_i18n_source_text", source)
                value = t(source)
            return original_set_text(self, value)

        set_text._luna_i18n_wrapped = True  # type: ignore[attr-defined]
        Text.set_text = set_text

    original_draw = FigureCanvasQTAgg.draw
    if not getattr(original_draw, "_luna_i18n_wrapped", False):

        def draw(self, *args, **kwargs):
            fig = getattr(self, "figure", None)
            if fig is not None:
                translate_figure(fig)
            return original_draw(self, *args, **kwargs)

        draw._luna_i18n_wrapped = True  # type: ignore[attr-defined]
        FigureCanvasQTAgg.draw = draw

    _mpl_hooks_installed = True


def _translate_text_artist(artist) -> None:
    try:
        text = artist.get_text()
    except Exception:
        return
    if not text:
        return
    source = getattr(artist, "_luna_i18n_source_text", None) or source_text(text)
    setattr(artist, "_luna_i18n_source_text", source)
    new_text = t(source)
    if new_text != text:
        try:
            artist.set_text(new_text)
        except Exception:
            pass


def translate_figure(fig) -> None:
    if fig is None:
        return
    try:
        for text in list(getattr(fig, "texts", []) or []):
            _translate_text_artist(text)
        for ax in list(getattr(fig, "axes", []) or []):
            for text in [
                ax.title,
                ax.xaxis.label,
                ax.yaxis.label,
                *ax.get_xticklabels(),
                *ax.get_yticklabels(),
                *list(getattr(ax, "texts", []) or []),
            ]:
                _translate_text_artist(text)
            legend = ax.get_legend()
            if legend is not None:
                for text in legend.get_texts():
                    _translate_text_artist(text)
                title = legend.get_title()
                if title is not None:
                    _translate_text_artist(title)
    except Exception:
        return


def retranslate_ui(root) -> None:
    """Retranslate a Qt widget tree and any embedded matplotlib figures."""
    try:
        from PyQt6.QtGui import QAction
        from PyQt6.QtWidgets import (
            QAbstractButton,
            QComboBox,
            QGroupBox,
            QLabel,
            QLineEdit,
            QMenu,
            QProgressBar,
            QTabWidget,
            QTableWidget,
        )
    except Exception:
        return

    widgets = [root] + list(root.findChildren(QObject))
    for widget in widgets:
        if isinstance(widget, QLabel):
            source = getattr(widget, "_luna_i18n_text", None) or source_text(widget.text())
            widget.setText(source)
        if isinstance(widget, QAbstractButton):
            source = getattr(widget, "_luna_i18n_text", None) or source_text(widget.text())
            widget.setText(source)
        if isinstance(widget, QGroupBox):
            source = getattr(widget, "_luna_i18n_title", None) or source_text(widget.title())
            widget.setTitle(source)
        if isinstance(widget, QLineEdit):
            source = getattr(widget, "_luna_i18n_placeholder", None) or source_text(widget.placeholderText())
            widget.setPlaceholderText(source)
        if isinstance(widget, QProgressBar):
            source = getattr(widget, "_luna_i18n_format", None) or source_text(widget.format())
            widget.setFormat(source)
        if isinstance(widget, QComboBox):
            for index in range(widget.count()):
                if widget.itemData(index) is not None:
                    widget.setItemText(index, source_text(widget.itemText(index)))
        if isinstance(widget, QTabWidget):
            for index in range(widget.count()):
                widget.setTabText(index, source_text(widget.tabText(index)))
        if isinstance(widget, QTableWidget):
            for index in range(widget.columnCount()):
                item = widget.horizontalHeaderItem(index)
                if item is not None:
                    item.setText(t(source_text(item.text())))
        try:
            tooltip = widget.toolTip()
        except Exception:
            tooltip = ""
        if tooltip:
            widget.setToolTip(source_text(tooltip))
        fig = getattr(widget, "figure", None)
        if fig is not None:
            translate_figure(fig)
            try:
                widget.draw_idle()
            except Exception:
                pass

    if hasattr(root, "menuBar"):
        for action in root.menuBar().actions():
            _retranslate_action(action, QAction, QMenu)


def _retranslate_action(action, action_cls, menu_cls) -> None:
    if isinstance(action, action_cls):
        source = getattr(action, "_luna_i18n_text", None) or source_text(action.text())
        action.setText(source)
        if action.toolTip():
            action.setToolTip(source_text(action.toolTip()))
    menu = action.menu()
    if isinstance(menu, menu_cls):
        source = getattr(menu, "_luna_i18n_title", None) or source_text(menu.title())
        menu.setTitle(source)
        for child_action in menu.actions():
            _retranslate_action(child_action, action_cls, menu_cls)
