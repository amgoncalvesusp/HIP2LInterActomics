"""Small runtime translation layer for HIP²LInterActomics.

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
    "Sobre o HIP²LInterActomics": {"en": "About HIP²LInterActomics", "es": "Acerca de HIP²LInterActomics"},
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
    "1. Setup": {"pt": "1. Inicio", "en": "1. Home", "es": "1. Inicio"},
    "1. Inicio": {"pt": "1. Inicio", "en": "1. Home", "es": "1. Inicio"},
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
    "Limpar detecção": {"en": "Clear detection", "es": "Limpiar detección"},
    "Atualizar": {"en": "Refresh", "es": "Actualizar"},
    "Remover toda a lista": {"en": "Clear entire list", "es": "Eliminar toda la lista"},
    "Remove todos os projetos do histórico sem apagar seus arquivos.": {"en": "Removes all projects from history without deleting their files.", "es": "Elimina todos los proyectos del historial sin borrar sus archivos."},
    "Limpar histórico": {"en": "Clear history", "es": "Limpiar historial"},
    "Remover todos os projetos da lista? Os arquivos dos projetos serão preservados.": {"en": "Remove every project from the list? Project files will be preserved.", "es": "¿Eliminar todos los proyectos de la lista? Los archivos de los proyectos se conservarán."},
    "Cancelar": {"en": "Cancel", "es": "Cancelar"},
    "Fechar": {"en": "Close", "es": "Cerrar"},
    "OK": {"en": "OK", "es": "OK"},
    "Erro": {"en": "Error", "es": "Error"},
    "Sem dados": {"en": "No data", "es": "Sin datos"},
    "Processando...": {"en": "Processing...", "es": "Procesando..."},
    "Processando... (pode levar alguns minutos)": {"en": "Processing... (this may take a few minutes)", "es": "Procesando... (puede tardar unos minutos)"},

    # Setup tab.
    "Verificando ambiente...": {"en": "Checking environment...", "es": "Verificando entorno..."},
    "O que é o HIP²LInterActomics": {"en": "What is HIP²LInterActomics", "es": "Qué es HIP²LInterActomics"},
    "Configuração do ambiente e instalação de pacotes": {
        "en": "Environment setup and package installation",
        "es": "Configuración del entorno e instalación de paquetes",
    },
    (
        "HIP²LInterActomics é uma interface gráfica para estabelecer um fluxo de trabalho "
        "eficiente na análise de interações intermoleculares protein-ligand e protein-protein.\n\n"
        "O software ajuda a avaliar virtual screening, inspecionar poses de docking ou frames "
        "de dinâmica molecular, calcular fingerprints interpretáveis, organizar dados para "
        "modelos de machine learning e extrair pharmacophoric features de um alvo específico.\n\n"
        "Autores do software: Daniel Andrés Grajales Ruiz e Adriano Marques Gonçalves."
    ): {
        "en": (
            "HIP²LInterActomics is a graphical interface for an efficient workflow in "
            "protein-ligand and protein-protein intermolecular interaction analysis. The software "
            "helps evaluate virtual screening, inspect docking poses or molecular dynamics frames, "
            "calculate interpretable fingerprints, organize data for machine-learning models and "
            "extract pharmacophoric features from a specific target.\n\n"
            "Software authors: Daniel Andrés Grajales Ruiz and Adriano Marques Gonçalves."
        ),
        "es": (
            "HIP²LInterActomics es una interfaz gráfica para establecer un flujo de trabajo "
            "eficiente en el análisis de interacciones intermoleculares protein-ligand y protein-protein. "
            "El software ayuda a evaluar virtual screening, inspeccionar poses de docking o frames "
            "de dinámica molecular, calcular fingerprints interpretables, organizar datos para modelos "
            "de machine learning y extraer pharmacophoric features de un blanco específico.\n\n"
            "Autores del software: Daniel Andrés Grajales Ruiz y Adriano Marques Gonçalves."
        ),
    },
    (
        "O HIP²LInterActomics é uma interface gráfica para preparar, executar e interpretar análises "
        "baseadas no LUNA. Ele ajuda a avaliar resultados de virtual screening, comparar "
        "poses ou trajetórias de dinâmica molecular, construir modelos de machine learning "
        "com fingerprints de interação e extrair pharmacophoric features de um alvo "
        "específico. Nesta aba você prepara o ambiente: verifica o Conda, instala ou atualiza "
        "o ambiente separado 'luna-env' e garante que as dependências analíticas, como "
        "scikit-learn, estejam disponíveis."
    ): {
        "en": (
            "HIP²LInterActomics is a graphical interface for preparing, running and interpreting analyses "
            "based on LUNA. It helps evaluate virtual screening results, compare poses or "
            "molecular dynamics trajectories, build machine-learning models from interaction "
            "fingerprints and extract pharmacophoric features from a specific target. In this "
            "tab you prepare the environment: check Conda, install or update the isolated "
            "'luna-env' environment and ensure analytical dependencies, such as scikit-learn, "
            "are available."
        ),
        "es": (
            "HIP²LInterActomics es una interfaz gráfica para preparar, ejecutar e interpretar análisis "
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
        "pré-processar arquivos de complexos de docking ou frames de dinâmica e então informa proteína e ligantes para a análise."
    ): {
        "en": "In this tab you define the project mode, choose the working directory, optionally pre-process docking complexes or dynamics frames, then provide the protein and ligands for the analysis.",
        "es": "En esta pestaña defines el modo del proyecto, eliges el directorio de trabajo, puedes preprocesar complejos de docking o frames de dinámica y luego indicas la proteína y los ligandos para el análisis.",
    },
    "Entradas": {"en": "Inputs", "es": "Entradas"},
    "Modo do projeto:": {"en": "Project mode:", "es": "Modo del proyecto:"},
    "Fork de projeto existente (desmarcado = projeto novo)": {"en": "Fork an existing project (unchecked = new project)", "es": "Fork de un proyecto existente (desmarcado = proyecto nuevo)"},
    "Projeto fonte:": {"en": "Source project:", "es": "Proyecto fuente:"},
    "Diretório de trabalho:": {"en": "Working directory:", "es": "Directorio de trabajo:"},
    "Pré-processamento:": {"en": "Pre-processing:", "es": "Preprocesamiento:"},
    "Preparar arquivos de complexos...": {"en": "Prepare complex files...", "es": "Preparar archivos de complejos..."},
    "Proteína (PDB):": {"en": "Protein (PDB):", "es": "Proteína (PDB):"},
    "Ligantes (MOL2/SDF):": {"en": "Ligands (MOL2/SDF):", "es": "Ligandos (MOL2/SDF):"},
    "Incluir águas (HOH) — análise hidratada": {"en": "Include waters (HOH) - hydrated analysis", "es": "Incluir aguas (HOH) - análisis hidratado"},
    "Águas detectadas nos inputs: 0": {"en": "Waters detected in inputs: 0", "es": "Aguas detectadas en los inputs: 0"},
    "Análise de trajetória de dinâmica molecular/poses de docking (entradas = frames/poses)": {
        "en": "Molecular dynamics trajectory/docking pose analysis (entries = frames/poses)",
        "es": "Análisis de trayectoria de dinámica molecular/poses de docking (entradas = frames/poses)",
    },
    (
        "Marque quando cada entrada do projeto representar um frame de dinâmica molecular "
        "ou uma pose de docking. A aba Resultados > Estatísticas mostrará gráficos por "
        "frame/pose e percentuais de entradas por resíduo/interação."
    ): {
        "en": (
            "Check this when each project entry represents a molecular dynamics frame "
            "or a docking pose. Results > Statistics will show plots per frame/pose "
            "and entry percentages by residue/interaction."
        ),
        "es": (
            "Márcalo cuando cada entrada del proyecto representa un frame de dinámica molecular "
            "o una pose de docking. Resultados > Estadísticas mostrará gráficos por frame/pose "
            "y porcentajes de entradas por residuo/interacción."
        ),
    },
    "Ligantes detectados": {"en": "Detected ligands", "es": "Ligandos detectados"},
    (
        "Depois de carregar os ligantes, marque apenas os que devem entrar na análise. "
        "O filtro ajuda a localizar nomes específicos sem apagar a lista."
    ): {
        "en": "After loading ligands, check only those that should enter the analysis. The filter helps locate specific names without clearing the list.",
        "es": "Después de cargar los ligandos, marca solo los que deben entrar en el análisis. El filtro ayuda a localizar nombres específicos sin borrar la lista.",
    },
    "Filtrar por nome (texto livre)...": {"en": "Filter by name (free text)...", "es": "Filtrar por nombre (texto libre)..."},
    "Remove todos os ligantes detectados da lista atual sem apagar os arquivos de entrada.": {
        "en": "Removes all detected ligands from the current list without deleting input files.",
        "es": "Elimina todos los ligandos detectados de la lista actual sin borrar los archivos de entrada.",
    },
    "0 ligantes": {"en": "0 ligands", "es": "0 ligandos"},
    "Pasta com proteínas PDB": {"en": "Folder with PDB proteins", "es": "Carpeta con proteínas PDB"},
    "Selecionar pasta com proteínas PDB": {"en": "Select folder with PDB proteins", "es": "Seleccionar carpeta con proteínas PDB"},
    "Selecionar proteína": {"en": "Select protein", "es": "Seleccionar proteína"},
    "Selecionar ligantes": {"en": "Select ligands", "es": "Seleccionar ligandos"},
    "Pasta com arquivos MOL2/SDF": {"en": "Folder with MOL2/SDF files", "es": "Carpeta con archivos MOL2/SDF"},
    "Consolidado": {"en": "Consolidated", "es": "Consolidado"},
    "Próximo passo": {"en": "Next step", "es": "Próximo paso"},
    "Pasta com arquivos de complexos": {"en": "Folder with complex files", "es": "Carpeta con archivos de complejos"},
    (
        "Selecione uma pasta com arquivos .mol2, .pdb ou .ent (proteína + ligante no mesmo arquivo).\n"
        "Serão geradas subpastas separadas para proteínas e ligantes compatíveis com LUNA: MOL2 quando a origem é MOL2, SDF quando a origem é PDB/ENT, com águas preservadas junto à proteína."
    ): {
        "en": (
            "Select a folder with .mol2, .pdb or .ent files (protein + ligand in the same file).\n"
            "Separate folders will be generated for proteins and LUNA-compatible ligands: MOL2 when the source is MOL2, SDF when the source is PDB/ENT, with waters preserved with the protein."
        ),
        "es": (
            "Selecciona una carpeta con archivos .mol2, .pdb o .ent (proteína + ligando en el mismo archivo).\n"
            "Se generarán carpetas separadas para proteínas y ligandos compatibles con LUNA: MOL2 cuando el origen es MOL2, SDF cuando el origen es PDB/ENT, con las aguas preservadas junto a la proteína."
        ),
    },
    "Divide complexos MOL2/PDB (proteína+ligante+águas) em entradas separadas.": {
        "en": "Splits MOL2/PDB complexes (protein+ligand+waters) into separate inputs.",
        "es": "Divide complejos MOL2/PDB (proteína+ligando+aguas) en entradas separadas.",
    },
    (
        "Arquivo PDB da proteína para todos os ligantes, ou pasta com um PDB por ligante/complexo.\n"
        "Quando uma pasta é usada, cada ligante será pareado com o PDB de mesmo nome-base."
    ): {
        "en": "Protein PDB file for all ligands, or a folder with one PDB per ligand/complex.\nWhen a folder is used, each ligand is paired with the PDB that has the same base name.",
        "es": "Archivo PDB de proteína para todos los ligandos, o carpeta con un PDB por ligando/complejo.\nCuando se usa una carpeta, cada ligando se empareja con el PDB del mismo nombre base.",
    },
    "Seleciona um arquivo PDB de proteína para todos os ligantes.": {
        "en": "Selects one protein PDB file for all ligands.",
        "es": "Selecciona un archivo PDB de proteína para todos los ligandos.",
    },
    "Seleciona uma pasta com um arquivo .pdb para cada ligante/complexo.": {
        "en": "Selects a folder with one .pdb file for each ligand/complex.",
        "es": "Selecciona una carpeta con un archivo .pdb para cada ligando/complejo.",
    },
    "Arquivo ou pasta com ligantes. Pode ser MOL2/SDF/PDB único, pasta de PDBs ou arquivo consolidado com vários ligantes.": {
        "en": "Ligand file or folder. It can be a single MOL2/SDF/PDB, a PDB folder, or a consolidated file with multiple ligands.",
        "es": "Archivo o carpeta con ligandos. Puede ser un MOL2/SDF/PDB único, una carpeta de PDBs o un archivo consolidado con varios ligandos.",
    },

    # Analyses tab.
    (
        "Aqui você escolhe quais saídas o LUNA vai gerar. Se estiver em dúvida, deixe os "
        "valores padrão e ative apenas o que pretende inspecionar depois: fingerprints, "
        "similaridade, sessões PyMOL ou filtros específicos."
    ): {
        "en": "Here you choose which outputs LUNA will generate. If unsure, keep the defaults and enable only what you plan to inspect later: fingerprints, similarity, PyMOL sessions or specific filters.",
        "es": "Aquí eliges qué salidas generará LUNA. Si tienes dudas, deja los valores por defecto y activa solo lo que quieras inspeccionar después: fingerprints, similitud, sesiones PyMOL o filtros específicos.",
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
    "Rótulos para importância de fingerprints": {"en": "Labels for fingerprint importance", "es": "Etiquetas para importancia de fingerprints"},
    (
        "Se marcado, a aba FP análises usa este CSV para treinar a importância das features. "
        "Informe o arquivo, a coluna do ID do ligante e a coluna que contém os rótulos/classes."
    ): {
        "en": "When checked, the FP analyses tab uses this CSV to train feature importance. Provide the file, the ligand ID column and the column containing labels/classes.",
        "es": "Cuando está marcado, la pestaña de análisis FP usa este CSV para entrenar la importancia de las features. Indica el archivo, la columna del ID del ligando y la columna con etiquetas/clases.",
    },
    "Arquivo CSV:": {"en": "CSV file:", "es": "Archivo CSV:"},
    "Coluna do ligand_id:": {"en": "ligand_id column:", "es": "Columna ligand_id:"},
    "Coluna de rótulo:": {"en": "Label column:", "es": "Columna de rótulo:"},
    "Tarefa:": {"en": "Task:", "es": "Tarea:"},
    "Regressão": {"en": "Regression", "es": "Regresión"},
    "Classificação": {"en": "Classification", "es": "Clasificación"},
    "Aplicar Otsu tambem a interacoes/residuos": {"pt": "Aplicar Otsu também a interações/resíduos", "en": "Apply Otsu also to interactions/residues", "es": "Aplicar Otsu también a interacciones/residuos"},
    "Matriz de similaridade (Tanimoto)": {"en": "Similarity matrix (Tanimoto)", "es": "Matriz de similitud (Tanimoto)"},
    (
        "A matriz de similaridade ajuda a ver quais ligantes se comportam de forma parecida. "
        "Ela também alimenta os gráficos e clusters da aba de resultados."
    ): {
        "en": "The similarity matrix helps identify ligands with similar interaction behavior. It also feeds charts and clusters in the results tab.",
        "es": "La matriz de similitud ayuda a ver qué ligandos se comportan de forma parecida. También alimenta los gráficos y clusters de la pestaña de resultados.",
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
        "es": "Usa esta pestaña para revisar lo que LUNA generó. Puedes cargar un workdir existente, visualizar tablas y gráficos, exportar figuras y agrupar ligandos por similitud.",
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
    "Mapa de calor da similaridade entre ligantes. Valores mais altos indicam fingerprints mais parecidos.": {
        "en": "Heatmap of ligand similarity. Higher values indicate more similar fingerprints.",
        "es": "Mapa de calor de la similitud entre ligandos. Valores más altos indican fingerprints más parecidos.",
    },
    "Estatísticas": {"en": "Statistics", "es": "Estadísticas"},
    "Resume quantas interações de cada tipo aparecem no conjunto analisado. É uma visão global do perfil químico observado.": {
        "en": "Summarizes how many interactions of each type appear in the analyzed set. It is a global view of the observed chemical profile.",
        "es": "Resume cuántas interacciones de cada tipo aparecen en el conjunto analizado. Es una visión global del perfil químico observado.",
    },
    "Calcular estatísticas (usa luna-env)": {"en": "Calculate statistics (uses luna-env)", "es": "Calcular estadísticas (usa luna-env)"},
    "Heatmap por tipo": {"pt": "Mapa de calor por tipo", "en": "Heatmap by type", "es": "Mapa de calor por tipo"},
    "Mapa de calor por tipo": {"en": "Heatmap by type", "es": "Mapa de calor por tipo"},
    "Mostra, para um tipo de interação escolhido, quais resíduos da proteína aparecem associados aos ligantes.": {
        "en": "Shows, for the selected interaction type, which protein residues are associated with ligands.",
        "es": "Muestra, para un tipo de interacción elegido, qué residuos de la proteína aparecen asociados a los ligandos.",
    },
    "Calcular heatmap (usa luna-env)": {"pt": "Calcular mapa de calor (usa luna-env)", "en": "Calculate heatmap (uses luna-env)", "es": "Calcular mapa de calor (usa luna-env)"},
    "Método:": {"en": "Method:", "es": "Método:"},
    "Clusters:": {"en": "Clusters:", "es": "Clusters:"},
    "Atualizar clusters": {"en": "Update clusters", "es": "Actualizar clusters"},
    "Agrupa ligantes com base na matriz de similaridade. O dendrograma mostra a relação hierárquica, e a tabela lista o cluster atribuído a cada ligante.": {
        "en": "Clusters ligands based on the similarity matrix. The dendrogram shows the hierarchical relationship, and the table lists the cluster assigned to each ligand.",
        "es": "Agrupa ligandos con base en la matriz de similitud. El dendrograma muestra la relación jerárquica y la tabla lista el cluster atribuido a cada ligando.",
    },
    "Exportar clusters CSV": {"en": "Export clusters CSV", "es": "Exportar clusters CSV"},
    "Abrir no PyMOL": {"en": "Open in PyMOL", "es": "Abrir en PyMOL"},
    "Lista as sessões PyMOL exportadas pelo LUNA. Abra um arquivo para inspecionar visualmente interações e poses.": {
        "en": "Lists the PyMOL sessions exported by LUNA. Open a file to visually inspect interactions and poses.",
        "es": "Lista las sesiones PyMOL exportadas por LUNA. Abre un archivo para inspeccionar visualmente interacciones y poses.",
    },
    "Estatisticas": {"pt": "Estatísticas", "en": "Statistics", "es": "Estadísticas"},
    "Heatmap completo ligantes x residuos": {"pt": "Mapa de calor completo ligantes x resíduos", "en": "Complete ligand x residue heatmap", "es": "Mapa de calor completo ligandos x residuos"},
    "Sessoes em PyMOL": {"pt": "Triagem virtual racional", "en": "Rational virtual screening", "es": "Cribado virtual racional"},
    "Triagem virtual racional": {"en": "Rational virtual screening", "es": "Cribado virtual racional"},
    "FP analises": {"pt": "Análises FP", "en": "FP analyses", "es": "Análisis FP"},
    "Análises FP": {"en": "FP analyses", "es": "Análisis FP"},
    "FP sessao em PyMOL": {"pt": "Sessão FP em PyMOL", "en": "FP PyMOL session", "es": "Sesión FP en PyMOL"},
    "Sessão FP em PyMOL": {"en": "FP PyMOL session", "es": "Sesión FP en PyMOL"},
    "Descrição do método classificação e seleção": {
        "en": "Classification and selection method description",
        "es": "Descripción del método de clasificación y selección",
    },
    "Descrição do método: classificação e seleção": {
        "en": "Classification and selection method description",
        "es": "Descripción del método: clasificación y selección",
    },
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
    "Full screen": {"pt": "Tela cheia", "en": "Full screen", "es": "Pantalla completa"},
    "Afastar": {"en": "Zoom out", "es": "Alejar"},
    "Aproximar": {"en": "Zoom in", "es": "Acercar"},
    "Ajustar a tela": {"pt": "Ajustar à tela", "en": "Fit to window", "es": "Ajustar a pantalla"},
    "Ligar/desligar todas": {"en": "Toggle all", "es": "Activar/desactivar todas"},
    "Ocultar desfav./repulsivas": {
        "en": "Hide unfavorable/repulsive",
        "es": "Ocultar desfav./repulsivas",
    },
    "Ocultar empilhamentos": {"en": "Hide stacking", "es": "Ocultar apilamientos"},
    "Fingerprints relevantes:": {"en": "Relevant fingerprints:", "es": "Fingerprints relevantes:"},
    "Aplicar selecao": {"pt": "Aplicar seleção", "en": "Apply selection", "es": "Aplicar selección"},
    "Usar linhas selecionadas": {"en": "Use selected rows", "es": "Usar filas seleccionadas"},
    "Restabelecer calculado": {"en": "Restore calculated", "es": "Restablecer calculado"},
    "Atualizar heatmaps": {"pt": "Atualizar mapas de calor", "en": "Update heatmaps", "es": "Actualizar mapas de calor"},
    "Atualizar mapas de calor": {"en": "Update heatmaps", "es": "Actualizar mapas de calor"},
    (
        "Mostra, em cada par ligante x residuo, todas as interacoes presentes em faixas coloridas. "
        "As cores indicam a classe da interacao, nao a contagem."
    ): {
        "en": "Shows all interactions present in each ligand x residue pair as colored bands. Colors indicate interaction class, not count.",
        "es": "Muestra, en cada par ligando x residuo, todas las interacciones presentes como franjas de color. Los colores indican la clase de interacción, no el conteo.",
    },
    "Carregar analises de FP": {"pt": "Carregar análises de FP", "en": "Load FP analyses", "es": "Cargar análisis FP"},
    "Carregar análises de FP": {"en": "Load FP analyses", "es": "Cargar análisis FP"},
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
    "Frequencia por classe": {"pt": "Frequência por classe", "en": "Frequency by class", "es": "Frecuencia por clase"},
    "Frequência por classe": {"en": "Frequency by class", "es": "Frecuencia por clase"},
    "Cobertura e importancia": {"pt": "Cobertura e importância", "en": "Coverage and importance", "es": "Cobertura e importancia"},
    "Cobertura e importância": {"en": "Coverage and importance", "es": "Cobertura e importancia"},
    "Heatmap importância": {"pt": "Mapa de calor de importância", "en": "Importance heatmap", "es": "Mapa de calor de importancia"},
    "Mapa de calor de importância": {"en": "Importance heatmap", "es": "Mapa de calor de importancia"},
    "Frequencia interacoes": {"pt": "Frequência de interações", "en": "Interaction frequency", "es": "Frecuencia de interacciones"},
    "Frequência de interações": {"en": "Interaction frequency", "es": "Frecuencia de interacciones"},
    "Interações prevalentes": {"en": "Prevalent interactions", "es": "Interacciones prevalentes"},
    "Heatmap interações": {"pt": "Mapa de calor de interações", "en": "Interaction heatmap", "es": "Mapa de calor de interacciones"},
    "Mapa de calor de interações": {"en": "Interaction heatmap", "es": "Mapa de calor de interacciones"},
    "Fingerprint:": {"en": "Fingerprint:", "es": "Fingerprint:"},
    "Molecula:": {"pt": "Molécula:", "en": "Molecule:", "es": "Molécula:"},
    "Molécula:": {"en": "Molecule:", "es": "Molécula:"},
    "Gerar sessao": {"pt": "Gerar sessão", "en": "Generate session", "es": "Generar sesión"},
    "Gerar sessão": {"en": "Generate session", "es": "Generar sesión"},
    (
        "Gera uma sessao PyMOL para o fingerprint escolhido, recuperando os shells "
        "que originaram aquele bit para um ligante especifico."
    ): {
        "en": "Generates a PyMOL session for the selected fingerprint, recovering the shells that originated that bit for a specific ligand.",
        "es": "Genera una sesión PyMOL para el fingerprint elegido, recuperando los shells que originaron ese bit para un ligando específico.",
    },
    "Abrir sessao": {"pt": "Abrir sessão", "en": "Open session", "es": "Abrir sesión"},
    "Abrir sessão": {"en": "Open session", "es": "Abrir sesión"},

    # Tables.
    "Feature": {"en": "Feature", "es": "Feature"},
    "Moleculas": {"pt": "Moléculas", "en": "Molecules", "es": "Moléculas"},
    "Moléculas": {"en": "Molecules", "es": "Moléculas"},
    "Cobertura (%)": {"en": "Coverage (%)", "es": "Cobertura (%)"},
    "Classe prevalente (%)": {"en": "Prevalent class (%)", "es": "Clase prevalente (%)"},
    "Z-score classe": {"en": "Class z-score", "es": "Z-score clase"},
    "Classe atribuida": {"pt": "Classe atribuída", "en": "Assigned class", "es": "Clase asignada"},
    "Classe atribuída": {"en": "Assigned class", "es": "Clase asignada"},
    "Importancia": {"pt": "Importância", "en": "Importance", "es": "Importancia"},
    "Importância": {"en": "Importance", "es": "Importancia"},
    "Z-score Importance": {"en": "Importance z-score", "es": "Z-score importancia"},
    "Colisoes": {"pt": "Colisões", "en": "Collisions", "es": "Colisiones"},
    "Colisões": {"en": "Collisions", "es": "Colisiones"},
    "Perfil da base": {"en": "Dataset profile", "es": "Perfil de la base"},
    "Nível assinado": {"en": "Assigned level", "es": "Nivel asignado"},
    "Níveis shell": {"en": "Shell levels", "es": "Niveles shell"},
    "Níveis colisão": {"en": "Collision levels", "es": "Niveles de colisión"},
    "Ligante": {"en": "Ligand", "es": "Ligando"},
    "Ordem": {"en": "Order", "es": "Orden"},

    # Binding mode editor and preparation dialog.
    "Editor de Binding Modes": {"en": "Binding Modes Editor", "es": "Editor de Binding Modes"},
    "Tipo de interação": {"en": "Interaction type", "es": "Tipo de interacción"},
    "accept_only (lista)": {"en": "accept_only (list)", "es": "accept_only (lista)"},
    "+ Adicionar tipo": {"en": "+ Add type", "es": "+ Agregar tipo"},
    "− Remover linha": {"en": "- Remove row", "es": "- Remover fila"},
    "- Remover linha": {"en": "- Remove row", "es": "- Remover fila"},
    "Carregar .cfg": {"en": "Load .cfg", "es": "Cargar .cfg"},
    "Salvar .cfg": {"en": "Save .cfg", "es": "Guardar .cfg"},
    "Preparar arquivos de complexos": {"en": "Prepare complex files", "es": "Preparar archivos de complejos"},
    "Pasta de origem:": {"en": "Source folder:", "es": "Carpeta de origen:"},
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
    "Similaridade de Tanimoto": {"en": "Tanimoto similarity", "es": "Similitud de Tanimoto"},
    "Sem dados de interação": {"en": "No interaction data", "es": "Sin datos de interacción"},
    "Sem interações para esta visão.": {"en": "No interactions for this view.", "es": "Sin interacciones para esta vista."},
    "Aminoácidos": {"en": "Amino acids", "es": "Aminoácidos"},
    "Átomos do ligante": {"en": "Ligand atoms", "es": "Átomos del ligando"},
    "Estrutura 2D\nIDs dos átomos": {"en": "2D structure\nAtom IDs", "es": "Estructura 2D\nIDs de átomos"},
    "Interações por átomo do ligante": {"en": "Interactions by ligand atom", "es": "Interacciones por átomo del ligando"},
    "Residuos": {"en": "Residues", "es": "Residuos"},
    "Resíduos": {"en": "Residues", "es": "Residuos"},
    "Ligantes": {"en": "Ligands", "es": "Ligandos"},
    "Frames": {"en": "Frames", "es": "Frames"},
    "Ausente": {"en": "Absent", "es": "Ausente"},
    "% de ligantes": {"en": "% of ligands", "es": "% de ligandos"},
    "% de frames": {"en": "% of frames", "es": "% de frames"},
    "Total (todas as entradas)": {"en": "Total (all entries)", "es": "Total (todas las entradas)"},
    "Total neste ligante": {"en": "Total for this ligand", "es": "Total en este ligando"},
    "Contagem por tipo de interação": {"en": "Count by interaction type", "es": "Conteo por tipo de interacción"},
    "Interações por tipo": {"en": "Interactions by type", "es": "Interacciones por tipo"},
    "Tipos de interacao por par ligante x residuo": {"pt": "Tipos de interação por par ligante x resíduo", "en": "Interaction types by ligand x residue pair", "es": "Tipos de interacción por par ligando x residuo"},
    "Tipos de interação por par ligante x resíduo": {"en": "Interaction types by ligand x residue pair", "es": "Tipos de interacción por par ligando x residuo"},
    "Clique para ocultar/mostrar": {"en": "Click to hide/show", "es": "Click para ocultar/mostrar"},
    "Sem matriz de residuos disponivel": {"pt": "Sem matriz de resíduos disponível", "en": "No residue matrix available", "es": "Sin matriz de residuos disponible"},
    "Sem matriz de resíduos disponível": {"en": "No residue matrix available", "es": "Sin matriz de residuos disponible"},
    "Sem matriz completa disponivel": {"pt": "Sem matriz completa disponível", "en": "No complete matrix available", "es": "Sin matriz completa disponible"},
    "Sem matriz completa disponível": {"en": "No complete matrix available", "es": "Sin matriz completa disponible"},
    "Sem dados de interacao para o heatmap completo": {"pt": "Sem dados de interação para o mapa de calor completo", "en": "No interaction data for the complete heatmap", "es": "Sin datos de interacción para el mapa de calor completo"},
    "Sem features importantes para resumir.": {"en": "No important features to summarize.", "es": "Sin features importantes para resumir."},
    "Sem features importantes para plotar.": {"en": "No important features to plot.", "es": "Sin features importantes para graficar."},
    "Sem features importantes para gerar o heatmap.": {"pt": "Sem features importantes para gerar o mapa de calor.", "en": "No important features to generate the heatmap.", "es": "Sin features importantes para generar el mapa de calor."},
    "Sem features importantes para gerar o mapa de calor.": {"en": "No important features to generate the heatmap.", "es": "Sin features importantes para generar el mapa de calor."},
    "Sem interacoes prevalentes confiaveis nas features importantes.": {"pt": "Sem interações prevalentes confiáveis nas features importantes.", "en": "No reliable prevalent interactions among important features.", "es": "Sin interacciones prevalentes confiables en las features importantes."},
    "Sem interações prevalentes confiáveis nas features importantes.": {"en": "No reliable prevalent interactions among important features.", "es": "Sin interacciones prevalentes confiables en las features importantes."},
    "Todas as interacoes estao ocultas pela legenda.": {"pt": "Todas as interações estão ocultas pela legenda.", "en": "All interactions are hidden by the legend.", "es": "Todas las interacciones están ocultas por la leyenda."},
    "Todas as interações estão ocultas pela legenda.": {"en": "All interactions are hidden by the legend.", "es": "Todas las interacciones están ocultas por la leyenda."},
    "Sem interacoes prevalentes para gerar o heatmap.": {"pt": "Sem interações prevalentes para gerar o mapa de calor.", "en": "No prevalent interactions to generate the heatmap.", "es": "Sin interacciones prevalentes para generar el mapa de calor."},
    "% Important features": {"pt": "% de features importantes", "en": "% Important features", "es": "% de features importantes"},
    "% de features importantes": {"en": "% Important features", "es": "% de features importantes"},
    "Distribuicao das classes entre as features mais importantes": {"pt": "Distribuição das classes entre as features mais importantes", "en": "Class distribution among the most important features", "es": "Distribución de clases entre las features más importantes"},
    "Distribuição das classes entre as features mais importantes": {"en": "Class distribution among the most important features", "es": "Distribución de clases entre las features más importantes"},
    "Assignment frequency of each class (%)": {"pt": "Frequência de atribuição de cada classe (%)", "en": "Assignment frequency of each class (%)", "es": "Frecuencia de asignación de cada clase (%)"},
    "Frequência de atribuição de cada classe (%)": {"en": "Assignment frequency of each class (%)", "es": "Frecuencia de asignación de cada clase (%)"},
    "Feature id": {"pt": "ID da feature", "en": "Feature id", "es": "ID de feature"},
    "ID da feature": {"en": "Feature id", "es": "ID de feature"},
    "Frequencia de atribuicao de classes nas features importantes": {"pt": "Frequência de atribuição de classes nas features importantes", "en": "Class assignment frequency in important features", "es": "Frecuencia de asignación de clases en features importantes"},
    "Frequência de atribuição de classes nas features importantes": {"en": "Class assignment frequency in important features", "es": "Frecuencia de asignación de clases en features importantes"},
    "% Fingerprints containing the feature": {"pt": "% de fingerprints contendo a feature", "en": "% Fingerprints containing the feature", "es": "% de fingerprints que contienen la feature"},
    "% de fingerprints contendo a feature": {"en": "% Fingerprints containing the feature", "es": "% de fingerprints que contienen la feature"},
    "Cobertura das features importantes e importancia do modelo": {"pt": "Cobertura das features importantes e importância do modelo", "en": "Important feature coverage and model importance", "es": "Cobertura de features importantes e importancia del modelo"},
    "Cobertura das features importantes e importância do modelo": {"en": "Important feature coverage and model importance", "es": "Cobertura de features importantes e importancia del modelo"},
    "Importância relativa": {"en": "Relative importance", "es": "Importancia relativa"},
    "Cobertura / importância relativa (%)": {"en": "Coverage / relative importance (%)", "es": "Cobertura / importancia relativa (%)"},
    "Mapa de presenca das features importantes por classe": {"pt": "Mapa de presença das features importantes por classe", "en": "Presence map of important features by class", "es": "Mapa de presencia de features importantes por clase"},
    "Mapa de presença das features importantes por classe": {"en": "Presence map of important features by class", "es": "Mapa de presencia de features importantes por clase"},
    "Assignment frequency of each interaction (%)": {"pt": "Frequência de atribuição de cada interação (%)", "en": "Assignment frequency of each interaction (%)", "es": "Frecuencia de asignación de cada interacción (%)"},
    "Frequência de atribuição de cada interação (%)": {"en": "Assignment frequency of each interaction (%)", "es": "Frecuencia de asignación de cada interacción (%)"},
    "Frequencia de atribuicao da interacao prevalente nas features importantes": {"pt": "Frequência de atribuição da interação prevalente nas features importantes", "en": "Prevalent interaction assignment frequency in important features", "es": "Frecuencia de asignación de la interacción prevalente en features importantes"},
    "Frequência de atribuição da interação prevalente nas features importantes": {"en": "Prevalent interaction assignment frequency in important features", "es": "Frecuencia de asignación de la interacción prevalente en features importantes"},
    "Numero de ligantes": {"pt": "Número de ligantes", "en": "Number of ligands", "es": "Número de ligandos"},
    "Número de ligantes": {"en": "Number of ligands", "es": "Número de ligandos"},
    "Interacao e residuo prevalentes nas features importantes": {"pt": "Interação e resíduo prevalentes nas features importantes", "en": "Prevalent interaction and residue in important features", "es": "Interacción y residuo prevalentes en features importantes"},
    "Interação e resíduo prevalentes nas features importantes": {"en": "Prevalent interaction and residue in important features", "es": "Interacción y residuo prevalentes en features importantes"},
    "Interacoes prevalentes das features importantes por ligante": {"pt": "Interações prevalentes das features importantes por ligante", "en": "Prevalent interactions of important features by ligand", "es": "Interacciones prevalentes de features importantes por ligando"},
    "Interações prevalentes das features importantes por ligante": {"en": "Prevalent interactions of important features by ligand", "es": "Interacciones prevalentes de features importantes por ligando"},
    "Mapa de calor de interações prevalentes": {"en": "Prevalent interaction heatmap", "es": "Mapa de calor de interacciones prevalentes"},
    "Importância da feature": {"en": "Feature importance", "es": "Importancia de la feature"},
    "ID da feature e nível": {"en": "Feature ID and level", "es": "ID de feature y nivel"},
    "Top 50 features": {"en": "Top 50 features", "es": "Top 50 features"},
    "Top 50 features por Extra Trees": {"en": "Top 50 features by Extra Trees", "es": "Top 50 features por Extra Trees"},
    "Top 50 features por Gradient Boosting": {"en": "Top 50 features by Gradient Boosting", "es": "Top 50 features por Gradient Boosting"},
    "Posição": {"en": "Rank", "es": "Posición"},
    "Coluna": {"en": "Column", "es": "Columna"},
    "Interpretação": {"en": "Interpretation", "es": "Interpretación"},
    "Como interpretar as análises de fingerprints": {"en": "How to interpret fingerprint analyses", "es": "Cómo interpretar los análisis de fingerprints"},
    "Guia das colunas de Análises FP": {"en": "FP Analyses column guide", "es": "Guía de columnas de Análisis FP"},
    "Consulte este dicionário ao interpretar as tabelas, os rankings e os mapas de fingerprints.": {
        "en": "Use this dictionary when interpreting fingerprint tables, rankings, and maps.",
        "es": "Consulte este diccionario al interpretar las tablas, los rankings y los mapas de fingerprints.",
    },
    "Resumo das análises de fingerprints": {"en": "Fingerprint analysis summary", "es": "Resumen de los análisis de fingerprints"},
    "A seção Análises FP transforma os bits dos fingerprints em variáveis interpretáveis. Ela combina cobertura, natureza química, nível de shell, colisões e importância preditiva para priorizar padrões que merecem inspeção estrutural.": {
        "en": "The FP Analyses section turns fingerprint bits into interpretable variables. It combines coverage, chemical nature, shell level, collisions and predictive importance to prioritize patterns that warrant structural inspection.",
        "es": "La sección Análisis FP transforma los bits de fingerprints en variables interpretables. Combina cobertura, naturaleza química, nivel de shell, colisiones e importancia predictiva para priorizar patrones que requieren inspección estructural.",
    },
    "Extra Trees e Gradient Boosting são ajustados para a tarefa ativa, como classificadores para rótulos discretos ou regressores para valores contínuos. Os rankings devem ser comparados: concordância entre os métodos reforça a estabilidade, enquanto divergências sinalizam dependência do modelo.": {
        "en": "Extra Trees and Gradient Boosting are fitted to the active task, as classifiers for discrete labels or regressors for continuous values. Compare their rankings: agreement supports stability, while divergence signals model dependence.",
        "es": "Extra Trees y Gradient Boosting se ajustan a la tarea activa, como clasificadores para etiquetas discretas o regresores para valores continuos. Compare sus rankings: la concordancia respalda la estabilidad y las divergencias señalan dependencia del modelo.",
    },
    "Importância não prova causalidade molecular. A priorização final deve considerar cobertura, colisões, p-value, interações prevalentes e a posição da feature na estrutura do ligante e do receptor.": {
        "en": "Importance does not prove molecular causality. Final prioritization should consider coverage, collisions, p-value, prevalent interactions and the feature position in ligand and receptor structures.",
        "es": "La importancia no demuestra causalidad molecular. La priorización final debe considerar cobertura, colisiones, p-value, interacciones prevalentes y la posición de la feature en las estructuras del ligando y receptor.",
    },
    "Identificador do bit/atributo do fingerprint usado para localizar a mesma feature em tabelas, gráficos e sessões estruturais.": {"en": "Fingerprint bit/attribute identifier used to locate the same feature in tables, charts and structural sessions.", "es": "Identificador del bit/atributo del fingerprint usado para localizar la misma feature en tablas, gráficos y sesiones estructurales."},
    "Número de ligantes em que a feature está presente.": {"en": "Number of ligands in which the feature is present.", "es": "Número de ligandos en los que está presente la feature."},
    "Percentual da base que contém a feature, calculado sobre todos os ligantes processados.": {"en": "Percentage of the dataset containing the feature, calculated over all processed ligands.", "es": "Porcentaje de la base que contiene la feature, calculado sobre todos los ligandos procesados."},
    "Maior participação percentual entre as classes químicas observadas para a feature.": {"en": "Largest percentage share among the chemical classes observed for the feature.", "es": "Mayor participación porcentual entre las clases químicas observadas para la feature."},
    "Distância padronizada entre a prevalência da classe e a distribuição das demais features.": {"en": "Standardized distance between class prevalence and the distribution of the remaining features.", "es": "Distancia estandarizada entre la prevalencia de la clase y la distribución de las demás features."},
    "Natureza química aceita para a feature após a aplicação do critério de confiabilidade.": {"en": "Chemical nature accepted for the feature after applying the reliability criterion.", "es": "Naturaleza química aceptada para la feature después de aplicar el criterio de confiabilidad."},
    "Peso fornecido pelo modelo supervisionado ou pelo fallback analítico para a tarefa configurada.": {"en": "Weight supplied by the supervised model or analytical fallback for the configured task.", "es": "Peso proporcionado por el modelo supervisado o el fallback analítico para la tarea configurada."},
    "Importância padronizada dentro do nível de fingerprint correspondente.": {"en": "Importance standardized within the corresponding fingerprint level.", "es": "Importancia estandarizada dentro del nivel de fingerprint correspondiente."},
    "Probabilidade de cauda derivada do Z-score de importância; valores menores indicam maior evidência de relevância.": {"en": "Tail probability derived from the importance Z-score; lower values indicate stronger evidence of relevance.", "es": "Probabilidad de cola derivada del Z-score de importancia; valores menores indican mayor evidencia de relevancia."},
    "Quantidade de ocorrências em que o mesmo bit agregou shells ou naturezas químicas distintas.": {"en": "Number of occurrences in which the same bit aggregated distinct shells or chemical natures.", "es": "Cantidad de ocurrencias en las que el mismo bit agregó shells o naturalezas químicas distintas."},
    "Nível do fingerprint atribuído à feature para separar modelos e interpretações por escala estrutural.": {"en": "Fingerprint level assigned to the feature to separate models and interpretations by structural scale.", "es": "Nivel del fingerprint asignado a la feature para separar modelos e interpretaciones por escala estructural."},
    "Distribuição dos níveis de shell efetivamente associados à feature.": {"en": "Distribution of shell levels effectively associated with the feature.", "es": "Distribución de los niveles de shell efectivamente asociados a la feature."},
    "Níveis de shell encontrados nas ocorrências classificadas como colisão.": {"en": "Shell levels found in occurrences classified as collisions.", "es": "Niveles de shell encontrados en las ocurrencias clasificadas como colisión."},
    "Resumo de contagens e percentuais de classes para todas as ocorrências da feature na base.": {"en": "Summary of class counts and percentages for every occurrence of the feature in the dataset.", "es": "Resumen de conteos y porcentajes de clases para todas las ocurrencias de la feature en la base."},
    "Clustering hierárquico": {"en": "Hierarchical clustering", "es": "Clustering jerárquico"},
    "Distância": {"en": "Distance", "es": "Distancia"},
    "Matriz reordenada por cluster": {"en": "Matrix reordered by cluster", "es": "Matriz reordenada por cluster"},
    "Parâmetro": {"en": "Parameter", "es": "Parámetro"},
    "Valor": {"en": "Value", "es": "Valor"},
    "Interações por aminoácido no conjunto de ligantes": {"en": "Interactions by amino acid across all ligands", "es": "Interacciones por aminoácido en el conjunto de ligandos"},
    "Interações por aminoácido ao longo da trajetória": {"en": "Interactions by amino acid across the trajectory", "es": "Interacciones por aminoácido a lo largo de la trayectoria"},
    "% de frames (entradas)": {"en": "% of frames (entries)", "es": "% de frames (entradas)"},
    "Número de interações no ligante": {"en": "Number of interactions in the ligand", "es": "Número de interacciones en el ligando"},
    "Número de interações no frame": {"en": "Number of interactions in the frame", "es": "Número de interacciones en el frame"},
    "Sem matriz resíduo x interação.\nCalcule o heatmap por tipo ou reexecute a análise.": {"pt": "Sem matriz resíduo x interação.\nCalcule o mapa de calor por tipo ou reexecute a análise.", "en": "No residue x interaction matrix.\nCalculate the heatmap by type or rerun the analysis.", "es": "Sin matriz residuo x interacción.\nCalcula el mapa de calor por tipo o reejecuta el análisis."},
    "Sem matriz resíduo x interação para trajetória.\nCalcule o heatmap por tipo ou reexecute a análise.": {"pt": "Sem matriz resíduo x interação para trajetória.\nCalcule o mapa de calor por tipo ou reexecute a análise.", "en": "No residue x interaction matrix for trajectory.\nCalculate the heatmap by type or rerun the analysis.", "es": "Sin matriz residuo x interacción para trayectoria.\nCalcula el mapa de calor por tipo o reejecuta el análisis."},
    "Heatmap completo ligantes x resíduos": {"pt": "Mapa de calor completo ligantes x resíduos", "en": "Complete ligand x residue heatmap", "es": "Mapa de calor completo ligandos x residuos"},
    "FP - classes importantes": {"en": "FP - important classes", "es": "FP - clases importantes"},
    "FP - frequência por classe": {"en": "FP - frequency by class", "es": "FP - frecuencia por clase"},
    "FP - cobertura e importância": {"en": "FP - coverage and importance", "es": "FP - cobertura e importancia"},
    "FP - heatmap de importância": {"pt": "FP - mapa de calor de importância", "en": "FP - importance heatmap", "es": "FP - mapa de calor de importancia"},
    "FP - frequência de interações": {"en": "FP - interaction frequency", "es": "FP - frecuencia de interacciones"},
    "FP - interações prevalentes": {"en": "FP - prevalent interactions", "es": "FP - interacciones prevalentes"},
    "FP - heatmap de interações": {"pt": "FP - mapa de calor de interações", "en": "FP - interaction heatmap", "es": "FP - mapa de calor de interacciones"},
    "Has noncovalent interactions with the protein": {
        "pt": "Possui interações não covalentes com a proteína",
        "en": "Has noncovalent interactions with the protein",
        "es": "Tiene interacciones no covalentes con la proteína",
    },
    "Ligand's level 0 features only": {
        "pt": "Apenas features de nível 0 do ligante",
        "en": "Ligand's level 0 features only",
        "es": "Solo features de nivel 0 del ligando",
    },
    "Protein's level 0 features only": {
        "pt": "Apenas features de nível 0 da proteína",
        "en": "Protein's level 0 features only",
        "es": "Solo features de nivel 0 de la proteína",
    },
    "Upper level with ligand atomic information only": {
        "pt": "Nível superior apenas com informação atômica do ligante",
        "en": "Upper level with ligand atomic information only",
        "es": "Nivel superior solo con información atómica del ligando",
    },
    "Upper level with protein atomic information only": {
        "pt": "Nível superior apenas com informação atômica da proteína",
        "en": "Upper level with protein atomic information only",
        "es": "Nivel superior solo con información atómica de la proteína",
    },
    "Intraligand interactions only": {
        "pt": "Apenas interações intraligante",
        "en": "Intraligand interactions only",
        "es": "Solo interacciones intraligando",
    },
    "Intraprotein interactions only": {
        "pt": "Apenas interações intraproteína",
        "en": "Intraprotein interactions only",
        "es": "Solo interacciones intraproteína",
    },
    "Features with collision in the same complex": {
        "pt": "Features com colisão no mesmo complexo",
        "en": "Features with collision in the same complex",
        "es": "Features con colisión en el mismo complejo",
    },
    "Unreliable feature": {
        "pt": "Feature não confiável",
        "en": "Unreliable feature",
        "es": "Feature no confiable",
    },
    "Ligantes selecionados": {"en": "Selected ligands", "es": "Ligandos seleccionados"},
    "Matriz de similaridade": {"en": "Similarity matrix", "es": "Matriz de similitud"},
    "Otsu fallback": {"en": "Otsu fallback", "es": "Fallback Otsu"},
    "Clusters": {"en": "Clusters", "es": "Clusters"},
    "Tipo": {"en": "Type", "es": "Tipo"},
    "Total": {"en": "Total", "es": "Total"},
    "total": {"en": "total", "es": "total"},
    "Contagem de interações": {"en": "Interaction count", "es": "Conteo de interacciones"},
    "Resumo de interações": {"en": "Interaction summary", "es": "Resumen de interacciones"},
    "Mapa de interação por resíduo": {"en": "Residue interaction map", "es": "Mapa de interacción por residuo"},
    "Mapa de calor completo ligantes x resíduos": {"en": "Complete ligands x residues heatmap", "es": "Mapa de calor completo ligandos x residuos"},
    "Cadeia/Resíduo/Num": {"en": "Chain/Residue/Num", "es": "Cadena/Residuo/Núm"},
    "Top 30 resíduos com mais interações": {"en": "Top 30 residues with most interactions", "es": "Top 30 residuos con más interacciones"},
    "Dendrograma e matriz clusterizada": {"en": "Dendrogram and clustered matrix", "es": "Dendrograma y matriz clusterizada"},
    "Relatório LUNA": {"en": "LUNA report", "es": "Reporte LUNA"},
    (
        "Cada célula pode conter várias faixas de cor, uma por classe de interação. Esse gráfico mostra quando um mesmo resíduo participa de mecanismos químicos diferentes entre ligantes."
    ): {
        "en": "Each cell can contain several color bands, one for each interaction class. This chart shows when the same residue participates in different chemical mechanisms across ligands.",
        "es": "Cada celda puede contener varias franjas de color, una por clase de interacción. Este gráfico muestra cuándo un mismo residuo participa en mecanismos químicos diferentes entre ligandos.",
    },
    (
        "Resume as classes atribuídas às features que passaram pelo corte de p-value. Classes dominantes indicam a natureza estrutural das features que mais influenciam o modelo."
    ): {
        "en": "Summarizes the classes assigned to features that passed the p-value cutoff. Dominant classes indicate the structural nature of the features that most influence the model.",
        "es": "Resume las clases asignadas a las features que pasaron el corte de p-value. Las clases dominantes indican la naturaleza estructural de las features que más influyen en el modelo.",
    },
    (
        "Mostra a composição de classe das features importantes. Barras mistas sugerem features com colisões ou natureza ambígua; barras concentradas sugerem assinatura mais interpretável."
    ): {
        "en": "Shows the class composition of important features. Mixed bars suggest features with collisions or ambiguous nature; concentrated bars suggest a more interpretable signature.",
        "es": "Muestra la composición de clase de las features importantes. Barras mixtas sugieren features con colisiones o naturaleza ambigua; barras concentradas sugieren una firma más interpretable.",
    },
    (
        "Compara quantos ligantes possuem cada feature com a importância estimada pelo modelo. Features com alta importância e cobertura moderada costumam ser boas candidatas para inspeção."
    ): {
        "en": "Compares how many ligands contain each feature with the importance estimated by the model. Features with high importance and moderate coverage are usually good candidates for inspection.",
        "es": "Compara cuántos ligandos poseen cada feature con la importancia estimada por el modelo. Features con alta importancia y cobertura moderada suelen ser buenas candidatas para inspección.",
    },
    (
        "Relaciona features importantes e classes de fingerprint. Tons mais intensos indicam maior peso relativo; use os IDs das features para gerar sessão PyMOL correspondente."
    ): {
        "en": "Relates important features and fingerprint classes. More intense tones indicate higher relative weight; use feature IDs to generate the corresponding PyMOL session.",
        "es": "Relaciona features importantes y clases de fingerprint. Tonos más intensos indican mayor peso relativo; usa los IDs de las features para generar la sesión PyMOL correspondiente.",
    },
    (
        "Mostra a distribuição das interações prevalentes nos shells das features importantes. Ajuda a ligar uma feature abstrata a contatos químicos observáveis."
    ): {
        "en": "Shows the distribution of prevalent interactions in the shells of important features. It helps connect an abstract feature to observable chemical contacts.",
        "es": "Muestra la distribución de las interacciones prevalentes en los shells de las features importantes. Ayuda a conectar una feature abstracta con contactos químicos observables.",
    },
    (
        "Resume quais tipos de interação aparecem como dominantes nas features importantes, após aplicar o limiar configurado por z-score ou Otsu."
    ): {
        "en": "Summarizes which interaction types appear as dominant in the important features after applying the configured z-score or Otsu threshold.",
        "es": "Resume qué tipos de interacción aparecen como dominantes en las features importantes, después de aplicar el umbral configurado por z-score u Otsu.",
    },
    (
        "Cruza features importantes com interações prevalentes dos shells reais do LUNA. Ele revela se diferentes features importantes apontam para a mesma família de contatos."
    ): {
        "en": "Crosses important features with prevalent interactions from the real LUNA shells. It reveals whether different important features point to the same contact family.",
        "es": "Cruza features importantes con interacciones prevalentes de los shells reales de LUNA. Revela si diferentes features importantes apuntan a la misma familia de contactos.",
    },
    "Sem dados de interação para o heatmap completo": {"pt": "Sem dados de interação para o mapa de calor completo", "en": "No interaction data for the complete heatmap", "es": "Sin datos de interacción para el mapa de calor completo"},
    (
        "Linhas representam ligantes e colunas representam resíduos. A intensidade da célula indica quantas vezes o tipo de interação selecionado aparece naquele par; colunas densas destacam resíduos recorrentes."
    ): {
        "en": "Rows represent ligands and columns represent residues. Cell intensity indicates how often the selected interaction type appears in that pair; dense columns highlight recurrent residues.",
        "es": "Las filas representan ligandos y las columnas residuos. La intensidad de la celda indica cuántas veces el tipo de interacción seleccionado aparece en ese par; columnas densas destacan residuos recurrentes.",
    },
    (
        "Heatmap completo ligantes x resíduos: usa faixas coloridas para mostrar múltiplos tipos de interação no mesmo par ligante-resíduo. Ele é útil para enxergar complementaridade química: o mesmo resíduo pode estabilizar ligantes por mecanismos diferentes."
    ): {
        "pt": "Mapa de calor completo ligantes x resíduos: usa faixas coloridas para mostrar múltiplos tipos de interação no mesmo par ligante-resíduo. Ele é útil para enxergar complementaridade química: o mesmo resíduo pode estabilizar ligantes por mecanismos diferentes.",
        "en": "Complete ligand x residue heatmap: uses colored bands to show multiple interaction types in the same ligand-residue pair. It helps visualize chemical complementarity: the same residue can stabilize ligands through different mechanisms.",
        "es": "Mapa de calor completo ligandos x residuos: usa franjas de color para mostrar múltiples tipos de interacción en el mismo par ligando-residuo. Es útil para ver complementariedad química: el mismo residuo puede estabilizar ligandos por mecanismos diferentes.",
    },
    (
        "Clusters: reorganizam a matriz de similaridade para revelar famílias de ligantes por comportamento no sítio. Use os grupos como hipótese para priorização e para escolher representantes para inspeção no PyMOL."
    ): {
        "en": "Clusters: reorganize the similarity matrix to reveal ligand families by site behavior. Use the groups as prioritization hypotheses and to choose representatives for PyMOL inspection.",
        "es": "Clusters: reorganizan la matriz de similitud para revelar familias de ligandos por comportamiento en el sitio. Usa los grupos como hipótesis de priorización y para elegir representantes para inspección en PyMOL.",
    },
    (
        "Fingerprints e FP análises: cada feature resume uma vizinhança de interação. A classe atribuída descreve a natureza dominante da feature; a importância do modelo estima quanto ela ajuda a separar classes/rótulos ou valores de atividade."
    ): {
        "en": "Fingerprints and FP analyses: each feature summarizes an interaction neighborhood. The assigned class describes the feature's dominant nature; model importance estimates how much it helps separate classes/labels or activity values.",
        "es": "Fingerprints y análisis FP: cada feature resume una vecindad de interacción. La clase asignada describe la naturaleza dominante de la feature; la importancia del modelo estima cuánto ayuda a separar clases/etiquetas o valores de actividad.",
    },
    (
        "Importância e p-value: o z-score de importância compara uma feature contra a distribuição de importâncias do conjunto. O p-value é calculado pela equação de Keiser and Hert [1], p = 1 - exp(-exp(((-z*pi)/sqrt(6)) - 0.577215665))."
    ): {
        "en": "Importance and p-value: the importance z-score compares a feature against the importance distribution of the set. The p-value is calculated with the Keiser and Hert equation [1], p = 1 - exp(-exp(((-z*pi)/sqrt(6)) - 0.577215665)).",
        "es": "Importancia y p-value: el z-score de importancia compara una feature contra la distribución de importancias del conjunto. El p-value se calcula con la ecuación de Keiser and Hert [1], p = 1 - exp(-exp(((-z*pi)/sqrt(6)) - 0.577215665)).",
    },
    (
        "Otsu's Thresholding: quando nenhuma feature passa pelo critério z-score > 1, Otsu define um limiar alternativo baseado na separação da distribuição de percentuais. Isso evita aceitar apenas casos 100% prevalentes quando a base tem padrões intermediários."
    ): {
        "en": "Otsu's Thresholding: when no feature passes the z-score > 1 criterion, Otsu defines an alternative threshold based on the separation of the percentage distribution. This avoids accepting only 100% prevalent cases when the dataset has intermediate patterns.",
        "es": "Otsu's Thresholding: cuando ninguna feature pasa el criterio z-score > 1, Otsu define un umbral alternativo basado en la separación de la distribución de porcentajes. Esto evita aceptar solo casos 100% prevalentes cuando la base tiene patrones intermedios.",
    },
    (
        "Sessões PyMOL: permitem validar visualmente se as features ou filtros representam contatos plausíveis no complexo. A filtragem dinâmica pode gerar sessões novas ou, se o projeto salvo não reabrir, copiar sessões existentes compatíveis com a matriz cacheada."
    ): {
        "en": "PyMOL sessions: visually validate whether features or filters represent plausible contacts in the complex. Dynamic filtering can generate new sessions or, if the saved project cannot be reopened, copy existing sessions compatible with the cached matrix.",
        "es": "Sesiones PyMOL: permiten validar visualmente si las features o filtros representan contactos plausibles en el complejo. El filtrado dinámico puede generar sesiones nuevas o, si el proyecto guardado no se reabre, copiar sesiones existentes compatibles con la matriz cacheada.",
    },
    (
        "O dendrograma mostra a distância entre perfis de interação e a matriz reordenada evidencia grupos. Clusters compactos sugerem ligantes que compartilham padrões de contato e podem ser priorizados em conjunto."
    ): {
        "en": "The dendrogram shows distances between interaction profiles and the reordered matrix highlights groups. Compact clusters suggest ligands that share contact patterns and can be prioritized together.",
        "es": "El dendrograma muestra la distancia entre perfiles de interacción y la matriz reordenada evidencia grupos. Clusters compactos sugieren ligandos que comparten patrones de contacto y pueden priorizarse en conjunto.",
    },
    (
        "Cada linha resume uma base de fingerprints carregada. Se o modelo aparecer como fallback ou indisponível, a interpretação deve ser tratada como exploratória e a interface mostra a causa no campo de método."
    ): {
        "en": "Each row summarizes a loaded fingerprint dataset. If the model appears as fallback or unavailable, treat the interpretation as exploratory; the interface shows the cause in the method field.",
        "es": "Cada línea resume una base de fingerprints cargada. Si el modelo aparece como fallback o indisponible, la interpretación debe tratarse como exploratoria; la interfaz muestra la causa en el campo de método.",
    },
    "Tabela com os primeiros ligantes e seus clusters hierárquicos.": {
        "en": "Table with the first ligands and their hierarchical clusters.",
        "es": "Tabla con los primeros ligandos y sus clusters jerárquicos.",
    },

    # Report text.
    "Relatório HIP²LInterActomics": {"en": "HIP²LInterActomics report", "es": "Reporte HIP²LInterActomics"},
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
    "Rótulos FP": {"en": "FP labels", "es": "Etiquetas FP"},
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
        "pt": "Mapa de calor por tipo: cruza ligantes e resíduos para uma interação escolhida. Células mais intensas indicam mais ocorrências daquele contato. Colunas densas sugerem resíduos-chave; linhas densas sugerem ligantes com muitos contatos daquele tipo.",
        "en": "Heatmap by type: crosses ligands and residues for a selected interaction. More intense cells indicate more occurrences of that contact. Dense columns suggest key residues; dense rows suggest ligands with many contacts of that type.",
        "es": "Mapa de calor por tipo: cruza ligandos y residuos para una interacción elegida. Celdas más intensas indican más ocurrencias de ese contacto. Columnas densas sugieren residuos clave; filas densas sugieren ligandos con muchos contactos de ese tipo.",
    },
    (
        "Matriz de similaridade: compara ligantes pelos fingerprints de interação. Valores próximos de 1 indicam perfis de interação semelhantes; valores baixos indicam modos de interação distintos, mesmo quando as moléculas parecem estruturalmente parecidas."
    ): {
        "en": "Similarity matrix: compares ligands by interaction fingerprints. Values close to 1 indicate similar interaction profiles; low values indicate distinct interaction modes, even when molecules look structurally similar.",
        "es": "Matriz de similitud: compara ligandos por fingerprints de interacción. Valores cercanos a 1 indican perfiles de interacción semejantes; valores bajos indican modos de interacción distintos, incluso cuando las moléculas parecen estructuralmente parecidas.",
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
        "es": "Cada celda compara dos ligandos por fingerprint de interacción. Tonos más intensos indican mayor similitud; bloques a lo largo de la diagonal sugieren familias con modos de interacción semejantes.",
    },

    # Standard interaction labels used in plots.  Data keys remain unchanged.
    "Hydrogen bond": {"pt": "Ligação de hidrogênio", "en": "Hydrogen bond", "es": "Puente de hidrógeno"},
    "Weak hydrogen bond": {"pt": "Ligação de hidrogênio fraca", "en": "Weak hydrogen bond", "es": "Puente de hidrógeno débil"},
    "Hydrophobic": {"pt": "Hidrofóbica", "en": "Hydrophobic", "es": "Hidrofóbica"},
    "Halogen bond": {"pt": "Ligação de halogênio", "en": "Halogen bond", "es": "Enlace de halógeno"},
    "Halogen-pi": {"pt": "Halogênio-pi", "en": "Halogen-pi", "es": "Halógeno-pi"},
    "Chalcogen bond": {"pt": "Ligação de calcogênio", "en": "Chalcogen bond", "es": "Enlace de calcógeno"},
    "Chalcogen-pi": {"pt": "Calcogênio-pi", "en": "Chalcogen-pi", "es": "Calcógeno-pi"},
    "Ionic": {"pt": "Iônica", "en": "Ionic", "es": "Iónica"},
    "Salt bridge": {"pt": "Ponte salina", "en": "Salt bridge", "es": "Puente salino"},
    "Cation-pi": {"pt": "Cátion-pi", "en": "Cation-pi", "es": "Catión-pi"},
    "Cation-nucleophile": {"pt": "Cátion-nucleófilo", "en": "Cation-nucleophile", "es": "Catión-nucleófilo"},
    "Anion-electrophile": {"pt": "Ânion-eletrófilo", "en": "Anion-electrophile", "es": "Anión-electrófilo"},
    "Anion-pi": {"pt": "Ânion-pi", "en": "Anion-pi", "es": "Anión-pi"},
    "Pi-stacking": {"pt": "Empilhamento pi", "en": "Pi-stacking", "es": "Apilamiento pi"},
    "Aromatic stacking": {"pt": "Empilhamento aromático", "en": "Aromatic stacking", "es": "Apilamiento aromático"},
    "Face-to-face": {"pt": "Face-a-face", "en": "Face-to-face", "es": "Cara a cara"},
    "Edge-to-face": {"pt": "Borda-a-face", "en": "Edge-to-face", "es": "Borde a cara"},
    "Face-to-edge pi-stacking": {"pt": "Empilhamento pi face-borda", "en": "Face-to-edge pi-stacking", "es": "Apilamiento pi cara-borde"},
    "Face-to-face pi-stacking": {"pt": "Empilhamento pi face-face", "en": "Face-to-face pi-stacking", "es": "Apilamiento pi cara-cara"},
    "Face-to-slope pi-stacking": {"pt": "Empilhamento pi face-inclinado", "en": "Face-to-slope pi-stacking", "es": "Apilamiento pi cara-inclinado"},
    "Displaced face-to-edge pi-stacking": {"pt": "Empilhamento pi deslocado face-borda", "en": "Displaced face-to-edge pi-stacking", "es": "Apilamiento pi desplazado cara-borde"},
    "Displaced face-to-face pi-stacking": {"pt": "Empilhamento pi deslocado face-face", "en": "Displaced face-to-face pi-stacking", "es": "Apilamiento pi desplazado cara-cara"},
    "Displaced face-to-slope pi-stacking": {"pt": "Empilhamento pi deslocado face-inclinado", "en": "Displaced face-to-slope pi-stacking", "es": "Apilamiento pi desplazado cara-inclinado"},
    "Parallel": {"pt": "Paralelo", "en": "Parallel", "es": "Paralelo"},
    "Parallel multipolar": {"pt": "Multipolar paralelo", "en": "Parallel multipolar", "es": "Multipolar paralelo"},
    "Antiparallel multipolar": {"pt": "Multipolar antiparalelo", "en": "Antiparallel multipolar", "es": "Multipolar antiparalelo"},
    "Orthogonal multipolar": {"pt": "Multipolar ortogonal", "en": "Orthogonal multipolar", "es": "Multipolar ortogonal"},
    "Tilted multipolar": {"pt": "Multipolar inclinado", "en": "Tilted multipolar", "es": "Multipolar inclinado"},
    "T-shaped": {"pt": "Forma de T", "en": "T-shaped", "es": "Forma de T"},
    "Amide-aromatic stacking": {"pt": "Empilhamento amida-aromático", "en": "Amide-aromatic stacking", "es": "Apilamiento amida-aromático"},
    "Charge-dipole interaction": {"pt": "Interação carga-dipolo", "en": "Charge-dipole interaction", "es": "Interacción carga-dipolo"},
    "Unfavorable charge-dipole interaction": {"pt": "Interação carga-dipolo desfavorável", "en": "Unfavorable charge-dipole interaction", "es": "Interacción carga-dipolo desfavorable"},
    "Unfavorable dipole interaction": {"pt": "Interação dipolar desfavorável", "en": "Unfavorable dipole interaction", "es": "Interacción dipolar desfavorable"},
    "Water-bridged hydrogen bond": {"pt": "Ponte de hidrogênio mediada por água", "en": "Water-bridged hydrogen bond", "es": "Puente de hidrógeno mediado por agua"},
    "Disulfide bond": {"pt": "Ponte dissulfeto", "en": "Disulfide bond", "es": "Puente disulfuro"},
    "Metal coordination": {"pt": "Coordenação metálica", "en": "Metal coordination", "es": "Coordinación metálica"},
    "Van der Waals": {"pt": "Van der Waals", "en": "Van der Waals", "es": "Van der Waals"},
    "Proximal": {"pt": "Proximal", "en": "Proximal", "es": "Proximal"},
    "Multipolar interaction": {"pt": "Interação multipolar", "en": "Multipolar interaction", "es": "Interacción multipolar"},
    "Multiple interactions": {"pt": "Interações múltiplas", "en": "Multiple interactions", "es": "Interacciones múltiples"},
    "Repulsive": {"pt": "Repulsiva", "en": "Repulsive", "es": "Repulsiva"},
    "Unfavorable": {"pt": "Desfavorável", "en": "Unfavorable", "es": "Desfavorable"},
    "Unfavorable anion-nucleophile": {"pt": "Ânion-nucleófilo desfavorável", "en": "Unfavorable anion-nucleophile", "es": "Anión-nucleófilo desfavorable"},
    "Unfavorable cation-electrophile": {"pt": "Cátion-eletrófilo desfavorável", "en": "Unfavorable cation-electrophile", "es": "Catión-electrófilo desfavorable"},
    "Unfavorable electrophile-electrophile": {"pt": "Eletrófilo-eletrófilo desfavorável", "en": "Unfavorable electrophile-electrophile", "es": "Electrófilo-electrófilo desfavorable"},
    "Unfavorable nucleophile-nucleophile": {"pt": "Nucleófilo-nucleófilo desfavorável", "en": "Unfavorable nucleophile-nucleophile", "es": "Nucleófilo-nucleófilo desfavorable"},
    "Aparência": {"en": "Appearance", "es": "Apariencia"},
    "Tema do sistema": {"en": "System theme", "es": "Tema del sistema"},
    "Claro": {"en": "Light", "es": "Claro"},
    "Escuro": {"en": "Dark", "es": "Oscuro"},
    "Tema atualizado": {"en": "Theme updated", "es": "Tema actualizado"},
    "Resultados sob demanda": {"en": "On-demand results", "es": "Resultados bajo demanda"},
    "Abrir resultados": {"en": "Open results", "es": "Abrir resultados"},
    (
        "Os módulos científicos e gráficos serão carregados somente quando necessários, "
        "reduzindo a memória e acelerando a abertura do aplicativo."
    ): {
        "en": "Scientific and charting modules load only when needed, reducing memory use and speeding up application startup.",
        "es": "Los módulos científicos y gráficos se cargan solo cuando son necesarios, reduciendo la memoria y acelerando el inicio de la aplicación.",
    },
    "Projeto aberto, mas a seleção de ligantes não pôde ser restaurada": {
        "en": "The project opened, but the ligand selection could not be restored",
        "es": "El proyecto se abrió, pero no se pudo restaurar la selección de ligandos",
    },
    "Projeto aberto, mas os resultados não puderam ser carregados": {
        "en": "The project opened, but its results could not be loaded",
        "es": "El proyecto se abrió, pero no se pudieron cargar sus resultados",
    },
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
            re.compile(r"^Sim, pH (.+)$"),
            lambda m: (
                f"Yes, pH {m.group(1)}"
                if lang == "en"
                else f"Sí, pH {m.group(1)}"
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
            re.compile(r"^Interações por tipo - (.+)$"),
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
            re.compile(r"^Interações por aminoácido - ligante (.+)$"),
            lambda m: (
                f"Interactions by amino acid - ligand {m.group(1)}"
                if lang == "en"
                else f"Interacciones por aminoácido - ligando {m.group(1)}"
            ),
        ),
        (
            re.compile(r"^Interações por aminoácido - frame (.+)$"),
            lambda m: (
                f"Interactions by amino acid - frame {m.group(1)}"
                if lang == "en"
                else f"Interacciones por aminoácido - frame {m.group(1)}"
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
            re.compile(r"^Sem ocorrências de '(.+)'$"),
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
        setattr(Text, "_luna_i18n_original_set_text", original_set_text)

        def set_text(self, value):
            if isinstance(value, str):
                previous_source = getattr(self, "_luna_i18n_source_text", None)
                source = source_text(value)
                if previous_source and value == t(previous_source):
                    source = previous_source
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
            original_set_text = getattr(type(artist), "_luna_i18n_original_set_text", None)
            if original_set_text is not None:
                original_set_text(artist, new_text)
            else:
                artist.set_text(new_text)
        except Exception:
            pass


def translate_figure(fig) -> None:
    if fig is None:
        return
    try:
        for text in list(getattr(fig, "texts", []) or []):
            _translate_text_artist(text)
        for legend in list(getattr(fig, "legends", []) or []):
            for text in legend.get_texts():
                _translate_text_artist(text)
            title = legend.get_title()
            if title is not None:
                _translate_text_artist(title)
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
        )
    except Exception:
        return

    widgets = [root] + list(root.findChildren(QObject))
    for widget in widgets:
        if isinstance(widget, QLabel):
            try:
                has_pixmap = widget.pixmap() is not None
            except Exception:
                has_pixmap = False
            if not (has_pixmap and not widget.text()):
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
