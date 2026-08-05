"""HTML and PDF reports for completed HIP2LInterActomics projects."""
from __future__ import annotations

import base64
import csv
import gc
import html
import json
import mimetypes
import multiprocessing
import os
import pickle
import shutil
import tempfile
import textwrap
import time
import traceback
from datetime import datetime
from pathlib import Path
from typing import Callable

from jinja2 import Environment, PackageLoader, select_autoescape

from .plot_manifest import load_manifest, resolve_plot_path
from .project import ProjectConfig
from ..i18n import language as active_language, set_language, t


_A4_LANDSCAPE = (11.69, 8.27)
_PDF_DPI = 300
_IMAGE_SUFFIXES = {".png", ".jpg", ".jpeg", ".tif", ".tiff", ".bmp"}
_REPORT_TEMP_PREFIXES = ("_report_", "_report_pdf_")
_IFP_ORDER = {"EIFP": 0, "FIFP": 1, "HIFP": 2}
_MODEL_ORDER = {"extra_trees": 0, "gradient_boosting": 1}
_REPORT_CATEGORY_ORDER = {
    "distribution": 10,
    "interaction_heatmap": 20,
    "complete_heatmap": 30,
    "similarity": 40,
    "clusters": 50,
    "fingerprint": 100,
    "appendix": 10000,
}

_REPORT_LANGUAGES = {"en", "pt", "es"}
_REPORT_COPY: dict[str, dict[str, str]] = {
    "en": {
        "report_title": "HIP2LInterActomics report",
        "generated": "Generated on",
        "summary": "Summary",
        "protein": "Protein",
        "ligands": "Ligands",
        "selected_total": "Total selected",
        "workdir": "Workdir",
        "processed_entries": "Processed entries",
        "configuration": "Configuration",
        "parameter": "Parameter",
        "value": "Value",
        "interaction_count": "Interaction count by type",
        "type": "Type",
        "total": "Total",
        "top_residues": "Top 30 residues with most interactions",
        "residue": "Chain/residue/number",
        "count": "Count",
        "source": "Source",
        "cluster_assignment": "Cluster assignment",
        "cluster": "Cluster",
        "ligand": "Ligand",
        "fp_interpret": "How to interpret fingerprint analyses",
        "fp_columns": "FP Analyses column guide",
        "fp_summary": "Fingerprint analyses summary",
        "column": "Column",
        "interpretation": "Interpretation",
        "top_features": "Top 50 features",
        "rank": "Rank",
        "feature": "Feature",
        "assigned_level": "Assigned level",
        "assigned_class": "Assigned class",
        "coverage": "Coverage (%)",
        "importance": "Importance",
        "appendix": "Appendix",
        "yes": "Yes",
        "no": "No",
        "not_reported": "not reported",
        "default_luna": "LUNA default",
        "no_interactions": "No counted interactions.",
        "no_residues": "No counted residues.",
        "analysis_interpretation": "How to interpret the analyses",
        "analysis_intro": "This report brings together the project parameters, numerical summary, and all plots available for the active report language. Visual patterns are hypotheses for structural and experimental validation.",
        "frequent_interactions": "Most frequent interactions",
        "frequent_residues": "Most frequent residues",
        "fp_guide_intro": "Use this guide when interpreting fingerprint tables, rankings, and heatmaps.",
        "fp_summary_intro": "Each row summarizes a fingerprint dataset, its feature count, and the selection threshold.",
        "top_features_intro": "Independent feature ranking for comparison between the two ensemble methods.",
        "cluster_table_intro": "Table of ligands and their hierarchical groups.",
        "image_unavailable": "Unavailable plot",
        "image_missing": "Image not found",
        "image_skipped": "The file could not be read and was skipped.",
        "fp_features": "features",
        "fp_important": "important",
        "fp_model": "model",
        "fp_threshold": "threshold",
        "cfg_selected_ligands": "Selected ligands",
        "cfg_include_waters": "Include waters",
        "cfg_add_h": "Add H",
        "cfg_count_fingerprint": "Count fingerprint",
        "cfg_similarity": "Similarity matrix",
        "cfg_pymol": "PyMOL sessions",
        "cfg_binding_modes": "Binding modes filter",
        "cfg_interactions": "Interaction config",
        "cfg_max_distance": "Global maximum distance",
        "cfg_fp_labels": "FP labels",
        "cfg_fp_task": "FP task",
        "cfg_otsu": "Otsu fallback",
        "cfg_cores": "Cores",
        "cfg_levels": "levels",
        "cfg_radius": "radius",
        "cfg_length": "length",
        "image_similarity_title": "Similarity matrix",
        "image_similarity_caption": "Each cell compares two ligands by their interaction fingerprints. Higher values indicate more similar profiles; contiguous blocks suggest families sharing molecular recognition modes.",
        "image_cluster_title": "Ligand clustering",
        "image_cluster_caption": "The dendrogram and reordered matrix group ligands by proximity between interaction profiles. Compact groups help select representatives and identify series with convergent binding-site behavior.",
        "image_heatmap_title": "Interaction heatmap",
        "image_heatmap_caption": "The heatmap summarizes interaction frequency or intensity between ligands, residues, and contact types. Stronger regions reveal recurring patterns and potential key residues for recognition.",
        "image_distribution_title": "Interaction distribution",
        "image_distribution_caption": "The distribution shows the relative abundance of contact types in the analyzed set. Dominant bars or regions indicate recurring chemical forces and help compare contact density and diversity.",
        "image_network_title": "Interaction network",
        "image_network_caption": "The network represents molecular entities as nodes and interactions as edges. Connectivity, communities, and central nodes help identify residues or ligands that organize the global contact pattern.",
        "image_fingerprint_title": "Fingerprint analysis",
        "image_fingerprint_caption": "The plot summarizes interaction-fingerprint features, their frequency, class, or importance. Recurrent features describe shared structural signatures; discriminant features help distinguish groups.",
        "image_model_title": "Model performance",
        "image_model_caption": "This visualization summarizes model behavior or selected variables. Interpret separation, error, coverage, and stability together before prioritizing a molecular hypothesis.",
        "image_ligand_title": "Ligand representation",
        "image_ligand_caption": "The representation relates ligand atom numbering to the descriptors and contacts used in the analyses, locating observed patterns on the molecular scaffold.",
        "image_generic_title": "Result plot",
        "image_generic_caption": "This visualization complements the project results. Interpret scales, labels, and legend together with the tables and structural inspection.",
    },
    "pt": {
        "report_title": "Relat\u00f3rio HIP2LInterActomics", "generated": "Gerado em", "summary": "Resumo", "protein": "Prote\u00edna", "ligands": "Ligantes", "selected_total": "Total selecionado", "workdir": "Diret\u00f3rio de trabalho", "processed_entries": "Entradas processadas", "configuration": "Configura\u00e7\u00e3o", "parameter": "Par\u00e2metro", "value": "Valor", "interaction_count": "Contagem por tipo de intera\u00e7\u00e3o", "type": "Tipo", "total": "Total", "top_residues": "Top 30 res\u00edduos com mais intera\u00e7\u00f5es", "residue": "Cadeia/res\u00edduo/n\u00famero", "count": "Contagem", "source": "Fonte", "cluster_assignment": "Atribui\u00e7\u00e3o de clusters", "cluster": "Cluster", "ligand": "Ligante", "fp_interpret": "Como interpretar as an\u00e1lises de fingerprints", "fp_columns": "Guia das colunas de An\u00e1lises FP", "fp_summary": "Resumo das an\u00e1lises de fingerprints", "column": "Coluna", "interpretation": "Interpreta\u00e7\u00e3o", "top_features": "Top 50 features", "rank": "Posi\u00e7\u00e3o", "feature": "Feature", "assigned_level": "N\u00edvel atribu\u00eddo", "assigned_class": "Classe atribu\u00edda", "coverage": "Cobertura (%)", "importance": "Import\u00e2ncia", "appendix": "Ap\u00eandice", "yes": "Sim", "no": "N\u00e3o", "not_reported": "n\u00e3o informado", "default_luna": "Padr\u00e3o LUNA", "no_interactions": "Sem intera\u00e7\u00f5es contabilizadas.", "no_residues": "Sem res\u00edduos contabilizados.", "analysis_interpretation": "Como interpretar as an\u00e1lises", "analysis_intro": "Este relat\u00f3rio re\u00fane os par\u00e2metros do projeto, o resumo num\u00e9rico e todos os gr\u00e1ficos dispon\u00edveis para o idioma ativo do relat\u00f3rio. Padr\u00f5es visuais s\u00e3o hip\u00f3teses para valida\u00e7\u00e3o estrutural e experimental.", "frequent_interactions": "Intera\u00e7\u00f5es mais frequentes", "frequent_residues": "Res\u00edduos mais frequentes", "fp_guide_intro": "Use este guia ao interpretar tabelas, rankings e mapas de calor de fingerprints.", "fp_summary_intro": "Cada linha resume uma base de fingerprints, a quantidade de features e o limiar de sele\u00e7\u00e3o.", "top_features_intro": "Ranking independente das features para compara\u00e7\u00e3o entre os dois m\u00e9todos de ensemble.", "cluster_table_intro": "Tabela dos ligantes e seus grupos hier\u00e1rquicos.", "image_unavailable": "Gr\u00e1fico indispon\u00edvel", "image_missing": "Imagem n\u00e3o encontrada", "image_skipped": "O arquivo n\u00e3o p\u00f4de ser lido e foi ignorado.", "fp_features": "features", "fp_important": "importantes", "fp_model": "modelo", "fp_threshold": "limiar", "cfg_selected_ligands": "Ligantes selecionados", "cfg_include_waters": "Incluir \u00e1guas", "cfg_add_h": "Adicionar H", "cfg_count_fingerprint": "Fingerprint de contagem", "cfg_similarity": "Matriz de similaridade", "cfg_pymol": "Sess\u00f5es PyMOL", "cfg_binding_modes": "Filtro de modos de liga\u00e7\u00e3o", "cfg_interactions": "Configura\u00e7\u00e3o de intera\u00e7\u00f5es", "cfg_max_distance": "Dist\u00e2ncia m\u00e1xima global", "cfg_fp_labels": "R\u00f3tulos FP", "cfg_fp_task": "Tarefa FP", "cfg_otsu": "Fallback Otsu", "cfg_cores": "N\u00facleos", "cfg_levels": "n\u00edveis", "cfg_radius": "raio", "cfg_length": "tamanho",
        "image_similarity_title": "Matriz de similaridade", "image_similarity_caption": "Cada c\u00e9lula compara dois ligantes por seus fingerprints de intera\u00e7\u00e3o. Valores mais altos indicam perfis mais semelhantes; blocos cont\u00edguos sugerem fam\u00edlias que compartilham modos de reconhecimento molecular.", "image_cluster_title": "Agrupamento de ligantes", "image_cluster_caption": "O dendrograma e a matriz reordenada agrupam ligantes pela proximidade entre perfis de intera\u00e7\u00e3o. Grupos compactos ajudam a selecionar representantes e identificar s\u00e9ries com comportamento convergente no s\u00edtio.", "image_heatmap_title": "Mapa de calor de intera\u00e7\u00f5es", "image_heatmap_caption": "O mapa de calor resume a frequ\u00eancia ou intensidade das intera\u00e7\u00f5es entre ligantes, res\u00edduos e tipos de contato. Regi\u00f5es mais intensas revelam padr\u00f5es recorrentes e poss\u00edveis res\u00edduos-chave para o reconhecimento.", "image_distribution_title": "Distribui\u00e7\u00e3o de intera\u00e7\u00f5es", "image_distribution_caption": "A distribui\u00e7\u00e3o mostra a abund\u00e2ncia relativa dos tipos de contato no conjunto analisado. Barras ou regi\u00f5es dominantes indicam as for\u00e7as qu\u00edmicas mais recorrentes e ajudam a comparar densidade e diversidade de contatos.", "image_network_title": "Rede de intera\u00e7\u00f5es", "image_network_caption": "A rede representa entidades moleculares como n\u00f3s e suas intera\u00e7\u00f5es como arestas. Conectividade, comunidades e n\u00f3s centrais ajudam a reconhecer res\u00edduos ou ligantes que organizam o padr\u00e3o global de contatos.", "image_fingerprint_title": "An\u00e1lise de fingerprints", "image_fingerprint_caption": "O gr\u00e1fico resume features dos fingerprints de intera\u00e7\u00e3o, sua frequ\u00eancia, classe ou import\u00e2ncia. Features recorrentes descrevem assinaturas estruturais compartilhadas; features discriminantes ajudam a separar grupos.", "image_model_title": "Desempenho do modelo", "image_model_caption": "Esta visualiza\u00e7\u00e3o resume o comportamento do modelo ou das vari\u00e1veis selecionadas. A interpreta\u00e7\u00e3o deve considerar separa\u00e7\u00e3o, erro, cobertura e estabilidade antes de priorizar uma hip\u00f3tese molecular.", "image_ligand_title": "Representa\u00e7\u00e3o do ligante", "image_ligand_caption": "A representa\u00e7\u00e3o relaciona a numera\u00e7\u00e3o at\u00f4mica do ligante aos descritores e contatos usados nas an\u00e1lises, localizando os padr\u00f5es observados no esqueleto molecular.", "image_generic_title": "Gr\u00e1fico de resultados", "image_generic_caption": "Esta visualiza\u00e7\u00e3o complementa os resultados do projeto. Interprete escalas, r\u00f3tulos e legenda em conjunto com as tabelas e a inspe\u00e7\u00e3o estrutural.",
    },
    "es": {
        "report_title": "Reporte HIP2LInterActomics", "generated": "Generado el", "summary": "Resumen", "protein": "Prote\u00edna", "ligands": "Ligandos", "selected_total": "Total seleccionado", "workdir": "Directorio de trabajo", "processed_entries": "Entradas procesadas", "configuration": "Configuraci\u00f3n", "parameter": "Par\u00e1metro", "value": "Valor", "interaction_count": "Conteo por tipo de interacci\u00f3n", "type": "Tipo", "total": "Total", "top_residues": "Top 30 residuos con m\u00e1s interacciones", "residue": "Cadena/residuo/n\u00famero", "count": "Conteo", "source": "Fuente", "cluster_assignment": "Asignaci\u00f3n de clusters", "cluster": "Cluster", "ligand": "Ligando", "fp_interpret": "C\u00f3mo interpretar los an\u00e1lisis de fingerprints", "fp_columns": "Gu\u00eda de columnas de An\u00e1lisis FP", "fp_summary": "Resumen de an\u00e1lisis de fingerprints", "column": "Columna", "interpretation": "Interpretaci\u00f3n", "top_features": "Top 50 features", "rank": "Posici\u00f3n", "feature": "Feature", "assigned_level": "Nivel asignado", "assigned_class": "Clase asignada", "coverage": "Cobertura (%)", "importance": "Importancia", "appendix": "Ap\u00e9ndice", "yes": "S\u00ed", "no": "No", "not_reported": "no informado", "default_luna": "Predeterminado de LUNA", "no_interactions": "No hay interacciones contabilizadas.", "no_residues": "No hay residuos contabilizadas.", "analysis_interpretation": "C\u00f3mo interpretar los an\u00e1lisis", "analysis_intro": "Este reporte re\u00fane los par\u00e1metros del proyecto, el resumen num\u00e9rico y todos los gr\u00e1ficos disponibles para el idioma activo del reporte. Los patrones visuales son hip\u00f3tesis para validaci\u00f3n estructural y experimental.", "frequent_interactions": "Interacciones m\u00e1s frecuentes", "frequent_residues": "Residuos m\u00e1s frecuentes", "fp_guide_intro": "Use esta gu\u00eda al interpretar tablas, rankings y mapas de calor de fingerprints.", "fp_summary_intro": "Cada fila resume una base de fingerprints, la cantidad de features y el umbral de selecci\u00f3n.", "top_features_intro": "Ranking independiente de features para comparar los dos m\u00e9todos de ensemble.", "cluster_table_intro": "Tabla de ligandos y sus grupos jer\u00e1rquicos.", "image_unavailable": "Gr\u00e1fico no disponible", "image_missing": "Imagen no encontrada", "image_skipped": "El archivo no se pudo leer y fue omitido.", "fp_features": "features", "fp_important": "importantes", "fp_model": "modelo", "fp_threshold": "umbral", "cfg_selected_ligands": "Ligandos seleccionados", "cfg_include_waters": "Incluir aguas", "cfg_add_h": "Agregar H", "cfg_count_fingerprint": "Fingerprint de conteo", "cfg_similarity": "Matriz de similitud", "cfg_pymol": "Sesiones PyMOL", "cfg_binding_modes": "Filtro de modos de uni\u00f3n", "cfg_interactions": "Configuraci\u00f3n de interacciones", "cfg_max_distance": "Distancia m\u00e1xima global", "cfg_fp_labels": "Etiquetas FP", "cfg_fp_task": "Tarea FP", "cfg_otsu": "Fallback Otsu", "cfg_cores": "N\u00facleos", "cfg_levels": "niveles", "cfg_radius": "radio", "cfg_length": "longitud",
        "image_similarity_title": "Matriz de similitud", "image_similarity_caption": "Cada celda compara dos ligandos por sus fingerprints de interacci\u00f3n. Valores m\u00e1s altos indican perfiles m\u00e1s similares; bloques contiguos sugieren familias que comparten modos de reconocimiento molecular.", "image_cluster_title": "Agrupamiento de ligandos", "image_cluster_caption": "El dendrograma y la matriz reordenada agrupan ligandos por proximidad entre perfiles de interacci\u00f3n. Los grupos compactos ayudan a seleccionar representantes e identificar series con comportamiento convergente en el sitio.", "image_heatmap_title": "Mapa de calor de interacciones", "image_heatmap_caption": "El mapa de calor resume la frecuencia o intensidad de las interacciones entre ligandos, residuos y tipos de contacto. Las regiones m\u00e1s intensas revelan patrones recurrentes y posibles residuos clave para el reconocimiento.", "image_distribution_title": "Distribuci\u00f3n de interacciones", "image_distribution_caption": "La distribuci\u00f3n muestra la abundancia relativa de los tipos de contacto en el conjunto analizado. Barras o regiones dominantes indican las fuerzas qu\u00edmicas m\u00e1s recurrentes y ayudan a comparar densidad y diversidad de contactos.", "image_network_title": "Red de interacciones", "image_network_caption": "La red representa entidades moleculares como nodos y sus interacciones como aristas. La conectividad, las comunidades y los nodos centrales ayudan a reconocer residuos o ligandos que organizan el patr\u00f3n global de contactos.", "image_fingerprint_title": "An\u00e1lisis de fingerprints", "image_fingerprint_caption": "El gr\u00e1fico resume features de los fingerprints de interacci\u00f3n, su frecuencia, clase o importancia. Las features recurrentes describen firmas estructurales compartidas; las discriminantes ayudan a separar grupos.", "image_model_title": "Desempe\u00f1o del modelo", "image_model_caption": "Esta visualizaci\u00f3n resume el comportamiento del modelo o de las variables seleccionadas. La interpretaci\u00f3n debe considerar separaci\u00f3n, error, cobertura y estabilidad antes de priorizar una hip\u00f3tesis molecular.", "image_ligand_title": "Representaci\u00f3n del ligando", "image_ligand_caption": "La representaci\u00f3n relaciona la numeraci\u00f3n at\u00f3mica del ligando con los descriptores y contactos usados en los an\u00e1lisis, localizando los patrones observados en el esqueleto molecular.", "image_generic_title": "Gr\u00e1fico de resultados", "image_generic_caption": "Esta visualizaci\u00f3n complementa los resultados del proyecto. Interprete escalas, etiquetas y leyenda junto con las tablas y la inspecci\u00f3n estructural.",
    },
}
_REPORT_COPY["es"]["no_residues"] = "No hay residuos contabilizados."

# Section titles are intentionally kept in the report catalog instead of the
# GUI catalog: report generation also runs from the terminal and PDF workers.
_REPORT_COPY["en"].update({
    "interaction_distribution": "Interaction distribution",
    "interaction_summary": "Interaction summary",
    "amino_acid_distribution": "Interactions by amino acid across all ligands",
    "ligand_atom_distribution": "Interactions by ligand atoms",
    "interaction_heatmaps": "Interaction heatmaps",
    "heatmaps_by_type": "Heatmaps by interaction type",
    "complete_interaction_heatmap": "Complete interaction heatmap",
    "fingerprint_analyses": "Fingerprint analyses",
    "fingerprint_definitions": "EIFP records interaction environments; FIFP records interaction features; HIFP records hybrid interaction features. Each fingerprint describes the molecular recognition pattern at a different structural granularity.",
    "similarity_matrix": "Similarity matrix",
    "hierarchical_clustering": "Hierarchical clustering",
    "matrix_reordered": "Matrix reordered by cluster",
    "brief_supervised_learning": "Brief supervised learning",
    "supervised_intro": "The supervised-learning panels rank fingerprint features with Extra Trees and Gradient Boosting. Interpret importance together with coverage, assigned class, interaction assignment, and the underlying structural evidence.",
    "report_appendix": "Additional available plots",
    "fp_top50_plot": "FP top 50 importance",
    "fp_class_summary_plot": "FP class summary",
    "fp_class_assignment_plot": "FP class assignment",
    "fp_coverage_importance_plot": "FP coverage importance",
    "fp_feature_presence_heatmap_plot": "FP feature presence heatmap",
    "fp_interaction_assignment_plot": "FP interaction assignment",
    "fp_prevalent_interactions_plot": "FP prevalent interactions",
    "fp_prevalent_interactions_heatmap_plot": "FP prevalent interactions heatmap",
})
_REPORT_COPY["pt"].update({
    "interaction_distribution": "Distribui\u00e7\u00e3o de intera\u00e7\u00f5es",
    "interaction_summary": "Resumo de intera\u00e7\u00f5es",
    "amino_acid_distribution": "Intera\u00e7\u00f5es por amino\u00e1cido em todos os ligantes",
    "ligand_atom_distribution": "Intera\u00e7\u00f5es por \u00e1tomos do ligante",
    "interaction_heatmaps": "Mapas de calor de intera\u00e7\u00f5es",
    "heatmaps_by_type": "Mapas de calor por tipo de intera\u00e7\u00e3o",
    "complete_interaction_heatmap": "Mapa de calor completo de intera\u00e7\u00f5es",
    "fingerprint_analyses": "An\u00e1lises de fingerprints",
    "fingerprint_definitions": "EIFP registra ambientes de intera\u00e7\u00e3o; FIFP registra features de intera\u00e7\u00e3o; HIFP registra features h\u00edbridas de intera\u00e7\u00e3o. Cada fingerprint descreve o padr\u00e3o de reconhecimento molecular em uma granularidade estrutural diferente.",
    "similarity_matrix": "Matriz de similaridade",
    "hierarchical_clustering": "Agrupamento hier\u00e1rquico",
    "matrix_reordered": "Matriz reordenada por cluster",
    "brief_supervised_learning": "Breve aprendizado supervisionado",
    "supervised_intro": "Os pain\u00e9is de aprendizado supervisionado classificam features de fingerprints com Extra Trees e Gradient Boosting. Interprete a import\u00e2ncia junto com cobertura, classe atribu\u00edda, atribui\u00e7\u00e3o de intera\u00e7\u00f5es e a evid\u00eancia estrutural subjacente.",
    "report_appendix": "Gr\u00e1ficos adicionais dispon\u00edveis",
    "fp_top50_plot": "Import\u00e2ncia das top 50 features de FP",
    "fp_class_summary_plot": "Resumo de classes de FP",
    "fp_class_assignment_plot": "Atribui\u00e7\u00e3o de classes de FP",
    "fp_coverage_importance_plot": "Cobertura e import\u00e2ncia de FP",
    "fp_feature_presence_heatmap_plot": "Mapa de calor de presen\u00e7a de features de FP",
    "fp_interaction_assignment_plot": "Atribui\u00e7\u00e3o de intera\u00e7\u00f5es de FP",
    "fp_prevalent_interactions_plot": "Intera\u00e7\u00f5es prevalentes de FP",
    "fp_prevalent_interactions_heatmap_plot": "Mapa de calor de intera\u00e7\u00f5es prevalentes de FP",
})
_REPORT_COPY["es"].update({
    "interaction_distribution": "Distribuci\u00f3n de interacciones",
    "interaction_summary": "Resumen de interacciones",
    "amino_acid_distribution": "Interacciones por amino\u00e1cido en todos los ligandos",
    "ligand_atom_distribution": "Interacciones por \u00e1tomos del ligando",
    "interaction_heatmaps": "Mapas de calor de interacciones",
    "heatmaps_by_type": "Mapas de calor por tipo de interacci\u00f3n",
    "complete_interaction_heatmap": "Mapa de calor completo de interacciones",
    "fingerprint_analyses": "An\u00e1lisis de fingerprints",
    "fingerprint_definitions": "EIFP registra entornos de interacci\u00f3n; FIFP registra features de interacci\u00f3n; HIFP registra features h\u00edbridas de interacci\u00f3n. Cada fingerprint describe el patr\u00f3n de reconocimiento molecular con una granularidad estructural diferente.",
    "similarity_matrix": "Matriz de similitud",
    "hierarchical_clustering": "Agrupamiento jer\u00e1rquico",
    "matrix_reordered": "Matriz reordenada por cluster",
    "brief_supervised_learning": "Breve aprendizaje supervisado",
    "supervised_intro": "Los paneles de aprendizaje supervisado clasifican features de fingerprints con Extra Trees y Gradient Boosting. Interprete la importancia junto con la cobertura, la clase asignada, la asignaci\u00f3n de interacciones y la evidencia estructural subyacente.",
    "report_appendix": "Gr\u00e1ficos adicionales disponibles",
    "fp_top50_plot": "Importancia de las top 50 features de FP",
    "fp_class_summary_plot": "Resumen de clases de FP",
    "fp_class_assignment_plot": "Asignaci\u00f3n de clases de FP",
    "fp_coverage_importance_plot": "Cobertura e importancia de FP",
    "fp_feature_presence_heatmap_plot": "Mapa de calor de presencia de features de FP",
    "fp_interaction_assignment_plot": "Asignaci\u00f3n de interacciones de FP",
    "fp_prevalent_interactions_plot": "Interacciones prevalentes de FP",
    "fp_prevalent_interactions_heatmap_plot": "Mapa de calor de interacciones prevalentes de FP",
})


def _report_language(value: str | None = None) -> str:
    code = str(value or active_language() or "en").lower()
    return code if code in _REPORT_LANGUAGES else "en"


def _rt(key: str, language: str | None = None) -> str:
    code = _report_language(language)
    return _REPORT_COPY[code].get(key, _REPORT_COPY["en"].get(key, key))


def _translate_report_data(value: object, language: str) -> str:
    """Translate known data labels while leaving IDs, paths, and values intact."""
    return t(value, lang=_report_language(language))


def _report_display_text(value: object, language: str | None = None) -> str:
    """Translate legacy report literals that may have been stored with bad encoding."""
    text = "" if value is None else str(value)
    code = _report_language(language)
    prefix_map = (
        ("Relat", "report_title"),
        ("Config", "configuration"),
        ("Contagem", "interaction_count"),
        ("Top 30", "top_residues"),
        ("Atribui", "cluster_assignment"),
        ("Ap", "appendix"),
        ("Par", "parameter"),
        ("Valor", "value"),
        ("Gerado", "generated"),
    )
    for prefix, key in prefix_map:
        if text.startswith(prefix):
            return _rt(key, code)
    if "Como interpretar" in text and "fingerprints" in text:
        prefix = text.split(":", 1)[0].strip()
        return f"{prefix}: {_rt('fp_interpret', code)}" if prefix else _rt("fp_interpret", code)
    if "Guia das colunas" in text and "FP" in text:
        prefix = text.split(":", 1)[0].strip()
        return f"{prefix}: {_rt('fp_columns', code)}" if prefix else _rt("fp_columns", code)
    if "Resumo das" in text and "fingerprints" in text:
        prefix = text.split(":", 1)[0].strip()
        return f"{prefix}: {_rt('fp_summary', code)}" if prefix else _rt("fp_summary", code)
    return _translate_report_data(text, code)


_CATEGORY_IMAGE_KIND = {
    "distribution": "distribution",
    "interaction_heatmap": "heatmap",
    "complete_heatmap": "heatmap",
    "similarity": "similarity",
    "clusters": "cluster",
    "fingerprint": "fingerprint",
}


def _localized_report_caption(category: str, caption: object, language: str) -> str:
    image_kind = _CATEGORY_IMAGE_KIND.get(str(category or ""))
    if image_kind:
        return _rt(f"image_{image_kind}_caption", language)
    return _translate_report_data(caption, language)

_FP_COLUMN_GUIDE = [
    ("Feature", "Identificador do bit/atributo do fingerprint usado para localizar a mesma feature em tabelas, gráficos e sessões estruturais."),
    ("Moléculas", "Número de ligantes em que a feature está presente."),
    ("Cobertura (%)", "Percentual da base que contém a feature, calculado sobre todos os ligantes processados."),
    ("Classe prevalente (%)", "Maior participação percentual entre as classes químicas observadas para a feature."),
    ("Z-score classe", "Distância padronizada entre a prevalência da classe e a distribuição das demais features."),
    ("Classe atribuída", "Natureza química aceita para a feature após a aplicação do critério de confiabilidade."),
    ("Importância", "Peso fornecido pelo modelo supervisionado ou pelo fallback analítico para a tarefa configurada."),
    ("Z-score Importance", "Importância padronizada dentro do nível de fingerprint correspondente."),
    ("p-value", "Probabilidade de cauda derivada do Z-score de importância; valores menores indicam maior evidência de relevância."),
    ("Colisões", "Quantidade de ocorrências em que o mesmo bit agregou shells ou naturezas químicas distintas."),
    ("Nível assinado", "Nível do fingerprint atribuído à feature para separar modelos e interpretações por escala estrutural."),
    ("Níveis shell", "Distribuição dos níveis de shell efetivamente associados à feature."),
    ("Níveis colisão", "Níveis de shell encontrados nas ocorrências classificadas como colisão."),
    ("Perfil da base", "Resumo de contagens e percentuais de classes para todas as ocorrências da feature na base."),
]

_FP_EDUCATION = [
    "A seção Análises FP transforma os bits dos fingerprints em variáveis interpretáveis. Ela combina cobertura, natureza química, nível de shell, colisões e importância preditiva para priorizar padrões que merecem inspeção estrutural.",
    "Extra Trees e Gradient Boosting são ajustados para a tarefa ativa, como classificadores para rótulos discretos ou regressores para valores contínuos. Os rankings devem ser comparados: concordância entre os métodos reforça a estabilidade, enquanto divergências sinalizam dependência do modelo.",
    "Importância não prova causalidade molecular. A priorização final deve considerar cobertura, colisões, p-value, interações prevalentes e a posição da feature na estrutura do ligante e do receptor.",
]

_EXPLANATIONS = {
    "similarity": (
        "Matriz de similaridade",
        "Cada célula compara dois ligantes por seus fingerprints de interação. Valores mais altos indicam perfis "
        "mais semelhantes; blocos contíguos sugerem famílias que compartilham modos de reconhecimento molecular.",
    ),
    "cluster": (
        "Agrupamento de ligantes",
        "O dendrograma e a matriz reordenada agrupam ligantes por proximidade entre perfis de interação. Grupos "
        "compactos ajudam a selecionar representantes e a identificar séries com comportamento convergente no sítio.",
    ),
    "heatmap": (
        "Mapa de calor de interacoes",
        "O mapa de calor resume a frequência ou intensidade das interações entre ligantes, resíduos e tipos de "
        "contato. Regiões mais intensas revelam padrões recorrentes e possíveis resíduos-chave para o reconhecimento.",
    ),
    "distribution": (
        "Distribuição de contatos",
        "A distribuição mostra a abundância relativa dos tipos de contato no conjunto analisado. Barras ou regiões "
        "dominantes indicam as forças químicas mais recorrentes e ajudam a comparar densidade e diversidade de contatos.",
    ),
    "network": (
        "Rede de interações",
        "A rede representa entidades moleculares como nós e suas interações como arestas. Conectividade, comunidades "
        "e nós centrais ajudam a reconhecer resíduos ou ligantes que organizam o padrão global de contatos.",
    ),
    "fingerprint": (
        "Análise de fingerprints",
        "O gráfico resume features dos fingerprints de interação, sua frequência, classe ou importância. Features "
        "recorrentes descrevem assinaturas estruturais compartilhadas; features discriminantes ajudam a separar grupos.",
    ),
    "model": (
        "Desempenho do modelo",
        "Esta visualização resume o comportamento do modelo ou das variáveis selecionadas. A interpretação deve "
        "considerar em conjunto separação, erro, cobertura e estabilidade antes de priorizar uma hipótese molecular.",
    ),
    "ligand": (
        "Representação do ligante",
        "A representação relaciona a numeração atômica do ligante aos descritores e contatos usados nas análises. "
        "Ela permite localizar no esqueleto molecular as regiões associadas aos padrões observados.",
    ),
    "generic": (
        "Gráfico de resultados",
        "Esta visualização complementa o conjunto de resultados do projeto. Interprete escalas, rótulos e legenda em "
        "conjunto com as tabelas e com a inspeção estrutural para evitar conclusões baseadas apenas na intensidade visual.",
    ),
}


class PdfReportError(RuntimeError):
    """Raised when the isolated PDF renderer cannot complete safely."""


def _esc(value) -> str:
    return html.escape(str(value))


def _img_b64(path: Path) -> str:
    if not path.exists():
        return ""
    mime = mimetypes.guess_type(path.name)[0] or "image/png"
    return f"data:{mime};base64," + base64.b64encode(path.read_bytes()).decode()


def _image_kind(path: Path) -> str:
    token = " ".join(part.lower() for part in path.parts[-4:])
    if any(value in token for value in ("similarity", "similaridade", "tanimoto")):
        return "similarity"
    if any(value in token for value in ("cluster", "dendrogram")):
        return "cluster"
    if any(value in token for value in ("heatmap", "heat_map", "mapa_calor", "residue_map", "matrix", "matriz")):
        return "heatmap"
    if any(value in token for value in ("network", "rede", "graph", "grafo")):
        return "network"
    if any(value in token for value in ("distribution", "distribu", "interaction_summary", "contact", "stat")):
        return "distribution"
    if any(value in token for value in ("fingerprint", "feature", "importance", "ifp", "fp_")):
        return "fingerprint"
    if any(value in token for value in ("roc", "confusion", "regression", "classification", "prediction", "model")):
        return "model"
    if any(value in token for value in ("ligand_atom", "ligante", "molecule", "molecula", "structure_2d")):
        return "ligand"
    return "generic"


def describe_report_image(path: str | Path, language: str | None = None) -> tuple[str, str]:
    """Return a scientific title and fixed interpretation paragraph for a plot."""
    image_path = Path(path)
    image_kind = _image_kind(image_path)
    base_title = _rt(f"image_{image_kind}_title", language)
    explanation = _rt(f"image_{image_kind}_caption", language)
    label = image_path.stem.replace("_", " ").replace("-", " ").strip()
    label = " ".join(label.split())
    title = f"{base_title}: {label}" if label and label.casefold() not in base_title.casefold() else base_title
    return title, explanation


def collect_result_images(
    workdir: str | Path,
    excluded_paths: list[str | Path] | tuple[str | Path, ...] | None = None,
    language: str | None = None,
) -> list[tuple[str, Path, str]]:
    """Collect every supported chart below the project's ``results`` directory."""
    wd = Path(workdir)
    excluded = {
        str(Path(path).resolve(strict=False)).casefold()
        for path in (excluded_paths or [])
        if path
    }
    candidates: list[Path] = []
    results_dir = wd / "results"
    if results_dir.exists():
        candidates.extend(path for path in results_dir.rglob("*") if path.is_file())

    pages: list[tuple[str, Path, str]] = []
    seen: set[str] = set()
    for path in sorted(candidates, key=lambda item: str(item).casefold()):
        if path.suffix.lower() not in _IMAGE_SUFFIXES:
            continue
        if path.name.lower().startswith(_REPORT_TEMP_PREFIXES):
            continue
        try:
            relative_parts = path.relative_to(results_dir).parts
        except ValueError:
            relative_parts = ()
        if relative_parts and relative_parts[0].casefold() == "plots":
            continue
        resolved = str(path.resolve(strict=False)).casefold()
        if resolved in excluded or resolved in seen:
            continue
        try:
            if path.stat().st_size <= 0:
                continue
        except OSError:
            continue
        seen.add(resolved)
        title, caption = describe_report_image(path, language)
        pages.append((title, path, caption))
    return pages


def _selected_images(
    heatmap_png: str | Path | None,
    interactions_png: str | Path | None,
    cluster_png: str | Path | None,
    extra_images: list[tuple[str, str | Path, str]] | None,
) -> list[tuple[str, Path, str]]:
    pages: list[tuple[str, Path, str]] = []
    if interactions_png:
        pages.append(("Distribuição de interações", Path(interactions_png), _EXPLANATIONS["distribution"][1]))
    if heatmap_png:
        pages.append(("Matriz de similaridade", Path(heatmap_png), _EXPLANATIONS["similarity"][1]))
    if cluster_png:
        pages.append(("Clusters", Path(cluster_png), _EXPLANATIONS["cluster"][1]))
    pages.extend((str(title), Path(image_path), str(caption)) for title, image_path, caption in (extra_images or []))
    return [(title, path, caption) for title, path, caption in pages if path.exists()]


def _selected_images_localized(
    heatmap_png: str | Path | None,
    interactions_png: str | Path | None,
    cluster_png: str | Path | None,
    extra_images: list[tuple[str, str | Path, str]] | None,
    language: str,
) -> list[tuple[str, Path, str]]:
    pages: list[tuple[str, Path, str]] = []
    if interactions_png:
        pages.append(
            (
                _rt("image_distribution_title", language),
                Path(interactions_png),
                _rt("image_distribution_caption", language),
            )
        )
    if heatmap_png:
        pages.append(
            (
                _rt("image_similarity_title", language),
                Path(heatmap_png),
                _rt("image_similarity_caption", language),
            )
        )
    if cluster_png:
        pages.append(
            (
                _rt("image_cluster_title", language),
                Path(cluster_png),
                _rt("image_cluster_caption", language),
            )
        )
    pages.extend(
        (
            _translate_report_data(title, language),
            Path(image_path),
            _translate_report_data(caption, language),
        )
        for title, image_path, caption in (extra_images or [])
    )
    return [(title, path, caption) for title, path, caption in pages if path.exists()]


def _all_report_images(
    cfg: ProjectConfig,
    heatmap_png: str | Path | None,
    interactions_png: str | Path | None,
    cluster_png: str | Path | None,
    extra_images: list[tuple[str, str | Path, str]] | None,
) -> list[tuple[str, Path, str]]:
    semantic = _semantic_report_images(
        cfg,
        heatmap_png,
        interactions_png,
        cluster_png,
        extra_images,
    )
    return [(row["title"], row["path"], row["caption"]) for row in semantic]


def _fallback_sequence(path: Path) -> tuple[int, str]:
    # Imported/legacy results do not have a manifest. Recognize the stable
    # English and Portuguese filenames before relying on broad image keywords.
    stem = path.stem.casefold()
    if "interactions_by_amino" in stem or "interacoes_por_amino" in stem:
        return 11, "distribution"
    if "interactions_by_ligand_atom" in stem or "interacoes_por_atomo" in stem:
        return 12, "distribution"
    if "hierarchical_clustering" in stem:
        return 50, "clusters"
    if "reordered_matrix_cluster" in stem:
        return 60, "clusters"
    if "similarity" in stem:
        return 40, "similarity"
    if "complete" in stem or "completo" in stem:
        return 30, "complete_heatmap"
    if "distribution" in stem or "distribui" in stem or "interaction_summary" in stem:
        return 10, "distribution"
    kind = _image_kind(path)
    if kind == "distribution":
        return 10, "distribution"
    if kind == "heatmap":
        is_complete = "complete" in stem or "completo" in stem
        return (30 if is_complete else 20), ("complete_heatmap" if is_complete else "interaction_heatmap")
    if kind == "similarity":
        return 40, "similarity"
    if kind == "cluster":
        return 50, "clusters"
    if kind == "fingerprint":
        return 100, "fingerprint"
    return 10000, "appendix"


def _path_fp_metadata(path: Path) -> tuple[str, str]:
    tokens = [part.casefold() for part in path.parts]
    ifp_type = next((value for value in _IFP_ORDER if value.casefold() in tokens), "")
    if not ifp_type and "fingerprints" in tokens:
        suffixes = {"e": "EIFP", "f": "FIFP", "h": "HIFP"}
        ifp_type = next((suffixes[token] for token in tokens if token in suffixes), "")
    model = ""
    if "extra_trees" in tokens or "extra trees" in tokens:
        model = "extra_trees"
    elif "gradient_boosting" in tokens or "gradient boosting" in tokens:
        model = "gradient_boosting"
    return ifp_type, model


def _discover_report_profile_images(
    cfg: ProjectConfig,
    language: str,
    excluded_paths: list[Path],
) -> list[dict]:
    """Supplement the manifest with report plots already present on disk.

    Rendering can be interrupted after an image is written but before the
    manifest is saved. A report must still include that valid active-language
    plot instead of silently dropping it.
    """
    root = Path(cfg.workdir) / "results" / "plots" / language / "report"
    if not root.exists():
        return []
    excluded = {str(path.resolve(strict=False)).casefold() for path in excluded_paths}
    rows: list[dict] = []
    for path in sorted(root.rglob("*"), key=lambda item: str(item).casefold()):
        if not path.is_file() or path.suffix.lower() not in _IMAGE_SUFFIXES:
            continue
        key = str(path.resolve(strict=False)).casefold()
        if key in excluded:
            continue
        sequence, category = _fallback_sequence(path)
        ifp_type, model = _path_fp_metadata(path)
        title, caption = describe_report_image(path, language)
        rows.append({
            "plot_id": path.stem,
            "title": title,
            "path": path,
            "caption": caption,
            "sequence": sequence,
            "category": category,
            "ifp_type": ifp_type,
            "model": model,
            "appendix": category == "appendix",
        })
    return rows


def _semantic_report_images(
    cfg: ProjectConfig,
    heatmap_png: str | Path | None,
    interactions_png: str | Path | None,
    cluster_png: str | Path | None,
    extra_images: list[tuple[str, str | Path, str]] | None,
) -> list[dict]:
    language = _report_language(getattr(cfg, "language", "en"))
    manifest = load_manifest(cfg.workdir) if str(cfg.workdir).strip() else None
    manifest_records = manifest.select(language=language, profile="report") if manifest else []
    rows: list[dict] = []
    excluded: list[Path] = []
    if manifest_records:
        for record in manifest_records:
            path = resolve_plot_path(record, cfg.workdir)
            if not path.exists():
                continue
            excluded.append(path)
            rows.append({
                "plot_id": record.plot_id,
                "title": _translate_report_data(record.title, language),
                "path": path,
                "caption": _localized_report_caption(record.category, record.caption, language),
                "sequence": int(record.sequence),
                "category": record.category,
                "ifp_type": record.ifp_type,
                "model": record.model,
                "appendix": record.category == "appendix",
            })
    else:
        selected = _selected_images_localized(
            heatmap_png,
            interactions_png,
            cluster_png,
            extra_images,
            language,
        )
        excluded.extend(item[1] for item in selected)
        for title, path, caption in selected:
            sequence, category = _fallback_sequence(path)
            ifp_type, model = _path_fp_metadata(path)
            rows.append({
                "plot_id": path.stem,
                "title": title,
                "path": path,
                "caption": caption,
                "sequence": sequence,
                "category": category,
                "ifp_type": ifp_type,
                "model": model,
                "appendix": category == "appendix",
            })

    rows.extend(_discover_report_profile_images(cfg, language, excluded))
    discovered = collect_result_images(cfg.workdir, excluded, language) if str(cfg.workdir).strip() else []
    for title, path, caption in discovered:
        sequence, category = _fallback_sequence(path)
        if manifest_records:
            sequence, category = 10000, "appendix"
        ifp_type, model = _path_fp_metadata(path)
        rows.append({
            "plot_id": path.stem,
            "title": title,
            "path": path,
            "caption": caption,
            "sequence": sequence,
            "category": category,
            "ifp_type": ifp_type,
            "model": model,
            "appendix": category == "appendix",
        })
    unique: dict[str, dict] = {}
    for row in rows:
        key = str(row["path"].resolve(strict=False)).casefold()
        unique.setdefault(key, row)
    return sorted(
        unique.values(),
        key=lambda row: (
            bool(row["appendix"]),
            _REPORT_CATEGORY_ORDER.get(row.get("category", ""), int(row["sequence"])),
            int(row["sequence"]),
            _IFP_ORDER.get(row["ifp_type"], 99),
            _MODEL_ORDER.get(row["model"], 99),
            str(row["plot_id"]),
        ),
    )


def _legacy_all_report_images(
    cfg: ProjectConfig,
    heatmap_png: str | Path | None,
    interactions_png: str | Path | None,
    cluster_png: str | Path | None,
    extra_images: list[tuple[str, str | Path, str]] | None,
) -> list[tuple[str, Path, str]]:
    selected = _selected_images(heatmap_png, interactions_png, cluster_png, extra_images)
    discovered = (
        collect_result_images(cfg.workdir, [item[1] for item in selected])
        if str(cfg.workdir).strip()
        else []
    )
    seen: set[str] = set()
    pages: list[tuple[str, Path, str]] = []
    for page in selected + discovered:
        key = str(page[1].resolve(strict=False)).casefold()
        if key not in seen:
            seen.add(key)
            pages.append(page)
    return pages


def _fp_model_tables(fp_dashboards: dict | None) -> list[dict]:
    tables: list[dict] = []
    seen: set[tuple[str, str]] = set()
    for dashboard in (fp_dashboards or {}).values():
        if not isinstance(dashboard, dict):
            continue
        ifp_type = str(dashboard.get("ifp_type") or dashboard.get("ifp_label") or "IFP")
        rankings = dashboard.get("top_features_by_model") or {}
        for model_key, model_title in (
            ("extra_trees", "Extra Trees"),
            ("gradient_boosting", "Gradient Boosting"),
        ):
            identity = (ifp_type, model_key)
            if identity in seen:
                continue
            rows = list(rankings.get(model_key, []) or [])[:50]
            if not rows:
                continue
            seen.add(identity)
            tables.append({
                "ifp_type": ifp_type,
                "model_key": model_key,
                "model_title": model_title,
                "rows": rows,
            })
    return sorted(
        tables,
        key=lambda table: (
            _IFP_ORDER.get(table["ifp_type"], 99),
            _MODEL_ORDER.get(table["model_key"], 99),
        ),
    )


def _build_report_legacy(
    cfg: ProjectConfig,
    analysis: dict,
    heatmap_png: Path | None = None,
    interactions_png: Path | None = None,
    cluster_png: Path | None = None,
    clusters: list[tuple[str, int]] | None = None,
    fp_dashboards: dict | None = None,
    extra_images: list[tuple[str, str | Path, str]] | None = None,
) -> str:
    set_language(str(getattr(cfg, "language", "en") or "en"))
    payload = _pdf_payload(
        cfg,
        analysis,
        heatmap_png,
        interactions_png,
        cluster_png,
        clusters,
        fp_dashboards,
        extra_images,
    )
    rows_cfg = "".join(
        f"<tr><td>{_esc(t(key))}</td><td>{_esc(t(value))}</td></tr>"
        for key, value in payload["cfg_rows"]
    )
    inter_counts = analysis.get("interaction_counts", {}) or {}
    rows_inter = "".join(
        f"<tr><td>{_esc(key)}</td><td class='number'>{_esc(value)}</td></tr>"
        for key, value in sorted(inter_counts.items(), key=lambda item: -item[1])
    ) or "<tr><td colspan='2'>-</td></tr>"
    res_counts = analysis.get("residue_counts", {}) or {}
    rows_res = "".join(
        f"<tr><td>{_esc(key)}</td><td class='number'>{_esc(value)}</td></tr>"
        for key, value in sorted(res_counts.items(), key=lambda item: -item[1])[:30]
    ) or "<tr><td colspan='2'>-</td></tr>"
    cluster_rows = "".join(
        f"<tr><td>{_esc(label)}</td><td class='number'>{cluster_id}</td></tr>"
        for label, cluster_id in payload["clusters"]
    )
    images_html = "".join(
        "<section class='plot'>"
        f"<h2>{_esc(title)}</h2>"
        f"<img src=\"{_img_b64(path)}\" alt=\"{_esc(title)}\">"
        f"<p>{_esc(caption)}</p>"
        f"<small>Fonte: {_esc(path)}</small>"
        "</section>"
        for title, path, caption in [
            (title, Path(path), caption) for title, path, caption in payload["images"]
        ]
    )
    cluster_table = (
        "<h2>Atribuição de clusters</h2>"
        f"<table><tr><th>Ligante</th><th>Cluster</th></tr>{cluster_rows}</table>"
        if cluster_rows else ""
    )
    fp_summary_rows = "".join(
        f"<tr><td>{_esc(key)}</td><td>{_esc(value)}</td></tr>"
        for key, value in payload["fp_rows"]
    )
    fp_summary_section = (
        f"<h2>{_esc(t('Resumo das análises de fingerprints'))}</h2>"
        f"<table><tbody>{fp_summary_rows}</tbody></table>"
        if fp_summary_rows
        else ""
    )
    fp_column_rows = "".join(
        f"<tr><td>{_esc(t(column))}</td><td>{_esc(t(description))}</td></tr>"
        for column, description in _FP_COLUMN_GUIDE
    )
    fp_education = "".join(f"<p>{_esc(t(paragraph))}</p>" for paragraph in _FP_EDUCATION)
    fp_model_html = ""
    for table in payload["fp_model_tables"]:
        model_rows = "".join(
            "<tr>"
            f"<td class='number'>{_esc(row.get('rank', '-'))}</td>"
            f"<td>{_esc(row.get('feature_id', '-'))}</td>"
            f"<td>{_esc(row.get('assigned_level') or '-')}</td>"
            f"<td>{_esc(t(row.get('assigned_class') or '-'))}</td>"
            f"<td class='number'>{float(row.get('coverage_pct', 0.0) or 0.0):.2f}</td>"
            f"<td class='number'>{float(row.get('importance_score', 0.0) or 0.0):.8f}</td>"
            "</tr>"
            for row in table["rows"]
        )
        fp_model_html += (
            "<section class='data-section page-break'>"
            f"<h2>{_esc(t('Top 50 features'))}: {_esc(table['ifp_type'])} / {_esc(table['model_title'])}</h2>"
            "<table><thead><tr>"
            f"<th>{_esc(t('Posição'))}</th><th>{_esc(t('Feature'))}</th>"
            f"<th>{_esc(t('Nível assinado'))}</th><th>{_esc(t('Classe atribuída'))}</th>"
            f"<th>{_esc(t('Cobertura (%)'))}</th><th>{_esc(t('Importância'))}</th>"
            f"</tr></thead><tbody>{model_rows}</tbody></table></section>"
        )

    html_language = {"pt": "pt-br", "en": "en", "es": "es"}.get(str(cfg.language), "en")

    return f"""<!doctype html>
<html lang="{html_language}"><head><meta charset="utf-8">
<title>Relatório HIP²LInterActomics</title>
<style>
@page{{size:A4 landscape;margin:14mm}}
*{{box-sizing:border-box}}
body{{font-family:'Segoe UI',Arial,sans-serif;max-width:1180px;margin:2em auto;padding:0 1.2em;color:#2a221d;background:#faf6f0;line-height:1.45}}
h1{{border-bottom:3px solid #0f766e;padding-bottom:.35em}}
h2{{color:#0f766e;margin-top:1.6em}}
table{{border-collapse:collapse;width:100%;table-layout:fixed;margin:.6em 0;background:#fffdfa}}
thead{{display:table-header-group}} tr{{break-inside:avoid;page-break-inside:avoid}}
th,td{{border:1px solid #ddd3c3;padding:7px 10px;font-size:13px;overflow-wrap:anywhere;vertical-align:top}}
th{{background:#efe7db;text-align:left}} .number{{text-align:right}}
.meta,small{{color:#6a5d52;font-size:12px}}
.plot{{break-before:page;page-break-before:always;break-inside:avoid;page-break-inside:avoid}}
.plot img{{display:block;max-width:100%;max-height:155mm;width:auto;height:auto;object-fit:contain;margin:.6em auto}}
.plot p{{line-height:1.45;background:#fffdfa;border:1px solid #e7d9c5;padding:10px 12px}}
.page-break{{break-before:page;page-break-before:always}} .data-section{{break-inside:auto}}
</style></head><body>
<h1>Relatório HIP²LInterActomics</h1>
<p class="meta">Gerado em {datetime.now().strftime('%Y-%m-%d %H:%M')}</p>
<h2>Resumo</h2>
<ul>
  <li><b>Proteína:</b> {_esc(cfg.protein_file)}</li>
  <li><b>Ligantes:</b> {_esc(cfg.ligand_file)}</li>
  <li><b>Total selecionado:</b> {len(cfg.selected_ligands)}</li>
  <li><b>Workdir:</b> {_esc(cfg.workdir)}</li>
  <li><b>Entradas processadas:</b> {_esc(analysis.get('entries', '-'))}</li>
</ul>
<h2>Configuração</h2>
<table><tr><th>Parâmetro</th><th>Valor</th></tr>{rows_cfg}</table>
<h2>Contagem por tipo de interação</h2>
<table><tr><th>Tipo</th><th>Total</th></tr>{rows_inter}</table>
<h2>Top 30 resíduos com mais interações</h2>
<table><tr><th>Cadeia/Resíduo/Número</th><th>Contagem</th></tr>{rows_res}</table>
<section class="page-break">
<h2>{_esc(t('Como interpretar as análises de fingerprints'))}</h2>
{fp_education}
<h2>{_esc(t('Guia das colunas de Análises FP'))}</h2>
<table><thead><tr><th>{_esc(t('Coluna'))}</th><th>{_esc(t('Interpretação'))}</th></tr></thead><tbody>{fp_column_rows}</tbody></table>
{fp_summary_section}
</section>
{images_html}
{fp_model_html}
{cluster_table}
</body></html>
"""


def _sorted_count_rows(values: dict, limit: int | None = None) -> list[tuple[str, int | float]]:
    rows = sorted(
        ((str(key), value) for key, value in (values or {}).items()),
        key=lambda item: (-float(item[1]), item[0].casefold(), item[0]),
    )
    return rows if limit is None else rows[: max(0, int(limit))]


def _localized_count_summary(
    rows: list[tuple[str, int | float]],
    empty_key: str,
    language: str,
) -> str:
    if not rows:
        return _rt(empty_key, language)
    return ", ".join(
        f"{_translate_report_data(key, language)}: {value}"
        for key, value in rows
    )


def _fp_report_sections(fp_dashboards: dict | None, images: list[dict], language: str) -> list[dict]:
    dashboards_by_type: dict[str, dict] = {}
    for dashboard in (fp_dashboards or {}).values():
        if not isinstance(dashboard, dict):
            continue
        ifp_type = str(dashboard.get("ifp_type") or dashboard.get("ifp_label") or "IFP")
        dashboards_by_type.setdefault(ifp_type, dashboard)
    for image in images:
        if image.get("category") != "fingerprint":
            continue
        ifp_type = str(image.get("ifp_type") or "").upper()
        if ifp_type:
            dashboards_by_type.setdefault(ifp_type, {"ifp_type": ifp_type})
    tables = {
        (table["ifp_type"], table["model_key"]): table
        for table in _fp_model_tables(fp_dashboards)
    }
    sections: list[dict] = []
    for ifp_type, dashboard in sorted(
        dashboards_by_type.items(),
        key=lambda item: (_IFP_ORDER.get(item[0], 99), item[0]),
    ):
        models: list[dict] = []
        for model_key, model_title in (
            ("extra_trees", "Extra Trees"),
            ("gradient_boosting", "Gradient Boosting"),
        ):
            table = tables.get((ifp_type, model_key))
            model_images = [
                image
                for image in images
                if image.get("ifp_type") == ifp_type and image.get("model") == model_key
            ]
            if table or model_images:
                models.append({
                    "key": model_key,
                    "title": model_title,
                    "rows": list(table.get("rows", []) if table else []),
                    "images": model_images,
                })
        sections.append({
            "ifp_type": ifp_type,
            "education": [t(paragraph) for paragraph in _FP_EDUCATION],
            "column_guide": [(t(column), t(description)) for column, description in _FP_COLUMN_GUIDE],
            "summary_rows": _localized_fp_rows({ifp_type: dashboard}, language),
            "models": models,
        })
    return sections


def _load_persisted_fp_dashboards(workdir: str | Path) -> dict[str, dict]:
    """Read terminal-exported FP dashboards without requiring an open GUI tab."""
    root = Path(workdir) / "results" / "terminal" / "fingerprints"
    dashboards: dict[str, dict] = {}
    suffixes = {"E": "EIFP", "F": "FIFP", "H": "HIFP"}
    for suffix, ifp_type in suffixes.items():
        path = root / suffix / "fp_dashboard.json"
        if not path.exists():
            continue
        try:
            data = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, ValueError, TypeError):
            continue
        if isinstance(data, dict):
            data.setdefault("ifp_type", ifp_type)
            dashboards[ifp_type] = data
    return dashboards


def _report_dashboards(workdir: str | Path, fp_dashboards: dict | None) -> dict[str, dict]:
    """Merge live GUI state with on-disk analysis artifacts, preferring live data."""
    merged = _load_persisted_fp_dashboards(workdir)
    for key, dashboard in (fp_dashboards or {}).items():
        if not isinstance(dashboard, dict):
            continue
        ifp_type = str(dashboard.get("ifp_type") or dashboard.get("ifp_label") or key).upper()
        if ifp_type:
            merged[ifp_type] = dashboard
    return merged


def _load_cluster_assignments(workdir: str | Path) -> dict[str, list[tuple[str, str]]]:
    """Load all terminal cluster CSVs so exports do not depend on GUI cache."""
    root = Path(workdir) / "results" / "terminal"
    assignments: dict[str, list[tuple[str, str]]] = {}
    for suffix, ifp_type in (("E", "EIFP"), ("F", "FIFP"), ("H", "HIFP")):
        path = root / f"clusters_{suffix}.csv"
        if not path.exists():
            continue
        try:
            with path.open("r", encoding="utf-8", newline="") as handle:
                rows = [
                    (str(row.get("ligand_id") or ""), str(row.get("cluster_id") or ""))
                    for row in csv.DictReader(handle)
                    if row.get("ligand_id") is not None and row.get("cluster_id") is not None
                ]
        except OSError:
            continue
        if rows:
            assignments[ifp_type] = rows
    return assignments


def _model_plot_order(image: dict) -> tuple[int, int, str]:
    stem = str(image.get("plot_id") or "").casefold()
    order = (
        ("fp_top50", 0),
        ("fp_class_summary", 1),
        ("fp_class_assignment", 2),
        ("fp_coverage_importance", 3),
        ("fp_feature_presence_heatmap", 4),
        ("fp_interaction_assignment", 5),
        ("fp_prevalent_interactions_heatmap", 7),
        ("fp_prevalent_interactions", 6),
    )
    for marker, item_order in order:
        if marker in stem:
            return item_order, int(image.get("sequence", 999)), stem
    return 99, int(image.get("sequence", 999)), stem


def _model_plot_title(image: dict, language: str) -> str:
    stem = str(image.get("plot_id") or "").casefold()
    title_key = next(
        (
            key
            for marker, key in (
                ("fp_top50", "fp_top50_plot"),
                ("fp_class_summary", "fp_class_summary_plot"),
                ("fp_class_assignment", "fp_class_assignment_plot"),
                ("fp_coverage_importance", "fp_coverage_importance_plot"),
                ("fp_feature_presence_heatmap", "fp_feature_presence_heatmap_plot"),
                ("fp_interaction_assignment", "fp_interaction_assignment_plot"),
                ("fp_prevalent_interactions_heatmap", "fp_prevalent_interactions_heatmap_plot"),
                ("fp_prevalent_interactions", "fp_prevalent_interactions_plot"),
            )
            if marker in stem
        ),
        "",
    )
    if not title_key:
        return str(image.get("title") or _rt("image_fingerprint_title", language))
    return f"{image.get('ifp_type') or 'IFP'} / {_MODEL_ORDER_LABELS.get(str(image.get('model') or ''), image.get('model') or '')}: {_rt(title_key, language)}"


_MODEL_ORDER_LABELS = {
    "extra_trees": "Extra Trees",
    "gradient_boosting": "Gradient Boosting",
}


def _numbered_title(number: str, title: str) -> str:
    return f"{number}. {title}" if number else title


def _image_report_section(number: str, title: str, image: dict) -> dict:
    item = dict(image)
    item["title"] = _numbered_title(number, title)
    return {
        "kind": "image",
        "number": number,
        "title": title,
        "image": item,
        "paragraphs": [],
        "table_headers": [],
        "table_rows": [],
        "pdf_rows": [],
    }


def _text_report_section(
    number: str,
    title: str,
    *,
    paragraphs: list[str] | None = None,
    table_headers: list[str] | None = None,
    table_rows: list[tuple] | None = None,
    pdf_rows: list[tuple[str, str]] | None = None,
    kind: str = "data",
) -> dict:
    return {
        "kind": kind,
        "number": number,
        "title": title,
        "paragraphs": list(paragraphs or []),
        "table_headers": list(table_headers or []),
        "table_rows": list(table_rows or []),
        "pdf_rows": list(pdf_rows or []),
    }


def _top_feature_table_rows(rows: list[dict], language: str) -> tuple[list[tuple], list[tuple[str, str]]]:
    table_rows: list[tuple] = []
    pdf_rows: list[tuple[str, str]] = []
    for row in rows[:50]:
        rank = str(row.get("rank", "-"))
        feature_id = str(row.get("feature_id", "-"))
        level = str(row.get("assigned_level") or "-")
        assigned_class = _translate_report_data(row.get("assigned_class") or "-", language)
        coverage = f"{float(row.get('coverage_pct', 0.0) or 0.0):.2f}"
        importance = f"{float(row.get('importance_score', 0.0) or 0.0):.8f}"
        table_rows.append((rank, feature_id, level, assigned_class, coverage, importance))
        pdf_rows.append((
            f"{rank}. {_rt('feature', language)} {feature_id}",
            "; ".join((
                f"{_rt('assigned_level', language)}={level}",
                f"{_rt('assigned_class', language)}={assigned_class}",
                f"{_rt('coverage', language)}={coverage}%",
                f"{_rt('importance', language)}={importance}",
            )),
        ))
    return table_rows, pdf_rows


def _build_report_sections(payload: dict) -> list[dict]:
    """Build the complete report outline shared by HTML and both PDF engines."""
    language = _report_language(payload.get("language"))
    images = [dict(image) for image in payload.get("semantic_images", [])]
    used: set[int] = set()

    def take(predicate) -> list[dict]:
        selected: list[dict] = []
        for index, image in enumerate(images):
            if index in used or not predicate(image):
                continue
            used.add(index)
            selected.append(image)
        return selected

    def display_title(base: str, image: dict) -> str:
        original = str(image.get("title") or "").strip()
        if not original or original.casefold() in base.casefold():
            return base
        return f"{base}: {original}"

    sections: list[dict] = [
        _text_report_section(
            "1",
            _rt("summary", language),
            paragraphs=[f"{_rt('generated', language)} {payload['generated_at']}.", _rt("analysis_intro", language)],
            table_headers=[_rt("parameter", language), _rt("value", language)],
            table_rows=payload["summary_rows"],
            pdf_rows=payload["summary_rows"],
        ),
        _text_report_section(
            "2",
            _rt("configuration", language),
            table_headers=[_rt("parameter", language), _rt("value", language)],
            table_rows=payload["cfg_rows"],
            pdf_rows=payload["cfg_rows"],
        ),
        _text_report_section(
            "3",
            _rt("interaction_count", language),
            table_headers=[_rt("type", language), _rt("total", language)],
            table_rows=payload["interaction_rows"],
            pdf_rows=payload["interaction_rows"],
        ),
        _text_report_section(
            "4",
            _rt("top_residues", language),
            table_headers=[_rt("residue", language), _rt("count", language)],
            table_rows=payload["top_res_rows"],
            pdf_rows=[(str(key), str(value)) for key, value in payload["top_res_rows"]],
        ),
    ]

    distribution_images = take(lambda image: image.get("category") == "distribution")
    if distribution_images:
        sections.append(_text_report_section("5", _rt("interaction_distribution", language), kind="heading"))
        distribution_titles = {
            "interaction_distribution": _rt("interaction_summary", language),
            "interactions_by_amino_acid": _rt("amino_acid_distribution", language),
            "interactions_by_ligand_atom": _rt("ligand_atom_distribution", language),
        }
        distribution_numbers = {
            "interaction_distribution": "5.1",
            "interactions_by_amino_acid": "5.2",
            "interactions_by_ligand_atom": "5.3",
        }
        next_number = 4
        for image in distribution_images:
            plot_id = str(image.get("plot_id") or "")
            number = distribution_numbers.get(plot_id, f"5.{next_number}")
            if plot_id not in distribution_numbers:
                next_number += 1
            sections.append(_image_report_section(
                number,
                display_title(distribution_titles.get(plot_id, str(image["title"])), image),
                image,
            ))

    heatmaps_by_type = take(lambda image: image.get("category") == "interaction_heatmap")
    complete_heatmaps = take(lambda image: image.get("category") == "complete_heatmap")
    if heatmaps_by_type or complete_heatmaps:
        sections.append(_text_report_section("6", _rt("interaction_heatmaps", language), kind="heading"))
        for index, image in enumerate(heatmaps_by_type, start=1):
            sections.append(_image_report_section(
                f"6.1.{index}",
                display_title(_rt("heatmaps_by_type", language), image),
                image,
            ))
        for index, image in enumerate(complete_heatmaps, start=1):
            number = "6.2" if index == 1 else f"6.2.{index}"
            sections.append(_image_report_section(
                number,
                display_title(_rt("complete_interaction_heatmap", language), image),
                image,
            ))

    similarity_images = take(lambda image: image.get("category") == "similarity")
    cluster_images = take(lambda image: image.get("category") == "clusters")
    fingerprint_images = take(lambda image: image.get("category") == "fingerprint" and bool(image.get("ifp_type")))
    dashboards = payload.get("fp_dashboards") or {}
    clusters_by_type = payload.get("clusters_by_type") or {}
    fp_types = {
        str(image.get("ifp_type") or "").upper()
        for image in similarity_images + cluster_images + fingerprint_images
        if image.get("ifp_type")
    }
    fp_types.update(str(key).upper() for key in dashboards if str(key).strip())
    fp_types.update(str(key).upper() for key in clusters_by_type if str(key).strip())
    fp_types.discard("")
    ordered_types = sorted(fp_types, key=lambda value: (_IFP_ORDER.get(value, 99), value))

    if ordered_types:
        sections.append(_text_report_section(
            "7",
            _rt("fingerprint_analyses", language),
            paragraphs=[_rt("fingerprint_definitions", language)],
            kind="heading",
        ))
        hierarchical_images = [image for image in cluster_images if "hierarchical_clustering" in str(image.get("plot_id", ""))]
        reordered_images = [image for image in cluster_images if "reordered_matrix_cluster" in str(image.get("plot_id", ""))]
        other_cluster_images = [image for image in cluster_images if image not in hierarchical_images and image not in reordered_images]
        for kind_number, label, images_for_kind in (
            ("7.1", _rt("similarity_matrix", language), similarity_images),
            ("7.2", _rt("hierarchical_clustering", language), hierarchical_images),
            ("7.3", _rt("matrix_reordered", language), reordered_images + other_cluster_images),
        ):
            for index, image in enumerate(images_for_kind, start=1):
                ifp_type = str(image.get("ifp_type") or "IFP")
                sections.append(_image_report_section(
                    f"{kind_number}.{index}",
                    display_title(f"{ifp_type}: {label}", image),
                    image,
                ))
        for index, ifp_type in enumerate(ordered_types, start=1):
            rows = clusters_by_type.get(ifp_type) or []
            if rows:
                sections.append(_text_report_section(
                    f"7.4.{index}",
                    f"{ifp_type}: {_rt('cluster_assignment', language)}",
                    paragraphs=[_rt("cluster_table_intro", language)],
                    table_headers=[_rt("ligand", language), _rt("cluster", language)],
                    table_rows=rows,
                    pdf_rows=rows,
                ))

    if ordered_types:
        sections.extend((
            _text_report_section(
                "8",
                _rt("brief_supervised_learning", language),
                paragraphs=[_rt("supervised_intro", language)],
                kind="heading",
            ),
            _text_report_section(
                "8.0",
                _rt("fp_interpret", language),
                paragraphs=[_translate_report_data(paragraph, language) for paragraph in _FP_EDUCATION],
                kind="narrative",
            ),
            _text_report_section(
                "8.1",
                _rt("fp_columns", language),
                paragraphs=[_rt("fp_guide_intro", language)],
                table_headers=[_rt("column", language), _rt("interpretation", language)],
                table_rows=[(_translate_report_data(column, language), _translate_report_data(description, language)) for column, description in _FP_COLUMN_GUIDE],
                pdf_rows=[(_translate_report_data(column, language), _translate_report_data(description, language)) for column, description in _FP_COLUMN_GUIDE],
            ),
            _text_report_section(
                "8.2",
                _rt("fp_summary", language),
                paragraphs=[_rt("fp_summary_intro", language)],
                table_headers=[_rt("type", language), _rt("value", language)],
                table_rows=payload.get("fp_rows", []),
                pdf_rows=payload.get("fp_rows", []),
            ),
        ))
        for ifp_index, ifp_type in enumerate(ordered_types, start=1):
            dashboard = dashboards.get(ifp_type) or {"ifp_type": ifp_type}
            model_tables = {
                table["model_key"]: table
                for table in _fp_model_tables({ifp_type: dashboard})
            }
            model_images = [
                image for image in fingerprint_images
                if str(image.get("ifp_type") or "").upper() == ifp_type
            ]
            sections.append(_text_report_section(f"8.2.{ifp_index}", ifp_type, kind="heading"))
            for model_index, (model_key, model_title) in enumerate(_MODEL_ORDER_LABELS.items(), start=1):
                table = model_tables.get(model_key)
                images_for_model = sorted(
                    [image for image in model_images if image.get("model") == model_key],
                    key=_model_plot_order,
                )
                if not table and not images_for_model:
                    continue
                model_number = f"8.2.{ifp_index}.{model_index}"
                if table:
                    table_rows, pdf_rows = _top_feature_table_rows(table["rows"], language)
                    sections.append(_text_report_section(
                        f"{model_number}.1",
                        f"{_rt('top_features', language)}: {ifp_type} / {model_title}",
                        paragraphs=[_rt("top_features_intro", language)],
                        table_headers=[_rt("rank", language), _rt("feature", language), _rt("assigned_level", language), _rt("assigned_class", language), _rt("coverage", language), _rt("importance", language)],
                        table_rows=table_rows,
                        pdf_rows=pdf_rows,
                    ))
                for plot_index, image in enumerate(images_for_model, start=2):
                    sections.append(_image_report_section(
                        f"{model_number}.{plot_index}",
                        display_title(_model_plot_title(image, language), image),
                        image,
                    ))

    appendix_images = [image for index, image in enumerate(images) if index not in used]
    if appendix_images:
        sections.append(_text_report_section("", _rt("report_appendix", language), kind="heading"))
        for index, image in enumerate(appendix_images, start=1):
            sections.append(_image_report_section(f"A.{index}", str(image["title"]), image))
    return sections


def build_report(
    cfg: ProjectConfig,
    analysis: dict,
    heatmap_png: Path | None = None,
    interactions_png: Path | None = None,
    cluster_png: Path | None = None,
    clusters: list[tuple[str, int]] | None = None,
    fp_dashboards: dict | None = None,
    extra_images: list[tuple[str, str | Path, str]] | None = None,
) -> str:
    language = _report_language(getattr(cfg, "language", "en"))
    set_language(language)
    payload = _pdf_payload(
        cfg,
        analysis,
        heatmap_png,
        interactions_png,
        cluster_png,
        clusters,
        fp_dashboards,
        extra_images,
    )
    report_sections = []
    for source in payload["report_sections"]:
        section = dict(source)
        if section.get("kind") == "image":
            image = dict(section["image"])
            image["data_uri"] = _img_b64(Path(image["path"]))
            image["path"] = str(image["path"])
            section["image"] = image
        report_sections.append(section)
    # Kept as empty compatibility fields while third-party templates migrate to
    # the unified outline below.
    general_images: list[dict] = []
    appendix_images: list[dict] = []
    fp_sections: list[dict] = []
    environment = Environment(
        loader=PackageLoader("luna_gui", "templates"),
        autoescape=select_autoescape(("html", "xml")),
    )
    template = environment.get_template("report.html.j2")
    context = {
        "html_language": {"pt": "pt-br", "en": "en", "es": "es"}.get(language, "en"),
        "report_title": t("Relatório HIP²LInterActomics"),
        "generated_label": t("Gerado em"),
        "generated_at": payload["generated_at"],
        "summary_title": t("Resumo"),
        "protein_label": t("Proteína"),
        "ligands_label": t("Ligantes"),
        "selected_label": t("Total selecionado"),
        "processed_label": t("Entradas processadas"),
        "configuration_title": t("Configuração"),
        "parameter_label": t("Parâmetro"),
        "value_label": t("Valor"),
        "interaction_count_title": t("Contagem por tipo de interação"),
        "type_label": t("Tipo"),
        "total_label": t("Total"),
        "top_residue_title": t("Top 30 resíduos com mais interações"),
        "residue_label": t("Cadeia/Resíduo/Número"),
        "count_label": t("Contagem"),
        "source_label": t("Fonte"),
        "cluster_assignment_title": t("Atribuição de clusters"),
        "ligand_label": t("Ligante"),
        "fp_interpret_title": t("Como interpretar as análises de fingerprints"),
        "fp_columns_title": t("Guia das colunas de Análises FP"),
        "fp_summary_title": t("Resumo das análises de fingerprints"),
        "column_label": t("Coluna"),
        "interpretation_label": t("Interpretação"),
        "top_features_title": t("Top 50 features"),
        "rank_label": t("Posição"),
        "assigned_level_label": t("Nível assinado"),
        "assigned_class_label": t("Classe atribuída"),
        "coverage_label": t("Cobertura (%)"),
        "importance_label": t("Importância"),
        "appendix_title": t("Apêndice de exceções"),
        "cfg": cfg,
        "analysis": analysis,
        "cfg_rows": [(t(key), t(value)) for key, value in payload["cfg_rows"]],
        "interaction_rows": [
            (_translate_report_data(key, language), value)
            for key, value in _sorted_count_rows(analysis.get("interaction_counts", {}) or {})
        ],
        "top_residue_rows": _sorted_count_rows(analysis.get("residue_counts", {}) or {}, 30),
        "general_images": general_images,
        "clusters": payload["clusters"],
        "fp_sections": fp_sections,
        "appendix_images": appendix_images,
    }
    context.update({
        "report_title": _rt("report_title", language),
        "generated_label": _rt("generated", language),
        "summary_title": _rt("summary", language),
        "protein_label": _rt("protein", language),
        "ligands_label": _rt("ligands", language),
        "selected_label": _rt("selected_total", language),
        "workdir_label": _rt("workdir", language),
        "processed_label": _rt("processed_entries", language),
        "configuration_title": _rt("configuration", language),
        "parameter_label": _rt("parameter", language),
        "value_label": _rt("value", language),
        "interaction_count_title": _rt("interaction_count", language),
        "type_label": _rt("type", language),
        "total_label": _rt("total", language),
        "top_residue_title": _rt("top_residues", language),
        "residue_label": _rt("residue", language),
        "count_label": _rt("count", language),
        "source_label": _rt("source", language),
        "cluster_assignment_title": _rt("cluster_assignment", language),
        "cluster_label": _rt("cluster", language),
        "ligand_label": _rt("ligand", language),
        "fp_interpret_title": _rt("fp_interpret", language),
        "fp_columns_title": _rt("fp_columns", language),
        "fp_summary_title": _rt("fp_summary", language),
        "column_label": _rt("column", language),
        "interpretation_label": _rt("interpretation", language),
        "top_features_title": _rt("top_features", language),
        "rank_label": _rt("rank", language),
        "feature_label": _rt("feature", language),
        "assigned_level_label": _rt("assigned_level", language),
        "assigned_class_label": _rt("assigned_class", language),
        "coverage_label": _rt("coverage", language),
        "importance_label": _rt("importance", language),
        "appendix_title": _rt("appendix", language),
        "cfg_rows": payload["cfg_rows"],
    })
    context["report_title"] = _rt("report_title", language)
    context["source_label"] = _rt("source", language)
    context["report_sections"] = report_sections
    return template.render(**context)


def save_report(path: str | Path, **kwargs) -> Path:
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_name(f".{output.name}.part")
    temporary.write_text(build_report(**kwargs), encoding="utf-8")
    temporary.replace(output)
    return output


def _cfg_rows_for_pdf(cfg: ProjectConfig) -> list[tuple[str, str]]:
    return [
        ("Proteína", cfg.protein_file or "-"),
        ("Ligantes", cfg.ligand_file or "-"),
        ("Workdir", cfg.workdir or "-"),
        ("Ligantes selecionados", str(len(cfg.selected_ligands))),
        ("Incluir águas", "Sim" if cfg.include_waters else "Não"),
        ("Adicionar H", f"Sim, pH {cfg.ph:g}" if cfg.add_h else "Não"),
        ("IFP", f"{cfg.ifp_type}; níveis={cfg.ifp_levels}; raio={cfg.ifp_radius:g}; tamanho={cfg.ifp_length}"),
        ("Fingerprint de contagem", "Não" if cfg.ifp_bit else "Sim"),
        ("Matriz de similaridade", "Sim" if cfg.sim_matrix else "Não"),
        ("Sessões PyMOL", "Sim" if cfg.out_pse else "Não"),
        ("Filtro binding modes", cfg.binding_modes_cfg or "-"),
        ("Config. interações", cfg.interaction_config_file or "Padrão LUNA"),
        ("Distância máxima global", f"{cfg.inter_max_distance_cap:g} A" if cfg.inter_max_distance_cap else "-"),
        ("Rótulos FP", cfg.fp_labels_csv or "-"),
        ("Tarefa FP", cfg.fp_label_task),
        ("Otsu fallback", "Sim" if getattr(cfg, "fp_use_otsu_threshold", False) else "Nao"),
        ("Núcleos", str(cfg.nproc)),
    ]


def _localized_cfg_rows(cfg: ProjectConfig, language: str) -> list[tuple[str, str]]:
    yes = _rt("yes", language)
    no = _rt("no", language)
    return [
        (_rt("protein", language), cfg.protein_file or "-"),
        (_rt("ligands", language), cfg.ligand_file or "-"),
        (_rt("workdir", language), cfg.workdir or "-"),
        (_rt("cfg_selected_ligands", language), str(len(cfg.selected_ligands))),
        (_rt("cfg_include_waters", language), yes if cfg.include_waters else no),
        (_rt("cfg_add_h", language), f"{yes}, pH {cfg.ph:g}" if cfg.add_h else no),
        (
            "IFP",
            (
                f"{cfg.ifp_type}; {_rt('cfg_levels', language)}={cfg.ifp_levels}; "
                f"{_rt('cfg_radius', language)}={cfg.ifp_radius:g}; "
                f"{_rt('cfg_length', language)}={cfg.ifp_length}"
            ),
        ),
        (_rt("cfg_count_fingerprint", language), no if cfg.ifp_bit else yes),
        (_rt("cfg_similarity", language), yes if cfg.sim_matrix else no),
        (_rt("cfg_pymol", language), yes if cfg.out_pse else no),
        (_rt("cfg_binding_modes", language), cfg.binding_modes_cfg or "-"),
        (_rt("cfg_interactions", language), cfg.interaction_config_file or _rt("default_luna", language)),
        (_rt("cfg_max_distance", language), f"{cfg.inter_max_distance_cap:g} A" if cfg.inter_max_distance_cap else "-"),
        (_rt("cfg_fp_labels", language), cfg.fp_labels_csv or "-"),
        (_rt("cfg_fp_task", language), _translate_report_data(cfg.fp_label_task, language)),
        (_rt("cfg_otsu", language), yes if getattr(cfg, "fp_use_otsu_threshold", False) else no),
        (_rt("cfg_cores", language), str(cfg.nproc)),
    ]


def _save_pdf_page(pdf, fig, page_state: list[int]) -> None:
    page_state[0] += 1
    fig.text(0.965, 0.025, str(page_state[0]), ha="right", va="bottom", fontsize=7.5, color="#6a5d52")
    fig.text(0.045, 0.025, "HIP2LInterActomics", ha="left", va="bottom", fontsize=7.5, color="#6a5d52")
    pdf.savefig(fig, dpi=_PDF_DPI, facecolor=fig.get_facecolor())
    fig.clear()


def _new_text_page(title: str):
    from matplotlib.figure import Figure
    from matplotlib.patches import Rectangle

    fig = Figure(figsize=_A4_LANDSCAPE, dpi=140)
    fig.patch.set_facecolor("#fbf7ef")
    ax = fig.add_subplot(111)
    ax.axis("off")
    ax.add_patch(Rectangle((0.0, 0.915), 1.0, 0.085, transform=ax.transAxes, facecolor="#145c58", edgecolor="none"))
    ax.text(0.045, 0.958, t(title), fontsize=16.5, weight="bold", va="top", color="white")
    return fig, ax, 0.875


def _add_text_page(
    pdf,
    title: str,
    paragraphs: list[str],
    rows: list[tuple[str, str]] | None,
    page_state: list[int],
) -> None:
    from matplotlib.patches import Rectangle

    fig, ax, y = _new_text_page(title)
    content_bottom = 0.075

    def save_and_continue() -> None:
        nonlocal fig, ax, y
        _save_pdf_page(pdf, fig, page_state)
        fig, ax, y = _new_text_page(title)

    def ensure_space(height: float) -> None:
        if y - height < content_bottom:
            save_and_continue()

    for paragraph in paragraphs:
        lines = textwrap.wrap(t(str(paragraph)), width=142, break_long_words=False, break_on_hyphens=False) or [""]
        line_height = 0.027
        ensure_space(line_height * len(lines) + 0.025)
        for line in lines:
            ax.text(0.05, y, line, fontsize=9.2, va="top", color="#2d251e")
            y -= line_height
        y -= 0.021

    if rows:
        def draw_table_header() -> None:
            nonlocal y
            ax.add_patch(Rectangle((0.047, y - 0.027), 0.906, 0.032, transform=ax.transAxes, facecolor="#e7d9c5", edgecolor="#d7c7b2", linewidth=0.35))
            ax.text(0.058, y - 0.004, t("Parâmetro"), fontsize=8.4, weight="bold", va="top", color="#145c58")
            ax.text(0.30, y - 0.004, t("Valor"), fontsize=8.4, weight="bold", va="top", color="#145c58")
            y -= 0.038

        ensure_space(0.05)
        draw_table_header()
        for row_index, (key, value) in enumerate(rows):
            key_lines = textwrap.wrap(t(str(key)), width=31, break_long_words=True, break_on_hyphens=False) or [""]
            value_lines = textwrap.wrap(t(str(value)), width=92, break_long_words=True, break_on_hyphens=False) or [""]
            line_count = max(len(key_lines), len(value_lines))
            row_height = 0.0185 * line_count + 0.010
            if y - row_height < content_bottom:
                save_and_continue()
                draw_table_header()
            ax.add_patch(Rectangle((0.047, y - row_height + 0.004), 0.906, row_height, transform=ax.transAxes, facecolor="#fffdf8" if row_index % 2 == 0 else "#f7efe4", edgecolor="#eadfce", linewidth=0.35))
            for line_index, line in enumerate(key_lines):
                ax.text(0.058, y - 0.005 - line_index * 0.0185, line, fontsize=7.9, weight="bold", va="top", color="#145c58")
            for line_index, line in enumerate(value_lines):
                ax.text(0.30, y - 0.005 - line_index * 0.0185, line, fontsize=7.9, va="top", color="#2d251e")
            y -= row_height + 0.003
    _save_pdf_page(pdf, fig, page_state)


def _fit_image_box(image_width: int, image_height: int) -> tuple[float, float, float, float]:
    """Fit image dimensions into the landscape page without changing aspect ratio."""
    page_width, page_height = _A4_LANDSCAPE
    available_width = 10.65
    available_height = 4.95
    aspect = max(float(image_width), 1.0) / max(float(image_height), 1.0)
    fitted_width = min(available_width, available_height * aspect)
    fitted_height = fitted_width / aspect
    left = (page_width - fitted_width) / 2.0
    bottom = 1.72 + (available_height - fitted_height) / 2.0
    return left / page_width, bottom / page_height, fitted_width / page_width, fitted_height / page_height


def _load_report_image(image_path: Path):
    from PIL import Image, ImageOps

    with Image.open(image_path) as source:
        oriented = ImageOps.exif_transpose(source)
        resampling = getattr(Image, "Resampling", Image)
        oriented.thumbnail((2600, 1800), resampling.LANCZOS)
        converted = oriented.convert("RGBA" if "A" in oriented.getbands() else "RGB")
        return converted.copy()


def _add_image_page(
    pdf,
    title: str,
    image_path: Path,
    caption: str,
    page_state: list[int],
) -> str | None:
    from matplotlib.figure import Figure
    from matplotlib.patches import Rectangle

    if not image_path.exists():
        return f"{_rt('image_missing')}: {image_path}"
    try:
        image = _load_report_image(image_path)
    except Exception as exc:
        _add_text_page(
            pdf,
            f"{_rt('image_unavailable')}: {_report_display_text(title)}",
            [f"{_rt('image_skipped')} {image_path}", f"{type(exc).__name__}: {exc}"],
            None,
            page_state,
        )
        return f"{image_path}: {type(exc).__name__}: {exc}"

    fig = Figure(figsize=_A4_LANDSCAPE, dpi=140)
    fig.patch.set_facecolor("#fbf7ef")
    fig.text(0.045, 0.952, t(str(title)), fontsize=16.0, weight="bold", color="#145c58", va="top")
    left, bottom, width, height = _fit_image_box(image.width, image.height)
    image_ax = fig.add_axes([left, bottom, width, height])
    image_ax.axis("off")
    image_ax.imshow(image, interpolation="antialiased", aspect="equal")

    caption_ax = fig.add_axes([0.045, 0.075, 0.91, 0.12])
    caption_ax.axis("off")
    caption_ax.add_patch(Rectangle((0.0, 0.0), 1.0, 1.0, transform=caption_ax.transAxes, facecolor="#fffdf8", edgecolor="#e7d9c5", linewidth=0.7))
    caption_lines = textwrap.wrap(t(str(caption)), width=168, break_long_words=False, break_on_hyphens=False)
    for line_index, line in enumerate(caption_lines[:5]):
        caption_ax.text(0.015, 0.80 - line_index * 0.19, line, fontsize=8.7, color="#2d251e", va="top")
    fig.text(0.955, 0.205, str(image_path), ha="right", va="bottom", fontsize=6.8, color="#6a5d52")
    _save_pdf_page(pdf, fig, page_state)
    image.close()
    return None


def _fp_rows(fp_dashboards: dict | None) -> list[tuple[str, str]]:
    rows: list[tuple[str, str]] = []
    for key, dashboard in (fp_dashboards or {}).items():
        if not isinstance(dashboard, dict):
            continue
        rows.append(
            (
                str(key),
                f"features={len(dashboard.get('features', []) or [])}; "
                f"importantes={len(dashboard.get('important_features', []) or [])}; "
                f"modelo={dashboard.get('model_name', '-')}; "
                f"limiar={float(dashboard.get('threshold_pct', 0.0) or 0.0):.2f}%",
            )
        )
    return rows[:30]


def _localized_fp_rows(fp_dashboards: dict | None, language: str) -> list[tuple[str, str]]:
    rows: list[tuple[str, str]] = []
    for key, dashboard in (fp_dashboards or {}).items():
        if not isinstance(dashboard, dict):
            continue
        rows.append(
            (
                str(key),
                "; ".join(
                    (
                        f"{_rt('fp_features', language)}={len(dashboard.get('features', []) or [])}",
                        f"{_rt('fp_important', language)}={len(dashboard.get('important_features', []) or [])}",
                        f"{_rt('fp_model', language)}={dashboard.get('model_name', '-')}",
                        f"{_rt('fp_threshold', language)}={float(dashboard.get('threshold_pct', 0.0) or 0.0):.2f}%",
                    )
                ),
            )
        )
    return rows[:30]


def _pdf_payload(
    cfg: ProjectConfig,
    analysis: dict,
    heatmap_png: str | Path | None,
    interactions_png: str | Path | None,
    cluster_png: str | Path | None,
    clusters: list[tuple[str, int]] | None,
    fp_dashboards: dict | None,
    extra_images: list[tuple[str, str | Path, str]] | None,
) -> dict:
    language = _report_language(getattr(cfg, "language", "en"))
    # PDF jobs can be spawned outside Qt, so establish the project locale before
    # building translated FP prose, table labels, and fallback plot descriptions.
    set_language(language)
    analysis = dict(analysis or {})
    inter_counts = analysis.get("interaction_counts", {}) or {}
    residue_counts = analysis.get("residue_counts", {}) or {}
    top_inter = _sorted_count_rows(inter_counts, 12)
    top_res = _sorted_count_rows(residue_counts, 30)
    semantic_images = _semantic_report_images(
        cfg,
        heatmap_png,
        interactions_png,
        cluster_png,
        extra_images,
    )
    dashboards = _report_dashboards(cfg.workdir, fp_dashboards)
    clusters_by_type = _load_cluster_assignments(cfg.workdir)
    if clusters:
        selected_types = list(cfg.selected_ifp_types()) if hasattr(cfg, "selected_ifp_types") else []
        current_type = str(selected_types[0] if selected_types else "IFP").upper()
        clusters_by_type[current_type] = [(str(label), str(cluster_id)) for label, cluster_id in clusters]
    payload = {
        "language": language,
        "analysis": analysis,
        "generated_at": datetime.now().strftime("%Y-%m-%d %H:%M"),
        "entries": analysis.get("entries", "não informado"),
        "summary_rows": [
            ("Proteína", cfg.protein_file or "-"),
            ("Ligantes", cfg.ligand_file or "-"),
            ("Total selecionado", str(len(cfg.selected_ligands))),
            ("Workdir", cfg.workdir or "-"),
            ("Entradas processadas", str(analysis.get("entries", "-"))),
        ],
        "cfg_rows": _cfg_rows_for_pdf(cfg),
        "top_inter": ", ".join(f"{key}: {value}" for key, value in top_inter) or "Sem interações contabilizadas.",
        "top_res": ", ".join(f"{key}: {value}" for key, value in top_res) or "Sem resíduos contabilizados.",
        "images": [
            (row["title"], str(row["path"]), row["caption"])
            for row in semantic_images
        ],
        "semantic_images": semantic_images,
        "fp_dashboards": dashboards,
        "fp_sections": _fp_report_sections(dashboards, semantic_images, language),
        "top_res_rows": top_res,
        "fp_rows": _fp_rows(dashboards),
        "fp_model_tables": _fp_model_tables(dashboards),
        "clusters_by_type": clusters_by_type,
        "clusters": [row for rows in clusters_by_type.values() for row in rows],
    }
    payload.update({
        "entries": analysis.get("entries", _rt("not_reported", language)),
        "summary_rows": [
            (_rt("protein", language), cfg.protein_file or "-"),
            (_rt("ligands", language), cfg.ligand_file or "-"),
            (_rt("selected_total", language), str(len(cfg.selected_ligands))),
            (_rt("workdir", language), cfg.workdir or "-"),
            (_rt("processed_entries", language), str(analysis.get("entries", "-"))),
        ],
        "cfg_rows": _localized_cfg_rows(cfg, language),
        "interaction_rows": [
            (_translate_report_data(key, language), str(value))
            for key, value in _sorted_count_rows(inter_counts)
        ],
        "top_inter": _localized_count_summary(top_inter, "no_interactions", language),
        "top_res": _localized_count_summary(top_res, "no_residues", language),
        "fp_rows": _localized_fp_rows(dashboards, language),
    })
    payload["report_sections"] = _build_report_sections(payload)
    return payload


def _reportlab_write_text_page(canvas, title: str, paragraphs: list[str], rows: list[tuple[str, str]] | None, page_state: list[int]) -> None:
    from reportlab.lib.pagesizes import A4, landscape
    from reportlab.pdfbase.pdfmetrics import stringWidth

    language = _report_language()
    t = lambda value: _report_display_text(value, language)
    page_width, page_height = landscape(A4)
    margin = 36.0
    line_height = 13.0

    def new_page() -> float:
        page_state[0] += 1
        canvas.setFillColorRGB(0.08, 0.36, 0.35)
        canvas.rect(0, page_height - 58, page_width, 58, fill=1, stroke=0)
        canvas.setFillColorRGB(1, 1, 1)
        canvas.setFont("Helvetica-Bold", 15)
        canvas.drawString(margin, page_height - 37, t(title))
        canvas.setFillColorRGB(0.38, 0.33, 0.29)
        canvas.setFont("Helvetica", 7.5)
        canvas.drawString(margin, 20, "HIP2LInterActomics")
        canvas.drawRightString(page_width - margin, 20, str(page_state[0]))
        return page_height - 78

    def wrapped_lines(value: str, width: float, font: str = "Helvetica", size: float = 9.0) -> list[str]:
        words = t(str(value)).split()
        lines: list[str] = []
        current = ""
        for word in words:
            candidate = f"{current} {word}".strip()
            if current and stringWidth(candidate, font, size) > width:
                lines.append(current)
                current = word
            else:
                current = candidate
        if current or not lines:
            lines.append(current)
        return lines

    y = new_page()
    canvas.setFillColorRGB(0.18, 0.15, 0.12)
    for paragraph in paragraphs:
        lines = wrapped_lines(paragraph, page_width - 2 * margin)
        if y - line_height * len(lines) < 42:
            canvas.showPage()
            y = new_page()
        canvas.setFont("Helvetica", 9)
        for line in lines:
            canvas.drawString(margin, y, line)
            y -= line_height
        y -= 8

    if rows:
        key_width = 210.0
        value_width = page_width - 2 * margin - key_width - 16
        for key, value in rows:
            key_lines = wrapped_lines(str(key), key_width, "Helvetica-Bold", 8)
            value_lines = wrapped_lines(str(value), value_width, "Helvetica", 8)
            row_lines = max(len(key_lines), len(value_lines))
            row_height = row_lines * 11 + 8
            if y - row_height < 42:
                canvas.showPage()
                y = new_page()
            canvas.setFillColorRGB(0.98, 0.96, 0.92)
            canvas.rect(margin, y - row_height + 3, page_width - 2 * margin, row_height, fill=1, stroke=0)
            canvas.setFillColorRGB(0.08, 0.36, 0.35)
            canvas.setFont("Helvetica-Bold", 8)
            for index, line in enumerate(key_lines):
                canvas.drawString(margin + 6, y - 8 - index * 11, line)
            canvas.setFillColorRGB(0.18, 0.15, 0.12)
            canvas.setFont("Helvetica", 8)
            for index, line in enumerate(value_lines):
                canvas.drawString(margin + key_width + 12, y - 8 - index * 11, line)
            y -= row_height + 3
    canvas.showPage()


def _reportlab_write_image_page(canvas, image: dict, page_state: list[int]) -> str | None:
    from reportlab.lib.pagesizes import A4, landscape
    from reportlab.lib.utils import ImageReader
    from reportlab.pdfbase.pdfmetrics import stringWidth
    from PIL import Image

    language = _report_language()
    path = Path(image["path"])
    if not path.exists():
        return f"Imagem não encontrada: {path}"
    page_width, page_height = landscape(A4)
    margin = 36.0
    try:
        with Image.open(path) as source:
            image_width, image_height = source.size
    except Exception as exc:
        return f"{path}: {type(exc).__name__}: {exc}"

    page_state[0] += 1
    canvas.setFillColorRGB(0.08, 0.36, 0.35)
    canvas.setFont("Helvetica-Bold", 14)
    canvas.drawString(margin, page_height - 38, _report_display_text(image["title"], language))
    available_width = page_width - 2 * margin
    available_height = page_height - 150
    scale = min(available_width / max(image_width, 1), available_height / max(image_height, 1))
    draw_width = image_width * scale
    draw_height = image_height * scale
    x = (page_width - draw_width) / 2
    y = 78 + (available_height - draw_height) / 2
    reader = ImageReader(str(path))
    canvas.drawImage(reader, x, y, width=draw_width, height=draw_height, preserveAspectRatio=True, anchor="c", mask="auto")

    caption = _report_display_text(image.get("caption") or "", language)
    words = caption.split()
    lines: list[str] = []
    current = ""
    for word in words:
        candidate = f"{current} {word}".strip()
        if current and stringWidth(candidate, "Helvetica", 8) > available_width:
            lines.append(current)
            current = word
        else:
            current = candidate
    if current:
        lines.append(current)
    canvas.setFillColorRGB(0.18, 0.15, 0.12)
    canvas.setFont("Helvetica", 8)
    for index, line in enumerate(lines[:4]):
        canvas.drawString(margin, 62 - index * 10, line)
    canvas.setFillColorRGB(0.38, 0.33, 0.29)
    canvas.setFont("Helvetica", 7)
    canvas.drawString(margin, 20, "HIP2LInterActomics")
    canvas.drawRightString(page_width - margin, 20, str(page_state[0]))
    canvas.showPage()
    del reader
    gc.collect()
    return None


def _write_reportlab_sections(canvas, payload: dict, page_state: list[int]) -> list[str]:
    """Render the ordered outline used by the HTML report, one section at a time."""
    warnings: list[str] = []
    for section in payload.get("report_sections", []):
        if section.get("kind") == "image":
            warning = _reportlab_write_image_page(canvas, section["image"], page_state)
            if warning:
                warnings.append(warning)
            continue
        title = _numbered_title(str(section.get("number") or ""), str(section.get("title") or ""))
        _reportlab_write_text_page(
            canvas,
            title,
            list(section.get("paragraphs") or []),
            list(section.get("pdf_rows") or []),
            page_state,
        )
    return warnings


def _write_matplotlib_sections(pdf, payload: dict, page_state: list[int]) -> list[str]:
    """Fallback equivalent of the ReportLab section renderer."""
    warnings: list[str] = []
    language = _report_language(payload.get("language"))
    for section in payload.get("report_sections", []):
        if section.get("kind") == "image":
            image = section["image"]
            warning = _add_image_page(
                pdf,
                str(image["title"]),
                Path(image["path"]),
                str(image.get("caption") or ""),
                page_state,
            )
            if warning:
                warnings.append(warning)
            continue
        _add_text_page(
            pdf,
            _numbered_title(str(section.get("number") or ""), str(section.get("title") or "")),
            list(section.get("paragraphs") or []),
            list(section.get("pdf_rows") or []),
            page_state,
        )
    return warnings


def _write_pdf_payload_reportlab_legacy(path: str | Path, payload: dict) -> tuple[Path, list[str]]:
    from reportlab.lib.pagesizes import A4, landscape
    from reportlab.pdfgen.canvas import Canvas

    set_language(str(payload.get("language") or "en"))
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_name(f".{output.name}.part")
    temporary.unlink(missing_ok=True)
    warnings: list[str] = []
    page_state = [0]
    canvas = Canvas(str(temporary), pagesize=landscape(A4), pageCompression=1)
    try:
        _reportlab_write_text_page(
            canvas,
            "Relatório HIP²LInterActomics",
            [
                f"{t('Gerado em')} {payload['generated_at']}.",
            ],
            payload["summary_rows"],
            page_state,
        )
        _reportlab_write_text_page(
            canvas,
            "Configuração",
            [],
            payload["cfg_rows"],
            page_state,
        )
        _reportlab_write_text_page(
            canvas,
            "Contagem por tipo de interação",
            [],
            [
                (str(key), str(value))
                for key, value in _sorted_count_rows(
                    (payload.get("analysis") or {}).get("interaction_counts", {})
                )
            ],
            page_state,
        )
        _reportlab_write_text_page(
            canvas,
            "Top 30 resíduos com mais interações",
            [],
            [(str(key), str(value)) for key, value in payload["top_res_rows"]],
            page_state,
        )
        general_images = [
            image for image in payload["semantic_images"]
            if not image.get("appendix") and image.get("category") != "fingerprint"
        ]
        for image in general_images:
            warning = _reportlab_write_image_page(canvas, image, page_state)
            if warning:
                warnings.append(warning)
        if payload["clusters"]:
            _reportlab_write_text_page(
                canvas,
                "Atribuição de clusters",
                ["Tabela dos ligantes e seus grupos hierárquicos."],
                payload["clusters"],
                page_state,
            )
        for section in payload["fp_sections"]:
            _reportlab_write_text_page(
                canvas,
                f"{section['ifp_type']}: Como interpretar as análises de fingerprints",
                section["education"],
                None,
                page_state,
            )
            _reportlab_write_text_page(
                canvas,
                f"{section['ifp_type']}: Guia das colunas de Análises FP",
                [],
                section["column_guide"],
                page_state,
            )
            _reportlab_write_text_page(
                canvas,
                f"{section['ifp_type']}: Resumo das análises de fingerprints",
                [],
                section["summary_rows"],
                page_state,
            )
            for model in section["models"]:
                model_rows = [
                    (
                        f"{row.get('rank', '-')}. feature {row.get('feature_id', '-')}",
                        f"nível={row.get('assigned_level') or '-'}; classe={row.get('assigned_class') or '-'}; cobertura={float(row.get('coverage_pct', 0.0) or 0.0):.2f}%; importância={float(row.get('importance_score', 0.0) or 0.0):.8f}",
                    )
                    for row in model["rows"]
                ]
                _reportlab_write_text_page(
                    canvas,
                    f"Top 50 features: {section['ifp_type']} / {model['title']}",
                    [],
                    model_rows,
                    page_state,
                )
                for image in model["images"]:
                    warning = _reportlab_write_image_page(canvas, image, page_state)
                    if warning:
                        warnings.append(warning)
        appendix = [
            image
            for image in payload["semantic_images"]
            if image.get("appendix") or (
                image.get("category") == "fingerprint"
                and (not image.get("ifp_type") or not image.get("model"))
            )
        ]
        if appendix:
            _reportlab_write_text_page(canvas, "Apêndice de exceções", [], None, page_state)
            for image in appendix:
                warning = _reportlab_write_image_page(canvas, image, page_state)
                if warning:
                    warnings.append(warning)
        canvas.save()
        temporary.replace(output)
    except BaseException:
        try:
            canvas.save()
        except Exception:
            pass
        temporary.unlink(missing_ok=True)
        raise
    return output, warnings


def _write_pdf_payload_reportlab(path: str | Path, payload: dict) -> tuple[Path, list[str]]:
    """Render the PDF from the exact report outline used by the HTML template."""
    from reportlab.lib.pagesizes import A4, landscape
    from reportlab.pdfgen.canvas import Canvas

    language = _report_language(payload.get("language"))
    set_language(language)
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_name(f".{output.name}.part")
    temporary.unlink(missing_ok=True)
    page_state = [0]
    canvas = Canvas(str(temporary), pagesize=landscape(A4), pageCompression=1)
    try:
        _reportlab_write_text_page(
            canvas,
            _rt("report_title", language),
            [f"{_rt('generated', language)} {payload['generated_at']}."],
            None,
            page_state,
        )
        warnings = _write_reportlab_sections(canvas, payload, page_state)
        canvas.save()
        temporary.replace(output)
    except BaseException:
        try:
            canvas.save()
        except Exception:
            pass
        temporary.unlink(missing_ok=True)
        raise
    return output, warnings


def _write_pdf_payload_matplotlib_legacy(path: str | Path, payload: dict) -> tuple[Path, list[str]]:
    from matplotlib.backends.backend_pdf import PdfPages

    set_language(str(payload.get("language") or "en"))
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_name(f".{output.name}.part")
    temporary.unlink(missing_ok=True)
    page_state = [0]
    warnings: list[str] = []
    try:
        with PdfPages(temporary) as pdf:
            _add_text_page(
                pdf,
                "Relatório HIP²LInterActomics",
                [
                    f"Gerado em {payload['generated_at']}.",
                    "Este relatório reúne os parâmetros do projeto, o resumo numérico e todos os gráficos encontrados na pasta de resultados. Padrões visuais devem ser tratados como hipóteses para validação estrutural e experimental.",
                    f"Entradas processadas: {payload['entries']}.",
                ],
                payload["cfg_rows"],
                page_state,
            )
            _add_text_page(
                pdf,
                "Como interpretar as análises",
                [
                    _EXPLANATIONS["distribution"][1],
                    _EXPLANATIONS["heatmap"][1],
                    _EXPLANATIONS["similarity"][1],
                    _EXPLANATIONS["cluster"][1],
                    _EXPLANATIONS["fingerprint"][1],
                    _EXPLANATIONS["network"][1],
                ],
                [("Interações mais frequentes", payload["top_inter"]), ("Resíduos mais frequentes", payload["top_res"])],
                page_state,
            )
            _add_text_page(
                pdf,
                "Como interpretar as análises de fingerprints",
                list(_FP_EDUCATION),
                None,
                page_state,
            )
            _add_text_page(
                pdf,
                "Guia das colunas de Análises FP",
                ["Consulte este dicionário ao interpretar as tabelas, os rankings e os mapas de fingerprints."],
                list(_FP_COLUMN_GUIDE),
                page_state,
            )
            if payload["fp_rows"]:
                _add_text_page(
                    pdf,
                    "Resumo das análises de fingerprints",
                    ["Cada linha resume uma base de fingerprints, a quantidade de features e o corte usado na seleção."],
                    payload["fp_rows"],
                    page_state,
                )
            for title, image_path, caption in payload["images"]:
                warning = _add_image_page(pdf, title, Path(image_path), caption, page_state)
                if warning:
                    warnings.append(warning)
            for table in payload["fp_model_tables"]:
                rows = [
                    (
                        f"{row.get('rank', '-')}. feature {row.get('feature_id', '-')} (L{row.get('assigned_level') or '-'})",
                        f"classe={row.get('assigned_class') or '-'}; cobertura={float(row.get('coverage_pct', 0.0) or 0.0):.2f}%; importância={float(row.get('importance_score', 0.0) or 0.0):.8f}",
                    )
                    for row in table["rows"]
                ]
                _add_text_page(
                    pdf,
                    f"Top 50 features: {table['ifp_type']} / {table['model_title']}",
                    ["Ranking independente das features para comparação entre os dois métodos de ensemble."],
                    rows,
                    page_state,
                )
            if payload["clusters"]:
                _add_text_page(pdf, "Atribuição de clusters", ["Tabela dos ligantes e seus grupos hierárquicos."], payload["clusters"], page_state)
        temporary.replace(output)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise
    return output, warnings


def _write_pdf_payload_matplotlib(path: str | Path, payload: dict) -> tuple[Path, list[str]]:
    """Fallback renderer that follows the same sections and order as HTML."""
    from matplotlib.backends.backend_pdf import PdfPages

    language = _report_language(payload.get("language"))
    set_language(language)
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_name(f".{output.name}.part")
    temporary.unlink(missing_ok=True)
    page_state = [0]
    try:
        with PdfPages(temporary) as pdf:
            _add_text_page(
                pdf,
                _rt("report_title", language),
                [f"{_rt('generated', language)} {payload['generated_at']}."],
                None,
                page_state,
            )
            warnings = _write_matplotlib_sections(pdf, payload, page_state)
        temporary.replace(output)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise
    return output, warnings


def _write_pdf_payload(path: str | Path, payload: dict) -> tuple[Path, list[str]]:
    """Prefer incremental ReportLab pages, retaining the legacy fallback."""
    try:
        import reportlab  # noqa: F401
    except ImportError:
        return _write_pdf_payload_matplotlib(path, payload)
    return _write_pdf_payload_reportlab(path, payload)


def save_pdf_report(
    path: str | Path,
    cfg: ProjectConfig,
    analysis: dict,
    heatmap_png: Path | None = None,
    interactions_png: Path | None = None,
    cluster_png: Path | None = None,
    clusters: list[tuple[str, int]] | None = None,
    fp_dashboards: dict | None = None,
    extra_images: list[tuple[str, str | Path, str]] | None = None,
) -> Path:
    """Generate the PDF in-process. Prefer ``save_pdf_report_isolated`` in GUIs."""
    payload = _pdf_payload(cfg, analysis, heatmap_png, interactions_png, cluster_png, clusters, fp_dashboards, extra_images)
    output, _warnings = _write_pdf_payload(path, payload)
    return output


def _write_worker_status(status_path: str, payload: dict) -> None:
    output = Path(status_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_name(f".{output.name}.part")
    temporary.write_text(json.dumps(payload, ensure_ascii=False), encoding="utf-8")
    temporary.replace(output)


def create_pdf_render_job(
    path: str | Path,
    cfg: ProjectConfig,
    analysis: dict,
    heatmap_png: Path | None = None,
    interactions_png: Path | None = None,
    cluster_png: Path | None = None,
    clusters: list[tuple[str, int]] | None = None,
    fp_dashboards: dict | None = None,
    extra_images: list[tuple[str, str | Path, str]] | None = None,
) -> tuple[Path, Path]:
    """Serialize a trusted local PDF job for a dedicated renderer process."""
    payload = _pdf_payload(
        cfg,
        analysis,
        heatmap_png,
        interactions_png,
        cluster_png,
        clusters,
        fp_dashboards,
        extra_images,
    )
    job_dir = Path(tempfile.mkdtemp(prefix="hip2l-pdf-job-"))
    job_path = job_dir / "job.pickle"
    status_path = job_dir / "status.json"
    temporary = job_dir / ".job.pickle.part"
    try:
        with temporary.open("wb") as handle:
            pickle.dump(
                {"output": str(Path(path)), "payload": payload},
                handle,
                protocol=pickle.HIGHEST_PROTOCOL,
            )
        temporary.replace(job_path)
    except BaseException:
        shutil.rmtree(job_dir, ignore_errors=True)
        raise
    return job_path, status_path


def execute_pdf_render_job(job_path: str | Path, status_path: str | Path) -> int:
    """Execute a serialized job without importing Qt in the worker process."""
    try:
        with Path(job_path).open("rb") as handle:
            job = pickle.load(handle)
        output, warnings = _write_pdf_payload(job["output"], job["payload"])
        _write_worker_status(
            str(status_path),
            {"ok": True, "path": str(output), "warnings": warnings},
        )
        return 0
    except BaseException as exc:
        try:
            _write_worker_status(
                str(status_path),
                {
                    "ok": False,
                    "error": f"{type(exc).__name__}: {exc}",
                    "traceback": traceback.format_exc(limit=12),
                },
            )
        except BaseException:
            pass
        return 1


def read_pdf_render_status(status_path: str | Path) -> dict | None:
    try:
        payload = json.loads(Path(status_path).read_text(encoding="utf-8"))
    except (OSError, ValueError, TypeError):
        return None
    return payload if isinstance(payload, dict) else None


def cleanup_pdf_render_job(job_path: str | Path) -> None:
    shutil.rmtree(Path(job_path).parent, ignore_errors=True)


def _pdf_process_worker(status_path: str, path: str, payload: dict) -> None:
    try:
        output, warnings = _write_pdf_payload(path, payload)
        _write_worker_status(status_path, {"ok": True, "path": str(output), "warnings": warnings})
    except BaseException as exc:
        _write_worker_status(status_path, {
            "ok": False,
            "error": f"{type(exc).__name__}: {exc}",
            "traceback": traceback.format_exc(limit=12),
        })


def save_pdf_report_isolated(
    path: str | Path,
    cfg: ProjectConfig,
    analysis: dict,
    heatmap_png: Path | None = None,
    interactions_png: Path | None = None,
    cluster_png: Path | None = None,
    clusters: list[tuple[str, int]] | None = None,
    fp_dashboards: dict | None = None,
    extra_images: list[tuple[str, str | Path, str]] | None = None,
    timeout: int = 900,
    progress_callback: Callable[[], None] | None = None,
) -> Path:
    """Render in a spawned process so native PDF failures cannot terminate Qt."""
    payload = _pdf_payload(cfg, analysis, heatmap_png, interactions_png, cluster_png, clusters, fp_dashboards, extra_images)
    context = multiprocessing.get_context("spawn")
    status_handle, status_path = tempfile.mkstemp(prefix="hip2l-pdf-", suffix=".json")
    os.close(status_handle)
    Path(status_path).unlink(missing_ok=True)
    process = context.Process(target=_pdf_process_worker, args=(status_path, str(path), payload), name="hip2l-pdf-renderer")
    try:
        process.start()
    except Exception as exc:
        Path(status_path).unlink(missing_ok=True)
        raise PdfReportError(f"Nao foi possivel iniciar o gerador de PDF: {exc}") from exc

    deadline = time.monotonic() + max(int(timeout), 1)
    while process.is_alive() and time.monotonic() < deadline:
        process.join(0.1)
        if progress_callback is not None:
            try:
                progress_callback()
            except Exception:
                pass

    if process.is_alive():
        process.terminate()
        process.join(5)
        if process.is_alive() and hasattr(process, "kill"):
            process.kill()
            process.join(2)
        Path(status_path).unlink(missing_ok=True)
        process.close()
        raise PdfReportError(f"A geracao do PDF excedeu {timeout} segundos e foi interrompida sem fechar o aplicativo.")

    message = None
    status_file = Path(status_path)
    if status_file.exists():
        try:
            message = json.loads(status_file.read_text(encoding="utf-8"))
        except (OSError, ValueError):
            message = None
    status_file.unlink(missing_ok=True)
    exit_code = process.exitcode
    process.close()

    if not message:
        raise PdfReportError(
            "O processo isolado de PDF foi encerrado inesperadamente "
            f"(codigo {exit_code}). Isso pode indicar falta de memoria ou falha de uma biblioteca nativa; o aplicativo permaneceu aberto."
        )
    if not message.get("ok"):
        detail = str(message.get("error") or "erro desconhecido")
        raise PdfReportError(f"Falha ao gerar o PDF: {detail}")
    output = Path(message.get("path") or path)
    if not output.exists():
        raise PdfReportError("O gerador terminou sem produzir o arquivo PDF.")
    return output
