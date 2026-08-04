"""HTML and PDF reports for completed HIP2LInterActomics projects."""
from __future__ import annotations

import base64
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
from ..i18n import set_language, t


_A4_LANDSCAPE = (11.69, 8.27)
_PDF_DPI = 300
_IMAGE_SUFFIXES = {".png", ".jpg", ".jpeg", ".tif", ".tiff", ".bmp"}
_REPORT_TEMP_PREFIXES = ("_report_", "_report_pdf_")
_IFP_ORDER = {"EIFP": 0, "FIFP": 1, "HIFP": 2}
_MODEL_ORDER = {"extra_trees": 0, "gradient_boosting": 1}

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


def describe_report_image(path: str | Path) -> tuple[str, str]:
    """Return a scientific title and fixed interpretation paragraph for a plot."""
    image_path = Path(path)
    base_title, explanation = _EXPLANATIONS[_image_kind(image_path)]
    label = image_path.stem.replace("_", " ").replace("-", " ").strip()
    label = " ".join(label.split())
    title = f"{base_title}: {label}" if label and label.casefold() not in base_title.casefold() else base_title
    return title, explanation


def collect_result_images(
    workdir: str | Path,
    excluded_paths: list[str | Path] | tuple[str | Path, ...] | None = None,
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
        title, caption = describe_report_image(path)
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
    kind = _image_kind(path)
    if kind == "distribution":
        return 10, "distribution"
    if kind == "heatmap":
        is_complete = "complete" in path.stem.casefold() or "completo" in path.stem.casefold()
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
    model = ""
    if "extra_trees" in tokens or "extra trees" in tokens:
        model = "extra_trees"
    elif "gradient_boosting" in tokens or "gradient boosting" in tokens:
        model = "gradient_boosting"
    return ifp_type, model


def _semantic_report_images(
    cfg: ProjectConfig,
    heatmap_png: str | Path | None,
    interactions_png: str | Path | None,
    cluster_png: str | Path | None,
    extra_images: list[tuple[str, str | Path, str]] | None,
) -> list[dict]:
    language = str(getattr(cfg, "language", "en") or "en")
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
                "title": record.title,
                "path": path,
                "caption": record.caption,
                "sequence": int(record.sequence),
                "category": record.category,
                "ifp_type": record.ifp_type,
                "model": record.model,
                "appendix": record.category == "appendix",
            })
    else:
        selected = _selected_images(heatmap_png, interactions_png, cluster_png, extra_images)
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

    discovered = collect_result_images(cfg.workdir, excluded) if str(cfg.workdir).strip() else []
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


def _fp_report_sections(fp_dashboards: dict | None, images: list[dict]) -> list[dict]:
    dashboards_by_type: dict[str, dict] = {}
    for dashboard in (fp_dashboards or {}).values():
        if not isinstance(dashboard, dict):
            continue
        ifp_type = str(dashboard.get("ifp_type") or dashboard.get("ifp_label") or "IFP")
        dashboards_by_type.setdefault(ifp_type, dashboard)
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
            "summary_rows": _fp_rows({ifp_type: dashboard}),
            "models": models,
        })
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
    language = str(getattr(cfg, "language", "en") or "en")
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
    images = [dict(row) for row in payload["semantic_images"]]
    for row in images:
        row["data_uri"] = _img_b64(Path(row["path"]))
        row["path"] = str(row["path"])
    general_images = [
        row for row in images
        if not row.get("appendix") and row.get("category") != "fingerprint"
    ]
    appendix_images = [
        row for row in images
        if row.get("appendix") or (
            row.get("category") == "fingerprint"
            and (not row.get("ifp_type") or not row.get("model"))
        )
    ]
    fp_sections = _fp_report_sections(fp_dashboards, images)
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
        "interaction_rows": _sorted_count_rows(analysis.get("interaction_counts", {}) or {}),
        "top_residue_rows": _sorted_count_rows(analysis.get("residue_counts", {}) or {}, 30),
        "general_images": general_images,
        "clusters": payload["clusters"],
        "fp_sections": fp_sections,
        "appendix_images": appendix_images,
    }
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
        return f"Imagem não encontrada: {image_path}"
    try:
        image = _load_report_image(image_path)
    except Exception as exc:
        _add_text_page(
            pdf,
            f"Gráfico indisponível: {title}",
            [f"O arquivo {image_path} não pode ser lido e foi ignorado.", f"{type(exc).__name__}: {exc}"],
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
    return {
        "language": str(getattr(cfg, "language", "en") or "en"),
        "generated_at": datetime.now().strftime("%Y-%m-%d %H:%M"),
        "entries": analysis.get("entries", "não informado"),
        "cfg_rows": _cfg_rows_for_pdf(cfg),
        "top_inter": ", ".join(f"{key}: {value}" for key, value in top_inter) or "Sem interações contabilizadas.",
        "top_res": ", ".join(f"{key}: {value}" for key, value in top_res) or "Sem resíduos contabilizados.",
        "images": [
            (row["title"], str(row["path"]), row["caption"])
            for row in semantic_images
        ],
        "semantic_images": semantic_images,
        "fp_sections": _fp_report_sections(fp_dashboards, semantic_images),
        "top_res_rows": top_res,
        "fp_rows": _fp_rows(fp_dashboards),
        "fp_model_tables": _fp_model_tables(fp_dashboards),
        "clusters": [(str(label), str(cluster_id)) for label, cluster_id in (clusters or [])],
    }


def _reportlab_write_text_page(canvas, title: str, paragraphs: list[str], rows: list[tuple[str, str]] | None, page_state: list[int]) -> None:
    from reportlab.lib.pagesizes import A4, landscape
    from reportlab.pdfbase.pdfmetrics import stringWidth

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
    canvas.drawString(margin, page_height - 38, t(str(image["title"])))
    available_width = page_width - 2 * margin
    available_height = page_height - 150
    scale = min(available_width / max(image_width, 1), available_height / max(image_height, 1))
    draw_width = image_width * scale
    draw_height = image_height * scale
    x = (page_width - draw_width) / 2
    y = 78 + (available_height - draw_height) / 2
    reader = ImageReader(str(path))
    canvas.drawImage(reader, x, y, width=draw_width, height=draw_height, preserveAspectRatio=True, anchor="c", mask="auto")

    caption = t(str(image.get("caption") or ""))
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


def _write_pdf_payload_reportlab(path: str | Path, payload: dict) -> tuple[Path, list[str]]:
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
                f"Gerado em {payload['generated_at']}.",
                f"Entradas processadas: {payload['entries']}.",
            ],
            payload["cfg_rows"],
            page_state,
        )
        _reportlab_write_text_page(
            canvas,
            "Top 30 resíduos com mais interações",
            [payload["top_inter"]],
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
        appendix = [image for image in payload["semantic_images"] if image.get("appendix")]
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


def _write_pdf_payload_matplotlib(path: str | Path, payload: dict) -> tuple[Path, list[str]]:
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
