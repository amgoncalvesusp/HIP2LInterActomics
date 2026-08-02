"""HTML and PDF reports for completed HIP2LInterActomics projects."""
from __future__ import annotations

import base64
import html
import json
import mimetypes
import multiprocessing
import os
import tempfile
import textwrap
import time
import traceback
from datetime import datetime
from pathlib import Path
from typing import Callable

from .project import ProjectConfig
from ..i18n import t


_A4_LANDSCAPE = (11.69, 8.27)
_PDF_DPI = 180
_IMAGE_SUFFIXES = {".png", ".jpg", ".jpeg", ".tif", ".tiff", ".bmp"}
_REPORT_TEMP_PREFIXES = ("_report_", "_report_pdf_")

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


def build_report(
    cfg: ProjectConfig,
    analysis: dict,
    heatmap_png: Path | None = None,
    interactions_png: Path | None = None,
    cluster_png: Path | None = None,
    clusters: list[tuple[str, int]] | None = None,
    extra_images: list[tuple[str, str | Path, str]] | None = None,
) -> str:
    rows_cfg = "".join(
        f"<tr><td>{_esc(key)}</td><td>{_esc(value)}</td></tr>"
        for key, value in vars(cfg).items()
        if not (isinstance(value, list) and not value) and value not in ("", False, None)
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
        for label, cluster_id in (clusters or [])
    )
    images_html = "".join(
        "<section class='plot'>"
        f"<h2>{_esc(title)}</h2>"
        f"<img src=\"{_img_b64(path)}\" alt=\"{_esc(title)}\">"
        f"<p>{_esc(caption)}</p>"
        f"<small>Fonte: {_esc(path)}</small>"
        "</section>"
        for title, path, caption in _all_report_images(
            cfg, heatmap_png, interactions_png, cluster_png, extra_images
        )
    )
    cluster_table = (
        "<h2>Atribuição de clusters</h2>"
        f"<table><tr><th>Ligante</th><th>Cluster</th></tr>{cluster_rows}</table>"
        if cluster_rows else ""
    )

    return f"""<!doctype html>
<html lang="pt-br"><head><meta charset="utf-8">
<title>Relatório HIP²LInterActomics</title>
<style>
@page{{size:A4 landscape;margin:14mm}}
*{{box-sizing:border-box}}
body{{font-family:'Segoe UI',Arial,sans-serif;max-width:1180px;margin:2em auto;padding:0 1.2em;color:#2a221d;background:#faf6f0}}
h1{{border-bottom:3px solid #0f766e;padding-bottom:.35em}}
h2{{color:#0f766e;margin-top:1.6em}}
table{{border-collapse:collapse;width:100%;margin:.6em 0;background:#fffdfa}}
th,td{{border:1px solid #ddd3c3;padding:7px 10px;font-size:13px}}
th{{background:#efe7db;text-align:left}} .number{{text-align:right}}
.meta,small{{color:#6a5d52;font-size:12px}}
.plot{{break-before:page;page-break-before:always;break-inside:avoid;page-break-inside:avoid}}
.plot img{{display:block;max-width:100%;max-height:155mm;width:auto;height:auto;object-fit:contain;margin:.6em auto}}
.plot p{{line-height:1.45;background:#fffdfa;border:1px solid #e7d9c5;padding:10px 12px}}
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
{images_html}
{cluster_table}
</body></html>
"""


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
    top_inter = sorted(inter_counts.items(), key=lambda item: -item[1])[:12]
    top_res = sorted(residue_counts.items(), key=lambda item: -item[1])[:12]
    return {
        "generated_at": datetime.now().strftime("%Y-%m-%d %H:%M"),
        "entries": analysis.get("entries", "não informado"),
        "cfg_rows": _cfg_rows_for_pdf(cfg),
        "top_inter": ", ".join(f"{key}: {value}" for key, value in top_inter) or "Sem interações contabilizadas.",
        "top_res": ", ".join(f"{key}: {value}" for key, value in top_res) or "Sem resíduos contabilizados.",
        "images": [(title, str(path), caption) for title, path, caption in _all_report_images(cfg, heatmap_png, interactions_png, cluster_png, extra_images)],
        "fp_rows": _fp_rows(fp_dashboards),
        "clusters": [(str(label), str(cluster_id)) for label, cluster_id in (clusters or [])[:80]],
    }


def _write_pdf_payload(path: str | Path, payload: dict) -> tuple[Path, list[str]]:
    from matplotlib.backends.backend_pdf import PdfPages

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
            for title, image_path, caption in payload["images"]:
                warning = _add_image_page(pdf, title, Path(image_path), caption, page_state)
                if warning:
                    warnings.append(warning)
            if payload["fp_rows"]:
                _add_text_page(pdf, "Resumo das análises de fingerprints", ["Cada linha resume uma base de fingerprints carregada e o corte usado na seleção de features."], payload["fp_rows"], page_state)
            if payload["clusters"]:
                _add_text_page(pdf, "Atribuição de clusters", ["Tabela dos ligantes e seus grupos hierárquicos."], payload["clusters"], page_state)
        temporary.replace(output)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise
    return output, warnings


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
    Path(status_path).write_text(json.dumps(payload, ensure_ascii=False), encoding="utf-8")


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
