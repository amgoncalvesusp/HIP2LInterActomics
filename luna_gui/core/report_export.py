"""Self-contained HTML report generator for a finished LUNA run."""
from __future__ import annotations

import base64
import html
import textwrap
from datetime import datetime
from pathlib import Path

from .project import ProjectConfig
from ..i18n import t


def _esc(value) -> str:
    return html.escape(str(value))


def _img_b64(path: Path) -> str:
    if not path.exists():
        return ""
    return "data:image/png;base64," + base64.b64encode(path.read_bytes()).decode()


def build_report(
    cfg: ProjectConfig,
    analysis: dict,
    heatmap_png: Path | None = None,
    interactions_png: Path | None = None,
    cluster_png: Path | None = None,
    clusters: list[tuple[str, int]] | None = None,
) -> str:
    rows_cfg = "".join(
        f"<tr><td>{_esc(k)}</td><td>{_esc(v)}</td></tr>"
        for k, v in vars(cfg).items()
        if not (isinstance(v, list) and not v) and v not in ("", False, None)
    )
    inter_counts = analysis.get("interaction_counts", {})
    rows_inter = "".join(
        f"<tr><td>{_esc(k)}</td><td style='text-align:right'>{_esc(v)}</td></tr>"
        for k, v in sorted(inter_counts.items(), key=lambda x: -x[1])
    ) or "<tr><td colspan=2>—</td></tr>"

    res_counts = analysis.get("residue_counts", {})
    top_res = sorted(res_counts.items(), key=lambda x: -x[1])[:30]
    rows_res = "".join(
        f"<tr><td>{_esc(k)}</td><td style='text-align:right'>{_esc(v)}</td></tr>"
        for k, v in top_res
    ) or "<tr><td colspan=2>—</td></tr>"

    cluster_rows = ""
    if clusters:
        cluster_rows = "".join(
            f"<tr><td>{_esc(label)}</td><td style='text-align:right'>{cluster_id}</td></tr>"
            for label, cluster_id in clusters
        )

    cluster_table = ""
    if cluster_rows:
        cluster_table = (
            "<h2>Clusters</h2>"
            "<table><tr><th>Ligante</th><th>Cluster</th></tr>"
            f"{cluster_rows}</table>"
        )

    cluster_html = ""
    if cluster_png and cluster_png.exists():
        cluster_html = (
            "<h2>Dendrograma e matriz clusterizada</h2>"
            f"<img src=\"{_img_b64(cluster_png)}\" style=\"max-width:100%\">"
        )

    heatmap_html = ""
    if heatmap_png and heatmap_png.exists():
        heatmap_html = (
            "<h2>Matriz de similaridade</h2>"
            f"<img src=\"{_img_b64(heatmap_png)}\" style=\"max-width:100%\">"
        )

    inter_html = ""
    if interactions_png and interactions_png.exists():
        inter_html = (
            "<h2>Distribuição de interações</h2>"
            f"<img src=\"{_img_b64(interactions_png)}\" style=\"max-width:100%\">"
        )

    return f"""<!doctype html>
<html lang="pt-br"><head><meta charset="utf-8">
<title>Relatório LUNA</title>
<style>
body{{font-family:'Segoe UI',Arial,sans-serif;max-width:1080px;margin:2em auto;padding:0 1.2em;color:#2a221d;background:#faf6f0}}
h1{{border-bottom:3px solid #0f766e;padding-bottom:.35em}}
h2{{color:#0f766e;margin-top:1.6em}}
table{{border-collapse:collapse;width:100%;margin:.6em 0;background:#fffdfa}}
th,td{{border:1px solid #ddd3c3;padding:7px 10px;font-size:13px}}
th{{background:#efe7db;text-align:left}}
.meta{{color:#6a5d52;font-size:12px}}
</style></head><body>
<h1>Relatório LUNA</h1>
<p class="meta">Gerado em {datetime.now().strftime('%Y-%m-%d %H:%M')}</p>

<h2>Resumo</h2>
<ul>
  <li><b>Proteína:</b> {_esc(cfg.protein_file)}</li>
  <li><b>Ligantes:</b> {_esc(cfg.ligand_file)}</li>
  <li><b>Total selecionado:</b> {len(cfg.selected_ligands)}</li>
  <li><b>Workdir:</b> {_esc(cfg.workdir)}</li>
  <li><b>Entries processadas:</b> {_esc(analysis.get('entries', '—'))}</li>
</ul>

<h2>Configuração</h2>
<table><tr><th>Parâmetro</th><th>Valor</th></tr>{rows_cfg}</table>

<h2>Contagem por tipo de interação</h2>
<table><tr><th>Tipo</th><th>Total</th></tr>{rows_inter}</table>

{inter_html}
{heatmap_html}
{cluster_html}
{cluster_table}

<h2>Top 30 resíduos com mais interações</h2>
<table><tr><th>Cadeia/Resíduo/Num</th><th>Contagem</th></tr>{rows_res}</table>

</body></html>
"""


def save_report(path: str | Path, **kwargs) -> Path:
    path = Path(path)
    path.write_text(build_report(**kwargs), encoding="utf-8")
    return path


def _cfg_rows_for_pdf(cfg: ProjectConfig) -> list[tuple[str, str]]:
    rows = [
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
        ("Otsu fallback", "Sim" if getattr(cfg, "fp_use_otsu_threshold", False) else "Não"),
        ("Núcleos", str(cfg.nproc)),
    ]
    return rows


def _add_text_page(pdf, title: str, paragraphs: list[str], rows: list[tuple[str, str]] | None = None) -> None:
    from matplotlib.figure import Figure
    from matplotlib.patches import Rectangle

    title = t(title)
    paragraphs = [t(paragraph) for paragraph in paragraphs]
    rows = [(t(key), t(value)) for key, value in (rows or [])] or None
    content_top = 0.885
    content_bottom = 0.075

    def new_page():
        fig = Figure(figsize=(8.27, 11.69), dpi=150)
        fig.patch.set_facecolor("#fbf7ef")
        ax = fig.add_subplot(111)
        ax.axis("off")
        ax.add_patch(
            Rectangle(
                (0.0, 0.925),
                1.0,
                0.075,
                transform=ax.transAxes,
                facecolor="#145c58",
                edgecolor="none",
            )
        )
        ax.text(0.06, 0.965, title, fontsize=16.5, weight="bold", va="top", color="white")
        return fig, ax, content_top

    fig, ax, y = new_page()

    def save_and_new_page():
        nonlocal fig, ax, y
        pdf.savefig(fig, dpi=220)
        fig, ax, y = new_page()

    def ensure_space(required_height: float) -> None:
        if y - required_height < content_bottom:
            save_and_new_page()

    def draw_paragraph(paragraph: str) -> None:
        nonlocal y
        lines = textwrap.wrap(str(paragraph), width=84, break_long_words=False, break_on_hyphens=False) or [""]
        line_height = 0.023
        block_height = (line_height * len(lines)) + 0.021
        ensure_space(block_height)
        for line in lines:
            ax.text(0.06, y, line, fontsize=9.3, va="top", color="#2d251e")
            y -= line_height
        y -= 0.021

    def draw_table(rows_to_draw: list[tuple[str, str]]) -> None:
        nonlocal y
        if not rows_to_draw:
            return
        ensure_space(0.055)
        ax.add_patch(
            Rectangle(
                (0.055, y - 0.032),
                0.89,
                0.036,
                transform=ax.transAxes,
                facecolor="#e7d9c5",
                edgecolor="#d7c7b2",
                linewidth=0.35,
            )
        )
        ax.text(0.065, y - 0.006, t("Parâmetro"), fontsize=8.9, weight="bold", va="top", color="#145c58")
        ax.text(0.325, y - 0.006, t("Valor"), fontsize=8.9, weight="bold", va="top", color="#145c58")
        y -= 0.046

        for row_index, (key, value) in enumerate(rows_to_draw):
            key_lines = textwrap.wrap(str(key), width=19, break_long_words=True, break_on_hyphens=False) or [""]
            value_lines = textwrap.wrap(str(value), width=61, break_long_words=True, break_on_hyphens=False) or [""]
            line_count = max(len(key_lines), len(value_lines))
            line_height = 0.0215
            row_height = (line_height * line_count) + 0.020
            ensure_space(row_height + 0.006)
            face = "#fffdf8" if row_index % 2 == 0 else "#f7efe4"
            ax.add_patch(
                Rectangle(
                    (0.055, y - row_height + 0.004),
                    0.89,
                    row_height,
                    transform=ax.transAxes,
                    facecolor=face,
                    edgecolor="#eadfce",
                    linewidth=0.35,
                )
            )
            text_y = y - 0.008
            for line_idx, line in enumerate(key_lines):
                ax.text(0.065, text_y - (line_idx * line_height), line, fontsize=8.4, weight="bold", va="top", color="#145c58")
            for line_idx, line in enumerate(value_lines):
                ax.text(0.325, text_y - (line_idx * line_height), line, fontsize=8.4, va="top", color="#2d251e")
            y -= row_height + 0.006

    for paragraph in paragraphs:
        draw_paragraph(paragraph)
    if rows:
        y -= 0.004
        draw_table(rows)
    pdf.savefig(fig, dpi=220)


def _add_image_page(pdf, title: str, image_path: Path, caption: str) -> None:
    from matplotlib.figure import Figure
    from matplotlib.patches import Rectangle
    import matplotlib.image as mpimg

    if not image_path.exists():
        return
    title = t(title)
    caption = t(caption)
    fig = Figure(figsize=(13.2, 9.3), dpi=220)
    fig.patch.set_facecolor("#fbf7ef")
    title_ax = fig.add_axes([0.045, 0.905, 0.91, 0.06])
    title_ax.axis("off")
    title_ax.text(0.0, 0.9, title, fontsize=16.5, weight="bold", color="#145c58", va="top")
    ax = fig.add_axes([0.045, 0.20, 0.91, 0.67])
    ax.axis("off")
    img = mpimg.imread(str(image_path))
    ax.imshow(img, interpolation="none", resample=False)
    caption_ax = fig.add_axes([0.045, 0.045, 0.91, 0.12])
    caption_ax.axis("off")
    caption_ax.add_patch(
        Rectangle(
            (0.0, 0.0),
            1.0,
            1.0,
            transform=caption_ax.transAxes,
            facecolor="#fffdf8",
            edgecolor="#e7d9c5",
            linewidth=0.7,
        )
    )
    caption_lines = textwrap.wrap(str(caption), width=150, break_long_words=False, break_on_hyphens=False)
    caption_y = 0.82
    for line in caption_lines[:5]:
        caption_ax.text(0.018, caption_y, line, fontsize=9.0, color="#2d251e", va="top")
        caption_y -= 0.18
    pdf.savefig(fig, dpi=600)


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
    """Generate a compact PDF report with parameters and interpretation notes."""
    from matplotlib.backends.backend_pdf import PdfPages

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    inter_counts = analysis.get("interaction_counts", {}) or {}
    top_inter = sorted(inter_counts.items(), key=lambda item: -item[1])[:12]
    residue_counts = analysis.get("residue_counts", {}) or {}
    top_res = sorted(residue_counts.items(), key=lambda item: -item[1])[:12]
    top_inter_text = ", ".join(f"{k}: {v}" for k, v in top_inter) or "Sem interações contabilizadas."
    top_res_text = ", ".join(f"{k}: {v}" for k, v in top_res) or "Sem resíduos contabilizados."

    with PdfPages(path) as pdf:
        _add_text_page(
            pdf,
            "Relatório HIP²LInterActomics",
            [
                f"Gerado em {datetime.now().strftime('%Y-%m-%d %H:%M')}.",
                "Este relatório resume os resultados carregados na aba 5.Resultados, os parâmetros usados no projeto e uma leitura guiada dos gráficos. Use-o como material de triagem: padrões fortes indicam hipóteses para inspeção estrutural, não uma conclusão automática de afinidade.",
                f"Entradas processadas: {analysis.get('entries', 'não informado')}.",
            ],
            _cfg_rows_for_pdf(cfg),
        )
        _add_text_page(
            pdf,
            "Como interpretar as análises",
            [
                "Estatísticas de interação: contam quantas vezes cada tipo de contato foi observado. Barras altas indicam mecanismos recorrentes, como ligações de hidrogênio, contatos hidrofóbicos ou interações iônicas. Compare tipos dominantes com resíduos frequentes para separar padrões químicos reais de ruído de pose.",
                "Mapa de calor por tipo: cruza ligantes e resíduos para uma interação escolhida. Células mais intensas indicam mais ocorrências daquele contato. Colunas densas sugerem resíduos-chave; linhas densas sugerem ligantes com muitos contatos daquele tipo.",
                "Mapa de calor completo ligantes x resíduos: usa faixas coloridas para mostrar múltiplos tipos de interação no mesmo par ligante-resíduo. Ele é útil para enxergar complementaridade química: o mesmo resíduo pode estabilizar ligantes por mecanismos diferentes.",
                "Matriz de similaridade: compara ligantes pelos fingerprints de interação. Valores próximos de 1 indicam perfis de interação semelhantes; valores baixos indicam modos de interação distintos, mesmo quando as moléculas parecem estruturalmente parecidas.",
                "Clusters: reorganizam a matriz de similaridade para revelar famílias de ligantes por comportamento no sítio. Use os grupos como hipótese para priorização e para escolher representantes para inspeção no PyMOL.",
                "Fingerprints e FP análises: cada feature resume uma vizinhança de interação. A classe atribuída descreve a natureza dominante da feature; a importância do modelo estima quanto ela ajuda a separar classes/rótulos ou valores de atividade.",
                "Importância e p-value: o z-score de importância compara uma feature contra a distribuição de importâncias do conjunto. O p-value é calculado pela equação de Keiser and Hert [1], p = 1 - exp(-exp(((-z*pi)/sqrt(6)) - 0.577215665)).",
                "Otsu's Thresholding: quando nenhuma feature passa pelo critério z-score > 1, Otsu define um limiar alternativo baseado na separação da distribuição de percentuais. Isso evita aceitar apenas casos 100% prevalentes quando a base tem padrões intermediários.",
                "Sessões PyMOL: permitem validar visualmente se as features ou filtros representam contatos plausíveis no complexo. A filtragem dinâmica pode gerar sessões novas ou, se o projeto salvo não reabrir, copiar sessões existentes compatíveis com a matriz cacheada.",
            ],
            [
                ("Interações mais frequentes", top_inter_text),
                ("Resíduos mais frequentes", top_res_text),
            ],
        )
        if interactions_png:
            _add_image_page(
                pdf,
                "Distribuição de interações",
                Path(interactions_png),
                "Cada barra representa a contagem de uma classe de interação. Use barras dominantes para identificar forças químicas recorrentes e barras raras para procurar contatos específicos que podem diferenciar poucos ligantes.",
            )
        if heatmap_png:
            _add_image_page(
                pdf,
                "Matriz de similaridade",
                Path(heatmap_png),
                "Cada célula compara dois ligantes pelo fingerprint de interação. Tons mais intensos indicam maior similaridade; blocos ao longo da diagonal sugerem famílias com modos de interação semelhantes.",
            )
        if cluster_png:
            _add_image_page(
                pdf,
                "Clusters",
                Path(cluster_png),
                "O dendrograma mostra a distância entre perfis de interação e a matriz reordenada evidencia grupos. Clusters compactos sugerem ligantes que compartilham padrões de contato e podem ser priorizados em conjunto.",
            )
        for title, image_path, caption in extra_images or []:
            _add_image_page(pdf, str(title), Path(image_path), str(caption))
        fp_rows: list[tuple[str, str]] = []
        for key, dashboard in (fp_dashboards or {}).items():
            if not isinstance(dashboard, dict):
                continue
            fp_rows.append(
                (
                    str(key),
                    (
                        f"features={len(dashboard.get('features', []) or [])}; "
                        f"importantes={len(dashboard.get('important_features', []) or [])}; "
                        f"modelo={dashboard.get('model_name', '-')}; "
                        f"limiar={float(dashboard.get('threshold_pct', 0.0) or 0.0):.2f}%"
                    ),
                )
            )
        if fp_rows:
            _add_text_page(
                pdf,
                "Resumo das FP análises",
                [
                    "Cada linha resume uma base de fingerprints carregada. Se o modelo aparecer como fallback ou indisponível, a interpretação deve ser tratada como exploratória e a interface mostra a causa no campo de método.",
                ],
                fp_rows[:30],
            )
        if clusters:
            _add_text_page(
                pdf,
                "Atribuição de clusters",
                ["Tabela com os primeiros ligantes e seus clusters hierárquicos."],
                [(label, str(cluster_id)) for label, cluster_id in clusters[:80]],
            )
    return path
