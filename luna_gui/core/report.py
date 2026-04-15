"""Self-contained HTML report generator for a finished LUNA run."""
from __future__ import annotations

import base64
import html
from datetime import datetime
from pathlib import Path

from .project import ProjectConfig


def _esc(x) -> str:
    return html.escape(str(x))


def _img_b64(path: Path) -> str:
    if not path.exists():
        return ""
    return "data:image/png;base64," + base64.b64encode(path.read_bytes()).decode()


def build_report(
    cfg: ProjectConfig,
    analysis: dict,
    heatmap_png: Path | None = None,
    interactions_png: Path | None = None,
) -> str:
    """Return the HTML report as a string."""
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

    heatmap_html = ""
    if heatmap_png and heatmap_png.exists():
        heatmap_html = f'<h2>Matriz de similaridade</h2><img src="{_img_b64(heatmap_png)}" style="max-width:100%">'

    inter_html = ""
    if interactions_png and interactions_png.exists():
        inter_html = f'<h2>Distribuição de interações</h2><img src="{_img_b64(interactions_png)}" style="max-width:100%">'

    return f"""<!doctype html>
<html lang="pt-br"><head><meta charset="utf-8">
<title>Relatório LUNA</title>
<style>
body{{font-family:Segoe UI,Arial,sans-serif;max-width:960px;margin:2em auto;padding:0 1em;color:#222}}
h1{{border-bottom:2px solid #2c5aa0;padding-bottom:.3em}}
h2{{color:#2c5aa0;margin-top:1.6em}}
table{{border-collapse:collapse;width:100%;margin:.5em 0}}
th,td{{border:1px solid #ddd;padding:6px 10px;font-size:13px}}
th{{background:#f4f6fa;text-align:left}}
.meta{{color:#666;font-size:12px}}
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

<h2>Top 30 resíduos com mais interações</h2>
<table><tr><th>Cadeia/Resíduo/Num</th><th>Contagem</th></tr>{rows_res}</table>

{heatmap_html}

</body></html>
"""


def save_report(path: str | Path, **kwargs) -> Path:
    p = Path(path)
    p.write_text(build_report(**kwargs), encoding="utf-8")
    return p
