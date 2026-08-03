"""Post-analysis exports for the HIP2LInterActomics terminal workflow.

The module writes portable artifacts from the same workdir consumed by the GUI.
Fingerprint dashboards, static figures, cluster assignments, a lightweight
interactive cluster explorer, and optional FP-PyMOL sessions are exported.
"""
from __future__ import annotations

import json
from datetime import datetime, timezone
from html import escape
from pathlib import Path
from typing import Any, Mapping

import numpy as np

from . import analysis_runtime, results_analysis
from .project import IFP_SUFFIXES, ProjectConfig, resolve_sim_matrix_output_paths
from ..i18n import t


_PLOT_DPI = 300


def _int_setting(
    settings: Mapping[str, Any],
    key: str,
    default: int,
    minimum: int = 1,
) -> int:
    try:
        return max(minimum, int(settings.get(key, default)))
    except (TypeError, ValueError):
        return default


def _relative_to_workdir(path: Path, workdir: Path) -> str:
    try:
        return path.resolve().relative_to(workdir.resolve()).as_posix()
    except (OSError, ValueError):
        return str(path)


def _write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")


def _safe_token(value: str) -> str:
    return "".join(char if char.isalnum() else "_" for char in value).strip("_") or "result"


def _get_pyplot():
    try:
        import matplotlib

        matplotlib.use("Agg", force=True)
        import matplotlib.pyplot as plt
    except Exception:
        return None
    return plt


def _rendered_entry_count(entry_count: int) -> int:
    return min(200, max(1, int(entry_count)))


def _save_plot(fig, output: Path, plt) -> Path:
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=_PLOT_DPI, bbox_inches="tight", pad_inches=0.18)
    plt.close(fig)
    return output


def _save_interaction_summary(summary: dict, output_dir: Path, language: str = "pt") -> Path | None:
    counts = summary.get("interaction_counts") if isinstance(summary, dict) else None
    if not isinstance(counts, dict):
        return None

    rows = [
        (str(name), float(value))
        for name, value in counts.items()
        if str(name).strip()
    ]
    if not rows:
        return None
    rows.sort(key=lambda item: item[1])

    plt = _get_pyplot()
    if plt is None:
        return None
    height = max(4.0, 0.28 * len(rows) + 1.8)
    fig, axis = plt.subplots(figsize=(9.0, height))
    axis.barh([name for name, _ in rows], [value for _, value in rows], color="#2f7f83")
    axis.set_xlabel(t("Contagem de interações", lang=language))
    axis.set_title(t("Resumo de interações", lang=language))
    axis.grid(axis="x", alpha=0.2)
    fig.tight_layout()
    output = output_dir / "interaction_summary.png"
    return _save_plot(fig, output, plt)


def _save_residue_heatmaps(residue_artifact: dict, output_dir: Path, language: str = "pt") -> list[Path]:
    entries = list(residue_artifact.get("entries", []) or [])
    residues = list(residue_artifact.get("residues", []) or [])
    matrices = residue_artifact.get("matrix")
    if not entries or not residues or not isinstance(matrices, dict):
        return []

    plt = _get_pyplot()
    if plt is None:
        return []

    outputs: list[Path] = []
    rendered_entries = _rendered_entry_count(len(entries))
    height = max(4.5, 3.0 + 0.035 * rendered_entries)
    width = max(8.0, 4.0 + 0.10 * min(len(residues), 60))
    for interaction_type, values in matrices.items():
        if not values:
            continue
        try:
            row_count = len(values)
            column_count = len(values[0]) if row_count else 0
        except TypeError:
            continue
        if row_count != len(entries) or column_count != len(residues):
            continue

        fig, axis = plt.subplots(figsize=(width, height))
        image = axis.imshow(values, cmap="viridis", aspect="auto")
        axis.set_title(f"{t('Mapa de interação por resíduo', lang=language)}: {t(interaction_type, lang=language)}")
        axis.set_xlabel(t("Resíduos", lang=language))
        axis.set_ylabel(f"{t('Ligantes', lang=language)} ({len(entries)} {t('total', lang=language)})")
        if len(residues) <= 40:
            axis.set_xticks(range(len(residues)))
            axis.set_xticklabels(residues, rotation=90, fontsize=7)
        if len(entries) <= 40:
            axis.set_yticks(range(len(entries)))
            axis.set_yticklabels(entries, fontsize=7)
        fig.colorbar(image, ax=axis, fraction=0.046, pad=0.04)
        fig.tight_layout()
        output = output_dir / f"residue_map_{_safe_token(str(interaction_type))}.png"
        outputs.append(_save_plot(fig, output, plt))
    return outputs


def _save_complete_residue_heatmap(
    residue_artifact: dict,
    output_dir: Path,
    language: str = "pt",
) -> Path | None:
    entries, residues, matrix, _interaction_types = results_analysis.build_complete_heatmap(residue_artifact)
    values = np.asarray(matrix, dtype=float)
    if not entries or not residues or values.size == 0:
        return None
    plt = _get_pyplot()
    if plt is None:
        return None
    fig, axis = plt.subplots(
        figsize=(
            max(9.0, 4.0 + 0.12 * min(len(residues), 90)),
            max(6.0, 3.2 + 0.04 * _rendered_entry_count(len(entries))),
        )
    )
    image = axis.imshow(values, cmap="viridis", aspect="auto")
    axis.set_title(t("Mapa de calor completo ligantes x resíduos", lang=language))
    axis.set_xlabel(t("Resíduos", lang=language))
    axis.set_ylabel(f"{t('Ligantes', lang=language)} ({len(entries)} {t('total', lang=language)})")
    if len(residues) <= 50:
        axis.set_xticks(range(len(residues)))
        axis.set_xticklabels(residues, rotation=90, fontsize=7)
    if len(entries) <= 40:
        axis.set_yticks(range(len(entries)))
        axis.set_yticklabels(entries, fontsize=7)
    fig.colorbar(image, ax=axis, fraction=0.035, pad=0.03)
    fig.tight_layout()
    return _save_plot(fig, output_dir / "complete_ligands_residues_heatmap.png", plt)


def _fp_plot_features(dashboard: dict, limit: int = 50) -> list[dict]:
    rows = list(dashboard.get("important_features", []) or [])
    if not rows:
        rows = list(dashboard.get("features", []) or [])
    return rows[: max(1, int(limit))]


def _save_fp_dashboard_figures(
    dashboard: dict,
    output_dir: Path,
    ifp_type: str,
    language: str,
) -> list[Path]:
    plt = _get_pyplot()
    if plt is None:
        return []
    suffix = IFP_SUFFIXES.get(ifp_type, ifp_type)
    target = output_dir / "fingerprints" / str(suffix)
    target.mkdir(parents=True, exist_ok=True)
    outputs: list[Path] = []
    features = _fp_plot_features(dashboard)
    colors = dashboard.get("class_colors", {}) or {}

    class_share = dashboard.get("class_share", {}) or {}
    class_labels = [label for label in results_analysis.FP_CLASS_ORDER if float(class_share.get(label, 0.0) or 0.0) > 0.0]
    if class_labels:
        fig, axis = plt.subplots(figsize=(max(9.0, 1.0 * len(class_labels) + 3.5), 6.0))
        values = [float(class_share[label]) for label in class_labels]
        axis.bar(range(len(class_labels)), values, color=[colors.get(label, "#2f7f83") for label in class_labels])
        axis.set_xticks(range(len(class_labels)))
        axis.set_xticklabels([t(label, lang=language) for label in class_labels], rotation=30, ha="right")
        axis.set_ylabel(t("% de features importantes", lang=language))
        axis.set_title(t("Distribuição das classes entre as features mais importantes", lang=language))
        fig.tight_layout()
        outputs.append(_save_plot(fig, target / "fp_class_summary.png", plt))

    if features:
        labels = [str(feature.get("feature_id", "-")) for feature in features]
        y = np.arange(len(features))
        fig, axis = plt.subplots(figsize=(11.5, max(6.5, 0.30 * len(features) + 2.5)))
        left = np.zeros(len(features), dtype=float)
        for class_name in results_analysis.FP_CLASS_ORDER:
            values = np.asarray(
                [float((feature.get("class_percentages") or {}).get(class_name, 0.0) or 0.0) for feature in features]
            )
            if not np.any(values > 0):
                continue
            axis.barh(y, values, left=left, label=t(class_name, lang=language), color=colors.get(class_name, "#6f9ec7"))
            left += values
        axis.set_yticks(y)
        axis.set_yticklabels(labels)
        axis.invert_yaxis()
        axis.set_xlim(0, 100)
        axis.set_xlabel(t("Frequência de atribuição de cada classe (%)", lang=language))
        axis.set_ylabel(t("ID da feature", lang=language))
        axis.set_title(t("Frequência de atribuição de classes nas features importantes", lang=language))
        axis.legend(loc="upper center", bbox_to_anchor=(0.5, -0.08), ncol=3, frameon=False)
        fig.tight_layout()
        outputs.append(_save_plot(fig, target / "fp_class_assignment.png", plt))

        fig, axis = plt.subplots(figsize=(11.5, max(6.5, 0.30 * len(features) + 2.5)))
        coverage = [float(feature.get("coverage_pct", 0.0) or 0.0) for feature in features]
        importance = [float(feature.get("importance_pct", 0.0) or 0.0) for feature in features]
        axis.barh(y, coverage, color=[colors.get(str(feature.get("assigned_class", "")), "#2f7f83") for feature in features])
        axis.scatter(importance, y, marker="*", color="#b42318", s=45, label=t("Importância relativa", lang=language))
        axis.set_yticks(y)
        axis.set_yticklabels(labels)
        axis.invert_yaxis()
        axis.set_xlim(0, 105)
        axis.set_xlabel(t("Cobertura / importância relativa (%)", lang=language))
        axis.set_ylabel(t("ID da feature", lang=language))
        axis.set_title(t("Cobertura das features importantes e importância do modelo", lang=language))
        axis.legend(frameon=False)
        fig.tight_layout()
        outputs.append(_save_plot(fig, target / "fp_coverage_importance.png", plt))

        entries = list(dashboard.get("entry_labels", []) or [])
        if entries:
            presence = np.asarray(
                [[1.0 if entry in (feature.get("entry_counts") or {}) else 0.0 for feature in features] for entry in entries],
                dtype=float,
            )
            fig, axis = plt.subplots(figsize=(max(9.0, 0.30 * len(features) + 4.0), max(6.5, 0.06 * _rendered_entry_count(len(entries)) + 3.0)))
            image = axis.imshow(presence, cmap="Blues", aspect="auto", vmin=0, vmax=1)
            axis.set_xticks(range(len(features)))
            axis.set_xticklabels(labels, rotation=90, fontsize=7)
            if len(entries) <= 40:
                axis.set_yticks(range(len(entries)))
                axis.set_yticklabels(entries, fontsize=7)
            axis.set_xlabel(t("ID da feature", lang=language))
            axis.set_ylabel(t("Ligantes", lang=language))
            axis.set_title(t("Mapa de presença das features importantes por classe", lang=language))
            fig.colorbar(image, ax=axis, fraction=0.03, pad=0.02)
            fig.tight_layout()
            outputs.append(_save_plot(fig, target / "fp_feature_presence_heatmap.png", plt))

        interaction_names = []
        for feature in features:
            for interaction_name, count in (feature.get("interaction_breakdown") or {}).items():
                if int(count or 0) > 0 and interaction_name not in interaction_names:
                    interaction_names.append(str(interaction_name))
        interaction_names = sorted(interaction_names, key=results_analysis.interaction_priority_key)[:12]
        if interaction_names:
            fig, axis = plt.subplots(figsize=(11.5, max(6.5, 0.30 * len(features) + 2.5)))
            left = np.zeros(len(features), dtype=float)
            for interaction_name in interaction_names:
                values = []
                for feature in features:
                    breakdown = feature.get("interaction_breakdown") or {}
                    total = float(sum(int(value or 0) for value in breakdown.values()))
                    values.append(100.0 * float(breakdown.get(interaction_name, 0) or 0) / total if total else 0.0)
                values_array = np.asarray(values, dtype=float)
                if not np.any(values_array > 0):
                    continue
                axis.barh(y, values_array, left=left, label=t(interaction_name, lang=language), color=results_analysis.get_interaction_color(interaction_name))
                left += values_array
            axis.set_yticks(y)
            axis.set_yticklabels(labels)
            axis.invert_yaxis()
            axis.set_xlim(0, 100)
            axis.set_xlabel(t("Frequência de atribuição de cada interação (%)", lang=language))
            axis.set_ylabel(t("ID da feature", lang=language))
            axis.set_title(t("Frequência de atribuição da interação prevalente nas features importantes", lang=language))
            axis.legend(loc="upper center", bbox_to_anchor=(0.5, -0.08), ncol=3, frameon=False)
            fig.tight_layout()
            outputs.append(_save_plot(fig, target / "fp_interaction_assignment.png", plt))

        prevalent = [feature for feature in features if str(feature.get("prevalent_interaction", "")) not in {"", results_analysis.CLASS_UNRELIABLE}]
        if prevalent:
            prevalent_labels = [str(feature.get("feature_id", "-")) for feature in prevalent]
            heights = [len(feature.get("prevalent_pair_entries", []) or []) for feature in prevalent]
            fig, axis = plt.subplots(figsize=(max(10.0, 0.48 * len(prevalent) + 4.0), 6.2))
            axis.bar(
                range(len(prevalent)),
                heights,
                color=[results_analysis.get_interaction_color(str(feature.get("prevalent_interaction", ""))) for feature in prevalent],
            )
            axis.set_xticks(range(len(prevalent)))
            axis.set_xticklabels(prevalent_labels, rotation=45, ha="right")
            axis.set_xlabel(t("ID da feature", lang=language))
            axis.set_ylabel(t("Número de ligantes", lang=language))
            axis.set_title(t("Interação e resíduo prevalentes nas features importantes", lang=language))
            fig.tight_layout()
            outputs.append(_save_plot(fig, target / "fp_prevalent_interactions.png", plt))

            if entries:
                prevalence_matrix = np.asarray(
                    [[1.0 if entry in (feature.get("prevalent_pair_entries") or []) else 0.0 for feature in prevalent] for entry in entries],
                    dtype=float,
                )
                fig, axis = plt.subplots(figsize=(max(9.0, 0.30 * len(prevalent) + 4.0), max(6.5, 0.06 * _rendered_entry_count(len(entries)) + 3.0)))
                image = axis.imshow(prevalence_matrix, cmap="YlGnBu", aspect="auto", vmin=0, vmax=1)
                axis.set_xticks(range(len(prevalent)))
                axis.set_xticklabels(prevalent_labels, rotation=90, fontsize=7)
                if len(entries) <= 40:
                    axis.set_yticks(range(len(entries)))
                    axis.set_yticklabels(entries, fontsize=7)
                axis.set_xlabel(t("ID da feature", lang=language))
                axis.set_ylabel(t("Ligantes", lang=language))
                axis.set_title(t("Mapa de calor de interações prevalentes", lang=language))
                fig.colorbar(image, ax=axis, fraction=0.03, pad=0.02)
                fig.tight_layout()
                outputs.append(_save_plot(fig, target / "fp_prevalent_interactions_heatmap.png", plt))

    for model_key, model_title in (("extra_trees", "Extra Trees"), ("gradient_boosting", "Gradient Boosting")):
        rows = list((dashboard.get("top_features_by_model") or {}).get(model_key, []) or [])[:50]
        if not rows:
            continue
        rows = list(reversed(rows))
        labels = [f"{row.get('feature_id', '-')} (L{row.get('assigned_level') or '-'})" for row in rows]
        values = [float(row.get("importance_score", 0.0) or 0.0) for row in rows]
        fig, axis = plt.subplots(figsize=(11.5, max(7.0, 0.28 * len(rows) + 2.7)))
        axis.barh(range(len(rows)), values, color="#2f7f83" if model_key == "extra_trees" else "#b45f3b")
        axis.set_yticks(range(len(rows)))
        axis.set_yticklabels(labels, fontsize=7.5)
        axis.set_xlabel(t("Importância da feature", lang=language))
        axis.set_ylabel(t("ID da feature e nível", lang=language))
        axis.set_title(t(f"Top 50 features por {model_title}", lang=language))
        fig.tight_layout()
        outputs.append(_save_plot(fig, target / f"fp_top50_{model_key}.png", plt))

    return outputs


def _similarity_candidates(workdir: Path, cfg: ProjectConfig, ifp_type: str) -> list[Path]:
    suffix = IFP_SUFFIXES.get(ifp_type, "E")
    configured = resolve_sim_matrix_output_paths(cfg).get(ifp_type, "")
    candidates: list[Path] = []
    if configured:
        candidate = Path(configured)
        if not candidate.is_absolute():
            candidate = workdir / candidate
        candidates.append(candidate)

    candidates.append(workdir / f"sim_matrix_{suffix}.csv")
    unique: list[Path] = []
    seen: set[str] = set()
    for candidate in candidates:
        square = candidate
        if not candidate.stem.endswith("_square"):
            square = candidate.with_name(f"{candidate.stem}_square{candidate.suffix}")
        for item in (square, candidate):
            key = str(item).lower()
            if key not in seen:
                seen.add(key)
                unique.append(item)
    return unique


def _square_matrix_size_hint(path: Path) -> int | None:
    try:
        with path.open("r", encoding="utf-8", errors="replace", newline="") as handle:
            for line in handle:
                if not line.strip():
                    continue
                cells = [cell.strip().lower() for cell in line.rstrip("\r\n").split(",")]
                if cells[:3] == ["entry1", "entry2", "similarity"]:
                    return None
                if len(cells) > 1:
                    return len(cells) - 1
                return None
    except OSError:
        return None
    return None


def _save_similarity_figure(
    labels: list[str],
    matrix,
    output_dir: Path,
    ifp_type: str,
    language: str = "pt",
) -> Path | None:
    plt = _get_pyplot()
    if plt is None:
        return None
    rendered_entries = _rendered_entry_count(len(labels))
    width = max(8.0, 4.0 + 0.035 * rendered_entries)
    height = max(7.0, 3.0 + 0.035 * rendered_entries)
    fig, axis = plt.subplots(figsize=(width, height))
    image = axis.imshow(matrix, cmap="magma", aspect="auto", vmin=0.0, vmax=1.0)
    axis.set_title(f"{t('Matriz de similaridade', lang=language)}: {ifp_type}")
    axis.set_xlabel(f"{t('Ligantes', lang=language)} ({len(labels)} {t('total', lang=language)})")
    axis.set_ylabel(f"{t('Ligantes', lang=language)} ({len(labels)} {t('total', lang=language)})")
    if len(labels) <= 40:
        axis.set_xticks(range(len(labels)))
        axis.set_xticklabels(labels, rotation=90, fontsize=7)
        axis.set_yticks(range(len(labels)))
        axis.set_yticklabels(labels, fontsize=7)
    fig.colorbar(image, ax=axis, fraction=0.046, pad=0.04)
    fig.tight_layout()
    output = output_dir / f"similarity_{IFP_SUFFIXES.get(ifp_type, ifp_type)}.png"
    return _save_plot(fig, output, plt)


def _save_cluster_figure(
    result,
    output_dir: Path,
    ifp_type: str,
    language: str = "pt",
) -> Path | None:
    plt = _get_pyplot()
    if plt is None:
        return None
    try:
        from scipy.cluster.hierarchy import dendrogram
    except Exception:
        return None

    rendered_entries = _rendered_entry_count(len(result.labels))
    width = max(9.0, 4.5 + 0.035 * rendered_entries)
    height = max(8.0, 4.2 + 0.035 * rendered_entries)
    fig = plt.figure(figsize=(width, height))
    grid = fig.add_gridspec(2, 1, height_ratios=[1.2, 2.3], hspace=0.35)
    tree_axis = fig.add_subplot(grid[0])
    labels = result.labels if len(result.labels) <= 40 else None
    dendrogram(result.linkage_matrix, labels=labels, leaf_rotation=90, leaf_font_size=7, ax=tree_axis)
    if labels is None:
        tree_axis.set_xticks([])
    tree_axis.set_title(
        f"{t('Clustering hierárquico', lang=language)}: {ifp_type} "
        f"({result.n_clusters} clusters, {result.method})"
    )
    tree_axis.set_ylabel(t("Distância", lang=language))

    matrix_axis = fig.add_subplot(grid[1])
    image = matrix_axis.imshow(result.ordered_matrix, cmap="magma", aspect="auto", vmin=0.0, vmax=1.0)
    matrix_axis.set_title(t("Matriz reordenada por cluster", lang=language))
    matrix_axis.set_xlabel(f"{t('Ligantes', lang=language)} ({len(result.labels)} {t('total', lang=language)})")
    matrix_axis.set_ylabel(f"{t('Ligantes', lang=language)} ({len(result.labels)} {t('total', lang=language)})")
    if len(result.ordered_labels) <= 40:
        matrix_axis.set_xticks(range(len(result.ordered_labels)))
        matrix_axis.set_xticklabels(result.ordered_labels, rotation=90, fontsize=7)
        matrix_axis.set_yticks(range(len(result.ordered_labels)))
        matrix_axis.set_yticklabels(result.ordered_labels, fontsize=7)
    fig.colorbar(image, ax=matrix_axis, fraction=0.046, pad=0.04)
    fig.subplots_adjust(left=0.09, right=0.90, top=0.93, bottom=0.08, hspace=0.38)
    output = output_dir / f"clusters_{IFP_SUFFIXES.get(ifp_type, ifp_type)}.png"
    return _save_plot(fig, output, plt)


def _write_cluster_explorer(
    result,
    output_dir: Path,
    ifp_type: str,
    include_matrix: bool,
    matrix_png: Path | None,
    assignments_csv: Path,
) -> Path:
    rows = []
    for leaf_order, item_index in enumerate(result.leaves, start=1):
        rows.append(
            {
                "ligand_id": result.labels[item_index],
                "cluster_id": int(result.cluster_ids[item_index]),
                "leaf_order": leaf_order,
            }
        )
    payload = {
        "ifp_type": ifp_type,
        "rows": rows,
        "matrix": result.ordered_matrix.tolist() if include_matrix else None,
        "matrix_available": bool(include_matrix),
        "matrix_png": matrix_png.name if matrix_png else "",
        "assignments_csv": assignments_csv.name,
    }
    payload_text = json.dumps(payload, ensure_ascii=False).replace("</", "<\\/")
    title = escape(f"{ifp_type} interactive cluster explorer")
    matrix_note = (
        "The complete matrix is embedded for hover inspection."
        if include_matrix
        else "The complete matrix remains in the workdir CSV; the explorer keeps all cluster assignments interactive."
    )
    html = """<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>__TITLE__</title>
<style>
body { font-family: Arial, sans-serif; margin: 24px; color: #17212b; }
h1 { margin: 0 0 6px; font-size: 22px; }
p { max-width: 1000px; }
.controls { display: flex; flex-wrap: wrap; gap: 12px; align-items: center; margin: 18px 0; }
input, select { padding: 7px; border: 1px solid #9aa8b5; border-radius: 4px; }
#cluster-summary { display: flex; flex-wrap: wrap; gap: 8px; margin: 12px 0; }
.cluster-chip { padding: 5px 9px; border: 1px solid #57748a; background: #eef5f8; border-radius: 3px; cursor: pointer; }
.canvas-wrap { max-width: 760px; overflow: auto; border: 1px solid #b9c5ce; padding: 10px; }
canvas { display: block; background: #ffffff; }
table { border-collapse: collapse; width: 100%; max-width: 1100px; margin-top: 14px; }
th, td { text-align: left; border-bottom: 1px solid #d9e1e7; padding: 7px; }
th { background: #edf3f6; }
.small { color: #52616d; font-size: 13px; }
</style>
</head>
<body>
<h1>__TITLE__</h1>
<p>__MATRIX_NOTE__</p>
<div class="controls">
<label>Search <input id="search" type="search" placeholder="Ligand ID"></label>
<label>Cluster <select id="cluster"><option value="">All clusters</option></select></label>
<a id="csv-link">Download cluster assignments</a>
<a id="png-link">Open static cluster figure</a>
</div>
<div id="cluster-summary"></div>
<div id="matrix-section" class="canvas-wrap" hidden>
<canvas id="matrix" width="700" height="700"></canvas>
<p id="matrix-hover" class="small">Hover over the matrix to inspect an ordered ligand pair.</p>
</div>
<table>
<thead><tr><th>Ligand ID</th><th>Cluster</th><th>Hierarchical order</th></tr></thead>
<tbody id="rows"></tbody>
</table>
<script>
const payload = __PAYLOAD__;
const search = document.getElementById("search");
const cluster = document.getElementById("cluster");
const tableBody = document.getElementById("rows");
const summary = document.getElementById("cluster-summary");
document.getElementById("csv-link").href = payload.assignments_csv;
document.getElementById("png-link").href = payload.matrix_png || payload.assignments_csv;

const clusterIds = [...new Set(payload.rows.map(row => row.cluster_id))].sort((a, b) => a - b);
for (const id of clusterIds) {
  const option = document.createElement("option");
  option.value = String(id);
  option.textContent = "Cluster " + id;
  cluster.appendChild(option);
  const count = payload.rows.filter(row => row.cluster_id === id).length;
  const chip = document.createElement("button");
  chip.className = "cluster-chip";
  chip.textContent = "Cluster " + id + ": " + count;
  chip.addEventListener("click", () => { cluster.value = String(id); renderRows(); });
  summary.appendChild(chip);
}

function renderRows() {
  const needle = search.value.trim().toLowerCase();
  const selected = cluster.value;
  tableBody.replaceChildren();
  for (const row of payload.rows) {
    if (selected && String(row.cluster_id) !== selected) continue;
    if (needle && !row.ligand_id.toLowerCase().includes(needle)) continue;
    const tr = document.createElement("tr");
    for (const value of [row.ligand_id, row.cluster_id, row.leaf_order]) {
      const td = document.createElement("td");
      td.textContent = String(value);
      tr.appendChild(td);
    }
    tableBody.appendChild(tr);
  }
}
search.addEventListener("input", renderRows);
cluster.addEventListener("change", renderRows);
renderRows();

if (payload.matrix_available && Array.isArray(payload.matrix)) {
  const section = document.getElementById("matrix-section");
  const canvas = document.getElementById("matrix");
  const hover = document.getElementById("matrix-hover");
  const ctx = canvas.getContext("2d");
  const n = payload.rows.length;
  const pixel = canvas.width / n;
  for (let row = 0; row < n; row += 1) {
    for (let column = 0; column < n; column += 1) {
      const value = Math.max(0, Math.min(1, Number(payload.matrix[row][column])));
      const red = Math.round(255 * Math.pow(value, 0.7));
      const green = Math.round(90 * value);
      const blue = Math.round(120 * (1 - value));
      ctx.fillStyle = "rgb(" + red + "," + green + "," + blue + ")";
      ctx.fillRect(column * pixel, row * pixel, Math.ceil(pixel), Math.ceil(pixel));
    }
  }
  canvas.addEventListener("mousemove", event => {
    const rect = canvas.getBoundingClientRect();
    const column = Math.min(n - 1, Math.max(0, Math.floor((event.clientX - rect.left) * n / rect.width)));
    const row = Math.min(n - 1, Math.max(0, Math.floor((event.clientY - rect.top) * n / rect.height)));
    const value = Number(payload.matrix[row][column]).toFixed(5);
    hover.textContent = payload.rows[row].ligand_id + " vs " + payload.rows[column].ligand_id + ": similarity " + value;
  });
  section.hidden = false;
}
</script>
</body>
</html>
"""
    html = (
        html.replace("__TITLE__", title)
        .replace("__MATRIX_NOTE__", escape(matrix_note))
        .replace("__PAYLOAD__", payload_text)
    )
    output = output_dir / f"clusters_{IFP_SUFFIXES.get(ifp_type, ifp_type)}.html"
    output.write_text(html, encoding="utf-8")
    return output


def _run_requested_fp_session(
    py_exe: str,
    workdir: Path,
    output_dir: Path,
    settings: Mapping[str, Any],
) -> tuple[Path | None, str | None]:
    request = settings.get("fp_session")
    if not isinstance(request, Mapping):
        return None, None

    ifp_type = str(request.get("ifp_type") or "").strip().upper()
    entry_name = str(request.get("entry_name") or "").strip()
    feature_id = request.get("feature_id")
    if not ifp_type or not entry_name or feature_id in (None, ""):
        return None, "fp_session requires ifp_type, entry_name, and feature_id."

    try:
        feature_id = int(feature_id)
    except (TypeError, ValueError):
        return None, "fp_session.feature_id must be an integer."

    requested_path = str(request.get("output_path") or "").strip()
    output = Path(requested_path) if requested_path else output_dir / (
        f"fp_session_{_safe_token(ifp_type)}_{_safe_token(entry_name)}_{feature_id}.pse"
    )
    if not output.is_absolute():
        output = workdir / output
    output.parent.mkdir(parents=True, exist_ok=True)
    response = analysis_runtime.generate_fp_session(
        py_exe,
        str(workdir),
        ifp_type,
        entry_name,
        feature_id,
        str(output),
    )
    if not isinstance(response, dict) or response.get("error"):
        message = response.get("error") if isinstance(response, dict) else "Unknown FP-PyMOL session error."
        return None, str(message)
    return output, None


def run_terminal_results(
    cfg: ProjectConfig,
    py_exe: str,
    settings: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Export all cached Results artifacts for a completed workdir."""
    settings = settings or {}
    workdir = Path(cfg.workdir)
    if not workdir.exists():
        raise FileNotFoundError(f"Workdir does not exist: {workdir}")

    output_dir = workdir / "results" / "terminal"
    output_dir.mkdir(parents=True, exist_ok=True)
    manifest: dict[str, Any] = {
        "schema_version": 1,
        "created_at": datetime.now(timezone.utc).isoformat(),
        "workdir": str(workdir),
        "outputs": {},
        "warnings": [],
        "errors": [],
        "gui_dashboard": "Fingerprint dashboard tables and high-resolution figures are exported when cached IFP artifacts are available.",
    }
    language = str(getattr(cfg, "language", "pt") or "pt")

    summary = analysis_runtime.run_analysis(py_exe, str(workdir))
    if isinstance(summary, dict) and not summary.get("error"):
        summary_path = workdir / "results" / "analysis_summary.json"
        if not summary_path.exists():
            _write_json(summary_path, summary)
        manifest["outputs"]["analysis_summary"] = _relative_to_workdir(summary_path, workdir)
        chart = _save_interaction_summary(summary, output_dir, language)
        if chart:
            manifest["outputs"]["interaction_summary_chart"] = _relative_to_workdir(chart, workdir)
    else:
        manifest["errors"].append(
            str(summary.get("error") if isinstance(summary, dict) else "Unable to load interaction summary.")
        )

    residue = analysis_runtime.run_residue_matrix(py_exe, str(workdir))
    if isinstance(residue, dict) and not residue.get("error"):
        residue_path = workdir / "results" / "residue_matrix.json"
        if not residue_path.exists():
            _write_json(residue_path, residue)
        manifest["outputs"]["residue_matrix"] = _relative_to_workdir(residue_path, workdir)
        heatmaps = _save_residue_heatmaps(residue, output_dir, language)
        if heatmaps:
            manifest["outputs"]["residue_heatmaps"] = [
                _relative_to_workdir(path, workdir) for path in heatmaps
            ]
        complete_heatmap = _save_complete_residue_heatmap(residue, output_dir, language)
        if complete_heatmap:
            manifest["outputs"]["complete_ligands_residues_heatmap"] = _relative_to_workdir(
                complete_heatmap, workdir
            )
    else:
        manifest["errors"].append(
            str(residue.get("error") if isinstance(residue, dict) else "Unable to load residue matrix.")
        )

    fp_artifacts = results_analysis.load_fp_analysis_artifacts(workdir)
    requested_types = set(cfg.selected_ifp_types())
    for ifp_type, artifact in sorted(fp_artifacts.items()):
        if requested_types and ifp_type not in requested_types:
            continue
        try:
            dashboard = results_analysis.build_fp_analysis_dashboard(
                workdir,
                artifact,
                labels_csv=cfg.fp_labels_csv,
                labels_id_column=cfg.fp_labels_id_column,
                labels_column=cfg.fp_labels_column,
                task_kind_preference=cfg.fp_label_task,
                use_otsu_threshold=cfg.fp_use_otsu_threshold,
            )
            dashboard_path = output_dir / "fingerprints" / IFP_SUFFIXES.get(ifp_type, ifp_type) / "fp_dashboard.json"
            _write_json(dashboard_path, dashboard)
            figures = _save_fp_dashboard_figures(
                dashboard,
                output_dir,
                ifp_type,
                str(getattr(cfg, "language", "pt") or "pt"),
            )
            manifest["outputs"][f"{ifp_type}_fp_dashboard"] = _relative_to_workdir(dashboard_path, workdir)
            manifest["outputs"][f"{ifp_type}_fp_charts"] = [
                _relative_to_workdir(path, workdir) for path in figures
            ]
        except Exception as exc:
            manifest["errors"].append(f"{ifp_type}: fingerprint dashboard export failed: {exc}")

    matrix_limit = _int_setting(settings, "terminal_matrix_max_entries", 5000)
    cluster_limit = _int_setting(settings, "terminal_cluster_max_entries", 5000)
    interactive_limit = _int_setting(settings, "terminal_interactive_max_entries", 2500)
    cluster_count = _int_setting(settings, "terminal_cluster_count", 4, minimum=2)
    cluster_method = str(settings.get("terminal_cluster_method") or "average").strip().lower()

    for ifp_type in cfg.selected_ifp_types():
        paths = _similarity_candidates(workdir, cfg, ifp_type)
        matrix_path = next((path for path in paths if path.exists()), None)
        if matrix_path is None:
            manifest["warnings"].append(f"{ifp_type}: no similarity matrix was found.")
            continue

        size_hint = _square_matrix_size_hint(matrix_path)
        if size_hint is None:
            manifest["warnings"].append(
                f"{ifp_type}: skipped terminal matrix loading because no square matrix header was available."
            )
            continue
        if size_hint > matrix_limit:
            manifest["warnings"].append(
                f"{ifp_type}: matrix with {size_hint} entries exceeds terminal_matrix_max_entries={matrix_limit}."
            )
            continue

        try:
            labels, matrix = results_analysis.load_similarity_matrix(matrix_path)
        except Exception as exc:
            manifest["errors"].append(f"{ifp_type}: could not load similarity matrix: {exc}")
            continue

        manifest["outputs"][f"{ifp_type}_similarity_matrix"] = _relative_to_workdir(matrix_path, workdir)
        chart = _save_similarity_figure(labels, matrix, output_dir, ifp_type, language)
        if chart:
            manifest["outputs"][f"{ifp_type}_similarity_chart"] = _relative_to_workdir(chart, workdir)

        if len(labels) < 2:
            manifest["warnings"].append(f"{ifp_type}: at least two ligands are required for clustering.")
            continue
        if len(labels) > cluster_limit:
            manifest["warnings"].append(
                f"{ifp_type}: clustering skipped for {len(labels)} entries; "
                f"terminal_cluster_max_entries={cluster_limit}."
            )
            continue

        try:
            result = results_analysis.cluster_similarity_matrix(
                labels,
                matrix,
                method=cluster_method,
                n_clusters=cluster_count,
            )
            suffix = IFP_SUFFIXES.get(ifp_type, ifp_type)
            assignments = output_dir / f"clusters_{suffix}.csv"
            results_analysis.export_cluster_assignments(assignments, result)
            cluster_chart = _save_cluster_figure(result, output_dir, ifp_type, language)
            explorer = _write_cluster_explorer(
                result,
                output_dir,
                ifp_type,
                include_matrix=len(labels) <= interactive_limit,
                matrix_png=cluster_chart,
                assignments_csv=assignments,
            )
            manifest["outputs"][f"{ifp_type}_cluster_assignments"] = _relative_to_workdir(assignments, workdir)
            manifest["outputs"][f"{ifp_type}_cluster_explorer"] = _relative_to_workdir(explorer, workdir)
            if cluster_chart:
                manifest["outputs"][f"{ifp_type}_cluster_chart"] = _relative_to_workdir(cluster_chart, workdir)
            if len(labels) > interactive_limit:
                manifest["warnings"].append(
                    f"{ifp_type}: the interactive explorer keeps every cluster assignment, "
                    f"but its matrix hover view is disabled above {interactive_limit} entries."
                )
        except Exception as exc:
            manifest["errors"].append(f"{ifp_type}: clustering failed: {exc}")

    session, session_error = _run_requested_fp_session(py_exe, workdir, output_dir, settings)
    if session:
        manifest["outputs"]["fp_pymol_session"] = _relative_to_workdir(session, workdir)
    if session_error:
        manifest["errors"].append(session_error)

    manifest_path = output_dir / "terminal_results_manifest.json"
    _write_json(manifest_path, manifest)
    manifest["outputs"]["manifest"] = _relative_to_workdir(manifest_path, workdir)
    _write_json(manifest_path, manifest)
    return manifest
