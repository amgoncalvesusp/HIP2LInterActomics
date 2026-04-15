"""Helpers for loading, clustering, and exporting result matrices."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import csv

import numpy as np
from scipy.cluster.hierarchy import fcluster, leaves_list, linkage
from scipy.spatial.distance import squareform


@dataclass
class ClusterResult:
    labels: list[str]
    matrix: np.ndarray
    linkage_matrix: np.ndarray
    cluster_ids: list[int]
    leaves: list[int]
    ordered_labels: list[str]
    ordered_cluster_ids: list[int]
    ordered_matrix: np.ndarray
    method: str
    requested_clusters: int

    @property
    def n_clusters(self) -> int:
        return len(set(self.cluster_ids))


def load_similarity_matrix(path: str | Path) -> tuple[list[str], np.ndarray]:
    """Read a square similarity matrix CSV and return labels plus a dense matrix."""
    path = Path(path)
    with path.open("r", encoding="utf-8", errors="replace", newline="") as fh:
        rows = [row for row in csv.reader(fh) if any(cell.strip() for cell in row)]

    if not rows:
        raise ValueError("A matriz de similaridade está vazia.")

    first = rows[0]
    has_header = len(first) > 1 and any(not _is_float(cell) for cell in first[1:])
    header_labels = [cell.strip() for cell in first[1:]] if has_header else []

    labels: list[str] = []
    data: list[list[float]] = []
    start_row = 1 if has_header else 0

    for row in rows[start_row:]:
        if not row:
            continue
        if not _is_float(row[0]):
            label = row[0].strip() or f"Ligante {len(labels) + 1}"
            labels.append(label)
            values = row[1:]
        else:
            values = row
        data.append([_safe_float(value) for value in values])

    if not data:
        raise ValueError("A matriz de similaridade não possui dados.")

    widths = {len(row) for row in data}
    if len(widths) != 1:
        raise ValueError("A matriz de similaridade possui linhas com tamanhos diferentes.")

    matrix = np.asarray(data, dtype=float)
    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
        raise ValueError(
            f"A matriz de similaridade deve ser quadrada. Recebido: {matrix.shape}."
        )

    if labels and len(labels) != matrix.shape[0]:
        raise ValueError("O número de rótulos não corresponde ao tamanho da matriz.")

    if not labels:
        if header_labels and len(header_labels) == matrix.shape[0]:
            labels = header_labels
        else:
            labels = [f"Ligante {i + 1}" for i in range(matrix.shape[0])]

    matrix = np.nan_to_num(matrix, nan=0.0, posinf=1.0, neginf=0.0)
    matrix = np.clip((matrix + matrix.T) / 2.0, 0.0, 1.0)
    np.fill_diagonal(matrix, 1.0)
    return labels, matrix


def cluster_similarity_matrix(
    labels: list[str],
    matrix: np.ndarray,
    method: str = "average",
    n_clusters: int = 4,
) -> ClusterResult:
    """Cluster a similarity matrix and return an ordering plus cluster assignments."""
    method = method.lower().strip()
    if method not in {"average", "complete", "single"}:
        raise ValueError(f"Método de clusterização não suportado: {method}")

    matrix = np.asarray(matrix, dtype=float)
    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
        raise ValueError("A matriz deve ser quadrada.")
    if matrix.shape[0] < 2:
        raise ValueError("São necessários pelo menos dois ligantes para clusterizar.")
    if len(labels) != matrix.shape[0]:
        raise ValueError("O número de rótulos não corresponde ao tamanho da matriz.")

    distance = 1.0 - np.clip((matrix + matrix.T) / 2.0, 0.0, 1.0)
    np.fill_diagonal(distance, 0.0)
    condensed = squareform(distance, checks=False)

    linkage_matrix = linkage(condensed, method=method)
    requested = max(2, min(int(n_clusters), matrix.shape[0]))
    cluster_ids = fcluster(linkage_matrix, t=requested, criterion="maxclust")
    leaves = leaves_list(linkage_matrix).tolist()

    ordered_labels = [labels[i] for i in leaves]
    ordered_cluster_ids = [int(cluster_ids[i]) for i in leaves]
    ordered_matrix = matrix[np.ix_(leaves, leaves)]

    return ClusterResult(
        labels=list(labels),
        matrix=matrix,
        linkage_matrix=linkage_matrix,
        cluster_ids=[int(x) for x in cluster_ids],
        leaves=list(leaves),
        ordered_labels=ordered_labels,
        ordered_cluster_ids=ordered_cluster_ids,
        ordered_matrix=ordered_matrix,
        method=method,
        requested_clusters=requested,
    )


def export_cluster_assignments(path: str | Path, result: ClusterResult) -> Path:
    """Save cluster assignments as CSV."""
    path = Path(path)
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["ligand_id", "cluster_id", "leaf_order"])
        for leaf_pos, item_index in enumerate(result.leaves, start=1):
            writer.writerow(
                [
                    result.labels[item_index],
                    result.cluster_ids[item_index],
                    leaf_pos,
                ]
            )
    return path


def cluster_rows(result: ClusterResult) -> list[tuple[str, int, int]]:
    """Return rows ready for tabular display/export."""
    rows: list[tuple[str, int, int]] = []
    for leaf_pos, item_index in enumerate(result.leaves, start=1):
        rows.append(
            (
                result.labels[item_index],
                result.cluster_ids[item_index],
                leaf_pos,
            )
        )
    return rows


def _is_float(value: str) -> bool:
    try:
        float(value)
        return True
    except Exception:
        return False


def _safe_float(value: str) -> float:
    try:
        return float(value)
    except Exception:
        return 0.0
