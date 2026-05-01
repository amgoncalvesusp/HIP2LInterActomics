"""Helpers for loading cached result artifacts, matrices, fingerprints, and clusters."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import csv
import json
import math
import re

import numpy as np
from scipy.cluster.hierarchy import fcluster, leaves_list, linkage
from scipy.spatial.distance import squareform

IFP_SUFFIX_TO_TYPE = {"E": "EIFP", "H": "HIFP", "F": "FIFP"}

CLASS_L0_LIGAND = "Ligand's level 0 features only"
CLASS_L0_PROTEIN = "Protein's level 0 features only"
CLASS_UPPER_LIGAND = "Upper level with ligand atomic information only"
CLASS_UPPER_PROTEIN = "Upper level with protein atomic information only"
CLASS_INTRALIGAND = "Intraligand interactions only"
CLASS_INTRAPROTEIN = "Intraprotein interactions only"
CLASS_NONCOVALENT = "Has noncovalent interactions with the protein"
CLASS_COLLISION = "Features with collision in the same complex"
CLASS_UNRELIABLE = "Unreliable feature"
EULER_MASCHERONI = 0.577215665
PREVALENCE_EXCLUDED_INTERACTIONS = {
    "weak hydrogen bond",
    "hydrophobic",
}

INTERACTION_COLORS = {
    "Hydrogen bond": "#05FF23",
    "Weak hydrogen bond": "#C5FFCC",
    "Halogen bond": "#00D6FA",
    "Halogen-pi": "#36B6C4",
    "Chalcogen bond": "#E7FF05",
    "Chalcogen-pi": "#B5C800",
    "Ionic": "#F69E00",
    "Salt bridge": "#FF9C25",
    "Cation-pi": "#AC5E00",
    "Cation-nucleophile": "#00AC3D",
    "Anion-electrophile": "#FFB7B7",
    "Anion-pi": "#FF5B5B",
    "Pi-stacking": "#6B6ECF",
    "Aromatic stacking": "#EA00C9",
    "Edge-to-face": "#EA00C9",
    "Face-to-face": "#EA00C9",
    "Face-to-edge pi-stacking": "#EA00C9",
    "Face-to-face pi-stacking": "#EA00C9",
    "Face-to-slope pi-stacking": "#EA00C9",
    "Displaced face-to-edge pi-stacking": "#EA00C9",
    "Displaced face-to-face pi-stacking": "#EA00C9",
    "Displaced face-to-slope pi-stacking": "#EA00C9",
    "Parallel": "#D5D5D5",
    "Parallel multipolar": "#A9A9A9",
    "Antiparallel multipolar": "#4E6DBA",
    "Orthogonal multipolar": "#324982",
    "Tilted multipolar": "#273863",
    "T-shaped": "#EA00C9",
    "Hydrophobic": "#FDD595",
    "Amide-aromatic stacking": "#BCBD22",
    "Charge-dipole interaction": "#9FBF00",
    "Unfavorable charge-dipole interaction": "#B3DE69",
    "Unfavorable dipole interaction": "#CCEBC5",
    "Water-bridged hydrogen bond": "#0701FF",
    "Disulfide bond": "#A8A858",
    "Metal coordination": "#5C01FF",
    "Van der Waals": "#808080",
    "Proximal": "#808080",
    "Multipolar interaction": "#8E63CE",
    "Multiple interactions": "#B2ABD2",
    "Repulsive": "#FF0101",
    "Unfavorable anion-nucleophile": "#C44E52",
    "Unfavorable cation-electrophile": "#DD8452",
    "Unfavorable electrophile-electrophile": "#937860",
    "Unfavorable nucleophile-nucleophile": "#8C8C8C",
}

INTERACTION_PRIORITY = [
    "Ionic",
    "Cation-nucleophile",
    "Cation-pi",
    "Anion-electrophile",
    "Anion-pi",
    "Hydrogen bond",
    "Halogen bond",
    "Chalcogen bond",
    "Chalcogen-pi",
    "Face-to-face pi-stacking",
    "Face-to-slope pi-stacking",
    "Face-to-edge pi-stacking",
    "Displaced face-to-face pi-stacking",
    "Displaced face-to-slope pi-stacking",
    "Displaced face-to-edge pi-stacking",
    "Weak hydrogen bond",
    "Hydrophobic",
    "Repulsive",
    "Unfavorable",
]

INTERACTION_PRIORITY_ALIASES = {
    "salt bridge": "Ionic",
    "face-to-face": "Face-to-face pi-stacking",
    "face-to-slope": "Face-to-slope pi-stacking",
    "face-to-edge": "Face-to-edge pi-stacking",
    "displaced face-to-face": "Displaced face-to-face pi-stacking",
    "displaced face-to-slope": "Displaced face-to-slope pi-stacking",
    "displaced face-to-edge": "Displaced face-to-edge pi-stacking",
}

_ENTRY_ID_SUFFIXES = (
    "_ligand",
    "-ligand",
    "_lig",
    "-lig",
    "_LIGAND",
    "-LIGAND",
    "_LIG",
    "-LIG",
)

FP_CLASS_ORDER = [
    CLASS_NONCOVALENT,
    CLASS_L0_LIGAND,
    CLASS_L0_PROTEIN,
    CLASS_UPPER_LIGAND,
    CLASS_UPPER_PROTEIN,
    CLASS_INTRALIGAND,
    CLASS_INTRAPROTEIN,
    CLASS_COLLISION,
    CLASS_UNRELIABLE,
]

FP_CLASS_COLORS = {
    CLASS_NONCOVALENT: "#f28e63",
    CLASS_L0_LIGAND: "#6f9ec7",
    CLASS_L0_PROTEIN: "#c39ac9",
    CLASS_UPPER_LIGAND: "#a24bd6",
    CLASS_UPPER_PROTEIN: "#f2cf44",
    CLASS_INTRALIGAND: "#2f9c61",
    CLASS_INTRAPROTEIN: "#a6c95f",
    CLASS_COLLISION: "#c8c8c8",
    CLASS_UNRELIABLE: "#7d7d7d",
}

FP_CLASS_ALIASES = {
    "ligand level 0 features only": CLASS_L0_LIGAND,
    "ligand's level 0 features": CLASS_L0_LIGAND,
    "ligand’s level 0 features only": CLASS_L0_LIGAND,
    "ligand’s level 0 features": CLASS_L0_LIGAND,
    "protein level 0 features only": CLASS_L0_PROTEIN,
    "protein's level 0 features": CLASS_L0_PROTEIN,
    "protein’s level 0 features only": CLASS_L0_PROTEIN,
    "protein’s level 0 features": CLASS_L0_PROTEIN,
    "upper level with ligand atomic information only": CLASS_UPPER_LIGAND,
    "upper level with ligand atomic information": CLASS_UPPER_LIGAND,
    "upper level with protein atomic information only": CLASS_UPPER_PROTEIN,
    "upper level with protein atomic information": CLASS_UPPER_PROTEIN,
    "intraligand interactions only": CLASS_INTRALIGAND,
    "only intra-ligand interactions": CLASS_INTRALIGAND,
    "intraprotein interactions only": CLASS_INTRAPROTEIN,
    "only intra-protein interactions": CLASS_INTRAPROTEIN,
    "has noncovalent interactions with the protein": CLASS_NONCOVALENT,
    "has non-covalent interactions with the protein": CLASS_NONCOVALENT,
    "features with collision in the same complex": CLASS_COLLISION,
    "with collision": CLASS_COLLISION,
    "unreliable feature": CLASS_UNRELIABLE,
}

FP_CLASS_RANK = {label: idx for idx, label in enumerate(FP_CLASS_ORDER)}

AA_THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "HID": "H", "HIE": "H", "HIP": "H", "HSD": "H", "HSE": "H", "HSP": "H",
    "CYM": "C", "CYX": "C", "MSE": "M", "SEC": "U",
}


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


def load_analysis_summary(workdir: str | Path) -> dict | None:
    """Load the cached interaction-summary JSON from a workdir, if available."""
    return _load_cached_json(workdir, ["results/analysis_summary.json", "analysis_summary.json"])


def load_residue_matrix_artifact(workdir: str | Path) -> dict | None:
    """Load the cached residue-matrix JSON from a workdir, if available."""
    return _load_cached_json(workdir, ["results/residue_matrix.json", "residue_matrix.json"])


def load_fp_analysis_artifacts(workdir: str | Path) -> dict[str, dict]:
    """Load cached fingerprint-analysis JSON artifacts keyed by IFP type."""
    fp_dir = Path(workdir) / "results" / "fingerprints"
    if not fp_dir.exists():
        return {}

    artifacts: dict[str, dict] = {}
    for candidate in sorted(fp_dir.glob("fp_analysis_*.json")):
        suffix = candidate.stem.rsplit("_", 1)[-1].upper()
        ifp_type = IFP_SUFFIX_TO_TYPE.get(suffix)
        if not ifp_type:
            continue
        with candidate.open("r", encoding="utf-8", errors="replace") as fh:
            artifacts[ifp_type] = json.load(fh)
    return artifacts


def load_fp_detail_artifact(workdir: str | Path, ifp_type: str | None) -> dict | None:
    if not ifp_type:
        return None
    suffix = {"EIFP": "E", "HIFP": "H", "FIFP": "F"}.get(str(ifp_type).upper())
    if suffix is None:
        return None
    path = Path(workdir) / "results" / "fingerprints" / f"fp_detail_{suffix}.json"
    if not path.exists():
        return None
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return None


def build_complete_heatmap(artifact: dict) -> tuple[list[str], list[str], np.ndarray, list[str]]:
    """Collapse the residue matrix into a categorical entry x residue map."""
    entries = list(artifact.get("entries", []) or [])
    residues = list(artifact.get("residues", []) or [])
    interaction_types = list(artifact.get("interaction_types", []) or [])

    if not entries or not residues or not interaction_types:
        return entries, residues, np.zeros((0, 0), dtype=int), interaction_types

    cube: list[np.ndarray] = []
    for interaction_type in interaction_types:
        values = artifact.get("matrix", {}).get(interaction_type) or []
        cube.append(np.asarray(values, dtype=float))

    if not cube:
        return entries, residues, np.zeros((0, 0), dtype=int), interaction_types

    stacked = np.stack(cube, axis=0)
    if stacked.ndim != 3:
        return entries, residues, np.zeros((0, 0), dtype=int), interaction_types

    any_signal = np.any(stacked > 0, axis=0)
    dominant = np.argmax(stacked, axis=0) + 1
    categorical = np.where(any_signal, dominant, 0)

    keep = np.any(categorical > 0, axis=0)
    residues = [residue for residue, include in zip(residues, keep) if include]
    categorical = categorical[:, keep]
    return entries, residues, categorical, interaction_types


def build_complete_heatmap_layers(
    artifact: dict,
) -> tuple[list[str], list[str], list[list[list[str]]], list[str]]:
    entries = list(artifact.get("entries", []) or [])
    residues = list(artifact.get("residues", []) or [])
    raw_interaction_types = [normalize_interaction_name(name) for name in (artifact.get("interaction_types", []) or [])]

    if not entries or not residues or not raw_interaction_types:
        return entries, residues, [], []

    present_matrix: dict[str, np.ndarray] = {}
    for raw_name, original_name in zip(raw_interaction_types, (artifact.get("interaction_types", []) or [])):
        values = artifact.get("matrix", {}).get(original_name) or []
        arr = np.asarray(values, dtype=float)
        if arr.ndim != 2 or arr.shape[0] != len(entries) or arr.shape[1] != len(residues):
            continue
        present_matrix[raw_name] = arr

    if not present_matrix:
        return entries, residues, [], []

    keep_mask = np.zeros(len(residues), dtype=bool)
    for arr in present_matrix.values():
        keep_mask |= np.any(arr > 0.0, axis=0)

    kept_residues = [residue for residue, keep in zip(residues, keep_mask) if keep]
    if not kept_residues:
        return entries, [], [[[] for _ in range(0)] for _ in entries], []

    ordered_interactions = sorted(present_matrix, key=interaction_priority_key)
    trimmed = {
        name: arr[:, keep_mask]
        for name, arr in present_matrix.items()
    }

    layered: list[list[list[str]]] = []
    present_types: set[str] = set()
    for row_idx in range(len(entries)):
        row: list[list[str]] = []
        for col_idx in range(len(kept_residues)):
            cell_types = [
                name
                for name in ordered_interactions
                if float(trimmed[name][row_idx, col_idx]) > 0.0
            ]
            for name in cell_types:
                present_types.add(name)
            row.append(cell_types)
        layered.append(row)

    ordered_present = [name for name in ordered_interactions if name in present_types]
    return entries, kept_residues, layered, ordered_present


def build_trajectory_frame_percentages(artifact: dict) -> tuple[list[str], list[str], np.ndarray]:
    """Return residue x interaction percentages across trajectory frames."""
    entries = list(artifact.get("entries", []) or [])
    residues = list(artifact.get("residues", []) or [])
    interaction_types = list(artifact.get("interaction_types", []) or [])
    if not entries or not residues or not interaction_types:
        return [], [], np.zeros((0, 0), dtype=float)

    n_frames = max(1, len(entries))
    columns: list[np.ndarray] = []
    kept_types: list[str] = []
    for interaction_type in interaction_types:
        raw = (artifact.get("matrix") or {}).get(interaction_type)
        if raw is None:
            continue
        arr = np.asarray(raw, dtype=float)
        if arr.ndim != 2 or arr.shape[0] != len(entries) or arr.shape[1] != len(residues):
            continue
        pct = 100.0 * np.count_nonzero(arr > 0.0, axis=0).astype(float) / float(n_frames)
        if np.any(pct > 0.0):
            columns.append(pct)
            kept_types.append(str(interaction_type))

    if not columns:
        return [], [], np.zeros((0, 0), dtype=float)
    matrix = np.vstack(columns).T
    keep_residue = np.sum(matrix, axis=1) > 0.0
    filtered_residues = [residue for residue, keep in zip(residues, keep_residue) if bool(keep)]
    return filtered_residues, kept_types, matrix[keep_residue, :]


def build_trajectory_entry_counts(
    artifact: dict,
    entry_name: str,
) -> tuple[list[str], list[str], np.ndarray]:
    """Return residue x interaction counts for one trajectory frame/entry."""
    entries = list(artifact.get("entries", []) or [])
    residues = list(artifact.get("residues", []) or [])
    interaction_types = list(artifact.get("interaction_types", []) or [])
    try:
        entry_index = entries.index(entry_name)
    except ValueError:
        return [], [], np.zeros((0, 0), dtype=float)

    columns: list[np.ndarray] = []
    kept_types: list[str] = []
    for interaction_type in interaction_types:
        raw = (artifact.get("matrix") or {}).get(interaction_type)
        if raw is None:
            continue
        arr = np.asarray(raw, dtype=float)
        if arr.ndim != 2 or arr.shape[0] <= entry_index or arr.shape[1] != len(residues):
            continue
        counts = np.asarray(arr[entry_index, :], dtype=float)
        if np.any(counts > 0.0):
            columns.append(counts)
            kept_types.append(str(interaction_type))

    if not columns:
        return [], [], np.zeros((0, 0), dtype=float)
    matrix = np.vstack(columns).T
    keep_residue = np.sum(matrix, axis=1) > 0.0
    filtered_residues = [residue for residue, keep in zip(residues, keep_residue) if bool(keep)]
    return filtered_residues, kept_types, matrix[keep_residue, :]


def format_residue_label(value: str | None) -> str:
    raw = str(value or "").strip()
    if not raw:
        return ""

    parts = raw.split("/")
    if len(parts) < 3:
        return raw

    chain = parts[0].strip() or "?"
    residue = parts[1].strip().upper()
    index = "/".join(parts[2:]).strip()
    if not index:
        return raw

    aa = AA_THREE_TO_ONE.get(residue[:3], residue[:3] or "?")
    match = re.match(r"^(-?\d+)([A-Za-z]?)$", index)
    if match:
        seq = match.group(1)
        insertion = match.group(2)
        return f"{chain}:{aa}{seq}{insertion}"
    return f"{chain}:{aa}{index}"


def trajectory_frame_number(entry_name: str) -> int | None:
    """Return the full trailing frame/pose number from an entry name."""
    text = str(entry_name)
    matches = re.findall(r"(?:frame|pose)?[_\-\s]*(\d+)(?=\D*$)", text, flags=re.IGNORECASE)
    if not matches:
        return None
    try:
        return int(matches[-1])
    except ValueError:
        return None


def normalize_interaction_name(name: str | None) -> str:
    raw = str(name or "").strip()
    if not raw:
        return ""
    alias = INTERACTION_PRIORITY_ALIASES.get(raw.lower())
    if alias:
        return alias
    if raw.lower().startswith("unfavorable"):
        return raw
    return raw


def interaction_priority_key(name: str | None) -> tuple[int, str]:
    normalized = normalize_interaction_name(name)
    lowered = normalized.lower()
    for idx, candidate in enumerate(INTERACTION_PRIORITY):
        if candidate == "Unfavorable":
            if lowered.startswith("unfavorable"):
                return idx, lowered
            continue
        if lowered == candidate.lower():
            return idx, lowered
    return len(INTERACTION_PRIORITY), lowered


def get_interaction_color(name: str | None) -> str:
    raw = str(name or "").strip()
    if not raw:
        return "#2f7f83"
    direct = INTERACTION_COLORS.get(raw)
    if direct:
        return direct
    for key, value in INTERACTION_COLORS.items():
        if key.lower() == raw.lower():
            return value
    normalized = normalize_interaction_name(raw)
    direct = INTERACTION_COLORS.get(normalized)
    if direct:
        return direct
    if normalized.lower().startswith("unfavorable"):
        return "#8c8c8c"
    for key, value in INTERACTION_COLORS.items():
        if key.lower() == normalized.lower():
            return value
    return "#2f7f83"


def load_similarity_matrix(path: str | Path) -> tuple[list[str], np.ndarray]:
    """Read a square similarity matrix CSV and return labels plus a dense matrix."""
    path = Path(path)
    with path.open("r", encoding="utf-8", errors="replace", newline="") as fh:
        rows = [row for row in csv.reader(fh) if any(cell.strip() for cell in row)]

    if not rows:
        raise ValueError("A matriz de similaridade está vazia.")

    first = rows[0]
    if _looks_like_similarity_edge_list(first):
        return _load_similarity_edge_list(rows[1:])

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


def _looks_like_similarity_edge_list(header: list[str]) -> bool:
    normalized = [cell.strip().lower() for cell in header[:3]]
    return normalized == ["entry1", "entry2", "similarity"]


def _load_similarity_edge_list(rows: list[list[str]]) -> tuple[list[str], np.ndarray]:
    labels: list[str] = []
    seen: set[str] = set()
    edges: dict[tuple[str, str], list[float]] = {}

    for row in rows:
        if len(row) < 3:
            continue
        entry1 = row[0].strip()
        entry2 = row[1].strip()
        if not entry1 or not entry2:
            continue
        for label in (entry1, entry2):
            if label not in seen:
                seen.add(label)
                labels.append(label)
        key = tuple(sorted((entry1, entry2)))
        edges.setdefault(key, []).append(_safe_float(row[2]))

    if not labels:
        raise ValueError("A matriz de similaridade n\u00e3o possui arestas.")

    size = len(labels)
    index = {label: pos for pos, label in enumerate(labels)}
    matrix = np.zeros((size, size), dtype=float)
    np.fill_diagonal(matrix, 1.0)

    for (entry1, entry2), values in edges.items():
        value = sum(values) / max(1, len(values))
        i = index[entry1]
        j = index[entry2]
        matrix[i, j] = value
        matrix[j, i] = value

    matrix = np.nan_to_num(matrix, nan=0.0, posinf=1.0, neginf=0.0)
    matrix = np.clip((matrix + matrix.T) / 2.0, 0.0, 1.0)
    np.fill_diagonal(matrix, 1.0)
    return labels, matrix


def _load_cached_json(workdir: str | Path, relatives: list[str]) -> dict | None:
    root = Path(workdir)
    for relative in relatives:
        candidate = root / relative
        if not candidate.exists():
            continue
        with candidate.open("r", encoding="utf-8", errors="replace") as fh:
            return json.load(fh)
    return None


def normalize_fp_class_name(name: str | None) -> str:
    raw = str(name or "").strip()
    if not raw:
        return CLASS_UNRELIABLE

    raw = raw.replace("’", "'")

    alias = FP_CLASS_ALIASES.get(raw.lower())
    if alias:
        return alias

    if raw in FP_CLASS_ORDER:
        return raw

    lowered = raw.lower()
    if "ligand-centered" in lowered or "protein-centered" in lowered or "unclassified" in lowered:
        return CLASS_UNRELIABLE
    return raw


def normalize_fp_breakdown(breakdown: dict | None) -> dict[str, int]:
    normalized: dict[str, int] = {}
    for label, value in (breakdown or {}).items():
        try:
            count = int(value)
        except Exception:
            count = 0
        if count <= 0:
            continue
        norm_label = normalize_fp_class_name(label)
        normalized[norm_label] = normalized.get(norm_label, 0) + count
    return normalized


def resolve_ifp_source_file(workdir: str | Path, artifact: dict) -> Path | None:
    source = str(artifact.get("source_file") or "").strip()
    candidates: list[Path] = []
    if source:
        candidates.append(Path(source))

    workdir = Path(workdir)
    fp_dir = workdir / "results" / "fingerprints"
    ifp_type = str(artifact.get("ifp_type") or "").upper()
    suffix = {"EIFP": "E", "HIFP": "H", "FIFP": "F"}.get(ifp_type)
    if suffix:
        candidates.append(fp_dir / f"ifp_{suffix}.csv")
    candidates.append(fp_dir / "ifp.csv")

    for candidate in candidates:
        if candidate.exists():
            return candidate
    return None


def resolve_fp_labels_file(workdir: str | Path, labels_csv: str | Path | None) -> Path | None:
    raw = str(labels_csv or "").strip()
    if not raw:
        return None

    candidate = Path(raw)
    if candidate.exists():
        return candidate

    workdir = Path(workdir)
    if not candidate.is_absolute():
        relative = workdir / candidate
        if relative.exists():
            return relative
    return None


def _normalize_label_key(value: str | None) -> str:
    return re.sub(r"\s+", "", str(value or "").strip()).lower()


def _strip_entry_suffix(value: str) -> str:
    current = str(value or "").strip()
    changed = True
    while changed and current:
        changed = False
        for suffix in _ENTRY_ID_SUFFIXES:
            if current.endswith(suffix):
                current = current[: -len(suffix)]
                changed = True
                break
    return current


def _candidate_label_ids(entry_name: str) -> list[str]:
    raw = str(entry_name or "").strip()
    if not raw:
        return []

    parts: list[str] = [raw]
    colon_tokens = [token.strip() for token in raw.split(":") if token.strip()] if ":" in raw else []
    slash_tokens = [token.strip() for token in raw.split("/") if token.strip()] if "/" in raw else []

    # Prioritize the ligand-side token from names such as PROTEIN:ID_LIG so that
    # external labels keyed only by ID match the ligand, not the protein prefix.
    if colon_tokens:
        parts.append(colon_tokens[-1])
    elif slash_tokens:
        parts.append(slash_tokens[-1])

    parts.extend(colon_tokens)
    parts.extend(slash_tokens)

    expanded: list[str] = []
    seen: set[str] = set()
    for part in parts:
        candidate = str(part).strip()
        while candidate:
            key = _normalize_label_key(candidate)
            if key and key not in seen:
                seen.add(key)
                expanded.append(candidate)
            stripped = _strip_entry_suffix(candidate)
            if stripped == candidate:
                break
            candidate = stripped
    return expanded


def _resolve_external_label_value(label_map: dict[str, str], entry_name: str) -> str:
    normalized = {_normalize_label_key(key): value for key, value in label_map.items()}
    for candidate in _candidate_label_ids(entry_name):
        value = normalized.get(_normalize_label_key(candidate))
        if value not in (None, ""):
            return str(value)
    return ""


def _safe_float_or_none(value: object) -> float | None:
    text = str(value or "").strip()
    if not text:
        return None
    try:
        return float(text.replace(",", "."))
    except Exception:
        return None


def _dominant_label(counts: dict[str, int]) -> tuple[str, int]:
    if not counts:
        return "", 0
    label = max(sorted(counts), key=lambda key: (counts[key], key))
    return label, int(counts[label])


def _prevalence_interaction_name(value: str) -> str:
    text = str(value or "").strip()
    if "||" in text:
        text = text.split("||", 1)[0].strip()
    return text


def _prevalence_residue_from_pair(value: str) -> str:
    text = str(value or "").strip()
    if "||" not in text:
        return ""
    return text.split("||", 1)[1].strip()


def _is_prevalence_interaction_allowed(value: str) -> bool:
    name = re.sub(r"\s+", " ", _prevalence_interaction_name(value)).casefold()
    return bool(name) and name not in PREVALENCE_EXCLUDED_INTERACTIONS


def _filtered_prevalence_counts(feature: dict) -> tuple[dict[str, int], dict[str, int], dict[str, int], dict[str, dict]]:
    interaction_breakdown = {
        str(key): int(value)
        for key, value in (feature.get("interaction_breakdown") or {}).items()
        if int(value) > 0 and _is_prevalence_interaction_allowed(str(key))
    }
    pair_breakdown = {
        str(key): int(value)
        for key, value in (feature.get("pair_breakdown") or {}).items()
        if int(value) > 0 and _is_prevalence_interaction_allowed(str(key))
    }
    residue_breakdown: dict[str, int] = {}
    for pair_name, count in pair_breakdown.items():
        residue_name = _prevalence_residue_from_pair(pair_name)
        if residue_name:
            residue_breakdown[residue_name] = residue_breakdown.get(residue_name, 0) + int(count)

    entry_details: dict[str, dict] = {}
    for entry_name, entry_info in (feature.get("entry_details") or {}).items():
        filtered_interactions = {
            str(key): int(value)
            for key, value in ((entry_info or {}).get("interaction_counts") or {}).items()
            if int(value) > 0 and _is_prevalence_interaction_allowed(str(key))
        }
        filtered_pairs = {
            str(key): int(value)
            for key, value in ((entry_info or {}).get("pair_counts") or {}).items()
            if int(value) > 0 and _is_prevalence_interaction_allowed(str(key))
        }
        filtered_residues: dict[str, int] = {}
        for pair_name, count in filtered_pairs.items():
            residue_name = _prevalence_residue_from_pair(pair_name)
            if residue_name:
                filtered_residues[residue_name] = filtered_residues.get(residue_name, 0) + int(count)
        entry_details[str(entry_name)] = {
            "shell_count": int((entry_info or {}).get("shell_count", 0) or 0),
            "interaction_counts": filtered_interactions,
            "residue_counts": filtered_residues,
            "pair_counts": filtered_pairs,
        }

    return interaction_breakdown, residue_breakdown, pair_breakdown, entry_details


def _zscore_threshold(values: list[float]) -> tuple[float, str, list[float]]:
    arr = np.asarray([float(value) for value in values if float(value) > 0.0], dtype=float)
    if arr.size == 0:
        return 100.0, "no_positive_values", []
    mean_value = float(arr.mean())
    std_value = float(arr.std())
    zscores: list[float] = []
    for value in values:
        raw = float(value)
        zscores.append((raw - mean_value) / std_value if raw > 0.0 and std_value > 1e-12 else 0.0)
    threshold_candidates = [float(value) for value, zscore in zip(values, zscores) if zscore > 1.0]
    if threshold_candidates:
        return float(min(threshold_candidates)), "zscore_gt_1", zscores
    return 100.0, "no_reference_single_class_only", zscores


def _otsu_threshold(values: list[float]) -> tuple[float, str]:
    """Return an Otsu threshold for positive percentages.

    Otsu is used only as a fallback when the z-score rule cannot define a
    reference set. Values at or above the returned threshold are retained.
    """
    arr = np.asarray([float(value) for value in values if float(value) > 0.0], dtype=float)
    if arr.size == 0:
        return 100.0, "no_positive_values"
    unique = np.unique(arr)
    if unique.size == 1:
        return float(unique[0]), "otsu_single_value"

    candidates = (unique[:-1] + unique[1:]) / 2.0
    best_threshold = float(candidates[0])
    best_score = -1.0
    for threshold in candidates:
        low = arr[arr <= threshold]
        high = arr[arr > threshold]
        if low.size == 0 or high.size == 0:
            continue
        weight_low = float(low.size) / float(arr.size)
        weight_high = float(high.size) / float(arr.size)
        score = weight_low * weight_high * (float(low.mean()) - float(high.mean())) ** 2
        if score > best_score:
            best_score = score
            best_threshold = float(threshold)
    return best_threshold, "otsu"


def _threshold_with_otsu_fallback(
    values: list[float],
    use_otsu_threshold: bool = False,
) -> tuple[float, str, list[float]]:
    threshold, source, zscores = _zscore_threshold(values)
    if source == "zscore_gt_1" or not use_otsu_threshold:
        return threshold, source, zscores
    otsu_threshold, otsu_source = _otsu_threshold(values)
    return otsu_threshold, otsu_source, zscores


def _gumbel_tail_p_value(zscore: float) -> float:
    z = float(zscore)
    value = 1.0 - math.exp(-math.exp(((-z * math.pi) / math.sqrt(6.0)) - EULER_MASCHERONI))
    return float(min(1.0, max(0.0, value)))


def _importance_p_value(zscore: float) -> float:
    return _gumbel_tail_p_value(zscore)


def load_external_fp_labels(
    path: str | Path,
    label_column: str | None = None,
    id_column: str | None = None,
) -> dict[str, str]:
    path = Path(path)
    with path.open("r", encoding="utf-8", errors="replace", newline="") as fh:
        sample = fh.read(4096)
        fh.seek(0)
        try:
            dialect = csv.Sniffer().sniff(sample, delimiters=",\t;")
        except Exception:
            dialect = csv.excel_tab if path.suffix.lower() == ".tsv" else csv.excel
        reader = csv.reader(fh, dialect=dialect)
        rows = [row for row in reader if any(str(cell).strip() for cell in row)]

    if not rows:
        raise ValueError("O CSV de rótulos está vazio.")

    header = [str(cell).strip() for cell in rows[0]]
    id_col = 0
    label_col = 1 if len(header) > 1 else 0
    has_header = False

    lowered = [cell.lower() for cell in header]
    id_candidates = {
        "ligand_id",
        "entry",
        "entry_id",
        "ligand",
        "molecule",
        "mol",
        "name",
        "id",
        "molecule_chembl_id",
        "chembl_id",
        "compound_id",
        "molecule_id",
    }
    label_candidates = {
        "label",
        "rótulo",
        "rotulo",
        "classe",
        "class",
        "target",
        "activity",
        "group",
        "grupo",
        "category",
        "categoria",
        "y",
        "value",
        "valor",
    }
    requested_id_column = str(id_column or "").strip()
    requested_label_column = str(label_column or "").strip()
    if requested_id_column:
        requested = requested_id_column.strip().lower()
        for idx, name in enumerate(lowered):
            if name == requested:
                id_col = idx
                has_header = True
                break
        else:
            raise ValueError(
                f"A coluna de ID '{requested_id_column}' não foi encontrada no CSV."
            )
    else:
        for idx, name in enumerate(lowered):
            if name in id_candidates:
                id_col = idx
                has_header = True
                break
    if requested_label_column:
        requested = requested_label_column.strip().lower()
        for idx, name in enumerate(lowered):
            if name == requested:
                label_col = idx
                has_header = True
                break
        else:
            raise ValueError(
                f"A coluna de rótulo '{requested_label_column}' não foi encontrada no CSV."
            )
    else:
        for idx, name in enumerate(lowered):
            if name in label_candidates:
                label_col = idx
                has_header = True
                break

    if len(header) < 2 and not has_header:
        raise ValueError("O CSV de rótulos precisa ter pelo menos duas colunas: identificador e classe.")

    start_row = 1 if has_header else 0
    mapping: dict[str, str] = {}
    for row in rows[start_row:]:
        if len(row) <= max(id_col, label_col):
            continue
        ligand_id = str(row[id_col]).strip()
        label = str(row[label_col]).strip()
        if not ligand_id or not label:
            continue
        mapping[ligand_id] = label

    if not mapping:
        raise ValueError("Nenhum par ligand_id/label válido foi encontrado no CSV de rótulos.")
    return mapping


def load_ifp_sparse_matrix(path: str | Path) -> tuple[list[str], list[int], np.ndarray]:
    path = Path(path)
    labels: list[str] = []
    rows: list[dict[int, float]] = []
    feature_ids: set[int] = set()

    with path.open("r", encoding="utf-8", errors="replace", newline="") as fh:
        reader = csv.DictReader(fh)
        for record in reader:
            label = str(record.get("ligand_id") or "").strip()
            if not label:
                continue

            raw_bits = [token.strip() for token in str(record.get("on_bits") or "").split("\t") if token.strip()]
            raw_counts = [token.strip() for token in str(record.get("count") or "").split("\t") if token.strip()]
            bits = [int(token) for token in raw_bits]

            if raw_counts and len(raw_counts) != len(bits):
                raise ValueError(
                    f"A linha '{label}' possui {len(bits)} bits, mas {len(raw_counts)} contagens."
                )

            counts = [float(token) for token in raw_counts] if raw_counts else [1.0] * len(bits)
            row_map: dict[int, float] = {}
            for bit, count in zip(bits, counts):
                row_map[bit] = row_map.get(bit, 0.0) + count
            labels.append(label)
            rows.append(row_map)
            feature_ids.update(row_map)

    ordered_features = sorted(feature_ids)
    matrix = np.zeros((len(labels), len(ordered_features)), dtype=float)
    feature_index = {feature_id: idx for idx, feature_id in enumerate(ordered_features)}
    for row_index, row_map in enumerate(rows):
        for feature_id, count in row_map.items():
            matrix[row_index, feature_index[feature_id]] = count
    return labels, ordered_features, matrix


def build_fp_analysis_dashboard(
    workdir: str | Path,
    artifact: dict,
    labels_csv: str | Path | None = None,
    labels_id_column: str | None = None,
    labels_column: str | None = None,
    algorithm_preference: str = "gradient_boosting",
    task_kind_preference: str | None = None,
    use_otsu_threshold: bool = False,
) -> dict:
    source_file = resolve_ifp_source_file(workdir, artifact)
    if source_file is None:
        raise FileNotFoundError("Nenhum arquivo de fingerprint foi encontrado para esta analise.")

    labels, feature_ids, matrix = load_ifp_sparse_matrix(source_file)
    feature_lookup = {int(feature_id): col for col, feature_id in enumerate(feature_ids)}

    raw_features = artifact.get("features", []) or []
    enriched: list[dict] = []
    prevalence_values: list[float] = []

    for row in raw_features:
        feature_id = int(row.get("feature_id", 0))
        breakdown = normalize_fp_breakdown(row.get("nature_breakdown"))
        if not breakdown:
            breakdown = {CLASS_UNRELIABLE: 1}

        total_class_hits = sum(breakdown.values())
        top_class = max(
            sorted(breakdown),
            key=lambda key: (breakdown[key], -FP_CLASS_RANK.get(key, len(FP_CLASS_ORDER)), key),
        )
        top_pct = (100.0 * breakdown[top_class] / total_class_hits) if total_class_hits else 0.0
        class_percentages = {
            class_name: (100.0 * breakdown.get(class_name, 0) / total_class_hits) if total_class_hits else 0.0
            for class_name in FP_CLASS_ORDER
            if breakdown.get(class_name, 0) > 0
        }

        coverage_pct = float(row.get("coverage_pct", 0.0) or 0.0)
        feature_matrix_col = feature_lookup.get(feature_id)
        coverage_from_matrix = 0.0
        if feature_matrix_col is not None and len(labels) > 0:
            coverage_from_matrix = 100.0 * float(np.count_nonzero(matrix[:, feature_matrix_col] > 0)) / float(len(labels))
        prevalence_values.append(top_pct)
        enriched.append(
            {
                "feature_id": feature_id,
                "molecule_hits": int(row.get("molecule_hits", 0) or 0),
                "coverage_pct": coverage_pct if coverage_pct > 0 else coverage_from_matrix,
                "top_class": top_class,
                "top_class_pct": top_pct,
                "collision_hits": int(row.get("collision_hits", 0) or 0),
                "total_count": int(row.get("total_count", 0) or 0),
                "class_breakdown": breakdown,
                "class_percentages": class_percentages,
                "missing_molecules": max(0, len(labels) - int(row.get("molecule_hits", 0) or 0)),
                "zscore": 0.0,
                "assigned_class": CLASS_UNRELIABLE,
                "reliable": False,
                "importance_eligible": False,
                "importance_score": 0.0,
                "importance_pct": 0.0,
                "importance_rank": 0,
            }
        )

    if not enriched:
        return {
            "ifp_type": artifact.get("ifp_type"),
            "ifp_label": artifact.get("ifp_label"),
            "source_file": str(source_file),
            "features": [],
            "important_features": [],
            "class_share": {},
            "class_counts": {},
            "threshold_pct": 100.0,
            "model_name": "Unavailable",
            "model_note": "Sem features para analisar.",
            "cluster_count": 0,
            "class_order": list(FP_CLASS_ORDER),
        }

    prevalence_source = [
        feature for feature in enriched
        if int(feature.get("molecule_hits", 0) or 0) > 0 and float(feature.get("top_class_pct", 0.0) or 0.0) > 0.0
    ]
    prevalence_arr = np.asarray(
        [float(feature.get("top_class_pct", 0.0) or 0.0) for feature in prevalence_source],
        dtype=float,
    )
    mean_prev = float(prevalence_arr.mean()) if prevalence_arr.size else 0.0
    std_prev = float(prevalence_arr.std()) if prevalence_arr.size else 0.0
    for feature in enriched:
        prevalence = float(feature.get("top_class_pct", 0.0) or 0.0)
        has_occurrence = int(feature.get("molecule_hits", 0) or 0) > 0 and prevalence > 0.0
        feature["zscore"] = (
            (prevalence - mean_prev) / std_prev
            if has_occurrence and std_prev > 1e-12
            else 0.0
        )

    threshold_candidates = [feature["top_class_pct"] for feature in enriched if feature["zscore"] > 1.0]
    if threshold_candidates:
        threshold_pct = float(min(threshold_candidates))
        threshold_source = "zscore_gt_1"
    else:
        threshold_pct = float(max(prevalence_values)) if prevalence_values else 100.0
        threshold_source = "fallback_max_prevalence"

    for feature in enriched:
        reliable = feature["top_class_pct"] >= threshold_pct
        feature["reliable"] = bool(reliable)
        feature["assigned_class"] = feature["top_class"] if reliable else CLASS_UNRELIABLE

    reliable_features = [feature for feature in enriched if feature["reliable"]]
    reliable_ids = [feature["feature_id"] for feature in reliable_features if feature["feature_id"] in feature_lookup]
    model_name = "Unavailable"
    model_note = "Nao foi possivel calcular importancias."
    cluster_count = 0
    label_source = "derived_clusters"
    raw_label_path = str(labels_csv or "").strip()
    label_path = resolve_fp_labels_file(workdir, labels_csv)
    label_warning = ""
    if raw_label_path and label_path is None:
        label_warning = (
            f"O CSV de rótulos '{raw_label_path}' não foi encontrado. "
            "Foi usado o fallback automático."
        )

    if len(labels) >= 2 and len(importance_ids) >= 1:
        X = np.asarray([matrix[:, feature_lookup[feature_id]] for feature_id in importance_ids], dtype=float).T
        X_train, y, cluster_count, label_source, resolved_warning = _resolve_training_labels(
            labels,
            X,
            label_path,
            labels_column=labels_column,
        )
        if resolved_warning:
            label_warning = resolved_warning
        if len(set(y.tolist())) >= 2:
            scores, model_name, model_note = _compute_feature_importances(
                X_train,
                y,
                prefer_extra_trees=(label_source != "external_csv"),
            )
            if label_source == "external_csv":
                model_note = (
                    f"Importancias calculadas com rótulos externos de '{Path(label_path).name}'. "
                    + model_note
                )
            else:
                model_note = "Importancias calculadas a partir de grupos derivados dos fingerprints. " + model_note
            if fallback_importance_rule:
                model_note += (
                    " Como nenhum bit apresentou z-score de classe > 1, a atribuicao de classe "
                    "foi mantida apenas para bits com uma unica classe no nature_breakdown, e a "
                    "analise de importancia foi restrita a features com 100% de classe prevalente "
                    "e cobertura >= 5%."
                )
            model_note += f" Features elegiveis para importancia: {importance_eligible_count}."
            if label_warning:
                model_note += " " + label_warning
            max_score = float(np.max(scores)) if scores.size else 0.0
            for rank, (feature_id, score) in enumerate(
                sorted(zip(reliable_ids, scores.tolist()), key=lambda item: (-item[1], item[0])),
                start=1,
            ):
                feature = next(item for item in enriched if item["feature_id"] == feature_id)
                feature["importance_score"] = float(score)
                feature["importance_pct"] = (100.0 * float(score) / max_score) if max_score > 1e-12 else 0.0
                feature["importance_rank"] = rank
        else:
            model_name = "Unavailable"
            if label_source == "external_csv":
                model_note = (
                    "Os rótulos externos não forneceram pelo menos duas classes distintas para treino."
                )
            else:
                model_note = "Os fingerprints confiaveis nao formaram pelo menos dois grupos distintos."
            if fallback_importance_rule:
                model_note += (
                    " Neste caso, a atribuicao de classe foi mantida apenas para bits com uma unica "
                    "classe no nature_breakdown, e a elegibilidade para importancia exigiu 100% de "
                    "classe prevalente e cobertura >= 5%."
                )
            model_note += f" Features elegiveis para importancia: {importance_eligible_count}."
            if label_warning:
                model_note += " " + label_warning
    else:
        model_name = "Unavailable"
        model_note = "Sao necessarios pelo menos dois ligantes e uma feature confiavel para calcular importancias."
        if fallback_importance_rule:
            model_note += (
                " Como nenhum bit apresentou z-score de classe > 1, a atribuicao de classe foi "
                "mantida apenas para bits com uma unica classe no nature_breakdown, e apenas "
                "features com 100% de classe prevalente e cobertura >= 5% foram consideradas "
                "elegiveis para importancia."
            )
        model_note += f" Features elegiveis para importancia: {importance_eligible_count}."
        if label_warning:
            model_note += " " + label_warning

    important_features = [feature for feature in enriched if feature["importance_score"] > 0.0]
    if not important_features:
        important_features = sorted(
            reliable_features,
            key=lambda row: (-float(row["coverage_pct"]), -float(row["zscore"]), int(row["feature_id"])),
        )[:30]
    else:
        important_features = sorted(
            important_features,
            key=lambda row: (-float(row["importance_score"]), -float(row["coverage_pct"]), int(row["feature_id"])),
        )[:30]

    class_counts: dict[str, int] = {}
    for feature in important_features:
        class_name = str(feature["assigned_class"])
        class_counts[class_name] = class_counts.get(class_name, 0) + 1

    class_share = {
        class_name: (100.0 * count / len(important_features)) if important_features else 0.0
        for class_name, count in class_counts.items()
    }

    all_features = sorted(
        enriched,
        key=lambda row: (
            -float(row["importance_score"]),
            -float(row["top_class_pct"]),
            -float(row["coverage_pct"]),
            int(row["feature_id"]),
        ),
    )

    return {
        "ifp_type": artifact.get("ifp_type"),
        "ifp_label": artifact.get("ifp_label"),
        "source_file": str(source_file),
        "threshold_pct": float(threshold_pct),
        "threshold_source": threshold_source,
        "class_order": list(FP_CLASS_ORDER),
        "class_colors": dict(FP_CLASS_COLORS),
        "features": all_features,
        "important_features": important_features,
        "class_counts": class_counts,
        "class_share": class_share,
        "model_name": model_name,
        "model_note": model_note,
        "cluster_count": int(cluster_count),
        "total_molecules": len(labels),
        "label_source": label_source,
        "labels_csv": str(label_path) if label_path else "",
        "labels_column": str(labels_column or "").strip(),
    }


def _resolve_training_labels(
    labels: list[str],
    matrix: np.ndarray,
    label_path: Path | None,
    labels_column: str | None = None,
) -> tuple[np.ndarray, np.ndarray, int, str, str]:
    if label_path is not None:
        try:
            label_map = load_external_fp_labels(label_path, label_column=labels_column)
            y_values = [label_map.get(label, "") for label in labels]
            matched = [idx for idx, value in enumerate(y_values) if value]
            distinct = {y_values[idx] for idx in matched}
            if len(matched) >= 2 and len(distinct) >= 2:
                X_train = np.asarray(matrix[matched, :], dtype=float)
                y_train = np.asarray([y_values[idx] for idx in matched], dtype=object)
                return X_train, y_train, len(distinct), "external_csv", ""
            warning = (
                f"O CSV de rótulos '{label_path.name}' não cobre amostras suficientes "
                "ou não possui pelo menos duas classes. Foi usado o fallback automático."
            )
        except Exception as exc:
            warning = (
                f"Não foi possível ler o CSV de rótulos '{label_path.name}' "
                f"({type(exc).__name__}). Foi usado o fallback automático."
            )
        y, cluster_count = _derive_cluster_labels(labels, matrix)
        return np.asarray(matrix, dtype=float), y, cluster_count, "derived_clusters", warning

    y, cluster_count = _derive_cluster_labels(labels, matrix)
    return np.asarray(matrix, dtype=float), y, cluster_count, "derived_clusters", ""


def _derive_cluster_labels(labels: list[str], matrix: np.ndarray) -> tuple[np.ndarray, int]:
    matrix = np.asarray(matrix, dtype=float)
    n_samples = len(labels)
    if n_samples < 2:
        return np.ones(n_samples, dtype=int), 1

    binary = (matrix > 0).astype(float)
    similarity = _tanimoto_similarity(binary)
    requested_clusters = max(2, min(4, int(round(np.sqrt(n_samples)))))
    cluster_result = cluster_similarity_matrix(labels, similarity, n_clusters=requested_clusters)
    return np.asarray(cluster_result.cluster_ids, dtype=int), int(cluster_result.n_clusters)


def _tanimoto_similarity(binary_matrix: np.ndarray) -> np.ndarray:
    binary = np.asarray(binary_matrix) > 0
    if binary.ndim != 2:
        raise ValueError("A matriz binaria deve ser bidimensional.")

    if binary.shape[1] == 0:
        matrix = np.zeros((binary.shape[0], binary.shape[0]), dtype=float)
        np.fill_diagonal(matrix, 1.0)
        return matrix

    counts = binary.sum(axis=1).astype(float)
    similarity = np.zeros((binary.shape[0], binary.shape[0]), dtype=float)
    # Avoid BLAS-backed matrix multiplication here. Some Windows conda builds
    # used with LUNA can terminate the process during ``binary @ binary.T``.
    block_size = 64
    for start in range(0, binary.shape[0], block_size):
        stop = min(start + block_size, binary.shape[0])
        block = binary[start:stop]
        intersections = np.count_nonzero(block[:, None, :] & binary[None, :, :], axis=2).astype(float)
        unions = counts[start:stop, None] + counts[None, :] - intersections
        similarity[start:stop, :] = np.divide(
            intersections,
            unions,
            out=np.zeros_like(intersections, dtype=float),
            where=unions > 0,
        )
    np.fill_diagonal(similarity, 1.0)
    return similarity


def count_tanimoto_similarity(count_matrix: np.ndarray) -> np.ndarray:
    matrix = np.asarray(count_matrix, dtype=float)
    if matrix.ndim != 2:
        raise ValueError("A matriz de contagens deve ser bidimensional.")

    n_rows = matrix.shape[0]
    if matrix.shape[1] == 0:
        result = np.zeros((n_rows, n_rows), dtype=float)
        np.fill_diagonal(result, 1.0)
        return result

    result = np.zeros((n_rows, n_rows), dtype=float)
    for i in range(n_rows):
        result[i, i] = 1.0
        for j in range(i + 1, n_rows):
            minima = np.minimum(matrix[i], matrix[j]).sum()
            maxima = np.maximum(matrix[i], matrix[j]).sum()
            value = float(minima / maxima) if maxima > 1e-12 else 1.0
            result[i, j] = value
            result[j, i] = value
    return result


def _compute_feature_importances(
    X: np.ndarray,
    y: np.ndarray,
    prefer_extra_trees: bool = False,
) -> tuple[np.ndarray, str, str]:
    classes, y_encoded = np.unique(np.asarray(y), return_inverse=True)
    primary_error = None
    if not prefer_extra_trees:
        try:
            from sklearn.ensemble import GradientBoostingClassifier

            model = GradientBoostingClassifier(
                n_estimators=200,
                learning_rate=0.06,
                max_depth=3,
                random_state=0,
            )
            model.fit(X, y_encoded)
            scores = np.asarray(getattr(model, "feature_importances_", np.zeros(X.shape[1])), dtype=float)
            return scores, "GradientBoosting", "Importancias calculadas com GradientBoostingClassifier sobre grupos derivados dos fingerprints."
        except Exception as exc:
            primary_error = exc

    try:
        from sklearn.ensemble import ExtraTreesClassifier

        model = ExtraTreesClassifier(
            n_estimators=256,
            random_state=0,
            class_weight="balanced",
        )
        model.fit(X, y_encoded)
        scores = np.asarray(getattr(model, "feature_importances_", np.zeros(X.shape[1])), dtype=float)
        if prefer_extra_trees:
            note = "Importancias calculadas com ExtraTrees conforme configuracao padrao sem rotulos externos."
        else:
            note = "Importancias calculadas com ExtraTrees porque o GradientBoosting nao ficou disponivel."
        if primary_error is not None and not prefer_extra_trees:
            note += f" Falha do GradientBoosting: {type(primary_error).__name__}."
        return scores, "ExtraTrees fallback", note
    except Exception as exc:
        if primary_error is not None:
            note = (
                f"Nao foi possivel calcular importancias com GradientBoosting ({type(primary_error).__name__}) "
                f"nem com ExtraTrees ({type(exc).__name__})."
            )
        else:
            note = f"Nao foi possivel calcular importancias ({type(exc).__name__})."
        return np.zeros(X.shape[1], dtype=float), "Unavailable", note


def build_fp_analysis_dashboard(
    workdir: str | Path,
    artifact: dict,
    labels_csv: str | Path | None = None,
    labels_id_column: str | None = None,
    labels_column: str | None = None,
    algorithm_preference: str = "gradient_boosting",
    task_kind_preference: str | None = None,
    use_otsu_threshold: bool = False,
) -> dict:
    source_file = resolve_ifp_source_file(workdir, artifact)
    if source_file is None:
        raise FileNotFoundError("Nenhum arquivo de fingerprint foi encontrado para esta analise.")

    labels, feature_ids, matrix = load_ifp_sparse_matrix(source_file)
    feature_lookup = {int(feature_id): col for col, feature_id in enumerate(feature_ids)}
    detail_artifact = load_fp_detail_artifact(workdir, artifact.get("ifp_type"))
    detail_lookup = (detail_artifact or {}).get("feature_details", {}) or {}

    raw_features = artifact.get("features", []) or []
    enriched: list[dict] = []
    prevalence_values: list[float] = []
    for row in raw_features:
        feature_id = int(row.get("feature_id", 0))
        breakdown = normalize_fp_breakdown(row.get("nature_breakdown"))
        if not breakdown:
            breakdown = {CLASS_UNRELIABLE: 1}

        total_class_hits = sum(breakdown.values())
        top_class = max(
            sorted(breakdown),
            key=lambda key: (breakdown[key], -FP_CLASS_RANK.get(key, len(FP_CLASS_ORDER)), key),
        )
        top_pct = (100.0 * breakdown[top_class] / total_class_hits) if total_class_hits else 0.0
        class_percentages = {
            class_name: (100.0 * breakdown.get(class_name, 0) / total_class_hits) if total_class_hits else 0.0
            for class_name in FP_CLASS_ORDER
            if breakdown.get(class_name, 0) > 0
        }

        feature_matrix_col = feature_lookup.get(feature_id)
        coverage_pct = float(row.get("coverage_pct", 0.0) or 0.0)
        molecule_hits = int(row.get("molecule_hits", 0) or 0)
        entry_counts: dict[str, float] = {}
        if feature_matrix_col is not None and len(labels) > 0:
            column = np.asarray(matrix[:, feature_matrix_col], dtype=float)
            entry_counts = {
                labels[idx]: float(value)
                for idx, value in enumerate(column.tolist())
                if float(value) > 0.0
            }
            molecule_hits = len(entry_counts)
            coverage_pct = 100.0 * float(molecule_hits) / float(len(labels))

        detail = detail_lookup.get(str(feature_id), {}) or {}
        entry_details = {
            str(entry_name): {
                "shell_count": int((entry_info or {}).get("shell_count", 0) or 0),
                "interaction_counts": {
                    str(key): int(value)
                    for key, value in ((entry_info or {}).get("interaction_counts") or {}).items()
                    if int(value) > 0
                },
                "residue_counts": {
                    str(key): int(value)
                    for key, value in ((entry_info or {}).get("residue_counts") or {}).items()
                    if int(value) > 0
                },
                "pair_counts": {
                    str(key): int(value)
                    for key, value in ((entry_info or {}).get("pair_counts") or {}).items()
                    if int(value) > 0
                },
            }
            for entry_name, entry_info in (detail.get("entries") or {}).items()
        }
        interaction_breakdown = {
            str(key): int(value)
            for key, value in (detail.get("interaction_counts") or {}).items()
            if int(value) > 0
        }
        residue_breakdown = {
            str(key): int(value)
            for key, value in (detail.get("residue_counts") or {}).items()
            if int(value) > 0
        }
        pair_breakdown = {
            str(key): int(value)
            for key, value in (detail.get("pair_counts") or {}).items()
            if int(value) > 0
        }
        prevalence_values.append(top_pct)
        enriched.append(
            {
                "feature_id": feature_id,
                "molecule_hits": int(molecule_hits),
                "coverage_pct": float(coverage_pct),
                "top_class": top_class,
                "top_class_pct": float(top_pct),
                "collision_hits": int(row.get("collision_hits", 0) or 0),
                "total_count": int(row.get("total_count", 0) or 0),
                "class_breakdown": breakdown,
                "class_percentages": class_percentages,
                "missing_molecules": max(0, len(labels) - int(molecule_hits)),
                "zscore": 0.0,
                "assigned_class": CLASS_UNRELIABLE,
                "reliable": False,
                "importance_score": 0.0,
                "importance_pct": 0.0,
                "importance_rank": 0,
                "importance_zscore": 0.0,
                "importance_pvalue": 1.0,
                "entry_counts": entry_counts,
                "interaction_breakdown": interaction_breakdown,
                "residue_breakdown": residue_breakdown,
                "pair_breakdown": pair_breakdown,
                "entry_details": entry_details,
                "prevalent_interaction": "",
                "prevalent_interaction_pct": 0.0,
                "prevalent_interaction_count": 0,
                "prevalent_interaction_zscore": 0.0,
                "prevalent_residue": "",
                "prevalent_residue_pct": 0.0,
                "prevalent_residue_count": 0,
                "prevalent_residue_zscore": 0.0,
                "prevalent_pair": "",
                "prevalent_pair_pct": 0.0,
                "prevalent_pair_count": 0,
                "prevalent_pair_zscore": 0.0,
                "prevalent_pair_entries": [],
                "prevalent_interaction_entries": [],
            }
        )

    if not enriched:
        return {
            "ifp_type": artifact.get("ifp_type"),
            "ifp_label": artifact.get("ifp_label"),
            "source_file": str(source_file),
            "features": [],
            "important_features": [],
            "class_share": {},
            "class_counts": {},
            "threshold_pct": 100.0,
            "model_name": "Unavailable",
            "model_note": "Sem features para analisar.",
            "cluster_count": 0,
            "class_order": list(FP_CLASS_ORDER),
            "entry_labels": list(labels),
            "detail_available": bool(detail_artifact),
        }

    class_prevalence_population = np.asarray(
        [float(value) for value in prevalence_values if float(value) > 0.0],
        dtype=float,
    )
    class_prevalence_mean = float(class_prevalence_population.mean()) if class_prevalence_population.size else 0.0
    class_prevalence_std = float(class_prevalence_population.std()) if class_prevalence_population.size else 0.0
    threshold_pct, threshold_source, prevalence_zscores = _threshold_with_otsu_fallback(
        prevalence_values,
        use_otsu_threshold=True,
    )
    for feature, zscore in zip(enriched, prevalence_zscores):
        feature["zscore"] = float(zscore)
        if threshold_source in {"zscore_gt_1", "otsu", "otsu_single_value"}:
            reliable = (
                int(feature.get("molecule_hits", 0) or 0) > 0
                and float(feature.get("top_class_pct", 0.0) or 0.0) >= threshold_pct
            )
        else:
            class_breakdown = dict(feature.get("class_breakdown", {}) or {})
            reliable = (
                int(feature.get("molecule_hits", 0) or 0) > 0
                and len(class_breakdown) == 1
                and str(feature.get("top_class", "") or "") != CLASS_UNRELIABLE
            )
        feature["reliable"] = bool(reliable)
        feature["assigned_class"] = feature["top_class"] if reliable else CLASS_UNRELIABLE

    fallback_importance_rule = threshold_source != "zscore_gt_1"
    for feature in enriched:
        importance_eligible = (
            feature["feature_id"] in feature_lookup
            and int(feature.get("molecule_hits", 0) or 0) > 0
            and bool(feature.get("reliable", False))
        )
        feature["importance_eligible"] = bool(importance_eligible)

    importance_features = [feature for feature in enriched if feature["importance_eligible"]]
    importance_ids = [feature["feature_id"] for feature in importance_features if feature["feature_id"] in feature_lookup]
    importance_eligible_count = len(importance_features)
    model_name = "Unavailable"
    model_note = "Nao foi possivel calcular importancias."
    cluster_count = 0
    label_source = "derived_clusters"
    label_kind = "classification"
    matched_molecules = 0
    raw_label_path = str(labels_csv or "").strip()
    label_path = resolve_fp_labels_file(workdir, labels_csv)
    label_warning = ""
    if raw_label_path and label_path is None:
        label_warning = (
            f"O CSV de rÃ³tulos '{raw_label_path}' nÃ£o foi encontrado. "
            "Foi usado o fallback automÃ¡tico."
        )

    if len(labels) >= 2 and len(importance_ids) >= 1:
        X = np.asarray([matrix[:, feature_lookup[feature_id]] for feature_id in importance_ids], dtype=float).T
        X_train, y, cluster_count, label_source, resolved_warning, label_kind, matched_molecules = _resolve_training_labels(
            labels,
            X,
            label_path,
            labels_id_column=labels_id_column,
            labels_column=labels_column,
            task_kind_preference=task_kind_preference,
        )
        if resolved_warning:
            label_warning = resolved_warning
        if len(set(np.asarray(y).tolist())) >= 2:
            scores, model_name, model_note = _compute_feature_importances(
                X_train,
                y,
                task_kind=label_kind,
                algorithm_preference=algorithm_preference,
            )
            if label_source == "external_csv":
                model_note = (
                    f"Importancias calculadas com rótulos externos de '{Path(label_path).name}'. "
                    + model_note
                )
            else:
                model_note = "Importancias calculadas a partir de grupos derivados dos fingerprints. " + model_note
            if fallback_importance_rule:
                if threshold_source in {"otsu", "otsu_single_value"}:
                    model_note += (
                        " Como nenhum bit apresentou z-score de classe > 1, a atribuicao de classe "
                        "usou Otsu's Thresholding como limiar alternativo, e a analise de importancia "
                        "usou apenas as features confiaveis."
                    )
                else:
                    model_note += (
                        " Como nenhum bit apresentou z-score de classe > 1, a atribuicao de classe "
                        "foi mantida apenas para bits com uma unica classe no nature_breakdown, mas a "
                        "analise de importancia usou apenas as features confiaveis."
                    )
            model_note += f" Features elegiveis para importancia: {importance_eligible_count}."
            if label_warning:
                model_note += " " + label_warning
            max_score = float(np.max(scores)) if scores.size else 0.0
            for rank, (feature_id, score) in enumerate(
                sorted(zip(importance_ids, scores.tolist()), key=lambda item: (-item[1], item[0])),
                start=1,
            ):
                feature = next(item for item in enriched if item["feature_id"] == feature_id)
                feature["importance_score"] = float(score)
                feature["importance_pct"] = (100.0 * float(score) / max_score) if max_score > 1e-12 else 0.0
                feature["importance_rank"] = rank
        else:
            model_name = "Unavailable"
            if label_source == "external_csv":
                model_note = "Os rÃ³tulos externos nÃ£o forneceram valores suficientes para treino."
            else:
                model_note = "Os fingerprints confiaveis nao formaram pelo menos dois grupos distintos."
            if fallback_importance_rule:
                if threshold_source in {"otsu", "otsu_single_value"}:
                    model_note += (
                        " Neste caso, a atribuicao de classe usou Otsu's Thresholding como limiar "
                        "alternativo, mas a importancia ficou restrita as features confiaveis."
                    )
                else:
                    model_note += (
                        " Neste caso, a atribuicao de classe foi mantida apenas para bits com uma unica "
                        "classe no nature_breakdown, mas a importancia ficou restrita as features confiaveis."
                    )
            model_note += f" Features elegiveis para importancia: {importance_eligible_count}."
            if label_warning:
                model_note += " " + label_warning
    else:
        model_name = "Unavailable"
        model_note = "Sao necessarios pelo menos dois ligantes e uma feature confiavel para calcular importancias."
        if fallback_importance_rule:
            if threshold_source in {"otsu", "otsu_single_value"}:
                model_note += (
                    " Como nenhum bit apresentou z-score de classe > 1, a atribuicao de classe usou "
                    "Otsu's Thresholding como limiar alternativo."
                )
            else:
                model_note += (
                    " Como nenhum bit apresentou z-score de classe > 1, a atribuicao de classe foi "
                    "mantida apenas para bits com uma unica classe no nature_breakdown."
                )
        model_note += f" Features elegiveis para importancia: {importance_eligible_count}."
        if label_warning:
            model_note += " " + label_warning

    importance_population = [
        feature
        for feature in enriched
        if feature["importance_eligible"]
    ]
    importance_scores = [float(feature.get("importance_score", 0.0) or 0.0) for feature in importance_population]
    importance_mean = float(np.mean(importance_scores)) if importance_scores else 0.0
    importance_std = float(np.std(importance_scores)) if importance_scores else 0.0
    for feature in enriched:
        score = float(feature.get("importance_score", 0.0) or 0.0)
        if feature["importance_eligible"] and importance_std > 1e-12:
            feature["importance_zscore"] = (score - importance_mean) / importance_std
        else:
            feature["importance_zscore"] = 0.0
        feature["importance_pvalue"] = _importance_p_value(feature.get("importance_zscore", 0.0) or 0.0)

    important_features = sorted(
        [
            feature for feature in enriched
            if float(feature.get("importance_pvalue", 1.0) or 1.0) < 0.01
        ],
        key=lambda row: (
            float(row.get("importance_pvalue", 1.0) or 1.0),
            -float(row.get("importance_score", 0.0) or 0.0),
            int(row.get("feature_id", 0) or 0),
        ),
    )
    important_selection = "pvalue_lt_0.01"

    class_counts: dict[str, int] = {}
    for feature in important_features:
        class_name = str(feature["assigned_class"])
        class_counts[class_name] = class_counts.get(class_name, 0) + 1
    class_share = {
        class_name: (100.0 * count / len(important_features)) if important_features else 0.0
        for class_name, count in class_counts.items()
    }

    interaction_features: list[dict] = []
    for feature in important_features:
        if str(feature.get("assigned_class", "")) != CLASS_NONCOVALENT:
            continue
        interaction_breakdown, residue_breakdown, pair_breakdown, entry_details = _filtered_prevalence_counts(feature)
        if sum(interaction_breakdown.values()) <= 0:
            continue
        feature["interaction_breakdown"] = interaction_breakdown
        feature["residue_breakdown"] = residue_breakdown
        feature["pair_breakdown"] = pair_breakdown
        feature["entry_details"] = entry_details
        interaction_features.append(feature)
    interaction_top_pcts: list[float] = []
    residue_top_pcts: list[float] = []
    pair_top_pcts: list[float] = []
    for feature in interaction_features:
        interaction_name, interaction_count = _dominant_label(feature.get("interaction_breakdown") or {})
        interaction_total = sum((feature.get("interaction_breakdown") or {}).values())
        interaction_pct = (100.0 * interaction_count / interaction_total) if interaction_total else 0.0
        residue_name, residue_count = _dominant_label(feature.get("residue_breakdown") or {})
        residue_total = sum((feature.get("residue_breakdown") or {}).values())
        residue_pct = (100.0 * residue_count / residue_total) if residue_total else 0.0
        pair_name, pair_count = _dominant_label(feature.get("pair_breakdown") or {})
        pair_total = sum((feature.get("pair_breakdown") or {}).values())
        pair_pct = (100.0 * pair_count / pair_total) if pair_total else 0.0
        feature["prevalent_interaction"] = interaction_name
        feature["prevalent_interaction_count"] = int(interaction_count)
        feature["prevalent_interaction_pct"] = float(interaction_pct)
        feature["prevalent_residue"] = residue_name
        feature["prevalent_residue_count"] = int(residue_count)
        feature["prevalent_residue_pct"] = float(residue_pct)
        feature["prevalent_pair"] = pair_name
        feature["prevalent_pair_count"] = int(pair_count)
        feature["prevalent_pair_pct"] = float(pair_pct)
        interaction_top_pcts.append(interaction_pct)
        residue_top_pcts.append(residue_pct)
        pair_top_pcts.append(pair_pct)

    interaction_threshold_pct, interaction_threshold_source, interaction_zscores = _threshold_with_otsu_fallback(
        interaction_top_pcts,
        use_otsu_threshold=use_otsu_threshold,
    )
    residue_threshold_pct, residue_threshold_source, residue_zscores = _threshold_with_otsu_fallback(
        residue_top_pcts,
        use_otsu_threshold=use_otsu_threshold,
    )
    pair_threshold_pct, pair_threshold_source, pair_zscores = _threshold_with_otsu_fallback(
        pair_top_pcts,
        use_otsu_threshold=use_otsu_threshold,
    )
    for feature, interaction_zscore, residue_zscore, pair_zscore in zip(
        interaction_features,
        interaction_zscores,
        residue_zscores,
        pair_zscores,
    ):
        feature["prevalent_interaction_zscore"] = float(interaction_zscore)
        feature["prevalent_residue_zscore"] = float(residue_zscore)
        feature["prevalent_pair_zscore"] = float(pair_zscore)
        if float(feature["prevalent_interaction_pct"]) < interaction_threshold_pct:
            feature["prevalent_interaction"] = CLASS_UNRELIABLE
        if float(feature["prevalent_residue_pct"]) < residue_threshold_pct:
            feature["prevalent_residue"] = CLASS_UNRELIABLE
        if float(feature["prevalent_pair_pct"]) < pair_threshold_pct:
            feature["prevalent_pair"] = CLASS_UNRELIABLE
        else:
            pair_name = str(feature["prevalent_pair"])
            interaction_name = pair_name.split("||", 1)[0]
            feature["prevalent_pair_entries"] = [
                entry_name
                for entry_name, entry_info in (feature.get("entry_details") or {}).items()
                if pair_name in ((entry_info or {}).get("pair_counts") or {})
            ]
        if str(feature.get("prevalent_interaction", "")) != CLASS_UNRELIABLE:
            interaction_name = str(feature["prevalent_interaction"])
            feature["prevalent_interaction_entries"] = [
                entry_name
                for entry_name, entry_info in (feature.get("entry_details") or {}).items()
                if interaction_name in ((entry_info or {}).get("interaction_counts") or {})
            ]

    all_features = sorted(
        enriched,
        key=lambda row: (
            float(row.get("importance_pvalue", 1.0) or 1.0),
            -float(row.get("importance_score", 0.0) or 0.0),
            -float(row.get("top_class_pct", 0.0) or 0.0),
            -float(row.get("coverage_pct", 0.0) or 0.0),
            int(row.get("feature_id", 0) or 0),
        ),
    )

    return {
        "ifp_type": artifact.get("ifp_type"),
        "ifp_label": artifact.get("ifp_label"),
        "source_file": str(source_file),
        "threshold_pct": float(threshold_pct),
        "threshold_source": threshold_source,
        "class_order": list(FP_CLASS_ORDER),
        "class_colors": dict(FP_CLASS_COLORS),
        "features": all_features,
        "important_features": important_features,
        "class_counts": class_counts,
        "class_share": class_share,
        "model_name": model_name,
        "model_note": model_note,
        "cluster_count": int(cluster_count),
        "total_molecules": len(labels),
        "label_source": label_source,
        "label_kind": label_kind,
        "labels_csv": str(label_path) if label_path else "",
        "labels_id_column": str(labels_id_column or "").strip(),
        "labels_column": str(labels_column or "").strip(),
        "matched_molecules": int(matched_molecules),
        "importance_eligible_count": int(importance_eligible_count),
        "class_zscore_mean": float(class_prevalence_mean),
        "class_zscore_std": float(class_prevalence_std),
        "importance_zscore_mean": float(importance_mean),
        "importance_zscore_std": float(importance_std),
        "entry_labels": list(labels),
        "interaction_threshold_pct": float(interaction_threshold_pct),
        "interaction_threshold_source": interaction_threshold_source,
        "residue_threshold_pct": float(residue_threshold_pct),
        "residue_threshold_source": residue_threshold_source,
        "pair_threshold_pct": float(pair_threshold_pct),
        "pair_threshold_source": pair_threshold_source,
        "detail_available": bool(detail_artifact),
        "important_selection": important_selection,
        "algorithm_preference": str(algorithm_preference or "gradient_boosting"),
        "task_kind_preference": str(task_kind_preference or ""),
        "use_otsu_threshold": bool(use_otsu_threshold or threshold_source in {"otsu", "otsu_single_value"}),
    }


def _resolve_training_labels(
    labels: list[str],
    matrix: np.ndarray,
    label_path: Path | None,
    labels_id_column: str | None = None,
    labels_column: str | None = None,
    task_kind_preference: str | None = None,
) -> tuple[np.ndarray, np.ndarray, int, str, str, str, int]:
    if label_path is not None:
        try:
            label_map = load_external_fp_labels(
                label_path,
                label_column=labels_column,
                id_column=labels_id_column,
            )
            matched: list[int] = []
            y_values: list[object] = []
            numeric_values: list[float] = []
            numeric_ok = True
            for idx, label in enumerate(labels):
                value = _resolve_external_label_value(label_map, label)
                if value == "":
                    continue
                matched.append(idx)
                y_values.append(value)
                parsed = _safe_float_or_none(value)
                if parsed is None:
                    numeric_ok = False
                else:
                    numeric_values.append(parsed)

            requested_task = str(task_kind_preference or "").strip().lower()
            if len(matched) >= 2:
                X_train = np.asarray(matrix[matched, :], dtype=float)
                if requested_task == "classification":
                    distinct = {str(value) for value in y_values}
                    if len(distinct) >= 2:
                        return (
                            X_train,
                            np.asarray([str(value) for value in y_values], dtype=object),
                            len(distinct),
                            "external_csv",
                            "",
                            "classification",
                            len(matched),
                        )
                elif requested_task == "regression":
                    if numeric_ok and len(set(numeric_values)) >= 2:
                        return (
                            X_train,
                            np.asarray(numeric_values, dtype=float),
                            len(set(numeric_values)),
                            "external_csv",
                            "",
                            "regression",
                            len(matched),
                        )
                elif numeric_ok and len(set(numeric_values)) >= 2:
                    return (
                        X_train,
                        np.asarray(numeric_values, dtype=float),
                        len(set(numeric_values)),
                        "external_csv",
                        "",
                        "regression",
                        len(matched),
                    )
                distinct = {str(value) for value in y_values}
                if len(distinct) >= 2:
                    return (
                        X_train,
                        np.asarray([str(value) for value in y_values], dtype=object),
                        len(distinct),
                        "external_csv",
                        "",
                        "classification",
                        len(matched),
                    )
            warning = (
                f"O CSV de rÃ³tulos '{label_path.name}' nÃ£o cobre amostras suficientes "
                "ou nÃ£o possui variacao suficiente. Foi usado o fallback automÃ¡tico."
            )
        except Exception as exc:
            warning = (
                f"NÃ£o foi possÃ­vel ler o CSV de rÃ³tulos '{label_path.name}' "
                f"({type(exc).__name__}). Foi usado o fallback automÃ¡tico."
            )
        y, cluster_count = _derive_cluster_labels(labels, matrix)
        return np.asarray(matrix, dtype=float), y, cluster_count, "derived_clusters", warning, "classification", len(labels)

    y, cluster_count = _derive_cluster_labels(labels, matrix)
    return np.asarray(matrix, dtype=float), y, cluster_count, "derived_clusters", "", "classification", len(labels)


def _compute_feature_importances(
    X: np.ndarray,
    y: np.ndarray,
    task_kind: str = "classification",
    algorithm_preference: str = "gradient_boosting",
) -> tuple[np.ndarray, str, str]:
    preference = str(algorithm_preference or "auto").strip().lower()
    if preference == "extra_trees":
        algorithm_order = ["extra_trees", "gradient_boosting"]
    else:
        algorithm_order = ["gradient_boosting", "extra_trees"]
    errors: list[tuple[str, Exception]] = []

    if str(task_kind) == "regression":
        for algorithm in algorithm_order:
            try:
                if algorithm == "gradient_boosting":
                    from sklearn.ensemble import GradientBoostingRegressor

                    model = GradientBoostingRegressor(
                        n_estimators=220,
                        learning_rate=0.05,
                        max_depth=3,
                        random_state=0,
                    )
                    model_label = "GradientBoosting"
                    class_label = "GradientBoostingRegressor"
                else:
                    from sklearn.ensemble import ExtraTreesRegressor

                    model = ExtraTreesRegressor(
                        n_estimators=256,
                        random_state=0,
                    )
                    model_label = "ExtraTrees"
                    class_label = "ExtraTreesRegressor"
                model.fit(X, np.asarray(y, dtype=float))
                scores = np.asarray(getattr(model, "feature_importances_", np.zeros(X.shape[1])), dtype=float)
                used_fallback = bool(errors)
                note = f"Importancias calculadas com {class_label} da biblioteca scikit-learn sobre variavel continua."
                if used_fallback:
                    failed_name, failed_exc = errors[0]
                    note += (
                        f" Metodo alternativo usado apos falha de {_algorithm_display_name(failed_name)} "
                        f"({type(failed_exc).__name__})."
                    )
                return scores, model_label if not used_fallback else f"{model_label} fallback", note
            except Exception as exc:
                errors.append((algorithm, exc))
        scores = _fallback_regression_importances(X, np.asarray(y, dtype=float))
        failures = "; ".join(
            f"{_algorithm_display_name(name)}={type(exc).__name__}" for name, exc in errors
        ) or "sklearn indisponivel"
        install_hint = (
            " Instale/atualize scikit-learn no luna-env pela aba Setup > Instalar LUNA."
            if errors and all(type(exc).__name__ == "ModuleNotFoundError" for _name, exc in errors)
            else ""
        )
        return scores, "Fallback regression", (
            f"Nao foi possivel calcular importancias com os algoritmos scikit-learn solicitados ({failures}). "
            "Foram usadas importancias exploratorias por correlacao absoluta."
            + install_hint
        )

    classes, y_encoded = np.unique(np.asarray(y), return_inverse=True)
    for algorithm in algorithm_order:
        try:
            if algorithm == "gradient_boosting":
                from sklearn.ensemble import GradientBoostingClassifier

                model = GradientBoostingClassifier(
                    n_estimators=200,
                    learning_rate=0.06,
                    max_depth=3,
                    random_state=0,
                )
                model_label = "GradientBoosting"
                class_label = "GradientBoostingClassifier"
            else:
                from sklearn.ensemble import ExtraTreesClassifier

                model = ExtraTreesClassifier(
                    n_estimators=256,
                    random_state=0,
                    class_weight="balanced",
                )
                model_label = "ExtraTrees"
                class_label = "ExtraTreesClassifier"
            model.fit(X, y_encoded)
            scores = np.asarray(getattr(model, "feature_importances_", np.zeros(X.shape[1])), dtype=float)
            used_fallback = bool(errors)
            note = f"Importancias calculadas com {class_label} da biblioteca scikit-learn."
            if used_fallback:
                failed_name, failed_exc = errors[0]
                note += (
                    f" Metodo alternativo usado apos falha de {_algorithm_display_name(failed_name)} "
                    f"({type(failed_exc).__name__})."
                )
            return scores, model_label if not used_fallback else f"{model_label} fallback", note
        except Exception as exc:
            errors.append((algorithm, exc))
    scores = _fallback_classification_importances(X, y_encoded)
    failures = "; ".join(
        f"{_algorithm_display_name(name)}={type(exc).__name__}" for name, exc in errors
    ) or "sklearn indisponivel"
    install_hint = (
        " Instale/atualize scikit-learn no luna-env pela aba Setup > Instalar LUNA."
        if errors and all(type(exc).__name__ == "ModuleNotFoundError" for _name, exc in errors)
        else ""
    )
    return scores, "Fallback classification", (
        f"Nao foi possivel calcular importancias com os algoritmos scikit-learn solicitados ({failures}). "
        "Foram usadas importancias exploratorias por separacao entre classes."
        + install_hint
    )


def _algorithm_display_name(name: str) -> str:
    if name == "gradient_boosting":
        return "GradientBoosting"
    if name == "extra_trees":
        return "ExtraTrees"
    return str(name)


def _fallback_regression_importances(X: np.ndarray, y: np.ndarray) -> np.ndarray:
    X = np.asarray(X, dtype=float)
    y = np.asarray(y, dtype=float)
    if X.ndim != 2 or X.shape[1] == 0 or len(y) != X.shape[0]:
        return np.zeros(X.shape[1] if X.ndim == 2 else 0, dtype=float)
    y_centered = y - float(np.mean(y))
    y_norm = float(np.linalg.norm(y_centered))
    if y_norm <= 1e-12:
        return np.zeros(X.shape[1], dtype=float)
    scores = np.zeros(X.shape[1], dtype=float)
    for idx in range(X.shape[1]):
        column = np.asarray(X[:, idx], dtype=float)
        column_centered = column - float(np.mean(column))
        denom = float(np.linalg.norm(column_centered) * y_norm)
        if denom > 1e-12:
            scores[idx] = abs(float(np.dot(column_centered, y_centered) / denom))
    return scores


def _fallback_classification_importances(X: np.ndarray, y_encoded: np.ndarray) -> np.ndarray:
    X = np.asarray(X, dtype=float)
    y_encoded = np.asarray(y_encoded, dtype=int)
    if X.ndim != 2 or X.shape[1] == 0 or len(y_encoded) != X.shape[0]:
        return np.zeros(X.shape[1] if X.ndim == 2 else 0, dtype=float)
    global_mean = np.mean(X, axis=0)
    classes = np.unique(y_encoded)
    between = np.zeros(X.shape[1], dtype=float)
    within = np.zeros(X.shape[1], dtype=float)
    for class_id in classes:
        mask = y_encoded == class_id
        class_rows = X[mask]
        if class_rows.size == 0:
            continue
        class_mean = np.mean(class_rows, axis=0)
        between += float(class_rows.shape[0]) * np.square(class_mean - global_mean)
        within += np.sum(np.square(class_rows - class_mean), axis=0)
    scores = np.divide(
        between,
        within + 1e-12,
        out=np.zeros_like(between, dtype=float),
        where=(within + 1e-12) > 0,
    )
    return scores
