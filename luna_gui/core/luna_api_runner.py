"""Run LUNA via its Python API (inside luna-env) instead of the CLI.

Used when the user enables any advanced option (water handling, custom
InteractionCalculator flags, DefaultInteractionConfig overrides, per-type
PSE filtering, project forks, or multi-IFP export). Writes a self-contained
Python script to the workdir and returns the argv needed to run it.
"""
from __future__ import annotations

import json
from pathlib import Path

from .project import ProjectConfig, resolve_ifp_output_paths, resolve_sim_matrix_output_paths


_LIGAND_SUFFIXES = ("_ligand", "-ligand", "_lig", "-lig")
_PREPARED_PROTEIN_MARKER = "REMARK   Separated Protein"


API_RUNNER_SCRIPT = r'''
"""Auto-generated LUNA runner. Reads params from argv[1] (JSON file)."""
import copy, csv, gzip, json, os, pickle, re, shutil, sys
from pathlib import Path

params_file = sys.argv[1]
with open(params_file, "r", encoding="utf-8") as fh:
    p = json.load(fh)

import luna
from luna.mol.entry import MolFileEntry
from luna.projects import LocalProject
from luna.interaction.config import DefaultInteractionConfig, InteractionConfig
from luna.interaction.filter import InteractionFilter
from luna.interaction.calc import InteractionCalculator
from luna.interaction.fp.type import IFPType
from luna.interaction.fp.shell import ShellGenerator

try:
    from rdkit import RDLogger
    RDLogger.DisableLog("rdApp.warning")
except Exception:
    pass


def _entry_key(entry):
    try:
        return entry.to_string()
    except Exception:
        return str(entry)


def _merge_entries(existing_entries, new_entries):
    merged = []
    seen = set()
    for entry in list(existing_entries) + list(new_entries):
        key = _entry_key(entry)
        if key in seen:
            continue
        seen.add(key)
        merged.append(entry)
    return merged


def _prepare_fork_tree(src_dir, dst_dir):
    if not src_dir:
        return False

    src = Path(src_dir)
    dst = Path(dst_dir)
    if not src.exists() or not src.is_dir():
        raise OSError(f"Projeto fonte inválido: {src}")

    try:
        same_location = src.resolve() == dst.resolve()
    except Exception:
        same_location = src == dst

    if same_location:
        return False

    shutil.copytree(
        src,
        dst,
        dirs_exist_ok=True,
        ignore=shutil.ignore_patterns(".luna_gui.json", "_luna_api_runner.py", "_luna_api_params.json"),
    )
    return True


def _load_existing_entries(workdir):
    project_file = LocalProject.get_project_file(workdir)
    if not os.path.exists(project_file):
        return []
    loaded = LocalProject.load(project_file, logging_enabled=True)
    return list(getattr(loaded, "entries", []) or [])


IFP_SUFFIX = {"EIFP": "E", "HIFP": "H", "FIFP": "F"}
IFP_LABELS = {"EIFP": "Extended", "HIFP": "Hybrid", "FIFP": "Functional"}


def _safe_name(value):
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value)).strip("_") or "entry"


def _group_has_water_residue(atm_grp):
    for atom in getattr(atm_grp, "atoms", []) or []:
        residue = getattr(atom, "parent", None)
        if residue is None:
            continue
        resname = getattr(residue, "resname", None)
        if resname is None:
            try:
                resname = residue.get_resname()
            except Exception:
                resname = ""
        if str(resname or "").strip().upper() in WATER_RESIDUES:
            return True
    return False


def _group_role(atm_grp):
    if atm_grp is None:
        return "unknown"
    try:
        if atm_grp.has_water():
            return "water"
    except Exception:
        pass
    try:
        if _group_has_water_residue(atm_grp):
            return "water"
    except Exception:
        pass
    try:
        if atm_grp.has_residue() or atm_grp.has_nucleotide():
            return "protein"
    except Exception:
        pass
    try:
        if atm_grp.has_hetatm():
            return "ligand"
    except Exception:
        pass
    return "other"


def _classify_shell_natures(shell):
    natures = set()
    center_role = _group_role(getattr(shell, "central_atm_grp", None))
    neighborhood = list(getattr(shell, "neighborhood", []) or [])
    neighborhood_roles = {_group_role(group) for group in neighborhood}
    interactions = list(getattr(shell, "interactions", []) or [])

    if getattr(shell, "level", 0) == 0:
        if center_role == "ligand":
            natures.add("Ligand's level 0 features only")
        elif center_role == "protein":
            natures.add("Protein's level 0 features only")
    else:
        if neighborhood_roles and neighborhood_roles <= {"ligand"}:
            natures.add("Upper level with ligand atomic information only")
        if neighborhood_roles and neighborhood_roles <= {"protein"}:
            natures.add("Upper level with protein atomic information only")

    has_ligand_protein = False
    has_ligand_water = False
    has_protein_water = False
    has_intraligand = False
    has_intraprotein = False
    for interaction in interactions:
        src_role = _group_role(getattr(interaction, "src_grp", None))
        trgt_role = _group_role(getattr(interaction, "trgt_grp", None))
        pair = {src_role, trgt_role}
        if pair == {"ligand", "protein"}:
            has_ligand_protein = True
        elif pair == {"ligand", "water"}:
            has_ligand_water = True
        elif pair == {"protein", "water"}:
            has_protein_water = True
        elif src_role == "ligand" and trgt_role == "ligand":
            has_intraligand = True
        elif src_role == "protein" and trgt_role == "protein":
            has_intraprotein = True

    has_water_mediated_protein_context = has_ligand_water and has_protein_water
    has_protein_context = has_ligand_protein or has_water_mediated_protein_context
    if has_intraligand and not has_protein_context and not has_intraprotein:
        natures.add("Intraligand interactions only")
    if has_intraprotein and not has_protein_context and not has_intraligand:
        natures.add("Intraprotein interactions only")
    if has_protein_context:
        natures.add("Has noncovalent interactions with the protein")

    if not natures:
        if center_role == "ligand":
            natures.add("Unreliable feature")
        elif center_role == "protein":
            natures.add("Unreliable feature")
        else:
            natures.add("Unreliable feature")
    return sorted(natures)


def _has_mixed_class_collision(shell_nature_sets):
    signatures = {
        tuple(sorted(str(nature) for nature in (natures or []) if str(nature).strip()))
        for natures in shell_nature_sets
        if natures
    }
    return len(signatures) > 1


DISTANCE_CAP_KEYS = [
    "max_da_dist_hb_inter",
    "max_ha_dist_hb_inter",
    "max_da_dist_whb_inter",
    "max_ha_dist_whb_inter",
    "max_dc_dist_whb_inter",
    "max_hc_dist_whb_inter",
    "max_dist_repuls_inter",
    "max_dist_attract_inter",
    "max_cc_dist_pi_pi_inter",
    "max_cc_dist_amide_pi_inter",
    "max_dist_hydrop_inter",
    "max_dist_cation_pi_inter",
    "max_xa_dist_xbond_inter",
    "max_xc_dist_xbond_inter",
    "max_ya_dist_ybond_inter",
    "max_yc_dist_ybond_inter",
    "max_ne_dist_multipolar_inter",
    "max_id_dist_ion_multipole_inter",
    "max_dist_proximal",
    "max_ma_dist_metal_coord",
    "bsite_cutoff",
    "cache_cutoff",
]


def _write_square_similarity_matrix(edge_list_path, output_path):
    edge_list_path = Path(edge_list_path)
    if not edge_list_path.exists():
        return None

    with edge_list_path.open("r", encoding="utf-8", errors="replace", newline="") as fh:
        rows = [row for row in csv.reader(fh) if any(cell.strip() for cell in row)]

    if not rows or [cell.strip().lower() for cell in rows[0][:3]] != ["entry1", "entry2", "similarity"]:
        return None

    labels = []
    seen = set()
    edges = {}
    for row in rows[1:]:
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
        edges.setdefault(key, []).append(float(row[2]))

    if not labels:
        return None

    size = len(labels)
    idx = {label: pos for pos, label in enumerate(labels)}
    matrix = [[0.0] * size for _ in range(size)]
    for pos in range(size):
        matrix[pos][pos] = 1.0

    for (entry1, entry2), values in edges.items():
        value = sum(values) / max(1, len(values))
        i = idx[entry1]
        j = idx[entry2]
        matrix[i][j] = value
        matrix[j][i] = value

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow([""] + labels)
        for label, values in zip(labels, matrix):
            writer.writerow([label] + [f"{value:.6f}" for value in values])
    return str(output_path)


def _load_ifp_sparse_rows(path):
    path = Path(path)
    if not path.exists():
        return [], []

    labels = []
    rows = []
    with path.open("r", encoding="utf-8", errors="replace", newline="") as fh:
        reader = csv.DictReader(fh)
        for record in reader:
            label = str(record.get("ligand_id") or "").strip()
            if not label:
                continue
            raw_bits = [token.strip() for token in str(record.get("on_bits") or "").split("\t") if token.strip()]
            raw_counts = [token.strip() for token in str(record.get("count") or "").split("\t") if token.strip()]
            if raw_counts and len(raw_counts) != len(raw_bits):
                continue
            counts = [float(token) for token in raw_counts] if raw_counts else [1.0] * len(raw_bits)
            row_map = {}
            for bit, count in zip(raw_bits, counts):
                row_map[int(bit)] = row_map.get(int(bit), 0.0) + float(count)
            labels.append(label)
            rows.append(row_map)
    return labels, rows


def _count_tanimoto_similarity(row_a, row_b):
    keys = set(row_a) | set(row_b)
    if not keys:
        return 1.0
    intersection = 0.0
    union = 0.0
    for key in keys:
        aval = float(row_a.get(key, 0.0) or 0.0)
        bval = float(row_b.get(key, 0.0) or 0.0)
        intersection += min(aval, bval)
        union += max(aval, bval)
    if union <= 1e-12:
        return 1.0
    return max(0.0, min(1.0, intersection / union))


def _write_similarity_outputs_from_ifp(ifp_csv_path, edge_output_path, square_output_path):
    labels, rows = _load_ifp_sparse_rows(ifp_csv_path)
    if not labels:
        return None, None

    size = len(labels)
    matrix = [[0.0] * size for _ in range(size)]
    for i in range(size):
        matrix[i][i] = 1.0
    for i in range(size):
        for j in range(i + 1, size):
            sim = _count_tanimoto_similarity(rows[i], rows[j])
            matrix[i][j] = sim
            matrix[j][i] = sim

    edge_output_path = Path(edge_output_path)
    edge_output_path.parent.mkdir(parents=True, exist_ok=True)
    with edge_output_path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["entry1", "entry2", "similarity"])
        for i in range(size):
            for j in range(i + 1, size):
                writer.writerow([labels[i], labels[j], f"{matrix[i][j]:.6f}"])

    square_output_path = Path(square_output_path)
    square_output_path.parent.mkdir(parents=True, exist_ok=True)
    with square_output_path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow([""] + labels)
        for label, values in zip(labels, matrix):
            writer.writerow([label] + [f"{value:.6f}" for value in values])
    return str(edge_output_path), str(square_output_path)


def _entry_meta(entry):
    meta = {
        "kind": entry.__class__.__name__,
        "entry_str": _entry_key(entry),
        "pdb_id": getattr(entry, "pdb_id", None),
        "sep": getattr(entry, "sep", ":"),
        "is_hetatm": bool(getattr(entry, "is_hetatm", True)),
    }
    if isinstance(entry, MolFileEntry):
        meta.update(
            {
                "mol_id": getattr(entry, "mol_id", None),
                "mol_file": getattr(entry, "mol_file", None),
                "mol_file_ext": getattr(entry, "mol_file_ext", None),
                "mol_obj_type": getattr(entry, "mol_obj_type", None),
                "overwrite_mol_name": bool(getattr(entry, "overwrite_mol_name", False)),
                "is_multimol_file": bool(getattr(entry, "is_multimol_file", False)),
            }
        )
    return meta


def _detach_shell_for_storage(shell):
    try:
        stored_shell = copy.copy(shell)
    except Exception:
        stored_shell = shell
    try:
        stored_shell._manager = None
    except Exception:
        pass

    groups = [getattr(stored_shell, "central_atm_grp", None)]
    try:
        groups.extend(list(getattr(stored_shell, "neighborhood", []) or []))
    except Exception:
        pass

    for group in groups:
        if group is None:
            continue
        try:
            group._manager = None
        except Exception:
            pass

    for interaction in list(getattr(stored_shell, "interactions", []) or []):
        for attr in ("src_grp", "trgt_grp"):
            group = getattr(interaction, attr, None)
            if group is None:
                continue
            try:
                group._manager = None
            except Exception:
                pass

    return stored_shell


def _save_feature_shell_payload(path, entry, feature_shells, pdb_dir):
    payload = {
        "entry_meta": _entry_meta(entry),
        "feature_shells": feature_shells,
        "pdb_dir": pdb_dir,
    }
    with gzip.open(path, "wb") as fh:
        pickle.dump(payload, fh, protocol=pickle.HIGHEST_PROTOCOL)


def _export_fp_artifacts(proj, params, type_name):
    workdir = Path(params["workdir"])
    suffix = IFP_SUFFIX.get(type_name, type_name)
    shell_dir = workdir / "results" / "fingerprints" / "_shells" / suffix
    shell_dir.mkdir(parents=True, exist_ok=True)

    total_entries = 0
    feature_summary = {}
    entry_index = {}

    for entry in list(getattr(proj, "entries", []) or []):
        entry_name = getattr(entry, "to_string", lambda: str(entry))()
        try:
            entry_results = proj.get_entry_results(entry)
            atm_grps_mngr = entry_results.atm_grps_mngr
        except Exception as ex:
            entry_index[entry_name] = [{"error": str(ex)}]
            continue

        sg = ShellGenerator(
            proj.ifp_num_levels,
            proj.ifp_radius_step,
            diff_comp_classes=getattr(proj, "ifp_diff_comp_classes", True),
            ifp_type=proj.ifp_type,
        )
        sm = sg.create_shells(atm_grps_mngr)
        unique_shells = not proj.ifp_count
        fp = sm.to_fingerprint(
            fold_to_length=proj.ifp_length,
            unique_shells=unique_shells,
            count_fp=proj.ifp_count,
        )
        feature_shells = {}

        if proj.ifp_count:
            bit_counts = {int(bit): int(count) for bit, count in fp.counts.items()}
        else:
            bit_counts = {int(bit): 1 for bit in fp.get_on_bits()}

        total_entries += 1
        entry_rows = []
        for feature_id, count in sorted(bit_counts.items()):
            traced = list(sm.trace_back_feature(feature_id, fp, unique_shells=unique_shells))
            raw_collision = len(traced) > 1 or any(len(found_shells) > 1 for _, found_shells in traced)
            original_features = []
            nature_tags = set()
            stored_shells = []
            shell_nature_sets = []

            for ori_feature, found_shells in traced:
                original_features.append(int(ori_feature))
                for shell in found_shells:
                    shell_natures = _classify_shell_natures(shell)
                    shell_nature_sets.append(shell_natures)
                    nature_tags.update(shell_natures)
                    stored_shells.append(_detach_shell_for_storage(shell))

            collision = bool(raw_collision and _has_mixed_class_collision(shell_nature_sets))
            if collision:
                nature_tags.add("Features with collision in the same complex")
            if not nature_tags:
                nature_tags.add("Unreliable feature")
            if stored_shells:
                feature_shells[str(int(feature_id))] = stored_shells

            feature_info = feature_summary.setdefault(
                int(feature_id),
                {"molecule_hits": 0, "total_count": 0, "collision_hits": 0, "nature_counts": {}},
            )
            feature_info["molecule_hits"] += 1
            feature_info["total_count"] += int(count)
            if collision:
                feature_info["collision_hits"] += 1
            for nature in nature_tags:
                feature_info["nature_counts"][nature] = feature_info["nature_counts"].get(nature, 0) + 1

            dominant_nature = sorted(nature_tags)[0]

            entry_rows.append(
                {
                    "feature_id": int(feature_id),
                    "count": int(count),
                    "collision": bool(collision),
                    "dominant_nature": dominant_nature,
                    "nature_tags": sorted(nature_tags),
                    "original_features": sorted(set(original_features)),
                }
            )
        entry_index[entry_name] = entry_rows
        try:
            _save_feature_shell_payload(
                shell_dir / f"{_safe_name(entry_name)}.pkl.gz",
                entry,
                feature_shells,
                params["pdb_dir"],
            )
        except Exception as ex:
            print(
                f"[luna-api] aviso: nao foi possivel salvar shells de {entry_name}: {ex}",
                flush=True,
            )

    values = [row["molecule_hits"] for row in feature_summary.values()]
    mean_hits = (sum(values) / len(values)) if values else 0.0
    std_hits = 0.0
    if values:
        variance = sum((value - mean_hits) ** 2 for value in values) / len(values)
        std_hits = variance ** 0.5

    features = []
    for feature_id in sorted(feature_summary):
        info = feature_summary[feature_id]
        coverage_pct = (100.0 * info["molecule_hits"] / total_entries) if total_entries else 0.0
        dominant_nature = "Unreliable feature"
        if info["nature_counts"]:
            dominant_nature = max(
                sorted(info["nature_counts"]),
                key=lambda key: (info["nature_counts"][key], key),
            )
        zscore = 0.0
        if std_hits > 1e-12:
            zscore = (info["molecule_hits"] - mean_hits) / std_hits
        features.append(
            {
                "feature_id": int(feature_id),
                "molecule_hits": int(info["molecule_hits"]),
                "coverage_pct": round(coverage_pct, 3),
                "cutoff_pct": round(coverage_pct, 3),
                "zscore": round(zscore, 6),
                "total_count": int(info["total_count"]),
                "collision_hits": int(info["collision_hits"]),
                "dominant_nature": dominant_nature,
                "nature_breakdown": info["nature_counts"],
            }
        )

    artifact = {
        "ifp_type": type_name,
        "ifp_label": IFP_LABELS.get(type_name, type_name),
        "fingerprint_length": int(proj.ifp_length),
        "count_fingerprint": bool(proj.ifp_count),
        "total_molecules": int(total_entries),
        "features": features,
        "entry_index": entry_index,
        "source_file": params.get("ifp_outputs", {}).get(type_name),
    }
    artifact_path = workdir / "results" / "fingerprints" / f"fp_analysis_{suffix}.json"
    _save_json(artifact_path, artifact)
    print(f"[luna-api] an\u00e1lise de fingerprints salva em {artifact_path}", flush=True)


def _generate_ifps(proj, params):
    ifp_types = list(params.get("ifp_types") or [])
    if not ifp_types:
        return

    outputs = params.get("ifp_outputs") or {}
    sim_outputs = params.get("sim_matrix_outputs") or {}
    wants_similarity = bool(params.get("sim_matrix", False))

    proj.calc_ifp = True
    proj.calc_mfp = False
    proj.ifp_num_levels = params.get("ifp_levels", 2)
    proj.ifp_radius_step = params.get("ifp_radius", 5.73171)
    proj.ifp_length = params.get("ifp_length", 4096)
    proj.ifp_count = not params.get("ifp_bit", False)

    total_ifp_types = max(1, len(ifp_types))
    for ifp_index, type_name in enumerate(ifp_types, start=1):
        proj.ifp_type = getattr(IFPType, type_name)
        proj.ifp_output = outputs.get(type_name)
        # Avoid LUNA's internal multiprocessing similarity calculation on Windows:
        # it can exhaust memory while pickling large jobs into worker queues.
        proj.ifp_sim_matrix_output = None
        start_pct = 76 + int(12 * (ifp_index - 1) / total_ifp_types)
        end_pct = 76 + int(12 * ifp_index / total_ifp_types)
        print(f"[luna-api-progress] {start_pct}% - gerando fingerprints {type_name}", flush=True)
        print(f"[luna-api] gerando fingerprints {type_name}...", flush=True)
        proj.generate_fps()
        if proj.ifp_output:
            print(f"[luna-api] IFP {type_name} salvo em {proj.ifp_output}", flush=True)
        print(f"[luna-api-progress] {end_pct}% - fingerprint {type_name} salvo", flush=True)
        sim_output_path = sim_outputs.get(type_name) if wants_similarity else None
        if sim_output_path:
            square_path = str(Path(sim_output_path).with_name(Path(sim_output_path).stem + "_square.csv"))
            saved_edge = None
            saved_square = None
            if proj.ifp_output and Path(proj.ifp_output).exists():
                saved_edge, saved_square = _write_similarity_outputs_from_ifp(
                    proj.ifp_output,
                    sim_output_path,
                    square_path,
                )
                if saved_edge:
                    print(f"[luna-api] Similaridade {type_name} reconstruida e salva em {saved_edge}", flush=True)
            if saved_square:
                print(f"[luna-api] Similaridade quadrada {type_name} salva em {saved_square}", flush=True)
        _export_fp_artifacts(proj, params, type_name)


AA_RESIDUES = {
    "ALA", "ARG", "ASN", "ASP", "ASH", "CYS", "CYM", "CYX", "GLN", "GLU", "GLH",
    "GLY", "HIS", "HID", "HIE", "HIP", "ILE", "LEU", "LYS", "LYN", "MET", "PHE",
    "PRO", "SER", "THR", "TRP", "TYR", "VAL",
}

WATER_RESIDUES = {"HOH", "WAT", "TIP", "SOL", "T3P", "H2O", "OH2", "DOD"}


def _interaction_type_name(interaction):
    value = getattr(interaction, "type", None)
    name = getattr(value, "name", None) or getattr(value, "value", None)
    if name:
        return str(name)
    return str(value or interaction.__class__.__name__)


def _residue_label(atom, include_water=False):
    residue = getattr(atom, "parent", None)
    if residue is None:
        return None
    resname = str(getattr(residue, "resname", "") or "").upper()
    if resname not in AA_RESIDUES and not (include_water and resname in WATER_RESIDUES):
        return None
    chain = getattr(getattr(residue, "parent", None), "id", "?")
    resid = getattr(residue, "id", None)
    if isinstance(resid, (tuple, list)) and len(resid) > 1:
        resid = resid[1]
    return f"{chain}/{resname}/{resid}"


def _labels_from_group(group, include_water=False):
    labels = set()
    for atom in getattr(group, "atoms", []) or []:
        label = _residue_label(atom, include_water=include_water)
        if label:
            labels.add(label)
    return labels


def _interaction_residue_labels(interaction):
    src_grp = getattr(interaction, "src_grp", None)
    trgt_grp = getattr(interaction, "trgt_grp", None)
    src_role = _group_role(src_grp)
    trgt_role = _group_role(trgt_grp)
    pair = {src_role, trgt_role}
    labels = set()
    if pair == {"ligand", "protein"}:
        labels.update(_labels_from_group(src_grp if src_role == "protein" else trgt_grp))
    elif pair == {"protein", "water"}:
        labels.update(_labels_from_group(src_grp if src_role == "protein" else trgt_grp))
    elif pair == {"ligand", "water"}:
        labels.update(_labels_from_group(src_grp if src_role == "water" else trgt_grp, include_water=True))
    return labels


def _save_json(path, data):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(data, fh, indent=2)


def _export_summary_artifacts(proj, workdir):
    summary = {
        "entries": 0,
        "interaction_counts": {},
        "entry_interaction_counts": {},
        "errors": [],
    }
    residue_counts = {}
    entries = []
    residues = set()
    interaction_types = set()

    for entry in list(getattr(proj, "entries", []) or []):
        entry_name = getattr(entry, "to_string", lambda: str(entry))()
        entries.append(entry_name)
        try:
            entry_results = proj.get_entry_results(entry)
            interactions = list(getattr(entry_results.interactions_mngr, "interactions", []) or [])
        except Exception as ex:
            summary["errors"].append(f"{entry_name}: {ex}")
            continue

        per_entry_counts = {}
        for interaction in interactions:
            itype = _interaction_type_name(interaction)
            interaction_types.add(itype)
            per_entry_counts[itype] = per_entry_counts.get(itype, 0) + 1
            summary["interaction_counts"][itype] = summary["interaction_counts"].get(itype, 0) + 1

            touched_residues = _interaction_residue_labels(interaction)
            for label in touched_residues:
                residues.add(label)
                by_type = residue_counts.setdefault(itype, {})
                by_entry = by_type.setdefault(entry_name, {})
                by_entry[label] = by_entry.get(label, 0) + 1

        if per_entry_counts:
            summary["entries"] += 1
        summary["entry_interaction_counts"][entry_name] = per_entry_counts

    residue_labels = sorted(residues)
    ordered_entries = list(entries)
    matrix = {}
    for itype in sorted(interaction_types):
        by_entry = residue_counts.get(itype, {})
        matrix[itype] = [
            [float(by_entry.get(entry_name, {}).get(label, 0.0)) for label in residue_labels]
            for entry_name in ordered_entries
        ]

    residue_artifact = {
        "interaction_types": sorted(interaction_types),
        "residues": residue_labels,
        "entries": ordered_entries,
        "matrix": matrix,
        "errors": list(summary["errors"]),
    }

    summary_path = os.path.join(workdir, "results", "analysis_summary.json")
    residue_path = os.path.join(workdir, "results", "residue_matrix.json")
    _save_json(summary_path, summary)
    _save_json(residue_path, residue_artifact)
    print(f"[luna-api] resumo salvo em {summary_path}", flush=True)
    print(f"[luna-api] matriz de resíduos salva em {residue_path}", flush=True)


workdir = p["workdir"]
pdb_dir = p["pdb_dir"]
lig_file = p["lig_file"]
lig_ext = Path(lig_file).suffix.lower()
lig_mol_obj_type = p.get("lig_mol_obj_type") or (
    "rdkit" if lig_ext in {".sdf", ".sd", ".mol"} else "openbabel"
)
fork_from = p.get("fork_from") or ""
add_h = bool(p.get("add_h", True))
amend_mol = bool(p.get("amend_mol", True))
ph = float(p.get("ph", 7.4))
entry_specs = p.get("entry_specs") or [
    {"pdb_id": p["pdb_id"], "ligand_name": name}
    for name in p.get("entries", [])
]

if add_h:
    print(f"[luna-api] add_h ativado em pH {ph}", flush=True)
else:
    print("[luna-api] add_h desativado", flush=True)
print(f"[luna-api] ligantes serao lidos com {lig_mol_obj_type}", flush=True)
if not amend_mol:
    print("[luna-api] PDB receptor sera preservado sem amend_mol", flush=True)

fork_copied = _prepare_fork_tree(fork_from, workdir)
if fork_copied:
    print(f"[luna-api] projeto fonte copiado de {fork_from} para {workdir}", flush=True)

# ----- Build MolFileEntry objects -----
entry_objs = []
for spec in entry_specs:
    name = spec["ligand_name"]
    pdb_id = spec["pdb_id"]
    try:
        e = MolFileEntry.from_mol_file(
            pdb_id=pdb_id,
            mol_id=name,
            mol_file=lig_file,
            is_multimol_file=True,
            mol_obj_type=lig_mol_obj_type,
        )
        entry_objs.append(e)
    except Exception as ex:
        print(f"[warn] pulando '{name}': {ex}", flush=True)
print(f"[luna-api] {len(entry_objs)} novas entries carregadas", flush=True)
print("[luna-api-progress] 10% - entradas carregadas", flush=True)

existing_entries = []
if fork_from:
    project_file = LocalProject.get_project_file(workdir)
    if not os.path.exists(project_file):
        raise OSError(f"Nenhum projeto LUNA válido foi encontrado em '{fork_from}'.")
    existing_entries = _load_existing_entries(workdir)
    print(f"[luna-api] fork com {len(existing_entries)} entries já existentes", flush=True)

project_entries = _merge_entries(existing_entries, entry_objs)
print(f"[luna-api] total de entries nesta execução: {len(project_entries)}", flush=True)
print("[luna-api-progress] 14% - projeto preparado", flush=True)

# ----- DefaultInteractionConfig / external .cfg with overrides -----
config_file = (p.get("interaction_config_file") or "").strip()
if config_file and os.path.exists(config_file):
    inter_config = InteractionConfig.from_config_file(config_file)
    print(f"[luna-api] inter_config carregado de {config_file}", flush=True)
else:
    inter_config = DefaultInteractionConfig()
    if config_file:
        print(f"[warn] arquivo de interações não encontrado: {config_file}", flush=True)
for k, v in (p.get("inter_config_overrides") or {}).items():
    try:
        inter_config[k] = v
        print(f"[luna-api] inter_config[{k}] = {v}", flush=True)
    except Exception as ex:
        print(f"[warn] inter_config[{k}] inválido: {ex}", flush=True)

inter_max_distance_cap = float(p.get("inter_max_distance_cap") or 0.0)
if inter_max_distance_cap > 0:
    for key in DISTANCE_CAP_KEYS:
        try:
            current = float(inter_config[key])
        except Exception:
            continue
        if current > inter_max_distance_cap:
            try:
                inter_config[key] = inter_max_distance_cap
                print(
                    f"[luna-api] inter_config[{key}] capped to {inter_max_distance_cap}",
                    flush=True,
                )
            except Exception as ex:
                print(f"[warn] cap global em {key} falhou: {ex}", flush=True)

# ----- InteractionFilter + Calculator -----
inter_filter = InteractionFilter.new_pli_filter(
    ignore_self_inter=p.get("ic_ignore_self_inter", True),
    ignore_any_h2o=not p.get("include_waters", False),
)
ic_kwargs = {
    "inter_filter": inter_filter,
    "inter_config": inter_config,
    "add_proximal": p.get("ic_add_proximal", False),
    "add_atom_atom": p.get("ic_add_atom_atom", False),
    "add_dependent_inter": p.get("ic_add_dependent_inter", True),
}
try:
    import inspect

    ic_sig = inspect.signature(InteractionCalculator.__init__)
    if "add_h2o_pairs_with_no_target" in ic_sig.parameters:
        ic_kwargs["add_h2o_pairs_with_no_target"] = p.get("ic_add_h2o_pairs_with_no_target", True)
except Exception:
    pass
inter_calc = InteractionCalculator(**ic_kwargs)

# ----- LocalProject -----
opts = dict(
    entries=project_entries,
    working_path=workdir,
    pdb_path=pdb_dir,
    overwrite_path=bool(p.get("overwrite", False) and not fork_from),
    add_h=add_h,
    ph=ph,
    amend_mol=amend_mol,
    calc_ifp=False,
    out_pse=False,
    nproc=p.get("nproc", 1),
    inter_calc=inter_calc,
    append_mode=bool(existing_entries),
)
proj = LocalProject(**opts)
print("[luna-api] iniciando proj.run()...", flush=True)
print("[luna-api-progress] 20% - calculando interacoes", flush=True)
proj.run()
print("[luna-api-progress] 70% - interacoes calculadas", flush=True)
_export_summary_artifacts(proj, workdir)
print("[luna-api-progress] 74% - resumo e matriz de residuos salvos", flush=True)

# ----- Export IFP CSV / similarity matrix -----
_generate_ifps(proj, p)

# ----- PSE sessions (with optional per-type filter) -----
if p.get("out_pse", False):
    print("[luna-api-progress] 92% - gerando sessoes PyMOL", flush=True)
    try:
        from luna.interaction.view import InteractionViewer
    except Exception as ex:
        print(f"[warn] InteractionViewer indisponível: {ex}", flush=True)
    else:
        pse_dir = p.get("pse_path") or os.path.join(workdir, "results", "pse")
        os.makedirs(pse_dir, exist_ok=True)
        types_filter = p.get("pse_interaction_types") or []
        for entry in proj.entries:
            er = proj.get_entry_results(entry)
            if types_filter:
                inter = list(er.interactions_mngr.filter_by_types(types_filter))
            else:
                inter = er.interactions_mngr
            safe_name = entry.to_string().replace(":", "_").replace("/", "_")
            pse_path = os.path.join(pse_dir, f"{safe_name}.pse")
            pdb_file = getattr(entry, "pdb_file", "")
            if not pdb_file:
                workdir_pdb = os.path.join(workdir, "pdbs", f"{entry.pdb_id}.pdb")
                pdb_file = workdir_pdb if os.path.exists(workdir_pdb) else os.path.join(pdb_dir, f"{entry.pdb_id}.pdb")
            try:
                viewer = InteractionViewer(show_hydrop_surface=False)
                viewer.new_session([(entry, inter, pdb_file)], pse_path)
            except Exception as ex:
                print(f"[warn] PSE para {safe_name} falhou: {ex}", flush=True)
        print(f"[luna-api] PSE salvos em {pse_dir}", flush=True)
        print("[luna-api-progress] 97% - sessoes PyMOL salvas", flush=True)

print("[luna-api] concluído", flush=True)
print("[luna-api-progress] 100% - concluido", flush=True)
'''


def _normalize_complex_name(name: str) -> str:
    stem = Path(name).stem.strip()
    lowered = stem.lower()
    for suffix in _LIGAND_SUFFIXES:
        if lowered.endswith(suffix):
            return stem[:-len(suffix)]
    return stem


def protein_has_explicit_hydrogens(protein_file: str | Path) -> bool:
    """Return True when the receptor PDB already contains explicit hydrogens."""
    path = Path(protein_file)
    if not path.exists():
        return False

    with path.open("r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            element = line[76:78].strip() if len(line) >= 78 else ""
            atom_name = line[12:16].strip() if len(line) >= 16 else ""
            if element.upper() == "H":
                return True
            if not element and atom_name.upper().startswith("H"):
                return True
    return False


def protein_is_gui_preprocessed(protein_file: str | Path) -> bool:
    """Return True when the PDB was generated by the docking pre-processor."""
    path = Path(protein_file)
    if not path.exists():
        return False

    with path.open("r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            stripped = line.strip()
            if not stripped:
                continue
            return line.startswith(_PREPARED_PROTEIN_MARKER)
    return False


def _candidate_protein_files(cfg: ProjectConfig) -> list[Path]:
    protein_path = Path(cfg.protein_file)
    if not protein_path.exists():
        return []
    if not cfg.include_waters:
        if protein_path.is_dir():
            return []
        return [protein_path]
    if protein_path.is_dir():
        return sorted(protein_path.glob("*.pdb"))
    return sorted(protein_path.parent.glob("*.pdb"))


def _iter_ligand_hydrogen_flags(ligand_file: str | Path) -> list[bool]:
    path = Path(ligand_file)
    if not path.exists():
        return []
    suffix = path.suffix.lower()
    if suffix != ".mol2":
        return []

    flags: list[bool] = []
    for block, _name in _split_mol2_blocks(path):
        flags.append(_mol2_block_has_hydrogens(block))
    return flags


def ligand_file_has_explicit_hydrogens(ligand_file: str | Path) -> bool:
    return any(_iter_ligand_hydrogen_flags(ligand_file))


def validate_hydrogen_inputs(cfg: ProjectConfig) -> list[str]:
    """Return user-facing warnings when Add_H is off but inputs are not protonated."""
    if cfg.add_h:
        return []

    messages: list[str] = []
    missing_protein_h = [
        path.name for path in _candidate_protein_files(cfg)
        if not protein_has_explicit_hydrogens(path)
    ]
    if missing_protein_h:
        if len(missing_protein_h) == 1:
            messages.append(
                "Add_H está desmarcado, mas a proteína "
                f"'{missing_protein_h[0]}' não contém hidrogênios explícitos."
            )
        else:
            messages.append(
                "Add_H está desmarcado, mas algumas proteínas não contêm hidrogênios explícitos: "
                + ", ".join(missing_protein_h[:5])
                + ("..." if len(missing_protein_h) > 5 else "")
            )

    if Path(cfg.ligand_file).suffix.lower() == ".mol2":
        ligand_flags = _iter_ligand_hydrogen_flags(cfg.ligand_file)
        if ligand_flags and not all(ligand_flags):
            messages.append(
                "Add_H está desmarcado, mas o arquivo de ligantes tem moléculas sem hidrogênios explícitos."
            )
        elif ligand_flags and not any(ligand_flags):
            messages.append(
                "Add_H está desmarcado, mas o arquivo de ligantes não contém hidrogênios explícitos."
            )

    if messages:
        messages.append(
            "Marque 'Adicionar hidrogênios antes da análise' na aba 3.Análises ou use arquivos já protonados."
        )
    return messages


def _split_mol2_blocks(path: Path) -> list[tuple[list[str], str]]:
    lines = path.read_text(errors="replace").splitlines(keepends=True)
    out: list[tuple[list[str], str]] = []
    cur: list[str] = []
    name = ""
    for line in lines:
        if line.startswith("@<TRIPOS>MOLECULE"):
            if cur:
                out.append((cur, name))
            cur = [line]
            name = ""
            continue
        cur.append(line)
        if name == "" and len(cur) == 2:
            name = line.strip()
    if cur:
        out.append((cur, name))
    return out


def _mol2_block_has_hydrogens(block: list[str]) -> bool:
    in_atoms = False
    for line in block:
        if line.startswith("@<TRIPOS>ATOM"):
            in_atoms = True
            continue
        if line.startswith("@<TRIPOS>") and in_atoms:
            break
        if not in_atoms:
            continue
        parts = line.split()
        if len(parts) < 6:
            continue
        atom_name = parts[1].strip().upper()
        atom_type = parts[5].strip().upper()
        if atom_type.startswith("H") or atom_name.lstrip("0123456789").startswith("H"):
            return True
    return False


def resolve_protein_processing_flags(cfg: ProjectConfig) -> dict[str, bool | float]:
    """Decide how the GUI should run LUNA without rewriting user inputs."""
    protein_files = _candidate_protein_files(cfg)
    has_h = any(protein_has_explicit_hydrogens(path) for path in protein_files)
    is_preprocessed = any(protein_is_gui_preprocessed(path) for path in protein_files)
    add_h = bool(cfg.add_h)
    ligand_has_h = ligand_file_has_explicit_hydrogens(cfg.ligand_file)
    return {
        "protein_has_hydrogens": has_h,
        "protein_is_preprocessed": is_preprocessed,
        "ligand_has_hydrogens": ligand_has_h,
        "add_h": add_h,
        "ph": float(cfg.ph),
        "amend_mol": not (has_h or is_preprocessed),
        "stage_protein_without_h": False,
        "stage_ligand_without_h": False,
    }


def should_use_api_runner(cfg: ProjectConfig) -> bool:
    """Return True when we need the Python API runner instead of the CLI."""
    flags = resolve_protein_processing_flags(cfg)
    return cfg.uses_python_api() or not flags["amend_mol"]


def _protein_index(protein_dir: Path) -> tuple[dict[str, str], set[str]]:
    index: dict[str, str] = {}
    duplicates: set[str] = set()
    for pdb_file in protein_dir.glob("*.pdb"):
        key = _normalize_complex_name(pdb_file.stem)
        if key in index and index[key] != pdb_file.stem:
            duplicates.add(key)
            continue
        index[key] = pdb_file.stem
    return index, duplicates


def build_entry_specs(cfg: ProjectConfig, entries: list[str]) -> list[dict[str, str]]:
    """Return the per-entry protein/ligand mapping for the API runner."""
    protein_path = Path(cfg.protein_file)
    if not cfg.include_waters:
        pdb_id = protein_path.stem
        return [{"pdb_id": pdb_id, "ligand_name": ligand_name} for ligand_name in entries]

    protein_dir = protein_path if protein_path.is_dir() else protein_path.parent
    protein_index, duplicates = _protein_index(protein_dir)

    errors: list[str] = []
    specs: list[dict[str, str]] = []
    for ligand_name in entries:
        key = _normalize_complex_name(ligand_name)
        if key in duplicates:
            errors.append(
                f"Mais de um PDB corresponde ao ligante '{ligand_name}' na pasta {protein_dir}."
            )
            continue
        pdb_id = protein_index.get(key)
        if pdb_id is None:
            errors.append(
                f"Nenhum PDB com o mesmo nome do ligante '{ligand_name}' foi encontrado em {protein_dir}."
            )
            continue
        specs.append({"pdb_id": pdb_id, "ligand_name": ligand_name})

    if errors:
        raise ValueError("\n".join(errors))
    return specs


def validate_entry_specs(cfg: ProjectConfig, entries: list[str]) -> list[str]:
    """Return user-facing validation errors for per-complex execution."""
    try:
        build_entry_specs(cfg, entries)
    except ValueError as exc:
        return [line for line in str(exc).splitlines() if line.strip()]
    return []


def write_runner(workdir: str | Path) -> str:
    """Write the API runner script to the workdir and return its path."""
    wd = Path(workdir)
    wd.mkdir(parents=True, exist_ok=True)
    path = wd / "_luna_api_runner.py"
    path.write_text(API_RUNNER_SCRIPT, encoding="utf-8")
    return str(path)


def write_params(workdir: str | Path, cfg: ProjectConfig, entries: list[str]) -> str:
    """Dump the relevant ProjectConfig subset as JSON next to the runner."""
    wd = Path(workdir)
    wd.mkdir(parents=True, exist_ok=True)
    entry_specs = build_entry_specs(cfg, entries)
    protein_flags = resolve_protein_processing_flags(cfg)
    protein_path = Path(cfg.protein_file)
    protein_dir = protein_path if protein_path.is_dir() else protein_path.parent
    protein_files = sorted(protein_dir.glob("*.pdb")) if protein_path.is_dir() else [protein_path]
    pdb_id = protein_path.stem
    if protein_path.is_dir() and protein_files:
        pdb_id = protein_files[0].stem
    params = {
        "workdir": cfg.workdir,
        "pdb_id": pdb_id,
        "pdb_dir": str(protein_dir),
        "lig_file": cfg.ligand_file,
        "lig_mol_obj_type": ligand_mol_obj_type(cfg.ligand_file),
        "entries": entries,
        "entry_specs": entry_specs,
        "fork_from": cfg.fork_from,
        "overwrite": cfg.overwrite,
        "nproc": max(1, cfg.nproc or 1),
        "add_h": protein_flags["add_h"],
        "ph": protein_flags["ph"],
        "amend_mol": protein_flags["amend_mol"],
        "protein_has_hydrogens": protein_flags["protein_has_hydrogens"],
        "protein_is_preprocessed": protein_flags["protein_is_preprocessed"],
        "stage_protein_without_h": protein_flags["stage_protein_without_h"],
        "stage_ligand_without_h": protein_flags["stage_ligand_without_h"],
        "include_waters": cfg.include_waters,
        "out_ifp": cfg.out_ifp,
        "ifp_type": cfg.ifp_type,
        "ifp_types": cfg.selected_ifp_types() if (cfg.out_ifp or cfg.sim_matrix) else [],
        "ifp_levels": cfg.ifp_levels,
        "ifp_radius": cfg.ifp_radius,
        "ifp_length": cfg.ifp_length,
        "ifp_bit": cfg.ifp_bit,
        "ifp_output": cfg.ifp_output,
        "ifp_outputs": resolve_ifp_output_paths(cfg) if (cfg.out_ifp or cfg.sim_matrix) else {},
        "sim_matrix": cfg.sim_matrix,
        "sim_matrix_output": cfg.sim_matrix_output,
        "sim_matrix_outputs": resolve_sim_matrix_output_paths(cfg),
        "out_pse": cfg.out_pse,
        "pse_path": cfg.pse_path,
        "pse_interaction_types": cfg.pse_interaction_types,
        "ic_add_proximal": cfg.ic_add_proximal,
        "ic_add_atom_atom": cfg.ic_add_atom_atom,
        "ic_add_dependent_inter": cfg.ic_add_dependent_inter,
        "ic_add_h2o_pairs_with_no_target": cfg.ic_add_h2o_pairs_with_no_target,
        "ic_ignore_self_inter": cfg.ic_ignore_self_inter,
        "interaction_config_file": cfg.interaction_config_file,
        "inter_max_distance_cap": cfg.inter_max_distance_cap,
        "inter_config_overrides": cfg.inter_config_overrides,
    }
    path = wd / "_luna_api_params.json"
    path.write_text(json.dumps(params, indent=2), encoding="utf-8")
    return str(path)


def ligand_mol_obj_type(ligand_file: str | Path) -> str:
    """Choose the LUNA molecule backend that can parse the ligand file."""
    suffix = Path(ligand_file).suffix.lower()
    if suffix in {".sdf", ".sd", ".mol"}:
        return "rdkit"
    return "openbabel"


def build_api_command(py_exe: str, cfg: ProjectConfig, entries: list[str]) -> list[str]:
    """Return argv to run the Python-API runner."""
    runner = write_runner(cfg.workdir)
    params = write_params(cfg.workdir, cfg, entries)
    return [py_exe, runner, params]
