"""Run LUNA via its Python API (inside luna-env) instead of the CLI.

Used when the user enables any advanced option (water handling, custom
InteractionCalculator flags, DefaultInteractionConfig overrides, per-type
PSE filtering, project forks, or multi-IFP export). Writes a self-contained
Python script to the workdir and returns the argv needed to run it.
"""
from __future__ import annotations

import json
import re
from pathlib import Path

from .luna_runner import safe_nproc
from .project import ProjectConfig, resolve_ifp_output_paths, resolve_sim_matrix_output_paths
from .results_analysis import INTERACTION_COLORS


_LIGAND_SUFFIXES = ("_ligand", "-ligand", "_lig", "-lig")
_LIGAND_FILE_SUFFIXES = {".mol2", ".sdf", ".sd", ".mol", ".pdb", ".ent"}
_PREPARED_PROTEIN_MARKER = "REMARK   Separated Protein"
_WATER_RESIDUE_NAMES = {"HOH", "WAT", "WTM", "TIP", "SOL", "T3P", "H2O", "OH2", "DOD"}


API_RUNNER_SCRIPT = r'''
"""Auto-generated LUNA runner. Reads params from argv[1] (JSON file)."""
import copy, csv, gzip, json, os, pickle, re, shutil, sys
from pathlib import Path

params_file = sys.argv[1]
with open(params_file, "r", encoding="utf-8") as fh:
    p = json.load(fh)

# A receptor PDB can intentionally contain non-water HETATM residues such as
# cofactors and coordinated metal ions.  The GUI keeps this opt-in so legacy
# projects retain their original interpretation.
INCLUDE_PROTEIN_HETEROATOMS = bool(p.get("include_protein_heteroatoms", False))
PROTEIN_HETEROATOM_RESIDUES = {
    str(value).strip()
    for value in (p.get("protein_heteroatom_residues") or [])
    if str(value).strip()
}

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


def _pymol_color_name(interaction_name):
    return "hip2l_" + re.sub(r"[^a-z0-9]+", "_", str(interaction_name).lower()).strip("_")


def _pymol_interaction_tokens(interaction_name):
    normalized = re.sub(r"[^a-z0-9]+", "", str(interaction_name).lower())
    return (normalized,)


def _pymol_interaction_markers(interaction_name):
    # Exact LUNA object-name segments for a displayed interaction.
    normalized = re.sub(r"[^a-z0-9]+", "", str(interaction_name).lower())
    short_names = {
        "hydrogenbond": ("hbond",),
        "weakhydrogenbond": ("weak_hbond",),
        "waterbridgedhydrogenbond": ("water_hbond",),
        "halogenbond": ("xbond",),
        "halogenpi": ("x-pi",),
        "chalcogenbond": ("ybond",),
        "chalcogenpi": ("y-pi",),
        "ionic": ("ionic",),
        "saltbridge": ("salt_bridge",),
        "cationpi": ("cation-pi",),
        "cationnucleophile": ("cat_nucleop",),
        "anionelectrophile": ("ani_electrop",),
        "pistacking": ("pi-stack",),
        "aromaticstacking": ("pi-stack", "aromatic-bond"),
        "edgetoface": ("edge-to-face_stack",),
        "facetoface": ("face-to-face_stack",),
        "facetoedgepistacking": ("face-to-edge_stack",),
        "facetofacepistacking": ("face-to-face_stack",),
        "facetoslopepistacking": ("face-to-slope_stack",),
        "displacedfacetoedgepistacking": ("disp_face-to-edge_stack",),
        "displacedfacetofacepistacking": ("disp_face-to-face_stack",),
        "displacedfacetoslopepistacking": ("disp_face-to-slope_stack",),
        "parallel": ("par_multipol",),
        "parallelmultipolar": ("par_multipol",),
        "antiparallelmultipolar": ("antipar_multipol",),
        "orthogonalmultipolar": ("ort_multipol",),
        "tiltedmultipolar": ("tilted_multipol",),
        "hydrophobic": ("hphobe",),
        "amidearomaticstacking": ("amide_stack",),
        "metalcoordination": ("metal-coord",),
        "vanderwaals": ("vdw",),
        "proximal": ("prox",),
        "multipolarinteraction": ("multipol",),
        "repulsive": ("repuls",),
        "unfavorableanionnucleophile": ("unf_ani_nucleop",),
        "unfavorablecationelectrophile": ("unf_cat_electrop",),
        "unfavorableelectrophileelectrophile": ("unf_electrop_electrop",),
        "unfavorablenucleophilenucleophile": ("unf_nucleop_nucleop",),
    }
    return tuple(f".{name}." for name in short_names.get(normalized, ()))


def _matches_pymol_interaction_object(object_name, normalized_name, interaction_name):
    markers = _pymol_interaction_markers(interaction_name)
    if markers:
        object_name = str(object_name).lower()
        return any(marker in object_name for marker in markers)
    return any(token and token in normalized_name for token in _pymol_interaction_tokens(interaction_name))


def _apply_pse_interaction_colors(pse_path, palette, *, save=True):
    """Apply the GUI's canonical palette to the active PyMOL session."""
    if not isinstance(palette, dict) or not palette:
        return 0
    try:
        from pymol import cmd
    except Exception as exc:
        print(f"[warn] paleta PSE ignorada: pymol.cmd indisponivel ({ex})", flush=True)
        return 0
    try:
        object_names = list(cmd.get_names("objects") or [])
    except Exception as exc:
        print(f"[warn] paleta PSE ignorada: nao foi possivel listar objetos ({ex})", flush=True)
        return 0

    normalized_names = {
        name: re.sub(r"[^a-z0-9]+", "", str(name).lower())
        for name in object_names
    }
    colored = 0
    for interaction_name, hex_color in palette.items():
        color = str(hex_color or "").lstrip("#")
        if len(color) != 6:
            continue
        try:
            rgb = [int(color[index:index + 2], 16) / 255.0 for index in (0, 2, 4)]
        except ValueError:
            continue
        color_name = _pymol_color_name(interaction_name)
        try:
            cmd.set_color(color_name, rgb)
        except Exception:
            continue
        for object_name, normalized_name in normalized_names.items():
            if not _matches_pymol_interaction_object(
                object_name,
                normalized_name,
                interaction_name,
            ):
                continue
            try:
                cmd.color(color_name, object_name)
                cmd.set("dash_color", color_name, object_name)
                cmd.set("label_color", color_name, object_name)
                colored += 1
            except Exception:
                continue
    if colored and save:
        try:
            cmd.save(pse_path)
        except Exception as exc:
            print(f"[warn] PSE colorido, mas nao foi salvo: {exc}", flush=True)
            return 0
    return colored


def _save_viewer_session_with_palette(viewer, session_rows, pse_path, palette):
    """Apply colors immediately before InteractionViewer writes the PSE.

    LUNA's InteractionViewer builds all dash and arrow objects in PyMOL and
    calls ``cmd.save`` internally.  Wrapping that single save call avoids a
    second load/save cycle for every generated session.
    """
    if not isinstance(palette, dict) or not palette:
        viewer.new_session(session_rows, pse_path)
        return 0
    try:
        from pymol import cmd
        original_save = cmd.save
    except Exception as exc:
        raise RuntimeError(f"pymol.cmd indisponivel para salvar a paleta: {exc}") from exc

    colored = [0]

    def _save_with_palette(*args, **kwargs):
        colored[0] = _apply_pse_interaction_colors(pse_path, palette, save=False)
        return original_save(*args, **kwargs)

    try:
        cmd.save = _save_with_palette
        viewer.new_session(session_rows, pse_path)
    finally:
        cmd.save = original_save
    return colored[0]


def _reset_pymol_session():
    """Discard the previous entry before creating the next standalone PSE."""
    try:
        from pymol import cmd
        cmd.reinitialize()
        return True
    except Exception:
        return False


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
LIGAND_STRUCTURE_SUFFIXES = {".mol2", ".sdf", ".sd", ".mol", ".pdb", ".ent"}


def _ifp_seed(params):
    try:
        return int(params.get("ifp_seed", 0) or 0)
    except Exception:
        return 0


def _write_ifp_seed(params, type_name):
    workdir = Path(params["workdir"])
    suffix = IFP_SUFFIX.get(type_name, type_name)
    seed = _ifp_seed(params)
    seed_dir = workdir / "results" / "fingerprints"
    seed_dir.mkdir(parents=True, exist_ok=True)
    seed_path = seed_dir / f"seed_ifp_{suffix}_importance.txt"
    seed_path.write_text(f"{seed}\n", encoding="utf-8")
    source = str(params.get("ifp_seed_file") or "").strip()
    if source:
        try:
            source_path = Path(source)
            if source_path.exists():
                source_path.write_text(f"{seed}\n", encoding="utf-8")
        except Exception as ex:
            print(f"[luna-api] aviso: nao foi possivel reescrever seed informado: {ex}", flush=True)
    return str(seed_path), seed


def _safe_name(value):
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value)).strip("_") or "entry"


def _shell_level_key(shell):
    value = getattr(shell, "level", None)
    try:
        return str(int(value))
    except Exception:
        return str(value if value is not None else "unknown")


def _sorted_level_keys(values):
    keys = {str(value) for value in (values or []) if str(value).strip()}
    return sorted(
        keys,
        key=lambda value: (
            not str(value).lstrip("-").isdigit(),
            int(value) if str(value).lstrip("-").isdigit() else str(value),
        ),
    )


def _feature_key(feature_id, level):
    text = str(level or "").strip()
    return f"{int(feature_id)}_{text}" if text else ""


def _legacy_level_breakdown_from_nature(nature_breakdown):
    level_map = {
        "Ligand's level 0 features only": "0",
        "Protein's level 0 features only": "0",
        "Intraligand interactions only": "0",
        "Intraprotein interactions only": "0",
        "Has noncovalent interactions with the protein": "0",
        "Upper level with ligand atomic information only": "1",
        "Upper level with protein atomic information only": "1",
    }
    inferred = {}
    unknown_count = 0
    for raw_name, raw_count in (nature_breakdown or {}).items():
        try:
            count = int(raw_count)
        except Exception:
            continue
        if count <= 0:
            continue
        name = str(raw_name or "").strip()
        level = level_map.get(name)
        if level is None:
            match = re.search(r"\blevel\s*(-?\d+)\b", name, flags=re.IGNORECASE)
            if match:
                level = match.group(1)
        if level is None:
            unknown_count += count
            continue
        inferred[level] = inferred.get(level, 0) + count
    if unknown_count > 0:
        return {}
    return inferred


def _level_threshold(values, use_otsu=False):
    positives = [float(value) for value in values if float(value) > 0.0]
    if not positives:
        return 100.0, "no_positive_values", [0.0 for _value in values]
    mean_value = sum(positives) / len(positives)
    variance = sum((value - mean_value) ** 2 for value in positives) / len(positives)
    std_value = variance ** 0.5
    zscores = [
        ((float(value) - mean_value) / std_value if float(value) > 0.0 and std_value > 1e-12 else 0.0)
        for value in values
    ]
    candidates = [float(value) for value, zscore in zip(values, zscores) if zscore > 1.0]
    if candidates:
        return min(candidates), "zscore_gt_1", zscores
    if not use_otsu:
        return 100.0, "no_reference", zscores
    unique = sorted(set(positives))
    if len(unique) == 1:
        return unique[0], "otsu_single_value", zscores
    best_threshold = (unique[0] + unique[1]) / 2.0
    best_score = -1.0
    for idx in range(len(unique) - 1):
        threshold = (unique[idx] + unique[idx + 1]) / 2.0
        low = [value for value in positives if value <= threshold]
        high = [value for value in positives if value > threshold]
        if not low or not high:
            continue
        score = (len(low) / len(positives)) * (len(high) / len(positives)) * ((sum(low) / len(low)) - (sum(high) / len(high))) ** 2
        if score > best_score:
            best_score = score
            best_threshold = threshold
    return best_threshold, "otsu", zscores


def _assign_feature_levels(features, use_otsu=False):
    collision_candidates = []
    for feature in features:
        breakdown = feature.get("shell_level_breakdown") or {}
        breakdown = {str(key): int(value) for key, value in breakdown.items() if int(value) > 0}
        feature["assigned_level"] = ""
        feature["assigned_level_pct"] = 0.0
        feature["assigned_level_source"] = "undetermined"
        feature["feature_key"] = ""
        if not breakdown:
            breakdown = _legacy_level_breakdown_from_nature(feature.get("nature_breakdown") or {})
        if not breakdown:
            continue
        level, count = max(sorted(breakdown.items()), key=lambda item: (int(item[1]), str(item[0])))
        total = sum(breakdown.values())
        pct = (100.0 * count / total) if total else 0.0
        feature["top_level"] = str(level)
        feature["top_level_pct"] = float(pct)
        if len(breakdown) == 1:
            feature["assigned_level"] = str(level)
            feature["assigned_level_pct"] = 100.0
            feature["assigned_level_source"] = "single_level"
            feature["feature_key"] = _feature_key(feature["feature_id"], level)
        else:
            collision_candidates.append(feature)
    threshold, source, zscores = _level_threshold([feature.get("top_level_pct", 0.0) for feature in collision_candidates], use_otsu)
    for feature, zscore in zip(collision_candidates, zscores):
        feature["assigned_level_zscore"] = float(zscore)
        if source in {"zscore_gt_1", "otsu", "otsu_single_value"} and float(feature.get("top_level_pct", 0.0)) >= threshold:
            level = str(feature.get("top_level", ""))
            feature["assigned_level"] = level
            feature["assigned_level_pct"] = float(feature.get("top_level_pct", 0.0))
            feature["assigned_level_source"] = source
            feature["feature_key"] = _feature_key(feature["feature_id"], level)
    return {"threshold_pct": threshold, "threshold_source": source}


def _rewrite_ifp_csv_with_feature_keys(path, features):
    if not path:
        return {"rewritten": False}
    path = Path(path)
    if not path.exists():
        return {"rewritten": False}
    mapping = {
        str(int(feature.get("feature_id", 0))): str(feature.get("feature_key", ""))
        for feature in features
        if str(feature.get("feature_key", "")).strip()
    }
    if not mapping:
        return {"rewritten": False}
    rows = []
    changed = False
    with path.open("r", encoding="utf-8", errors="replace", newline="") as fh:
        reader = csv.DictReader(fh)
        fieldnames = list(reader.fieldnames or [])
        if "count" not in fieldnames:
            fieldnames.append("count")
        for record in reader:
            bits = [token.strip() for token in str(record.get("on_bits") or "").split("\t") if token.strip()]
            raw_counts = [token.strip() for token in str(record.get("count") or "").split("\t") if token.strip()]
            counts = [float(token) for token in raw_counts] if len(raw_counts) == len(bits) else [1.0] * len(bits)
            row_counts = {}
            for bit, count in zip(bits, counts):
                new_bit = bit if "_" in bit else mapping.get(str(int(bit)) if str(bit).lstrip("-").isdigit() else bit, "")
                if not new_bit:
                    changed = True
                    continue
                if new_bit != bit:
                    changed = True
                row_counts[new_bit] = row_counts.get(new_bit, 0.0) + float(count)
            ordered = sorted(row_counts, key=lambda value: (int(str(value).split("_", 1)[0]), str(value)))
            updated = dict(record)
            updated["on_bits"] = "\t".join(ordered)
            updated["count"] = "\t".join(str(int(row_counts[bit])) if float(row_counts[bit]).is_integer() else f"{row_counts[bit]:.10g}" for bit in ordered)
            rows.append(updated)
    if not changed:
        return {"rewritten": False}
    backup = path.with_suffix(path.suffix + ".pre_level_assignment.bak")
    if not backup.exists():
        shutil.copy2(path, backup)
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    return {"rewritten": True, "backup": str(backup)}


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


def _residue_name_from_atom(atom):
    residue = getattr(atom, "parent", None)
    if residue is None:
        return ""
    resname = getattr(residue, "resname", None)
    if resname is None:
        try:
            resname = residue.get_resname()
        except Exception:
            resname = ""
    return str(resname or "").strip().upper()


def _residue_key_from_atom(atom):
    residue = getattr(atom, "parent", None)
    if residue is None:
        return ""
    resname = _residue_name_from_atom(atom)
    if not resname:
        return ""
    chain = getattr(getattr(residue, "parent", None), "id", "?")
    resid = getattr(residue, "id", None)
    if isinstance(resid, (tuple, list)) and len(resid) > 1:
        resid = resid[1]
    if resid in (None, ""):
        return ""
    return f"{chain}/{resname}/{resid}"


def _group_has_ligand_residue(atm_grp):
    for atom in getattr(atm_grp, "atoms", []) or []:
        resname = _residue_name_from_atom(atom)
        if resname and resname not in AA_RESIDUES and resname not in WATER_RESIDUES:
            return True
    return False


def _group_is_configured_protein_heteroatom(atm_grp):
    return any(
        _residue_key_from_atom(atom) in PROTEIN_HETEROATOM_RESIDUES
        for atom in (getattr(atm_grp, "atoms", []) or [])
    )


def _group_has_structural_residue(atm_grp):
    try:
        return bool(atm_grp.has_residue() or atm_grp.has_nucleotide())
    except Exception:
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
    if INCLUDE_PROTEIN_HETEROATOMS and _group_is_configured_protein_heteroatom(atm_grp):
        return "protein"
    if _group_has_structural_residue(atm_grp):
        try:
            is_hetatm = bool(atm_grp.has_hetatm())
        except Exception:
            is_hetatm = False
        # Keep amino-acid variants that happen to be encoded as HETATM (for
        # example MSE) on their historical protein path.  The opt-in controls
        # only non-standard residues such as cofactors and metal ions.
        if (
            not is_hetatm
            or not _group_has_ligand_residue(atm_grp)
            or INCLUDE_PROTEIN_HETEROATOMS
        ):
            return "protein"
    try:
        if atm_grp.has_hetatm():
            return "ligand"
    except Exception:
        pass
    if _group_has_ligand_residue(atm_grp):
        return "ligand"
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
                key = int(bit) if str(bit).lstrip("-").isdigit() else str(bit)
                row_map[key] = row_map.get(key, 0.0) + float(count)
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


def _portable_atom_record(atom):
    """Convert an ExtendedAtom/BioPython atom to JSON-safe identity data."""
    try:
        full_id = tuple(atom.get_full_id())
    except Exception:
        full_id = tuple(getattr(atom, "full_id", ()) or ())
    if len(full_id) < 5:
        return None
    residue_id = full_id[3]
    if isinstance(residue_id, (list, tuple)):
        residue_id = list(residue_id[:3])
    else:
        residue_id = [" ", residue_id, " "]
    atom_id = full_id[4]
    atom_name = atom_id[0] if isinstance(atom_id, (list, tuple)) else atom_id
    altloc = atom_id[1] if isinstance(atom_id, (list, tuple)) and len(atom_id) > 1 else " "
    residue = getattr(atom, "parent", None)
    resname = str(getattr(residue, "resname", "") or "").strip()
    chain = str(full_id[2])
    hetflag = str(residue_id[0])
    is_ligand = bool(getattr(atom, "is_ligand", False))
    if not is_ligand:
        is_ligand = chain.lower() == "z" or resname.upper() == "LIG" or residue_id[0] == "H_LIG"
    if is_ligand:
        object_role = "ligand"
    elif resname.upper() in {"HOH", "WAT", "WTM", "TIP", "SOL", "T3P", "H2O", "OH2", "DOD"}:
        object_role = "water"
    elif str(getattr(atom, "element", "") or "").strip().upper() in {"ZN", "MG", "CA", "FE", "MN", "CU", "CO", "NI", "NA", "K", "CL"}:
        object_role = "metal"
    elif hetflag.strip() and not hetflag.startswith(" "):
        object_role = "heteroatom"
    else:
        object_role = "protein"
    try:
        coord = [float(value) for value in atom.get_coord()]
    except Exception:
        coord = None
    return {
        "model": full_id[1],
        "chain": chain,
        "resname": resname,
        "hetflag": hetflag,
        "resseq": residue_id[1],
        "icode": str(residue_id[2]),
        "name": str(atom_name).strip(),
        "altloc": str(altloc),
        "element": str(getattr(atom, "element", "") or "").strip(),
        "coord": coord,
        "object_role": object_role,
        "is_ligand": is_ligand,
        "is_hydrogen": str(atom_name).strip().upper() in {"H", "D"}
        or str(getattr(atom, "element", "") or "").strip().upper() in {"H", "D"},
    }


def _portable_group_atoms(group):
    records = []
    for atom in list(getattr(group, "atoms", []) or []):
        record = _portable_atom_record(atom)
        if record is not None:
            records.append(record)
    return records


def _portable_shell_record(shell):
    def _bool_method(obj, name):
        try:
            value = getattr(obj, name, False)
            return bool(value() if callable(value) else value)
        except Exception:
            return False

    central_group = getattr(shell, "central_atm_grp", None)
    central_centroid = None
    try:
        central_centroid = [float(value) for value in list(central_group.centroid)[:3]]
    except Exception:
        pass
    interactions = []
    for interaction in list(getattr(shell, "interactions", []) or []):
        src_centroid = None
        trgt_centroid = None
        for attr, target in (("src_centroid", "src_centroid"), ("trgt_centroid", "trgt_centroid")):
            try:
                value = [float(item) for item in list(getattr(interaction, attr))[:3]]
            except Exception:
                value = None
            if target == "src_centroid":
                src_centroid = value
            else:
                trgt_centroid = value
        interactions.append(
            {
                "type": str(getattr(interaction, "type", type(interaction).__name__)),
                "src_group_atoms": _portable_group_atoms(getattr(interaction, "src_grp", None)),
                "trgt_group_atoms": _portable_group_atoms(getattr(interaction, "trgt_grp", None)),
                "src_interacting_atoms": _portable_group_atoms(getattr(interaction, "src_interacting_atms", None)),
                "trgt_interacting_atoms": _portable_group_atoms(getattr(interaction, "trgt_interacting_atms", None)),
                "src_centroid": src_centroid,
                "trgt_centroid": trgt_centroid,
                "directional": _bool_method(interaction, "is_directional"),
                "intramolecular": _bool_method(interaction, "is_intramol_interaction"),
                "unfavorable": "unfavorable" in str(getattr(interaction, "type", "")).lower(),
            }
        )
    return {
        "identifier": getattr(shell, "identifier", None),
        "level": getattr(shell, "level", None),
        "radius": getattr(shell, "radius", None),
        "valid": bool(getattr(shell, "valid", True)),
        "central_atoms": _portable_group_atoms(central_group),
        "central_centroid": central_centroid,
        "neighborhood_atoms": [
            atom
            for group in list(getattr(shell, "neighborhood", []) or [])
            for atom in _portable_group_atoms(group)
        ],
        "interactions": interactions,
    }


def _save_portable_shell_payload(path, entry, feature_shells, pdb_dir, proj, type_name):
    payload = {
        "schema": "hip2l.fp-shells",
        "schema_version": 2,
        "entry_meta": _entry_meta(entry),
        "feature_shells": {
            str(feature_id): [_portable_shell_record(shell) for shell in shells]
            for feature_id, shells in feature_shells.items()
        },
        "pdb_dir": pdb_dir,
        "ifp": {
            "type": str(type_name),
            "num_levels": int(getattr(proj, "ifp_num_levels", 0) or 0),
            "radius_step": float(getattr(proj, "ifp_radius_step", 0) or 0),
            "length": int(getattr(proj, "ifp_length", 0) or 0),
            "count": bool(getattr(proj, "ifp_count", False)),
            "diff_comp_classes": bool(getattr(proj, "ifp_diff_comp_classes", True)),
        },
    }
    temporary = path.with_name(f".{path.name}.part")
    try:
        with gzip.open(temporary, "wt", encoding="utf-8") as fh:
            json.dump(payload, fh, ensure_ascii=False, separators=(",", ":"))
        os.replace(temporary, path)
    except Exception:
        try:
            temporary.unlink()
        except OSError:
            pass
        raise


def _save_feature_shell_payload(path, entry, feature_shells, pdb_dir):
    payload = {
        "entry_meta": _entry_meta(entry),
        "feature_shells": feature_shells,
        "pdb_dir": pdb_dir,
    }
    with gzip.open(path, "wb") as fh:
        pickle.dump(payload, fh, protocol=pickle.HIGHEST_PROTOCOL)


def _fp_detail_interaction_name(interaction):
    value = getattr(interaction, "type", None)
    if value is None:
        return type(interaction).__name__
    return str(value).strip() or type(interaction).__name__


def _fp_detail_residue_label_from_group(atm_grp):
    labels = {}
    for atom in getattr(atm_grp, "atoms", []) or []:
        residue = getattr(atom, "parent", None)
        if residue is None:
            continue
        chain = "?"
        full_id = getattr(residue, "full_id", None)
        if isinstance(full_id, (list, tuple)) and len(full_id) >= 3:
            chain = str(full_id[2] or "?").strip() or "?"
        else:
            parent = getattr(residue, "parent", None)
            chain = str(getattr(parent, "id", "?") or "?").strip() or "?"
        resname = str(getattr(residue, "resname", None) or getattr(residue, "get_resname", lambda: "")()).strip() or "UNK"
        resid = getattr(residue, "id", None)
        seq = ""
        if isinstance(resid, (list, tuple)) and len(resid) >= 3:
            seq = f"{resid[1]}{str(resid[2] or '').strip()}".strip()
        elif resid is not None:
            seq = str(resid).strip()
        if not seq:
            continue
        label = f"{chain}/{resname}/{seq}"
        labels[label] = labels.get(label, 0) + 1
    if not labels:
        return ""
    return max(sorted(labels), key=lambda key: labels[key])


def _fp_detail_interaction_residue_group(src_grp, trgt_grp, src_role, trgt_role):
    pair = {src_role, trgt_role}
    if pair == {"ligand", "protein"}:
        return src_grp if src_role == "protein" else trgt_grp
    if pair == {"protein", "water"}:
        return src_grp if src_role == "protein" else trgt_grp
    if pair == {"ligand", "water"}:
        return src_grp if src_role == "water" else trgt_grp
    return None


def _increment_fp_detail_counter(container, key, amount=1):
    text = str(key or "").strip()
    if not text:
        return
    container[text] = container.get(text, 0) + int(amount)


def _add_fp_detail_shell(feature_details, feature_id, entry_name, shell):
    feature_key = str(int(feature_id))
    detail = feature_details.setdefault(
        feature_key,
        {
            "interaction_counts": {},
            "residue_counts": {},
            "pair_counts": {},
            "shell_level_counts": {},
            "entries": {},
        },
    )
    entry_info = detail["entries"].setdefault(
        entry_name,
        {
            "shell_count": 0,
            "interaction_counts": {},
            "residue_counts": {},
            "pair_counts": {},
            "shell_level_counts": {},
        },
    )
    entry_info["shell_count"] += 1
    level_key = _shell_level_key(shell)
    _increment_fp_detail_counter(detail["shell_level_counts"], level_key)
    _increment_fp_detail_counter(entry_info["shell_level_counts"], level_key)

    for interaction in list(getattr(shell, "interactions", []) or []):
        src_grp = getattr(interaction, "src_grp", None)
        trgt_grp = getattr(interaction, "trgt_grp", None)
        src_role = _group_role(src_grp)
        trgt_role = _group_role(trgt_grp)
        residue_group = _fp_detail_interaction_residue_group(src_grp, trgt_grp, src_role, trgt_role)
        if residue_group is None:
            continue
        interaction_name = _fp_detail_interaction_name(interaction)
        residue_label = _fp_detail_residue_label_from_group(residue_group)
        pair_key = interaction_name if not residue_label else f"{interaction_name}||{residue_label}"

        _increment_fp_detail_counter(detail["interaction_counts"], interaction_name)
        _increment_fp_detail_counter(entry_info["interaction_counts"], interaction_name)
        if residue_label:
            _increment_fp_detail_counter(detail["residue_counts"], residue_label)
            _increment_fp_detail_counter(entry_info["residue_counts"], residue_label)
        _increment_fp_detail_counter(detail["pair_counts"], pair_key)
        _increment_fp_detail_counter(entry_info["pair_counts"], pair_key)


def _export_fp_artifacts(proj, params, type_name):
    workdir = Path(params["workdir"])
    suffix = IFP_SUFFIX.get(type_name, type_name)
    shell_dir = workdir / "results" / "fingerprints" / "_shells" / suffix
    shell_dir.mkdir(parents=True, exist_ok=True)

    total_entries = 0
    feature_summary = {}
    entry_index = {}
    feature_details = {}

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
            shell_levels = []
            shell_level_counts = {}

            for ori_feature, found_shells in traced:
                original_features.append(int(ori_feature))
                for shell in found_shells:
                    level_key = _shell_level_key(shell)
                    shell_levels.append(level_key)
                    shell_level_counts[level_key] = shell_level_counts.get(level_key, 0) + 1
                    shell_natures = _classify_shell_natures(shell)
                    shell_nature_sets.append(shell_natures)
                    nature_tags.update(shell_natures)
                    _add_fp_detail_shell(feature_details, int(feature_id), entry_name, shell)
                    stored_shells.append(_detach_shell_for_storage(shell))

            collision = bool(raw_collision and _has_mixed_class_collision(shell_nature_sets))
            collision_level_counts = dict(shell_level_counts) if raw_collision else {}
            if collision:
                nature_tags.add("Features with collision in the same complex")
            if not nature_tags:
                nature_tags.add("Unreliable feature")
            if stored_shells:
                feature_shells[str(int(feature_id))] = stored_shells

            feature_info = feature_summary.setdefault(
                int(feature_id),
                {
                    "molecule_hits": 0,
                    "total_count": 0,
                    "collision_hits": 0,
                    "raw_collision_hits": 0,
                    "nature_counts": {},
                    "level_counts": {},
                    "collision_level_counts": {},
                },
            )
            feature_info["molecule_hits"] += 1
            feature_info["total_count"] += int(count)
            if raw_collision:
                feature_info["raw_collision_hits"] += 1
            if collision:
                feature_info["collision_hits"] += 1
            for level_key, level_count in shell_level_counts.items():
                feature_info["level_counts"][level_key] = feature_info["level_counts"].get(level_key, 0) + int(level_count)
            for level_key, level_count in collision_level_counts.items():
                feature_info["collision_level_counts"][level_key] = (
                    feature_info["collision_level_counts"].get(level_key, 0) + int(level_count)
                )
            for nature in nature_tags:
                feature_info["nature_counts"][nature] = feature_info["nature_counts"].get(nature, 0) + 1

            dominant_nature = sorted(nature_tags)[0]

            entry_rows.append(
                {
                    "feature_id": int(feature_id),
                    "count": int(count),
                    "collision": bool(collision),
                    "raw_collision": bool(raw_collision),
                    "mixed_class_collision": bool(collision),
                    "dominant_nature": dominant_nature,
                    "nature_tags": sorted(nature_tags),
                    "original_features": sorted(set(original_features)),
                    "shell_levels": _sorted_level_keys(shell_levels),
                    "shell_level_breakdown": shell_level_counts,
                    "collision_shell_levels": _sorted_level_keys(collision_level_counts.keys()),
                    "collision_level_breakdown": collision_level_counts,
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
            _save_portable_shell_payload(
                shell_dir / f"{_safe_name(entry_name)}.shells.json.gz",
                entry,
                feature_shells,
                params["pdb_dir"],
                proj,
                type_name,
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
                "raw_collision_hits": int(info.get("raw_collision_hits", 0)),
                "dominant_nature": dominant_nature,
                "nature_breakdown": info["nature_counts"],
                "shell_levels": _sorted_level_keys(info.get("level_counts", {}).keys()),
                "shell_level_breakdown": info.get("level_counts", {}),
                "collision_shell_levels": _sorted_level_keys(info.get("collision_level_counts", {}).keys()),
                "collision_level_breakdown": info.get("collision_level_counts", {}),
            }
        )

    level_assignment = _assign_feature_levels(features, bool(params.get("fp_use_otsu_threshold", False)))
    assigned_matrix = _rewrite_ifp_csv_with_feature_keys(params.get("ifp_outputs", {}).get(type_name), features)

    artifact = {
        "ifp_type": type_name,
        "ifp_label": IFP_LABELS.get(type_name, type_name),
        "random_seed": _ifp_seed(params),
        "seed_file": str(workdir / "results" / "fingerprints" / f"seed_ifp_{suffix}_importance.txt"),
        "use_otsu_threshold": bool(params.get("fp_use_otsu_threshold", False)),
        "level_assignment": level_assignment,
        "assigned_matrix": assigned_matrix,
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

    detail_artifact = {
        "ifp_type": type_name,
        "feature_details": feature_details,
        "protein_heteroatom_residues": sorted(PROTEIN_HETEROATOM_RESIDUES),
        "source": "luna_api_live_shells",
    }
    detail_path = workdir / "results" / "fingerprints" / f"fp_detail_{suffix}.json"
    _save_json(detail_path, detail_artifact)
    print(f"[luna-api] detalhes de fingerprints salvos em {detail_path}", flush=True)


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
        seed_path, ifp_seed = _write_ifp_seed(params, type_name)
        print(f"[luna-api] seed IFP {type_name}: {ifp_seed} salvo em {seed_path}", flush=True)
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
        _export_fp_artifacts(proj, params, type_name)
        if sim_output_path and proj.ifp_output and Path(proj.ifp_output).exists():
            square_path = str(Path(sim_output_path).with_name(Path(sim_output_path).stem + "_square.csv"))
            saved_edge, saved_square = _write_similarity_outputs_from_ifp(
                proj.ifp_output,
                sim_output_path,
                square_path,
            )
            if saved_edge:
                print(f"[luna-api] Similaridade {type_name} salva em {saved_edge}", flush=True)
            if saved_square:
                print(f"[luna-api] Similaridade quadrada {type_name} atualizada com niveis assinados em {saved_square}", flush=True)


AA_RESIDUES = {
    "ALA", "ARG", "ASN", "ASP", "ASH", "CYS", "CYM", "CYX", "GLN", "GLU", "GLH",
    "GLY", "HIS", "HID", "HIE", "HIP", "ILE", "LEU", "LYS", "LYN", "MET", "PHE",
    "PRO", "SER", "THR", "TRP", "TYR", "VAL",
}

WATER_RESIDUES = {"HOH", "WAT", "WTM", "TIP", "SOL", "T3P", "H2O", "OH2", "DOD"}


def _interaction_type_name(interaction):
    value = getattr(interaction, "type", None)
    name = getattr(value, "name", None) or getattr(value, "value", None)
    if name:
        return str(name)
    return str(value or interaction.__class__.__name__)


def _residue_label(atom, include_water=False, include_protein_heteroatoms=False):
    residue = getattr(atom, "parent", None)
    if residue is None:
        return None
    resname = str(getattr(residue, "resname", "") or "").upper()
    is_water = resname in WATER_RESIDUES
    if (
        resname not in AA_RESIDUES
        and not (include_water and is_water)
        and not (include_protein_heteroatoms and resname and not is_water)
    ):
        return None
    chain = getattr(getattr(residue, "parent", None), "id", "?")
    resid = getattr(residue, "id", None)
    if isinstance(resid, (tuple, list)) and len(resid) > 1:
        resid = resid[1]
    return f"{chain}/{resname}/{resid}"


def _labels_from_group(group, include_water=False, include_protein_heteroatoms=False):
    labels = set()
    for atom in getattr(group, "atoms", []) or []:
        label = _residue_label(
            atom,
            include_water=include_water,
            include_protein_heteroatoms=include_protein_heteroatoms,
        )
        if label:
            labels.add(label)
    return labels


def _atom_label(atom):
    value = getattr(atom, "name", None)
    if not value:
        try:
            value = atom.get_name()
        except Exception:
            value = ""
    if not value:
        value = getattr(atom, "id", None) or getattr(atom, "element", None)
    text = str(value or "").strip()
    if text:
        return text
    serial = getattr(atom, "serial_number", None)
    if serial is None:
        try:
            serial = atom.get_serial_number()
        except Exception:
            serial = ""
    return str(serial or "atom").strip()


def _atom_sort_key(label):
    text = str(label or "")
    parts = re.split(r"(\d+)", text)
    key = []
    for part in parts:
        if part.isdigit():
            key.append((0, int(part)))
        elif part:
            key.append((1, part.lower()))
    return key or [(1, text.lower())]


def _interaction_ligand_atom_labels(interaction):
    explicit = _explicit_interaction_atoms(interaction)
    if explicit:
        return {
            label
            for atom in explicit
            for label in [_atom_label(atom)]
            if label
        }
    labels = set()
    for group in (getattr(interaction, "src_grp", None), getattr(interaction, "trgt_grp", None)):
        if _group_role(group) != "ligand":
            continue
        for atom in getattr(group, "atoms", []) or []:
            label = _atom_label(atom)
            if label:
                labels.add(label)
    return labels


def _explicit_interaction_atoms(interaction):
    atoms = []
    for attr in (
        "src_atom", "trgt_atom", "src_atm", "trgt_atm", "atom",
        "src_atoms", "trgt_atoms", "atoms", "atom_pairs", "pairs",
    ):
        value = getattr(interaction, attr, None)
        if value is None:
            continue
        values = value if isinstance(value, (list, tuple, set)) else [value]
        for item in values:
            pair = item if isinstance(item, (list, tuple, set)) else [item]
            atoms.extend(atom for atom in pair if hasattr(atom, "name") or hasattr(atom, "get_name"))
    return atoms


def _atom_serial(atom):
    value = getattr(atom, "serial_number", None)
    if value is None:
        try:
            value = atom.get_serial_number()
        except Exception:
            value = None
    try:
        return int(value)
    except Exception:
        return None


def _interaction_ligand_atom_records(interaction):
    records = []
    for group in (getattr(interaction, "src_grp", None), getattr(interaction, "trgt_grp", None)):
        if _group_role(group) != "ligand":
            continue
        for atom in getattr(group, "atoms", []) or []:
            label = _atom_label(atom)
            if not label:
                continue
            records.append({
                "matrix_label": label,
                "serial": _atom_serial(atom),
                "element": str(getattr(atom, "element", "") or "").strip().upper(),
            })
    return records


def _interaction_residue_payload(interaction):
    src_grp = getattr(interaction, "src_grp", None)
    trgt_grp = getattr(interaction, "trgt_grp", None)
    src_role = _group_role(src_grp)
    trgt_role = _group_role(trgt_grp)
    pair = {src_role, trgt_role}
    group = None
    group_role = "unknown"
    include_water = False
    if pair == {"ligand", "protein"}:
        group = src_grp if src_role == "protein" else trgt_grp
        group_role = "protein"
    elif pair == {"protein", "water"}:
        group = src_grp if src_role == "protein" else trgt_grp
        group_role = "protein"
    elif pair == {"ligand", "water"}:
        group = src_grp if src_role == "water" else trgt_grp
        group_role = "water"
        include_water = True
    if group is None:
        return {}
    is_protein_heteroatom = (
        group_role == "protein"
        and INCLUDE_PROTEIN_HETEROATOMS
        and _group_has_ligand_residue(group)
    )
    residue_kind = "protein_heteroatom" if is_protein_heteroatom else group_role
    return {
        label: residue_kind
        for label in _labels_from_group(
            group,
            include_water=include_water,
            include_protein_heteroatoms=(group_role == "protein" and INCLUDE_PROTEIN_HETEROATOMS),
        )
    }


def _interaction_residue_labels(interaction):
    return set(_interaction_residue_payload(interaction))


def _save_json(path, data):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(data, fh, indent=2)


def _iter_ligand_source_files(ligand_path):
    path = Path(str(ligand_path or ""))
    if not path.exists():
        return []
    if path.is_file() and path.suffix.lower() in LIGAND_STRUCTURE_SUFFIXES:
        return [path]
    if not path.is_dir():
        return []
    return [
        candidate
        for candidate in sorted(path.iterdir())
        if candidate.is_file() and candidate.suffix.lower() in LIGAND_STRUCTURE_SUFFIXES
    ]


def _first_mol2_block_atom_names(path):
    names = []
    in_atoms = False
    try:
        with Path(path).open("r", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if line.startswith("@<TRIPOS>ATOM"):
                    in_atoms = True
                    continue
                if line.startswith("@<TRIPOS>") and in_atoms:
                    break
                if not in_atoms:
                    continue
                parts = line.split()
                if len(parts) >= 2:
                    names.append(parts[1].strip())
    except Exception:
        return []
    return names


def _rdkit_mol_from_source(path):
    suffix = Path(path).suffix.lower()
    try:
        from rdkit import Chem
    except Exception:
        return None
    try:
        if suffix in {".sdf", ".sd"}:
            supplier = Chem.SDMolSupplier(str(path), removeHs=False)
            for mol in supplier:
                if mol is not None:
                    return mol
            return None
        if suffix == ".mol":
            return Chem.MolFromMolFile(str(path), removeHs=False)
        if suffix == ".mol2":
            return Chem.MolFromMol2File(str(path), sanitize=True, removeHs=False)
        if suffix in {".pdb", ".ent"}:
            return Chem.MolFromPDBFile(str(path), sanitize=True, removeHs=False)
    except Exception:
        return None
    return None


def _rdkit_atom_labels_for_source(path, mol):
    labels = []
    if Path(path).suffix.lower() == ".mol2":
        labels = _first_mol2_block_atom_names(path)
    if len(labels) < int(mol.GetNumAtoms()):
        labels = []
        for atom in mol.GetAtoms():
            label = ""
            for prop_name in ("molAtomMapNumber", "_TriposAtomName", "atomLabel", "name"):
                try:
                    if atom.HasProp(prop_name):
                        label = str(atom.GetProp(prop_name)).strip()
                        if label:
                            break
                except Exception:
                    pass
            if not label:
                label = f"{atom.GetSymbol()}{atom.GetIdx() + 1}"
            labels.append(label)
    return [str(label or idx + 1) for idx, label in enumerate(labels[: int(mol.GetNumAtoms())])]


def _normalized_atom_label(value):
    return re.sub(r"[^A-Za-z0-9]+", "", str(value or "")).upper()


def _numeric_atom_token(value):
    values = re.findall(r"\d+", str(value or ""))
    return int(values[-1]) if values else None


def _align_atom_map_labels(source_labels, matrix_labels, source_records=None, matrix_metadata=None, return_mapping=False):
    """Use the matrix IDs on the drawing while preserving the molecule atom order."""
    source = [str(label) for label in source_labels]
    preferred = [str(label) for label in (matrix_labels or []) if str(label).strip()]
    if not preferred:
        mapping = [
            {"source_index": index, "source_label": label, "display_label": label, "matrix_label": "", "matched": False}
            for index, label in enumerate(source)
        ]
        return (source, mapping) if return_mapping else source

    source_records = list(source_records or [])
    metadata = dict(matrix_metadata or {})
    while len(source_records) < len(source):
        source_records.append({})

    preferred_by_key = {}
    for label in preferred:
        preferred_by_key.setdefault(_normalized_atom_label(label), []).append(label)
    aligned = []
    used = set()
    for label in source:
        matches = preferred_by_key.get(_normalized_atom_label(label), [])
        match = next((candidate for candidate in matches if candidate not in used), None)
        aligned.append(match or "")
        if match:
            used.add(match)

    def assign(source_index, candidate):
        if not candidate or candidate in used or aligned[source_index]:
            return False
        aligned[source_index] = candidate
        used.add(candidate)
        return True

    # LUNA and RDKit may expose different names for the same atom. Serial
    # numbers are the stable bridge for PDB and MOL2 inputs.
    for source_index, record in enumerate(source_records[:len(source)]):
        if aligned[source_index]:
            continue
        source_serial = record.get("serial")
        if source_serial is None:
            continue
        candidates = [
            label for label in preferred
            if label not in used and (metadata.get(label) or {}).get("serial") == source_serial
        ]
        if len(candidates) == 1:
            assign(source_index, candidates[0])

    # Some parsers drop serial metadata but retain it as the numeric suffix of
    # the matrix identifier. Only accept unambiguous matches.
    for source_index, record in enumerate(source_records[:len(source)]):
        if aligned[source_index]:
            continue
        source_tokens = {
            value for value in (
                record.get("serial"),
                record.get("source_index", source_index) + 1,
                _numeric_atom_token(source[source_index]),
            )
            if value is not None
        }
        candidates = [
            label for label in preferred
            if label not in used and _numeric_atom_token(label) in source_tokens
        ]
        if len(candidates) == 1:
            assign(source_index, candidates[0])

    remaining = [label for label in preferred if label not in used]
    if len(preferred) == len(source):
        remaining_iter = iter(remaining)
        aligned = [label or next(remaining_iter) for label in aligned]
    else:
        aligned = [label or source[index] for index, label in enumerate(aligned)]
    mapping = [
        {
            "source_index": index,
            "source_label": source[index],
            "display_label": aligned[index],
            "matrix_label": aligned[index] if aligned[index] in preferred else "",
            "matched": aligned[index] in preferred,
        }
        for index in range(len(source))
    ]
    return (aligned, mapping) if return_mapping else aligned


def _rdkit_source_atom_records(path, mol, labels):
    records = []
    mol2_serials = []
    if Path(path).suffix.lower() == ".mol2":
        try:
            in_atoms = False
            for line in Path(path).read_text(encoding="utf-8", errors="replace").splitlines():
                if line.startswith("@<TRIPOS>ATOM"):
                    in_atoms = True
                    continue
                if line.startswith("@<TRIPOS>") and in_atoms:
                    break
                if in_atoms:
                    parts = line.split()
                    if parts:
                        mol2_serials.append(int(parts[0]))
        except Exception:
            mol2_serials = []
    for atom in mol.GetAtoms():
        index = int(atom.GetIdx())
        serial = mol2_serials[index] if index < len(mol2_serials) else None
        try:
            pdb_info = atom.GetPDBResidueInfo()
            if pdb_info is not None:
                serial = int(pdb_info.GetSerialNumber())
        except Exception:
            pass
        records.append({
            "source_index": index,
            "source_label": labels[index] if index < len(labels) else f"{atom.GetSymbol()}{index + 1}",
            "serial": serial if serial is not None else index + 1,
            "element": str(atom.GetSymbol() or "").strip().upper(),
        })
    return records


def _rdkit_heavy_atom_mol_and_labels(path, mol, matrix_labels=None, matrix_metadata=None):
    try:
        from rdkit import Chem
    except Exception:
        return mol, _rdkit_atom_labels_for_source(path, mol)
    source_mol = Chem.Mol(mol)
    source_labels = _rdkit_atom_labels_for_source(path, source_mol)
    source_records = _rdkit_source_atom_records(path, source_mol, source_labels)
    for atom in source_mol.GetAtoms():
        idx = int(atom.GetIdx())
        label = source_labels[idx] if idx < len(source_labels) else f"{atom.GetSymbol()}{idx + 1}"
        atom.SetProp("_hip2l_atom_label", str(label))
        atom.SetIntProp("_hip2l_source_index", idx)
    try:
        heavy_mol = Chem.RemoveHs(source_mol, sanitize=False)
    except Exception:
        heavy_mol = source_mol
    labels = []
    heavy_records = []
    for atom in heavy_mol.GetAtoms():
        try:
            label = atom.GetProp("_hip2l_atom_label")
        except Exception:
            label = f"{atom.GetSymbol()}{atom.GetIdx() + 1}"
        labels.append(str(label))
        try:
            source_index = int(atom.GetIntProp("_hip2l_source_index"))
        except Exception:
            source_index = int(atom.GetIdx())
        heavy_records.append(source_records[source_index] if source_index < len(source_records) else {"source_index": source_index})
    aligned, mapping = _align_atom_map_labels(
        labels,
        matrix_labels,
        source_records=heavy_records,
        matrix_metadata=matrix_metadata,
        return_mapping=True,
    )
    return heavy_mol, aligned, mapping


def _export_ligand_atom_map(params, residue_artifact):
    if not bool(params.get("trajectory_analysis", False)):
        return ""
    workdir = Path(params["workdir"])
    output_path = workdir / "results" / "ligand_atom_map.png"
    meta_path = workdir / "results" / "ligand_atom_map.json"
    for source in _iter_ligand_source_files(params.get("lig_file")):
        mol = _rdkit_mol_from_source(source)
        if mol is None:
            continue
        try:
            from rdkit import Chem
            from rdkit.Chem import AllChem, Draw

            matrix_labels = list(residue_artifact.get("ligand_atoms", []) or [])
            matrix_metadata = dict(residue_artifact.get("ligand_atom_metadata", {}) or {})
            draw_mol, labels, label_mapping = _rdkit_heavy_atom_mol_and_labels(
                source,
                mol,
                matrix_labels,
                matrix_metadata,
            )
            try:
                AllChem.Compute2DCoords(draw_mol)
            except Exception:
                pass
            options = Draw.MolDrawOptions()
            options.addAtomIndices = False
            for idx, label in enumerate(labels):
                options.atomLabels[idx] = str(label)
            drawer = Draw.MolDraw2DCairo(900, 650)
            drawer.SetDrawOptions(options)
            drawer.DrawMolecule(draw_mol)
            drawer.FinishDrawing()
            output_path.parent.mkdir(parents=True, exist_ok=True)
            output_path.write_bytes(drawer.GetDrawingText())
            meta = {
                "image": str(output_path),
                "source_file": str(source),
                "atom_labels": labels,
                "matrix_atom_labels": matrix_labels,
                "atom_label_mapping": label_mapping,
                "labels_match_matrix": all(label in labels for label in matrix_labels),
                "heavy_atoms_only": True,
            }
            _save_json(meta_path, meta)
            residue_artifact["ligand_atom_map"] = str(output_path)
            residue_artifact["ligand_atom_map_source"] = str(source)
            print(f"[luna-api] mapa 2D do ligante salvo em {output_path}", flush=True)
            return str(output_path)
        except Exception as ex:
            print(f"[luna-api] aviso: mapa 2D do ligante falhou para {source}: {ex}", flush=True)
            continue
    print("[luna-api] aviso: nao foi possivel gerar mapa 2D do ligante para estatisticas.", flush=True)
    return ""


def _export_summary_artifacts(proj, workdir, params=None):
    summary = {
        "entries": 0,
        "interaction_counts": {},
        "entry_interaction_counts": {},
        "errors": [],
    }
    residue_counts = {}
    ligand_atom_counts = {}
    entries = []
    residues = set()
    residue_roles = {}
    ligand_atoms = set()
    ligand_atom_metadata = {}
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

            touched_residues = _interaction_residue_payload(interaction)
            for label, residue_role in touched_residues.items():
                residues.add(label)
                residue_roles[label] = residue_role
                by_type = residue_counts.setdefault(itype, {})
                by_entry = by_type.setdefault(entry_name, {})
                by_entry[label] = by_entry.get(label, 0) + 1

            touched_ligand_atoms = _interaction_ligand_atom_labels(interaction)
            for label in touched_ligand_atoms:
                ligand_atoms.add(label)
                by_type = ligand_atom_counts.setdefault(itype, {})
                by_entry = by_type.setdefault(entry_name, {})
                by_entry[label] = by_entry.get(label, 0) + 1
            for record in _interaction_ligand_atom_records(interaction):
                label = str(record.get("matrix_label") or "")
                if label:
                    ligand_atom_metadata.setdefault(label, record)

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

    ligand_atom_labels = sorted(ligand_atoms, key=_atom_sort_key)
    ligand_atom_matrix = {}
    for itype in sorted(interaction_types):
        by_entry = ligand_atom_counts.get(itype, {})
        ligand_atom_matrix[itype] = [
            [float(by_entry.get(entry_name, {}).get(label, 0.0)) for label in ligand_atom_labels]
            for entry_name in ordered_entries
        ]

    residue_artifact = {
        "interaction_types": sorted(interaction_types),
        "residues": residue_labels,
        "protein_heteroatom_residues": sorted(PROTEIN_HETEROATOM_RESIDUES),
        "residue_roles": {
            label: residue_roles.get(label, "protein")
            for label in residue_labels
        },
        "ligand_atoms": ligand_atom_labels,
        "ligand_atom_metadata": {
            label: ligand_atom_metadata.get(label, {"matrix_label": label})
            for label in ligand_atom_labels
        },
        "entries": ordered_entries,
        "matrix": matrix,
        "ligand_atom_matrix": ligand_atom_matrix,
        "errors": list(summary["errors"]),
    }
    _export_ligand_atom_map(params or {}, residue_artifact)

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
    mol_file = spec.get("mol_file") or lig_file
    spec_ext = Path(mol_file).suffix.lower()
    spec_mol_obj_type = spec.get("mol_obj_type") or (
        "rdkit" if spec_ext in {".sdf", ".sd", ".mol"} else lig_mol_obj_type
    )
    is_multimol_file = bool(spec.get("is_multimol_file", True))
    try:
        e = MolFileEntry.from_mol_file(
            pdb_id=pdb_id,
            mol_id=name,
            mol_file=mol_file,
            is_multimol_file=is_multimol_file,
            mol_obj_type=spec_mol_obj_type,
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
    # A protein-residue/HETATM contact is the geometric basis for a ligand
    # contact with a receptor cofactor or coordinated ion.  Make LUNA's
    # permissive default explicit so it remains enabled across versions.
    ignore_res_hetatm=False,
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
_export_summary_artifacts(proj, workdir, p)
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
                # InteractionViewer reuses PyMOL's process-global object store. Reset it
                # for every entry so a .pse cannot inherit prior ligands/interactions.
                _reset_pymol_session()
                viewer = InteractionViewer(show_hydrop_surface=False)
                colored = _save_viewer_session_with_palette(
                    viewer,
                    [(entry, inter, pdb_file)],
                    pse_path,
                    p.get("pse_interaction_colors") or {},
                )
            except Exception as ex:
                print(f"[warn] PSE para {safe_name} falhou: {ex}", flush=True)
            finally:
                _reset_pymol_session()
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
    if protein_path.is_dir():
        return sorted(protein_path.glob("*.pdb"))
    return [protein_path]


def protein_heteroatom_residue_keys(protein_file: str | Path) -> list[str]:
    """Return non-water HETATM residue keys declared in receptor PDB input(s).

    LUNA deliberately classes every PDB ``HETATM`` as a hetero group.  The
    entry therefore does not retain enough role information to tell a zinc in
    the receptor PDB from a zinc supplied with a ligand.  Persisting these
    keys in the run parameters preserves the source-file decision throughout
    the API runner and all post-processing helpers.
    """
    keys: set[str] = set()
    for pdb_path in _candidate_protein_files(
        ProjectConfig(protein_file=str(protein_file))
    ):
        try:
            with pdb_path.open("r", encoding="utf-8", errors="replace") as fh:
                for line in fh:
                    if not line.startswith("HETATM"):
                        continue
                    resname = line[17:20].strip().upper()
                    if not resname or resname in _WATER_RESIDUE_NAMES:
                        continue
                    chain = line[21:22].strip() or "?"
                    resid = line[22:26].strip()
                    if resid:
                        keys.add(f"{chain}/{resname}/{resid}")
        except OSError:
            continue
    return sorted(keys)


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
    protein_is_dir = bool(cfg.protein_file) and Path(cfg.protein_file).is_dir()
    ligand_is_dir = bool(cfg.ligand_file) and Path(cfg.ligand_file).is_dir()
    return (
        cfg.uses_python_api()
        or protein_is_dir
        or ligand_is_dir
        or not flags["amend_mol"]
    )


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


def _ligand_file_index(ligand_dir: Path) -> tuple[dict[str, Path], set[str]]:
    index: dict[str, Path] = {}
    duplicates: set[str] = set()
    for ligand_file in sorted(ligand_dir.iterdir()):
        if not ligand_file.is_file() or ligand_file.suffix.lower() not in _LIGAND_FILE_SUFFIXES:
            continue
        key = _normalize_complex_name(ligand_file.stem)
        if key in index and index[key] != ligand_file:
            duplicates.add(key)
            continue
        index[key] = ligand_file
    return index, duplicates


def build_entry_specs(cfg: ProjectConfig, entries: list[str]) -> list[dict[str, str]]:
    """Return the per-entry protein/ligand mapping for the API runner."""
    protein_path = Path(cfg.protein_file)
    ligand_path = Path(cfg.ligand_file) if cfg.ligand_file else None
    ligand_index: dict[str, Path] = {}
    ligand_duplicates: set[str] = set()
    ligand_is_dir = ligand_path is not None and ligand_path.is_dir()
    if ligand_is_dir:
        ligand_index, ligand_duplicates = _ligand_file_index(ligand_path)

    protein_dir = protein_path if protein_path.is_dir() else protein_path.parent
    protein_index: dict[str, str] = {}
    protein_duplicates: set[str] = set()
    if protein_path.is_dir():
        protein_index, protein_duplicates = _protein_index(protein_dir)

    single_pdb_id = protein_path.stem
    errors: list[str] = []
    specs: list[dict[str, str]] = []
    for ligand_name in entries:
        key = _normalize_complex_name(ligand_name)
        if key in ligand_duplicates:
            errors.append(
                f"Mais de um arquivo de ligante corresponde a '{ligand_name}' na pasta {ligand_path}."
            )
            continue

        spec: dict[str, str] = {"ligand_name": ligand_name}
        if ligand_is_dir:
            mol_file = ligand_index.get(key)
            if mol_file is None:
                errors.append(
                    f"Nenhum arquivo de ligante com o mesmo nome de '{ligand_name}' foi encontrado em {ligand_path}."
                )
                continue
            spec["mol_file"] = str(mol_file)
            spec["mol_obj_type"] = ligand_mol_obj_type(mol_file)
            spec["is_multimol_file"] = False

        if not protein_path.is_dir():
            spec["pdb_id"] = single_pdb_id
            specs.append(spec)
            continue

        if key in protein_duplicates:
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
        spec["pdb_id"] = pdb_id
        specs.append(spec)

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


def read_ifp_seed_file(path: str | Path | None) -> int:
    """Read the first integer seed from a user-provided text file."""
    raw_path = str(path or "").strip()
    if not raw_path:
        return 0
    seed_path = Path(raw_path)
    if not seed_path.exists():
        raise ValueError(f"Arquivo seed IFP nao encontrado: {seed_path}")
    text = seed_path.read_text(encoding="utf-8", errors="replace")
    match = re.search(r"-?\d+", text)
    if match is None:
        raise ValueError(f"Arquivo seed IFP nao contem um inteiro: {seed_path}")
    return int(match.group(0))


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
        "nproc": safe_nproc(cfg.nproc),
        "add_h": protein_flags["add_h"],
        "ph": protein_flags["ph"],
        "amend_mol": protein_flags["amend_mol"],
        "protein_has_hydrogens": protein_flags["protein_has_hydrogens"],
        "protein_is_preprocessed": protein_flags["protein_is_preprocessed"],
        "stage_protein_without_h": protein_flags["stage_protein_without_h"],
        "stage_ligand_without_h": protein_flags["stage_ligand_without_h"],
        "trajectory_analysis": cfg.trajectory_analysis,
        "include_waters": cfg.include_waters,
        "include_protein_heteroatoms": cfg.include_protein_heteroatoms,
        "protein_heteroatom_residues": (
            protein_heteroatom_residue_keys(cfg.protein_file)
            if cfg.include_protein_heteroatoms else []
        ),
        "out_ifp": cfg.out_ifp,
        "ifp_type": cfg.ifp_type,
        "ifp_types": cfg.selected_ifp_types() if (cfg.out_ifp or cfg.sim_matrix) else [],
        "ifp_levels": cfg.ifp_levels,
        "ifp_radius": cfg.ifp_radius,
        "ifp_length": cfg.ifp_length,
        "ifp_bit": cfg.ifp_bit,
        "ifp_output": cfg.ifp_output,
        "ifp_seed_file": cfg.ifp_seed_file,
        "ifp_seed": read_ifp_seed_file(cfg.ifp_seed_file),
        "fp_use_otsu_threshold": cfg.fp_use_otsu_threshold,
        "ifp_outputs": resolve_ifp_output_paths(cfg) if (cfg.out_ifp or cfg.sim_matrix) else {},
        "sim_matrix": cfg.sim_matrix,
        "sim_matrix_output": cfg.sim_matrix_output,
        "sim_matrix_outputs": resolve_sim_matrix_output_paths(cfg),
        "out_pse": cfg.out_pse,
        "pse_path": cfg.pse_path,
        "pse_interaction_types": cfg.pse_interaction_types,
        "pse_interaction_colors": dict(INTERACTION_COLORS),
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
