"""Run LUNA's post-analysis API inside the luna-env via subprocess.

The GUI itself runs in a clean Python with no LUNA installed, so we
serialize the helper script as a string and feed it to `python -c`
through the env's interpreter.
"""
from __future__ import annotations

import json
import subprocess
from pathlib import Path

from .env_manager import python_process_env

# Helper script executed inside luna-env. It walks the workdir, loads each
# EntryResults pickle, and aggregates interaction-type counts and per-residue
# counts. Output is a single JSON document on stdout.
HELPER_SCRIPT = r"""
import sys, os, json, glob
workdir = sys.argv[1]
_project_options = {}
for _options_path in ("_luna_api_params.json", ".luna_gui.json"):
    try:
        with open(os.path.join(workdir, _options_path), "r", encoding="utf-8") as fh:
            _project_options = json.load(fh)
        break
    except Exception:
        pass
INCLUDE_PROTEIN_HETEROATOMS = bool(
    _project_options.get("include_protein_heteroatoms", False)
)
INCLUDE_WATERS = bool(_project_options.get("include_waters", False))
PROTEIN_HETEROATOM_RESIDUES = {
    str(value).strip()
    for value in (_project_options.get("protein_heteroatom_residues") or [])
    if str(value).strip()
}
out = {
    "entries": 0,
    "interaction_counts": {},
    "residue_counts": {},
    "residue_roles": {},
    "errors": [],
}
try:
    from luna.projects import EntryResults
except Exception as e:
    print(json.dumps({"error": "luna import failed: %s" % e}))
    sys.exit(0)
try:
    from luna.analysis.summary import count_interaction_types
except Exception:
    count_interaction_types = None

AA_RESIDUES = {
    "ALA", "ARG", "ASN", "ASP", "ASH", "CYS", "CYM", "CYX", "GLN", "GLU", "GLH",
    "GLY", "HIS", "HID", "HIE", "HIP", "HSD", "HSE", "HSP", "ILE", "LEU", "LYS",
    "LYN", "MET", "MSE", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL", "SEC",
}
WATER_RESIDUES = {"HOH", "WAT", "WTM", "TIP", "SOL", "T3P", "H2O", "OH2", "DOD"}


def _residue_name(atom):
    residue = getattr(atom, "parent", None)
    if residue is None:
        return ""
    value = getattr(residue, "resname", None)
    if value is None:
        try:
            value = residue.get_resname()
        except Exception:
            value = ""
    return str(value or "").strip().upper()


def _group_has_water_residue(group):
    return any(
        _residue_name(atom) in WATER_RESIDUES
        for atom in (getattr(group, "atoms", []) or [])
    )


def _group_has_structural_residue(group):
    try:
        return bool(group.has_residue() or group.has_nucleotide())
    except Exception:
        return False


def _group_has_ligand_residue(group):
    return any(
        _residue_name(atom) not in (AA_RESIDUES | WATER_RESIDUES)
        for atom in (getattr(group, "atoms", []) or [])
    )


def _residue_key(atom):
    residue = getattr(atom, "parent", None)
    if residue is None:
        return ""
    resname = _residue_name(atom)
    residue_id = getattr(residue, "id", None)
    if isinstance(residue_id, (tuple, list)) and len(residue_id) > 1:
        residue_id = residue_id[1]
    chain = getattr(getattr(residue, "parent", None), "id", "?")
    if not resname or residue_id in (None, ""):
        return ""
    return "%s/%s/%s" % (chain, resname, residue_id)


def _group_is_configured_protein_heteroatom(group):
    return any(
        _residue_key(atom) in PROTEIN_HETEROATOM_RESIDUES
        for atom in (getattr(group, "atoms", []) or [])
    )


def _group_role(group):
    if group is None:
        return "unknown"
    try:
        if group.has_water():
            return "water"
    except Exception:
        pass
    if _group_has_water_residue(group):
        return "water"
    if INCLUDE_PROTEIN_HETEROATOMS and _group_is_configured_protein_heteroatom(group):
        return "protein"
    if _group_has_structural_residue(group):
        try:
            is_hetatm = bool(group.has_hetatm())
        except Exception:
            is_hetatm = False
        if (
            not is_hetatm
            or not _group_has_ligand_residue(group)
            or INCLUDE_PROTEIN_HETEROATOMS
        ):
            return "protein"
    try:
        if group.has_hetatm():
            return "ligand"
    except Exception:
        pass
    return "other"


def _residue_label(atom, *, include_water=False, include_protein_heteroatoms=False):
    residue = getattr(atom, "parent", None)
    if residue is None:
        return None
    resname = _residue_name(atom)
    if (
        resname not in AA_RESIDUES
        and not (include_water and resname in WATER_RESIDUES)
        and not (include_protein_heteroatoms and resname and resname not in WATER_RESIDUES)
    ):
        return None
    residue_id = getattr(residue, "id", None)
    if isinstance(residue_id, (tuple, list)) and len(residue_id) > 1:
        residue_id = residue_id[1]
    chain = getattr(getattr(residue, "parent", None), "id", "?")
    return "%s/%s/%s" % (chain, resname, residue_id)


def _interaction_residue_payload(interaction):
    src_group = getattr(interaction, "src_grp", None)
    trgt_group = getattr(interaction, "trgt_grp", None)
    src_role = _group_role(src_group)
    trgt_role = _group_role(trgt_group)
    pair = {src_role, trgt_role}
    group = None
    group_role = "unknown"
    include_water = False
    if pair == {"ligand", "protein"} or pair == {"protein", "water"}:
        group = src_group if src_role == "protein" else trgt_group
        group_role = "protein"
    elif pair == {"ligand", "water"} and INCLUDE_WATERS:
        group = src_group if src_role == "water" else trgt_group
        group_role = "water"
        include_water = True
    if group is None:
        return {}
    role = group_role
    if group_role == "protein" and INCLUDE_PROTEIN_HETEROATOMS:
        has_nonstandard = any(
            _residue_name(atom) not in (AA_RESIDUES | WATER_RESIDUES)
            for atom in (getattr(group, "atoms", []) or [])
        )
        if has_nonstandard:
            role = "protein_heteroatom"
    payload = {}
    for atom in (getattr(group, "atoms", []) or []):
        label = _residue_label(
            atom,
            include_water=include_water,
            include_protein_heteroatoms=(group_role == "protein" and INCLUDE_PROTEIN_HETEROATOMS),
        )
        if label:
            payload[label] = role
    return payload

# LUNA stores per-entry pickles under <workdir>/results/ (compressed).
patterns = [
    os.path.join(workdir, "results", "**", "*.pkl.gz"),
    os.path.join(workdir, "results", "**", "*.pkl"),
    os.path.join(workdir, "**", "*.pkl.gz"),
]
seen = set()
files = []
for pat in patterns:
    for f in glob.glob(pat, recursive=True):
        if f in seen:
            continue
        seen.add(f); files.append(f)

for f in files:
    try:
        er = EntryResults.load(f)
    except Exception:
        continue
    im = getattr(er, "interactions_mngr", None)
    if im is None:
        continue
    interactions = list(getattr(im, "interactions", []) or [])
    if not interactions:
        continue
    out["entries"] += 1
    # Tally interaction types
    if count_interaction_types is not None:
        try:
            counts = count_interaction_types(interactions, must_have_target=True)
            for k, v in counts.items():
                out["interaction_counts"][k] = out["interaction_counts"].get(k, 0) + int(v)
        except Exception as e:
            out["errors"].append("count_interaction_types: %s" % e)
    # Tally the same receptor-side residues exported by the API runner.  This
    # is a compatibility path for older projects without cached artifacts.
    for it in interactions:
        interaction_residues = _interaction_residue_payload(it)
        try:
            for name, role in interaction_residues.items():
                out["residue_roles"][name] = role
        except Exception:
            pass
        for name in interaction_residues:
            out["residue_counts"][name] = out["residue_counts"].get(name, 0) + 1
print(json.dumps(out))
"""


RESIDUE_MATRIX_SCRIPT = r"""
import sys, os, json, glob, re
workdir = sys.argv[1]
_project_options = {}
for _options_path in ("_luna_api_params.json", ".luna_gui.json"):
    try:
        with open(os.path.join(workdir, _options_path), "r", encoding="utf-8") as fh:
            _project_options = json.load(fh)
        break
    except Exception:
        pass
try:
    INCLUDE_PROTEIN_HETEROATOMS = bool(
        _project_options.get("include_protein_heteroatoms", False)
    )
except Exception:
    INCLUDE_PROTEIN_HETEROATOMS = False
PROTEIN_HETEROATOM_RESIDUES = {
    str(value).strip()
    for value in (_project_options.get("protein_heteroatom_residues") or [])
    if str(value).strip()
}
out = {
    "interaction_types": [],
    "residues": [],
    "ligand_atoms": [],
    "entries": [],
    "matrix": {},
    "ligand_atom_matrix": {},
    "protein_heteroatom_residues": sorted(PROTEIN_HETEROATOM_RESIDUES),
    "residue_roles": {},
    "errors": [],
}
try:
    from luna.projects import EntryResults
except Exception as e:
    print(json.dumps({"error": "luna import failed: %s" % e}))
    sys.exit(0)

AA_RESIDUES = {
    "ALA", "ARG", "ASN", "ASP", "ASH", "CYS", "CYM", "CYX", "GLN", "GLU", "GLH",
    "GLY", "HIS", "HID", "HIE", "HIP", "HSD", "HSE", "HSP", "ILE", "LEU", "LYS",
    "LYN", "MET", "MSE", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL", "SEC",
}
WATER_RESIDUES = {"HOH", "WAT", "WTM", "TIP", "SOL", "T3P", "H2O", "OH2", "DOD"}
try:
    INCLUDE_WATERS = bool(_project_options.get("include_waters", False))
except Exception:
    INCLUDE_WATERS = False


def _residue_name(atom):
    residue = getattr(atom, "parent", None)
    if residue is None:
        return ""
    value = getattr(residue, "resname", None)
    if value is None:
        try:
            value = residue.get_resname()
        except Exception:
            value = ""
    return str(value or "").strip().upper()


def _group_has_water_residue(atm_grp):
    for atom in getattr(atm_grp, "atoms", []) or []:
        if _residue_name(atom) in WATER_RESIDUES:
            return True
    return False


def _group_has_ligand_residue(atm_grp):
    for atom in getattr(atm_grp, "atoms", []) or []:
        resname = _residue_name(atom)
        if resname and resname not in AA_RESIDUES and resname not in WATER_RESIDUES:
            return True
    return False


def _group_has_structural_residue(atm_grp):
    try:
        return bool(atm_grp.has_residue() or atm_grp.has_nucleotide())
    except Exception:
        return False


def _residue_key(atom):
    residue = getattr(atom, "parent", None)
    if residue is None:
        return ""
    resname = _residue_name(atom)
    residue_id = getattr(residue, "id", None)
    if isinstance(residue_id, (tuple, list)) and len(residue_id) > 1:
        residue_id = residue_id[1]
    chain = getattr(getattr(residue, "parent", None), "id", "?")
    if not resname or residue_id in (None, ""):
        return ""
    return "%s/%s/%s" % (chain, resname, residue_id)


def _group_is_configured_protein_heteroatom(atm_grp):
    return any(
        _residue_key(atom) in PROTEIN_HETEROATOM_RESIDUES
        for atom in (getattr(atm_grp, "atoms", []) or [])
    )


def _group_role(atm_grp):
    if atm_grp is None:
        return "unknown"
    try:
        if atm_grp.has_water():
            return "water"
    except Exception:
        pass
    if _group_has_water_residue(atm_grp):
        return "water"
    if INCLUDE_PROTEIN_HETEROATOMS and _group_is_configured_protein_heteroatom(atm_grp):
        return "protein"
    is_structural = _group_has_structural_residue(atm_grp)
    try:
        is_hetatm = bool(atm_grp.has_hetatm())
    except Exception:
        is_hetatm = False
    if (
        is_structural
        and (
            not is_hetatm
            or not _group_has_ligand_residue(atm_grp)
            or INCLUDE_PROTEIN_HETEROATOMS
        )
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


def _interaction_type_name(interaction):
    value = getattr(interaction, "type", None)
    name = getattr(value, "name", None) or getattr(value, "value", None)
    if name:
        return str(name)
    return str(value or interaction.__class__.__name__)


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


def _residue_label(atom, include_water=False, include_protein_heteroatoms=False):
    residue = getattr(atom, "parent", None)
    if residue is None:
        return None
    resname = _residue_name(atom)
    if (
        resname not in AA_RESIDUES
        and not (include_water and resname in WATER_RESIDUES)
        and not (include_protein_heteroatoms and resname and resname not in WATER_RESIDUES)
    ):
        return None
    residue_id = getattr(residue, "id", None)
    if isinstance(residue_id, (tuple, list)) and len(residue_id) > 1:
        residue_id = residue_id[1]
    chain = getattr(getattr(residue, "parent", None), "id", "?")
    return "%s/%s/%s" % (chain, resname, residue_id)


def _interaction_residue_payload(interaction):
    src_group = getattr(interaction, "src_grp", None)
    trgt_group = getattr(interaction, "trgt_grp", None)
    src_role = _group_role(src_group)
    trgt_role = _group_role(trgt_group)
    pair = {src_role, trgt_role}
    group = None
    group_role = "unknown"
    include_water = False
    if pair == {"ligand", "protein"} or pair == {"protein", "water"}:
        group = src_group if src_role == "protein" else trgt_group
        group_role = "protein"
    elif pair == {"ligand", "water"} and INCLUDE_WATERS:
        group = src_group if src_role == "water" else trgt_group
        group_role = "water"
        include_water = True
    if group is None:
        return {}
    residue_role = group_role
    if group_role == "protein" and INCLUDE_PROTEIN_HETEROATOMS:
        if _group_has_ligand_residue(group):
            residue_role = "protein_heteroatom"
    payload = {}
    for atom in (getattr(group, "atoms", []) or []):
        label = _residue_label(
            atom,
            include_water=include_water,
            include_protein_heteroatoms=(group_role == "protein" and INCLUDE_PROTEIN_HETEROATOMS),
        )
        if label:
            payload[label] = residue_role
    return payload


def _interaction_ligand_atom_labels(interaction):
    # Prefer the atom references attached to the interaction itself.  Group
    # atoms may contain an entire aromatic/cationic system (for example in a
    # cation-pi contact), while the interaction object can identify the exact
    # participating atom(s).
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

patterns = [
    os.path.join(workdir, "results", "**", "*.pkl.gz"),
    os.path.join(workdir, "results", "**", "*.pkl"),
    os.path.join(workdir, "**", "*.pkl.gz"),
]
seen = set(); files = []
for pat in patterns:
    for f in glob.glob(pat, recursive=True):
        if f in seen: continue
        seen.add(f); files.append(f)

entry_names = []
ligand_atoms = set()
ligand_atom_values = {}
residue_values = {}
residue_roles = {}
for f in files:
    try:
        er = EntryResults.load(f)
    except Exception:
        continue
    im = getattr(er, "interactions_mngr", None)
    if im is None: continue
    name = ""
    try:
        name = er.entry.to_string()
    except Exception:
        name = os.path.basename(f)
    entry_names.append(name)
    for interaction in list(getattr(im, "interactions", []) or []):
        itype = _interaction_type_name(interaction)
        for label, residue_role in _interaction_residue_payload(interaction).items():
            residue_roles[label] = residue_role
            by_entry = residue_values.setdefault(itype, {}).setdefault(name, {})
            by_entry[label] = by_entry.get(label, 0) + 1
        labels = _interaction_ligand_atom_labels(interaction)
        if not labels:
            continue
        by_entry = ligand_atom_values.setdefault(itype, {}).setdefault(name, {})
        for label in labels:
            ligand_atoms.add(label)
            by_entry[label] = by_entry.get(label, 0) + 1

if not entry_names:
    print(json.dumps(out)); sys.exit(0)

try:
    residues = sorted(residue_roles)
    entries_in_order = entry_names
    interaction_types = sorted(set(residue_values) | set(ligand_atom_values))
    matrix = {
        itype: [
            [float(residue_values.get(itype, {}).get(entry, {}).get(label, 0.0)) for label in residues]
            for entry in entries_in_order
        ]
        for itype in interaction_types
    }

    ligand_atom_labels = sorted(ligand_atoms, key=_atom_sort_key)
    ligand_atom_matrix = {}
    for it in sorted(ligand_atom_values):
        by_entry = ligand_atom_values.get(it, {})
        ligand_atom_matrix[it] = [
            [float(by_entry.get(e, {}).get(label, 0.0)) for label in ligand_atom_labels]
            for e in entries_in_order
        ]
    out["interaction_types"] = interaction_types
    out["residues"] = residues
    out["ligand_atoms"] = ligand_atom_labels
    out["entries"] = entries_in_order
    out["matrix"] = matrix
    out["ligand_atom_matrix"] = ligand_atom_matrix
    out["residue_roles"] = {label: residue_roles[label] for label in residues}
except Exception as e:
    out["errors"].append("matrix build: %s" % e)

print(json.dumps(out))
"""


def run_residue_matrix(py_exe: str, workdir: str, timeout: int = 600) -> dict:
    """Load residue-interaction matrix from a completed workdir."""
    if not py_exe or not Path(py_exe).exists():
        return {"error": "Python do luna-env não encontrado."}
    if not Path(workdir).exists():
        return {"error": f"Workdir não existe: {workdir}"}
    try:
        r = subprocess.run(
            [py_exe, "-c", RESIDUE_MATRIX_SCRIPT, workdir],
            capture_output=True, text=True, timeout=timeout,
        )
    except subprocess.TimeoutExpired:
        return {"error": "Residue-matrix excedeu o tempo limite."}
    except Exception as e:
        return {"error": str(e)}
    if r.returncode != 0:
        return {"error": f"Helper falhou:\n{r.stderr.strip()}"}
    try:
        result = json.loads(r.stdout.strip().splitlines()[-1])
    except Exception:
        return {"error": f"Saída inválida do helper:\n{r.stdout[:500]}"}
    if "error" not in result:
        try:
            out_path = Path(workdir) / "results" / "residue_matrix.json"
            out_path.parent.mkdir(parents=True, exist_ok=True)
            out_path.write_text(json.dumps(result, indent=2), encoding="utf-8")
        except Exception:
            pass
    return result


def run_analysis(py_exe: str, workdir: str, timeout: int = 900) -> dict:
    """Execute the helper inside luna-env. Returns parsed JSON or {'error': ...}."""
    if not py_exe or not Path(py_exe).exists():
        return {"error": "Python do luna-env não encontrado."}
    if not Path(workdir).exists():
        return {"error": f"Workdir não existe: {workdir}"}
    try:
        r = subprocess.run(
            [py_exe, "-c", HELPER_SCRIPT, workdir],
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
            timeout=timeout,
            cwd=workdir,
            env=python_process_env(py_exe),
        )
    except subprocess.TimeoutExpired:
        return {
            "error": (
                f"A análise excedeu o limite de {timeout} segundos. "
                "O processo foi interrompido com segurança; a interface permanece aberta."
            )
        }
    except Exception as e:
        return {"error": str(e)}
    if r.returncode != 0:
        return {"error": f"Helper falhou:\n{r.stderr.strip()}"}
    try:
        return json.loads(r.stdout.strip().splitlines()[-1])
    except Exception:
        return {"error": f"Saída inválida do helper:\n{r.stdout[:500]}"}
