"""Run LUNA's post-analysis API inside the luna-env via subprocess.

The GUI itself runs in a clean Python with no LUNA installed, so we
serialize the helper script as a string and feed it to `python -c`
through the env's interpreter.
"""
from __future__ import annotations

import json
import subprocess
from pathlib import Path

# Helper script executed inside luna-env. It walks the workdir, loads each
# EntryResults pickle, and aggregates interaction-type counts and per-residue
# counts. Output is a single JSON document on stdout.
HELPER_SCRIPT = r"""
import sys, os, json, glob
workdir = sys.argv[1]
out = {"entries": 0, "interaction_counts": {}, "residue_counts": {}, "errors": []}
try:
    from luna.projects import EntryResults
except Exception as e:
    print(json.dumps({"error": "luna import failed: %s" % e}))
    sys.exit(0)
try:
    from luna.analysis.summary import count_interaction_types
except Exception:
    count_interaction_types = None

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
    # Tally per-residue
    for it in interactions:
        try:
            for grp in (it.src_grp, it.trgt_grp):
                for atm in getattr(grp, "atoms", []) or []:
                    res = atm.parent
                    name = "%s/%s/%s" % (res.parent.id, res.resname, res.id[1])
                    out["residue_counts"][name] = out["residue_counts"].get(name, 0) + 1
        except Exception:
            pass
print(json.dumps(out))
"""


RESIDUE_MATRIX_SCRIPT = r"""
import sys, os, json, glob, re
workdir = sys.argv[1]
out = {
    "interaction_types": [],
    "residues": [],
    "ligand_atoms": [],
    "entries": [],
    "matrix": {},
    "ligand_atom_matrix": {},
    "errors": [],
}
try:
    from luna.projects import EntryResults
    from luna.analysis.residues import generate_residue_matrix
    from luna.interaction.calc import InteractionsManager
except Exception as e:
    print(json.dumps({"error": "luna import failed: %s" % e}))
    sys.exit(0)

AA_RESIDUES = {
    "ALA", "ARG", "ASN", "ASP", "ASH", "CYS", "CYM", "CYX", "GLN", "GLU", "GLH",
    "GLY", "HIS", "HID", "HIE", "HIP", "ILE", "LEU", "LYS", "LYN", "MET", "PHE",
    "PRO", "SER", "THR", "TRP", "TYR", "VAL",
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
    try:
        if atm_grp.has_hetatm() and _group_has_ligand_residue(atm_grp):
            return "ligand"
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


def _interaction_ligand_atom_labels(interaction):
    labels = set()
    for group in (getattr(interaction, "src_grp", None), getattr(interaction, "trgt_grp", None)):
        if _group_role(group) != "ligand":
            continue
        for atom in getattr(group, "atoms", []) or []:
            label = _atom_label(atom)
            if label:
                labels.add(label)
    return labels

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

managers = []
entry_names = []
ligand_atoms = set()
ligand_atom_values = {}
for f in files:
    try:
        er = EntryResults.load(f)
    except Exception:
        continue
    im = getattr(er, "interactions_mngr", None)
    if im is None: continue
    managers.append(im)
    name = ""
    try:
        name = er.entry.to_string()
    except Exception:
        name = os.path.basename(f)
    entry_names.append(name)
    for interaction in list(getattr(im, "interactions", []) or []):
        labels = _interaction_ligand_atom_labels(interaction)
        if not labels:
            continue
        itype = _interaction_type_name(interaction)
        by_entry = ligand_atom_values.setdefault(itype, {}).setdefault(name, {})
        for label in labels:
            ligand_atoms.add(label)
            by_entry[label] = by_entry.get(label, 0) + 1

if not managers:
    print(json.dumps(out)); sys.exit(0)

try:
    df = generate_residue_matrix(managers, by_interaction=True)
except Exception as e:
    out["errors"].append("generate_residue_matrix: %s" % e)
    print(json.dumps(out)); sys.exit(0)

# df.index is a MultiIndex (entry, interaction); df.columns are residue labels
try:
    rows = []
    interaction_types = sorted({idx[1] for idx in df.index})
    residues = list(df.columns)
    # Build {interaction_type: {residue: [values-per-entry]}} where entry order is stable
    entries_in_order = []
    seen_e = set()
    for (ent, _inter) in df.index:
        if ent not in seen_e:
            seen_e.add(ent); entries_in_order.append(ent)

    matrix = {}
    for it in interaction_types:
        sub = df.loc[[i for i in df.index if i[1] == it]]
        # Reindex rows so missing (entry, it) pairs become zeros
        sub_entries = [i[0] for i in sub.index]
        per_entry = {e: [0.0] * len(residues) for e in entries_in_order}
        for (e, _), vals in zip(sub.index, sub.values):
            per_entry[e] = [float(v) for v in vals]
        matrix[it] = [per_entry[e] for e in entries_in_order]

    ligand_atom_labels = sorted(ligand_atoms, key=_atom_sort_key)
    ligand_atom_matrix = {}
    for it in sorted(ligand_atom_values):
        by_entry = ligand_atom_values.get(it, {})
        ligand_atom_matrix[it] = [
            [float(by_entry.get(e, {}).get(label, 0.0)) for label in ligand_atom_labels]
            for e in entries_in_order
        ]
        if it not in interaction_types:
            interaction_types.append(it)

    out["interaction_types"] = interaction_types
    out["residues"] = residues
    out["ligand_atoms"] = ligand_atom_labels
    out["entries"] = entries_in_order
    out["matrix"] = matrix
    out["ligand_atom_matrix"] = ligand_atom_matrix
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


def run_analysis(py_exe: str, workdir: str, timeout: int = 300) -> dict:
    """Execute the helper inside luna-env. Returns parsed JSON or {'error': ...}."""
    if not py_exe or not Path(py_exe).exists():
        return {"error": "Python do luna-env não encontrado."}
    if not Path(workdir).exists():
        return {"error": f"Workdir não existe: {workdir}"}
    try:
        r = subprocess.run(
            [py_exe, "-c", HELPER_SCRIPT, workdir],
            capture_output=True, text=True, timeout=timeout,
        )
    except subprocess.TimeoutExpired:
        return {"error": "Análise excedeu o tempo limite."}
    except Exception as e:
        return {"error": str(e)}
    if r.returncode != 0:
        return {"error": f"Helper falhou:\n{r.stderr.strip()}"}
    try:
        return json.loads(r.stdout.strip().splitlines()[-1])
    except Exception:
        return {"error": f"Saída inválida do helper:\n{r.stdout[:500]}"}
