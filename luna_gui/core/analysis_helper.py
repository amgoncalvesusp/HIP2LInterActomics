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
import sys, os, json, glob
workdir = sys.argv[1]
out = {"interaction_types": [], "residues": [], "entries": [], "matrix": {}, "errors": []}
try:
    from luna.projects import EntryResults
    from luna.analysis.residues import generate_residue_matrix
    from luna.interaction.calc import InteractionsManager
except Exception as e:
    print(json.dumps({"error": "luna import failed: %s" % e}))
    sys.exit(0)

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

    out["interaction_types"] = interaction_types
    out["residues"] = residues
    out["entries"] = entries_in_order
    out["matrix"] = matrix
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
        return json.loads(r.stdout.strip().splitlines()[-1])
    except Exception:
        return {"error": f"Saída inválida do helper:\n{r.stdout[:500]}"}


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
