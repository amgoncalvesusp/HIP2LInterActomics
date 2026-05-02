"""Cached post-analysis loaders for the Results tab.

Prefer JSON artifacts generated during the LUNA run; fall back to the
legacy subprocess helpers only when those artifacts are missing.
"""
from __future__ import annotations

import json
import subprocess
from pathlib import Path

from . import analysis_helper
from .env_manager import python_process_env

IFP_SUFFIX_TO_TYPE = {"E": "EIFP", "H": "HIFP", "F": "FIFP"}

_FP_SESSION_SCRIPT = r"""
import gzip, json, pickle, re, sys
from pathlib import Path

from luna.interaction.fp.view import ShellViewer
from luna.mol.entry import Entry, MolFileEntry


def _safe_name(value):
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value)).strip("_") or "entry"


def _restore_entry(meta):
    kind = str(meta.get("kind") or "")
    if kind == "MolFileEntry":
        preferred_type = meta.get("mol_obj_type") or "rdkit"
        candidate_types = [preferred_type]
        if preferred_type != "rdkit":
            candidate_types.append("rdkit")
        if preferred_type != "openbabel":
            candidate_types.append("openbabel")
        errors = []
        for mol_obj_type in candidate_types:
            try:
                entry = MolFileEntry.from_mol_file(
                    meta["pdb_id"],
                    meta["mol_id"],
                    meta["mol_file"],
                    bool(meta.get("is_multimol_file", False)),
                    mol_file_ext=meta.get("mol_file_ext"),
                    mol_obj_type=mol_obj_type,
                    autoload=False,
                    overwrite_mol_name=bool(meta.get("overwrite_mol_name", False)),
                    sep=meta.get("sep") or ":",
                )
                entry.get_biopython_structure()
                return entry
            except Exception as exc:
                errors.append(f"{mol_obj_type}: {type(exc).__name__}")
        raise RuntimeError("Nao foi possivel restaurar MolFileEntry: " + "; ".join(errors))
    return Entry.from_string(
        meta["entry_str"],
        is_hetatm=bool(meta.get("is_hetatm", True)),
        sep=meta.get("sep") or ":",
    )


workdir = Path(sys.argv[1])
ifp_type = sys.argv[2]
entry_name = sys.argv[3]
feature_id = int(sys.argv[4])
output_path = Path(sys.argv[5])

suffix = {"EIFP": "E", "HIFP": "H", "FIFP": "F"}.get(ifp_type)
if suffix is None:
    print(json.dumps({"error": f"Tipo IFP inv\u00e1lido: {ifp_type}"}))
    sys.exit(0)

payload_path = workdir / "results" / "fingerprints" / "_shells" / suffix / f"{_safe_name(entry_name)}.pkl.gz"
if not payload_path.exists():
    print(json.dumps({"error": f"Shell artifact n\u00e3o encontrado: {payload_path}"}))
    sys.exit(0)

with gzip.open(payload_path, "rb") as fh:
    payload = pickle.load(fh)

if "entry_meta" in payload and "feature_shells" in payload:
    entry = _restore_entry(payload["entry_meta"])
    shells = payload.get("feature_shells", {}).get(str(feature_id)) or []
else:
    sm = payload["shell_manager"]
    entry = payload["entry"]
    ifp_length = int(payload["ifp_length"])
    ifp_count = bool(payload["ifp_count"])
    unique_shells = not ifp_count
    fp = sm.to_fingerprint(
        fold_to_length=ifp_length,
        unique_shells=unique_shells,
        count_fp=ifp_count,
    )

    recovered = list(sm.trace_back_feature(feature_id, fp, unique_shells=unique_shells))
    shells = []
    for _ori_feature, found_shells in recovered:
        shells.extend(found_shells)

if not shells:
    print(json.dumps({"error": f"Nenhum shell foi encontrado para o fingerprint {feature_id} em {entry_name}."}))
    sys.exit(0)

output_path.parent.mkdir(parents=True, exist_ok=True)
viewer = ShellViewer()
viewer.new_session([(entry, shells, payload["pdb_dir"])], str(output_path))
print(json.dumps({"ok": True, "output": str(output_path), "shells": len(shells)}))
"""

_FP_DETAIL_SCRIPT = r"""
import gzip, json, pickle, sys
from collections import Counter, defaultdict
from pathlib import Path

WATER_RESIDUES = {"HOH", "WAT", "TIP", "SOL", "T3P", "H2O", "OH2", "DOD"}


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


def _residue_label_from_group(atm_grp):
    atoms = list(getattr(atm_grp, "atoms", []) or [])
    labels = Counter()
    for atom in atoms:
        residue = getattr(atom, "parent", None)
        if residue is None:
            continue
        chain = "?"
        full_id = getattr(residue, "full_id", None)
        if isinstance(full_id, (list, tuple)) and len(full_id) >= 3:
            chain = str(full_id[2] or "?").strip() or "?"
        resname = str(getattr(residue, "resname", None) or getattr(residue, "get_resname", lambda: "")()).strip() or "UNK"
        resid = getattr(residue, "id", None)
        seq = ""
        if isinstance(resid, (list, tuple)) and len(resid) >= 3:
            seq = f"{resid[1]}{str(resid[2] or '').strip()}".strip()
        elif resid is not None:
            seq = str(resid).strip()
        if not seq:
            continue
        labels[f"{chain}/{resname}/{seq}"] += 1
    if not labels:
        return ""
    return labels.most_common(1)[0][0]


def _interaction_name(interaction):
    value = getattr(interaction, "type", None)
    if value is None:
        return type(interaction).__name__
    return str(value).strip() or type(interaction).__name__


def _interaction_residue_group(src_grp, trgt_grp, src_role, trgt_role):
    pair = {src_role, trgt_role}
    if pair == {"ligand", "protein"}:
        return src_grp if src_role == "protein" else trgt_grp
    if pair == {"protein", "water"}:
        return src_grp if src_role == "protein" else trgt_grp
    if pair == {"ligand", "water"}:
        return src_grp if src_role == "water" else trgt_grp
    return None


workdir = Path(sys.argv[1])
ifp_type = sys.argv[2]
output_path = Path(sys.argv[3])
suffix = {"EIFP": "E", "HIFP": "H", "FIFP": "F"}.get(ifp_type)
if suffix is None:
    print(json.dumps({"error": f"Tipo IFP invalido: {ifp_type}"}))
    sys.exit(0)

shell_dir = workdir / "results" / "fingerprints" / "_shells" / suffix
if not shell_dir.exists():
    print(json.dumps({"error": f"Diretorio de shells nao encontrado: {shell_dir}"}))
    sys.exit(0)

feature_details = {}
files = sorted(shell_dir.glob("*.pkl.gz"))
for payload_path in files:
    with gzip.open(payload_path, "rb") as fh:
        payload = pickle.load(fh)
    entry_meta = payload.get("entry_meta") or {}
    entry_name = str(entry_meta.get("entry_str") or payload_path.stem)
    feature_shells = payload.get("feature_shells") or {}
    for feature_key, shells in feature_shells.items():
        feature_id = str(feature_key)
        detail = feature_details.setdefault(
            feature_id,
            {
                "interaction_counts": {},
                "residue_counts": {},
                "pair_counts": {},
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
            },
        )
        entry_info["shell_count"] += len(shells)
        for shell in list(shells or []):
            interactions = list(getattr(shell, "interactions", []) or [])
            for interaction in interactions:
                src_grp = getattr(interaction, "src_grp", None)
                trgt_grp = getattr(interaction, "trgt_grp", None)
                src_role = _group_role(src_grp)
                trgt_role = _group_role(trgt_grp)
                residue_group = _interaction_residue_group(src_grp, trgt_grp, src_role, trgt_role)
                if residue_group is None:
                    continue
                interaction_name = _interaction_name(interaction)
                residue_label = _residue_label_from_group(residue_group)
                pair_key = interaction_name if not residue_label else f"{interaction_name}||{residue_label}"

                detail["interaction_counts"][interaction_name] = detail["interaction_counts"].get(interaction_name, 0) + 1
                entry_info["interaction_counts"][interaction_name] = entry_info["interaction_counts"].get(interaction_name, 0) + 1

                if residue_label:
                    detail["residue_counts"][residue_label] = detail["residue_counts"].get(residue_label, 0) + 1
                    entry_info["residue_counts"][residue_label] = entry_info["residue_counts"].get(residue_label, 0) + 1

                detail["pair_counts"][pair_key] = detail["pair_counts"].get(pair_key, 0) + 1
                entry_info["pair_counts"][pair_key] = entry_info["pair_counts"].get(pair_key, 0) + 1

payload = {
    "ifp_type": ifp_type,
    "feature_details": feature_details,
}
output_path.parent.mkdir(parents=True, exist_ok=True)
output_path.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
print(json.dumps({"ok": True, "output": str(output_path), "features": len(feature_details), "entries": len(files)}))
"""

_FP_DASHBOARD_SCRIPT = r"""
import json, sys
from pathlib import Path

package_root = Path(sys.argv[1])
workdir = Path(sys.argv[2])
ifp_type = str(sys.argv[3]).upper()
labels_csv = sys.argv[4]
labels_id_column = sys.argv[5]
labels_column = sys.argv[6]
algorithm_preference = sys.argv[7]
task_kind_preference = sys.argv[8]
use_otsu_threshold = str(sys.argv[9]).strip().lower() in {"1", "true", "yes", "sim"}

sys.path.insert(0, str(package_root))

from luna_gui.core.results_analysis import build_fp_analysis_dashboard, load_fp_analysis_artifacts

artifacts = load_fp_analysis_artifacts(workdir)
artifact = artifacts.get(ifp_type)
if artifact is None:
    print(json.dumps({"error": f"Artefato de fingerprint nao encontrado para {ifp_type}."}))
    sys.exit(0)

dashboard = build_fp_analysis_dashboard(
    workdir,
    artifact,
    labels_csv=labels_csv or None,
    labels_id_column=labels_id_column or None,
    labels_column=labels_column or None,
    algorithm_preference=algorithm_preference or "gradient_boosting",
    task_kind_preference=task_kind_preference or None,
    use_otsu_threshold=use_otsu_threshold,
)
print(json.dumps(dashboard, ensure_ascii=False))
"""

_PSE_FILTER_SCRIPT = r"""
import ast, configparser, fnmatch, gzip, json, pickle, re, shutil, sys
from pathlib import Path

from luna.interaction.view import InteractionViewer


def _safe_name(value):
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value)).strip("_") or "filtered"


def _load_project(workdir):
    candidates = sorted(Path(workdir).glob("project_v*.pkl.gz"), reverse=True)
    if not candidates:
        return None, f"Nenhum project_v*.pkl.gz encontrado em {workdir}."
    errors = []
    for candidate in candidates:
        try:
            with gzip.open(candidate, "rb") as fh:
                return pickle.load(fh), ""
        except Exception as exc:
            errors.append(f"{candidate.name}: {type(exc).__name__}: {exc}")
    return None, "Falha ao reabrir o projeto salvo pelo LUNA (" + "; ".join(errors[:3]) + ")."


def _load_residue_matrix(workdir):
    path = Path(workdir) / "results" / "residue_matrix.json"
    if not path.exists():
        return None
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return None


def _parse_rules(cfg_path):
    parser = configparser.ConfigParser()
    parser.optionxform = str
    parser.read(cfg_path, encoding="utf-8")
    rules = {}
    for section in parser.sections():
        accept_all = parser.get(section, "accept_all", fallback="False").strip().lower() == "true"
        raw_accept_only = parser.get(section, "accept_only", fallback="").strip()
        selectors = []
        if raw_accept_only:
            try:
                parsed = ast.literal_eval(raw_accept_only)
                if isinstance(parsed, (list, tuple)):
                    selectors = [str(item).strip() for item in parsed if str(item).strip()]
                else:
                    selectors = [str(parsed).strip()]
            except Exception:
                selectors = [part.strip().strip("'\"") for part in raw_accept_only.split(",") if part.strip()]
        rules[section.lower()] = {
            "name": section,
            "accept_all": accept_all,
            "selectors": selectors,
        }
    return rules


def _interaction_name(interaction):
    value = getattr(interaction, "type", None)
    if value is None:
        return type(interaction).__name__
    return str(value).strip() or type(interaction).__name__


def _atoms_from_group(group):
    return list(getattr(group, "atoms", []) or [])


def _atom_parts(atom):
    residue = getattr(atom, "parent", None)
    chain = "*"
    resname = "*"
    seq = "*"
    if residue is not None:
        full_id = getattr(residue, "full_id", None)
        if isinstance(full_id, (list, tuple)) and len(full_id) >= 3:
            chain = str(full_id[2] or "*").strip() or "*"
        resname = str(getattr(residue, "resname", None) or getattr(residue, "get_resname", lambda: "")()).strip() or "*"
        resid = getattr(residue, "id", None)
        if isinstance(resid, (list, tuple)) and len(resid) >= 2:
            seq = f"{resid[1]}{str(resid[2] or '').strip()}".strip() or "*"
        elif resid is not None:
            seq = str(resid).strip() or "*"
    atom_name = str(getattr(atom, "name", None) or getattr(atom, "get_name", lambda: "")()).strip() or "*"
    return [chain.upper(), resname.upper(), seq.upper(), atom_name.upper()]


def _selector_matches_atom(selector, atom):
    parts = [part.strip() or "*" for part in str(selector).split("/")]
    parts = (parts + ["*"] * 4)[:4]
    for pattern, value in zip(parts, _atom_parts(atom)):
        if pattern == "*":
            continue
        if not fnmatch.fnmatchcase(value.upper(), pattern.upper()):
            return False
    return True


def _selector_matches_residue_label(selector, residue_label):
    parts = [part.strip() or "*" for part in str(selector).split("/")]
    parts = (parts + ["*"] * 4)[:4]
    residue_parts = [part.strip() or "*" for part in str(residue_label).split("/")]
    residue_parts = (residue_parts + ["*"] * 3)[:3]
    for pattern, value in zip(parts[:3], residue_parts):
        if pattern == "*":
            continue
        if not fnmatch.fnmatchcase(value.upper(), pattern.upper()):
            return False
    return True


def _interaction_matches_rule(interaction, rule):
    if rule.get("accept_all"):
        return True
    selectors = list(rule.get("selectors") or [])
    if not selectors:
        return False
    atoms = []
    atoms.extend(_atoms_from_group(getattr(interaction, "src_grp", None)))
    atoms.extend(_atoms_from_group(getattr(interaction, "trgt_grp", None)))
    return any(_selector_matches_atom(selector, atom) for selector in selectors for atom in atoms)


def _iter_interactions(manager):
    try:
        return list(manager)
    except Exception:
        pass
    for attr in ("interactions", "interactions_list", "_interactions"):
        value = getattr(manager, attr, None)
        if value is not None:
            try:
                return list(value)
            except Exception:
                pass
    return []


def _find_existing_pse_files(workdir):
    files = []
    for folder in (Path(workdir) / "results" / "pse", Path(workdir) / "pse"):
        if folder.exists() and folder.is_dir():
            files.extend(folder.glob("*.pse"))
    return sorted(files)


def _pse_for_entry(pse_files, entry_name):
    safe_variants = {
        _safe_name(entry_name).lower(),
        str(entry_name).replace(":", "_").replace("/", "_").lower(),
        Path(str(entry_name)).stem.lower(),
    }
    for path in pse_files:
        if path.stem.lower() in safe_variants:
            return path
    normalized_entry = _safe_name(entry_name).lower()
    for path in pse_files:
        stem = path.stem.lower()
        if normalized_entry in stem or stem in normalized_entry:
            return path
    return None


def _entry_matches_cached_rules(matrix_artifact, entry_index, rules):
    entries = list(matrix_artifact.get("entries") or [])
    residues = list(matrix_artifact.get("residues") or [])
    matrices = matrix_artifact.get("matrix") or {}
    wildcard = rules.get("*")
    matched = 0

    for interaction_name, rows in matrices.items():
        rule = rules.get(str(interaction_name).lower()) or wildcard
        if not rule:
            continue
        try:
            row = list(rows[entry_index])
        except Exception:
            continue
        positive_residues = [
            residue
            for residue, value in zip(residues, row)
            if float(value or 0.0) > 0.0
        ]
        if not positive_residues:
            continue
        if rule.get("accept_all"):
            matched += len(positive_residues)
            continue
        selectors = list(rule.get("selectors") or [])
        if selectors and any(
            _selector_matches_residue_label(selector, residue)
            for selector in selectors
            for residue in positive_residues
        ):
            matched += 1
    return matched


def _copy_cached_pse_sessions(workdir, rules, output_dir, original_error):
    matrix_artifact = _load_residue_matrix(workdir)
    if not matrix_artifact:
        return {
            "error": original_error
            + " Também não foi encontrado results/residue_matrix.json para aplicar filtragem por arquivo.",
        }
    pse_files = _find_existing_pse_files(workdir)
    if not pse_files:
        return {
            "error": original_error
            + " Também não foram encontradas sessões .pse existentes para copiar.",
        }

    output_dir.mkdir(parents=True, exist_ok=True)
    created = 0
    matched_interactions = 0
    warnings = [
        original_error,
        "Fallback usado: o projeto salvo não reabriu, então a GUI filtrou arquivos .pse existentes usando residue_matrix.json. "
        "As sessões copiadas preservam o conteúdo visual original; a filtragem seleciona quais arquivos entram na subpasta.",
    ]
    for entry_index, entry_name in enumerate(list(matrix_artifact.get("entries") or [])):
        matched = _entry_matches_cached_rules(matrix_artifact, entry_index, rules)
        if matched <= 0:
            continue
        source = _pse_for_entry(pse_files, entry_name)
        if source is None:
            warnings.append(f"PSE não encontrado para {entry_name}.")
            continue
        target = output_dir / source.name
        shutil.copy2(source, target)
        created += 1
        matched_interactions += matched

    return {
        "ok": True,
        "fallback": True,
        "output_dir": str(output_dir),
        "created": created,
        "matched_interactions": matched_interactions,
        "warnings": warnings[:20],
    }


workdir = Path(sys.argv[1])
cfg_path = Path(sys.argv[2])
output_dir = Path(sys.argv[3])

if not cfg_path.exists():
    print(json.dumps({"error": f"Arquivo de regras nao encontrado: {cfg_path}"}))
    sys.exit(0)

rules = _parse_rules(cfg_path)
if not rules:
    print(json.dumps({"error": "O arquivo .cfg nao contem regras."}))
    sys.exit(0)

project, error = _load_project(workdir)
if error:
    print(json.dumps(_copy_cached_pse_sessions(workdir, rules, output_dir, error), ensure_ascii=False))
    sys.exit(0)

output_dir.mkdir(parents=True, exist_ok=True)
created = 0
matched_interactions = 0
warnings = []
wildcard = rules.get("*")
viewer = InteractionViewer(show_hydrop_surface=False)

for entry in list(getattr(project, "entries", []) or []):
    try:
        er = project.get_entry_results(entry)
        selected = []
        for interaction in _iter_interactions(er.interactions_mngr):
            name = _interaction_name(interaction)
            rule = rules.get(name.lower()) or wildcard
            if rule and _interaction_matches_rule(interaction, rule):
                selected.append(interaction)
        if not selected:
            continue
        matched_interactions += len(selected)
        safe_entry = _safe_name(entry.to_string() if hasattr(entry, "to_string") else str(entry))
        pse_path = output_dir / f"{safe_entry}.pse"
        pdb_file = getattr(entry, "pdb_file", "")
        if not pdb_file:
            pdb_id = getattr(entry, "pdb_id", "")
            pdb_file = workdir / "pdbs" / f"{pdb_id}.pdb"
        viewer.new_session([(entry, selected, str(pdb_file))], str(pse_path))
        created += 1
    except Exception as exc:
        warnings.append(f"{entry}: {type(exc).__name__}: {exc}")

print(json.dumps({
    "ok": True,
    "output_dir": str(output_dir),
    "created": created,
    "matched_interactions": matched_interactions,
    "warnings": warnings[:20],
}, ensure_ascii=False))
"""


def run_analysis(py_exe: str, workdir: str, timeout: int = 300) -> dict:
    wd = Path(workdir)
    cached = _load_cached_json(wd / "results" / "analysis_summary.json")
    if cached is not None:
        return cached
    return analysis_helper.run_analysis(py_exe, workdir, timeout=timeout)


def run_residue_matrix(py_exe: str, workdir: str, timeout: int = 600) -> dict:
    wd = Path(workdir)
    cached = _load_cached_json(wd / "results" / "residue_matrix.json")
    if cached is not None:
        return cached
    return analysis_helper.run_residue_matrix(py_exe, workdir, timeout=timeout)


def generate_fp_session(
    py_exe: str,
    workdir: str,
    ifp_type: str,
    entry_name: str,
    feature_id: int,
    output_path: str,
    timeout: int = 600,
) -> dict:
    """Generate a PyMOL session for one fingerprint feature using cached shells."""
    if not py_exe or not Path(py_exe).exists():
        return {"error": "Python do luna-env n\u00e3o encontrado."}
    if not Path(workdir).exists():
        return {"error": f"Workdir n\u00e3o existe: {workdir}"}

    try:
        result = subprocess.run(
            [
                py_exe,
                "-c",
                _FP_SESSION_SCRIPT,
                workdir,
                ifp_type,
                entry_name,
                str(int(feature_id)),
                output_path,
            ],
            capture_output=True,
            text=True,
            timeout=timeout,
            env=python_process_env(py_exe),
        )
    except subprocess.TimeoutExpired:
        return {"error": "A gera\u00e7\u00e3o da sess\u00e3o de fingerprint excedeu o tempo limite."}
    except Exception as exc:
        return {"error": str(exc)}

    if result.returncode != 0:
        return {"error": f"Helper falhou:\n{result.stderr.strip()}"}
    try:
        return json.loads(result.stdout.strip().splitlines()[-1])
    except Exception:
        return {"error": f"Sa\u00edda inv\u00e1lida do helper:\n{result.stdout[:500]}"}


def run_fp_detail_analysis(
    py_exe: str,
    workdir: str,
    ifp_type: str,
    timeout: int = 900,
) -> dict:
    """Generate or load cached interaction/residue summaries for one IFP type."""
    if not py_exe or not Path(py_exe).exists():
        return {"error": "Python do luna-env n\u00e3o encontrado."}
    wd = Path(workdir)
    if not wd.exists():
        return {"error": f"Workdir n\u00e3o existe: {workdir}"}

    suffix = {value: key for key, value in IFP_SUFFIX_TO_TYPE.items()}.get(str(ifp_type).upper())
    if suffix is None:
        return {"error": f"Tipo IFP inv\u00e1lido: {ifp_type}"}

    output_path = wd / "results" / "fingerprints" / f"fp_detail_{suffix}.json"
    cached = _load_cached_json(output_path)
    if cached is not None:
        return cached

    try:
        result = subprocess.run(
            [
                py_exe,
                "-c",
                _FP_DETAIL_SCRIPT,
                str(wd),
                str(ifp_type).upper(),
                str(output_path),
            ],
            capture_output=True,
            text=True,
            timeout=timeout,
            env=python_process_env(py_exe),
        )
    except subprocess.TimeoutExpired:
        return {"error": "A gera\u00e7\u00e3o do resumo de fingerprints excedeu o tempo limite."}
    except Exception as exc:
        return {"error": str(exc)}

    if result.returncode != 0:
        return {"error": f"Helper falhou:\n{result.stderr.strip()}"}
    try:
        parsed = json.loads(result.stdout.strip().splitlines()[-1])
    except Exception:
        parsed = {"ok": False, "error": f"Sa\u00edda inv\u00e1lida do helper:\n{result.stdout[:500]}"}
    if "error" in parsed:
        return parsed
    return _load_cached_json(output_path) or parsed


def run_fp_dashboard_analysis(
    py_exe: str,
    workdir: str,
    ifp_type: str,
    labels_csv: str = "",
    labels_id_column: str = "",
    labels_column: str = "",
    algorithm_preference: str = "gradient_boosting",
    task_kind_preference: str = "",
    use_otsu_threshold: bool = False,
    timeout: int = 900,
) -> dict:
    if not py_exe or not Path(py_exe).exists():
        return {"error": "Python do luna-env não encontrado."}
    wd = Path(workdir)
    if not wd.exists():
        return {"error": f"Workdir não existe: {workdir}"}

    package_root = Path(__file__).resolve().parents[2]
    try:
        result = subprocess.run(
            [
                py_exe,
                "-c",
                _FP_DASHBOARD_SCRIPT,
                str(package_root),
                str(wd),
                str(ifp_type).upper(),
                str(labels_csv or ""),
                str(labels_id_column or ""),
                str(labels_column or ""),
                str(algorithm_preference or "gradient_boosting"),
                str(task_kind_preference or ""),
                "1" if use_otsu_threshold else "0",
            ],
            capture_output=True,
            text=True,
            timeout=timeout,
            env=python_process_env(py_exe),
        )
    except subprocess.TimeoutExpired:
        return {"error": "A geração do dashboard de fingerprints excedeu o tempo limite."}
    except Exception as exc:
        return {"error": str(exc)}

    if result.returncode != 0:
        detail = (result.stderr.strip() or result.stdout.strip() or "sem saída de erro")
        return {"error": f"Helper falhou (exit {result.returncode:#x}):\n{detail[:2000]}"}
    try:
        return json.loads(result.stdout.strip().splitlines()[-1])
    except Exception:
        return {"error": f"Saída inválida do helper:\n{result.stdout[:500]}"}


def generate_filtered_pse_sessions(
    py_exe: str,
    workdir: str,
    binding_modes_cfg: str,
    output_dir: str,
    timeout: int = 1800,
) -> dict:
    """Generate filtered PyMOL sessions from a saved LUNA project."""
    if not py_exe or not Path(py_exe).exists():
        return {"error": "Python do luna-env não encontrado."}
    if not Path(workdir).exists():
        return {"error": f"Workdir não existe: {workdir}"}
    if not Path(binding_modes_cfg).exists():
        return {"error": f"Arquivo .cfg não existe: {binding_modes_cfg}"}

    try:
        result = subprocess.run(
            [
                py_exe,
                "-c",
                _PSE_FILTER_SCRIPT,
                str(workdir),
                str(binding_modes_cfg),
                str(output_dir),
            ],
            capture_output=True,
            text=True,
            timeout=timeout,
            env=python_process_env(py_exe),
        )
    except subprocess.TimeoutExpired:
        return {"error": "A filtragem de sessões PyMOL excedeu o tempo limite."}
    except Exception as exc:
        return {"error": str(exc)}

    if result.returncode != 0:
        detail = (result.stderr.strip() or result.stdout.strip() or "sem saída de erro")
        return {"error": f"Helper falhou (exit {result.returncode:#x}):\n{detail[:2000]}"}
    try:
        return json.loads(result.stdout.strip().splitlines()[-1])
    except Exception:
        return {"error": f"Saída inválida do helper:\n{result.stdout[:500]}"}


def _load_cached_json(path: Path) -> dict | None:
    if not path.exists():
        return None
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return None
