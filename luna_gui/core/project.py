"""Project state — what the user configured in the GUI for one analysis run."""
from __future__ import annotations

from dataclasses import dataclass, field, asdict
from pathlib import Path
import json
import re

PROJECT_FILENAME = ".luna_gui.json"
HISTORY_FILE = Path.home() / ".luna_gui_history.json"

IFP_ALL = "ALL"
IFP_SINGLE_TYPES = ("EIFP", "HIFP", "FIFP")
IFP_ALL_TYPES = ("HIFP", "EIFP", "FIFP")
IFP_SUFFIXES = {"HIFP": "H", "EIFP": "E", "FIFP": "F"}

_RELOCATABLE_CONFIG_FIELDS = (
    "protein_file",
    "ligand_file",
    "ifp_output",
    "ifp_seed_file",
    "fp_labels_csv",
    "sim_matrix_output",
    "pse_path",
    "binding_modes_cfg",
    "fork_from",
    "interaction_config_file",
)
_RELOCATABLE_PARAMS_FIELDS = (
    "workdir",
    "pdb_dir",
    "lig_file",
    "fork_from",
    "ifp_output",
    "ifp_seed_file",
    "sim_matrix_output",
    "pse_path",
    "interaction_config_file",
)


def _normal_path_text(value: str | Path) -> str:
    """Normalize separators for safe lexical path-boundary comparisons."""
    text = str(value).strip().replace("\\", "/")
    while len(text) > 1 and text.endswith("/"):
        text = text[:-1]
    return text


def _path_casefold(value: str) -> str:
    text = _normal_path_text(value)
    return text.casefold() if re.match(r"^(?:[A-Za-z]:/|//)", text) else text


def _relative_under(path: str, root: str) -> str | None:
    path_text = _path_casefold(path)
    root_text = _path_casefold(root)
    if not path_text or not root_text or path_text == root_text:
        return "" if path_text == root_text else None
    prefix = root_text + "/"
    if not path_text.startswith(prefix):
        return None
    return path_text[len(prefix):]


def _is_absolute_path(value: str) -> bool:
    text = _normal_path_text(value)
    return bool(re.match(r"^(?:[A-Za-z]:/|/|//)", text))


def _join_relocated(root: Path, relative: str) -> Path:
    return root.joinpath(*[part for part in relative.split("/") if part])


def _relocate_path_value(
    value: str,
    old_workdir: str,
    new_workdir: Path,
) -> tuple[str, bool, bool]:
    """Return (value, changed, unresolved) for one path-valued field."""
    raw = str(value or "").strip()
    if not raw or not _is_absolute_path(raw) or not old_workdir:
        return value, False, False

    old_root = _normal_path_text(old_workdir)
    new_root = new_workdir
    roots = [(old_root, new_root)]
    old_parent = str(Path(old_workdir).parent)
    new_parent = new_workdir.parent
    if old_parent and _normal_path_text(old_parent) != old_root:
        roots.append((old_parent, new_parent))

    for source_root, target_root in roots:
        relative = _relative_under(raw, source_root)
        if relative is None:
            continue
        candidate = target_root if not relative else _join_relocated(target_root, relative)
        if candidate.exists():
            return str(candidate), str(candidate) != raw, False
        return value, False, True
    return value, False, False


def is_project_workdir(path: str | Path) -> bool:
    """Return whether a directory contains enough markers to be a saved project."""
    root = Path(path)
    if not root.is_dir() or not (root / PROJECT_FILENAME).is_file():
        return False
    return bool(
        list(root.glob("project_v*.pkl.gz"))
        or (root / "_luna_api_params.json").is_file()
        or (root / "results").is_dir()
    )


def relocation_candidate(config_path: str | Path, cfg: "ProjectConfig") -> tuple[str, Path] | None:
    """Return (saved workdir, config parent) when a project was relocated."""
    selected_path = Path(config_path).resolve()
    if selected_path.name != PROJECT_FILENAME:
        return None
    selected = selected_path.parent
    if not is_project_workdir(selected):
        return None
    saved = str(getattr(cfg, "workdir", "") or "").strip()
    if saved and _path_casefold(saved) == _path_casefold(str(selected)):
        return None
    return saved, selected


def relocate_config_paths(cfg: "ProjectConfig", old_workdir: str, new_workdir: str | Path) -> dict:
    """Rebase existing project-local config paths, leaving external paths intact."""
    new_root = Path(new_workdir).resolve()
    report = {"changed_paths": [], "unresolved_paths": []}
    for field_name in _RELOCATABLE_CONFIG_FIELDS:
        current = getattr(cfg, field_name, "")
        updated, changed, unresolved = _relocate_path_value(current, old_workdir, new_root)
        if changed:
            setattr(cfg, field_name, updated)
            report["changed_paths"].append(field_name)
        elif unresolved:
            report["unresolved_paths"].append(field_name)
    if _path_casefold(str(getattr(cfg, "workdir", "") or "")) != _path_casefold(str(new_root)):
        cfg.workdir = str(new_root)
        report["changed_paths"].append("workdir")
    return report


def relocate_params_file(
    params_path: str | Path,
    old_workdir: str,
    new_workdir: str | Path,
) -> dict:
    """Rebase known path fields in the runner parameters JSON."""
    path = Path(params_path)
    report = {"changed_paths": [], "unresolved_paths": []}
    if not path.is_file():
        return report
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        report["unresolved_paths"].append(str(path))
        return report

    new_root = Path(new_workdir).resolve()
    for field_name in _RELOCATABLE_PARAMS_FIELDS:
        current = data.get(field_name, "")
        updated, changed, unresolved = _relocate_path_value(current, old_workdir, new_root)
        if changed:
            data[field_name] = updated
            report["changed_paths"].append(field_name)
        elif unresolved:
            report["unresolved_paths"].append(field_name)

    specs = data.get("entry_specs")
    if isinstance(specs, list):
        for index, spec in enumerate(specs):
            if not isinstance(spec, dict) or "mol_file" not in spec:
                continue
            updated, changed, unresolved = _relocate_path_value(spec["mol_file"], old_workdir, new_root)
            if changed:
                spec["mol_file"] = updated
                report["changed_paths"].append(f"entry_specs[{index}].mol_file")
            elif unresolved:
                report["unresolved_paths"].append(f"entry_specs[{index}].mol_file")

    for field_name in ("ifp_outputs", "sim_matrix_outputs"):
        outputs = data.get(field_name)
        if not isinstance(outputs, dict):
            continue
        for key, current in list(outputs.items()):
            updated, changed, unresolved = _relocate_path_value(current, old_workdir, new_root)
            if changed:
                outputs[key] = updated
                report["changed_paths"].append(f"{field_name}.{key}")
            elif unresolved:
                report["unresolved_paths"].append(f"{field_name}.{key}")

    if _path_casefold(str(data.get("workdir", "") or "")) != _path_casefold(str(new_root)):
        data["workdir"] = str(new_root)
        report["changed_paths"].append("workdir")
    path.write_text(json.dumps(data, indent=2, ensure_ascii=False), encoding="utf-8")
    return report


@dataclass
class ProjectConfig:
    # GUI
    language: str = "en"          # en | pt | es

    # Inputs
    protein_file: str = ""
    ligand_file: str = ""
    selected_ligands: list[str] = field(default_factory=list)
    trajectory_analysis: bool = False  # entradas representam frames de dinamica molecular

    # Workdir
    workdir: str = ""

    # Analyses
    out_ifp: bool = True
    ifp_type: str = "EIFP"        # EIFP | HIFP | FIFP | ALL
    ifp_levels: int = 2
    ifp_radius: float = 5.73171
    ifp_length: int = 4096
    ifp_bit: bool = False         # False = count fingerprints
    ifp_output: str = ""          # default: <workdir>/results/fingerprints/ifp.csv
    ifp_seed_file: str = ""       # optional file containing the random seed for FP analyses
    fp_labels_csv: str = ""       # optional supervised labels for FP analyses
    fp_labels_id_column: str = "" # optional column name containing ligand IDs
    fp_labels_column: str = ""    # column name containing labels/classes
    fp_label_task: str = "regression"  # regression or classification
    fp_use_otsu_threshold: bool = False # use Otsu fallback when z-score cutoff is unavailable

    sim_matrix: bool = False
    sim_matrix_output: str = ""   # default: <workdir>/sim_matrix.csv

    out_pse: bool = False
    pse_path: str = ""

    add_h: bool = True             # let LUNA add hydrogens
    ph: float = 7.4                # pH used when add_h=True

    filter_binding_modes: bool = False
    binding_modes_cfg: str = ""   # path to .cfg file

    fork_from: str = ""           # path to existing project to fork

    # Execution
    nproc: int = 1
    overwrite: bool = False

    # --- Advanced (used only by the Python-API runner) ---
    # When any of the fields below differ from defaults, the GUI switches
    # from the LUNA CLI to the Python-API runner (core/luna_api_runner.py).
    include_waters: bool = False          # keep HOH; passes ignore_any_h2o=False
    # Treat non-water HETATM residues present in the receptor PDB as protein
    # components (for example cofactors and coordinated metal ions).
    include_protein_heteroatoms: bool = False
    # InteractionCalculator flags
    ic_add_proximal: bool = False
    ic_add_atom_atom: bool = False
    ic_add_dependent_inter: bool = True
    ic_add_h2o_pairs_with_no_target: bool = True
    ic_ignore_self_inter: bool = True
    # DefaultInteractionConfig overrides — {key: value}; empty means "defaults"
    inter_config_overrides: dict = field(default_factory=dict)
    inter_max_distance_cap: float = 0.0
    # Optional full InteractionConfig file (.cfg) used by the Python API runner
    interaction_config_file: str = ""
    # PSE filter — only generate PSE for these interaction types; empty = all
    pse_interaction_types: list[str] = field(default_factory=list)
    # Force the Python-API runner even if no advanced field changed
    force_python_api: bool = False

    def selected_ifp_types(self) -> list[str]:
        """Return the concrete IFP types requested by the UI."""
        if self.ifp_type == IFP_ALL:
            return list(IFP_ALL_TYPES)
        if self.ifp_type in IFP_SINGLE_TYPES:
            return [self.ifp_type]
        return ["EIFP"]

    def uses_python_api(self) -> bool:
        """True when any advanced knob requires the Python-API runner."""
        if self.force_python_api:
            return True
        if self.fork_from:
            return True
        if self.ifp_type == IFP_ALL:
            return True
        if self.ifp_seed_file:
            return True
        if not self.add_h:
            return True
        if abs(float(self.ph) - 7.4) > 1e-9:
            return True
        if self.include_waters:
            return True
        if self.include_protein_heteroatoms:
            return True
        if self.interaction_config_file:
            return True
        if self.inter_config_overrides:
            return True
        if float(self.inter_max_distance_cap or 0.0) > 0.0:
            return True
        if self.pse_interaction_types:
            return True
        # Non-default InteractionCalculator flags
        defaults = dict(
            ic_add_proximal=False, ic_add_atom_atom=False,
            ic_add_dependent_inter=True, ic_add_h2o_pairs_with_no_target=True,
            ic_ignore_self_inter=True,
        )
        for k, v in defaults.items():
            if getattr(self, k) != v:
                return True
        return False

    def to_json(self) -> str:
        return json.dumps(asdict(self), indent=2)

    def save(self, path: str | Path) -> None:
        Path(path).write_text(self.to_json(), encoding="utf-8")

    @classmethod
    def load(cls, path: str | Path) -> "ProjectConfig":
        data = json.loads(Path(path).read_text(encoding="utf-8"))
        return cls(**data)


def save_to_workdir(cfg: ProjectConfig) -> None:
    if not cfg.workdir:
        return
    wd = Path(cfg.workdir)
    wd.mkdir(parents=True, exist_ok=True)
    cfg.save(wd / PROJECT_FILENAME)
    add_to_history(str(wd))


def add_to_history(workdir: str) -> None:
    items = load_history()
    items = [w for w in items if w != workdir]
    items.insert(0, workdir)
    items = items[:50]
    HISTORY_FILE.write_text(json.dumps(items, indent=2), encoding="utf-8")


def load_history() -> list[str]:
    if not HISTORY_FILE.exists():
        return []
    try:
        return json.loads(HISTORY_FILE.read_text(encoding="utf-8"))
    except Exception:
        return []


def remove_from_history(workdir: str) -> None:
    items = [item for item in load_history() if item != workdir]
    HISTORY_FILE.write_text(json.dumps(items, indent=2), encoding="utf-8")


def clear_history() -> None:
    HISTORY_FILE.write_text("[]\n", encoding="utf-8")


def resolve_ifp_output_paths(cfg: ProjectConfig) -> dict[str, str]:
    """Return the concrete IFP CSV path(s) implied by the current config."""
    types = cfg.selected_ifp_types()
    if not types:
        return {}

    root = Path(cfg.workdir or ".")
    if len(types) == 1:
        output = cfg.ifp_output.strip() if cfg.ifp_output else str(root / "results" / "fingerprints" / "ifp.csv")
        return {types[0]: output}

    base_dir = Path(cfg.ifp_output).parent if cfg.ifp_output else root / "results" / "fingerprints"
    return {
        ifp_type: str(base_dir / f"ifp_{IFP_SUFFIXES[ifp_type]}.csv")
        for ifp_type in types
    }


def resolve_sim_matrix_output_paths(cfg: ProjectConfig) -> dict[str, str]:
    """Return the concrete similarity-matrix path(s) implied by the current config."""
    if not cfg.sim_matrix:
        return {}

    types = cfg.selected_ifp_types()
    root = Path(cfg.workdir or ".")
    if len(types) == 1:
        if cfg.sim_matrix_output:
            output = cfg.sim_matrix_output.strip()
        else:
            output = str(root / f"sim_matrix_{IFP_SUFFIXES[types[0]]}.csv")
        return {types[0]: output}

    base_dir = Path(cfg.sim_matrix_output).parent if cfg.sim_matrix_output else root
    return {
        ifp_type: str(base_dir / f"sim_matrix_{IFP_SUFFIXES[ifp_type]}.csv")
        for ifp_type in types
    }
