"""Project state — what the user configured in the GUI for one analysis run."""
from __future__ import annotations

from dataclasses import dataclass, field, asdict
from pathlib import Path
import json

PROJECT_FILENAME = ".luna_gui.json"
HISTORY_FILE = Path.home() / ".luna_gui_history.json"


@dataclass
class ProjectConfig:
    # Inputs
    protein_file: str = ""
    ligand_file: str = ""
    selected_ligands: list[str] = field(default_factory=list)

    # Workdir
    workdir: str = ""

    # Analyses
    out_ifp: bool = True
    ifp_type: str = "EIFP"        # EIFP | HIFP | FIFP
    ifp_levels: int = 2
    ifp_radius: float = 5.73171
    ifp_length: int = 4096
    ifp_bit: bool = False         # False = count fingerprints
    ifp_output: str = ""          # default: <workdir>/results/fingerprints/ifp.csv

    sim_matrix: bool = False
    sim_matrix_output: str = ""   # default: <workdir>/sim_matrix.csv

    out_pse: bool = False
    pse_path: str = ""

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
    # InteractionCalculator flags
    ic_add_proximal: bool = False
    ic_add_atom_atom: bool = False
    ic_add_dependent_inter: bool = True
    ic_add_h2o_pairs_with_no_target: bool = True
    ic_ignore_self_inter: bool = True
    # DefaultInteractionConfig overrides — {key: value}; empty means "defaults"
    inter_config_overrides: dict = field(default_factory=dict)
    # PSE filter — only generate PSE for these interaction types; empty = all
    pse_interaction_types: list[str] = field(default_factory=list)
    # Force the Python-API runner even if no advanced field changed
    force_python_api: bool = False

    def uses_python_api(self) -> bool:
        """True when any advanced knob requires the Python-API runner."""
        if self.force_python_api:
            return True
        if self.include_waters:
            return True
        if self.inter_config_overrides:
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
