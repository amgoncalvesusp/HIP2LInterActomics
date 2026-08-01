"""Headless HIP²LInterActomics/LUNA command-line entry point.

Settings may come from JSON (or a legacy Python literal dictionary), direct
command-line options, or both. Direct options override the configuration file.
This module never imports the Qt UI package or creates a GUI event loop.
"""
from __future__ import annotations

import argparse
import ast
import json
import os
import subprocess
import sys
from dataclasses import fields
from pathlib import Path
from typing import Any


def _configure_headless_environment() -> None:
    """Force non-interactive rendering before scientific modules are imported."""
    os.environ.setdefault("MPLBACKEND", "Agg")
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")


_configure_headless_environment()

from luna_gui.core import env_manager as em
from luna_gui.core import ligand_io, luna_api_runner, luna_runner, terminal_results
from luna_gui.core.process_control import TerminationController, signal_exit_code
from luna_gui.core.project import ProjectConfig, save_to_workdir
from luna_gui.core.runtime_resources import detect_cpu_allocation, effective_nproc


PROJECT_FIELDS = {field.name for field in fields(ProjectConfig)}
BOOL_FIELDS = {
    "trajectory_analysis",
    "out_ifp",
    "ifp_bit",
    "fp_use_otsu_threshold",
    "sim_matrix",
    "out_pse",
    "add_h",
    "filter_binding_modes",
    "overwrite",
    "include_waters",
    "ic_add_proximal",
    "ic_add_atom_atom",
    "ic_add_dependent_inter",
    "ic_add_h2o_pairs_with_no_target",
    "ic_ignore_self_inter",
    "force_python_api",
}
INT_FIELDS = {"ifp_levels", "ifp_length", "nproc"}
FLOAT_FIELDS = {"ifp_radius", "ph", "inter_max_distance_cap"}
LIST_FIELDS = {"selected_ligands", "pse_interaction_types"}
DICT_FIELDS = {"inter_config_overrides"}
TERMINAL_KEYS = {
    "allow_hydrogen_warnings",
    "conda",
    "conda_exe",
    "dry_run",
    "entries_file",
    "env_name",
    "luna_env",
    "luna_python",
    "py_exe",
    "python",
    "python_exe",
    "run_py",
    "selected_ligands_file",
    "terminal_cluster_count",
    "terminal_cluster_max_entries",
    "terminal_cluster_method",
    "terminal_interactive_max_entries",
    "terminal_matrix_max_entries",
    "terminal_results",
    "fp_session",
}


EXAMPLE_CONFIG = """{
  "protein_file": "/caminho/para/proteinas_ou_receptor.pdb",
  "ligand_file": "/caminho/para/ligantes_ou_pasta",
  "selected_ligands": "ALL",
  "trajectory_analysis": false,
  "workdir": "/caminho/para/workdir",

  "out_ifp": true,
  "ifp_type": "ALL",
  "ifp_levels": 6,
  "ifp_radius": 2.0,
  "ifp_length": 4096,
  "ifp_bit": true,
  "sim_matrix": true,
  "out_pse": false,

  "include_waters": true,
  "add_h": true,
  "ph": 7.4,
  "nproc": 1,
  "overwrite": true,

  "fp_labels_csv": "",
  "fp_labels_id_column": "",
  "fp_labels_column": "",
  "fp_label_task": "regression",
  "fp_use_otsu_threshold": true,

  "python_exe": "",
  "conda_exe": "",
  "allow_hydrogen_warnings": false,

  "terminal_results": true,
  "terminal_cluster_method": "average",
  "terminal_cluster_count": 4,
  "terminal_matrix_max_entries": 5000,
  "terminal_cluster_max_entries": 5000,
  "terminal_interactive_max_entries": 2500,
  "fp_session": {
    "ifp_type": "EIFP",
    "entry_name": "",
    "feature_id": "",
    "output_path": ""
  }
}
"""


def _read_dict_file(path: Path) -> dict[str, Any]:
    text = path.read_text(encoding="utf-8", errors="replace").strip()
    if not text:
        raise ValueError(f"Arquivo de configuracao vazio: {path}")

    try:
        data = json.loads(text)
    except Exception:
        data = None
    if isinstance(data, dict):
        return data

    try:
        module = ast.parse(text, filename=str(path), mode="exec")
    except SyntaxError as exc:
        raise ValueError(f"Configuracao invalida em {path}: {exc}") from exc

    if len(module.body) == 1 and isinstance(module.body[0], ast.Expr):
        value = ast.literal_eval(module.body[0].value)
        if isinstance(value, dict):
            return value

    for node in module.body:
        if not isinstance(node, ast.Assign):
            continue
        for target in node.targets:
            if isinstance(target, ast.Name) and target.id in {"config", "cfg", "settings"}:
                value = ast.literal_eval(node.value)
                if isinstance(value, dict):
                    return value

    raise ValueError(
        "O arquivo deve conter JSON, um dicionario Python literal, ou uma atribuicao "
        "config = {...}."
    )


def _as_bool(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    if isinstance(value, (int, float)):
        return bool(value)
    text = str(value).strip().lower()
    if text in {"1", "true", "t", "yes", "y", "sim", "s", "on"}:
        return True
    if text in {"0", "false", "f", "no", "n", "nao", "não", "off", ""}:
        return False
    raise ValueError(f"Valor booleano invalido: {value!r}")


def _as_list(value: Any) -> list[str]:
    if value is None:
        return []
    if isinstance(value, (list, tuple, set)):
        return [str(item).strip() for item in value if str(item).strip()]
    text = str(value).strip()
    if not text:
        return []
    if text.upper() == "ALL":
        return ["ALL"]
    separators = [",", "\n", ";"]
    parts = [text]
    for separator in separators:
        if separator in text:
            parts = text.replace(";", ",").replace("\n", ",").split(",")
            break
    return [part.strip() for part in parts if part.strip()]


def _normalize_project_data(raw: dict[str, Any]) -> tuple[dict[str, Any], dict[str, Any]]:
    if isinstance(raw.get("project"), dict):
        project_data = dict(raw["project"])
        terminal_data = {key: value for key, value in raw.items() if key != "project"}
        if isinstance(raw.get("terminal"), dict):
            terminal_data.update(raw["terminal"])
        if isinstance(raw.get("execution"), dict):
            terminal_data.update(raw["execution"])
    else:
        project_data = {
            key: value
            for key, value in raw.items()
            if key in PROJECT_FIELDS
        }
        terminal_data = {
            key: value
            for key, value in raw.items()
            if key not in PROJECT_FIELDS
        }

    unknown = sorted(
        key
        for key in terminal_data
        if key not in TERMINAL_KEYS and key not in {"terminal", "execution"}
    )
    if unknown:
        raise ValueError(
            "Chaves desconhecidas no arquivo de configuracao: " + ", ".join(unknown)
        )

    normalized: dict[str, Any] = {}
    for key, value in project_data.items():
        if key not in PROJECT_FIELDS:
            raise ValueError(f"Campo de ProjectConfig desconhecido: {key}")
        if key in BOOL_FIELDS:
            normalized[key] = _as_bool(value)
        elif key in INT_FIELDS:
            normalized[key] = int(value)
        elif key in FLOAT_FIELDS:
            normalized[key] = float(value)
        elif key in LIST_FIELDS:
            normalized[key] = _as_list(value)
        elif key in DICT_FIELDS:
            if value in (None, ""):
                normalized[key] = {}
            elif isinstance(value, dict):
                normalized[key] = value
            else:
                raise ValueError(f"O campo {key} deve ser um dicionario.")
        else:
            normalized[key] = str(value) if isinstance(value, Path) else value
    return normalized, terminal_data


def _load_entries_from_file(path: str | os.PathLike[str]) -> list[str]:
    entries_path = Path(path)
    if not entries_path.exists():
        raise FileNotFoundError(f"Arquivo de entradas nao encontrado: {entries_path}")
    return [
        line.strip()
        for line in entries_path.read_text(encoding="utf-8", errors="replace").splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ]


def _resolve_selected_ligands(cfg: ProjectConfig, terminal_data: dict[str, Any]) -> None:
    entries_file = terminal_data.get("entries_file") or terminal_data.get("selected_ligands_file")
    if entries_file:
        cfg.selected_ligands = _load_entries_from_file(str(entries_file))
        return

    selected = _as_list(cfg.selected_ligands)
    if selected and selected != ["ALL"]:
        cfg.selected_ligands = selected
        return

    if not cfg.ligand_file:
        cfg.selected_ligands = []
        return
    cfg.selected_ligands = ligand_io.parse_ligand_file(cfg.ligand_file)


def _resolve_luna_python(terminal_data: dict[str, Any]) -> Path:
    explicit = (
        terminal_data.get("python_exe")
        or terminal_data.get("luna_python")
        or terminal_data.get("py_exe")
        or terminal_data.get("python")
    )
    if explicit:
        py = Path(str(explicit)).expanduser()
        if not py.exists():
            raise FileNotFoundError(f"Python informado nao existe: {py}")
        return py

    conda = terminal_data.get("conda_exe") or terminal_data.get("conda")
    if not conda:
        conda = em.find_conda()
    if not conda:
        raise RuntimeError(
            "Conda nao encontrado. Informe 'python_exe' apontando para o python do luna-env "
            "ou 'conda_exe' apontando para conda/mamba."
        )

    env_name = str(
        terminal_data.get("env_name")
        or terminal_data.get("luna_env")
        or em.ENV_NAME
    )
    py = em.env_python(str(conda), env_name)
    if py is None or not py.exists():
        raise RuntimeError(f"Ambiente conda '{env_name}' nao encontrado ou sem python.")
    return py


def _write_entries_and_project(cfg: ProjectConfig) -> Path:
    workdir = Path(cfg.workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    entries_file = workdir / "entries.txt"
    ligand_io.write_entries_file(entries_file, cfg.selected_ligands)
    save_to_workdir(cfg)
    return entries_file


def _run_command(cmd: list[str], env: dict[str, str], log_path: Path) -> int:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    termination = TerminationController()
    with termination.installed():
        with log_path.open(
            "w",
            encoding="utf-8",
            errors="replace",
            buffering=64 * 1024,
        ) as log:
            log.write("=== HIP2LInterActomics terminal ===\n")
            log.write("$ " + " ".join(cmd) + "\n\n")

            process = subprocess.Popen(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                encoding="utf-8",
                errors="replace",
                bufsize=1,
                env=env,
            )
            termination.attach(process)
            try:
                assert process.stdout is not None
                for line in process.stdout:
                    print(line, end="")
                    log.write(line)
                returncode = process.wait()
            finally:
                if process.stdout is not None:
                    process.stdout.close()
                termination.detach(process)

    if termination.received_signal is not None:
        return signal_exit_code(termination.received_signal)
    return returncode


def _effective_settings(
    config_path: Path | None,
    project_overrides: dict[str, Any] | None = None,
    terminal_overrides: dict[str, Any] | None = None,
) -> tuple[dict[str, Any], dict[str, Any]]:
    raw = _read_dict_file(config_path) if config_path is not None else {}
    project_data, terminal_data = _normalize_project_data(raw)
    project_data.update(project_overrides or {})

    supplied_terminal = dict(terminal_overrides or {})
    fp_session = supplied_terminal.pop("fp_session", None)
    terminal_data.update(supplied_terminal)
    if fp_session:
        merged_session = dict(terminal_data.get("fp_session") or {})
        merged_session.update(fp_session)
        terminal_data["fp_session"] = merged_session
    return project_data, terminal_data


def run_from_config(
    config_path: Path | None,
    dry_run: bool = False,
    project_overrides: dict[str, Any] | None = None,
    terminal_overrides: dict[str, Any] | None = None,
) -> int:
    project_data, terminal_data = _effective_settings(
        config_path,
        project_overrides,
        terminal_overrides,
    )
    if dry_run:
        terminal_data["dry_run"] = True

    cfg = ProjectConfig(**project_data)
    _resolve_selected_ligands(cfg, terminal_data)
    cfg.force_python_api = True
    allocation = detect_cpu_allocation()
    selected_nproc = effective_nproc(cfg.nproc)
    safe_nproc = luna_runner.safe_nproc(selected_nproc)
    print(
        f"[terminal] nproc={safe_nproc} "
        f"(resource source: {allocation.source})"
    )
    if selected_nproc != safe_nproc:
        print(
            "[terminal] Aviso: no Windows nativo o LUNA foi limitado para nproc=1. "
            "Use Linux ou WSL2 para paralelismo real."
        )
    cfg.nproc = safe_nproc

    py_exe = _resolve_luna_python(terminal_data)
    if not em.luna_installed(py_exe):
        raise RuntimeError(f"LUNA nao esta instalado no Python informado: {py_exe}")

    run_py = em.luna_run_py_path(py_exe)
    if run_py:
        print(f"[terminal] LUNA run.py: {run_py}")
    else:
        print("[terminal] Aviso: run.py nao encontrado; usando Python API runner.")

    errors = luna_runner.validate(cfg)
    if errors:
        raise ValueError("Configuracao invalida:\n" + "\n".join(errors))

    hydrogen_warnings = luna_api_runner.validate_hydrogen_inputs(cfg)
    if hydrogen_warnings:
        message = "Alerta de hidrogenios (Add_H):\n" + "\n".join(hydrogen_warnings)
        if _as_bool(terminal_data.get("allow_hydrogen_warnings", False)):
            print("[terminal] " + message.replace("\n", "\n[terminal] "))
        else:
            raise ValueError(message + "\nDefina allow_hydrogen_warnings=true para continuar mesmo assim.")

    entry_errors = luna_api_runner.validate_entry_specs(cfg, cfg.selected_ligands)
    if entry_errors:
        raise ValueError("Configuracao invalida:\n" + "\n".join(entry_errors))

    entries_file = _write_entries_and_project(cfg)
    cmd = luna_api_runner.build_api_command(str(py_exe), cfg, cfg.selected_ligands)
    log_path = Path(cfg.workdir) / "hipplinteractomics_terminal.log"

    print(f"[terminal] Configuracao salva em: {Path(cfg.workdir) / '.luna_gui.json'}")
    print(f"[terminal] Entries salvas em: {entries_file}")
    print(f"[terminal] Log terminal: {log_path}")
    if cfg.trajectory_analysis:
        print("[terminal] Modo trajetoria/poses ativo; graficos com frames/poses manterao a ordem temporal/numerica.")
    elif str(cfg.fp_labels_csv or "").strip():
        if str(cfg.fp_label_task or "").strip().lower() == "classification":
            print("[terminal] Graficos da GUI com moleculas serao agrupados pelos rotulos do CSV/TSV.")
        else:
            print("[terminal] Graficos da GUI com moleculas serao ordenados por rotulo em ordem decrescente.")
    if cfg.trajectory_analysis:
        print(f"[terminal] Modo trajetoria/poses ativo; mapa 2D do ligante sera salvo em {Path(cfg.workdir) / 'results' / 'ligand_atom_map.png'} quando possivel.")
    print("[terminal] Comando:")
    print(" ".join(cmd))

    if _as_bool(terminal_data.get("dry_run", False)):
        print("[terminal] dry_run=true; nada foi executado.")
        return 0

    env = em.python_process_env(py_exe)
    env["PYTHONIOENCODING"] = "utf-8"
    env["PYTHONUTF8"] = "1"
    code = _run_command(cmd, env, log_path)
    if code == 0:
        print(f"\n[terminal] Analise concluida. Resultados em: {cfg.workdir}")
        if _as_bool(terminal_data.get("terminal_results", True)):
            manifest = terminal_results.run_terminal_results(cfg, str(py_exe), terminal_data)
            manifest_path = manifest.get("outputs", {}).get("manifest", "")
            print(f"[terminal] Artefatos de Resultados: {manifest_path or 'results/terminal'}")
            if manifest.get("errors"):
                print("[terminal] Alguns artefatos de Resultados falharam; consulte o manifesto.")
    else:
        print(f"\n[terminal] LUNA finalizou com exit code {code}. Veja o log: {log_path}")
    return code


def run_results_from_config(
    config_path: Path | None,
    project_overrides: dict[str, Any] | None = None,
    terminal_overrides: dict[str, Any] | None = None,
) -> int:
    project_data, terminal_data = _effective_settings(
        config_path,
        project_overrides,
        terminal_overrides,
    )
    cfg = ProjectConfig(**project_data)
    if not cfg.workdir:
        raise ValueError("workdir e obrigatorio para --results-only.")

    try:
        py_exe = _resolve_luna_python(terminal_data)
    except Exception as exc:
        py_exe = Path(sys.executable)
        print(f"[terminal] Aviso: Python do luna-env nao resolvido ({exc}); usando caches disponiveis.")

    manifest = terminal_results.run_terminal_results(cfg, str(py_exe), terminal_data)
    manifest_path = manifest.get("outputs", {}).get("manifest", "")
    print(f"[terminal] Artefatos de Resultados: {manifest_path or 'results/terminal'}")
    if manifest.get("errors"):
        print("[terminal] Alguns artefatos falharam; consulte o manifesto.")
    return 0 if manifest.get("outputs") else 1


def _add_boolean_option(
    group: argparse._ArgumentGroup,
    option: str,
    *,
    dest: str | None = None,
    help_text: str,
) -> None:
    group.add_argument(
        option,
        dest=dest,
        action=argparse.BooleanOptionalAction,
        default=None,
        help=help_text,
    )


def _json_object_argument(value: str) -> dict[str, Any]:
    candidate = Path(value).expanduser()
    try:
        text = candidate.read_text(encoding="utf-8") if candidate.is_file() else value
        parsed = json.loads(text)
    except (OSError, json.JSONDecodeError) as exc:
        raise argparse.ArgumentTypeError(
            "expected a JSON object or a path to a JSON object"
        ) from exc
    if not isinstance(parsed, dict):
        raise argparse.ArgumentTypeError("the value must resolve to a JSON object")
    return parsed


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="hipplinteractomics-terminal",
        description=(
            "Run HIP²LInterActomics/LUNA in fully headless mode. Options supplied "
            "on the command line override values from the optional JSON config."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "config_path",
        nargs="?",
        help="JSON or legacy Python-literal configuration file.",
    )
    parser.add_argument(
        "--config",
        dest="config_option",
        metavar="PATH",
        help="Explicit configuration path (alternative to the positional argument).",
    )
    parser.add_argument(
        "--write-template",
        metavar="PATH",
        help="Write a complete JSON configuration template and exit.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Validate and print the generated LUNA command without executing it.",
    )
    parser.add_argument(
        "--results-only",
        action="store_true",
        help="Regenerate cached result artifacts without running LUNA.",
    )

    inputs = parser.add_argument_group("project inputs and outputs")
    inputs.add_argument("--protein-file", metavar="PATH")
    inputs.add_argument("--ligand-file", metavar="PATH")
    inputs.add_argument("--workdir", metavar="PATH")
    inputs.add_argument(
        "--selected-ligand",
        dest="selected_ligands",
        action="append",
        default=None,
        metavar="ID",
        help="Select one ligand; repeat for multiple IDs. Omit to auto-detect all.",
    )
    inputs.add_argument("--entries-file", metavar="PATH")
    inputs.add_argument("--fork-from", metavar="PATH")
    inputs.add_argument("--language", choices=("pt", "en", "es"))
    _add_boolean_option(
        inputs,
        "--trajectory-analysis",
        help_text="Treat entries as ordered trajectory frames or docking poses.",
    )
    _add_boolean_option(
        inputs,
        "--include-waters",
        help_text="Preserve water-mediated interactions.",
    )

    fingerprints = parser.add_argument_group("fingerprints and supervised labels")
    _add_boolean_option(
        fingerprints,
        "--out-ifp",
        help_text="Generate interaction fingerprints.",
    )
    fingerprints.add_argument(
        "--ifp-type",
        choices=("EIFP", "HIFP", "FIFP", "ALL"),
    )
    fingerprints.add_argument("--ifp-levels", type=int, metavar="N")
    fingerprints.add_argument(
        "--ifp-radius",
        type=float,
        metavar="VALUE",
        help="LUNA shell growth ratio/radius.",
    )
    fingerprints.add_argument("--ifp-length", type=int, metavar="BITS")
    fingerprints.add_argument(
        "--ifp-format",
        choices=("binary", "bin", "count", "cnt"),
        help="Convenience alias that sets --ifp-bit/--no-ifp-bit.",
    )
    _add_boolean_option(
        fingerprints,
        "--ifp-bit",
        help_text="Use binary fingerprints; disable for count fingerprints.",
    )
    fingerprints.add_argument("--ifp-output", metavar="PATH")
    fingerprints.add_argument("--ifp-seed-file", metavar="PATH")
    fingerprints.add_argument("--fp-labels-csv", metavar="PATH")
    fingerprints.add_argument("--fp-labels-id-column", metavar="NAME")
    fingerprints.add_argument("--fp-labels-column", metavar="NAME")
    fingerprints.add_argument(
        "--fp-label-task",
        choices=("regression", "classification"),
    )
    _add_boolean_option(
        fingerprints,
        "--fp-use-otsu-threshold",
        help_text="Use Otsu as a fallback feature-importance threshold.",
    )

    analyses = parser.add_argument_group("analyses and LUNA execution")
    _add_boolean_option(
        analyses,
        "--sim-matrix",
        help_text="Generate Tanimoto similarity matrices.",
    )
    analyses.add_argument("--sim-matrix-output", metavar="PATH")
    _add_boolean_option(
        analyses,
        "--out-pse",
        help_text="Generate PyMOL sessions.",
    )
    analyses.add_argument("--pse-path", metavar="PATH")
    analyses.add_argument(
        "--pse-interaction-type",
        dest="pse_interaction_types",
        action="append",
        default=None,
        metavar="TYPE",
        help="Interaction type retained in PSE output; repeat as needed.",
    )
    _add_boolean_option(
        analyses,
        "--add-h",
        help_text="Allow LUNA to add hydrogens.",
    )
    analyses.add_argument("--ph", type=float)
    _add_boolean_option(
        analyses,
        "--filter-binding-modes",
        help_text="Enable the binding-mode filter configuration.",
    )
    analyses.add_argument("--binding-modes-cfg", metavar="PATH")
    analyses.add_argument("--interaction-config-file", metavar="PATH")
    analyses.add_argument(
        "--inter-config-overrides",
        type=_json_object_argument,
        metavar="JSON_OR_PATH",
    )
    analyses.add_argument("--inter-max-distance-cap", type=float)
    analyses.add_argument("--nproc", type=int)
    _add_boolean_option(
        analyses,
        "--overwrite",
        help_text="Allow an existing project workdir to be overwritten.",
    )
    _add_boolean_option(
        analyses,
        "--force-python-api",
        help_text="Force the LUNA Python API runner.",
    )

    advanced = parser.add_argument_group("advanced InteractionCalculator flags")
    for option, help_text in (
        ("--ic-add-proximal", "Enable proximal contacts."),
        ("--ic-add-atom-atom", "Enable atom-atom contacts."),
        ("--ic-add-dependent-inter", "Enable dependent interactions."),
        (
            "--ic-add-h2o-pairs-with-no-target",
            "Include water pairs that do not have a target interaction.",
        ),
        ("--ic-ignore-self-inter", "Ignore self interactions."),
    ):
        _add_boolean_option(advanced, option, help_text=help_text)

    runtime = parser.add_argument_group("headless runtime and result export")
    runtime.add_argument(
        "--python-exe",
        "--luna-python",
        dest="python_exe",
        metavar="PATH",
        help="Python executable inside luna-env.",
    )
    runtime.add_argument(
        "--conda-exe",
        "--conda",
        dest="conda_exe",
        metavar="PATH",
    )
    runtime.add_argument("--env-name", default=None)
    _add_boolean_option(
        runtime,
        "--allow-hydrogen-warnings",
        help_text="Continue after recoverable hydrogen validation warnings.",
    )
    _add_boolean_option(
        runtime,
        "--terminal-results",
        help_text="Export terminal result tables, plots, and manifests.",
    )
    runtime.add_argument("--terminal-cluster-method")
    runtime.add_argument("--terminal-cluster-count", type=int)
    runtime.add_argument("--terminal-matrix-max-entries", type=int)
    runtime.add_argument("--terminal-cluster-max-entries", type=int)
    runtime.add_argument("--terminal-interactive-max-entries", type=int)
    runtime.add_argument("--fp-session-ifp-type", choices=("EIFP", "HIFP", "FIFP"))
    runtime.add_argument("--fp-session-entry-name")
    runtime.add_argument("--fp-session-feature-id")
    runtime.add_argument("--fp-session-output-path", metavar="PATH")
    return parser


_PROJECT_OVERRIDE_NAMES = (
    "language",
    "protein_file",
    "ligand_file",
    "selected_ligands",
    "trajectory_analysis",
    "workdir",
    "out_ifp",
    "ifp_type",
    "ifp_levels",
    "ifp_radius",
    "ifp_length",
    "ifp_bit",
    "ifp_output",
    "ifp_seed_file",
    "fp_labels_csv",
    "fp_labels_id_column",
    "fp_labels_column",
    "fp_label_task",
    "fp_use_otsu_threshold",
    "sim_matrix",
    "sim_matrix_output",
    "out_pse",
    "pse_path",
    "add_h",
    "ph",
    "filter_binding_modes",
    "binding_modes_cfg",
    "fork_from",
    "nproc",
    "overwrite",
    "include_waters",
    "ic_add_proximal",
    "ic_add_atom_atom",
    "ic_add_dependent_inter",
    "ic_add_h2o_pairs_with_no_target",
    "ic_ignore_self_inter",
    "inter_config_overrides",
    "inter_max_distance_cap",
    "interaction_config_file",
    "pse_interaction_types",
    "force_python_api",
)

_TERMINAL_OVERRIDE_NAMES = (
    "entries_file",
    "python_exe",
    "conda_exe",
    "env_name",
    "allow_hydrogen_warnings",
    "terminal_results",
    "terminal_cluster_method",
    "terminal_cluster_count",
    "terminal_matrix_max_entries",
    "terminal_cluster_max_entries",
    "terminal_interactive_max_entries",
)


def _collect_cli_overrides(
    args: argparse.Namespace,
    parser: argparse.ArgumentParser,
) -> tuple[dict[str, Any], dict[str, Any]]:
    project = {
        name: getattr(args, name)
        for name in _PROJECT_OVERRIDE_NAMES
        if getattr(args, name, None) is not None
    }
    terminal = {
        name: getattr(args, name)
        for name in _TERMINAL_OVERRIDE_NAMES
        if getattr(args, name, None) is not None
    }

    if args.ifp_format:
        format_is_binary = args.ifp_format in {"binary", "bin"}
        if args.ifp_bit is not None and args.ifp_bit != format_is_binary:
            parser.error("--ifp-format conflicts with --ifp-bit/--no-ifp-bit")
        project["ifp_bit"] = format_is_binary

    session_names = {
        "ifp_type": "fp_session_ifp_type",
        "entry_name": "fp_session_entry_name",
        "feature_id": "fp_session_feature_id",
        "output_path": "fp_session_output_path",
    }
    fp_session = {
        key: getattr(args, attribute)
        for key, attribute in session_names.items()
        if getattr(args, attribute) is not None
    }
    if fp_session:
        terminal["fp_session"] = fp_session
    return project, terminal


def main(argv: list[str] | None = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)

    if args.write_template:
        output = Path(args.write_template).expanduser()
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(EXAMPLE_CONFIG, encoding="utf-8")
        print(f"Template salvo em: {output}")
        return 0

    if args.config_path and args.config_option:
        parser.error("use either the positional config or --config, not both")
    config_value = args.config_option or args.config_path
    config_path = Path(config_value).expanduser() if config_value else None
    project_overrides, terminal_overrides = _collect_cli_overrides(args, parser)

    if args.results_only and args.dry_run:
        parser.error("--results-only cannot be combined with --dry-run")
    if config_path is None:
        required = ("workdir",) if args.results_only else (
            "protein_file",
            "ligand_file",
            "workdir",
        )
        missing = [name for name in required if not project_overrides.get(name)]
        if missing:
            options = ", ".join("--" + name.replace("_", "-") for name in missing)
            parser.error(
                "without a config file the following options are required: " + options
            )

    try:
        if args.results_only:
            return run_results_from_config(
                config_path,
                project_overrides,
                terminal_overrides,
            )
        return run_from_config(
            config_path,
            dry_run=bool(args.dry_run),
            project_overrides=project_overrides,
            terminal_overrides=terminal_overrides,
        )
    except Exception as exc:
        print(f"[terminal] ERRO: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
