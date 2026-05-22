"""Run HIP2LInterActomics/LUNA from a text configuration file.

The input file may be JSON or a Python literal dictionary. The dictionary keys
are the same fields stored by the GUI in ``.luna_gui.json``. Extra terminal-only
keys are accepted for environment discovery, for example ``python_exe`` or
``conda_exe``.
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

from luna_gui.core import env_manager as em
from luna_gui.core import ligand_io, luna_api_runner, luna_runner
from luna_gui.core.project import ProjectConfig, save_to_workdir


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
  "out_pse": true,

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
  "allow_hydrogen_warnings": false
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
    with log_path.open("w", encoding="utf-8", errors="replace") as log:
        log.write("=== HIP2LInterActomics terminal ===\n")
        log.write("$ " + " ".join(cmd) + "\n\n")
        log.flush()

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
        assert process.stdout is not None
        for line in process.stdout:
            print(line, end="")
            log.write(line)
            log.flush()
        return process.wait()


def run_from_config(config_path: Path, dry_run: bool = False) -> int:
    raw = _read_dict_file(config_path)
    project_data, terminal_data = _normalize_project_data(raw)
    if dry_run:
        terminal_data["dry_run"] = True

    cfg = ProjectConfig(**project_data)
    _resolve_selected_ligands(cfg, terminal_data)
    cfg.force_python_api = True

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
    else:
        print(f"\n[terminal] LUNA finalizou com exit code {code}. Veja o log: {log_path}")
    return code


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Executa HIP2LInterActomics/LUNA sem abrir a interface grafica.",
    )
    parser.add_argument(
        "config",
        nargs="?",
        help="Arquivo .txt/.json contendo um dicionario de configuracao.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Valida a configuracao, escreve os temporarios e mostra o comando sem executar.",
    )
    parser.add_argument(
        "--write-template",
        metavar="ARQUIVO",
        help="Grava um exemplo de configuracao e sai.",
    )
    args = parser.parse_args(argv)

    if args.write_template:
        out = Path(args.write_template)
        out.write_text(EXAMPLE_CONFIG, encoding="utf-8")
        print(f"Template salvo em: {out}")
        return 0

    if not args.config:
        parser.error("informe um arquivo de configuracao ou use --write-template")

    try:
        return run_from_config(Path(args.config), dry_run=bool(args.dry_run))
    except Exception as exc:
        print(f"[terminal] ERRO: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
