"""Conda environment management for LUNA.

The GUI runs in a regular Python interpreter; the heavy LUNA stack
(RDKit / OpenBabel / PyMOL / Biopython) lives in a dedicated conda env.
"""
from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

ENV_NAME = "luna-env"
GUI_ENV_NAME = "luna-gui"

# conda-forge packages required by LUNA (see LUNA docs §1.2.1).
# Biopython is pinned per LUNA's requirement.
CONDA_PACKAGES = [
    "python=3.9",
    "pip",
    "rdkit",
    "openbabel",
    "pymol-open-source",
    "biopython=1.79",  # LUNA uses Bio.Data.SCOPData, removed in >=1.80
    "numpy",
    "pandas",
    "scipy",
    "scikit-learn",
    "matplotlib",
    "seaborn",
    "networkx",
]

# Pure-pip dependencies with version constraints for LUNA compatibility.
# These are installed via pip because conda-forge's versions may be incompatible.
PIP_PREREQS = [
    "pdbecif",        # Required by LUNA at setup.py import time
    "mmh3<4",         # LUNA passes numpy arrays to mmh3.hash; mmh3>=4 rejects them
    "xopen",
    "colorlog",
]

RUNTIME_MODULES = [
    ("sklearn", "scikit-learn"),
    ("numpy", "numpy"),
    ("scipy", "scipy"),
    ("matplotlib", "matplotlib"),
]


# Path to the post-install patch script (runs inside luna-env)
_LUNA_PATCH_SCRIPT = str(Path(__file__).parent / "_luna_patch.py")


def find_conda() -> str | None:
    """Return the best conda executable for the current GUI session."""
    for candidate in _iter_conda_candidates():
        if candidate.exists():
            return str(candidate)
    return None


def _iter_conda_candidates() -> list[Path]:
    """Yield candidate conda executables in priority order."""
    candidates: list[Path] = []

    override = os.environ.get("HIP2LINTERACTOMICS_GUI_CONDA") or os.environ.get("LUNA_GUI_CONDA")
    if override:
        candidates.extend(_conda_candidates_from_path(Path(override)))

    # Highest priority: conda installation related to the current GUI session.
    for env_name in ("CONDA_EXE", "MAMBA_EXE"):
        value = os.environ.get(env_name)
        if value:
            candidates.extend(_conda_candidates_from_path(Path(value)))

    related_paths = [Path(sys.executable), Path(sys.prefix)]
    for env_name in ("CONDA_PREFIX", "MAMBA_ROOT_PREFIX"):
        value = os.environ.get(env_name)
        if value:
            related_paths.append(Path(value))

    for path in related_paths:
        candidates.extend(_conda_candidates_from_path(path))

    # Then try PATH-resolved executables.
    for name in ("conda", "mamba"):
        exe = shutil.which(name)
        if exe:
            candidates.extend(_conda_candidates_from_path(Path(exe)))

    # Finally, look through common installation roots.
    for root in _common_conda_roots():
        candidates.extend(_conda_executables_for_root(root))

    return _dedupe_paths(candidates)


def _conda_candidates_from_path(path: Path) -> list[Path]:
    """Derive possible conda executables from a related file or directory."""
    candidates: list[Path] = []

    if not path:
        return candidates

    roots: list[Path] = []
    if path.suffix.lower() in (".exe", ".bat", ".cmd"):
        parent_name = path.parent.name.lower()
        if parent_name in ("scripts", "condabin", "bin"):
            roots.append(path.parent.parent)
        roots.append(path.parent)
    else:
        roots.append(path)

    # If the path points inside ".../envs/<env_name>/...", prefer that base root.
    try:
        parts = list(path.parts)
        lower_parts = [part.lower() for part in parts]
        env_idx = lower_parts.index("envs")
    except ValueError:
        env_idx = -1
    if env_idx > 0:
        roots.append(Path(*parts[:env_idx]))

    # If the path already looks like "<base>/envs/luna-gui", derive "<base>".
    if path.name.lower() in ("python", "python.exe") and path.parent.name.lower() == GUI_ENV_NAME:
        envs_dir = path.parent.parent
        if envs_dir.name.lower() == "envs":
            roots.append(envs_dir.parent)
    elif path.name.lower() == GUI_ENV_NAME and path.parent.name.lower() == "envs":
        roots.append(path.parent.parent)

    expanded: list[Path] = []
    for root in roots:
        expanded.extend(_conda_executables_for_root(root))
    return expanded


def _conda_executables_for_root(root: Path) -> list[Path]:
    root = Path(root)
    if sys.platform == "win32":
        return [
            root / "Scripts" / "conda.exe",
            root / "Scripts" / "mamba.exe",
        ]
    return [
        root / "bin" / "conda",
        root / "bin" / "mamba",
    ]


def _common_conda_roots() -> list[Path]:
    home = Path.home()
    roots = [
        home / ".conda",
        home / "miniconda3",
        home / "anaconda3",
        home / "miniforge3",
        home / "mambaforge",
    ]
    if sys.platform == "win32":
        localapp = os.environ.get("LOCALAPPDATA")
        program_data = os.environ.get("ProgramData", "C:/ProgramData")
        if localapp:
            roots.extend([
                Path(localapp) / "miniconda3",
                Path(localapp) / "anaconda3",
                Path(localapp) / "miniforge3",
                Path(localapp) / "mambaforge",
            ])
        roots.extend([
            Path(program_data) / "miniconda3",
            Path(program_data) / "Miniconda3",
            Path(program_data) / "anaconda3",
            Path(program_data) / "Anaconda3",
            Path(program_data) / "miniforge3",
            Path(program_data) / "mambaforge",
        ])
    else:
        roots.extend([
            Path("/opt/miniconda3"),
            Path("/opt/anaconda3"),
            Path("/opt/miniforge3"),
            Path("/opt/mambaforge"),
        ])
    return _dedupe_paths(roots)


def _dedupe_paths(paths: list[Path]) -> list[Path]:
    seen: set[str] = set()
    unique: list[Path] = []
    for path in paths:
        norm = os.path.normcase(os.path.normpath(str(path)))
        if norm in seen:
            continue
        seen.add(norm)
        unique.append(path)
    return unique


def conda_root(conda: str | Path) -> Path:
    """Return the installation root for a conda executable."""
    conda_path = Path(conda)
    parent_name = conda_path.parent.name.lower()
    if parent_name in ("scripts", "bin", "condabin"):
        return conda_path.parent.parent
    return conda_path.parent


def conda_process_env(conda: str | Path) -> dict[str, str]:
    """Return a sanitized environment for running Conda commands."""
    env = os.environ.copy()
    root = conda_root(conda)

    root_paths: list[str]
    if sys.platform == "win32":
        root_paths = [
            str(root),
            str(root / "Library" / "mingw-w64" / "bin"),
            str(root / "Library" / "usr" / "bin"),
            str(root / "Library" / "bin"),
            str(root / "Scripts"),
            str(root / "bin"),
        ]
    else:
        root_paths = [str(root / "bin")]

    original_path = [p for p in env.get("PATH", "").split(os.pathsep) if p]
    merged_path: list[str] = []
    seen: set[str] = set()
    for path in root_paths + original_path:
        norm = os.path.normcase(os.path.normpath(path))
        if norm in seen:
            continue
        seen.add(norm)
        merged_path.append(path)
    env["PATH"] = os.pathsep.join(merged_path)

    for key in (
        "CONDA_PREFIX",
        "CONDA_DEFAULT_ENV",
        "CONDA_PROMPT_MODIFIER",
        "CONDA_SHLVL",
        "CONDA_PYTHON_EXE",
        "_CE_CONDA",
        "_CE_M",
        "PYTHONHOME",
        "PYTHONPATH",
        "PYTHONUSERBASE",
    ):
        env.pop(key, None)
    env["PYTHONNOUSERSITE"] = "1"

    conda_str = str(conda)
    if Path(conda_str).name.lower().startswith("mamba"):
        env["MAMBA_EXE"] = conda_str
        env.pop("CONDA_EXE", None)
    else:
        env["CONDA_EXE"] = conda_str
        env.pop("MAMBA_EXE", None)

    return env


def python_prefix(py: str | Path) -> Path:
    """Return the environment prefix for a Python executable."""
    py_path = Path(py)
    if sys.platform == "win32":
        if py_path.name.lower() == "python.exe":
            return py_path.parent
    elif py_path.parent.name == "bin" and py_path.name == "python":
        return py_path.parent.parent
    return py_path.parent


def python_process_env(py: str | Path) -> dict[str, str]:
    """Return a sanitized environment for running an env-owned Python."""
    env = os.environ.copy()
    prefix = python_prefix(py)

    if sys.platform == "win32":
        preferred_paths = [
            str(prefix),
            str(prefix / "Library" / "mingw-w64" / "bin"),
            str(prefix / "Library" / "usr" / "bin"),
            str(prefix / "Library" / "bin"),
            str(prefix / "Scripts"),
        ]
    else:
        preferred_paths = [str(prefix / "bin")]

    original_path = [p for p in env.get("PATH", "").split(os.pathsep) if p]
    merged_path: list[str] = []
    seen: set[str] = set()
    for path in preferred_paths + original_path:
        norm = os.path.normcase(os.path.normpath(path))
        if norm in seen:
            continue
        seen.add(norm)
        merged_path.append(path)
    env["PATH"] = os.pathsep.join(merged_path)

    for key in ("PYTHONHOME", "PYTHONPATH", "PYTHONUSERBASE"):
        env.pop(key, None)
    env["PYTHONNOUSERSITE"] = "1"
    return env


def conda_info(conda: str) -> dict[str, object]:
    """Return `conda info --json` parsed output, or an empty dict on failure."""
    try:
        out = subprocess.check_output(
            [conda, "info", "--json"],
            text=True,
            stderr=subprocess.STDOUT,
            env=conda_process_env(conda),
        )
    except Exception:
        return {}
    try:
        info = json.loads(out)
    except json.JSONDecodeError:
        return {}
    return info if isinstance(info, dict) else {}


def _listed_env_prefixes(conda: str) -> list[Path]:
    info = conda_info(conda)
    envs = info.get("envs")
    if not isinstance(envs, list):
        return []
    return [Path(prefix) for prefix in envs if isinstance(prefix, str) and prefix]


def _default_env_prefix(conda: str, name: str = ENV_NAME) -> Path:
    info = conda_info(conda)
    envs_dirs = info.get("envs_dirs")
    if isinstance(envs_dirs, list):
        for envs_dir in envs_dirs:
            if isinstance(envs_dir, str) and envs_dir:
                return Path(envs_dir) / name
    return conda_root(conda) / "envs" / name


def env_prefix(conda: str, name: str = ENV_NAME) -> Path:
    """Return the absolute prefix where the target env lives or should be created."""
    lowered_name = name.casefold()
    for prefix in _listed_env_prefixes(conda):
        if prefix.name.casefold() == lowered_name:
            return prefix
    return _default_env_prefix(conda, name)


def _env_python_path(prefix: Path) -> Path:
    if sys.platform == "win32":
        return prefix / "python.exe"
    return prefix / "bin" / "python"


def env_is_valid(prefix: str | Path) -> bool:
    """Return True when the target prefix looks like a usable Conda env."""
    prefix = Path(prefix)
    if not prefix.exists():
        return False

    py = _env_python_path(prefix)
    meta_dir = prefix / "conda-meta"
    if not py.exists() or not meta_dir.is_dir():
        return False

    if (meta_dir / "history").exists():
        return True

    try:
        return any(path.suffix.lower() == ".json" for path in meta_dir.iterdir())
    except OSError:
        return False


def env_is_partial(conda: str, name: str = ENV_NAME) -> bool:
    prefix = env_prefix(conda, name)
    return prefix.exists() and not env_is_valid(prefix)


def _is_expected_env_prefix(prefix: Path, name: str = ENV_NAME) -> bool:
    try:
        resolved = prefix.resolve()
    except OSError:
        resolved = prefix
    return (
        resolved.name.casefold() == name.casefold()
        and resolved.parent.name.casefold() == "envs"
    )


def cleanup_partial_env(conda: str, name: str = ENV_NAME) -> Path | None:
    """Remove an incomplete env directory so Conda can recreate it cleanly."""
    prefix = env_prefix(conda, name)
    if not prefix.exists() or env_is_valid(prefix):
        return None
    if not _is_expected_env_prefix(prefix, name):
        raise ValueError(f"Refusing to remove unexpected env prefix: {prefix}")
    shutil.rmtree(prefix)
    return prefix


def list_envs(conda: str) -> list[str]:
    return [prefix.name for prefix in _listed_env_prefixes(conda)]


def env_exists(conda: str, name: str = ENV_NAME) -> bool:
    return env_is_valid(env_prefix(conda, name))


def env_python(conda: str, name: str = ENV_NAME) -> Path | None:
    """Return path to the python executable inside the env."""
    py = _env_python_path(env_prefix(conda, name))
    return py if py.exists() else None


def luna_installed(py: Path) -> bool:
    try:
        r = subprocess.run(
            [str(py), "-c", "import luna; print(luna.__file__)"],
            capture_output=True, text=True, timeout=15, env=python_process_env(py),
        )
        return r.returncode == 0
    except Exception:
        return False


def luna_run_py_path(py: Path) -> Path | None:
    """Locate luna/run.py inside the env's site-packages."""
    try:
        r = subprocess.run(
            [str(py), "-c", "import luna, os; print(os.path.dirname(luna.__file__))"],
            capture_output=True, text=True, timeout=15, env=python_process_env(py),
        )
        if r.returncode != 0:
            return None
        run_py = Path(r.stdout.strip()) / "run.py"
        return run_py if run_py.exists() else None
    except Exception:
        return None


def missing_runtime_packages(py: Path) -> list[str]:
    """Return conda package names missing from the runtime Python."""
    script = (
        "import importlib.util, json; "
        f"mods={json.dumps(RUNTIME_MODULES)}; "
        "print(json.dumps([pkg for mod, pkg in mods if importlib.util.find_spec(mod) is None]))"
    )
    try:
        r = subprocess.run(
            [str(py), "-c", script],
            capture_output=True,
            text=True,
            timeout=20,
            env=python_process_env(py),
        )
        if r.returncode != 0:
            return [pkg for _mod, pkg in RUNTIME_MODULES]
        return list(json.loads(r.stdout.strip() or "[]"))
    except Exception:
        return [pkg for _mod, pkg in RUNTIME_MODULES]


def install_commands(conda: str, name: str = ENV_NAME) -> list[list[str]]:
    """Return the sequence of commands to create env + install LUNA.

    The caller (UI) is responsible for streaming each command's output.
    """
    prefix = env_prefix(conda, name)
    create_or_update_cmd = "install" if env_is_valid(prefix) else "create"
    create = [
        conda,
        create_or_update_cmd,
        "-p",
        str(prefix),
        "--override-channels",
        "-c",
        "conda-forge",
        "-y",
        *CONDA_PACKAGES,
    ]
    # LUNA's setup.py imports the package to read its version, which requires
    # several pure-pip deps to be present BEFORE `pip install luna`.
    pip_prereqs = [
        conda, "run", "-p", str(prefix), "python", "-s", "-m", "pip", "install",
        *PIP_PREREQS,
    ]
    pip_luna = [
        conda, "run", "-p", str(prefix), "python", "-s", "-m", "pip", "install",
        "--no-build-isolation", "-U", "luna",
    ]
    # Patch LUNA source after install to fix the Windows np.int_ overflow bug.
    patch = [
        conda, "run", "-p", str(prefix), "python", "-s", _LUNA_PATCH_SCRIPT,
    ]
    return [create, pip_prereqs, pip_luna, patch]


def miniconda_download_url() -> str:
    if sys.platform == "win32":
        return "https://repo.anaconda.com/miniconda/Miniconda3-latest-Windows-x86_64.exe"
    if sys.platform == "darwin":
        return "https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh"
    return "https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh"
