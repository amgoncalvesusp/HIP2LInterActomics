"""Conda environment management for LUNA.

The GUI runs in a regular Python interpreter; the heavy LUNA stack
(RDKit / OpenBabel / PyMOL / Biopython) lives in a dedicated conda env.
"""
from __future__ import annotations

import os
import shutil
import subprocess
import sys
from pathlib import Path

ENV_NAME = "luna-env"

# conda-forge packages required by LUNA (see LUNA docs §1.2.1).
# Biopython is pinned per LUNA's requirement.
CONDA_PACKAGES = [
    "python=3.9",
    "rdkit",
    "openbabel",
    "pymol-open-source",
    "biopython=1.79",  # LUNA uses Bio.Data.SCOPData, removed in >=1.80
    "numpy",
    "pandas",
    "scipy",
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


# Path to the post-install patch script (runs inside luna-env)
_LUNA_PATCH_SCRIPT = str(Path(__file__).parent / "_luna_patch.py")


def find_conda() -> str | None:
    """Return path to a conda executable, or None."""
    exe = shutil.which("conda") or shutil.which("mamba")
    if exe:
        return exe
    # Common install locations
    candidates = []
    if sys.platform == "win32":
        home = Path.home()
        candidates += [
            home / "miniconda3" / "Scripts" / "conda.exe",
            home / "anaconda3" / "Scripts" / "conda.exe",
            Path("C:/ProgramData/miniconda3/Scripts/conda.exe"),
            Path("C:/ProgramData/Anaconda3/Scripts/conda.exe"),
        ]
    else:
        home = Path.home()
        candidates += [
            home / "miniconda3" / "bin" / "conda",
            home / "anaconda3" / "bin" / "conda",
            Path("/opt/miniconda3/bin/conda"),
            Path("/opt/anaconda3/bin/conda"),
        ]
    for c in candidates:
        if c.exists():
            return str(c)
    return None


def list_envs(conda: str) -> list[str]:
    try:
        out = subprocess.check_output(
            [conda, "env", "list"], text=True, stderr=subprocess.STDOUT
        )
    except Exception:
        return []
    envs = []
    for line in out.splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        envs.append(line.split()[0])
    return envs


def env_exists(conda: str, name: str = ENV_NAME) -> bool:
    return name in list_envs(conda)


def env_python(conda: str, name: str = ENV_NAME) -> Path | None:
    """Return path to the python executable inside the env."""
    try:
        out = subprocess.check_output(
            [conda, "env", "list"], text=True, stderr=subprocess.STDOUT
        )
    except Exception:
        return None
    for line in out.splitlines():
        parts = line.split()
        if parts and parts[0] == name:
            env_path = Path(parts[-1])
            if sys.platform == "win32":
                py = env_path / "python.exe"
            else:
                py = env_path / "bin" / "python"
            return py if py.exists() else None
    return None


def luna_installed(py: Path) -> bool:
    try:
        r = subprocess.run(
            [str(py), "-c", "import luna; print(luna.__file__)"],
            capture_output=True, text=True, timeout=15,
        )
        return r.returncode == 0
    except Exception:
        return False


def luna_run_py_path(py: Path) -> Path | None:
    """Locate luna/run.py inside the env's site-packages."""
    try:
        r = subprocess.run(
            [str(py), "-c", "import luna, os; print(os.path.dirname(luna.__file__))"],
            capture_output=True, text=True, timeout=15,
        )
        if r.returncode != 0:
            return None
        run_py = Path(r.stdout.strip()) / "run.py"
        return run_py if run_py.exists() else None
    except Exception:
        return None


def install_commands(conda: str, name: str = ENV_NAME) -> list[list[str]]:
    """Return the sequence of commands to create env + install LUNA.

    The caller (UI) is responsible for streaming each command's output.
    """
    create = [conda, "create", "-n", name, "--override-channels",
              "-c", "conda-forge", "-y", *CONDA_PACKAGES]
    # LUNA's setup.py imports the package to read its version, which requires
    # several pure-pip deps to be present BEFORE `pip install luna`.
    pip_prereqs = [
        conda, "run", "-n", name, "python", "-m", "pip", "install",
        *PIP_PREREQS,
    ]
    pip_luna = [
        conda, "run", "-n", name, "python", "-m", "pip", "install",
        "--no-build-isolation", "-U", "luna",
    ]
    # Patch LUNA source after install to fix the Windows np.int_ overflow bug.
    patch = [
        conda, "run", "-n", name, "python", _LUNA_PATCH_SCRIPT,
    ]
    return [create, pip_prereqs, pip_luna, patch]


def miniconda_download_url() -> str:
    if sys.platform == "win32":
        return "https://repo.anaconda.com/miniconda/Miniconda3-latest-Windows-x86_64.exe"
    if sys.platform == "darwin":
        return "https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh"
    return "https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh"
