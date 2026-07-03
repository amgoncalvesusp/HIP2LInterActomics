"""Helpers for opening PyMOL session files from the GUI."""
from __future__ import annotations

import os
import shutil
import subprocess
import sys
from pathlib import Path

from .env_manager import python_prefix, python_process_env


class PymolLaunchError(RuntimeError):
    """Raised when a PyMOL session cannot be opened."""


def _dedupe(paths: list[Path]) -> list[Path]:
    seen: set[str] = set()
    result: list[Path] = []
    for path in paths:
        key = os.path.normcase(os.path.abspath(str(path)))
        if key in seen:
            continue
        seen.add(key)
        result.append(path)
    return result


def iter_pymol_candidates(py_exe: str | Path | None = None) -> list[Path]:
    """Return likely PyMOL executables, preferring the selected LUNA env."""
    candidates: list[Path] = []

    if py_exe:
        py_path = Path(py_exe)
        prefixes = [python_prefix(py_exe)]
        if py_path.parent.name.lower() in {"bin", "scripts"}:
            prefixes.append(py_path.parent.parent)
        for prefix in _dedupe(prefixes):
            candidates.extend(
                [
                    prefix / "Scripts" / "pymol.exe",
                    prefix / "Scripts" / "pymol.bat",
                    prefix / "Scripts" / "pymol.cmd",
                    prefix / "Library" / "bin" / "pymol.exe",
                    prefix / "bin" / "pymol",
                    prefix / "bin" / "pymol.exe",
                ]
            )

    for name in ("pymol", "pymol.exe", "pymol.bat", "pymol.cmd"):
        resolved = shutil.which(name)
        if resolved:
            candidates.append(Path(resolved))

    return _dedupe(candidates)


def find_pymol_executable(py_exe: str | Path | None = None) -> Path | None:
    """Find an installed PyMOL executable."""
    for candidate in iter_pymol_candidates(py_exe):
        try:
            if candidate.exists() and candidate.is_file():
                return candidate
        except OSError:
            continue
    return None


def _is_inside(child: Path, parent: Path) -> bool:
    try:
        child_resolved = child.resolve()
        parent_resolved = parent.resolve()
    except OSError:
        return False
    return child_resolved == parent_resolved or parent_resolved in child_resolved.parents


def _detached_popen(args: list[str], env: dict[str, str] | None = None) -> subprocess.Popen:
    kwargs: dict[str, object] = {"env": env} if env is not None else {}
    if sys.platform == "win32":
        flags = 0
        flags |= getattr(subprocess, "CREATE_NEW_PROCESS_GROUP", 0)
        flags |= getattr(subprocess, "DETACHED_PROCESS", 0)
        if flags:
            kwargs["creationflags"] = flags
    else:
        kwargs["start_new_session"] = True
    return subprocess.Popen(args, **kwargs)


def _association_open(path: Path) -> str:
    if sys.platform == "win32":
        os.startfile(str(path))  # type: ignore[attr-defined]
        return "associacao do Windows"
    if sys.platform == "darwin":
        _detached_popen(["open", str(path)])
        return "associacao do macOS"
    opener = shutil.which("xdg-open")
    if opener:
        _detached_popen([opener, str(path)])
        return "associacao do Linux"
    raise PymolLaunchError(
        "Nao encontrei PyMOL no ambiente do LUNA/PATH e tambem nao encontrei xdg-open."
    )


def launch_pse_session(path: str | Path, py_exe: str | Path | None = None) -> str:
    """Open a PyMOL session/script file, falling back to the OS association."""
    pse_path = Path(path)
    if not pse_path.exists():
        raise PymolLaunchError(f"Arquivo de sessao PyMOL nao encontrado:\n{pse_path}")
    if not pse_path.is_file():
        raise PymolLaunchError(f"O caminho selecionado nao e um arquivo:\n{pse_path}")

    pymol = find_pymol_executable(py_exe)
    if pymol:
        env = None
        if py_exe:
            try:
                prefix = python_prefix(py_exe)
                if _is_inside(pymol, prefix):
                    env = python_process_env(py_exe)
            except Exception:
                env = None
        _detached_popen([str(pymol), str(pse_path)], env=env)
        return f"PyMOL: {pymol}"

    try:
        return _association_open(pse_path)
    except Exception as exc:
        raise PymolLaunchError(
            "Nao consegui abrir a sessao PyMOL. Instale/ative o PyMOL no ambiente "
            "do LUNA ou associe arquivos .pse ao PyMOL no sistema.\n\n"
            f"Arquivo: {pse_path}\nErro original: {exc}"
        ) from exc
