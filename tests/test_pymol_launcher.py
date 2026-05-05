from __future__ import annotations

import os
import sys
from pathlib import Path

from luna_gui.core import pymol_launcher


def test_find_pymol_prefers_luna_env_windows_layout(tmp_path, monkeypatch):
    monkeypatch.setattr(pymol_launcher.shutil, "which", lambda _name: None)
    prefix = tmp_path / "luna-env"
    scripts = prefix / "Scripts"
    scripts.mkdir(parents=True)
    py_exe = prefix / "python.exe"
    py_exe.write_text("", encoding="utf-8")
    pymol_exe = scripts / "pymol.exe"
    pymol_exe.write_text("", encoding="utf-8")

    assert pymol_launcher.find_pymol_executable(py_exe) == pymol_exe


def test_find_pymol_prefers_luna_env_unix_layout(tmp_path, monkeypatch):
    monkeypatch.setattr(pymol_launcher.shutil, "which", lambda _name: None)
    prefix = tmp_path / "luna-env"
    bin_dir = prefix / "bin"
    bin_dir.mkdir(parents=True)
    py_exe = bin_dir / "python"
    py_exe.write_text("", encoding="utf-8")
    pymol_exe = bin_dir / "pymol"
    pymol_exe.write_text("", encoding="utf-8")

    assert pymol_launcher.find_pymol_executable(py_exe) == pymol_exe


def test_launch_pse_uses_luna_env_for_env_owned_pymol(tmp_path, monkeypatch):
    monkeypatch.setattr(pymol_launcher.shutil, "which", lambda _name: None)
    prefix = tmp_path / "luna-env"
    if sys.platform == "win32":
        bin_dir = prefix / "Scripts"
        py_exe = prefix / "python.exe"
        pymol_exe = bin_dir / "pymol.exe"
    else:
        bin_dir = prefix / "bin"
        py_exe = bin_dir / "python"
        pymol_exe = bin_dir / "pymol"
    bin_dir.mkdir(parents=True)
    py_exe.write_text("", encoding="utf-8")
    pymol_exe.write_text("", encoding="utf-8")
    pse_path = tmp_path / "session.pse"
    pse_path.write_text("pse", encoding="utf-8")

    calls: list[tuple[list[str], dict[str, str] | None]] = []

    class DummyProcess:
        pass

    def fake_popen(args, env=None):
        calls.append((args, env))
        return DummyProcess()

    monkeypatch.setattr(pymol_launcher, "_detached_popen", fake_popen)

    launcher = pymol_launcher.launch_pse_session(pse_path, py_exe)

    assert launcher == f"PyMOL: {pymol_exe}"
    assert len(calls) == 1
    assert calls[0][0] == [str(pymol_exe), str(pse_path)]
    env = calls[0][1]
    assert env is not None
    path_parts = [os.path.normcase(part) for part in env["PATH"].split(os.pathsep)]
    assert os.path.normcase(str(bin_dir)) in path_parts
