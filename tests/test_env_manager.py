from __future__ import annotations

import tempfile
import unittest
from pathlib import Path
from unittest import mock

from luna_gui.core import env_manager


class EnvManagerTests(unittest.TestCase):
    def test_env_prefix_uses_conda_envs_dirs_when_env_is_missing(self) -> None:
        conda = r"C:\ProgramData\Anaconda3\Scripts\conda.exe"
        info = {
            "envs": [r"C:\Users\danie\.conda\envs\luna-gui"],
            "envs_dirs": [
                r"C:\Users\danie\.conda\envs",
                r"C:\ProgramData\Anaconda3\envs",
            ],
        }
        with mock.patch.object(env_manager, "conda_info", return_value=info):
            prefix = env_manager.env_prefix(conda)

        self.assertEqual(prefix, Path(r"C:\Users\danie\.conda\envs\luna-env"))

    def test_find_conda_honors_manual_override(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            override_exe = root / "custom" / "Scripts" / "conda.exe"
            override_exe.parent.mkdir(parents=True, exist_ok=True)
            override_exe.write_text("", encoding="utf-8")

            with mock.patch.dict(
                env_manager.os.environ,
                {
                    "LUNA_GUI_CONDA": str(override_exe),
                    "USERPROFILE": str(root),
                },
                clear=True,
            ):
                found = env_manager.find_conda()

            self.assertEqual(found, str(override_exe))

    def test_conda_process_env_prefers_conda_root_bins_and_clears_active_env(self) -> None:
        conda = r"C:\ProgramData\Anaconda3\Scripts\conda.exe"
        original_path = r"C:\Users\danie\.conda\envs\luna-gui;C:\Windows\System32"
        with mock.patch.dict(
            env_manager.os.environ,
            {
                "PATH": original_path,
                "CONDA_PREFIX": r"C:\Users\danie\.conda\envs\luna-gui",
                "CONDA_DEFAULT_ENV": "luna-gui",
                "PYTHONHOME": "bad-home",
                "PYTHONPATH": "bad-path",
            },
            clear=True,
        ):
            env = env_manager.conda_process_env(conda)

        path_parts = env["PATH"].split(";")
        self.assertEqual(path_parts[0], r"C:\ProgramData\Anaconda3")
        self.assertEqual(path_parts[3], r"C:\ProgramData\Anaconda3\Library\bin")
        self.assertEqual(path_parts[4], r"C:\ProgramData\Anaconda3\Scripts")
        self.assertEqual(env["CONDA_EXE"], conda)
        self.assertEqual(env["PYTHONNOUSERSITE"], "1")
        self.assertNotIn("CONDA_PREFIX", env)
        self.assertNotIn("CONDA_DEFAULT_ENV", env)
        self.assertNotIn("PYTHONHOME", env)
        self.assertNotIn("PYTHONPATH", env)

    def test_python_process_env_prefers_env_bins_and_disables_user_site(self) -> None:
        py = r"C:\Users\danie\.conda\envs\luna-env\python.exe"
        original_path = r"C:\Users\danie\AppData\Roaming\Python\Python39\site-packages;C:\Windows\System32"
        with mock.patch.dict(
            env_manager.os.environ,
            {
                "PATH": original_path,
                "PYTHONHOME": "bad-home",
                "PYTHONPATH": "bad-path",
                "PYTHONUSERBASE": r"C:\Users\danie\AppData\Roaming\Python",
            },
            clear=True,
        ):
            env = env_manager.python_process_env(py)

        path_parts = env["PATH"].split(";")
        self.assertEqual(path_parts[0], r"C:\Users\danie\.conda\envs\luna-env")
        self.assertEqual(path_parts[3], r"C:\Users\danie\.conda\envs\luna-env\Library\bin")
        self.assertEqual(path_parts[4], r"C:\Users\danie\.conda\envs\luna-env\Scripts")
        self.assertEqual(env["PYTHONNOUSERSITE"], "1")
        self.assertNotIn("PYTHONHOME", env)
        self.assertNotIn("PYTHONPATH", env)
        self.assertNotIn("PYTHONUSERBASE", env)

    def test_find_conda_prefers_current_gui_env_over_path(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            preferred_root = root / "user-home" / ".conda"
            preferred_exe = preferred_root / "Scripts" / "conda.exe"
            preferred_python = preferred_root / "envs" / "luna-gui" / "python.exe"
            wrong_exe = root / "ProgramData" / "Anaconda3" / "Scripts" / "conda.exe"

            preferred_exe.parent.mkdir(parents=True, exist_ok=True)
            preferred_python.parent.mkdir(parents=True, exist_ok=True)
            wrong_exe.parent.mkdir(parents=True, exist_ok=True)
            preferred_exe.write_text("", encoding="utf-8")
            preferred_python.write_text("", encoding="utf-8")
            wrong_exe.write_text("", encoding="utf-8")

            with mock.patch.dict(env_manager.os.environ, {}, clear=True):
                with mock.patch.object(env_manager.Path, "home", return_value=root / "user-home"):
                    with mock.patch.object(env_manager.sys, "executable", str(preferred_python)):
                        with mock.patch.object(env_manager.sys, "prefix", str(preferred_python.parent)):
                            with mock.patch.object(
                                env_manager.shutil,
                                "which",
                                side_effect=lambda name: str(wrong_exe) if name == "conda" else None,
                            ):
                                found = env_manager.find_conda()

            self.assertEqual(found, str(preferred_exe))

    def test_find_conda_discovers_hidden_dot_conda_install(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            home = root / "home"
            hidden_root = home / ".conda"
            hidden_exe = hidden_root / "Scripts" / "conda.exe"
            hidden_exe.parent.mkdir(parents=True, exist_ok=True)
            hidden_exe.write_text("", encoding="utf-8")

            unrelated = root / "plain-python" / "python.exe"
            unrelated.parent.mkdir(parents=True, exist_ok=True)
            unrelated.write_text("", encoding="utf-8")

            with mock.patch.dict(env_manager.os.environ, {}, clear=True):
                with mock.patch.object(env_manager.Path, "home", return_value=home):
                    with mock.patch.object(env_manager.sys, "executable", str(unrelated)):
                        with mock.patch.object(env_manager.sys, "prefix", str(unrelated.parent)):
                            with mock.patch.object(env_manager.shutil, "which", return_value=None):
                                found = env_manager.find_conda()

            self.assertEqual(found, str(hidden_exe))

    def test_cleanup_partial_env_removes_incomplete_target_directory(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            prefix = root / "envs" / "luna-env"
            (prefix / ".condatmp").mkdir(parents=True, exist_ok=True)
            (prefix / "conda-meta").mkdir(parents=True, exist_ok=True)
            conda = str(root / "Scripts" / "conda.exe")
            info = {"envs": [], "envs_dirs": [str(prefix.parent)]}

            with mock.patch.object(env_manager, "conda_info", return_value=info):
                removed = env_manager.cleanup_partial_env(conda)

            self.assertEqual(removed, prefix)
            self.assertFalse(prefix.exists())

    def test_install_commands_use_absolute_prefix_for_create_and_run(self) -> None:
        conda = r"C:\ProgramData\Anaconda3\Scripts\conda.exe"
        prefix = Path(r"C:\Users\danie\.conda\envs\luna-env")

        with mock.patch.object(env_manager, "env_prefix", return_value=prefix):
            with mock.patch.object(env_manager, "env_is_valid", return_value=False):
                commands = env_manager.install_commands(conda)

        self.assertEqual(commands[0][:4], [conda, "create", "-p", str(prefix)])
        self.assertEqual(commands[1][:4], [conda, "run", "-p", str(prefix)])
        self.assertEqual(commands[2][:4], [conda, "run", "-p", str(prefix)])
        self.assertEqual(commands[3][:4], [conda, "run", "-p", str(prefix)])
        self.assertEqual(commands[1][4:9], ["python", "-s", "-m", "pip", "install"])
        self.assertEqual(commands[2][4:9], ["python", "-s", "-m", "pip", "install"])
        self.assertEqual(commands[3][4:6], ["python", "-s"])
        self.assertIn("scikit-learn", commands[0])


if __name__ == "__main__":
    unittest.main()
