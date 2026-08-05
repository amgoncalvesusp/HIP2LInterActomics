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
        with (
            mock.patch.object(env_manager, "conda_info", return_value=info),
            mock.patch.object(env_manager.Path, "home", return_value=Path(r"C:\\Users\\danie")),
        ):
            prefix = env_manager.env_prefix(conda)

        self.assertEqual(prefix, Path(r"C:\Users\danie\.conda\envs\luna-env"))

    def test_env_prefix_finds_user_env_when_frozen_conda_info_fails(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            home = root / "user-home"
            prefix = home / ".conda" / "envs" / "luna-env"
            (prefix / "conda-meta").mkdir(parents=True)
            (prefix / "conda-meta" / "history").write_text("created", encoding="utf-8")
            (prefix / "python.exe").write_text("", encoding="utf-8")
            conda = str(root / "ProgramData" / "Anaconda3" / "Scripts" / "conda.exe")

            with mock.patch.object(env_manager, "conda_info", return_value={}):
                with mock.patch.object(env_manager.Path, "home", return_value=home):
                    with mock.patch.object(env_manager.sys, "platform", "win32"):
                        found = env_manager.env_prefix(conda)

        self.assertEqual(found, prefix)

    def test_new_env_uses_user_prefix_when_conda_is_installed_system_wide(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            home = root / "regular-user"
            system_root = root / "ProgramData" / "miniconda3"
            (system_root / "envs").mkdir(parents=True)
            conda = str(system_root / "Scripts" / "conda.exe")
            info = {
                "envs": [],
                "envs_dirs": [
                    str(system_root / "envs"),
                    str(home / ".conda" / "envs"),
                ],
            }

            with mock.patch.object(env_manager, "conda_info", return_value=info):
                with mock.patch.object(env_manager.Path, "home", return_value=home):
                    prefix = env_manager.env_prefix(conda)

        self.assertEqual(prefix, home / ".conda" / "envs" / "luna-env")

    def test_new_env_honors_explicit_luna_prefix_override(self) -> None:
        conda = r"C:\ProgramData\miniconda3\Scripts\conda.exe"
        override = r"D:\chemistry-envs\luna-env"
        with mock.patch.dict(
            env_manager.os.environ,
            {"HIP2LINTERACTOMICS_LUNA_ENV": override},
            clear=True,
        ):
            with mock.patch.object(env_manager, "conda_info", return_value={}):
                prefix = env_manager.env_prefix(conda)

        self.assertEqual(prefix, Path(override))

    def test_luna_runtime_is_located_without_starting_external_python(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            prefix = Path(tmp) / "luna-env"
            py = prefix / "python.exe"
            run_py = prefix / "Lib" / "site-packages" / "luna" / "run.py"
            py.parent.mkdir(parents=True, exist_ok=True)
            run_py.parent.mkdir(parents=True, exist_ok=True)
            py.write_text("", encoding="utf-8")
            run_py.write_text("", encoding="utf-8")

            with mock.patch.object(env_manager.sys, "platform", "win32"):
                with mock.patch.object(
                    env_manager.subprocess,
                    "run",
                    side_effect=AssertionError("external Python should not be started"),
                ):
                    self.assertTrue(env_manager.luna_installed(py))
                    self.assertEqual(env_manager.luna_run_py_path(py), run_py)

    def test_find_luna_runtime_uses_user_conda_env_without_conda_process(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            home = Path(tmp) / "home"
            prefix = home / ".conda" / "envs" / "luna-env"
            py = prefix / "python.exe"
            run_py = prefix / "Lib" / "site-packages" / "luna" / "run.py"
            py.parent.mkdir(parents=True, exist_ok=True)
            run_py.parent.mkdir(parents=True, exist_ok=True)
            py.write_text("", encoding="utf-8")
            run_py.write_text("", encoding="utf-8")

            with mock.patch.object(env_manager.Path, "home", return_value=home):
                with mock.patch.object(env_manager.sys, "platform", "win32"):
                    with mock.patch.object(
                        env_manager.subprocess,
                        "run",
                        side_effect=AssertionError("Conda subprocess should not be needed"),
                    ):
                        runtime = env_manager.find_luna_runtime()

        self.assertEqual(runtime, (py, run_py))

    def test_find_conda_honors_manual_override(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            override_exe = root / "custom" / "Scripts" / "conda.exe"
            override_exe.parent.mkdir(parents=True, exist_ok=True)
            override_exe.write_text("", encoding="utf-8")

            with mock.patch.dict(
                env_manager.os.environ,
                {
                    "HIP2LINTERACTOMICS_GUI_CONDA": str(override_exe),
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
        self.assertEqual(env["CONDA_PREFIX"], r"C:\Users\danie\.conda\envs\luna-env")
        self.assertEqual(env["CONDA_DEFAULT_ENV"], "luna-env")
        self.assertEqual(env["PYTHONNOUSERSITE"], "1")
        self.assertNotIn("PYTHONHOME", env)
        self.assertNotIn("PYTHONPATH", env)
        self.assertNotIn("PYTHONUSERBASE", env)

    def test_chemistry_process_env_limits_native_thread_pools(self) -> None:
        with mock.patch.dict(
            env_manager.os.environ,
            {"PATH": r"C:\Windows\System32", "OPENBLAS_NUM_THREADS": "32"},
            clear=True,
        ):
            env = env_manager.chemistry_process_env(
                r"C:\Users\danie\.conda\envs\luna-env\python.exe"
            )

        for key in (
            "OMP_NUM_THREADS",
            "OPENBLAS_NUM_THREADS",
            "MKL_NUM_THREADS",
            "NUMEXPR_NUM_THREADS",
            "VECLIB_MAXIMUM_THREADS",
        ):
            self.assertEqual(env[key], "1")

    def test_python_process_env_restores_pre_bundle_linux_library_path(self) -> None:
        with mock.patch.dict(
            env_manager.os.environ,
            {
                "PATH": "/app/bin:/usr/bin",
                "APPIMAGE": "/downloads/HIP2LInterActomics.AppImage",
                "LD_LIBRARY_PATH": "/tmp/_MEI/lib",
                "LD_LIBRARY_PATH_ORIG": "/usr/local/lib",
                "QT_PLUGIN_PATH": "/tmp/_MEI/PyQt6/plugins",
            },
            clear=True,
        ):
            env = env_manager.python_process_env("/opt/conda/envs/luna-env/bin/python")

        self.assertEqual(env["LD_LIBRARY_PATH"], "/usr/local/lib")
        self.assertNotIn("LD_LIBRARY_PATH_ORIG", env)
        self.assertNotIn("QT_PLUGIN_PATH", env)

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

    def test_find_conda_accepts_linux_conda_exe_without_extension(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            home = root / "home" / "laqmedsom"
            conda_root = home / "softwares_laqmedsomm" / "anaconda3"
            conda_exe = conda_root / "bin" / "conda"
            gui_python = home / ".conda" / "envs" / "luna-gui" / "bin" / "python"
            conda_exe.parent.mkdir(parents=True, exist_ok=True)
            gui_python.parent.mkdir(parents=True, exist_ok=True)
            conda_exe.write_text("", encoding="utf-8")
            gui_python.write_text("", encoding="utf-8")

            with mock.patch.dict(env_manager.os.environ, {"CONDA_EXE": str(conda_exe)}, clear=True):
                with mock.patch.object(env_manager.sys, "platform", "linux"):
                    with mock.patch.object(env_manager.Path, "home", return_value=home):
                        with mock.patch.object(env_manager.sys, "executable", str(gui_python)):
                            with mock.patch.object(env_manager.sys, "prefix", str(gui_python.parent.parent)):
                                with mock.patch.object(env_manager.shutil, "which", return_value=None):
                                    found = env_manager.find_conda()

            self.assertEqual(found, str(conda_exe))

    def test_find_conda_discovers_nested_linux_home_install(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            home = root / "home" / "laqmedsom"
            conda_root = home / "softwares_laqmedsomm" / "anaconda3"
            conda_exe = conda_root / "bin" / "conda"
            unrelated = home / "plain-python" / "bin" / "python"
            conda_exe.parent.mkdir(parents=True, exist_ok=True)
            unrelated.parent.mkdir(parents=True, exist_ok=True)
            conda_exe.write_text("", encoding="utf-8")
            unrelated.write_text("", encoding="utf-8")

            with mock.patch.dict(env_manager.os.environ, {}, clear=True):
                with mock.patch.object(env_manager.sys, "platform", "linux"):
                    with mock.patch.object(env_manager.Path, "home", return_value=home):
                        with mock.patch.object(env_manager.sys, "executable", str(unrelated)):
                            with mock.patch.object(env_manager.sys, "prefix", str(unrelated.parent.parent)):
                                with mock.patch.object(env_manager.shutil, "which", return_value=None):
                                    found = env_manager.find_conda()

            self.assertEqual(found, str(conda_exe))

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


    def test_snap_qt_paths_are_removed_from_external_conda_environment(self) -> None:
        with mock.patch.dict(
            env_manager.os.environ,
            {
                "PATH": r"C:\Windows\System32",
                "HIP2LINTERACTOMICS_SNAP": "1",
                "LD_LIBRARY_PATH": "/snap/lib",
                "QML2_IMPORT_PATH": "/snap/qml",
                "QT_PLUGIN_PATH": "/snap/plugins",
                "QT_QPA_PLATFORM_PLUGIN_PATH": "/snap/platforms",
            },
            clear=True,
        ):
            env = env_manager.python_process_env(r"C:\envs\luna-env\python.exe")

        self.assertEqual(env["HIP2LINTERACTOMICS_SNAP"], "1")
        self.assertNotIn("LD_LIBRARY_PATH", env)
        self.assertNotIn("QML2_IMPORT_PATH", env)
        self.assertNotIn("QT_PLUGIN_PATH", env)
        self.assertNotIn("QT_QPA_PLATFORM_PLUGIN_PATH", env)

    def test_miniconda_download_url_uses_linux_arm64_build(self) -> None:
        with (
            mock.patch.object(env_manager.sys, "platform", "linux"),
            mock.patch.object(env_manager.platform, "machine", return_value="aarch64"),
        ):
            url = env_manager.miniconda_download_url()

        self.assertTrue(url.endswith("Linux-aarch64.sh"))

    def test_miniconda_download_url_rejects_unknown_linux_architecture(self) -> None:
        with (
            mock.patch.object(env_manager.sys, "platform", "linux"),
            mock.patch.object(env_manager.platform, "machine", return_value="riscv64"),
        ):
            with self.assertRaisesRegex(RuntimeError, "riscv64"):
                env_manager.miniconda_download_url()


if __name__ == "__main__":
    unittest.main()
