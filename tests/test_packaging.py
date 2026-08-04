from __future__ import annotations

import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def _read(relative: str) -> str:
    return (ROOT / relative).read_text(encoding="utf-8")


def test_native_windows_distribution_has_required_build_artifacts() -> None:
    spec = _read("HIP2LInterActomics.spec")
    inno = _read("installer/HIP2LInterActomics.iss")
    hook = _read("rthook_qt_bundled.py")

    assert '("luna_gui", "luna_gui")' in spec
    assert 'icon="luna_gui/assets/hip2l_interactomics_icon.ico"' in spec
    assert 'contents_directory="."' in spec
    assert "os.add_dll_directory" in hook
    assert "PrivilegesRequired=lowest" in inno
    assert "ChangesEnvironment=yes" in inno
    assert "preservestringtype" not in inno
    assert "uninsneveruninstall" not in inno
    assert "autodesktop" in inno
    assert "Tasks: desktopicon" not in inno
    assert "NeedsAddPath" in inno
    assert "{app}\\bin" in inno
    assert (ROOT / "installer/hipplinteractomics-terminal.cmd").is_file()
    assert (ROOT / "installer/hipplinteractomics-multiple-run.cmd").is_file()
    assert (ROOT / "luna_gui/assets/hip2l_interactomics_icon.ico").is_file()


def test_release_version_is_consistent_across_all_packagers() -> None:
    pyproject = _read("pyproject.toml")
    inno = _read("installer/HIP2LInterActomics.iss")
    linux = _read("installer/build_linux.sh")

    project_version = re.search(r'^version = "([^"]+)"$', pyproject, re.MULTILINE)
    windows_version = re.search(r'^#define MyAppVersion "([^"]+)"$', inno, re.MULTILINE)
    linux_version = re.search(r'^APP_VERSION="\$\{APP_VERSION:-([^}]+)\}"$', linux, re.MULTILINE)

    assert project_version is not None
    assert windows_version is not None
    assert linux_version is not None
    assert project_version.group(1) == windows_version.group(1) == linux_version.group(1)


def test_native_build_scripts_cover_windows_and_linux() -> None:
    windows = _read("installer/build_windows.ps1")
    linux = _read("installer/build_linux.sh")
    workflow = _read(".github/workflows/build-installers.yml")
    app_run = _read("installer/AppRun")
    shortcut_installer = _read("installer/install_linux_shortcuts.sh")
    desktop = _read("installer/hip2linteractomics.desktop")

    assert "pip install -r requirements.txt" in windows
    assert "import PyQt6, jinja2, matplotlib, numpy, reportlab" in windows
    assert "PyInstaller" in windows
    assert "HIP2LInterActomics.exe" in windows
    assert "pip install -r" in linux
    assert "import PyQt6, jinja2, matplotlib, numpy, reportlab" in linux
    assert "PyInstaller" in linux
    assert "linux-" in linux and ".tar.gz" in linux
    assert "appimagetool" in linux
    assert ".AppImage" in linux
    assert "HIP2LInterActomics/HIP2LInterActomics" in app_run
    assert "install-linux-shortcuts" in app_run
    assert "install_linux_shortcuts.sh" in linux
    assert "XDG_DATA_HOME" in shortcut_installer
    assert "xdg-user-dir DESKTOP" in shortcut_installer
    assert "HIP2LInterActomics.desktop" in shortcut_installer
    assert '.local/bin' in shortcut_installer
    assert 'hipplinteractomics-terminal' in shortcut_installer
    assert 'hipplinteractomics-multiple-run' in shortcut_installer
    assert '--terminal' in shortcut_installer
    assert '--multiple-run' in shortcut_installer
    assert "Type=Application" in desktop
    assert "Icon=hip2linteractomics" in desktop
    assert "runs-on: windows-latest" in workflow
    assert "runs-on: ubuntu-22.04" in workflow
    assert "gh release create" in workflow
    assert '--repo "$GITHUB_REPOSITORY"' in workflow
    assert "HIP2LInterActomics-Setup.exe" in workflow
    assert "HIP2LInterActomics-*.AppImage" in workflow


def test_pyproject_registers_gui_and_headless_commands() -> None:
    metadata = _read("pyproject.toml")
    assert '[project.scripts]' in metadata
    assert 'hipplinteractomics-terminal = "hipplinteractomics_terminal:main"' in metadata
    assert 'hipplinteractomics-multiple-run = "hipplinteractomics_multiple_run:main"' in metadata
    assert 'hip2linteractomics = "luna_gui.main:main"' in metadata


def test_package_manifest_includes_cli_modules_and_assets() -> None:
    metadata = _read("pyproject.toml")
    assert '"hipplinteractomics_terminal"' in metadata
    assert '"hipplinteractomics_multiple_run"' in metadata
    assert '"assets/*"' in metadata
    assert '"examples/*"' in metadata
    assert '"share/hip2linteractomics" = ["environment.yml"]' in metadata


def test_terminal_only_environment_keeps_pymol_without_app_gui_dependencies() -> None:
    environment = _read("environment.yml")
    lowered = environment.lower()
    assert "python=3.9" in environment
    assert "matplotlib-base" in environment
    assert "seaborn-base" in environment
    assert "biopython=1.79" in environment
    assert "pymol-open-source" in environment
    assert "mmh3<4" in environment
    assert "- pdbecif" in environment
    assert "pip install --no-build-isolation -U luna" in _read("README.md")
    assert "pip install --no-deps -e ." in _read("README.md")
    for forbidden in ("pyqt6", "pyside", "tkinter"):
        assert forbidden not in lowered


def test_all_distribution_paths_include_environment_yml() -> None:
    manifest = _read("MANIFEST.in")
    spec = _read("HIP2LInterActomics.spec")
    inno = _read("installer/HIP2LInterActomics.iss")
    windows = _read("installer/build_windows.ps1")
    linux = _read("installer/build_linux.sh")

    assert "include environment.yml" in manifest
    assert '("environment.yml", ".")' in spec
    assert "recursesubdirs" in inno
    assert "environment.yml" in windows
    assert "environment.yml" in linux


def test_desktop_bundles_include_headless_cli_sources() -> None:
    spec = _read("HIP2LInterActomics.spec")
    windows = _read("installer/build_windows.ps1")
    linux = _read("installer/build_linux.sh")

    for cli in (
        "hipplinteractomics_terminal.py",
        "hipplinteractomics_multiple_run.py",
    ):
        assert f'("{cli}", ".")' in spec
        assert cli in windows
        assert cli in linux


def test_packaged_launcher_dispatches_both_cli_commands_before_loading_qt() -> None:
    launcher = _read("run.py")
    gui_import = launcher.index("from luna_gui.main import main")
    assert launcher.index('command == "--terminal"') < gui_import
    assert launcher.index('command == "--multiple-run"') < gui_import
    assert launcher.index('command == "--render-pdf-job"') < gui_import
    assert "report_worker" in launcher
    assert "AttachConsole" in launcher
