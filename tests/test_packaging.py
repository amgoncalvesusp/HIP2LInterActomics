from __future__ import annotations

from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def _read(relative: str) -> str:
    return (ROOT / relative).read_text(encoding="utf-8")


def test_snap_uses_qt6_runtime_with_the_current_release_version() -> None:
    snapcraft = _read("snap/snapcraft.yaml")
    assert "base: core24" in snapcraft
    assert 'version: "1.0.0"' in snapcraft
    assert "python3-pyqt6" in snapcraft
    assert "qt6-qpa-plugins" in snapcraft
    assert "qt6-wayland" in snapcraft


def test_linux_installer_validates_qt_and_generates_a_quoted_launcher() -> None:
    installer = _read("dist/linux/install_hip2linteractomics.sh")
    assert "check_qt_runtime" in installer
    assert "ldd" in installer
    assert "printf 'exec %q" in installer
    assert "--conda-root requer um diretorio" in installer


def test_windows_conda_installer_writes_a_utf8_launcher() -> None:
    installer = _read("dist/windows/install_hip2linteractomics.ps1")
    assert "UTF8Encoding" in installer
    assert "-Encoding ASCII" not in installer
    assert "chcp 65001" in installer
    assert "-Wait -PassThru -WindowStyle Hidden" in installer


def test_native_windows_distribution_has_the_required_build_artifacts() -> None:
    spec = _read("HIP2LInterActomics.spec")
    inno = _read("installer/HIP2LInterActomics.iss")
    hook = _read("rthook_qt_bundled.py")

    assert '("luna_gui", "luna_gui")' in spec
    assert 'icon="luna_gui/assets/hip2l_interactomics_icon.ico"' in spec
    assert "os.add_dll_directory" in hook
    assert "PrivilegesRequired=lowest" in inno
    assert '#define MyAppVersion "1.0.0"' in inno
    assert (ROOT / "luna_gui/assets/hip2l_interactomics_icon.ico").is_file()
