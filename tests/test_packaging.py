from __future__ import annotations

from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def _read(relative: str) -> str:
    return (ROOT / relative).read_text(encoding="utf-8")


def test_snap_uses_available_pyqt6_base_and_wayland_plugin() -> None:
    snapcraft = _read("snap/snapcraft.yaml")
    assert "base: core24" in snapcraft
    assert 'version: "1.0.0"' in snapcraft
    assert "python3-pyqt6" in snapcraft
    assert "qt6-qpa-plugins" in snapcraft
    assert "qt6-wayland" in snapcraft


def test_linux_installer_validates_qt_and_generates_quoted_wrapper() -> None:
    installer = _read("dist/linux/install_hip2linteractomics.sh")
    assert "check_qt_runtime" in installer
    assert "ldd" in installer
    assert "printf 'exec %q" in installer
    assert '--conda-root requer um diretorio' in installer


def test_windows_launcher_is_written_as_utf8_without_bom() -> None:
    installer = _read("dist/windows/install_hip2linteractomics.ps1")
    assert "UTF8Encoding" in installer
    assert "-Encoding ASCII" not in installer
    assert "chcp 65001" in installer
    assert "Start-Process -FilePath $installer" in installer
    assert "-Wait -PassThru" in installer


def test_installer_versions_are_consistent() -> None:
    snapcraft = _read("snap/snapcraft.yaml")
    inno = _read("installer/HIP2LInterActomics.iss")
    assert 'version: "1.0.0"' in snapcraft
    assert '#define MyAppVersion "1.0.0"' in inno
