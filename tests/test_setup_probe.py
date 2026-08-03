from __future__ import annotations

from pathlib import Path
import unittest
from unittest.mock import patch

try:
    import PyQt6  # noqa: F401
except ImportError:
    HAS_PYQT6 = False
else:
    HAS_PYQT6 = True
    from luna_gui.ui.tab_setup import _probe_environment


@unittest.skipUnless(HAS_PYQT6, "PyQt6 is required for setup probe tests")
class SetupEnvironmentProbeTests(unittest.TestCase):
    @patch("luna_gui.ui.tab_setup.em.find_conda", return_value=None)
    def test_reports_missing_conda(self, _find_conda) -> None:
        result = _probe_environment()

        self.assertIn("Conda não encontrado", result["status"])
        self.assertFalse(result["install_enabled"])
        self.assertIsNone(result["ready"])

    @patch("luna_gui.ui.tab_setup.em.missing_runtime_packages", return_value=[])
    @patch("luna_gui.ui.tab_setup.em.find_luna_runtime", return_value=None)
    @patch(
        "luna_gui.ui.tab_setup.em.luna_run_py_path",
        return_value=Path("C:/runtime/luna/run.py"),
    )
    @patch("luna_gui.ui.tab_setup.em.luna_installed", return_value=True)
    @patch(
        "luna_gui.ui.tab_setup.em.env_python",
        return_value=Path("C:/runtime/python.exe"),
    )
    @patch("luna_gui.ui.tab_setup.em.env_exists", return_value=True)
    @patch("luna_gui.ui.tab_setup.em.env_is_partial", return_value=False)
    @patch("luna_gui.ui.tab_setup.em.env_prefix", return_value=Path("C:/runtime"))
    @patch("luna_gui.ui.tab_setup.em.find_conda", return_value="conda")
    def test_reports_ready_runtime(self, *_mocks) -> None:
        result = _probe_environment()

        self.assertEqual(
            result["ready"],
            (str(Path("C:/runtime/python.exe")), str(Path("C:/runtime/luna/run.py"))),
        )
        self.assertIn("LUNA pronto", result["status"])
        self.assertTrue(result["install_enabled"])


if __name__ == "__main__":
    unittest.main()
