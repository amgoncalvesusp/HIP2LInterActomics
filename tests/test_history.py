from __future__ import annotations

import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from luna_gui.core import project


class HistoryTests(unittest.TestCase):
    def test_remove_and_clear_history_preserve_project_directories(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            first = root / "first"
            second = root / "second"
            first.mkdir()
            second.mkdir()
            history_file = root / "history.json"
            with patch.object(project, "HISTORY_FILE", history_file):
                project.add_to_history(str(first))
                project.add_to_history(str(second))
                project.remove_from_history(str(first))
                self.assertEqual(project.load_history(), [str(second)])
                project.clear_history()
                self.assertEqual(project.load_history(), [])
            self.assertTrue(first.exists())
            self.assertTrue(second.exists())


if __name__ == "__main__":
    unittest.main()
