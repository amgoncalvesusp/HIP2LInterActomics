from __future__ import annotations

import ast
from pathlib import Path
import unittest


class MainWindowScrollTests(unittest.TestCase):
    def test_main_tabs_are_wrapped_in_scroll_areas(self) -> None:
        source = (Path(__file__).resolve().parents[1] / "luna_gui" / "ui" / "main_window.py").read_text(
            encoding="utf-8"
        )
        tree = ast.parse(source)
        main_window = next(
            node for node in tree.body if isinstance(node, ast.ClassDef) and node.name == "MainWindow"
        )
        methods = {
            node.name: node
            for node in main_window.body
            if isinstance(node, ast.FunctionDef)
        }

        self.assertIn("_add_scrollable_tab", methods)
        self.assertIn("_set_current_tab", methods)
        helper_source = ast.get_source_segment(source, methods["_add_scrollable_tab"]) or ""
        self.assertIn("QScrollArea", helper_source)
        self.assertIn("ScrollBarAsNeeded", helper_source)

        init_source = ast.get_source_segment(source, methods["__init__"]) or ""
        self.assertEqual(init_source.count("self._add_scrollable_tab("), 6)


if __name__ == "__main__":
    unittest.main()
