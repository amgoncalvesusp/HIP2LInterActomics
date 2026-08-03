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
        self.assertIn("_create_scroll_page", methods)
        self.assertIn("_set_current_tab", methods)
        helper_source = ast.get_source_segment(source, methods["_create_scroll_page"]) or ""
        self.assertIn("QScrollArea", helper_source)
        self.assertIn("ScrollBarAlwaysOff", helper_source)
        self.assertIn("QSizePolicy.Policy.Ignored", helper_source)

        init_source = ast.get_source_segment(source, methods["__init__"]) or ""
        self.assertEqual(init_source.count("self._add_scrollable_tab("), 6)

    def test_results_module_is_loaded_only_when_the_tab_is_opened(self) -> None:
        source = (Path(__file__).resolve().parents[1] / "luna_gui" / "ui" / "main_window.py").read_text(
            encoding="utf-8"
        )
        tree = ast.parse(source)
        runtime_imports = [
            node
            for node in tree.body
            if isinstance(node, ast.ImportFrom) and node.module == "tab_results_advanced"
        ]
        self.assertEqual(runtime_imports, [])
        self.assertIn("def _ensure_results_tab", source)
        self.assertIn("QTimer.singleShot(0, self._ensure_results_tab)", source)

    def test_setup_probe_starts_only_after_runtime_signals_are_connected(self) -> None:
        source = (Path(__file__).resolve().parents[1] / "luna_gui" / "ui" / "main_window.py").read_text(
            encoding="utf-8"
        )
        run_connection = source.index("self.tab_setup.luna_ready.connect(self.tab_run.set_luna)")
        probe_start = source.index("QTimer.singleShot(0, self.tab_setup.detect)")
        self.assertLess(run_connection, probe_start)

    def test_fp_plot_canvases_are_created_lazily(self) -> None:
        source = (
            Path(__file__).resolve().parents[1]
            / "luna_gui"
            / "ui"
            / "tab_results_advanced.py"
        ).read_text(encoding="utf-8")
        tree = ast.parse(source)
        results_tab = next(
            node for node in tree.body if isinstance(node, ast.ClassDef) and node.name == "ResultsTab"
        )
        methods = {
            node.name: node
            for node in results_tab.body
            if isinstance(node, ast.FunctionDef)
        }
        installer_source = ast.get_source_segment(source, methods["_install_fp_analysis_tab"]) or ""
        ensure_source = ast.get_source_segment(source, methods["_ensure_fp_plot_canvases"]) or ""

        self.assertNotIn("Figure(figsize=", installer_source)
        self.assertIn("Figure(figsize=size)", ensure_source)
        self.assertIn("self._ensure_fp_plot_canvases()", source)

    def test_trajectory_statistics_use_vertical_chart_legend_sequence(self) -> None:
        source = (
            Path(__file__).resolve().parents[1]
            / "luna_gui"
            / "ui"
            / "tab_results_advanced.py"
        ).read_text(encoding="utf-8")
        self.assertIn("height_ratios.extend([chart_height_in, residue_legend_height])", source)
        self.assertIn("height_ratios.extend([chart_height_in, atom_legend_height])", source)
        self.assertIn("chart_height_in = min(width_in * 2.0", source)
        self.assertIn("legend_gap_in = 2.0 / 2.54", source)
        self.assertIn("self.stats_scroll.setMinimumHeight(560)", source)


if __name__ == "__main__":
    unittest.main()
