from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from luna_gui.core.project import ProjectConfig
from luna_gui.core.plot_manifest import load_manifest
from luna_gui.core.terminal_results import _plot_worker_count, run_terminal_results


class TerminalResultsTests(unittest.TestCase):
    def test_exports_cached_results_and_interactive_cluster_explorer(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            workdir = Path(tmp)
            results = workdir / "results"
            results.mkdir()
            (results / "analysis_summary.json").write_text(
                json.dumps({"entries": 3, "interaction_counts": {"Hydrogen bond": 4, "Hydrophobic": 2}}),
                encoding="utf-8",
            )
            (results / "residue_matrix.json").write_text(
                json.dumps(
                    {
                        "entries": ["ligA", "ligB", "ligC"],
                        "residues": ["A/GLY/1", "A/TYR/2"],
                        "interaction_types": ["Hydrogen bond"],
                        "matrix": {"Hydrogen bond": [[1, 0], [0, 2], [3, 1]]},
                    }
                ),
                encoding="utf-8",
            )
            (workdir / "sim_matrix_E_square.csv").write_text(
                ",ligA,ligB,ligC\n"
                "ligA,1,0.9,0.2\n"
                "ligB,0.9,1,0.3\n"
                "ligC,0.2,0.3,1\n",
                encoding="utf-8",
            )
            cfg = ProjectConfig(workdir=str(workdir), ifp_type="EIFP", sim_matrix=True)

            manifest = run_terminal_results(
                cfg,
                sys.executable,
                {
                    "terminal_matrix_max_entries": 10,
                    "terminal_cluster_max_entries": 10,
                    "terminal_interactive_max_entries": 10,
                },
            )

            output_dir = results / "terminal"
            self.assertTrue((output_dir / "terminal_results_manifest.json").exists())
            self.assertTrue((output_dir / "clusters_E.csv").exists())
            explorer = output_dir / "clusters_E.html"
            self.assertTrue(explorer.exists())
            self.assertIn("matrix_available", explorer.read_text(encoding="utf-8"))
            self.assertIn("EIFP_cluster_explorer", manifest["outputs"])
            plot_manifest = load_manifest(workdir)
            self.assertEqual(manifest["outputs"]["plot_count"], len(plot_manifest.records))
            for language in ("en", "pt", "es"):
                for profile in ("screen", "report"):
                    complete = (
                        results
                        / "plots"
                        / language
                        / profile
                        / "heatmaps"
                        / "complete_ligands_residues_heatmap.png"
                    )
                    self.assertTrue(complete.exists())
                    selected = plot_manifest.select(language=language, profile=profile)
                    self.assertTrue(any(row.plot_id == "complete_interaction_heatmap" for row in selected))

    def test_extreme_memory_pressure_forces_sequential_language_workers(self) -> None:
        with patch(
            "luna_gui.core.terminal_results._available_memory_bytes",
            return_value=2 * 1024 ** 3,
        ):
            self.assertEqual(_plot_worker_count(), 1)


if __name__ == "__main__":
    unittest.main()
