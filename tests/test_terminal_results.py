from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path

from luna_gui.core.project import ProjectConfig
from luna_gui.core.terminal_results import run_terminal_results


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
            self.assertTrue((output_dir / "complete_ligands_residues_heatmap.png").exists())
            self.assertIn("complete_ligands_residues_heatmap", manifest["outputs"])


if __name__ == "__main__":
    unittest.main()
