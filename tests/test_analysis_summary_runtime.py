from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from luna_gui.core import analysis_runtime


class AnalysisSummaryRuntimeTests(unittest.TestCase):
    def test_rebuilds_and_caches_summary_from_residue_matrix(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            workdir = Path(temporary)
            results = workdir / "results"
            results.mkdir()
            (results / "residue_matrix.json").write_text(
                json.dumps(
                    {
                        "entries": ["lig1", "lig2"],
                        "residues": ["A/ASP/42", "A/PHE/99"],
                        "matrix": {
                            "Hydrogen bond": [[2, 0], [1, 0]],
                            "Hydrophobic": [[0, 3], [0, 4]],
                        },
                    }
                ),
                encoding="utf-8",
            )

            with patch.object(analysis_runtime.analysis_helper, "run_analysis") as legacy_helper:
                result = analysis_runtime.run_analysis("missing-python", str(workdir))

            legacy_helper.assert_not_called()
            self.assertEqual(result["entries"], 2)
            self.assertEqual(result["interaction_counts"]["Hydrogen bond"], 3)
            self.assertEqual(result["interaction_counts"]["Hydrophobic"], 7)
            self.assertEqual(result["residue_counts"]["A/PHE/99"], 7)
            cached = json.loads((results / "analysis_summary.json").read_text(encoding="utf-8"))
            self.assertEqual(cached["source"], "residue_matrix.json")
            self.assertEqual(cached["schema_version"], analysis_runtime.ANALYSIS_SUMMARY_SCHEMA_VERSION)
            self.assertEqual(len(cached["source_signature"]), 64)

    def test_invalidates_summary_when_residue_matrix_changes(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            workdir = Path(temporary)
            results = workdir / "results"
            results.mkdir()
            matrix_path = results / "residue_matrix.json"
            matrix_path.write_text(
                json.dumps({
                    "entries": ["lig1"],
                    "residues": ["A/ASP/42"],
                    "matrix": {"Hydrogen bond": [[1]]},
                }),
                encoding="utf-8",
            )
            first = analysis_runtime.run_analysis("missing-python", str(workdir))
            matrix_path.write_text(
                json.dumps({
                    "entries": ["lig1"],
                    "residues": ["A/ASP/42"],
                    "matrix": {"Hydrogen bond": [[9]]},
                }),
                encoding="utf-8",
            )

            second = analysis_runtime.run_analysis("missing-python", str(workdir))

            self.assertEqual(first["residue_counts"]["A/ASP/42"], 1)
            self.assertEqual(second["residue_counts"]["A/ASP/42"], 9)


if __name__ == "__main__":
    unittest.main()
