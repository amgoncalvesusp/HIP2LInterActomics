from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import hipplinteractomics_multiple_run as multiple_run
import hipplinteractomics_terminal as terminal
from luna_gui.core.mol2_prep import PrepResult
from luna_gui.core.project import ProjectConfig


class TerminalComplexPreparationTests(unittest.TestCase):
    def test_terminal_config_accepts_automatic_complex_preparation_fields(self) -> None:
        project, terminal_data = terminal._normalize_project_data({
            "project": {
                "workdir": "workdir",
                "protein_file": "old.pdb",
                "ligand_file": "old.sdf",
            },
            "terminal": {
                "prepare_complexes": True,
                "complex_folder": "complexes",
                "prepare_output": "prepared",
                "last_protein_atom": 42,
            },
        })

        self.assertEqual(project["workdir"], "workdir")
        self.assertTrue(terminal_data["prepare_complexes"])
        self.assertEqual(terminal_data["last_protein_atom"], 42)

    def test_preparation_replaces_inputs_before_ligand_resolution(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            complexes = root / "complexes"
            complexes.mkdir()
            prepared = root / "prepared"
            protein_dir = prepared / "proteinas_pdb"
            ligand_dir = prepared / "ligantes_sdf"
            result = PrepResult(
                files_processed=2,
                proteins_written=2,
                ligands_written=2,
                water_molecules_detected=0,
                protein_dir=str(protein_dir),
                ligand_dir=str(ligand_dir),
                errors=[],
            )
            cfg = ProjectConfig(
                workdir=str(root / "work"),
                protein_file="old.pdb",
                ligand_file="old.sdf",
                selected_ligands=["old"],
            )

            with patch.object(terminal.mol2_prep, "split_complex_folder", return_value=result) as split:
                terminal._prepare_configured_complexes(
                    cfg,
                    {
                        "prepare_complexes": True,
                        "complex_folder": str(complexes),
                        "prepare_output": str(prepared),
                        "last_protein_atom": "17",
                    },
                    root / "luna-python",
                )

            split.assert_called_once_with(
                complexes.resolve(),
                last_pa=17,
                out_folder=prepared.resolve(),
                chemistry_python=root / "luna-python",
            )
            self.assertEqual(cfg.protein_file, str(protein_dir))
            self.assertEqual(cfg.ligand_file, str(ligand_dir))
            self.assertEqual(cfg.selected_ligands, [])

    def test_multiple_run_keeps_preparation_per_combination(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            base = {
                "project": {
                    "workdir": str(root / "runs"),
                    "ligand_file": str(root / "will-be-created-by-preparation.sdf"),
                    "selected_ligands": [],
                },
                "terminal": {
                    "prepare_complexes": True,
                    "complex_folder": str(root / "complexes"),
                    "prepare_output": str(root / "shared-output-must-not-be-used"),
                },
            }

            prepared = multiple_run.prepare_static_inputs(base)
            generated = multiple_run.generate_configurations(
                base_config=prepared,
                bits=[1024],
                formats=["cnt"],
                levels_growth=[(2, 5.0)],
                config_dir=root / "configs",
            )

            document = json.loads(generated[0].config_path.read_text(encoding="utf-8"))
            self.assertEqual(document["project"]["selected_ligands"], [])
            self.assertEqual(document["terminal"]["prepare_output"], "")
            self.assertEqual(
                document["project"]["workdir"],
                str((root / "runs" / generated[0].run_id).resolve()),
            )


if __name__ == "__main__":
    unittest.main()
