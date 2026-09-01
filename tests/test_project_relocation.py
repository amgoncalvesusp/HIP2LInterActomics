from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

from luna_gui.core.project import (
    PROJECT_FILENAME,
    ProjectConfig,
    is_project_workdir,
    relocate_config_paths,
    relocate_params_file,
    relocation_candidate,
)


class ProjectRelocationTests(unittest.TestCase):
    def _make_project_tree(self, root: Path) -> tuple[Path, Path]:
        bundle = root / "bundle"
        workdir = bundle / "work_dir"
        (bundle / "protein_file").mkdir(parents=True)
        (bundle / "ligand_file").mkdir()
        workdir.mkdir()
        (workdir / "results").mkdir()
        (bundle / "protein_file" / "receptor.pdb").write_text("ATOM\n", encoding="utf-8")
        (bundle / "ligand_file" / "library.mol2").write_text("@<TRIPOS>MOLECULE\n", encoding="utf-8")
        return bundle, workdir

    def test_relocation_rebases_workdir_and_sibling_inputs(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            old_bundle, old_workdir = self._make_project_tree(root / "old")
            new_bundle, new_workdir = self._make_project_tree(root / "new")
            old_output = old_workdir / "results" / "ifp.csv"
            new_output = new_workdir / "results" / "ifp.csv"
            old_output.write_text("old", encoding="utf-8")
            new_output.write_text("new", encoding="utf-8")

            cfg = ProjectConfig(
                workdir=str(old_workdir),
                protein_file=str(old_bundle / "protein_file" / "receptor.pdb"),
                ligand_file=str(old_bundle / "ligand_file" / "library.mol2"),
                ifp_output=str(old_output),
                fp_labels_csv=str(root / "external" / "labels.csv"),
            )
            config_path = new_workdir / PROJECT_FILENAME
            cfg.save(config_path)
            (root / "external").mkdir()
            (root / "external" / "labels.csv").write_text("name,label\n", encoding="utf-8")

            loaded = ProjectConfig.load(config_path)
            candidate = relocation_candidate(config_path, loaded)
            self.assertEqual(candidate, (str(old_workdir), new_workdir.resolve()))

            report = relocate_config_paths(loaded, str(old_workdir), new_workdir)
            self.assertEqual(Path(loaded.workdir), new_workdir.resolve())
            self.assertEqual(Path(loaded.protein_file), (new_bundle / "protein_file" / "receptor.pdb").resolve())
            self.assertEqual(Path(loaded.ligand_file), (new_bundle / "ligand_file" / "library.mol2").resolve())
            self.assertEqual(Path(loaded.ifp_output), new_output.resolve())
            self.assertEqual(loaded.fp_labels_csv, str(root / "external" / "labels.csv"))
            self.assertIn("workdir", report["changed_paths"])

    def test_params_relocation_updates_nested_entry_specs_without_touching_unresolved(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            old_bundle, old_workdir = self._make_project_tree(root / "old")
            _new_bundle, new_workdir = self._make_project_tree(root / "new")
            params_path = new_workdir / "_luna_api_params.json"
            data = {
                "workdir": str(old_workdir),
                "pdb_dir": str(old_bundle / "protein_file"),
                "entry_specs": [
                    {
                        "ligand_name": "ligA",
                        "mol_file": str(old_bundle / "ligand_file" / "library.mol2"),
                    }
                ],
                "ifp_outputs": {"EIFP": str(old_workdir / "results" / "ifp.csv")},
            }
            (new_workdir / "results" / "ifp.csv").write_text("new", encoding="utf-8")
            params_path.write_text(json.dumps(data), encoding="utf-8")

            report = relocate_params_file(params_path, str(old_workdir), new_workdir)
            updated = json.loads(params_path.read_text(encoding="utf-8"))
            self.assertEqual(Path(updated["workdir"]), new_workdir.resolve())
            self.assertEqual(Path(updated["pdb_dir"]), (new_workdir.parent / "protein_file").resolve())
            self.assertEqual(
                Path(updated["entry_specs"][0]["mol_file"]),
                (new_workdir.parent / "ligand_file" / "library.mol2").resolve(),
            )
            self.assertEqual(Path(updated["ifp_outputs"]["EIFP"]), (new_workdir / "results" / "ifp.csv").resolve())
            self.assertFalse(report["unresolved_paths"])

    def test_relocation_candidate_requires_project_markers(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            workdir = Path(tmp) / "work_dir"
            workdir.mkdir()
            config_path = workdir / PROJECT_FILENAME
            ProjectConfig(workdir=str(Path(tmp) / "old_workdir")).save(config_path)
            self.assertFalse(is_project_workdir(workdir))
            self.assertIsNone(relocation_candidate(config_path, ProjectConfig.load(config_path)))


if __name__ == "__main__":
    unittest.main()
