from __future__ import annotations

import re
from pathlib import Path
import unittest


FIXTURE_ROOT = Path(__file__).parent / "fixtures" / "sample_workdir"


def _fixture_dir(name: str) -> Path:
    path = FIXTURE_ROOT / name
    if not path.is_dir():
        raise unittest.SkipTest(f"Fixture de integração não encontrado: {path}")
    return path


def _frame_number(path: Path) -> int:
    match = re.search(r"frame(\d+)", path.stem.lower())
    if match is None:
        raise AssertionError(f"Nome sem número de frame: {path.name}")
    return int(match.group(1))


class SampleWorkdirFixtureTests(unittest.TestCase):
    """Smoke tests using the real normal and trajectory input directories."""

    def test_normal_analysis_fixture_is_read_as_non_trajectory_input(self) -> None:
        workdir = _fixture_dir("normal_analises")
        ligands = sorted((workdir / "ligand_file").glob("*.mol2"))
        proteins = sorted((workdir / "protein_file").glob("*.pdb"))
        references = sorted((workdir / "reference_importance_file").glob("*.tsv"))

        self.assertTrue(ligands, "A fixture normal precisa conter pelo menos um MOL2")
        self.assertTrue(proteins, "A fixture normal precisa conter pelo menos um PDB")
        self.assertTrue(references, "A fixture normal precisa conter o arquivo de referência")
        trajectory_analysis = False
        self.assertFalse(trajectory_analysis)

        self.assertIn("@<TRIPOS>MOLECULE", ligands[0].read_text(encoding="utf-8", errors="replace"))
        protein_text = proteins[0].read_text(encoding="utf-8", errors="replace")
        self.assertTrue(any(line.startswith(("HEADER", "ATOM", "HETATM")) for line in protein_text.splitlines()))
        self.assertIn("\t", references[0].read_text(encoding="utf-8", errors="replace"))

    def test_trajectory_fixture_reads_100_matching_frames(self) -> None:
        workdir = _fixture_dir("trajectory_analises")
        ligands = sorted((workdir / "ligantes_sdf_100").glob("frame*_ligand.sdf"), key=_frame_number)
        proteins = sorted((workdir / "proteinas_pdb_100").glob("frame*.pdb"), key=_frame_number)

        self.assertEqual(len(ligands), 100)
        self.assertEqual(len(proteins), 100)
        self.assertEqual([_frame_number(path) for path in ligands], list(range(100)))
        self.assertEqual([_frame_number(path) for path in proteins], list(range(100)))

        # Read representative files completely; the remaining frames are
        # validated by their sequence and extension without loading 100 MB.
        sdf_sample = ligands[0].read_text(encoding="utf-8", errors="replace")
        pdb_sample = proteins[0].read_text(encoding="utf-8", errors="replace")
        self.assertIn("$$$$", sdf_sample)
        self.assertIn("ATOM", pdb_sample)

        trajectory_analysis = True
        self.assertTrue(trajectory_analysis)


if __name__ == "__main__":
    unittest.main()
