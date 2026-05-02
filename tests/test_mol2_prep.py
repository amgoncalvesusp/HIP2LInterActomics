from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from luna_gui.core.mol2_prep import (
    count_water_molecules,
    count_water_molecules_in_inputs,
    detect_last_protein_atom,
    split_docking_folder,
)


MOL2_COMPLEX = """@<TRIPOS>MOLECULE
complex_AB
 4 1 3 0 0
BIOPOLYMER
USER_CHARGES

@<TRIPOS>ATOM
      1 N   0.0 0.0 0.0 N.3 1 GLU1 0.0
      2 CA  1.0 0.0 0.0 C.3 2 ALA1 0.0
      3 C1  2.0 0.0 0.0 C.3 3 LIG1 0.0
      4 H1  3.0 0.0 0.0 H   3 LIG1 0.0
@<TRIPOS>BOND
     1    3    4 1
@<TRIPOS>SUBSTRUCTURE
     1 GLU1 1 RESIDUE 0 A GLU 0
     2 ALA1 2 RESIDUE 0 B ALA 0
     3 LIG1 3 RESIDUE 0 L LIG 0
"""


MOL2_COMPLEX_WITH_WATER_AFTER_LIGAND = """@<TRIPOS>MOLECULE
complex_water
 7 3 4 0 0
BIOPOLYMER
USER_CHARGES

@<TRIPOS>ATOM
      1 N   0.0 0.0 0.0 N.3 1 GLU1 0.0
      2 CA  1.0 0.0 0.0 C.3 2 ALA1 0.0
      3 C1  2.0 0.0 0.0 C.3 3 LIG1 0.0
      4 H1  3.0 0.0 0.0 H   3 LIG1 0.0
      5 OW  4.0 0.0 0.0 O.3 4 HOH10 0.0
      6 LP  4.1 0.0 0.0 LP  4 HOH10 0.0
      7 HW1 4.0 0.9 0.0 H   4 HOH10 0.0
@<TRIPOS>BOND
     1    3    4 1
     2    5    6 1
     3    5    7 1
@<TRIPOS>SUBSTRUCTURE
     1 GLU1 1 RESIDUE 0 A GLU 0
     2 ALA1 2 RESIDUE 0 A ALA 0
     3 LIG1 3 RESIDUE 0 L LIG 0
     4 HOH10 5 RESIDUE 0 W HOH 0
"""


class Mol2PrepTests(unittest.TestCase):
    def test_split_docking_folder_preserves_multiple_protein_chains(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            src = Path(tmp) / "src"
            src.mkdir()
            (src / "complex_AB.mol2").write_text(MOL2_COMPLEX, encoding="utf-8")

            result = split_docking_folder(src, last_pa=2)

            protein_pdb = Path(result.protein_dir) / "complex_AB.pdb"
            ligand_mol2 = Path(result.ligand_dir) / "complex_AB_ligand.mol2"
            self.assertTrue(protein_pdb.exists())
            self.assertTrue(ligand_mol2.exists())

            pdb_lines = protein_pdb.read_text(encoding="utf-8").splitlines()
            atom_lines = [line for line in pdb_lines if line.startswith("ATOM")]

            self.assertEqual(len(atom_lines), 2)
            self.assertEqual(atom_lines[0][21], "A")
            self.assertEqual(atom_lines[1][21], "B")
            self.assertEqual(atom_lines[0][76:78].strip(), "N")
            self.assertEqual(atom_lines[1][76:78].strip(), "C")

    def test_split_docking_folder_reports_progress(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            src = Path(tmp) / "src"
            src.mkdir()
            (src / "complex_1.mol2").write_text(MOL2_COMPLEX, encoding="utf-8")
            (src / "complex_2.mol2").write_text(MOL2_COMPLEX.replace("complex_AB", "complex_CD"), encoding="utf-8")

            progress_events: list[tuple[int, int, str, bool, str]] = []

            result = split_docking_folder(
                src,
                last_pa=2,
                progress_cb=lambda processed, total, filename, ok, error: progress_events.append(
                    (processed, total, filename, ok, error)
                ),
            )

            self.assertEqual(result.files_processed, 2)
            self.assertGreaterEqual(len(progress_events), 3)
            self.assertEqual(progress_events[0][:3], (0, 2, ""))
            self.assertEqual(progress_events[-1][:3], (2, 2, "complex_2.mol2"))
            self.assertTrue(progress_events[-1][3])

    def test_split_docking_folder_moves_water_to_protein_and_drops_water_lp(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            src = Path(tmp) / "src"
            src.mkdir()
            mol2 = src / "complex_water.mol2"
            mol2.write_text(MOL2_COMPLEX_WITH_WATER_AFTER_LIGAND, encoding="utf-8")

            detected = detect_last_protein_atom(mol2)
            self.assertEqual(detected.last_pa, 2)
            self.assertEqual(detected.n_protein_atoms, 4)
            self.assertEqual(detected.n_ligand_atoms, 2)
            self.assertEqual(detected.ligand_names, ["LIG1"])

            result = split_docking_folder(src, last_pa=detected.last_pa)
            self.assertEqual(result.water_molecules_detected, 1)

            protein_pdb = Path(result.protein_dir) / "complex_water.pdb"
            ligand_mol2 = Path(result.ligand_dir) / "complex_water_ligand.mol2"
            pdb_text = protein_pdb.read_text(encoding="utf-8")
            ligand_text = ligand_mol2.read_text(encoding="utf-8")

            self.assertIn("HOH", pdb_text)
            self.assertIn(" OW ", pdb_text)
            self.assertIn("HW1", pdb_text)
            water_lines = [line for line in pdb_text.splitlines() if " HOH " in line]
            self.assertTrue(water_lines)
            self.assertTrue(all(line.startswith("HETATM") for line in water_lines))
            self.assertNotIn(" LP ", pdb_text)
            self.assertNotIn("HOH", ligand_text)
            self.assertIn("2 1 1 0 0", ligand_text)
            self.assertIn("  1 1 2    1", ligand_text)

    def test_count_water_molecules_in_pdb_and_mol2_inputs(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            pdb_dir = root / "proteins"
            pdb_dir.mkdir()
            (pdb_dir / "complex_a.pdb").write_text(
                "ATOM      1  O   HOH A  10      0.000   0.000   0.000  1.00  0.00           O  \n"
                "ATOM      2  H1  HOH A  10      0.100   0.000   0.000  1.00  0.00           H  \n"
                "ATOM      3  O   WAT A  11      1.000   0.000   0.000  1.00  0.00           O  \n"
                "END\n",
                encoding="utf-8",
            )
            mol2 = root / "complex_water.mol2"
            mol2.write_text(MOL2_COMPLEX_WITH_WATER_AFTER_LIGAND, encoding="utf-8")

            self.assertEqual(count_water_molecules(pdb_dir), 2)
            self.assertEqual(count_water_molecules(mol2), 1)
            self.assertEqual(count_water_molecules_in_inputs(pdb_dir, mol2), 3)


if __name__ == "__main__":
    unittest.main()
