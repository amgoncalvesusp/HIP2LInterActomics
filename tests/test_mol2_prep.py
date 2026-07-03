from __future__ import annotations

import math
import tempfile
import unittest
from pathlib import Path

from luna_gui.core.mol2_prep import (
    count_water_molecules,
    count_water_molecules_in_inputs,
    detect_last_protein_atom,
    split_complex_folder,
    split_docking_folder,
)
from luna_gui.core.ligand_io import parse_ligand_file

try:
    from rdkit import Chem
except Exception:  # pragma: no cover - depends on optional chemistry runtime
    Chem = None


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


PDB_COMPLEX_WITH_WATER_AND_LP = """HEADER    TEST COMPLEX
ATOM      1  N   GLY A   1       0.000   0.000   0.000  1.00  0.00           N
ATOM      2  CA  GLY A   1       1.000   0.000   0.000  1.00  0.00           C
HETATM    3  O   HOH W  10       2.000   0.000   0.000  1.00  0.00           O
HETATM    4  LP  HOH W  10       2.100   0.000   0.000  1.00  0.00          LP
HETATM    5  C1  LIG L   1       3.000   0.000   0.000  1.00  0.00           C
HETATM    6  O1  LIG L   1       4.000   0.000   0.000  1.00  0.00           O
END
"""

PDB_COMPLEX_WITH_LIGAND_AS_ATOM = """HEADER    TEST COMPLEX
ATOM      1  N   GLY A   1       0.000   0.000   0.000  1.00  0.00           N
ATOM      2  CA  GLY A   1       1.000   0.000   0.000  1.00  0.00           C
ATOM      3  C1  LIG L   1       3.000   0.000   0.000  1.00  0.00           C
ATOM      4  O1  LIG L   1       4.000   0.000   0.000  1.00  0.00           O
END
"""

PDB_COMPLEX_WITH_WTM_WATER = """HEADER    TEST COMPLEX
ATOM      1  N   GLY A   7       0.000   0.000   0.000  1.00  0.00           N
ATOM      2  CA  GLY A   7       1.000   0.000   0.000  1.00  0.00           C
HETATM    3  O   WTM W  42       2.000   0.000   0.000  1.00  0.00           O
HETATM    4  C1  LIG L   5       3.000   0.000   0.000  1.00  0.00           C
HETATM    5  O1  LIG L   5       4.000   0.000   0.000  1.00  0.00           O
END
"""


def _pdb_atom(
    serial: int,
    atom_name: str,
    resname: str,
    chain: str,
    resseq: int,
    x: float,
    y: float,
    z: float,
    element: str,
    record: str = "HETATM",
    charge: str = "",
) -> str:
    return (
        f"{record:<6}{serial:5d} {atom_name:>4s} {resname:>3s} {chain:1s}{resseq:4d}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}{1.00:6.2f}{0.00:6.2f}          {element:>2s}{charge:>2s}\n"
    )


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

    def test_split_complex_folder_accepts_pdb_and_drops_water_lp(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            src = Path(tmp) / "src"
            src.mkdir()
            (src / "frame_100.pdb").write_text(PDB_COMPLEX_WITH_WATER_AND_LP, encoding="utf-8")

            result = split_complex_folder(src)

            self.assertEqual(result.files_processed, 1)
            self.assertEqual(result.proteins_written, 1)
            self.assertEqual(result.ligands_written, 1)
            self.assertEqual(result.water_molecules_detected, 1)
            protein_pdb = Path(result.protein_dir) / "frame_100.pdb"
            ligand_sdf = Path(result.ligand_dir) / "frame_100_ligand.sdf"
            self.assertTrue(protein_pdb.exists())
            self.assertTrue(ligand_sdf.exists())
            self.assertFalse((Path(result.ligand_dir) / "frame_100_ligand.pdb").exists())

            protein_text = protein_pdb.read_text(encoding="utf-8")
            ligand_text = ligand_sdf.read_text(encoding="utf-8")
            self.assertIn("HOH", protein_text)
            self.assertNotIn(" LP ", protein_text)
            self.assertNotIn("LIG", protein_text)
            self.assertEqual(parse_ligand_file(ligand_sdf), ["frame_100_LIG"])
            self.assertIn("M  END", ligand_text)
            self.assertIn("$$$$", ligand_text)
            self.assertNotIn("HOH", ligand_text)

    def test_split_complex_folder_detects_pdb_ligand_even_when_written_as_atom(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            src = Path(tmp) / "src"
            src.mkdir()
            (src / "pose_19.pdb").write_text(PDB_COMPLEX_WITH_LIGAND_AS_ATOM, encoding="utf-8")

            result = split_complex_folder(src)

            self.assertEqual(result.proteins_written, 1)
            self.assertEqual(result.ligands_written, 1)
            protein_text = (Path(result.protein_dir) / "pose_19.pdb").read_text(encoding="utf-8")
            ligand_sdf = Path(result.ligand_dir) / "pose_19_ligand.sdf"
            ligand_text = ligand_sdf.read_text(encoding="utf-8")
            self.assertIn(" GLY ", protein_text)
            self.assertNotIn(" LIG ", protein_text)
            self.assertEqual(parse_ligand_file(ligand_sdf), ["pose_19_LIG"])
            self.assertIn("M  END", ligand_text)
            self.assertTrue(any("tratado como ligante" in msg for msg in result.errors))
            log_text = Path(result.log_file).read_text(encoding="utf-8")
            self.assertIn("tratado como ligante", log_text)

    def test_split_complex_folder_treats_wtm_as_water_and_preserves_pdb_fields(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            src = Path(tmp) / "src"
            src.mkdir()
            (src / "pose_wtm.pdb").write_text(PDB_COMPLEX_WITH_WTM_WATER, encoding="utf-8")

            result = split_complex_folder(src)

            self.assertEqual(result.water_molecules_detected, 1)
            protein_text = (Path(result.protein_dir) / "pose_wtm.pdb").read_text(encoding="utf-8")
            ligand_sdf = Path(result.ligand_dir) / "pose_wtm_ligand.sdf"
            ligand_text = ligand_sdf.read_text(encoding="utf-8")
            self.assertIn(" GLY A   7", protein_text)
            self.assertIn(" WTM W  42", protein_text)
            self.assertNotIn(" WTM ", ligand_text)
            self.assertTrue(Path(result.log_file).exists())

    @unittest.skipIf(Chem is None, "RDKit is required to inspect SDF chemistry")
    def test_split_complex_folder_preserves_pdb_formal_charge_in_sdf(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            src = Path(tmp) / "src"
            src.mkdir()
            lines = [
                "HEADER    TEST COMPLEX\n",
                _pdb_atom(1, "CA", "GLY", "A", 1, 0.0, 0.0, 0.0, "C", record="ATOM"),
                _pdb_atom(10, "N1", "AMM", "L", 1, 5.0, 0.0, 0.0, "N", charge="1+"),
                _pdb_atom(11, "H1", "AMM", "L", 1, 6.0, 0.0, 0.0, "H"),
                _pdb_atom(12, "H2", "AMM", "L", 1, 4.0, 0.0, 0.0, "H"),
                _pdb_atom(13, "H3", "AMM", "L", 1, 5.0, 1.0, 0.0, "H"),
                _pdb_atom(14, "H4", "AMM", "L", 1, 5.0, 0.0, 1.0, "H"),
                "END\n",
            ]
            (src / "ammonium.pdb").write_text("".join(lines), encoding="utf-8")

            result = split_complex_folder(src)

            ligand_sdf = Path(result.ligand_dir) / "ammonium_ligand.sdf"
            mol = Chem.SDMolSupplier(str(ligand_sdf), removeHs=False)[0]
            self.assertIsNotNone(mol)
            self.assertEqual(sum(atom.GetFormalCharge() for atom in mol.GetAtoms()), 1)
            self.assertEqual(mol.GetProp("HIP2L_TotalFormalCharge"), "1")
            self.assertTrue(mol.HasProp("HIP2L_GasteigerCharges"))

    @unittest.skipIf(Chem is None, "RDKit is required to inspect SDF chemistry")
    def test_split_complex_folder_perceives_aromaticity_in_pdb_ligand(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            src = Path(tmp) / "src"
            src.mkdir()
            lines = [
                "HEADER    TEST COMPLEX\n",
                _pdb_atom(1, "CA", "GLY", "A", 1, 0.0, 0.0, 0.0, "C", record="ATOM"),
            ]
            for idx in range(6):
                angle = 2.0 * math.pi * idx / 6.0
                lines.append(
                    _pdb_atom(
                        10 + idx,
                        f"C{idx + 1}",
                        "BEN",
                        "L",
                        1,
                        5.0 + 1.397 * math.cos(angle),
                        1.397 * math.sin(angle),
                        0.0,
                        "C",
                    )
                )
            for idx in range(6):
                angle = 2.0 * math.pi * idx / 6.0
                lines.append(
                    _pdb_atom(
                        20 + idx,
                        f"H{idx + 1}",
                        "BEN",
                        "L",
                        1,
                        5.0 + 2.48 * math.cos(angle),
                        2.48 * math.sin(angle),
                        0.0,
                        "H",
                    )
                )
            lines.append("END\n")
            (src / "benzene.pdb").write_text("".join(lines), encoding="utf-8")

            result = split_complex_folder(src)

            ligand_sdf = Path(result.ligand_dir) / "benzene_ligand.sdf"
            mol = Chem.SDMolSupplier(str(ligand_sdf), removeHs=False)[0]
            self.assertIsNotNone(mol)
            self.assertEqual(sum(1 for bond in mol.GetBonds() if bond.GetIsAromatic()), 6)
            self.assertTrue(mol.HasProp("HIP2L_GasteigerCharges"))


if __name__ == "__main__":
    unittest.main()
