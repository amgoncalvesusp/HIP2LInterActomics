from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from luna_gui.core.ligand_io import (
    consolidate_folder,
    consolidate_folder_clean,
    consolidate_ligand_folder,
    consolidate_sdf_folder,
    parse_ligand_file,
    strip_hydrogens_from_mol2_file,
)


MOL2_TEMPLATE = """@<TRIPOS>MOLECULE
{mol_name}
 2 1 1 0 0
SMALL
USER_CHARGES

@<TRIPOS>ATOM
      1 C1 0.0 0.0 0.0 C.3 1 LIG 0.0
      2 H1 0.0 0.0 1.0 H   1 LIG 0.0
@<TRIPOS>BOND
     1    1    2 1
"""


SDF_TEMPLATE = """{mol_name}
  LUNA GUI

  1  0  0  0  0  0            999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
M  END
$$$$
"""


PDB_LIGAND_TEMPLATE = """COMPND    {mol_name}
HETATM    1  C1  LIG L   1       0.000   0.000   0.000  1.00  0.00           C
END
"""


class LigandIoTests(unittest.TestCase):
    def test_consolidate_folder_clean_uses_filename_stems_as_names(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            folder = Path(tmp)
            (folder / "ligand_A.mol2").write_text(
                MOL2_TEMPLATE.format(mol_name="protein_prepared"),
                encoding="utf-8",
            )
            (folder / "ligand_B.mol2").write_text(
                MOL2_TEMPLATE.format(mol_name="protein_prepared"),
                encoding="utf-8",
            )
            out = folder / "_consolidated_ligands.mol2"

            total, names = consolidate_folder_clean(
                folder,
                out,
                drop_lp=False,
                use_file_stem_as_name=True,
            )

            self.assertEqual(total, 2)
            self.assertEqual(names, ["ligand_A", "ligand_B"])
            self.assertEqual(parse_ligand_file(out), ["ligand_A", "ligand_B"])

    def test_consolidate_folder_clean_overwrites_existing_output_without_self_including(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            folder = Path(tmp)
            (folder / "ligand_A.mol2").write_text(
                MOL2_TEMPLATE.format(mol_name="ligand_A"),
                encoding="utf-8",
            )
            (folder / "_consolidated_ligands.mol2").write_text(
                MOL2_TEMPLATE.format(mol_name="stale_output"),
                encoding="utf-8",
            )
            out = folder / "_consolidated_ligands.mol2"

            total, names = consolidate_folder_clean(
                folder,
                out,
                drop_lp=False,
                use_file_stem_as_name=True,
            )

            self.assertEqual(total, 1)
            self.assertEqual(names, ["ligand_A"])
            self.assertEqual(parse_ligand_file(out), ["ligand_A"])

    def test_consolidate_folder_uses_unique_names_for_multimol_source_file(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            folder = Path(tmp)
            (folder / "bundle.mol2").write_text(
                MOL2_TEMPLATE.format(mol_name="complex_a")
                + "\n"
                + MOL2_TEMPLATE.format(mol_name="complex_b"),
                encoding="utf-8",
            )
            out = folder / "_consolidated_ligands.mol2"

            total = consolidate_folder(
                folder,
                out,
                use_file_stem_as_name=True,
            )

            self.assertEqual(total, 2)
            self.assertEqual(parse_ligand_file(out), ["bundle__1", "bundle__2"])

    def test_consolidate_folder_can_be_run_twice_without_duplicating_previous_output(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            folder = Path(tmp)
            (folder / "ligand_A.mol2").write_text(
                MOL2_TEMPLATE.format(mol_name="ligand_A"),
                encoding="utf-8",
            )
            (folder / "ligand_B.mol2").write_text(
                MOL2_TEMPLATE.format(mol_name="ligand_B"),
                encoding="utf-8",
            )
            out = folder / "_consolidated_ligands.mol2"

            first_total = consolidate_folder(folder, out, use_file_stem_as_name=False)
            second_total = consolidate_folder(folder, out, use_file_stem_as_name=False)

            self.assertEqual(first_total, 2)
            self.assertEqual(second_total, 2)
            self.assertEqual(parse_ligand_file(out), ["ligand_A", "ligand_B"])

    def test_consolidate_sdf_folder_uses_filename_stems_as_names(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            folder = Path(tmp)
            (folder / "ligand_A.sdf").write_text(
                SDF_TEMPLATE.format(mol_name="old_name"),
                encoding="utf-8",
            )
            (folder / "ligand_B.sdf").write_text(
                SDF_TEMPLATE.format(mol_name="old_name"),
                encoding="utf-8",
            )
            out = folder / "_consolidated_ligands.sdf"

            total, names = consolidate_sdf_folder(
                folder,
                out,
                use_file_stem_as_name=True,
            )

            self.assertEqual(total, 2)
            self.assertEqual(names, ["ligand_A", "ligand_B"])
            self.assertEqual(parse_ligand_file(out), ["ligand_A", "ligand_B"])

    def test_consolidate_ligand_folder_accepts_sdf_only_folder(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            folder = Path(tmp)
            (folder / "ligand_A.sdf").write_text(
                SDF_TEMPLATE.format(mol_name="original_A"),
                encoding="utf-8",
            )

            total, names, out = consolidate_ligand_folder(folder)

            self.assertEqual(total, 1)
            self.assertEqual(names, ["ligand_A"])
            self.assertEqual(out.name, "_consolidated_ligands.sdf")
            self.assertEqual(parse_ligand_file(out), ["ligand_A"])

    def test_parse_ligand_folder_accepts_pdb_ligands(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            folder = Path(tmp)
            (folder / "frame_19_ligand.pdb").write_text(
                PDB_LIGAND_TEMPLATE.format(mol_name="frame_19_LIG"),
                encoding="utf-8",
            )
            (folder / "frame_100_ligand.pdb").write_text(
                PDB_LIGAND_TEMPLATE.format(mol_name="frame_100_LIG"),
                encoding="utf-8",
            )

            self.assertEqual(
                parse_ligand_file(folder),
                ["frame_100_LIG", "frame_19_LIG"],
            )

    def test_consolidate_ligand_folder_rejects_mixed_mol2_sdf_folder(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            folder = Path(tmp)
            (folder / "ligand_A.mol2").write_text(
                MOL2_TEMPLATE.format(mol_name="ligand_A"),
                encoding="utf-8",
            )
            (folder / "ligand_B.sdf").write_text(
                SDF_TEMPLATE.format(mol_name="ligand_B"),
                encoding="utf-8",
            )

            with self.assertRaises(ValueError):
                consolidate_ligand_folder(folder)

    def test_strip_hydrogens_from_mol2_file_preserves_names_and_counts(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            folder = Path(tmp)
            src = folder / "ligands.mol2"
            src.write_text(
                MOL2_TEMPLATE.format(mol_name="ligand_A")
                + "\n"
                + MOL2_TEMPLATE.format(mol_name="ligand_B"),
                encoding="utf-8",
            )
            out = folder / "ligands_no_h.mol2"

            total, names = strip_hydrogens_from_mol2_file(src, out, drop_lp=False)

            self.assertEqual(total, 2)
            self.assertEqual(names, ["ligand_A", "ligand_B"])
            self.assertEqual(parse_ligand_file(out), ["ligand_A", "ligand_B"])
            content = out.read_text(encoding="utf-8")
            self.assertNotIn(" H ", content)
            self.assertIn("1 0 1 0 0", content)


if __name__ == "__main__":
    unittest.main()
