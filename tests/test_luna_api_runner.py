from __future__ import annotations

import ast
import inspect
import re
import tempfile
import unittest
import json
from pathlib import Path

from luna_gui.core.analysis_runtime import (
    _FP_DETAIL_SCRIPT,
    _FP_SESSION_SCRIPT,
    _PSE_FILTER_SCRIPT,
    _clean_helper_text,
)
from luna_gui.core.luna_api_runner import (
    API_RUNNER_SCRIPT,
    build_entry_specs,
    ligand_mol_obj_type,
    protein_has_explicit_hydrogens,
    protein_is_gui_preprocessed,
    read_ifp_seed_file,
    resolve_protein_processing_flags,
    should_use_api_runner,
    validate_entry_specs,
    validate_hydrogen_inputs,
    write_params,
)
from luna_gui.core.luna_runner import validate as validate_cli_inputs
from luna_gui.core.project import IFP_ALL, ProjectConfig


class LunaApiRunnerTests(unittest.TestCase):
    def test_ligand_atom_map_reuses_exact_matrix_identifiers(self) -> None:
        self.assertIn("def _align_atom_map_labels", API_RUNNER_SCRIPT)
        self.assertIn(
            "matrix_metadata=matrix_metadata",
            API_RUNNER_SCRIPT,
        )
        self.assertIn('"labels_match_matrix": all(label in labels for label in matrix_labels)', API_RUNNER_SCRIPT)
        self.assertIn('"atom_label_mapping": label_mapping', API_RUNNER_SCRIPT)

    def test_ligand_atom_alignment_uses_serials_for_interacting_subset(self) -> None:
        tree = ast.parse(API_RUNNER_SCRIPT)
        names = {"_normalized_atom_label", "_numeric_atom_token", "_align_atom_map_labels"}
        functions = [
            node for node in tree.body
            if isinstance(node, ast.FunctionDef) and node.name in names
        ]
        namespace = {"re": __import__("re")}
        exec(compile(ast.Module(body=functions, type_ignores=[]), "<atom-alignment>", "exec"), namespace)

        labels, mapping = namespace["_align_atom_map_labels"](
            ["C1", "C2", "N3", "O4"],
            ["LIG-CARBON", "N3"],
            source_records=[
                {"source_index": 0, "serial": 1},
                {"source_index": 1, "serial": 2},
                {"source_index": 2, "serial": 3},
                {"source_index": 3, "serial": 4},
            ],
            matrix_metadata={"LIG-CARBON": {"serial": 2}, "N3": {"serial": 3}},
            return_mapping=True,
        )

        self.assertEqual(labels, ["C1", "LIG-CARBON", "N3", "O4"])
        self.assertTrue(mapping[1]["matched"])
        self.assertEqual(mapping[1]["matrix_label"], "LIG-CARBON")

    def test_project_config_uses_python_api_for_fork_and_multi_ifp(self) -> None:
        self.assertTrue(ProjectConfig(fork_from="D:/old_project").uses_python_api())
        self.assertTrue(ProjectConfig(ifp_type=IFP_ALL).uses_python_api())
        self.assertTrue(ProjectConfig(ifp_seed_file="D:/seed.txt").uses_python_api())
        self.assertTrue(ProjectConfig(interaction_config_file="D:/config.cfg").uses_python_api())
        self.assertTrue(ProjectConfig(inter_max_distance_cap=4.5).uses_python_api())

    def test_protein_heteroatom_option_defaults_to_off_and_uses_api_runner(self) -> None:
        self.assertFalse(ProjectConfig().include_protein_heteroatoms)
        self.assertTrue(
            ProjectConfig(include_protein_heteroatoms=True).uses_python_api()
        )

    def test_receptor_heteroatom_role_is_opt_in_and_keeps_chain_residue_id(self) -> None:
        tree = ast.parse(API_RUNNER_SCRIPT)
        required = {
            "_group_has_water_residue",
            "_residue_name_from_atom",
            "_residue_key_from_atom",
            "_group_has_ligand_residue",
            "_group_is_configured_protein_heteroatom",
            "_group_has_structural_residue",
            "_group_role",
            "_residue_label",
            "_labels_from_group",
            "_interaction_residue_payload",
        }
        functions = [
            node for node in tree.body
            if isinstance(node, ast.FunctionDef) and node.name in required
        ]
        namespace = {
            "AA_RESIDUES": {"ALA", "GLY"},
            "WATER_RESIDUES": {"HOH"},
            "INCLUDE_PROTEIN_HETEROATOMS": False,
            "PROTEIN_HETEROATOM_RESIDUES": {"A/ZN/228"},
        }
        exec(compile(ast.Module(body=functions, type_ignores=[]), "<roles>", "exec"), namespace)

        class Chain:
            id = "A"

        class Residue:
            parent = Chain()
            resname = "ZN"
            id = ("H_ZN", 228, " ")

        class Atom:
            parent = Residue()

        class Group:
            atoms = [Atom()]

            def has_water(self):
                return False

            def has_residue(self):
                return False

            def has_nucleotide(self):
                return False

            def has_hetatm(self):
                return True

        class LigandGroup:
            atoms = []

            def has_water(self):
                return False

            def has_residue(self):
                return False

            def has_nucleotide(self):
                return False

            def has_hetatm(self):
                return True

        class Interaction:
            src_grp = LigandGroup()
            trgt_grp = Group()

        zinc = Group()
        self.assertEqual(namespace["_group_role"](zinc), "ligand")
        self.assertEqual(namespace["_interaction_residue_payload"](Interaction()), {})

        namespace["INCLUDE_PROTEIN_HETEROATOMS"] = True
        namespace["PROTEIN_HETEROATOM_RESIDUES"] = set()
        self.assertEqual(namespace["_group_role"](zinc), "ligand")
        namespace["PROTEIN_HETEROATOM_RESIDUES"] = {"A/ZN/228"}
        self.assertEqual(namespace["_group_role"](zinc), "protein")
        self.assertEqual(
            namespace["_interaction_residue_payload"](Interaction()),
            {"A/ZN/228": "protein_heteroatom"},
        )

    def test_intermediate_artifact_records_heteroatom_roles(self) -> None:
        self.assertIn('"residue_roles"', API_RUNNER_SCRIPT)
        self.assertIn("ignore_res_hetatm=False", API_RUNNER_SCRIPT)

    def test_api_runner_rebuilds_similarity_from_ifp_instead_of_using_luna_internal_queue(self) -> None:
        self.assertIn("proj.ifp_sim_matrix_output = None", API_RUNNER_SCRIPT)
        self.assertIn("_write_similarity_outputs_from_ifp(", API_RUNNER_SCRIPT)

    def test_api_runner_calculates_similarity_once_per_ifp_type(self) -> None:
        tree = ast.parse(API_RUNNER_SCRIPT)
        generator = next(node for node in tree.body if isinstance(node, ast.FunctionDef) and node.name == "_generate_ifps")
        calls = [
            node for node in ast.walk(generator)
            if isinstance(node, ast.Call)
            and isinstance(node.func, ast.Name)
            and node.func.id == "_write_similarity_outputs_from_ifp"
        ]
        self.assertEqual(len(calls), 1)

    def test_api_runner_saves_ifp_seed_files(self) -> None:
        self.assertIn("seed_ifp_{suffix}_importance.txt", API_RUNNER_SCRIPT)
        self.assertIn("_write_ifp_seed", API_RUNNER_SCRIPT)

    def test_read_ifp_seed_file_reads_first_integer(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "seed.txt"
            path.write_text("seed = 123\n", encoding="utf-8")

            self.assertEqual(read_ifp_seed_file(path), 123)

    def test_ligand_backend_uses_rdkit_for_sdf_and_openbabel_for_mol2(self) -> None:
        self.assertEqual(ligand_mol_obj_type("ligands.sdf"), "rdkit")
        self.assertEqual(ligand_mol_obj_type("ligands.sd"), "rdkit")
        self.assertEqual(ligand_mol_obj_type("ligands.mol"), "rdkit")
        self.assertEqual(ligand_mol_obj_type("ligands.mol2"), "openbabel")
        self.assertIn("lig_mol_obj_type", API_RUNNER_SCRIPT)

    def test_api_runner_suppresses_rdkit_deprecation_warnings(self) -> None:
        self.assertIn('RDLogger.DisableLog("rdApp.warning")', API_RUNNER_SCRIPT)

    def test_each_exported_pse_resets_the_pymol_object_store(self) -> None:
        self.assertIn("def _reset_pymol_session", API_RUNNER_SCRIPT)
        self.assertIn("cmd.reinitialize()", API_RUNNER_SCRIPT)
        self.assertIn("_reset_pymol_session()\n                viewer = InteractionViewer", API_RUNNER_SCRIPT)
        self.assertIn("finally:\n                _reset_pymol_session()", API_RUNNER_SCRIPT)
        self.assertIn("def _reset_pymol_session", _PSE_FILTER_SCRIPT)
        self.assertIn("viewer = InteractionViewer(show_hydrop_surface=False)", _PSE_FILTER_SCRIPT)

    def test_pse_palette_is_applied_in_the_active_session_without_reloading(self) -> None:
        palette_function = next(
            node
            for node in ast.parse(API_RUNNER_SCRIPT).body
            if isinstance(node, ast.FunctionDef)
            and node.name == "_apply_pse_interaction_colors"
        )
        calls_load = [
            node
            for node in ast.walk(palette_function)
            if isinstance(node, ast.Call)
            and isinstance(node.func, ast.Attribute)
            and node.func.attr == "load"
        ]
        self.assertEqual(calls_load, [])
        self.assertNotIn("objetos coloridos", API_RUNNER_SCRIPT)

    def test_prepared_protein_flags_keep_add_h_without_staging_inputs(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            protein_file = Path(tmp) / "protein_prepared.pdb"
            protein_file.write_text(
                "REMARK   Separated Protein\n"
                "ATOM      1  H   GLY A   1      0.000   0.000   0.000  1.00  0.00           H  \n"
                "END\n",
                encoding="utf-8",
            )

            cfg = ProjectConfig(protein_file=str(protein_file))
            flags = resolve_protein_processing_flags(cfg)

            self.assertTrue(protein_has_explicit_hydrogens(protein_file))
            self.assertTrue(protein_is_gui_preprocessed(protein_file))
            self.assertTrue(flags["add_h"])
            self.assertFalse(flags["amend_mol"])
            self.assertFalse(flags["stage_protein_without_h"])
            self.assertTrue(should_use_api_runner(cfg))

    def test_plain_protein_flags_keep_default_luna_processing(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            protein_file = Path(tmp) / "protein_plain.pdb"
            protein_file.write_text(
                "HEADER    TEST\n"
                "ATOM      1  CA  GLY A   1      0.000   0.000   0.000  1.00  0.00           C  \n"
                "END\n",
                encoding="utf-8",
            )

            cfg = ProjectConfig(protein_file=str(protein_file))
            flags = resolve_protein_processing_flags(cfg)

            self.assertFalse(protein_has_explicit_hydrogens(protein_file))
            self.assertFalse(protein_is_gui_preprocessed(protein_file))
            self.assertTrue(flags["add_h"])
            self.assertTrue(flags["amend_mol"])
            self.assertFalse(flags["stage_protein_without_h"])

    def test_build_entry_specs_uses_selected_protein_when_waters_are_off(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            protein_file = Path(tmp) / "reference_complex.pdb"
            protein_file.write_text("HEADER\n", encoding="utf-8")

            cfg = ProjectConfig(protein_file=str(protein_file), include_waters=False)
            specs = build_entry_specs(cfg, ["ligA", "ligB"])

            self.assertEqual(
                specs,
                [
                    {"pdb_id": "reference_complex", "ligand_name": "ligA"},
                    {"pdb_id": "reference_complex", "ligand_name": "ligB"},
                ],
            )

    def test_build_entry_specs_pairs_each_ligand_with_matching_pdb(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            protein_dir = Path(tmp)
            for name in ("CHEMBL112640", "CHEMBL2177539_ncw"):
                (protein_dir / f"{name}.pdb").write_text("HEADER\n", encoding="utf-8")

            cfg = ProjectConfig(
                protein_file=str(protein_dir),
                include_waters=True,
            )
            specs = build_entry_specs(cfg, ["CHEMBL112640_LIG", "CHEMBL2177539_ncw_LIG"])

            self.assertEqual(
                specs,
                [
                    {"pdb_id": "CHEMBL112640", "ligand_name": "CHEMBL112640_LIG"},
                    {"pdb_id": "CHEMBL2177539_ncw", "ligand_name": "CHEMBL2177539_ncw_LIG"},
                ],
            )

    def test_build_entry_specs_accepts_protein_folder_when_waters_are_on(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            protein_dir = Path(tmp) / "proteins"
            protein_dir.mkdir()
            for name in ("frame_19", "frame_100"):
                (protein_dir / f"{name}.pdb").write_text("HEADER\n", encoding="utf-8")

            cfg = ProjectConfig(
                protein_file=str(protein_dir),
                include_waters=True,
            )
            specs = build_entry_specs(cfg, ["frame_100_LIG", "frame_19_LIG"])

            self.assertEqual(
                specs,
                [
                    {"pdb_id": "frame_100", "ligand_name": "frame_100_LIG"},
                    {"pdb_id": "frame_19", "ligand_name": "frame_19_LIG"},
                ],
            )

    def test_build_entry_specs_accepts_protein_folder_without_waters(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            protein_dir = Path(tmp) / "proteins"
            protein_dir.mkdir()
            for name in ("frame_19", "frame_100"):
                (protein_dir / f"{name}.pdb").write_text("HEADER\n", encoding="utf-8")

            cfg = ProjectConfig(
                protein_file=str(protein_dir),
                include_waters=False,
            )
            specs = build_entry_specs(cfg, ["frame_100_LIG", "frame_19_LIG"])

            self.assertEqual(
                specs,
                [
                    {"pdb_id": "frame_100", "ligand_name": "frame_100_LIG"},
                    {"pdb_id": "frame_19", "ligand_name": "frame_19_LIG"},
                ],
            )
            self.assertTrue(should_use_api_runner(cfg))

    def test_build_entry_specs_accepts_ligand_folder_with_single_files(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            protein_dir = root / "proteins"
            ligand_dir = root / "ligands"
            protein_dir.mkdir()
            ligand_dir.mkdir()
            for name in ("frame_19", "frame_100"):
                (protein_dir / f"{name}.pdb").write_text("HEADER\n", encoding="utf-8")
                (ligand_dir / f"{name}_ligand.pdb").write_text(f"COMPND    {name}_LIG\nEND\n", encoding="utf-8")

            cfg = ProjectConfig(
                protein_file=str(protein_dir),
                ligand_file=str(ligand_dir),
            )
            specs = build_entry_specs(cfg, ["frame_100_LIG", "frame_19_LIG"])

            self.assertEqual(specs[0]["pdb_id"], "frame_100")
            self.assertEqual(specs[0]["mol_file"], str(ligand_dir / "frame_100_ligand.pdb"))
            self.assertEqual(specs[0]["mol_obj_type"], "openbabel")
            self.assertFalse(specs[0]["is_multimol_file"])

    def test_validate_entry_specs_reports_missing_matching_pdb(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            protein_dir = Path(tmp)
            (protein_dir / "CHEMBL112640.pdb").write_text("HEADER\n", encoding="utf-8")

            cfg = ProjectConfig(
                protein_file=str(protein_dir),
                include_waters=True,
            )
            errs = validate_entry_specs(cfg, ["CHEMBL999999_LIG"])

            self.assertEqual(len(errs), 1)
            self.assertIn("CHEMBL999999_LIG", errs[0])

    def test_write_params_expands_all_ifp_outputs_and_keeps_fork_source(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            workdir = root / "forked_project"
            protein_file = root / "proteins" / "complex_A.pdb"
            protein_file.parent.mkdir(parents=True, exist_ok=True)
            protein_file.write_text("HEADER\n", encoding="utf-8")
            ligand_file = root / "ligands.mol2"
            ligand_file.write_text("@<TRIPOS>MOLECULE\nligA\n", encoding="utf-8")
            seed_file = root / "seed_ifp.txt"
            seed_file.write_text("321\n", encoding="utf-8")

            cfg = ProjectConfig(
                protein_file=str(protein_file),
                ligand_file=str(ligand_file),
                workdir=str(workdir),
                fork_from=str(root / "source_project"),
                interaction_config_file=str(root / "custom_interactions.cfg"),
                inter_max_distance_cap=4.5,
                out_ifp=True,
                ifp_type=IFP_ALL,
                ifp_seed_file=str(seed_file),
                sim_matrix=True,
            )

            params_path = write_params(workdir, cfg, ["ligA"])
            params = json.loads(Path(params_path).read_text(encoding="utf-8"))

            self.assertEqual(params["fork_from"], str(root / "source_project"))
            self.assertEqual(params["interaction_config_file"], str(root / "custom_interactions.cfg"))
            self.assertEqual(params["inter_max_distance_cap"], 4.5)
            self.assertFalse(params["include_protein_heteroatoms"])
            self.assertTrue(params["add_h"])
            self.assertTrue(params["amend_mol"])
            self.assertEqual(params["ifp_seed"], 321)
            self.assertEqual(params["ifp_seed_file"], str(seed_file))
            self.assertEqual(params["ifp_types"], ["HIFP", "EIFP", "FIFP"])
            self.assertEqual(
                params["ifp_outputs"],
                {
                    "HIFP": str(workdir / "results" / "fingerprints" / "ifp_H.csv"),
                    "EIFP": str(workdir / "results" / "fingerprints" / "ifp_E.csv"),
                    "FIFP": str(workdir / "results" / "fingerprints" / "ifp_F.csv"),
                },
            )
            self.assertEqual(
                params["sim_matrix_outputs"],
                {
                    "HIFP": str(workdir / "sim_matrix_H.csv"),
                    "EIFP": str(workdir / "sim_matrix_E.csv"),
                    "FIFP": str(workdir / "sim_matrix_F.csv"),
                },
            )

    def test_write_params_persists_protein_heteroatom_option(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            protein_file = root / "protein.pdb"
            ligand_file = root / "ligands.mol2"
            protein_file.write_text(
                "HEADER\n"
                "HETATM 3594  ZN   ZN A 228       3.425   8.551   3.392  1.00  0.00          ZN2+  \n",
                encoding="utf-8",
            )
            ligand_file.write_text("@<TRIPOS>MOLECULE\nligA\n", encoding="utf-8")
            cfg = ProjectConfig(
                protein_file=str(protein_file),
                ligand_file=str(ligand_file),
                workdir=str(root / "workdir"),
                include_protein_heteroatoms=True,
            )

            params = json.loads(Path(write_params(cfg.workdir, cfg, ["ligA"])).read_text(encoding="utf-8"))

        self.assertTrue(params["include_protein_heteroatoms"])
        self.assertEqual(params["protein_heteroatom_residues"], ["A/ZN/228"])

    def test_write_params_uses_typed_similarity_output_for_single_ifp(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            workdir = root / "project"
            protein_file = root / "protein.pdb"
            protein_file.write_text("HEADER\n", encoding="utf-8")
            ligand_file = root / "ligands.mol2"
            ligand_file.write_text("@<TRIPOS>MOLECULE\nligA\n", encoding="utf-8")

            cfg = ProjectConfig(
                protein_file=str(protein_file),
                ligand_file=str(ligand_file),
                workdir=str(workdir),
                out_ifp=True,
                ifp_type="EIFP",
                sim_matrix=True,
            )

            params_path = write_params(workdir, cfg, ["ligA"])
            params = json.loads(Path(params_path).read_text(encoding="utf-8"))

            self.assertEqual(
                params["sim_matrix_outputs"],
                {"EIFP": str(workdir / "sim_matrix_E.csv")},
            )

    def test_write_params_selects_rdkit_backend_for_sdf_ligands(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            workdir = root / "project"
            protein_file = root / "protein.pdb"
            protein_file.write_text("HEADER\n", encoding="utf-8")
            ligand_file = root / "ligands.sdf"
            ligand_file.write_text("ligA\n  LUNA GUI\n\n  0  0  0  0  0  0            999 V2000\nM  END\n$$$$\n", encoding="utf-8")

            cfg = ProjectConfig(
                protein_file=str(protein_file),
                ligand_file=str(ligand_file),
                workdir=str(workdir),
            )

            params_path = write_params(workdir, cfg, ["ligA"])
            params = json.loads(Path(params_path).read_text(encoding="utf-8"))

            self.assertEqual(params["lig_mol_obj_type"], "rdkit")

    def test_write_params_uses_selected_protein_folder_for_hydrated_inputs(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            workdir = root / "project"
            protein_dir = root / "proteins"
            protein_dir.mkdir()
            (protein_dir / "ligA.pdb").write_text("HEADER\n", encoding="utf-8")
            ligand_file = root / "ligands.mol2"
            ligand_file.write_text("@<TRIPOS>MOLECULE\nligA_LIG\n", encoding="utf-8")

            cfg = ProjectConfig(
                protein_file=str(protein_dir),
                ligand_file=str(ligand_file),
                workdir=str(workdir),
                include_waters=True,
            )

            params_path = write_params(workdir, cfg, ["ligA_LIG"])
            params = json.loads(Path(params_path).read_text(encoding="utf-8"))

            self.assertEqual(params["pdb_dir"], str(protein_dir))
            self.assertEqual(params["entry_specs"], [{"pdb_id": "ligA", "ligand_name": "ligA_LIG"}])

    def test_validate_accepts_protein_folder_for_per_complex_mode(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            protein_dir = root / "proteins"
            protein_dir.mkdir()
            (protein_dir / "ligA.pdb").write_text("HEADER\n", encoding="utf-8")
            ligand_file = root / "ligands.mol2"
            ligand_file.write_text("@<TRIPOS>MOLECULE\nligA_LIG\n", encoding="utf-8")
            cfg = ProjectConfig(
                protein_file=str(protein_dir),
                ligand_file=str(ligand_file),
                workdir=str(root / "project"),
                selected_ligands=["ligA_LIG"],
            )

            self.assertEqual(validate_cli_inputs(cfg), [])

    def test_write_params_keeps_original_receptor_and_requested_ph(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            workdir = root / "project"
            protein_file = root / "proteins" / "complex_A.pdb"
            protein_file.parent.mkdir(parents=True, exist_ok=True)
            protein_file.write_text(
                "REMARK   Separated Protein\n"
                "ATOM      1  H1  GLY A   1      0.000   0.000   0.000  1.00  0.00           H  \n"
                "ATOM      2  CA  GLY A   1      1.000   0.000   0.000  1.00  0.00           C  \n"
                "END\n",
                encoding="utf-8",
            )
            ligand_file = root / "ligands.mol2"
            ligand_file.write_text("@<TRIPOS>MOLECULE\nligA\n", encoding="utf-8")
            cfg = ProjectConfig(
                protein_file=str(protein_file),
                ligand_file=str(ligand_file),
                workdir=str(workdir),
                ph=6.8,
            )

            params_path = write_params(workdir, cfg, ["ligA"])
            params = json.loads(Path(params_path).read_text(encoding="utf-8"))

            self.assertTrue(params["add_h"])
            self.assertEqual(params["ph"], 6.8)
            self.assertFalse(params["stage_protein_without_h"])
            self.assertEqual(params["pdb_id"], "complex_A")
            self.assertEqual(params["pdb_dir"], str(protein_file.parent))
            self.assertEqual(params["lig_file"], str(ligand_file))

    def test_write_params_preserves_inputs_when_add_h_is_disabled(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            workdir = root / "project"
            protein_file = root / "proteins" / "complex_A.pdb"
            protein_file.parent.mkdir(parents=True, exist_ok=True)
            protein_file.write_text(
                "HEADER    TEST\n"
                "ATOM      1  H1  GLY A   1      0.000   0.000   0.000  1.00  0.00           H  \n"
                "ATOM      2  CA  GLY A   1      0.000   0.000   0.000  1.00  0.00           C  \n"
                "END\n",
                encoding="utf-8",
            )
            ligand_file = root / "ligands.mol2"
            ligand_file.write_text(
                "@<TRIPOS>MOLECULE\n"
                "ligA\n"
                " 2 1 1 0 0\n"
                "SMALL\n"
                "USER_CHARGES\n\n"
                "@<TRIPOS>ATOM\n"
                "      1 C1 0.0 0.0 0.0 C.3 1 LIG 0.0\n"
                "      2 H1 0.0 0.0 1.0 H   1 LIG 0.0\n"
                "@<TRIPOS>BOND\n"
                "     1    1    2 1\n",
                encoding="utf-8",
            )
            cfg = ProjectConfig(
                protein_file=str(protein_file),
                ligand_file=str(ligand_file),
                workdir=str(workdir),
                add_h=False,
            )

            params_path = write_params(workdir, cfg, ["ligA"])
            params = json.loads(Path(params_path).read_text(encoding="utf-8"))

            self.assertFalse(params["add_h"])
            self.assertFalse(params["stage_ligand_without_h"])
            self.assertEqual(params["lig_file"], str(ligand_file))
            self.assertEqual(params["pdb_dir"], str(protein_file.parent))

    def test_validate_hydrogen_inputs_warns_when_add_h_is_disabled_without_hydrogens(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            protein_file = root / "proteins" / "complex_A.pdb"
            protein_file.parent.mkdir(parents=True, exist_ok=True)
            protein_file.write_text(
                "HEADER    TEST\n"
                "ATOM      1  CA  GLY A   1      0.000   0.000   0.000  1.00  0.00           C  \n"
                "END\n",
                encoding="utf-8",
            )
            ligand_file = root / "ligands.mol2"
            ligand_file.write_text(
                "@<TRIPOS>MOLECULE\n"
                "ligA\n"
                " 1 0 1 0 0\n"
                "SMALL\n"
                "USER_CHARGES\n\n"
                "@<TRIPOS>ATOM\n"
                "      1 C1 0.0 0.0 0.0 C.3 1 LIG 0.0\n",
                encoding="utf-8",
            )
            cfg = ProjectConfig(
                protein_file=str(protein_file),
                ligand_file=str(ligand_file),
                add_h=False,
            )

            errs = validate_hydrogen_inputs(cfg)

            self.assertGreaterEqual(len(errs), 2)
            self.assertIn("Add_H", errs[0])
            self.assertIn("ligantes", "\n".join(errs))

    def test_api_runner_uses_lightweight_fp_shell_payloads(self) -> None:
        self.assertIn("entry_meta", API_RUNNER_SCRIPT)
        self.assertIn("feature_shells", API_RUNNER_SCRIPT)
        self.assertIn("_save_feature_shell_payload", API_RUNNER_SCRIPT)
        self.assertNotIn('"shell_manager": sm', API_RUNNER_SCRIPT)

    def test_api_runner_exposes_canonical_fp_class_labels(self) -> None:
        self.assertIn("Ligand's level 0 features only", API_RUNNER_SCRIPT)
        self.assertIn("Protein's level 0 features only", API_RUNNER_SCRIPT)
        self.assertIn("Upper level with ligand atomic information only", API_RUNNER_SCRIPT)
        self.assertIn("Upper level with protein atomic information only", API_RUNNER_SCRIPT)
        self.assertIn("Intraligand interactions only", API_RUNNER_SCRIPT)
        self.assertIn("Intraprotein interactions only", API_RUNNER_SCRIPT)
        self.assertIn("Has noncovalent interactions with the protein", API_RUNNER_SCRIPT)
        self.assertIn("Features with collision in the same complex", API_RUNNER_SCRIPT)
        self.assertIn("Unreliable feature", API_RUNNER_SCRIPT)

    def test_api_runner_keeps_water_mediated_protein_context(self) -> None:
        self.assertIn("_group_has_water_residue", API_RUNNER_SCRIPT)
        self.assertIn("_group_has_water_residue", _FP_DETAIL_SCRIPT)
        self.assertIn('"WTM"', _FP_DETAIL_SCRIPT)
        self.assertIn("has_ligand_water", API_RUNNER_SCRIPT)
        self.assertIn("has_protein_water", API_RUNNER_SCRIPT)
        self.assertIn("has_water_mediated_protein_context", API_RUNNER_SCRIPT)
        self.assertIn('{"ligand", "water"}', _FP_DETAIL_SCRIPT)
        self.assertIn('{"protein", "water"}', _FP_DETAIL_SCRIPT)

    def test_api_runner_marks_collision_only_for_mixed_class_signatures(self) -> None:
        self.assertIn("def _has_mixed_class_collision", API_RUNNER_SCRIPT)
        self.assertIn("raw_collision and _has_mixed_class_collision", API_RUNNER_SCRIPT)
        self.assertNotIn("_classify_shell_natures(shell, collision=collision)", API_RUNNER_SCRIPT)

    def test_api_runner_exports_fp_shell_levels(self) -> None:
        self.assertIn("def _shell_level_key", API_RUNNER_SCRIPT)
        self.assertIn('"shell_levels"', API_RUNNER_SCRIPT)
        self.assertIn('"shell_level_breakdown"', API_RUNNER_SCRIPT)
        self.assertIn('"collision_shell_levels"', API_RUNNER_SCRIPT)
        self.assertIn('"collision_level_breakdown"', API_RUNNER_SCRIPT)
        self.assertIn('"raw_collision"', API_RUNNER_SCRIPT)

    def test_helper_output_filter_keeps_real_errors(self) -> None:
        noisy = (
            "/home/user/.conda/envs/luna-env/lib/python3.9/site-packages/openbabel/__init__.py:14: "
            'UserWarning: "import openbabel" is deprecated\n'
            "Module 'simplejson' not available. Built-in module 'json' will be imported.\n"
            "RuntimeError: real failure\n"
        )

        cleaned = _clean_helper_text(noisy)

        self.assertNotIn("openbabel", cleaned)
        self.assertNotIn("simplejson", cleaned)
        self.assertIn("real failure", cleaned)

    def test_fp_session_helper_restores_entries_from_metadata(self) -> None:
        self.assertIn("_restore_entry", _FP_SESSION_SCRIPT)
        self.assertIn("entry_meta", _FP_SESSION_SCRIPT)
        self.assertIn("feature_shells", _FP_SESSION_SCRIPT)

    def test_fp_session_helper_prefers_live_shell_regeneration(self) -> None:
        self.assertIn("_regenerate_shells_from_project", _FP_SESSION_SCRIPT)
        self.assertIn("_load_cached_shell_payload", _FP_SESSION_SCRIPT)
        self.assertLess(
            _FP_SESSION_SCRIPT.index("_regenerate_shells_from_project"),
            _FP_SESSION_SCRIPT.index("_load_cached_shell_payload"),
        )

    def test_fp_session_helper_adds_shell_number_labels(self) -> None:
        self.assertIn("def _add_shell_number_labels", _FP_SESSION_SCRIPT)
        self.assertIn("Shell {shell_index} | L{level}", _FP_SESSION_SCRIPT)
        self.assertIn('"shell_labels"', _FP_SESSION_SCRIPT)
        self.assertIn("cmd.pseudoatom", _FP_SESSION_SCRIPT)

    def test_generated_pse_sessions_use_the_canonical_interaction_palette(self) -> None:
        self.assertIn(
            '"pse_interaction_colors": dict(INTERACTION_COLORS)',
            inspect.getsource(write_params),
        )
        self.assertIn("def _apply_pse_interaction_colors", API_RUNNER_SCRIPT)
        self.assertIn('cmd.set("dash_color", color_name, object_name)', API_RUNNER_SCRIPT)
        self.assertIn("def _save_viewer_session_with_palette", API_RUNNER_SCRIPT)
        self.assertIn("cmd.save = _save_with_palette", API_RUNNER_SCRIPT)
        self.assertIn("pse_interaction_colors", _FP_SESSION_SCRIPT)
        self.assertIn("def _apply_pse_interaction_colors", _FP_SESSION_SCRIPT)
        self.assertIn("def _save_viewer_session_with_palette", _FP_SESSION_SCRIPT)
        self.assertIn("HIP2L_PSE_INTERACTION_COLORS", _PSE_FILTER_SCRIPT)
        self.assertIn('cmd.set("label_color", color_name, object_name)', _PSE_FILTER_SCRIPT)
        self.assertIn("def _save_viewer_session_with_palette", _PSE_FILTER_SCRIPT)

    def test_pse_palette_matches_exact_luna_interaction_object_segments(self) -> None:
        tree = ast.parse(API_RUNNER_SCRIPT)
        helper_names = {
            "_pymol_interaction_tokens",
            "_pymol_interaction_markers",
            "_matches_pymol_interaction_object",
        }
        helpers = [
            node for node in tree.body
            if isinstance(node, ast.FunctionDef) and node.name in helper_names
        ]
        namespace = {"re": re}
        exec(compile(ast.Module(body=helpers, type_ignores=[]), "<pse-palette>", "exec"), namespace)
        matches = namespace["_matches_pymol_interaction_object"]

        self.assertTrue(matches("entry.all_inters.inter.hbond.i0.line", "entryallintersinterhb...", "Hydrogen bond"))
        self.assertTrue(matches("entry.all_inters.inter.xbond.i0.line", "entryallintersinterxbondi0line", "Halogen bond"))
        self.assertTrue(matches("entry.all_inters.inter.disp_face-to-face_stack.i0.line", "entryallintersinterdispfacetofacestacki0line", "Displaced face-to-face pi-stacking"))
        self.assertFalse(matches("entry.all_inters.inter.face-to-face_stack.i0.line", "entryallintersinterfacetofacestacki0line", "Displaced face-to-face pi-stacking"))


if __name__ == "__main__":
    unittest.main()
