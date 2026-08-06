from __future__ import annotations

import json
import os
import subprocess
import sys
import tempfile
import threading
import unittest
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from unittest import mock

import hipplinteractomics_multiple_run as batch
import hipplinteractomics_terminal as terminal


ROOT = Path(__file__).resolve().parents[1]


class TerminalCliTests(unittest.TestCase):
    def test_import_does_not_load_qt(self) -> None:
        probe = (
            "import sys, hipplinteractomics_terminal; "
            "assert not any(name.startswith('PyQt') for name in sys.modules); "
            "assert __import__('os').environ['MPLBACKEND'] == 'Agg'"
        )
        completed = subprocess.run(
            [sys.executable, "-c", probe],
            cwd=ROOT,
            capture_output=True,
            text=True,
            check=False,
        )
        self.assertEqual(completed.returncode, 0, completed.stderr)

    def test_direct_arguments_override_nested_json(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            config_path = Path(temporary) / "project.json"
            config_path.write_text(
                json.dumps(
                    {
                        "project": {
                            "protein_file": "protein.pdb",
                            "ligand_file": "ligands",
                            "workdir": "old-workdir",
                            "ifp_length": 1024,
                        },
                        "terminal": {
                            "env_name": "old-env",
                            "fp_session": {"entry_name": "old-entry"},
                        },
                    }
                ),
                encoding="utf-8",
            )
            project, runtime = terminal._effective_settings(
                config_path,
                {"workdir": "new-workdir", "ifp_length": 2048},
                {
                    "env_name": "new-env",
                    "fp_session": {"feature_id": "42"},
                },
            )

        self.assertEqual(project["workdir"], "new-workdir")
        self.assertEqual(project["ifp_length"], 2048)
        self.assertEqual(runtime["env_name"], "new-env")
        self.assertEqual(runtime["fp_session"]["entry_name"], "old-entry")
        self.assertEqual(runtime["fp_session"]["feature_id"], "42")

    def test_parser_maps_format_and_repeated_values(self) -> None:
        parser = terminal._build_parser()
        args = parser.parse_args(
            [
                "--protein-file",
                "protein.pdb",
                "--ligand-file",
                "ligands",
                "--workdir",
                "output",
                "--ifp-format",
                "cnt",
                "--selected-ligand",
                "ligand-a",
                "--selected-ligand",
                "ligand-b",
            ]
        )
        project, runtime = terminal._collect_cli_overrides(args, parser)
        self.assertFalse(project["ifp_bit"])
        self.assertEqual(project["selected_ligands"], ["ligand-a", "ligand-b"])
        self.assertEqual(runtime, {})

    def test_parser_exposes_protein_heteroatom_option(self) -> None:
        parser = terminal._build_parser()
        args = parser.parse_args(["--include-protein-heteroatoms"])
        project, runtime = terminal._collect_cli_overrides(args, parser)

        self.assertTrue(project["include_protein_heteroatoms"])
        self.assertEqual(runtime, {})

    def test_template_is_valid_json(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "template.json"
            code = terminal.main(["--write-template", str(output)])
            payload = json.loads(output.read_text(encoding="utf-8"))

        self.assertEqual(code, 0)
        self.assertIn("protein_file", payload)
        self.assertIn("workdir", payload)

    def test_results_only_persists_the_workdir_configuration(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            workdir = root / "workdir"
            config_path = root / "project.json"
            config_path.write_text(
                json.dumps({"workdir": str(workdir)}),
                encoding="utf-8",
            )
            with mock.patch.object(
                terminal,
                "_resolve_luna_python",
                return_value=Path(sys.executable),
            ), mock.patch(
                "luna_gui.core.project.HISTORY_FILE",
                root / "history.json",
            ), mock.patch.object(
                terminal.terminal_results,
                "run_terminal_results",
                return_value={"outputs": {}, "errors": []},
            ):
                code = terminal.run_results_from_config(config_path)

            saved = workdir / ".luna_gui.json"
            self.assertTrue(saved.is_file())
            self.assertEqual(json.loads(saved.read_text(encoding="utf-8"))["workdir"], str(workdir))
        self.assertEqual(code, 1)

    def test_prepare_complexes_mode_processes_the_folder_and_exits(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "complexes"
            output = root / "prepared"
            source.mkdir()
            (source / "pose.mol2").write_text("@<TRIPOS>MOLECULE\n", encoding="utf-8")
            result = mock.Mock(
                files_processed=1,
                proteins_written=1,
                ligands_written=1,
                protein_dir=str(output / "proteinas_pdb"),
                ligand_dir=str(output / "ligantes_mol2"),
                log_file=str(output / "preprocess.log"),
                errors=[],
            )
            with mock.patch.object(
                terminal.mol2_prep,
                "split_complex_folder",
                return_value=result,
            ) as split:
                code = terminal.main(
                    [
                        "--prepare-complexes",
                        str(source),
                        "--prepare-output",
                        str(output),
                        "--last-protein-atom",
                        "42",
                    ]
                )

        self.assertEqual(code, 0)
        split.assert_called_once_with(
            source.resolve(),
            last_pa=42,
            out_folder=output.resolve(),
            chemistry_python=None,
        )


class MultipleRunTests(unittest.TestCase):
    def test_default_terminal_command_uses_the_shipped_terminal_script(self) -> None:
        with mock.patch.object(batch.shutil, "which", return_value=None):
            command = batch.resolve_terminal_command()

        self.assertEqual(command[0], sys.executable)
        self.assertEqual(Path(command[1]).name, "hipplinteractomics_terminal.py")

    def test_list_syntaxes_are_supported(self) -> None:
        self.assertEqual(batch.parse_bits(["[1024, 2048]"]), [1024, 2048])
        self.assertEqual(batch.parse_formats(["binary", "cnt"]), ["bin", "cnt"])
        self.assertEqual(
            batch.parse_levels_growth(["[(2,10), (3,5), (6,2)]"]),
            [(2, 10.0), (3, 5.0), (6, 2.0)],
        )
        self.assertEqual(
            batch.parse_levels_growth(["2:10", "3:5"]),
            [(2, 10.0), (3, 5.0)],
        )

    def test_cartesian_matrix_writes_unique_json_configs(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            generated = batch.generate_configurations(
                base_config={
                    "project": {
                        "protein_file": "protein.pdb",
                        "ligand_file": "ligands",
                        "workdir": str(root / "runs"),
                        "include_protein_heteroatoms": True,
                    },
                    "terminal": {"env_name": "luna-env"},
                },
                bits=[1024, 2048],
                formats=["bin", "cnt"],
                levels_growth=[(2, 10.0), (3, 5.0)],
                config_dir=root / "configs",
            )
            documents = [
                json.loads(item.config_path.read_text(encoding="utf-8"))
                for item in generated
            ]

        self.assertEqual(len(generated), 8)
        self.assertEqual(len({item.run_id for item in generated}), 8)
        self.assertTrue(documents[0]["project"]["ifp_bit"])
        self.assertFalse(documents[2]["project"]["ifp_bit"])
        self.assertEqual(documents[0]["project"]["ifp_levels"], 2)
        self.assertEqual(documents[0]["project"]["ifp_radius"], 10.0)
        self.assertTrue(documents[0]["project"]["include_protein_heteroatoms"])
        self.assertTrue(
            documents[0]["project"]["workdir"].endswith("B1024_L2_G10_bin")
        )

    def test_atomic_json_writes_do_not_share_temporary_names(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            destination = root / "summary.json"
            gate = threading.Barrier(12)

            def write(index: int) -> None:
                gate.wait()
                batch._write_json_atomic(destination, {"writer": index})

            with ThreadPoolExecutor(max_workers=12) as executor:
                list(executor.map(write, range(12)))

            payload = json.loads(destination.read_text(encoding="utf-8"))
            leftovers = list(root.glob("*.tmp"))

        self.assertIn(payload["writer"], range(12))
        self.assertEqual(leftovers, [])

    def test_static_ligands_are_materialized_once_for_all_configs(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            ligand_path = root / "ligands.sdf"
            ligand_path.write_text("ligand-a\n$$$$\n", encoding="utf-8")
            base = {
                "project": {
                    "ligand_file": str(ligand_path),
                    "selected_ligands": "ALL",
                    "workdir": str(root / "runs"),
                },
                "terminal": {},
            }
            with mock.patch.object(
                batch.ligand_io,
                "parse_ligand_file",
                return_value=["ligand-a", "ligand-b"],
            ) as parse:
                prepared = batch.prepare_static_inputs(base)
                generated = batch.generate_configurations(
                    base_config=prepared,
                    bits=[1024, 2048],
                    formats=["bin", "cnt"],
                    levels_growth=[(2, 10.0)],
                    config_dir=root / "configs",
                )

            documents = [
                json.loads(item.config_path.read_text(encoding="utf-8"))
                for item in generated
            ]

        parse.assert_called_once_with(str(ligand_path))
        self.assertTrue(
            all(
                document["project"]["selected_ligands"]
                == ["ligand-a", "ligand-b"]
                for document in documents
            )
        )

    def test_main_executes_all_generated_configs_serially(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            base_config = root / "base.json"
            base_config.write_text(
                json.dumps(
                    {
                        "protein_file": "protein.pdb",
                        "ligand_file": "ligands",
                        "workdir": str(root / "runs"),
                    }
                ),
                encoding="utf-8",
            )
            order_file = root / "order.txt"
            fake_terminal = root / "fake_terminal.py"
            fake_terminal.write_text(
                "import os, pathlib, sys\n"
                "config = pathlib.Path(sys.argv[1])\n"
                "with open(os.environ['ORDER_FILE'], 'a', encoding='utf-8') as out:\n"
                "    out.write(config.stem + '\\n')\n"
                "print('processed', config.name)\n",
                encoding="utf-8",
            )
            config_dir = root / "configs"
            with mock.patch.dict(
                os.environ,
                {"ORDER_FILE": str(order_file)},
                clear=False,
            ), mock.patch.object(
                batch,
                "_write_json_atomic",
                wraps=batch._write_json_atomic,
            ) as atomic_write:
                code = batch.main(
                    [
                        "--base-config",
                        str(base_config),
                        "--bits",
                        "1024",
                        "2048",
                        "--formats",
                        "bin",
                        "cnt",
                        "--levels-growth",
                        "2:10",
                        "3:5",
                        "--config-dir",
                        str(config_dir),
                        "--terminal-executable",
                        str(fake_terminal),
                    ]
                )

            observed = order_file.read_text(encoding="utf-8").splitlines()
            summary = json.loads(
                (config_dir / "pipeline_summary.json").read_text(encoding="utf-8")
            )
            logs = list((config_dir / "logs").glob("*.log"))
            summary_writes = [
                call
                for call in atomic_write.call_args_list
                if Path(call.args[0]).name == "pipeline_summary.json"
            ]

        self.assertEqual(code, 0)
        self.assertEqual(len(observed), 8)
        self.assertEqual(len(logs), 8)
        self.assertEqual(len(summary_writes), 10)
        self.assertEqual(summary["status"], "completed")
        self.assertEqual(
            observed[0],
            "hipplinteractomics_B1024_L2_G10_bin",
        )
        self.assertEqual(
            observed[-1],
            "hipplinteractomics_B2048_L3_G5_cnt",
        )

    def test_resume_skips_only_completed_matching_configurations(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            base_config = root / "base.json"
            base_config.write_text(
                json.dumps(
                    {
                        "protein_file": "protein.pdb",
                        "ligand_file": "ligands",
                        "selected_ligands": ["ligand-a"],
                        "workdir": str(root / "runs"),
                    }
                ),
                encoding="utf-8",
            )
            order_file = root / "order.txt"
            fake_terminal = root / "fake_terminal.py"
            fake_terminal.write_text(
                "import os, pathlib, sys\n"
                "config = pathlib.Path(sys.argv[1])\n"
                "with open(os.environ['ORDER_FILE'], 'a', encoding='utf-8') as out:\n"
                "    out.write(config.stem + '\\n')\n"
                "token = os.environ.get('FAIL_TOKEN')\n"
                "raise SystemExit(9 if token and token in config.stem else 0)\n",
                encoding="utf-8",
            )
            config_dir = root / "configs"
            arguments = [
                "--base-config",
                str(base_config),
                "--bits",
                "1024",
                "2048",
                "--formats",
                "bin",
                "--levels-growth",
                "2:10",
                "--config-dir",
                str(config_dir),
                "--terminal-executable",
                str(fake_terminal),
            ]

            with mock.patch.dict(
                os.environ,
                {"ORDER_FILE": str(order_file), "FAIL_TOKEN": "B2048"},
                clear=False,
            ):
                first_code = batch.main(arguments)
            first_summary = json.loads(
                (config_dir / "pipeline_summary.json").read_text(encoding="utf-8")
            )

            with mock.patch.dict(
                os.environ,
                {"ORDER_FILE": str(order_file), "FAIL_TOKEN": ""},
                clear=False,
            ):
                second_code = batch.main(arguments)

            observed = order_file.read_text(encoding="utf-8").splitlines()
            final_summary = json.loads(
                (config_dir / "pipeline_summary.json").read_text(encoding="utf-8")
            )

        self.assertEqual(first_code, 1)
        self.assertEqual(
            [run["status"] for run in first_summary["runs"]],
            ["completed", "failed"],
        )
        self.assertEqual(second_code, 0)
        self.assertEqual(
            observed.count("hipplinteractomics_B1024_L2_G10_bin"),
            1,
        )
        self.assertEqual(
            observed.count("hipplinteractomics_B2048_L2_G10_bin"),
            2,
        )
        self.assertEqual(final_summary["status"], "completed")
        self.assertTrue(
            all(run["status"] == "completed" for run in final_summary["runs"])
        )


if __name__ == "__main__":
    unittest.main()
