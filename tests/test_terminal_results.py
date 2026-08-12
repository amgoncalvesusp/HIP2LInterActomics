from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from PIL import Image

from luna_gui.core import results_analysis
from luna_gui.core.project import ProjectConfig
from luna_gui.core.plot_manifest import load_manifest
from luna_gui.core.terminal_results import (
    _annotate_bar_segments,
    _plot_worker_count,
    _save_complete_residue_heatmap,
    _save_residue_heatmaps,
    _save_similarity_figure,
    _save_stacked_interaction_distribution,
    run_terminal_results,
)


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
            plot_manifest = load_manifest(workdir)
            self.assertEqual(manifest["outputs"]["plot_count"], len(plot_manifest.records))
            for language in ("en", "pt", "es"):
                for profile in ("screen", "report"):
                    complete = (
                        results
                        / "plots"
                        / language
                        / profile
                        / "heatmaps"
                        / "complete_ligands_residues_heatmap.png"
                    )
                    self.assertTrue(complete.exists())
                    selected = plot_manifest.select(language=language, profile=profile)
                    self.assertTrue(any(row.plot_id == "complete_interaction_heatmap" for row in selected))
                    self.assertTrue(any(row.plot_id == "interactions_by_amino_acid" for row in selected))
                    self.assertTrue(
                        (results / "plots" / language / profile / "distribution" / "interactions_by_amino_acid.png").exists()
                    )

    def test_extreme_memory_pressure_forces_sequential_language_workers(self) -> None:
        with patch(
            "luna_gui.core.terminal_results._available_memory_bytes",
            return_value=2 * 1024 ** 3,
        ):
            self.assertEqual(_plot_worker_count(), 1)

    def test_complete_heatmap_export_keeps_each_interaction_color(self) -> None:
        artifact = {
            "entries": ["ligA"],
            "residues": ["A/TYR/2"],
            "interaction_types": ["Hydrogen bond", "Hydrophobic", "Ionic"],
            "matrix": {
                "Hydrogen bond": [[1]],
                "Hydrophobic": [[1]],
                "Ionic": [[1]],
            },
        }
        with tempfile.TemporaryDirectory() as tmp:
            output = _save_complete_residue_heatmap(artifact, Path(tmp))
            self.assertIsNotNone(output)
            with Image.open(output).convert("RGB") as image:
                pixels = set(image.getdata())
                self.assertEqual(image.getpixel((0, 0)), (255, 255, 255))

        for interaction_name in artifact["interaction_types"]:
            color = results_analysis.get_interaction_color(interaction_name).lstrip("#")
            expected = tuple(int(color[index : index + 2], 16) for index in range(0, 6, 2))
            self.assertIn(expected, pixels)

    def test_trajectory_similarity_plot_places_larger_frames_at_the_top(self) -> None:
        captured: dict[str, object] = {}

        def capture_plot(figure, output, _plt, **_kwargs):
            axis = figure.axes[0]
            captured["matrix"] = axis.images[0].get_array().tolist()
            captured["labels"] = [label.get_text() for label in axis.get_yticklabels()]
            _plt.close(figure)
            return Path(output)

        with tempfile.TemporaryDirectory() as tmp, patch(
            "luna_gui.core.terminal_results._save_plot",
            side_effect=capture_plot,
        ):
            _save_similarity_figure(
                ["frame0_LIG", "frame1000_LIG", "frame1001_LIG"],
                [[1.0, 0.1, 0.2], [0.1, 1.0, 0.3], [0.2, 0.3, 1.0]],
                Path(tmp),
                "EIFP",
                trajectory_analysis=True,
            )

        self.assertEqual(captured["matrix"], [[1.0, 0.3, 0.2], [0.3, 1.0, 0.1], [0.2, 0.1, 1.0]])
        self.assertEqual(captured["labels"], ["frame1001_LIG", "frame1000_LIG", "frame0000_LIG"])

    def test_trajectory_residue_heatmap_places_larger_frames_at_the_top(self) -> None:
        captured: dict[str, object] = {}

        def capture_plot(figure, output, _plt, **_kwargs):
            axis = figure.axes[0]
            captured["matrix"] = axis.images[0].get_array().tolist()
            captured["labels"] = [label.get_text() for label in axis.get_yticklabels()]
            _plt.close(figure)
            return Path(output)

        artifact = {
            "entries": ["frame0_LIG", "frame1000_LIG", "frame1001_LIG"],
            "residues": ["A/GLY/1"],
            "matrix": {"Hydrogen bond": [[0.0], [1000.0], [1001.0]]},
        }
        with tempfile.TemporaryDirectory() as tmp, patch(
            "luna_gui.core.terminal_results._save_plot",
            side_effect=capture_plot,
        ):
            _save_residue_heatmaps(artifact, Path(tmp), trajectory_analysis=True)

        self.assertEqual(captured["matrix"], [[1001.0], [1000.0], [0.0]])
        self.assertEqual(
            captured["labels"],
            ["Frame: frame1001_LIG", "Frame: frame1000_LIG", "Frame: frame0000_LIG"],
        )

    def test_bar_segment_annotations_start_at_five_percent(self) -> None:
        from luna_gui.core.terminal_results import _get_pyplot

        plt = _get_pyplot()
        if plt is None:
            self.skipTest("matplotlib is unavailable")
        figure, axis = plt.subplots()
        try:
            bars = axis.bar([0, 1], [5.0, 4.0], color="#2f7f83")
            _annotate_bar_segments(axis, bars, [5.0, 4.0], percentage=True)
            self.assertEqual([text.get_text() for text in axis.texts], ["5.0%"])
        finally:
            plt.close(figure)

    def test_count_segments_can_be_labeled_as_percentages_of_their_total(self) -> None:
        from luna_gui.core.terminal_results import _get_pyplot

        plt = _get_pyplot()
        if plt is None:
            self.skipTest("matplotlib is unavailable")
        figure, axis = plt.subplots()
        try:
            bars = axis.bar([0, 1], [10.0, 30.0], color="#2f7f83")
            _annotate_bar_segments(
                axis,
                bars,
                [10.0, 30.0],
                references=[40.0, 40.0],
                percentage=False,
                label_percentages=True,
            )
            self.assertEqual([text.get_text() for text in axis.texts], ["25%", "75%"])
        finally:
            plt.close(figure)

    def test_ligand_atom_distribution_embeds_the_atom_map_and_keeps_over_100_percent(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            output_dir = Path(tmp)
            atom_map_path = output_dir / "ligand_atom_map.png"
            Image.new("RGB", (12, 12), (255, 0, 0)).save(atom_map_path)
            output = _save_stacked_interaction_distribution(
                ["C1"],
                ["Hydrogen bond", "Hydrophobic"],
                values=[[100.0, 100.0]],
                output=output_dir / "interactions_by_ligand_atom.png",
                title="Interactions by ligand atoms",
                xlabel="Ligand atoms",
                ylabel="% of frames (entries)",
                language="en",
                profile=None,
                percentage=True,
                ligand_atom_map_path=atom_map_path,
            )
            self.assertIsNotNone(output)
            with Image.open(output).convert("RGB") as image:
                self.assertIn((255, 0, 0), set(image.getdata()))

    def test_trajectory_exports_ligand_atom_distribution_in_all_languages(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            workdir = Path(tmp)
            results = workdir / "results"
            results.mkdir()
            (results / "analysis_summary.json").write_text(
                json.dumps({"entries": 2, "interaction_counts": {"Hydrogen bond": 2}}),
                encoding="utf-8",
            )
            (results / "residue_matrix.json").write_text(
                json.dumps({
                    "entries": ["frame_000", "frame_001"],
                    "residues": ["A/GLY/1"],
                    "interaction_types": ["Hydrogen bond"],
                    "matrix": {"Hydrogen bond": [[1], [1]]},
                    "ligand_atoms": ["C1", "N2"],
                    "ligand_atom_matrix": {
                        "Hydrogen bond": [[1, 0], [0, 1]],
                    },
                }),
                encoding="utf-8",
            )
            cfg = ProjectConfig(workdir=str(workdir), trajectory_analysis=True)
            run_terminal_results(cfg, sys.executable, {"terminal_matrix_max_entries": 10})
            manifest = load_manifest(workdir)
            for language in ("en", "pt", "es"):
                selected = manifest.select(language=language, profile="report")
                self.assertTrue(any(row.plot_id == "interactions_by_ligand_atom" for row in selected))


if __name__ == "__main__":
    unittest.main()
