from __future__ import annotations

import ast
import inspect
import json
from pathlib import Path
import unittest
from unittest.mock import patch

import numpy as np

import luna_gui.core.results_analysis as results_analysis

from luna_gui.core.results_analysis import (
    CLASS_INTRAPROTEIN,
    CLASS_L0_LIGAND,
    CLASS_NONCOVALENT,
    CLASS_UNRELIABLE,
    CLASS_UNRELIABLE_BY_CLASS,
    _gumbel_tail_p_value,
    _resolve_external_label_value,
    _tanimoto_similarity,
    build_complete_heatmap,
    build_complete_heatmap_layers,
    build_fp_analysis_dashboard,
    build_ligand_atom_entry_counts,
    build_ligand_atom_frame_percentages,
    build_trajectory_entry_counts,
    build_trajectory_frame_percentages,
    cluster_rows,
    cluster_similarity_matrix,
    export_cluster_assignments,
    format_trajectory_entry_name,
    format_residue_label,
    get_interaction_color,
    is_pi_stacking_interaction,
    load_analysis_summary,
    load_external_fp_labels,
    load_fp_analysis_artifacts,
    load_ifp_sparse_matrix,
    load_residue_matrix_artifact,
    load_similarity_matrix,
    normalize_fp_breakdown,
    normalize_fp_class_name,
    trajectory_entry_order,
    trajectory_frame_count,
    trajectory_frame_number,
)


class ResultsAnalysisTests(unittest.TestCase):
    def test_fp_analysis_uses_one_canonical_implementation_per_helper(self) -> None:
        tree = ast.parse(inspect.getsource(results_analysis))
        for name in ("build_fp_analysis_dashboard", "_resolve_training_labels", "_compute_feature_importances"):
            definitions = [node for node in tree.body if isinstance(node, ast.FunctionDef) and node.name == name]
            self.assertEqual(len(definitions), 1, name)

    def test_trajectory_frame_number_uses_full_trailing_number(self) -> None:
        entries = [
            "gold_soln_pose_19",
            "gold_soln_pose_100",
            "gold_soln_pose_2",
            "frame_11",
            "frame_3",
        ]
        ordered = sorted(
            entries,
            key=lambda name: -(trajectory_frame_number(name) or -1),
        )

        self.assertEqual(
            ordered,
            [
                "gold_soln_pose_100",
                "gold_soln_pose_19",
                "frame_11",
                "frame_3",
                "gold_soln_pose_2",
            ],
        )

    def test_trajectory_entry_order_and_display_padding_keep_frame_zero_last(self) -> None:
        entries = ["frame0_LIG", "frame1000_LIG", "frame1001_LIG", "frame1_LIG"]
        order = trajectory_entry_order(entries)

        self.assertEqual(order, [2, 1, 3, 0])
        self.assertEqual(trajectory_frame_count(entries), 1002)
        self.assertEqual(
            [format_trajectory_entry_name(entries[index], 30000) for index in order],
            ["frame01001_LIG", "frame01000_LIG", "frame00001_LIG", "frame00000_LIG"],
        )

    @classmethod
    def setUpClass(cls) -> None:
        cls.tmp_root = Path("tests/.tmp")
        cls.tmp_root.mkdir(parents=True, exist_ok=True)

    def test_load_external_fp_labels_uses_maximum_for_duplicate_numeric_scores(self) -> None:
        path = self.tmp_root / "duplicate_scores.tsv"
        path.write_text(
            "ligand_id\tscore\n"
            "ligA\t11.2\n"
            "ligA\t14.7\n"
            "ligB\t3.0\n",
            encoding="utf-8",
        )

        labels = load_external_fp_labels(path, label_column="score", id_column="ligand_id")

        self.assertEqual(labels["ligA"], "14.7")
        self.assertEqual(labels["ligB"], "3.0")

    def test_load_similarity_matrix_normalizes_labels_and_symmetry(self) -> None:
        path = self.tmp_root / "sim_matrix.csv"
        path.write_text(
            ",ligA,ligB,ligC\n"
            "ligA,1,0.8,0.2\n"
            "ligB,0.7,1,0.3\n"
            "ligC,0.2,0.3,1\n",
            encoding="utf-8",
        )

        labels, matrix = load_similarity_matrix(path)

        self.assertEqual(labels, ["ligA", "ligB", "ligC"])
        self.assertEqual(matrix.shape, (3, 3))
        self.assertTrue(np.allclose(matrix, matrix.T))
        self.assertTrue(np.allclose(np.diag(matrix), [1.0, 1.0, 1.0]))
        self.assertAlmostEqual(matrix[0, 1], 0.75)

    def test_cluster_similarity_matrix_groups_close_pairs(self) -> None:
        labels = ["A", "B", "C", "D"]
        matrix = np.array(
            [
                [1.0, 0.91, 0.21, 0.15],
                [0.91, 1.0, 0.18, 0.10],
                [0.21, 0.18, 1.0, 0.87],
                [0.15, 0.10, 0.87, 1.0],
            ]
        )

        result = cluster_similarity_matrix(labels, matrix, method="average", n_clusters=2)
        groups = {}
        for label, cluster_id in zip(result.labels, result.cluster_ids):
            groups.setdefault(cluster_id, []).append(label)
        grouped = {tuple(sorted(v)) for v in groups.values()}

        self.assertEqual(result.n_clusters, 2)
        self.assertEqual(grouped, {("A", "B"), ("C", "D")})
        self.assertEqual(result.ordered_matrix.shape, (4, 4))

    def test_export_cluster_assignments_writes_leaf_order(self) -> None:
        labels = ["A", "B", "C"]
        matrix = np.array(
            [
                [1.0, 0.95, 0.30],
                [0.95, 1.0, 0.25],
                [0.30, 0.25, 1.0],
            ]
        )
        result = cluster_similarity_matrix(labels, matrix, n_clusters=2)

        out = self.tmp_root / "clusters.csv"
        export_cluster_assignments(out, result)
        content = out.read_text(encoding="utf-8").splitlines()

        self.assertEqual(content[0], "ligand_id,cluster_id,leaf_order")
        self.assertEqual(len(content), 4)
        self.assertEqual(len(cluster_rows(result)), 3)

    def test_load_similarity_matrix_accepts_edge_list_format(self) -> None:
        path = self.tmp_root / "sim_edges.csv"
        path.write_text(
            "entry1,entry2,similarity\n"
            "ligA,ligB,0.80\n"
            "ligA,ligC,0.25\n"
            "ligB,ligC,0.60\n",
            encoding="utf-8",
        )

        labels, matrix = load_similarity_matrix(path)

        self.assertEqual(labels, ["ligA", "ligB", "ligC"])
        self.assertTrue(np.allclose(np.diag(matrix), [1.0, 1.0, 1.0]))
        self.assertAlmostEqual(matrix[0, 1], 0.8)
        self.assertAlmostEqual(matrix[1, 2], 0.6)

    def test_load_analysis_summary_reads_cached_artifact(self) -> None:
        workdir = self.tmp_root / "summary_workdir"
        (workdir / "results").mkdir(parents=True, exist_ok=True)
        (workdir / "results" / "analysis_summary.json").write_text(
            '{"entries": 2, "interaction_counts": {"Hydrogen bond": 5}}',
            encoding="utf-8",
        )

        data = load_analysis_summary(workdir)

        self.assertEqual(data["entries"], 2)
        self.assertEqual(data["interaction_counts"]["Hydrogen bond"], 5)

    def test_load_residue_matrix_artifact_reads_cached_artifact(self) -> None:
        workdir = self.tmp_root / "matrix_workdir"
        (workdir / "results").mkdir(parents=True, exist_ok=True)
        (workdir / "results" / "residue_matrix.json").write_text(
            '{"interaction_types":["Hydrogen bond"],"residues":["A/GLY/1"],"entries":["ligA"],"matrix":{"Hydrogen bond":[[2.0]]}}',
            encoding="utf-8",
        )

        data = load_residue_matrix_artifact(workdir)

        self.assertEqual(data["entries"], ["ligA"])
        self.assertEqual(data["matrix"]["Hydrogen bond"], [[2.0]])

    def test_load_fp_analysis_artifacts_reads_all_ifp_types(self) -> None:
        workdir = self.tmp_root / "fp_workdir"
        fp_dir = workdir / "results" / "fingerprints"
        fp_dir.mkdir(parents=True, exist_ok=True)
        (fp_dir / "fp_analysis_E.json").write_text('{"ifp_type":"EIFP","features":[]}', encoding="utf-8")
        (fp_dir / "fp_analysis_H.json").write_text('{"ifp_type":"HIFP","features":[]}', encoding="utf-8")

        artifacts = load_fp_analysis_artifacts(workdir)

        self.assertEqual(sorted(artifacts), ["EIFP", "HIFP"])
        self.assertEqual(artifacts["EIFP"]["ifp_type"], "EIFP")

    def test_build_complete_heatmap_collapses_dominant_interaction_type(self) -> None:
        artifact = {
            "entries": ["ligA", "ligB"],
            "residues": ["A/GLY/1", "A/TYR/2"],
            "interaction_types": ["Hydrogen bond", "Hydrophobic"],
            "matrix": {
                "Hydrogen bond": [[2.0, 0.0], [0.0, 1.0]],
                "Hydrophobic": [[1.0, 3.0], [0.0, 0.0]],
            },
        }

        entries, residues, categorical, interaction_types = build_complete_heatmap(artifact)

        self.assertEqual(entries, ["ligA", "ligB"])
        self.assertEqual(residues, ["A/GLY/1", "A/TYR/2"])
        self.assertEqual(interaction_types, ["Hydrogen bond", "Hydrophobic"])
        self.assertEqual(categorical.tolist(), [[1, 2], [0, 1]])

    def test_build_complete_heatmap_layers_keeps_all_interactions_in_priority_order(self) -> None:
        artifact = {
            "entries": ["ligA"],
            "residues": ["A/GLY/1"],
            "interaction_types": ["Hydrophobic", "Hydrogen bond", "Ionic"],
            "matrix": {
                "Hydrophobic": [[1.0]],
                "Hydrogen bond": [[2.0]],
                "Ionic": [[1.0]],
            },
        }

        entries, residues, layered, interaction_types = build_complete_heatmap_layers(artifact)

        self.assertEqual(entries, ["ligA"])
        self.assertEqual(residues, ["A/GLY/1"])
        self.assertEqual(layered, [[["Ionic", "Hydrogen bond", "Hydrophobic"]]])
        self.assertEqual(interaction_types, ["Ionic", "Hydrogen bond", "Hydrophobic"])

    def test_build_trajectory_frame_percentages_and_entry_counts(self) -> None:
        artifact = {
            "entries": ["frame1", "frame2", "frame3", "frame4"],
            "residues": ["A/ASP/1", "A/TYR/2", "A/GLY/3"],
            "ligand_atoms": ["C1", "N2", "O3"],
            "interaction_types": ["Hydrogen bond", "Ionic"],
            "matrix": {
                "Hydrogen bond": [
                    [1, 0, 0],
                    [0, 2, 0],
                    [1, 0, 0],
                    [0, 0, 0],
                ],
                "Ionic": [
                    [0, 1, 0],
                    [0, 0, 0],
                    [1, 1, 0],
                    [0, 0, 0],
                ],
            },
            "ligand_atom_matrix": {
                "Hydrogen bond": [
                    [1, 0, 0],
                    [0, 1, 0],
                    [1, 0, 0],
                    [0, 0, 0],
                ],
                "Ionic": [
                    [0, 1, 0],
                    [0, 0, 0],
                    [1, 1, 0],
                    [0, 0, 0],
                ],
            },
        }

        residues, interaction_types, percentages = build_trajectory_frame_percentages(artifact)
        self.assertEqual(residues, ["A/ASP/1", "A/TYR/2"])
        self.assertEqual(interaction_types, ["Hydrogen bond", "Ionic"])
        self.assertEqual(percentages.tolist(), [[50.0, 25.0], [25.0, 50.0]])

        residues, interaction_types, counts = build_trajectory_entry_counts(artifact, "frame3")
        self.assertEqual(residues, ["A/ASP/1", "A/TYR/2"])
        self.assertEqual(interaction_types, ["Hydrogen bond", "Ionic"])
        self.assertEqual(counts.tolist(), [[1.0, 1.0], [0.0, 1.0]])

        atoms, interaction_types, percentages = build_ligand_atom_frame_percentages(artifact)
        self.assertEqual(atoms, ["C1", "N2"])
        self.assertEqual(interaction_types, ["Hydrogen bond", "Ionic"])
        self.assertEqual(percentages.tolist(), [[50.0, 25.0], [25.0, 50.0]])

        atoms, interaction_types, counts = build_ligand_atom_entry_counts(artifact, "frame3")
        self.assertEqual(atoms, ["C1", "N2"])
        self.assertEqual(interaction_types, ["Hydrogen bond", "Ionic"])
        self.assertEqual(counts.tolist(), [[1.0, 1.0], [0.0, 1.0]])

    def test_format_residue_label_converts_three_letter_code(self) -> None:
        self.assertEqual(format_residue_label("A/GLY/12"), "A:G12")
        self.assertEqual(format_residue_label("B/TYR/104A"), "B:Y104A")
        self.assertEqual(format_residue_label("A/ZN/228"), "A/ZN/228")

    def test_pi_stacking_interactions_use_magenta_pink_palette(self) -> None:
        names = [
            "Displaced face-to-face pi stacking",
            "Displaced face-to-edge pi stacking",
            "Displaced face-to-slope pi stacking",
            "Face-to-face pi stacking",
            "Face-to-edge pi stacking",
            "Face-to-slope pi stacking",
            "pi stacking",
            "T-shape",
        ]
        colors = [get_interaction_color(name) for name in names]

        self.assertTrue(all(is_pi_stacking_interaction(name) for name in names))
        self.assertEqual(len(set(colors)), len(colors))
        self.assertTrue(all(color.lower().startswith("#") for color in colors))

    def test_normalize_fp_class_name_accepts_legacy_and_curly_labels(self) -> None:
        self.assertEqual(
            normalize_fp_class_name("Ligand’s level 0 features"),
            CLASS_L0_LIGAND,
        )
        self.assertEqual(
            normalize_fp_class_name("Only intra-protein interactions"),
            CLASS_INTRAPROTEIN,
        )

    def test_normalize_fp_breakdown_merges_aliases(self) -> None:
        breakdown = normalize_fp_breakdown(
            {
                "Ligand’s level 0 features": 2,
                "Ligand's level 0 features only": 3,
                "Has non-covalent interactions with the protein": 1,
            }
        )

        self.assertEqual(breakdown[CLASS_L0_LIGAND], 5)
        self.assertEqual(breakdown[CLASS_NONCOVALENT], 1)

    def test_load_ifp_sparse_matrix_reads_sparse_counts(self) -> None:
        path = self.tmp_root / "ifp_sparse.csv"
        path.write_text(
            "ligand_id,on_bits,count\n"
            "ligA,\"10\t20\t20\",\"1\t2\t3\"\n"
            "ligB,\"10\t30\",\"5\t1\"\n",
            encoding="utf-8",
        )

        labels, feature_ids, matrix = load_ifp_sparse_matrix(path)

        self.assertEqual(labels, ["ligA", "ligB"])
        self.assertEqual(feature_ids, [10, 20, 30])
        self.assertEqual(matrix.tolist(), [[1.0, 5.0, 0.0], [5.0, 0.0, 1.0]])

    def test_tanimoto_similarity_does_not_require_blas_matmul(self) -> None:
        matrix = np.zeros((5, 11), dtype=float)
        matrix[0, [0, 1, 2]] = 1
        matrix[1, [1, 2, 3]] = 1
        matrix[2, [8]] = 1

        similarity = _tanimoto_similarity(matrix)

        self.assertEqual(similarity.shape, (5, 5))
        self.assertTrue(np.allclose(similarity, similarity.T))
        self.assertAlmostEqual(similarity[0, 1], 2 / 4)
        self.assertAlmostEqual(similarity[0, 2], 0.0)

    def test_load_external_fp_labels_reads_headered_csv(self) -> None:
        path = self.tmp_root / "fp_labels.csv"
        path.write_text(
            "ligand_id,label\n"
            "ligA,active\n"
            "ligB,inactive\n",
            encoding="utf-8",
        )

        labels = load_external_fp_labels(path)

        self.assertEqual(labels, {"ligA": "active", "ligB": "inactive"})

    def test_load_external_fp_labels_auto_detects_rotulo_column(self) -> None:
        path = self.tmp_root / "fp_labels_rotulo.csv"
        path.write_text(
            "ligand_id,rótulo\n"
            "ligA,ativo\n"
            "ligB,inativo\n",
            encoding="utf-8",
        )

        labels = load_external_fp_labels(path)

        self.assertEqual(labels, {"ligA": "ativo", "ligB": "inativo"})

    def test_load_external_fp_labels_accepts_explicit_label_column(self) -> None:
        path = self.tmp_root / "fp_labels_activity.csv"
        path.write_text(
            "ligand_id,activity\n"
            "ligA,active\n"
            "ligB,inactive\n",
            encoding="utf-8",
        )

        labels = load_external_fp_labels(path, label_column="activity")

        self.assertEqual(labels, {"ligA": "active", "ligB": "inactive"})

    def test_load_external_fp_labels_accepts_tsv(self) -> None:
        path = self.tmp_root / "fp_labels.tsv"
        path.write_text(
            "ligand_id\tactivity\n"
            "ligA\tactive\n"
            "ligB\tinactive\n",
            encoding="utf-8",
        )

        labels = load_external_fp_labels(path, label_column="activity")

        self.assertEqual(labels, {"ligA": "active", "ligB": "inactive"})

    def test_load_external_fp_labels_accepts_explicit_id_column(self) -> None:
        path = self.tmp_root / "fp_labels_idcol.tsv"
        path.write_text(
            "molecule_chembl_id\tpKi_value\n"
            "CHEMBL1\t7.1\n"
            "CHEMBL2\t8.4\n",
            encoding="utf-8",
        )

        labels = load_external_fp_labels(path, label_column="pKi_value", id_column="molecule_chembl_id")

        self.assertEqual(labels, {"CHEMBL1": "7.1", "CHEMBL2": "8.4"})

    def test_build_fp_analysis_dashboard_applies_threshold_and_importance(self) -> None:
        workdir = self.tmp_root / "fp_dashboard_workdir"
        fp_dir = workdir / "results" / "fingerprints"
        fp_dir.mkdir(parents=True, exist_ok=True)
        (fp_dir / "ifp_E.csv").write_text(
            "ligand_id,on_bits,count\n"
            "ligA,\"10\t20\t40\",\"1\t1\t1\"\n"
            "ligB,\"10\t20\t30\",\"1\t1\t1\"\n"
            "ligC,\"10\t30\",\"1\t1\"\n"
            "ligD,\"10\t30\t50\",\"1\t1\t1\"\n",
            encoding="utf-8",
        )
        (fp_dir / "seed_ifp_E_importance.txt").write_text("77\n", encoding="utf-8")
        artifact = {
            "ifp_type": "EIFP",
            "ifp_label": "Extended",
            "features": [
                {
                    "feature_id": 10,
                    "molecule_hits": 4,
                    "coverage_pct": 100.0,
                    "collision_hits": 1,
                    "total_count": 10,
                    "nature_breakdown": {
                        "Has noncovalent interactions with the protein": 1,
                        "Features with collision in the same complex": 1,
                    },
                    "shell_levels": ["0", "2"],
                    "shell_level_breakdown": {"0": 1, "2": 3},
                    "collision_shell_levels": ["0", "2"],
                    "collision_level_breakdown": {"0": 1, "2": 3},
                },
                {
                    "feature_id": 20,
                    "molecule_hits": 2,
                    "coverage_pct": 50.0,
                    "collision_hits": 0,
                    "total_count": 4,
                    "nature_breakdown": {"Ligand's level 0 features only": 4},
                },
                {
                    "feature_id": 30,
                    "molecule_hits": 3,
                    "coverage_pct": 50.0,
                    "collision_hits": 0,
                    "total_count": 4,
                    "nature_breakdown": {"Intraprotein interactions only": 4},
                },
                {
                    "feature_id": 40,
                    "molecule_hits": 1,
                    "coverage_pct": 25.0,
                    "collision_hits": 1,
                    "total_count": 6,
                    "nature_breakdown": {
                        "Ligand's level 0 features only": 3,
                        "Protein's level 0 features only": 3,
                    },
                },
                {
                    "feature_id": 50,
                    "molecule_hits": 1,
                    "coverage_pct": 25.0,
                    "collision_hits": 1,
                    "total_count": 2,
                    "nature_breakdown": {
                        "Ligand's level 0 features only": 1,
                        "Protein's level 0 features only": 1,
                    },
                },
                {
                    "feature_id": 60,
                    "molecule_hits": 0,
                    "coverage_pct": 0.0,
                    "collision_hits": 0,
                    "total_count": 0,
                    "nature_breakdown": {},
                },
            ],
        }

        dashboard = build_fp_analysis_dashboard(workdir, artifact)

        self.assertEqual(dashboard["threshold_source"], "zscore_gt_1")
        self.assertAlmostEqual(dashboard["threshold_pct"], 100.0)
        self.assertEqual(dashboard["class_assignment"]["collision_count"], 3)
        self.assertEqual(dashboard["total_molecules"], 4)
        self.assertEqual(dashboard["cluster_count"], 2)
        features = {row["feature_id"]: row for row in dashboard["features"]}
        self.assertEqual(features[10]["assigned_class"], CLASS_UNRELIABLE_BY_CLASS)
        self.assertEqual(features[10]["shell_levels"], ["0", "2"])
        self.assertEqual(features[10]["collision_shell_levels"], ["0", "2"])
        self.assertEqual(features[10]["shell_level_breakdown"], {"0": 1, "2": 3})
        self.assertEqual(features[10]["assigned_level_label"], CLASS_UNRELIABLE_BY_CLASS)
        self.assertEqual(features[20]["assigned_class"], CLASS_L0_LIGAND)
        self.assertEqual(features[30]["assigned_class"], CLASS_INTRAPROTEIN)
        self.assertEqual(features[40]["assigned_class"], CLASS_UNRELIABLE_BY_CLASS)
        self.assertEqual(features[50]["assigned_class"], CLASS_UNRELIABLE_BY_CLASS)
        self.assertGreaterEqual(features[60]["zscore"], 0.0)
        self.assertGreater(features[20]["importance_score"], 0.0)
        self.assertGreater(features[30]["importance_score"], 0.0)
        self.assertIn("importance_pvalue", features[20])
        self.assertIn("importance_zscore", features[20])
        self.assertAlmostEqual(
            features[20]["importance_pvalue"],
            _gumbel_tail_p_value(features[20]["importance_zscore"]),
        )
        self.assertEqual(dashboard["label_source"], "derived_clusters")
        self.assertIn(dashboard["model_name"], {"GradientBoosting", "ExtraTrees fallback"})
        self.assertEqual(dashboard["important_selection"], "pvalue_lt_0.01")
        self.assertEqual(dashboard["random_seed"], 77)
        self.assertEqual(
            set(dashboard["top_features_by_model"]),
            {"extra_trees", "gradient_boosting"},
        )
        self.assertTrue(dashboard["top_features_by_model"]["extra_trees"])
        self.assertTrue(dashboard["top_features_by_model"]["gradient_boosting"])
        self.assertLessEqual(len(dashboard["top_features_by_model"]["extra_trees"]), 50)

    def test_build_fp_analysis_dashboard_uses_external_labels_csv_when_available(self) -> None:
        workdir = self.tmp_root / "fp_dashboard_labels_workdir"
        fp_dir = workdir / "results" / "fingerprints"
        fp_dir.mkdir(parents=True, exist_ok=True)
        (fp_dir / "ifp_E.csv").write_text(
            "ligand_id,on_bits,count\n"
            "ligA,\"10\t20\",\"1\t1\"\n"
            "ligB,\"10\t20\",\"1\t1\"\n"
            "ligC,\"30\t40\",\"1\t1\"\n"
            "ligD,\"30\t40\",\"1\t1\"\n",
            encoding="utf-8",
        )
        labels_csv = workdir / "fp_labels.csv"
        labels_csv.write_text(
            "ligand_id,activity\n"
            "ligA,active\n"
            "ligB,active\n"
            "ligC,inactive\n"
            "ligD,inactive\n",
            encoding="utf-8",
        )
        artifact = {
            "ifp_type": "EIFP",
            "ifp_label": "Extended",
            "features": [
                {
                    "feature_id": 10,
                    "molecule_hits": 2,
                    "coverage_pct": 50.0,
                    "collision_hits": 0,
                    "total_count": 2,
                    "nature_breakdown": {"Ligand's level 0 features only": 4},
                },
                {
                    "feature_id": 20,
                    "molecule_hits": 2,
                    "coverage_pct": 50.0,
                    "collision_hits": 0,
                    "total_count": 2,
                    "nature_breakdown": {"Has noncovalent interactions with the protein": 4},
                },
                {
                    "feature_id": 30,
                    "molecule_hits": 2,
                    "coverage_pct": 50.0,
                    "collision_hits": 0,
                    "total_count": 2,
                    "nature_breakdown": {"Intraprotein interactions only": 4},
                },
                {
                    "feature_id": 40,
                    "molecule_hits": 2,
                    "coverage_pct": 50.0,
                    "collision_hits": 0,
                    "total_count": 2,
                    "nature_breakdown": {"Upper level with protein atomic information only": 4},
                },
            ],
        }

        dashboard = build_fp_analysis_dashboard(
            workdir,
            artifact,
            labels_csv=labels_csv,
            labels_column="activity",
        )

        self.assertEqual(dashboard["label_source"], "external_csv")
        self.assertEqual(dashboard["labels_csv"], str(labels_csv))
        self.assertEqual(dashboard["labels_column"], "activity")
        self.assertIn("rótulos externos", dashboard["model_note"])
        self.assertIn(dashboard["model_name"], {"GradientBoosting", "ExtraTrees fallback"})
        self.assertTrue(any("importance_score" in row for row in dashboard["features"]))

    def test_build_fp_analysis_dashboard_uses_otsu_and_reliable_features_when_no_zscore_gt_one(self) -> None:
        workdir = self.tmp_root / "fp_dashboard_fallback_rule_workdir"
        fp_dir = workdir / "results" / "fingerprints"
        fp_dir.mkdir(parents=True, exist_ok=True)

        rows = ["ligand_id,on_bits,count"]
        sep = "\t"
        for idx in range(1, 26):
            ligand = f"lig{idx:02d}"
            bits: list[str] = []
            counts: list[str] = []
            if idx == 1:
                bits.append("10"); counts.append("1")
            if idx in {1, 2}:
                bits.append("20"); counts.append("1")
            if idx in {3, 4}:
                bits.append("30"); counts.append("1")
            if idx == 5:
                bits.append("40"); counts.append("1")
            rows.append(f'{ligand},"{sep.join(bits)}","{sep.join(counts)}"')
        (fp_dir / "ifp_E.csv").write_text("\n".join(rows) + "\n", encoding="utf-8")

        artifact = {
            "ifp_type": "EIFP",
            "ifp_label": "Extended",
            "features": [
                {
                    "feature_id": 10,
                    "molecule_hits": 1,
                    "coverage_pct": 4.0,
                    "collision_hits": 0,
                    "total_count": 1,
                    "nature_breakdown": {"Ligand's level 0 features only": 1},
                },
                {
                    "feature_id": 20,
                    "molecule_hits": 2,
                    "coverage_pct": 8.0,
                    "collision_hits": 0,
                    "total_count": 2,
                    "nature_breakdown": {"Has noncovalent interactions with the protein": 2},
                },
                {
                    "feature_id": 30,
                    "molecule_hits": 2,
                    "coverage_pct": 8.0,
                    "collision_hits": 0,
                    "total_count": 2,
                    "nature_breakdown": {"Intraprotein interactions only": 2},
                },
                {
                    "feature_id": 40,
                    "molecule_hits": 1,
                    "coverage_pct": 4.0,
                    "collision_hits": 0,
                    "total_count": 2,
                    "nature_breakdown": {
                        "Ligand's level 0 features only": 1,
                        "Protein's level 0 features only": 1,
                    },
                },
            ],
        }

        dashboard = build_fp_analysis_dashboard(workdir, artifact)

        self.assertEqual(dashboard["threshold_source"], "otsu")
        self.assertIn("Otsu's Thresholding", dashboard["model_note"])
        self.assertIn("features confiaveis", dashboard["model_note"])
        self.assertEqual(dashboard["importance_eligible_count"], 3)
        features = {row["feature_id"]: row for row in dashboard["features"]}
        self.assertEqual(features[10]["assigned_class"], CLASS_L0_LIGAND)
        self.assertEqual(features[20]["assigned_class"], CLASS_NONCOVALENT)
        self.assertEqual(features[30]["assigned_class"], CLASS_INTRAPROTEIN)
        self.assertEqual(features[40]["assigned_class"], CLASS_UNRELIABLE_BY_CLASS)
        self.assertTrue(features[10]["reliable"])
        self.assertTrue(features[10]["importance_eligible"])
        self.assertTrue(features[20]["importance_eligible"])
        self.assertTrue(features[30]["importance_eligible"])
        self.assertFalse(features[40]["importance_eligible"])
        self.assertIn("importance_zscore", features[10])
        self.assertAlmostEqual(
            features[10]["importance_pvalue"],
            _gumbel_tail_p_value(features[10]["importance_zscore"]),
        )


    def test_build_fp_analysis_dashboard_matches_external_numeric_labels_by_entry_id(self) -> None:
        workdir = self.tmp_root / "fp_dashboard_regression_workdir"
        fp_dir = workdir / "results" / "fingerprints"
        fp_dir.mkdir(parents=True, exist_ok=True)
        (fp_dir / "ifp_E.csv").write_text(
            "ligand_id,on_bits,count\n"
            "PROT:CHEMBL1_LIG,\"10\t20\",\"2\t1\"\n"
            "PROT:CHEMBL2_LIG,\"10\t30\",\"1\t3\"\n"
            "PROT:CHEMBL3_LIG,\"20\t30\",\"4\t1\"\n"
            "PROT:CHEMBL4_LIG,\"30\t40\",\"2\t5\"\n",
            encoding="utf-8",
        )
        labels_csv = workdir / "labels.tsv"
        labels_csv.write_text(
            "molecule_chembl_id\tvalue\n"
            "CHEMBL1\t1.2\n"
            "CHEMBL2\t2.4\n"
            "CHEMBL3\t3.6\n"
            "CHEMBL4\t4.8\n"
            "CHEMBL999\t9.9\n",
            encoding="utf-8",
        )
        artifact = {
            "ifp_type": "EIFP",
            "ifp_label": "Extended",
            "features": [
                {"feature_id": 10, "molecule_hits": 2, "coverage_pct": 50.0, "collision_hits": 0, "total_count": 3, "nature_breakdown": {"Ligand's level 0 features only": 2}},
                {"feature_id": 20, "molecule_hits": 2, "coverage_pct": 50.0, "collision_hits": 0, "total_count": 5, "nature_breakdown": {"Has noncovalent interactions with the protein": 2}},
                {"feature_id": 30, "molecule_hits": 3, "coverage_pct": 75.0, "collision_hits": 0, "total_count": 6, "nature_breakdown": {"Intraprotein interactions only": 3}},
                {"feature_id": 40, "molecule_hits": 1, "coverage_pct": 25.0, "collision_hits": 0, "total_count": 5, "nature_breakdown": {"Upper level with protein atomic information only": 1}},
            ],
        }

        dashboard = build_fp_analysis_dashboard(
            workdir,
            artifact,
            labels_csv=labels_csv,
            labels_id_column="molecule_chembl_id",
            labels_column="value",
        )

        self.assertEqual(dashboard["label_source"], "external_csv")
        self.assertEqual(dashboard["label_kind"], "regression")
        self.assertEqual(dashboard["matched_molecules"], 4)
        self.assertEqual(dashboard["labels_id_column"], "molecule_chembl_id")
        self.assertIn(dashboard["model_name"], {"GradientBoosting", "ExtraTrees fallback"})
        features = {row["feature_id"]: row for row in dashboard["features"]}
        self.assertTrue(all("importance_pvalue" in row for row in features.values()))

    def test_interaction_prevalence_uses_only_noncovalent_features_and_excludes_weak_hbond_hydrophobic(self) -> None:
        workdir = self.tmp_root / "fp_dashboard_interaction_prevalence_workdir"
        fp_dir = workdir / "results" / "fingerprints"
        fp_dir.mkdir(parents=True, exist_ok=True)
        bit_ids = [str(bit_id) for bit_id in range(10, 260, 10)]
        rows = ["ligand_id,on_bits,count"]
        sep = "\t"
        bits_text = sep.join(bit_ids)
        counts_text = sep.join(["1"] * len(bit_ids))
        for ligand in ("ligA", "ligB", "ligC", "ligD"):
            rows.append(f'{ligand},"{bits_text}","{counts_text}"')
        (fp_dir / "ifp_E.csv").write_text("\n".join(rows) + "\n", encoding="utf-8")
        labels_csv = workdir / "labels.csv"
        labels_csv.write_text(
            "ligand_id,label\n"
            "ligA,active\n"
            "ligB,active\n"
            "ligC,inactive\n"
            "ligD,inactive\n",
            encoding="utf-8",
        )
        (fp_dir / "fp_detail_E.json").write_text(
            json.dumps(
                {
                    "ifp_type": "EIFP",
                    "feature_details": {
                        "10": {
                            "interaction_counts": {
                                "Hydrophobic": 20,
                                "Weak hydrogen bond": 10,
                                "Ionic": 5,
                            },
                            "residue_counts": {
                                "A/LEU/2": 20,
                                "A/ASN/3": 10,
                                "A/ASP/1": 5,
                            },
                            "pair_counts": {
                                "Hydrophobic||A/LEU/2": 20,
                                "Weak hydrogen bond||A/ASN/3": 10,
                                "Ionic||A/ASP/1": 5,
                            },
                            "entries": {
                                "ligA": {
                                    "shell_count": 4,
                                    "interaction_counts": {"Ionic": 5, "Hydrophobic": 4},
                                    "residue_counts": {"A/ASP/1": 5, "A/LEU/2": 4},
                                    "pair_counts": {"Ionic||A/ASP/1": 5, "Hydrophobic||A/LEU/2": 4},
                                },
                                "ligB": {
                                    "shell_count": 4,
                                    "interaction_counts": {"Hydrophobic": 20, "Weak hydrogen bond": 10},
                                    "residue_counts": {"A/LEU/2": 20, "A/ASN/3": 10},
                                    "pair_counts": {
                                        "Hydrophobic||A/LEU/2": 20,
                                        "Weak hydrogen bond||A/ASN/3": 10,
                                    },
                                },
                            },
                        },
                        "20": {
                            "interaction_counts": {"Ionic": 12},
                            "residue_counts": {"A/ASP/9": 12},
                            "pair_counts": {"Ionic||A/ASP/9": 12},
                            "entries": {"ligC": {"shell_count": 1, "interaction_counts": {"Ionic": 12}, "residue_counts": {"A/ASP/9": 12}, "pair_counts": {"Ionic||A/ASP/9": 12}}},
                        },
                    },
                },
                ensure_ascii=False,
            ),
            encoding="utf-8",
        )
        artifact = {
            "ifp_type": "EIFP",
            "ifp_label": "Extended",
            "features": [
                {
                    "feature_id": 10,
                    "molecule_hits": 4,
                    "coverage_pct": 100.0,
                    "collision_hits": 0,
                    "total_count": 4,
                    "nature_breakdown": {"Has noncovalent interactions with the protein": 4},
                },
                {
                    "feature_id": 20,
                    "molecule_hits": 4,
                    "coverage_pct": 100.0,
                    "collision_hits": 0,
                    "total_count": 4,
                    "nature_breakdown": {"Ligand's level 0 features only": 4},
                },
            ]
            + [
                {
                    "feature_id": bit_id,
                    "molecule_hits": 4,
                    "coverage_pct": 100.0,
                    "collision_hits": 0,
                    "total_count": 4,
                    "nature_breakdown": {"Protein's level 0 features only": 4},
                }
                for bit_id in range(30, 260, 10)
            ],
        }

        scores = np.zeros(len(bit_ids), dtype=float)
        scores[0] = 1000.0
        with patch(
            "luna_gui.core.results_analysis._compute_feature_importances",
            return_value=(scores, "GradientBoosting", "mocked"),
        ):
            dashboard = build_fp_analysis_dashboard(
                workdir,
                artifact,
                labels_csv=labels_csv,
                labels_column="label",
                task_kind_preference="classification",
            )

        important = {row["feature_id"]: row for row in dashboard["important_features"]}
        self.assertIn(10, important)
        self.assertNotIn(20, [row["feature_id"] for row in dashboard["important_features"] if row.get("prevalent_interaction")])
        feature = important[10]
        self.assertEqual(feature["assigned_class"], CLASS_NONCOVALENT)
        self.assertEqual(feature["interaction_breakdown"], {"Ionic": 5})
        self.assertEqual(feature["residue_breakdown"], {"A/ASP/1": 5})
        self.assertEqual(feature["pair_breakdown"], {"Ionic||A/ASP/1": 5})
        self.assertEqual(feature["prevalent_interaction"], "Ionic")
        self.assertEqual(feature["prevalent_interaction_pct"], 100.0)
        self.assertEqual(feature["prevalent_interaction_entries"], ["ligA"])

    def test_interaction_prevalence_thresholds_exact_interaction_residue_pair(self) -> None:
        workdir = self.tmp_root / "fp_dashboard_pair_prevalence_workdir"
        fp_dir = workdir / "results" / "fingerprints"
        fp_dir.mkdir(parents=True, exist_ok=True)
        bit_ids = [str(bit_id) for bit_id in range(10, 260, 10)]
        bits_text = "\t".join(bit_ids)
        counts_text = "\t".join(["1"] * len(bit_ids))
        (fp_dir / "ifp_E.csv").write_text(
            "ligand_id,on_bits,count\n"
            f'ligA,"{bits_text}","{counts_text}"\n'
            f'ligB,"{bits_text}","{counts_text}"\n'
            f'ligC,"{bits_text}","{counts_text}"\n'
            f'ligD,"{bits_text}","{counts_text}"\n',
            encoding="utf-8",
        )
        labels_csv = workdir / "labels.csv"
        labels_csv.write_text(
            "ligand_id,label\nligA,active\nligB,active\nligC,inactive\nligD,inactive\n",
            encoding="utf-8",
        )
        (fp_dir / "fp_detail_E.json").write_text(
            json.dumps(
                {
                    "ifp_type": "EIFP",
                    "feature_details": {
                        "10": {
                            "interaction_counts": {"Hydrogen bond": 9, "Ionic": 6},
                            "residue_counts": {"A/ASP/1": 10, "A/TYR/2": 5},
                            "pair_counts": {
                                "Ionic||A/ASP/1": 6,
                                "Hydrogen bond||A/TYR/2": 5,
                                "Hydrogen bond||A/ASP/1": 4,
                            },
                            "entries": {
                                "ligA": {
                                    "shell_count": 1,
                                    "interaction_counts": {"Ionic": 6, "Hydrogen bond": 4},
                                    "residue_counts": {"A/ASP/1": 10},
                                    "pair_counts": {"Ionic||A/ASP/1": 6, "Hydrogen bond||A/ASP/1": 4},
                                },
                                "ligB": {
                                    "shell_count": 1,
                                    "interaction_counts": {"Hydrogen bond": 5},
                                    "residue_counts": {"A/TYR/2": 5},
                                    "pair_counts": {"Hydrogen bond||A/TYR/2": 5},
                                },
                            },
                        },
                    },
                },
                ensure_ascii=False,
            ),
            encoding="utf-8",
        )
        artifact = {
            "ifp_type": "EIFP",
            "ifp_label": "Extended",
            "features": [
                {
                    "feature_id": 10,
                    "molecule_hits": 4,
                    "coverage_pct": 100.0,
                    "collision_hits": 0,
                    "total_count": 4,
                    "nature_breakdown": {"Has noncovalent interactions with the protein": 4},
                },
            ]
            + [
                {
                    "feature_id": bit_id,
                    "molecule_hits": 4,
                    "coverage_pct": 100.0,
                    "collision_hits": 0,
                    "total_count": 4,
                    "nature_breakdown": {"Protein's level 0 features only": 4},
                }
                for bit_id in range(20, 260, 10)
            ],
        }
        scores = np.zeros(len(bit_ids), dtype=float)
        scores[0] = 1000.0
        with patch(
            "luna_gui.core.results_analysis._compute_feature_importances",
            return_value=(scores, "GradientBoosting", "mocked"),
        ):
            dashboard = build_fp_analysis_dashboard(
                workdir,
                artifact,
                labels_csv=labels_csv,
                labels_column="label",
                task_kind_preference="classification",
                use_otsu_threshold=True,
            )

        feature = {row["feature_id"]: row for row in dashboard["important_features"]}[10]
        self.assertEqual(feature["prevalent_pair"], "Ionic||A/ASP/1")
        self.assertEqual(feature["prevalent_interaction"], "Ionic")
        self.assertEqual(feature["prevalent_residue"], "A/ASP/1")
        self.assertEqual(feature["prevalent_pair_entries"], ["ligA"])

    def test_build_fp_analysis_dashboard_uses_eligible_features_when_model_is_unavailable(self) -> None:
        workdir = self.tmp_root / "fp_dashboard_model_unavailable_workdir"
        fp_dir = workdir / "results" / "fingerprints"
        fp_dir.mkdir(parents=True, exist_ok=True)
        (fp_dir / "ifp_E.csv").write_text(
            "ligand_id,on_bits,count\n"
            "ligA,\"10\t20\",\"1\t1\"\n"
            "ligB,\"10\",\"1\"\n"
            "ligC,\"20\",\"1\"\n"
            "ligD,\"30\",\"1\"\n",
            encoding="utf-8",
        )
        artifact = {
            "ifp_type": "EIFP",
            "ifp_label": "Extended",
            "features": [
                {
                    "feature_id": 10,
                    "molecule_hits": 2,
                    "coverage_pct": 50.0,
                    "collision_hits": 0,
                    "total_count": 2,
                    "nature_breakdown": {"Ligand's level 0 features only": 2},
                },
                {
                    "feature_id": 20,
                    "molecule_hits": 2,
                    "coverage_pct": 50.0,
                    "collision_hits": 0,
                    "total_count": 2,
                    "nature_breakdown": {"Has noncovalent interactions with the protein": 2},
                },
                {
                    "feature_id": 30,
                    "molecule_hits": 1,
                    "coverage_pct": 25.0,
                    "collision_hits": 0,
                    "total_count": 1,
                    "nature_breakdown": {"Protein's level 0 features only": 1},
                },
            ],
        }

        with patch(
            "luna_gui.core.results_analysis._compute_feature_importances",
            return_value=(np.zeros(3, dtype=float), "Unavailable", "ModuleNotFoundError"),
        ):
            dashboard = build_fp_analysis_dashboard(workdir, artifact)

        self.assertEqual(dashboard["important_selection"], "pvalue_lt_0.01")
        self.assertEqual(dashboard["important_features"], [])

    def test_external_label_matching_prefers_ligand_id_over_protein_prefix(self) -> None:
        label_map = {
            "CHEMBL112640": "protein_like_value",
            "CHEMBL113076": "ligand_value",
        }
        self.assertEqual(
            _resolve_external_label_value(label_map, "CHEMBL112640:CHEMBL113076_LIG"),
            "ligand_value",
        )

if __name__ == "__main__":
    unittest.main()
