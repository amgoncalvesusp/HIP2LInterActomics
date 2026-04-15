from __future__ import annotations

from pathlib import Path
import unittest

import numpy as np

from luna_gui.core.results_analysis import (
    cluster_rows,
    cluster_similarity_matrix,
    export_cluster_assignments,
    load_similarity_matrix,
)


class ResultsAnalysisTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.tmp_root = Path("tests/.tmp")
        cls.tmp_root.mkdir(parents=True, exist_ok=True)

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


if __name__ == "__main__":
    unittest.main()
