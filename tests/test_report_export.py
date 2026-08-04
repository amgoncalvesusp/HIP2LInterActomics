from __future__ import annotations

import re
import tempfile
import unittest
from pathlib import Path

from PIL import Image

from luna_gui.core.project import ProjectConfig
from luna_gui.core.plot_manifest import PlotManifest, PlotRecord, manifest_path
from luna_gui.core.report_export import (
    PdfReportError,
    _fit_image_box,
    build_report,
    collect_result_images,
    save_pdf_report,
    save_pdf_report_isolated,
)


class ReportExportTests(unittest.TestCase):
    def setUp(self) -> None:
        self.tempdir = tempfile.TemporaryDirectory()
        self.root = Path(self.tempdir.name)
        self.results = self.root / "results" / "terminal"
        self.results.mkdir(parents=True)
        self.cfg = ProjectConfig(
            workdir=str(self.root),
            protein_file="protein.pdb",
            ligand_file="ligands.sdf",
            language="pt",
        )
        self.analysis = {
            "entries": 3,
            "interaction_counts": {"Hydrogen bond": 8, "Hydrophobic": 5},
            "residue_counts": {"A/ASP/42": 4},
        }

    def tearDown(self) -> None:
        self.tempdir.cleanup()

    @staticmethod
    def _image(path: Path, size: tuple[int, int], color: str) -> Path:
        path.parent.mkdir(parents=True, exist_ok=True)
        Image.new("RGB", size, color).save(path)
        return path

    def test_collects_all_supported_result_charts_and_skips_report_cache(self) -> None:
        heatmap = self._image(self.results / "residue_heatmap.png", (900, 300), "#174f4b")
        network = self._image(self.results / "networks" / "contact_network.jpg", (400, 700), "#c8693a")
        self._image(self.root / "_report_pdf_similarity.png", (100, 100), "white")
        (self.results / "table.csv").write_text("a,b\n1,2\n", encoding="utf-8")

        pages = collect_result_images(self.root)

        self.assertEqual([item[1] for item in pages], [network, heatmap])
        self.assertIn("Rede", pages[0][0])
        self.assertIn("Mapa de calor", pages[1][0])

    def test_html_embeds_all_result_charts_with_landscape_and_proportional_css(self) -> None:
        self._image(self.results / "interaction_summary.png", (800, 400), "#174f4b")
        self._image(self.results / "residue_map_Hbond.png", (1200, 500), "#c8693a")

        document = build_report(cfg=self.cfg, analysis=self.analysis)

        self.assertEqual(document.count("data:image/png;base64,"), 2)
        self.assertIn("@page{size:A4 landscape", document)
        self.assertIn("width:auto;height:auto;object-fit:contain", document)
        self.assertIn("frequência ou intensidade das interações", document)
        self.assertIn("abundância relativa dos tipos de contato", document)

    def test_reports_explain_fp_columns_and_include_both_top50_rankings(self) -> None:
        dashboard = {
            "ifp_type": "EIFP",
            "features": [{"feature_id": 10}],
            "important_features": [{"feature_id": 10}],
            "model_name": "GradientBoosting",
            "threshold_pct": 90.0,
            "top_features_by_model": {
                "extra_trees": [{
                    "rank": 1,
                    "feature_id": 10,
                    "assigned_level": "2",
                    "assigned_class": "Has noncovalent interactions with the protein",
                    "coverage_pct": 75.0,
                    "importance_score": 0.61,
                }],
                "gradient_boosting": [{
                    "rank": 1,
                    "feature_id": 20,
                    "assigned_level": "1",
                    "assigned_class": "Ligand's level 0 features only",
                    "coverage_pct": 50.0,
                    "importance_score": 0.49,
                }],
            },
        }

        document = build_report(
            cfg=self.cfg,
            analysis=self.analysis,
            fp_dashboards={"EIFP": dashboard},
        )

        self.assertIn("Guia das colunas de Análises FP", document)
        self.assertIn("Níveis colisão", document)
        self.assertIn("EIFP / Extra Trees", document)
        self.assertIn("EIFP / Gradient Boosting", document)
        self.assertIn("table-layout:fixed", document)
        self.assertIn("overflow-wrap:anywhere", document)

    def test_pdf_is_landscape_and_image_box_keeps_source_aspect_ratio(self) -> None:
        wide = self._image(self.results / "wide_heatmap.png", (1600, 400), "#174f4b")
        self._image(self.results / "vertical_distribution.png", (300, 1200), "#c8693a")
        output = self.root / "report.pdf"

        save_pdf_report(output, cfg=self.cfg, analysis=self.analysis, heatmap_png=wide)

        data = output.read_bytes()
        media_boxes = re.findall(rb"/MediaBox\s*\[\s*0\s+0\s+([0-9.]+)\s+([0-9.]+)\s*\]", data)
        self.assertTrue(media_boxes)
        self.assertTrue(all(float(width) > float(height) for width, height in media_boxes))
        left, bottom, width, height = _fit_image_box(1600, 400)
        del left, bottom
        physical_ratio = (width * 11.69) / (height * 8.27)
        self.assertAlmostEqual(physical_ratio, 4.0, places=6)

    def test_isolated_renderer_returns_a_clean_error_instead_of_killing_parent(self) -> None:
        blocked_parent = self.root / "blocked"
        blocked_parent.write_text("not a directory", encoding="utf-8")

        with self.assertRaisesRegex(PdfReportError, "Falha ao gerar o PDF"):
            save_pdf_report_isolated(
                blocked_parent / "report.pdf",
                cfg=self.cfg,
                analysis=self.analysis,
                timeout=30,
            )

    def test_isolated_renderer_creates_a_valid_pdf(self) -> None:
        self._image(self.results / "interaction_network.png", (900, 500), "#174f4b")
        output = self.root / "isolated-report.pdf"

        result = save_pdf_report_isolated(
            output,
            cfg=self.cfg,
            analysis=self.analysis,
            timeout=30,
        )

        self.assertEqual(result, output)
        self.assertTrue(output.read_bytes().startswith(b"%PDF"))

    def test_manifest_filters_exactly_one_language_and_profile(self) -> None:
        paths = {
            "pt_report": self._image(self.root / "results/plots/pt/report/a.png", (200, 100), "red"),
            "pt_screen": self._image(self.root / "results/plots/pt/screen/a.png", (200, 100), "green"),
            "en_report": self._image(self.root / "results/plots/en/report/a.png", (200, 100), "blue"),
        }
        manifest = PlotManifest()
        manifest.add(PlotRecord("distribution", "pt", "report", str(paths["pt_report"]), "PT report", "", 10, category="distribution"))
        manifest.add(PlotRecord("distribution", "pt", "screen", str(paths["pt_screen"]), "PT screen", "", 10, category="distribution"))
        manifest.add(PlotRecord("distribution", "en", "report", str(paths["en_report"]), "EN report", "", 10, category="distribution"))
        manifest.save(manifest_path(self.root))

        document = build_report(cfg=self.cfg, analysis=self.analysis)

        self.assertEqual(document.count("data:image/png;base64,"), 1)
        self.assertIn("PT report", document)
        self.assertNotIn("PT screen", document)
        self.assertNotIn("EN report", document)

    def test_report_uses_scientific_order_and_deterministic_top30(self) -> None:
        manifest = PlotManifest()
        specs = [
            ("clusters", "05 Clusters", 50, "clusters", "", ""),
            ("complete", "03 Complete", 30, "complete_heatmap", "", ""),
            ("distribution", "01 Distribution", 10, "distribution", "", ""),
            ("by_type", "02 By type", 20, "interaction_heatmap", "", ""),
            ("similarity", "04 Similarity", 40, "similarity", "", ""),
            ("et_plot", "ET chart", 104, "fingerprint", "EIFP", "extra_trees"),
            ("gb_plot", "GB chart", 144, "fingerprint", "EIFP", "gradient_boosting"),
        ]
        for plot_id, title, sequence, category, ifp_type, model in specs:
            path = self._image(
                self.root / f"results/plots/pt/report/{plot_id}.png",
                (200, 100),
                "#174f4b",
            )
            manifest.add(PlotRecord(
                plot_id,
                "pt",
                "report",
                str(path),
                title,
                "caption",
                sequence,
                ifp_type=ifp_type,
                model=model,
                category=category,
            ))
        manifest.save(manifest_path(self.root))
        ranking = [{
            "rank": 1,
            "feature_id": 1,
            "assigned_level": "1",
            "assigned_class": "class",
            "coverage_pct": 50.0,
            "importance_score": 0.5,
        }]
        dashboard = {
            "ifp_type": "EIFP",
            "features": [{"feature_id": 1}],
            "important_features": [{"feature_id": 1}],
            "top_features_by_model": {
                "extra_trees": ranking,
                "gradient_boosting": ranking,
            },
        }
        analysis = dict(self.analysis)
        analysis["residue_counts"] = {f"A/RES/{index:02d}": 1 for index in range(35)}

        document = build_report(
            cfg=self.cfg,
            analysis=analysis,
            clusters=[("ligA", 1)],
            fp_dashboards={"EIFP": dashboard},
        )

        ordered_tokens = [
            "01 Distribution",
            "02 By type",
            "03 Complete",
            "04 Similarity",
            "05 Clusters",
            "Atribuição de clusters",
            "Como interpretar as análises de fingerprints",
            "Top 50 features: EIFP / Extra Trees",
            "ET chart",
            "Top 50 features: EIFP / Gradient Boosting",
            "GB chart",
        ]
        positions = [document.index(token) for token in ordered_tokens]
        self.assertEqual(positions, sorted(positions))
        top_section = document.split("Top 30 resíduos com mais interações", 1)[1].split("</table>", 1)[0]
        self.assertEqual(top_section.count("<tr><td>"), 30)
        self.assertIn("A/RES/00", top_section)
        self.assertIn("A/RES/29", top_section)
        self.assertNotIn("A/RES/30", top_section)


if __name__ == "__main__":
    unittest.main()
