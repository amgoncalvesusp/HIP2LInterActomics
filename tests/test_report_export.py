from __future__ import annotations

import json
import re
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from PIL import Image

import luna_gui.core.report_export as report_export
from luna_gui.core.project import ProjectConfig
from luna_gui.core.plot_manifest import PlotManifest, PlotRecord, manifest_path
from luna_gui.core.report_export import (
    PdfReportError,
    _fit_image_box,
    _pdf_payload,
    build_report,
    cleanup_pdf_render_job,
    collect_result_images,
    create_pdf_render_job,
    execute_pdf_render_job,
    read_pdf_render_status,
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

        pages = collect_result_images(self.root, language="pt")

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

    def test_html_and_pdf_payload_follow_the_configured_export_language(self) -> None:
        expectations = {
            "en": ("HIP2LInterActomics report", "Interaction count by type", "Hydrogen bond", "levels="),
            "pt": ("Relat\u00f3rio HIP2LInterActomics", "Contagem por tipo de intera\u00e7\u00e3o", "Liga\u00e7\u00e3o de hidrog\u00eanio", "n\u00edveis="),
            "es": ("Reporte HIP2LInterActomics", "Conteo por tipo de interacci\u00f3n", "Puente de hidr\u00f3geno", "niveles="),
        }
        for language, (report_title, count_title, interaction_label, levels_label) in expectations.items():
            with self.subTest(language=language):
                cfg = ProjectConfig(
                    workdir=str(self.root),
                    protein_file="protein.pdb",
                    ligand_file="ligands.sdf",
                    language=language,
                )
                document = build_report(cfg=cfg, analysis=self.analysis)
                payload = _pdf_payload(cfg, self.analysis, None, None, None, None, None, None)
                captured_titles: list[str] = []
                output = self.root / f"{language}-localized.pdf"
                with patch.object(
                    report_export,
                    "_reportlab_write_text_page",
                    side_effect=lambda _canvas, title, *_args: captured_titles.append(title),
                ):
                    report_export._write_pdf_payload_reportlab(output, payload)

                self.assertIn(report_title, document)
                self.assertIn(count_title, document)
                self.assertIn(interaction_label, document)
                self.assertEqual(payload["interaction_rows"][0][0], interaction_label)
                self.assertIn(levels_label, next(value for key, value in payload["cfg_rows"] if key == "IFP"))
                self.assertTrue(any(report_title in title for title in captured_titles))
                self.assertTrue(any(count_title in title for title in captured_titles))

    def test_html_embeds_every_legacy_per_type_heatmap_passed_by_gui(self) -> None:
        heatmaps = [
            self._image(
                self.root / "results" / "plots" / "heatmaps" / name,
                (900, 500),
                "#174f4b",
            )
            for name in (
                "interaction_hydrogen_bond.png",
                "interaction_cation_pi.png",
                "interaction_hydrophobic.png",
            )
        ]
        extra = [
            (
                f"Heatmap por tipo: {path.stem}",
                path,
                "Explicação científica do tipo de interação.",
            )
            for path in heatmaps
        ]

        document = build_report(cfg=self.cfg, analysis=self.analysis, extra_images=extra)

        self.assertEqual(document.count("data:image/png;base64,"), len(heatmaps))
        for path in heatmaps:
            self.assertIn(path.stem, document)

    def test_legacy_english_complete_heatmap_stays_in_the_complete_section(self) -> None:
        self.cfg.language = "en"
        distribution = self._image(self.root / "interaction_distribution.png", (800, 400), "#174f4b")
        complete = self._image(self.root / "complete_interaction_heatmap.png", (800, 400), "#c8693a")

        document = build_report(
            cfg=self.cfg,
            analysis=self.analysis,
            extra_images=[
                ("Interaction distribution", distribution, "distribution"),
                ("Complete interaction heatmap", complete, "complete"),
            ],
        )

        self.assertLess(
            document.index("Interaction distribution"),
            document.index("Complete interaction heatmap"),
        )

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

    def test_external_worker_job_uses_atomic_status_and_valid_pdf(self) -> None:
        output = self.root / "worker-report.pdf"
        job_path, status_path = create_pdf_render_job(
            output,
            cfg=self.cfg,
            analysis=self.analysis,
        )
        try:
            exit_code = execute_pdf_render_job(job_path, status_path)
            status = read_pdf_render_status(status_path)

            self.assertEqual(exit_code, 0)
            self.assertIsNotNone(status)
            self.assertTrue(status["ok"])
            self.assertTrue(output.read_bytes().startswith(b"%PDF"))
            self.assertFalse(output.with_name(f".{output.name}.part").exists())
        finally:
            cleanup_pdf_render_job(job_path)

    def test_external_worker_job_reports_failure_without_raising_in_parent(self) -> None:
        blocked_parent = self.root / "blocked-worker"
        blocked_parent.write_text("not a directory", encoding="utf-8")
        job_path, status_path = create_pdf_render_job(
            blocked_parent / "report.pdf",
            cfg=self.cfg,
            analysis=self.analysis,
        )
        try:
            self.assertEqual(execute_pdf_render_job(job_path, status_path), 1)
            status = read_pdf_render_status(status_path)
            self.assertIsNotNone(status)
            self.assertFalse(status["ok"])
        finally:
            cleanup_pdf_render_job(job_path)

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

    def test_report_includes_every_per_type_heatmap_from_active_language(self) -> None:
        manifest = PlotManifest()
        for index, interaction_type in enumerate(("Hydrogen bond", "Hydrophobic", "Ionic")):
            path = self._image(
                self.root / f"results/plots/pt/report/heatmaps/residue_map_{index}.png",
                (600, 900),
                "#174f4b",
            )
            manifest.add(PlotRecord(
                f"interaction_heatmap_{index}",
                "pt",
                "report",
                str(path),
                f"Heatmap por tipo: {interaction_type}",
                "Explicação científica do padrão de frequência e permanência.",
                20 + index,
                category="interaction_heatmap",
            ))
        manifest.save(manifest_path(self.root))

        document = build_report(cfg=self.cfg, analysis=self.analysis)

        self.assertEqual(document.count("data:image/png;base64,"), 3)
        for interaction_type in ("Hydrogen bond", "Hydrophobic", "Ionic"):
            self.assertIn(f"Heatmap por tipo: {interaction_type}", document)
        self.assertEqual(document.count("O mapa de calor resume"), 3)

    def test_report_uses_scientific_order_and_deterministic_top30(self) -> None:
        manifest = PlotManifest()
        specs = [
            ("clusters", "05 Clusters", 50, "clusters", "", ""),
            ("complete", "03 Complete", 30, "complete_heatmap", "", ""),
            ("distribution", "01 Distribution", 10, "distribution", "", ""),
            ("amino_acids", "01.1 Amino acid distribution", 11, "distribution", "", ""),
            ("ligand_atoms", "01.2 Ligand atom distribution", 12, "distribution", "", ""),
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
            "01.1 Amino acid distribution",
            "01.2 Ligand atom distribution",
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

    def test_report_outline_scans_persisted_fp_assets_and_is_shared_by_pdf(self) -> None:
        manifest = PlotManifest()
        specs = [
            ("interaction_distribution", "interaction summary", 10, "distribution", "", ""),
            ("interactions_by_amino_acid", "amino acids", 11, "distribution", "", ""),
            ("interactions_by_ligand_atom", "ligand atoms", 12, "distribution", "", ""),
            ("interaction_heatmap_hbond", "hbond heatmap", 20, "interaction_heatmap", "", ""),
            ("interaction_heatmap_ionic", "ionic heatmap", 21, "interaction_heatmap", "", ""),
            ("complete_interaction_heatmap", "complete heatmap", 30, "complete_heatmap", "", ""),
            ("similarity_matrix", "similarity", 40, "similarity", "EIFP", ""),
            ("hierarchical_clustering", "hierarchical", 50, "clusters", "EIFP", ""),
            ("reordered_matrix_cluster", "reordered", 60, "clusters", "EIFP", ""),
        ]
        fp_plot_ids = [
            "fp_top50_extra_trees",
            "fp_class_summary",
            "fp_class_assignment",
            "fp_coverage_importance",
            "fp_feature_presence_heatmap",
            "fp_interaction_assignment",
            "fp_prevalent_interactions",
            "fp_prevalent_interactions_heatmap",
        ]
        for index, plot_id in enumerate(fp_plot_ids):
            specs.append((plot_id, plot_id, 100 + index, "fingerprint", "EIFP", "extra_trees"))
        for plot_id, title, sequence, category, ifp_type, model in specs:
            path = self._image(
                self.root / "results" / "plots" / "pt" / "report" / category / f"{plot_id}.png",
                (300, 160),
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

        dashboard_dir = self.root / "results" / "terminal" / "fingerprints" / "E"
        dashboard_dir.mkdir(parents=True)
        dashboard_dir.joinpath("fp_dashboard.json").write_text(json.dumps({
            "ifp_type": "EIFP",
            "features": [{"feature_id": 10}],
            "important_features": [{"feature_id": 10}],
            "top_features_by_model": {
                "extra_trees": [{
                    "rank": 1,
                    "feature_id": 10,
                    "assigned_level": "2",
                    "assigned_class": "class",
                    "coverage_pct": 75.0,
                    "importance_score": 0.5,
                }],
            },
        }), encoding="utf-8")
        (self.root / "results" / "terminal" / "clusters_E.csv").write_text(
            "ligand_id,cluster_id,leaf_order\nligA,1,0\n",
            encoding="utf-8",
        )

        payload = _pdf_payload(self.cfg, self.analysis, None, None, None, None, None, None)
        section_titles = [
            report_export._numbered_title(section["number"], section["title"])
            if section["kind"] != "image" else section["image"]["title"]
            for section in payload["report_sections"]
        ]
        expected = [
            "1. Resumo",
            "2. Configura\u00e7\u00e3o",
            "3. Contagem por tipo de intera\u00e7\u00e3o",
            "4. Top 30 res\u00edduos com mais intera\u00e7\u00f5es",
            "5. Distribui\u00e7\u00e3o de intera\u00e7\u00f5es",
            "6. Mapas de calor de intera\u00e7\u00f5es",
            "7. An\u00e1lises de fingerprints",
            "7.4.1. EIFP: Atribui\u00e7\u00e3o de clusters",
            "8. Breve aprendizado supervisionado",
            "8.2.1.1.1. Top 50 features: EIFP / Extra Trees",
        ]
        positions = [next(index for index, title in enumerate(section_titles) if expected_title in title) for expected_title in expected]
        self.assertEqual(positions, sorted(positions))
        self.assertTrue(any(title.startswith("8.2.1.1.9.") for title in section_titles))

        document = build_report(cfg=self.cfg, analysis=self.analysis)
        for expected_title in expected:
            self.assertIn(expected_title, document)

        captured: list[str] = []
        output = self.root / "outline.pdf"
        with patch.object(
            report_export,
            "_reportlab_write_text_page",
            side_effect=lambda _canvas, title, *_args: captured.append(title),
        ), patch.object(
            report_export,
            "_reportlab_write_image_page",
            side_effect=lambda _canvas, image, *_args: captured.append(image["title"]),
        ):
            report_export._write_pdf_payload_reportlab(output, payload)

        self.assertEqual(captured[1:], section_titles)


if __name__ == "__main__":
    unittest.main()
