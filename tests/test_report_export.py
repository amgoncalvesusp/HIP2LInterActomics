from __future__ import annotations

import re
import tempfile
import unittest
from pathlib import Path

from PIL import Image

from luna_gui.core.project import ProjectConfig
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
        self.cfg = ProjectConfig(workdir=str(self.root), protein_file="protein.pdb", ligand_file="ligands.sdf")
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


if __name__ == "__main__":
    unittest.main()
