from __future__ import annotations

import tempfile
from pathlib import Path

from luna_gui.core.plot_manifest import (
    PlotManifest,
    PlotRecord,
    load_manifest,
    manifest_path,
    plot_output_path,
)


def test_manifest_filters_strictly_by_language_and_profile() -> None:
    manifest = PlotManifest()
    manifest.add(PlotRecord("distribution", "en", "report", "en.png", "EN", "", 10))
    manifest.add(PlotRecord("distribution", "pt", "report", "pt.png", "PT", "", 10))
    manifest.add(PlotRecord("distribution", "en", "screen", "screen.png", "EN", "", 10))

    selected = manifest.select(language="en", profile="report")

    assert [item.path for item in selected] == ["en.png"]


def test_manifest_round_trip_and_semantic_sorting() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        manifest = PlotManifest()
        manifest.add(PlotRecord("clusters", "en", "report", "clusters.png", "Clusters", "", 50))
        manifest.add(PlotRecord("heatmap", "en", "report", "heatmap.png", "Heatmap", "", 20))
        manifest.save(manifest_path(tmp))

        loaded = load_manifest(tmp)

        assert [row.plot_id for row in loaded.select(language="en", profile="report")] == [
            "heatmap",
            "clusters",
        ]
        assert plot_output_path(tmp, "es", "screen", "heatmaps", "a.png") == (
            Path(tmp) / "results" / "plots" / "es" / "screen" / "heatmaps" / "a.png"
        )
