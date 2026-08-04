from __future__ import annotations

import tempfile
import zipfile
from pathlib import Path

import pytest

from luna_gui.core.results_archive import ResultsArchiveError, extract_results_archive


def test_extract_results_archive_finds_nested_workdir() -> None:
    with tempfile.TemporaryDirectory() as temporary:
        root = Path(temporary)
        archive = root / "portable.zip"
        with zipfile.ZipFile(archive, "w") as handle:
            handle.writestr("project/.luna_gui.json", "{}")
            handle.writestr("project/results/analysis_summary.json", "{}")

        workdir = extract_results_archive(archive, root / "extracted")

        assert workdir == root / "extracted" / "project"
        assert (workdir / "results" / "analysis_summary.json").is_file()


def test_extract_results_archive_rejects_path_traversal_and_cleans_output() -> None:
    with tempfile.TemporaryDirectory() as temporary:
        root = Path(temporary)
        archive = root / "unsafe.zip"
        destination = root / "extracted"
        with zipfile.ZipFile(archive, "w") as handle:
            handle.writestr("../outside.txt", "unsafe")

        with pytest.raises(ResultsArchiveError, match="Unsafe archive path"):
            extract_results_archive(archive, destination)

        assert not destination.exists()
        assert not (root / "outside.txt").exists()
