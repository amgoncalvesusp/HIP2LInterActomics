from __future__ import annotations

from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_run_tab_generates_plots_before_emitting_completion() -> None:
    source = (ROOT / "luna_gui" / "ui" / "tab_run.py").read_text(encoding="utf-8")

    assert "def _start_plot_generation" in source
    assert '"--results-only"' in source
    assert 'self._phase = "plots"' in source
    plot_start = source.split("def _start_plot_generation", 1)[1]
    assert plot_start.index("save_to_workdir(self.cfg)") < plot_start.index('Path(self.cfg.workdir) / ".luna_gui.json"')
    plot_branch = source.split('if phase == "plots":', 1)[1]
    assert plot_branch.index("self.finished_ok.emit()") < plot_branch.index("def _start_plot_generation")


def test_plot_backend_emits_progress_and_renders_three_languages_at_180_300_dpi() -> None:
    terminal_source = (ROOT / "luna_gui" / "core" / "terminal_results.py").read_text(encoding="utf-8")
    profile_source = (ROOT / "luna_gui" / "core" / "plot_profiles.py").read_text(encoding="utf-8")

    assert "[plots-progress]" in terminal_source
    assert "screen 180 DPI + report 300 DPI" in terminal_source
    assert 'SCREEN_PROFILE = PlotProfile("screen", 2480, 3508, 180)' in profile_source
    assert 'REPORT_PROFILE = PlotProfile("report", 2480, 3508, 300)' in profile_source


def test_gui_pdf_export_never_pumps_qt_events_while_waiting() -> None:
    for name in ("tab_results_enhanced.py", "tab_results_advanced.py"):
        source = (ROOT / "luna_gui" / "ui" / name).read_text(encoding="utf-8")
        assert "QApplication.processEvents" not in source
    launcher = (ROOT / "run.py").read_text(encoding="utf-8")
    assert launcher.index('command == "--render-pdf-job"') < launcher.index("from luna_gui.main import main")
