from __future__ import annotations

import threading
from types import SimpleNamespace
from unittest import mock

from PyQt6.QtCore import QTimer

from luna_gui.ui.async_task import AsyncTask
from luna_gui.core.project import ProjectConfig
from luna_gui.ui.tab_results_advanced import ResultsTab, _FpDashboardRequest
from luna_gui.ui.dialog_prep import DockingPrepDialog
from luna_gui.ui.tab_setup import SetupTab


def test_async_task_keeps_qt_event_loop_responsive(qtbot) -> None:
    release = threading.Event()
    ticks: list[bool] = []
    timer = QTimer()
    timer.timeout.connect(lambda: ticks.append(True))
    timer.start(5)

    task = AsyncTask(lambda: release.wait(1.0) or "done")
    task.start()

    qtbot.waitUntil(lambda: len(ticks) >= 2, timeout=500)
    with qtbot.waitSignal(task.succeeded, timeout=1000) as blocker:
        release.set()

    assert blocker.args == [True]
    qtbot.waitUntil(lambda: not task.is_running, timeout=1000)
    timer.stop()


def test_async_task_reports_exception_and_always_finishes(qtbot) -> None:
    def fail() -> None:
        raise RuntimeError("boom")

    task = AsyncTask(fail)
    errors: list[str] = []
    finished: list[bool] = []
    task.failed.connect(errors.append)
    task.finished.connect(lambda: finished.append(True))

    task.start()
    qtbot.waitUntil(lambda: bool(finished), timeout=1000)

    assert errors == ["boom"]
    assert finished == [True]


def test_results_statistics_run_outside_gui_thread(qtbot, tmp_path) -> None:
    release = threading.Event()

    def run_analysis(_python: str, _workdir: str) -> dict:
        release.wait(1.0)
        return {"entry_interaction_counts": {}}

    cfg = ProjectConfig(workdir=str(tmp_path))
    tab = ResultsTab(cfg)
    qtbot.addWidget(tab)
    tab.py_exe = "python"

    with (
        mock.patch("luna_gui.ui.tab_results_advanced.load_analysis_summary", return_value=None),
        mock.patch("luna_gui.ui.tab_results_advanced.run_analysis", side_effect=run_analysis),
        mock.patch.object(tab, "_handle_stats_result") as handle_result,
    ):
        tab.compute_stats()
        assert tab._async_tasks["statistics"].is_running
        release.set()
        qtbot.waitUntil(lambda: handle_result.called, timeout=1000)
        qtbot.waitUntil(lambda: "statistics" not in tab._async_tasks, timeout=1000)


def test_complex_preparation_runs_outside_gui_thread(qtbot, tmp_path) -> None:
    source = tmp_path / "source"
    source.mkdir()
    (source / "complex.mol2").write_text("@<TRIPOS>MOLECULE\n", encoding="utf-8")
    release = threading.Event()
    result = SimpleNamespace(
        files_processed=1,
        proteins_written=1,
        ligands_written=1,
        water_molecules_detected=0,
        protein_dir=str(tmp_path / "proteins"),
        ligand_dir=str(tmp_path / "ligands"),
        log_file="",
        errors=[],
    )

    def prepare(*_args, **_kwargs):
        release.wait(1.0)
        return result

    dialog = DockingPrepDialog()
    qtbot.addWidget(dialog)
    dialog.src_edit.setText(str(source))
    dialog.rb_manual.setChecked(True)

    with mock.patch("luna_gui.ui.dialog_prep.split_complex_folder", side_effect=prepare):
        dialog._run()
        assert dialog._prep_task is not None
        assert dialog._prep_task.is_running
        release.set()
        qtbot.waitUntil(lambda: dialog._prep_task is None, timeout=1000)

    assert dialog.result_ligand_dir == result.ligand_dir


def test_setup_detection_does_not_block_qt_event_loop(qtbot) -> None:
    release = threading.Event()
    ticks: list[bool] = []

    def find_conda():
        release.wait(1.0)
        return None

    timer = QTimer()
    timer.timeout.connect(lambda: ticks.append(True))
    timer.start(5)
    with mock.patch("luna_gui.ui.tab_setup.em.find_conda", side_effect=find_conda):
        tab = SetupTab()
        qtbot.addWidget(tab)
        qtbot.waitUntil(lambda: len(ticks) >= 2, timeout=500)
        release.set()
        qtbot.waitUntil(lambda: tab._detect_task is None, timeout=1000)

    timer.stop()
    assert "Conda não encontrado" in tab.status_label.text()


def test_results_discards_async_result_after_workdir_changes(qtbot, tmp_path) -> None:
    first = tmp_path / "first"
    second = tmp_path / "second"
    first.mkdir()
    second.mkdir()
    release = threading.Event()
    applied: list[dict] = []
    tab = ResultsTab(ProjectConfig(workdir=str(first)))
    qtbot.addWidget(tab)
    tab.wd_edit.setText(str(first))

    tab._run_async(
        "context_test",
        lambda: release.wait(1.0) or {"source": "first"},
        applied.append,
        lambda _message: None,
    )
    tab.wd_edit.setText(str(second))
    release.set()
    qtbot.waitUntil(lambda: "context_test" not in tab._async_tasks, timeout=1000)

    assert applied == []


def test_fp_dashboard_computation_runs_outside_gui_thread(qtbot, tmp_path) -> None:
    release = threading.Event()
    ticks: list[bool] = []
    request = _FpDashboardRequest(
        ifp_type="EIFP",
        workdir=tmp_path,
        artifact={},
        labels_csv="",
        labels_id_column="",
        labels_column="",
        task_kind_preference="classification",
        algorithm_preference="gradient_boosting",
        use_otsu_threshold=False,
        random_seed=7,
        python_executable="python",
        cache_key=("EIFP", "test"),
    )
    tab = ResultsTab(ProjectConfig(workdir=str(tmp_path)))
    qtbot.addWidget(tab)
    tab.cb_fp_analysis_type.addItem("EIFP", "EIFP")
    timer = QTimer()
    timer.timeout.connect(lambda: ticks.append(True))
    timer.start(5)

    def compute(_request):
        release.wait(1.0)
        return {"model_name": "GradientBoosting"}

    with (
        mock.patch.object(tab, "_fp_dashboard_request", return_value=request),
        mock.patch("luna_gui.ui.tab_results_advanced._compute_fp_dashboard", side_effect=compute),
        mock.patch.object(tab, "_apply_fp_dashboard") as apply_dashboard,
    ):
        tab._render_fp_analysis_table()
        qtbot.waitUntil(lambda: len(ticks) >= 2, timeout=500)
        release.set()
        qtbot.waitUntil(lambda: apply_dashboard.called, timeout=1000)

    timer.stop()
