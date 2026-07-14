from __future__ import annotations

from unittest import mock

from PyQt6.QtCore import QProcess

from luna_gui.core import env_manager
from luna_gui.core.project import ProjectConfig
from luna_gui.ui.tab_run import RunTab
from luna_gui.ui.tab_setup import SetupTab


def test_run_tab_recovers_from_failed_start(qtbot) -> None:
    tab = RunTab(ProjectConfig())
    qtbot.addWidget(tab)
    tab.proc = QProcess(tab)
    tab.btn_run.setEnabled(False)
    tab.btn_cancel.setEnabled(True)
    tab.proc.errorOccurred.connect(tab._on_process_error)

    tab.proc.start("definitely-missing-hip2l-executable", [])
    qtbot.waitUntil(lambda: tab.proc is None, timeout=1000)

    assert tab.proc is None
    assert tab.btn_run.isEnabled()
    assert not tab.btn_cancel.isEnabled()
    assert "erro ao iniciar" in tab.log.toPlainText().lower()


def test_setup_tab_recovers_from_failed_start(qtbot) -> None:
    with mock.patch("luna_gui.ui.tab_setup.SetupTab.detect"):
        tab = SetupTab()
    qtbot.addWidget(tab)
    tab.proc = QProcess(tab)
    tab._cmd_queue = [["remaining"]]
    tab.btn_install_luna.setEnabled(False)
    tab.proc.errorOccurred.connect(tab._on_process_error)

    tab.proc.start("definitely-missing-hip2l-executable", [])
    qtbot.waitUntil(lambda: tab.proc is None, timeout=1000)

    assert tab.proc is None
    assert tab.btn_install_luna.isEnabled()
    assert tab._cmd_queue == []


def test_conda_info_has_a_finite_timeout() -> None:
    with (
        mock.patch.object(env_manager, "conda_process_env", return_value={}),
        mock.patch.object(env_manager.subprocess, "check_output", return_value="{}") as check_output,
    ):
        assert env_manager.conda_info("conda") == {}

    assert check_output.call_args.kwargs["timeout"] == 10


def test_stale_cancel_timer_never_kills_new_process(qtbot) -> None:
    tab = RunTab(ProjectConfig())
    qtbot.addWidget(tab)
    old_process = mock.Mock()
    new_process = mock.Mock()
    tab.proc = new_process

    tab._kill_if_running(old_process)

    new_process.kill.assert_not_called()


def test_stale_taskkill_callback_does_not_clear_new_helper(qtbot) -> None:
    tab = RunTab(ProjectConfig())
    qtbot.addWidget(tab)
    process = mock.Mock()
    stale_helper = mock.Mock()
    current_helper = mock.Mock()
    tab._cancel_helper = current_helper

    tab._taskkill_finished(process, stale_helper, 0)

    assert tab._cancel_helper is current_helper
    stale_helper.deleteLater.assert_called_once()
    current_helper.deleteLater.assert_not_called()
