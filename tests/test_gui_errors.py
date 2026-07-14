from __future__ import annotations

from pathlib import Path
from unittest import mock

from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import QListWidgetItem

from luna_gui.ui.tab_history import HistoryTab


def test_history_remove_reports_write_error(qtbot, tmp_path: Path) -> None:
    with mock.patch("luna_gui.ui.tab_history.load_history", return_value=["workdir"]):
        tab = HistoryTab()
    qtbot.addWidget(tab)
    tab.list.clear()
    item = QListWidgetItem("workdir")
    item.setData(Qt.ItemDataRole.UserRole, "workdir")
    tab.list.addItem(item)
    tab.list.setCurrentItem(item)

    with (
        mock.patch("luna_gui.ui.tab_history.load_history", return_value=["workdir"]),
        mock.patch("luna_gui.core.project.HISTORY_FILE", tmp_path / "missing" / "history.json"),
        mock.patch("luna_gui.ui.tab_history.QMessageBox.critical") as critical,
    ):
        tab.remove_selected()

    critical.assert_called_once()
