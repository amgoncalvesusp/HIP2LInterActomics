from __future__ import annotations

from unittest import mock

from PyQt6.QtCore import QRect
from PyQt6.QtGui import QFontDatabase
from PyQt6.QtWidgets import QGridLayout, QMainWindow, QMessageBox, QPushButton, QTabWidget

from luna_gui.ui.binding_modes_editor import BindingModesEditor
from luna_gui.ui.dialog_prep import DockingPrepDialog
from luna_gui.ui.info import InfoButton
from luna_gui.ui.main_window import MainWindow
from luna_gui.ui.tab_results_enhanced import ResultsTab
from luna_gui.ui.theme import apply_theme
from luna_gui.core.project import ProjectConfig


class _Screen:
    def availableGeometry(self) -> QRect:
        return QRect(0, 0, 800, 600)


def test_fit_to_screen_never_exceeds_available_geometry(qtbot) -> None:
    window = QMainWindow()
    qtbot.addWidget(window)
    with mock.patch("luna_gui.ui.main_window.QGuiApplication.primaryScreen", return_value=_Screen()):
        MainWindow._fit_to_screen(window)

    assert window.width() <= 800
    assert window.height() <= 600
    assert window.minimumWidth() <= 800
    assert window.minimumHeight() <= 600


def test_theme_uses_platform_general_font(qapp) -> None:
    apply_theme(qapp)
    expected = QFontDatabase.systemFont(QFontDatabase.SystemFont.GeneralFont).family()
    assert qapp.font().family() == expected


def test_results_actions_wrap_to_two_grid_rows(qtbot) -> None:
    tab = ResultsTab(ProjectConfig())
    qtbot.addWidget(tab)
    controls = tab.layout().itemAt(1).layout()
    assert isinstance(controls, QGridLayout)
    assert controls.rowCount() >= 2


def test_dialogs_have_bounded_minimum_width(qtbot) -> None:
    with (
        mock.patch("luna_gui.ui.dialog_prep.QGuiApplication.primaryScreen", return_value=_Screen()),
        mock.patch("luna_gui.ui.binding_modes_editor.QGuiApplication.primaryScreen", return_value=_Screen()),
    ):
        prep = DockingPrepDialog()
        editor = BindingModesEditor()
    qtbot.addWidget(prep)
    qtbot.addWidget(editor)
    assert prep.minimumSizeHint().width() <= 900
    assert editor.minimumSizeHint().width() <= 900
    assert prep.width() <= 800
    assert prep.height() <= 600
    assert editor.width() <= 800
    assert editor.height() <= 600


def test_info_button_is_keyboard_accessible(qtbot) -> None:
    button = InfoButton("Ajuda contextual")
    qtbot.addWidget(button)
    assert button.width() >= 28
    assert button.accessibleName()
    assert button.accessibleDescription() == "Ajuda contextual"
    with mock.patch("luna_gui.ui.info.QToolTip.showText") as show_help:
        button.click()
    show_help.assert_called_once()


def test_all_main_window_buttons_are_wired_and_horizontally_accessible(qtbot) -> None:
    with mock.patch("luna_gui.ui.tab_setup.SetupTab.detect"):
        window = MainWindow()
    qtbot.addWidget(window)
    window.resize(900, 620)
    window.show()

    issues: list[str] = []
    for outer_index in range(window.tabs.count()):
        window.tabs.setCurrentIndex(outer_index)
        nested_tabs = window.tabs.currentWidget().findChildren(QTabWidget)
        inner_indices = range(nested_tabs[0].count()) if nested_tabs else [None]
        for inner_index in inner_indices:
            if inner_index is not None:
                nested_tabs[0].setCurrentIndex(inner_index)
            for button in window.findChildren(QPushButton):
                if button.objectName().startswith("qt_"):
                    continue
                if button.receivers(button.clicked) == 0:
                    issues.append(f"sem ação: {button.text()}")
                if not button.isVisibleTo(window):
                    continue
                left = button.mapTo(window, button.rect().topLeft()).x()
                right = button.mapTo(window, button.rect().bottomRight()).x()
                if left < 0 or right >= window.width():
                    issues.append(f"fora da largura: {button.text()}")

    with mock.patch(
        "luna_gui.ui.main_window.QMessageBox.question",
        return_value=QMessageBox.StandardButton.Yes,
    ):
        window.close()
    assert issues == []
