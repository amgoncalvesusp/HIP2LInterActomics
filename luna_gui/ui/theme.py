"""Shared Qt theme for the LUNA GUI."""
from __future__ import annotations

from PyQt6.QtGui import QColor, QFont, QPalette
from PyQt6.QtWidgets import QApplication


def apply_theme(app: QApplication) -> None:
    """Apply a consistent palette and stylesheet across the application."""
    app.setStyle("Fusion")
    app.setFont(QFont("Segoe UI", 10))

    palette = QPalette()
    palette.setColor(QPalette.ColorRole.Window, QColor("#ebe7df"))
    palette.setColor(QPalette.ColorRole.WindowText, QColor("#061a32"))
    palette.setColor(QPalette.ColorRole.Base, QColor("#fbfcff"))
    palette.setColor(QPalette.ColorRole.AlternateBase, QColor("#e7e2d9"))
    palette.setColor(QPalette.ColorRole.ToolTipBase, QColor("#fbfcff"))
    palette.setColor(QPalette.ColorRole.ToolTipText, QColor("#061a32"))
    palette.setColor(QPalette.ColorRole.Text, QColor("#061a32"))
    palette.setColor(QPalette.ColorRole.Button, QColor("#236893"))
    palette.setColor(QPalette.ColorRole.ButtonText, QColor("#ffffff"))
    palette.setColor(QPalette.ColorRole.Highlight, QColor("#236893"))
    palette.setColor(QPalette.ColorRole.HighlightedText, QColor("#ffffff"))
    palette.setColor(QPalette.ColorRole.PlaceholderText, QColor("#7c8ca0"))
    app.setPalette(palette)

    app.setStyleSheet(
        """
        QWidget {
            color: #061a32;
        }
        QMainWindow, QDialog {
            background: #ebe7df;
        }
        QStatusBar {
            background: #e3ded5;
            color: #263a52;
        }
        QMenuBar {
            background: #ebe7df;
            border-bottom: 1px solid #b6c7d8;
        }
        QMenuBar::item {
            padding: 6px 10px;
            background: transparent;
        }
        QMenuBar::item:selected,
        QMenu::item:selected {
            background: #d8e7f2;
        }
        QMenu {
            background: #fbfcff;
            border: 1px solid #94acc4;
        }
        QTabWidget::pane {
            border: 1px solid #9fb7cf;
            background: #fbfcff;
            border-radius: 12px;
            margin-top: 8px;
        }
        QTabBar::tab {
            background: #e0d9ce;
            border: 1px solid #b9c7d6;
            border-bottom: none;
            color: #061a32;
            padding: 9px 16px;
            margin-right: 6px;
            border-top-left-radius: 9px;
            border-top-right-radius: 9px;
        }
        QTabBar::tab:selected {
            background: #fbfcff;
            color: #061a32;
            border-top: 3px solid #e4b44f;
        }
        QGroupBox {
            background: #fbfcff;
            border: 1px solid #a8bed4;
            border-radius: 14px;
            margin-top: 16px;
            padding: 10px;
            font-weight: 600;
        }
        QGroupBox::title {
            subcontrol-origin: margin;
            left: 14px;
            padding: 0 6px;
            color: #061a32;
        }
        QPushButton {
            background: #236893;
            color: #ffffff;
            border: none;
            border-radius: 10px;
            padding: 6px 12px;
            min-height: 24px;
            font-weight: 600;
        }
        QPushButton:hover {
            background: #1d5b82;
        }
        QPushButton:pressed {
            background: #164969;
        }
        QPushButton:disabled {
            background: #cbd5df;
            color: #738398;
        }
        QLineEdit, QPlainTextEdit, QListWidget, QTableWidget, QComboBox, QSpinBox, QDoubleSpinBox {
            background: #fbfcff;
            border: 1px solid #88a3bd;
            border-radius: 10px;
            padding: 5px 8px;
            min-height: 24px;
            selection-background-color: #236893;
            selection-color: #ffffff;
        }
        QScrollArea {
            background: #ebe7df;
            border: none;
        }
        QLineEdit:focus, QPlainTextEdit:focus, QListWidget:focus, QTableWidget:focus,
        QComboBox:focus, QSpinBox:focus, QDoubleSpinBox:focus {
            border: 1px solid #236893;
        }
        QPlainTextEdit {
            background: #f7f9fc;
            color: #061a32;
            border: 1px solid #88a3bd;
        }
        QHeaderView::section {
            background: #e7e2d9;
            color: #061a32;
            border: none;
            border-bottom: 1px solid #a8bed4;
            padding: 6px;
            font-weight: 600;
        }
        QScrollBar:vertical {
            background: #e7e2d9;
            width: 12px;
            margin: 4px 0 4px 0;
            border-radius: 6px;
        }
        QScrollBar::handle:vertical {
            background: #8aa4bc;
            min-height: 24px;
            border-radius: 6px;
        }
        QScrollBar:horizontal {
            background: #e7e2d9;
            height: 12px;
            margin: 0 4px 0 4px;
            border-radius: 6px;
        }
        QScrollBar::handle:horizontal {
            background: #8aa4bc;
            min-width: 24px;
            border-radius: 6px;
        }
        QLabel[muted="true"] {
            color: #314c6b;
        }
        """
    )
