"""Shared Qt theme for the LUNA GUI."""
from __future__ import annotations

from PyQt6.QtGui import QColor, QFont, QPalette
from PyQt6.QtWidgets import QApplication


def apply_theme(app: QApplication) -> None:
    """Apply a consistent palette and stylesheet across the application."""
    app.setStyle("Fusion")
    app.setFont(QFont("Segoe UI", 10))

    palette = QPalette()
    palette.setColor(QPalette.ColorRole.Window, QColor("#f4efe7"))
    palette.setColor(QPalette.ColorRole.WindowText, QColor("#2e241f"))
    palette.setColor(QPalette.ColorRole.Base, QColor("#fffdfa"))
    palette.setColor(QPalette.ColorRole.AlternateBase, QColor("#efe8dd"))
    palette.setColor(QPalette.ColorRole.ToolTipBase, QColor("#fffdfa"))
    palette.setColor(QPalette.ColorRole.ToolTipText, QColor("#2e241f"))
    palette.setColor(QPalette.ColorRole.Text, QColor("#2e241f"))
    palette.setColor(QPalette.ColorRole.Button, QColor("#d8cdbd"))
    palette.setColor(QPalette.ColorRole.ButtonText, QColor("#2e241f"))
    palette.setColor(QPalette.ColorRole.Highlight, QColor("#0f766e"))
    palette.setColor(QPalette.ColorRole.HighlightedText, QColor("#fffdfa"))
    palette.setColor(QPalette.ColorRole.PlaceholderText, QColor("#8a7e71"))
    app.setPalette(palette)

    app.setStyleSheet(
        """
        QWidget {
            color: #2e241f;
        }
        QMainWindow, QDialog {
            background: #f4efe7;
        }
        QStatusBar {
            background: #eadfce;
            color: #4e4037;
        }
        QMenuBar {
            background: #eadfce;
            border-bottom: 1px solid #c9baa7;
        }
        QMenuBar::item {
            padding: 6px 10px;
            background: transparent;
        }
        QMenuBar::item:selected,
        QMenu::item:selected {
            background: #d6c7b6;
        }
        QMenu {
            background: #fffdfa;
            border: 1px solid #c9baa7;
        }
        QTabWidget::pane {
            border: 1px solid #d7ccbd;
            background: #fbf8f2;
            border-radius: 14px;
            margin-top: 8px;
        }
        QTabBar::tab {
            background: #ddd2c4;
            border: 1px solid #c8bba7;
            border-bottom: none;
            color: #4b3c31;
            padding: 10px 16px;
            margin-right: 6px;
            border-top-left-radius: 10px;
            border-top-right-radius: 10px;
        }
        QTabBar::tab:selected {
            background: #fbf8f2;
            color: #173f3d;
        }
        QGroupBox {
            background: #fffdfa;
            border: 1px solid #ddd2c3;
            border-radius: 14px;
            margin-top: 16px;
            padding: 14px;
            font-weight: 600;
        }
        QGroupBox::title {
            subcontrol-origin: margin;
            left: 14px;
            padding: 0 6px;
            color: #6d5f52;
        }
        QPushButton {
            background: #c8693a;
            color: #fffdfa;
            border: none;
            border-radius: 10px;
            padding: 8px 14px;
            font-weight: 600;
        }
        QPushButton:hover {
            background: #b85f35;
        }
        QPushButton:pressed {
            background: #9f522d;
        }
        QPushButton:disabled {
            background: #d8d0c6;
            color: #8d8174;
        }
        QLineEdit, QPlainTextEdit, QListWidget, QTableWidget, QComboBox, QSpinBox, QDoubleSpinBox {
            background: #fffdfa;
            border: 1px solid #d4c8b9;
            border-radius: 10px;
            padding: 6px 8px;
            selection-background-color: #0f766e;
            selection-color: #fffdfa;
        }
        QLineEdit:focus, QPlainTextEdit:focus, QListWidget:focus, QTableWidget:focus,
        QComboBox:focus, QSpinBox:focus, QDoubleSpinBox:focus {
            border: 1px solid #0f766e;
        }
        QPlainTextEdit {
            background: #1f2527;
            color: #e6edf0;
            border: 1px solid #30393d;
        }
        QHeaderView::section {
            background: #efe7db;
            color: #5a4b40;
            border: none;
            border-bottom: 1px solid #d4c8b9;
            padding: 6px;
            font-weight: 600;
        }
        QScrollBar:vertical {
            background: #efe7db;
            width: 12px;
            margin: 4px 0 4px 0;
            border-radius: 6px;
        }
        QScrollBar::handle:vertical {
            background: #c8b8a4;
            min-height: 24px;
            border-radius: 6px;
        }
        QScrollBar:horizontal {
            background: #efe7db;
            height: 12px;
            margin: 0 4px 0 4px;
            border-radius: 6px;
        }
        QScrollBar::handle:horizontal {
            background: #c8b8a4;
            min-width: 24px;
            border-radius: 6px;
        }
        QLabel[muted="true"] {
            color: #7d6e60;
        }
        """
    )
