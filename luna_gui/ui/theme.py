"""Shared light/dark Qt theme for HIP²LInterActomics."""
from __future__ import annotations

from PyQt6.QtCore import Qt
from PyQt6.QtGui import QColor, QPalette
from PyQt6.QtWidgets import QApplication


THEME_MODES = ("system", "light", "dark")

_LIGHT = {
    "window": "#F3F6F9",
    "surface": "#FFFFFF",
    "surface_alt": "#EAF0F5",
    "text": "#172033",
    "muted": "#52657A",
    "border": "#B8C7D6",
    "border_strong": "#849AAF",
    "primary": "#126B78",
    "primary_hover": "#0E5A66",
    "primary_pressed": "#0A4852",
    "primary_text": "#FFFFFF",
    "accent": "#C88719",
    "focus": "#172033",
    "selection": "#126B78",
    "disabled": "#DCE4EB",
    "disabled_text": "#718196",
    "danger": "#B42318",
}

_DARK = {
    "window": "#0D1520",
    "surface": "#162231",
    "surface_alt": "#1D2C3D",
    "text": "#EDF3F8",
    "muted": "#B3C1D0",
    "border": "#3B5065",
    "border_strong": "#58718A",
    "primary": "#197A89",
    "primary_hover": "#218B9B",
    "primary_pressed": "#126673",
    "primary_text": "#FFFFFF",
    "accent": "#F0B94D",
    "focus": "#FFFFFF",
    "selection": "#197A89",
    "disabled": "#2A3949",
    "disabled_text": "#8798AA",
    "danger": "#FF8A80",
}


def _resolved_mode(app: QApplication, mode: str) -> str:
    mode = str(mode or "system").strip().lower()
    if mode not in THEME_MODES:
        mode = "system"
    if mode != "system":
        return mode
    try:
        return (
            "dark"
            if app.styleHints().colorScheme() == Qt.ColorScheme.Dark
            else "light"
        )
    except (AttributeError, RuntimeError):
        return "light"


def apply_theme(app: QApplication, mode: str = "system") -> str:
    """Apply the requested theme and return the resolved light/dark mode."""
    resolved = _resolved_mode(app, mode)
    colors = _DARK if resolved == "dark" else _LIGHT
    app.setStyle("Fusion")
    font = app.font()
    font.setPointSize(10)
    app.setFont(font)

    palette = QPalette()
    palette.setColor(QPalette.ColorRole.Window, QColor(colors["window"]))
    palette.setColor(QPalette.ColorRole.WindowText, QColor(colors["text"]))
    palette.setColor(QPalette.ColorRole.Base, QColor(colors["surface"]))
    palette.setColor(QPalette.ColorRole.AlternateBase, QColor(colors["surface_alt"]))
    palette.setColor(QPalette.ColorRole.ToolTipBase, QColor(colors["surface"]))
    palette.setColor(QPalette.ColorRole.ToolTipText, QColor(colors["text"]))
    palette.setColor(QPalette.ColorRole.Text, QColor(colors["text"]))
    palette.setColor(QPalette.ColorRole.Button, QColor(colors["primary"]))
    palette.setColor(QPalette.ColorRole.ButtonText, QColor(colors["primary_text"]))
    palette.setColor(QPalette.ColorRole.Highlight, QColor(colors["selection"]))
    palette.setColor(QPalette.ColorRole.HighlightedText, QColor(colors["primary_text"]))
    palette.setColor(QPalette.ColorRole.PlaceholderText, QColor(colors["muted"]))
    palette.setColor(
        QPalette.ColorGroup.Disabled,
        QPalette.ColorRole.Text,
        QColor(colors["disabled_text"]),
    )
    app.setPalette(palette)
    app.setProperty("resolvedTheme", resolved)

    app.setStyleSheet(
        f"""
        QWidget {{
            color: {colors['text']};
        }}
        QMainWindow, QDialog {{
            background: {colors['window']};
        }}
        QToolTip {{
            background: {colors['surface']};
            color: {colors['text']};
            border: 1px solid {colors['border_strong']};
            padding: 5px 7px;
        }}
        QStatusBar {{
            background: {colors['surface_alt']};
            color: {colors['muted']};
            border-top: 1px solid {colors['border']};
            padding: 2px 8px;
        }}
        QMenuBar {{
            background: {colors['window']};
            border-bottom: 1px solid {colors['border']};
            padding: 2px 6px;
        }}
        QMenuBar::item {{
            padding: 7px 10px;
            border-radius: 6px;
            background: transparent;
        }}
        QMenuBar::item:selected, QMenu::item:selected {{
            background: {colors['surface_alt']};
        }}
        QMenu {{
            background: {colors['surface']};
            border: 1px solid {colors['border']};
            padding: 5px;
        }}
        QMenu::item {{
            padding: 7px 24px 7px 10px;
            border-radius: 5px;
        }}
        QTabWidget::pane {{
            border: 1px solid {colors['border']};
            background: {colors['surface']};
            border-radius: 14px;
            margin-top: 7px;
        }}
        QTabBar::tab {{
            background: {colors['surface_alt']};
            border: 1px solid transparent;
            color: {colors['muted']};
            padding: 9px 15px;
            margin-right: 5px;
            border-radius: 8px;
        }}
        QTabBar::tab:hover {{
            color: {colors['text']};
            border-color: {colors['border']};
        }}
        QTabBar::tab:selected {{
            background: {colors['surface']};
            color: {colors['text']};
            border-color: {colors['border_strong']};
            border-bottom: 3px solid {colors['accent']};
            font-weight: 600;
        }}
        QGroupBox {{
            background: {colors['surface']};
            border: 1px solid {colors['border']};
            border-radius: 13px;
            margin-top: 18px;
            padding: 13px;
            font-weight: 600;
        }}
        QGroupBox::title {{
            subcontrol-origin: margin;
            left: 13px;
            padding: 0 7px;
            color: {colors['text']};
            background: {colors['surface']};
        }}
        QPushButton {{
            background: {colors['primary']};
            color: {colors['primary_text']};
            border: 1px solid {colors['primary']};
            border-radius: 8px;
            padding: 5px 13px;
            min-height: 28px;
            font-weight: 600;
        }}
        QPushButton:hover {{
            background: {colors['primary_hover']};
            border-color: {colors['primary_hover']};
        }}
        QPushButton:pressed {{
            background: {colors['primary_pressed']};
            border-color: {colors['primary_pressed']};
        }}
        QPushButton:focus {{
            border: 2px solid {colors['focus']};
        }}
        QPushButton:disabled {{
            background: {colors['disabled']};
            border-color: {colors['disabled']};
            color: {colors['disabled_text']};
        }}
        QToolButton[infoButton="true"] {{
            background: {colors['surface_alt']};
            color: {colors['text']};
            border: 1px solid {colors['border_strong']};
            border-radius: 16px;
            font-weight: 700;
        }}
        QToolButton[infoButton="true"]:hover {{
            background: {colors['surface']};
            border-color: {colors['text']};
        }}
        QToolButton[infoButton="true"]:focus {{
            border: 2px solid {colors['focus']};
        }}
        QLineEdit, QPlainTextEdit, QTextEdit, QListWidget, QTableWidget,
        QComboBox, QSpinBox, QDoubleSpinBox {{
            background: {colors['surface']};
            border: 1px solid {colors['border_strong']};
            border-radius: 8px;
            padding: 5px 8px;
            min-height: 27px;
            selection-background-color: {colors['selection']};
            selection-color: {colors['primary_text']};
        }}
        QLineEdit:focus, QPlainTextEdit:focus, QTextEdit:focus, QListWidget:focus,
        QTableWidget:focus, QComboBox:focus, QSpinBox:focus, QDoubleSpinBox:focus {{
            border: 2px solid {colors['focus']};
        }}
        QPlainTextEdit, QTextEdit {{
            background: {colors['surface_alt']};
        }}
        QTableWidget {{
            gridline-color: {colors['border']};
            alternate-background-color: {colors['surface_alt']};
        }}
        QHeaderView::section {{
            background: {colors['surface_alt']};
            color: {colors['text']};
            border: none;
            border-bottom: 1px solid {colors['border']};
            padding: 7px;
            font-weight: 600;
        }}
        QCheckBox, QRadioButton {{
            spacing: 7px;
            min-height: 24px;
        }}
        QScrollArea {{
            background: {colors['window']};
            border: none;
        }}
        QScrollBar:vertical {{
            background: {colors['surface_alt']};
            width: 11px;
            margin: 3px;
            border-radius: 5px;
        }}
        QScrollBar::handle:vertical {{
            background: {colors['border_strong']};
            min-height: 28px;
            border-radius: 4px;
        }}
        QScrollBar:horizontal {{
            background: {colors['surface_alt']};
            height: 11px;
            margin: 3px;
            border-radius: 5px;
        }}
        QScrollBar::handle:horizontal {{
            background: {colors['border_strong']};
            min-width: 28px;
            border-radius: 4px;
        }}
        QScrollBar::add-line, QScrollBar::sub-line {{
            width: 0;
            height: 0;
        }}
        QLabel[muted="true"] {{
            color: {colors['muted']};
        }}
        QLabel#workflowSummary {{
            background: {colors['surface_alt']};
            color: {colors['text']};
            border: 1px solid {colors['border']};
            border-radius: 8px;
            padding: 10px;
            font-weight: 600;
        }}
        QLabel[error="true"] {{
            color: {colors['danger']};
        }}
        """
    )
    return resolved
