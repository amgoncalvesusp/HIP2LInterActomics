"""Small visible help affordance used beside GUI options."""
from __future__ import annotations

from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import QToolButton, QToolTip


class InfoButton(QToolButton):
    def __init__(self, tooltip: str, parent=None) -> None:
        super().__init__(parent)
        self.setText("!")
        self.setToolTip(tooltip)
        self.setAccessibleName("Ajuda")
        self.setAccessibleDescription(tooltip)
        self.setFocusPolicy(Qt.FocusPolicy.StrongFocus)
        self.clicked.connect(self._show_help)
        self.setAutoRaise(True)
        self.setFixedSize(28, 28)
        self.setStyleSheet(
            "QToolButton {"
            "border: 1px solid #8aa5bd; border-radius: 14px;"
            "background: #f4f8fc; color: #0b2b45; font-weight: 700;"
            "}"
            "QToolButton:hover { background: #dcebf6; }"
        )

    def _show_help(self) -> None:
        QToolTip.showText(self.mapToGlobal(self.rect().bottomLeft()), self.toolTip(), self)
