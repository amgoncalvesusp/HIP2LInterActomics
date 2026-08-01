"""Small visible help affordance used beside GUI options."""
from __future__ import annotations

from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import QToolButton


class InfoButton(QToolButton):
    def __init__(self, tooltip: str, parent=None) -> None:
        super().__init__(parent)
        self.setText("i")
        self.setToolTip(tooltip)
        self.setAccessibleName("Informação")
        self.setProperty("infoButton", True)
        self.setAutoRaise(False)
        self.setFocusPolicy(Qt.FocusPolicy.StrongFocus)
        self.setFixedSize(34, 34)
