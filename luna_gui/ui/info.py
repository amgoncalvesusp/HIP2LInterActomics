"""Small visible help affordance used beside GUI options."""
from __future__ import annotations

from PyQt6.QtWidgets import QToolButton


class InfoButton(QToolButton):
    def __init__(self, tooltip: str, parent=None) -> None:
        super().__init__(parent)
        self.setText("!")
        self.setToolTip(tooltip)
        self.setAutoRaise(True)
        self.setFixedSize(22, 22)
        self.setStyleSheet(
            "QToolButton {"
            "border: 1px solid #8aa5bd; border-radius: 11px;"
            "background: #f4f8fc; color: #0b2b45; font-weight: 700;"
            "}"
            "QToolButton:hover { background: #dcebf6; }"
        )
