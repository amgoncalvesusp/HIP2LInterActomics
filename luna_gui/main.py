"""Entry point: `python -m luna_gui.main`"""
from __future__ import annotations

import sys

from PyQt6.QtWidgets import QApplication

from .ui.main_window import MainWindow
from .ui.theme import apply_theme


def main() -> int:
    app = QApplication(sys.argv)
    app.setApplicationName("LUNA GUI")
    apply_theme(app)
    w = MainWindow()
    w.show()
    return app.exec()


if __name__ == "__main__":
    sys.exit(main())
