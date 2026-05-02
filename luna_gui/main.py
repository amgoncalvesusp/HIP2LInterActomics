"""Entry point: `python -m luna_gui.main`"""
from __future__ import annotations

import sys
from pathlib import Path

from PyQt6.QtGui import QIcon
from PyQt6.QtWidgets import QApplication

from .ui.main_window import MainWindow
from .ui.theme import apply_theme


def main() -> int:
    app = QApplication(sys.argv)
    app.setApplicationName("HIP2LInterActomics_GUI")
    icon_path = Path(__file__).resolve().parent / "assets" / "hip2l_interactomics_icon.svg"
    if icon_path.exists():
        app.setWindowIcon(QIcon(str(icon_path)))
    apply_theme(app)
    w = MainWindow()
    w.show()
    return app.exec()


if __name__ == "__main__":
    sys.exit(main())
