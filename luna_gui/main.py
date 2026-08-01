"""Entry point: `python -m luna_gui.main`"""
from __future__ import annotations

import sys
import traceback
from pathlib import Path

from PyQt6.QtCore import QSettings
from PyQt6.QtGui import QIcon
from PyQt6.QtWidgets import QApplication, QMessageBox

from .ui.main_window import MainWindow
from .ui.theme import apply_theme


def _exception_hook(exc_type, exc_value, exc_traceback) -> None:
    """Expose otherwise invisible callback failures in the windowed bundle."""
    if issubclass(exc_type, KeyboardInterrupt):
        sys.__excepthook__(exc_type, exc_value, exc_traceback)
        return
    traceback.print_exception(exc_type, exc_value, exc_traceback)
    QMessageBox.critical(
        None,
        "Erro inesperado",
        "O aplicativo encontrou um erro inesperado.\n\n"
        f"{exc_type.__name__}: {exc_value}",
    )


def main() -> int:
    app = QApplication(sys.argv)
    sys.excepthook = _exception_hook
    app.setApplicationName("HIP²LInterActomics")
    app.setOrganizationName("HIP2LInterActomics")
    icon_path = Path(__file__).resolve().parent / "assets" / "hip2l_interactomics_icon.png"
    if icon_path.exists():
        app.setWindowIcon(QIcon(str(icon_path)))
    preferred_theme = str(QSettings().value("appearance/theme", "system") or "system")
    apply_theme(app, preferred_theme)
    w = MainWindow()
    w.show()
    return app.exec()


if __name__ == "__main__":
    sys.exit(main())
