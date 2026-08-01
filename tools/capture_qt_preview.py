"""Capture deterministic offscreen previews of the native Qt interface."""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("output", type=Path)
    parser.add_argument("--theme", choices=("light", "dark"), default="light")
    parser.add_argument("--width", type=int, default=1280)
    parser.add_argument("--height", type=int, default=820)
    args = parser.parse_args()

    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    from PyQt6.QtGui import QFont, QFontDatabase
    from PyQt6.QtWidgets import QApplication

    from luna_gui.ui.main_window import MainWindow
    from luna_gui.ui.theme import apply_theme

    app = QApplication([])
    if sys.platform == "win32" and not QFontDatabase.families():
        system_font = Path(os.environ.get("WINDIR", "C:/Windows")) / "Fonts" / "segoeui.ttf"
        if system_font.exists() and QFontDatabase.addApplicationFont(str(system_font)) >= 0:
            app.setFont(QFont("Segoe UI", 10))
    apply_theme(app, args.theme)
    window = MainWindow()
    window.resize(args.width, args.height)
    window.show()

    probe_thread = window.tab_setup._probe_thread
    if probe_thread is not None:
        probe_thread.wait(60_000)
    app.processEvents()

    args.output.parent.mkdir(parents=True, exist_ok=True)
    if not window.grab().save(str(args.output)):
        raise RuntimeError(f"Could not save preview to {args.output}")
    print(args.output.resolve())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
