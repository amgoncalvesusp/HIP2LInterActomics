"""GUI launcher and packaged command dispatcher."""
from __future__ import annotations

from multiprocessing import freeze_support
from pathlib import Path

import os
import sys

freeze_support()


def _attach_windows_console() -> None:
    if os.name != "nt" or not getattr(sys, "frozen", False):
        return
    try:
        import ctypes

        if ctypes.windll.kernel32.AttachConsole(-1):
            sys.stdout = open("CONOUT$", "w", encoding="utf-8", buffering=1)
            sys.stderr = open("CONOUT$", "w", encoding="utf-8", buffering=1)
            sys.stdin = open("CONIN$", "r", encoding="utf-8")
    except Exception:
        pass


def _dispatch_cli(argv: list[str]) -> int | None:
    executable = Path(sys.argv[0]).stem.casefold().replace("_", "-")
    command = argv[0].casefold() if argv else ""
    if command == "--render-pdf-job":
        from luna_gui.core.report_worker import main

        return int(main(argv[1:]))
    if command == "--terminal" or executable == "hipplinteractomics-terminal":
        _attach_windows_console()
        from hipplinteractomics_terminal import main

        return int(main(argv[1:] if command == "--terminal" else argv))
    if command == "--multiple-run" or executable == "hipplinteractomics-multiple-run":
        _attach_windows_console()
        from hipplinteractomics_multiple_run import main

        return int(main(argv[1:] if command == "--multiple-run" else argv))
    return None


exit_code = _dispatch_cli(sys.argv[1:])
if exit_code is None:
    from luna_gui.main import main

    exit_code = int(main())
sys.exit(exit_code)
