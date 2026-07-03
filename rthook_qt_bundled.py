"""PyInstaller runtime hook: pin Qt DLL resolution to the bundled Qt.

End-user machines often have a different Qt on ``PATH`` (conda, Schrodinger
PyMOL, system installs). Because the bundled Qt DLLs live in the subdirectory
``PyQt6/Qt6/bin`` rather than the bundle root, Windows can resolve the
transitive ``Qt6Core.dll`` dependency of ``Qt6Gui.dll`` to a *foreign,
older* copy found on ``PATH``. The version mismatch then fails as:

    ImportError: DLL load failed while importing QtGui:
    The specified procedure could not be found.

Running before any ``PyQt6`` import, this hook makes the bundled
``PyQt6/Qt6/bin`` the first directory searched for DLLs, so the matching
Qt6Core/Qt6Gui pair is always used regardless of the user's PATH.
"""

import os
import sys

if getattr(sys, "frozen", False):
    base = getattr(sys, "_MEIPASS", None)
    if base:
        qt_bin = os.path.join(base, "PyQt6", "Qt6", "bin")
        if os.path.isdir(qt_bin):
            try:
                os.add_dll_directory(qt_bin)
            except OSError:
                pass
            os.environ["PATH"] = qt_bin + os.pathsep + os.environ.get("PATH", "")
