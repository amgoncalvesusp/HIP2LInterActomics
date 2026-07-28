# -*- mode: python ; coding: utf-8 -*-
"""PyInstaller spec for the HIP2LInterActomics GUI (native Windows .exe).

Builds a self-contained native Windows application so lay users no longer need
WSL/WSLg, an X11 server, or the Qt ``xcb`` Linux platform plugin. On Windows the
bundle ships Qt's native ``qwindows.dll`` platform plugin (collected automatically
by PyInstaller's PyQt6 hook), which is what failed under WSLg.

The GUI still drives the LUNA calculation engine through a separate conda
environment (``luna-env``) via subprocess; this bundle packages only the
front-end. The ``luna_gui`` source tree is shipped as on-disk data because the
external luna-env interpreter imports ``luna_gui.core.results_analysis`` for the
fingerprint dashboard (see ``analysis_runtime._app_root``).
"""

import os

from PyInstaller.utils.hooks import collect_data_files

block_cipher = None

# Ship the package source on disk (not only frozen bytecode) so the external
# luna-env interpreter can import luna_gui.core.results_analysis, and so the GUI
# can read its bundled assets/examples via Path(__file__).
datas = [("luna_gui", "luna_gui")]
datas += collect_data_files("matplotlib", subdir="mpl-data")

hiddenimports = [
    "matplotlib.backends.backend_qtagg",
]

a = Analysis(
    ["run.py"],
    pathex=[],
    binaries=[],
    datas=datas,
    hiddenimports=hiddenimports,
    hookspath=[],
    hooksconfig={},
    runtime_hooks=["rthook_qt_bundled.py"],
    excludes=["PyQt5", "PySide2", "PySide6", "tkinter"],
    win_no_prefer_redirects=False,
    win_private_assemblies=False,
    cipher=block_cipher,
    noarchive=False,
)

# Drop the stray ICU DLLs (icuuc/icudt/icuin) that PyInstaller pulls in from a
# foreign Qt/conda install on PATH to satisfy Qt6Core's optional ICU import.
# That copy exports versioned symbols (ucnv_open_75) instead of the unversioned
# names Qt6Core links against, so its presence breaks the QtGui DLL load with
# "procedure not found". PyQt6 6.11's Qt delay-loads ICU and runs without it.
a.binaries = [
    b for b in a.binaries
    if not os.path.basename(b[0]).lower().startswith(("icuuc", "icudt", "icuin"))
]

pyz = PYZ(a.pure, a.zipped_data, cipher=block_cipher)

exe = EXE(
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,
    name="HIP2LInterActomics",
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=False,
    console=False,
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
    icon="luna_gui/assets/hip2l_interactomics_icon.ico",
)

coll = COLLECT(
    exe,
    a.binaries,
    a.datas,
    strip=False,
    upx=False,
    upx_exclude=[],
    name="HIP2LInterActomics",
)
