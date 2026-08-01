"""Measure native startup memory and eager Matplotlib canvas creation."""
from __future__ import annotations

import argparse
import ctypes
import json
import os
import sys
import time
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))


def _rss_mb() -> float:
    if sys.platform == "win32":
        from ctypes import wintypes

        class ProcessMemoryCounters(ctypes.Structure):
            _fields_ = [
                ("cb", wintypes.DWORD),
                ("PageFaultCount", wintypes.DWORD),
                ("PeakWorkingSetSize", ctypes.c_size_t),
                ("WorkingSetSize", ctypes.c_size_t),
                ("QuotaPeakPagedPoolUsage", ctypes.c_size_t),
                ("QuotaPagedPoolUsage", ctypes.c_size_t),
                ("QuotaPeakNonPagedPoolUsage", ctypes.c_size_t),
                ("QuotaNonPagedPoolUsage", ctypes.c_size_t),
                ("PagefileUsage", ctypes.c_size_t),
                ("PeakPagefileUsage", ctypes.c_size_t),
            ]

        counters = ProcessMemoryCounters()
        counters.cb = ctypes.sizeof(counters)
        ctypes.windll.kernel32.GetCurrentProcess.restype = wintypes.HANDLE
        ctypes.windll.psapi.GetProcessMemoryInfo.argtypes = (
            wintypes.HANDLE,
            ctypes.POINTER(ProcessMemoryCounters),
            wintypes.DWORD,
        )
        if not ctypes.windll.psapi.GetProcessMemoryInfo(
            ctypes.windll.kernel32.GetCurrentProcess(),
            ctypes.byref(counters),
            counters.cb,
        ):
            raise ctypes.WinError()
        return counters.WorkingSetSize / (1024 * 1024)

    import resource

    maximum = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    divisor = 1024 if sys.platform.startswith("linux") else 1024 * 1024
    return maximum / divisor


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--open-results",
        action="store_true",
        help="Include construction of the on-demand results workspace.",
    )
    args = parser.parse_args()
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    started = time.perf_counter()

    from PyQt6.QtWidgets import QApplication

    from luna_gui.ui.main_window import MainWindow

    app = QApplication([])
    window = MainWindow()
    if args.open_results:
        window._ensure_results_tab()
    app.processEvents()
    startup_seconds = time.perf_counter() - started
    rss_mb = _rss_mb()
    canvases = sum(
        widget.__class__.__name__ == "FigureCanvasQTAgg"
        for widget in QApplication.allWidgets()
    )
    print(
        json.dumps(
            {
                "startup_seconds": round(startup_seconds, 3),
                "rss_mb": round(rss_mb, 1),
                "matplotlib_canvases": canvases,
            },
            sort_keys=True,
        ),
        flush=True,
    )
    probe_thread = window.tab_setup._probe_thread
    if probe_thread is not None:
        probe_thread.wait(60_000)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
