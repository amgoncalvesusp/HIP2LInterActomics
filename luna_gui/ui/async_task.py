"""Small Qt worker for running blocking callables outside the GUI thread."""
from __future__ import annotations

from collections.abc import Callable
from typing import Any

from PyQt6.QtCore import QObject, QThread, pyqtSignal, pyqtSlot


class _Worker(QObject):
    succeeded = pyqtSignal(object)
    failed = pyqtSignal(str)
    completed = pyqtSignal()

    def __init__(self, function: Callable[[], Any]) -> None:
        super().__init__()
        self._function = function

    @pyqtSlot()
    def run(self) -> None:
        try:
            result = self._function()
        except Exception as exc:
            self.failed.emit(str(exc))
        else:
            self.succeeded.emit(result)
        finally:
            self.completed.emit()


class AsyncTask(QObject):
    """Own one worker thread and relay its result back to the Qt main thread."""

    succeeded = pyqtSignal(object)
    failed = pyqtSignal(str)
    progress = pyqtSignal(object)
    finished = pyqtSignal()

    def __init__(self, function: Callable[[], Any], parent: QObject | None = None) -> None:
        super().__init__(parent)
        self._function = function
        self._thread: QThread | None = None
        self._worker: _Worker | None = None

    @property
    def is_running(self) -> bool:
        return self._thread is not None and self._thread.isRunning()

    def start(self) -> None:
        if self.is_running:
            raise RuntimeError("A tarefa já está em execução.")

        thread = QThread(self)
        worker = _Worker(self._function)
        worker.moveToThread(thread)

        thread.started.connect(worker.run)
        worker.succeeded.connect(self.succeeded)
        worker.failed.connect(self.failed)
        worker.completed.connect(thread.quit)
        worker.completed.connect(worker.deleteLater)
        thread.finished.connect(self._on_finished)

        self._thread = thread
        self._worker = worker
        thread.start()

    @pyqtSlot()
    def _on_finished(self) -> None:
        thread = self._thread
        self._worker = None
        self._thread = None
        if thread is not None:
            thread.deleteLater()
        self.finished.emit()
