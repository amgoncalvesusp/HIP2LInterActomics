"""Graceful signal handling for synchronous child processes."""

from __future__ import annotations

import signal
import subprocess
from contextlib import contextmanager
from types import FrameType
from typing import Iterator


class TerminationController:
    """Record termination requests and relay them to the active child."""

    def __init__(self) -> None:
        self.received_signal: int | None = None
        self._process: subprocess.Popen[str] | None = None

    def attach(self, process: subprocess.Popen[str]) -> None:
        self._process = process
        if self.received_signal is not None:
            self._terminate(process)

    def detach(self, process: subprocess.Popen[str]) -> None:
        if self._process is process:
            self._process = None

    def handle_signal(self, signum: int, _frame: FrameType | None) -> None:
        if self.received_signal is None:
            self.received_signal = signum
        process = self._process
        if process is not None:
            self._terminate(process)

    @staticmethod
    def _terminate(process: subprocess.Popen[str]) -> None:
        try:
            if process.poll() is None:
                process.terminate()
        except OSError:
            pass

    @contextmanager
    def installed(self) -> Iterator["TerminationController"]:
        previous: dict[int, object] = {}
        supported = [signal.SIGINT]
        if hasattr(signal, "SIGTERM"):
            supported.append(signal.SIGTERM)
        try:
            for signum in supported:
                previous[signum] = signal.getsignal(signum)
                signal.signal(signum, self.handle_signal)
            yield self
        finally:
            for signum, handler in previous.items():
                signal.signal(signum, handler)


def signal_exit_code(signum: int | None) -> int:
    return 128 + signum if signum is not None else 0
