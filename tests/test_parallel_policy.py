from unittest import mock

from luna_gui.core import luna_runner


def test_windows_native_nproc_is_clamped_to_one():
    with mock.patch.object(luna_runner.sys, "platform", "win32"):
        assert luna_runner.safe_nproc(8) == 1


def test_linux_nproc_keeps_requested_parallelism():
    with mock.patch.object(luna_runner.sys, "platform", "linux"):
        assert luna_runner.safe_nproc(8) == 8


def test_invalid_nproc_falls_back_to_one():
    with mock.patch.object(luna_runner.sys, "platform", "linux"):
        assert luna_runner.safe_nproc(0) == 1
