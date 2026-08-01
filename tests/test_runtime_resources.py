from __future__ import annotations

import os
import signal
import tempfile
from pathlib import Path
from unittest import mock

from luna_gui.core.process_control import TerminationController
from luna_gui.core.runtime_resources import detect_cpu_allocation, effective_nproc


def test_slurm_cpu_allocation_has_priority() -> None:
    allocation = detect_cpu_allocation(
        {
            "SLURM_CPUS_PER_TASK": "12",
            "SLURM_JOB_CPUS_PER_NODE": "24(x2)",
            "PBS_NUM_PPN": "8",
        },
        cpu_count=64,
    )

    assert allocation.cpus == 12
    assert allocation.source == "SLURM_CPUS_PER_TASK"
    assert allocation.scheduler == "slurm"


def test_slurm_compact_node_value_is_parsed() -> None:
    allocation = detect_cpu_allocation(
        {"SLURM_JOB_CPUS_PER_NODE": "16(x3)"},
        cpu_count=64,
    )

    assert allocation.cpus == 16
    assert allocation.source == "SLURM_JOB_CPUS_PER_NODE"


def test_pbs_nodefile_uses_slots_on_one_node() -> None:
    with tempfile.TemporaryDirectory() as temporary:
        nodefile = Path(temporary) / "nodes"
        nodefile.write_text("node-a\nnode-a\nnode-b\nnode-b\n", encoding="utf-8")
        allocation = detect_cpu_allocation(
            {"PBS_NODEFILE": str(nodefile)},
            cpu_count=64,
        )

    assert allocation.cpus == 2
    assert allocation.source == "PBS_NODEFILE"
    assert allocation.scheduler == "pbs"


def test_scheduler_allocation_drives_nproc_and_local_request_is_capped() -> None:
    assert effective_nproc(1, {"PBS_NUM_PPN": "10"}, cpu_count=64) == 10
    assert effective_nproc(32, {}, cpu_count=8) == 8
    assert effective_nproc(None, {}, cpu_count=6) == 6


def test_signal_controller_terminates_active_child() -> None:
    process = mock.Mock()
    process.poll.return_value = None
    controller = TerminationController()
    controller.attach(process)

    controller.handle_signal(signal.SIGTERM, None)

    process.terminate.assert_called_once_with()
    assert controller.received_signal == signal.SIGTERM
    controller.detach(process)


def test_signal_controller_installs_only_supported_signals() -> None:
    controller = TerminationController()
    with mock.patch.object(signal, "signal", wraps=signal.signal) as install:
        with controller.installed():
            assert controller.received_signal is None

    installed = {call.args[0] for call in install.call_args_list}
    assert signal.SIGINT in installed
    if hasattr(signal, "SIGTERM"):
        assert signal.SIGTERM in installed
