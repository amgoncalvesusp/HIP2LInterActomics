"""Resolve CPU allocations for local, Slurm, and PBS/Torque execution."""

from __future__ import annotations

import os
import re
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Mapping


@dataclass(frozen=True)
class CpuAllocation:
    cpus: int
    source: str
    scheduler: str | None = None


def _positive_integer(value: object) -> int | None:
    match = re.match(r"\s*(\d+)", str(value or ""))
    if not match:
        return None
    parsed = int(match.group(1))
    return parsed if parsed > 0 else None


def _pbs_node_slots(path: str | os.PathLike[str]) -> int | None:
    try:
        hosts = [
            line.strip()
            for line in Path(path).read_text(
                encoding="utf-8", errors="replace"
            ).splitlines()
            if line.strip()
        ]
    except OSError:
        return None
    if not hosts:
        return None
    # Python multiprocessing is node-local. The largest per-host slot count is
    # therefore safer than treating a multi-node allocation as shared memory.
    return max(Counter(hosts).values())


def detect_cpu_allocation(
    environ: Mapping[str, str] | None = None,
    *,
    cpu_count: int | None = None,
) -> CpuAllocation:
    """Return the CPU budget advertised by a scheduler or the local host."""

    values = os.environ if environ is None else environ
    for variable, scheduler in (
        ("SLURM_CPUS_PER_TASK", "slurm"),
        ("SLURM_JOB_CPUS_PER_NODE", "slurm"),
        ("PBS_NUM_PPN", "pbs"),
    ):
        parsed = _positive_integer(values.get(variable))
        if parsed is not None:
            return CpuAllocation(parsed, variable, scheduler)

    nodefile = values.get("PBS_NODEFILE")
    if nodefile:
        parsed = _pbs_node_slots(nodefile)
        if parsed is not None:
            return CpuAllocation(parsed, "PBS_NODEFILE", "pbs")

    available = cpu_count if cpu_count is not None else os.cpu_count()
    return CpuAllocation(max(1, int(available or 1)), "os.cpu_count", None)


def effective_nproc(
    requested: int | None,
    environ: Mapping[str, str] | None = None,
    *,
    cpu_count: int | None = None,
) -> int:
    """Select nproc without exceeding local resources or scheduler limits.

    Scheduler allocations deliberately drive the value so a legacy config with
    ``nproc=1`` does not silently underuse the resources requested from Slurm or
    PBS. Outside a scheduler, an explicit value is capped to the host CPU count.
    """

    allocation = detect_cpu_allocation(environ, cpu_count=cpu_count)
    if allocation.scheduler:
        return allocation.cpus
    parsed = _positive_integer(requested)
    if parsed is None:
        return allocation.cpus
    return min(parsed, allocation.cpus)
