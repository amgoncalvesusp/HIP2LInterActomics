"""Generate configuration matrices and run HIP2LInterActomics serially.

This module is fully headless: it imports neither PyQt nor any GUI entry point.
It always writes all JSON configurations before launching the first child.
"""

from __future__ import annotations

import argparse
import ast
import copy
import hashlib
import itertools
import json
import os
import shlex
import shutil
import subprocess
import sys
import tempfile
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Callable, Iterable, Sequence

from luna_gui.core import ligand_io
from luna_gui.core.process_control import TerminationController, signal_exit_code
from luna_gui.core.runtime_resources import detect_cpu_allocation, effective_nproc


_FORMAT_ALIASES = {
    "binary": "bin",
    "bin": "bin",
    "count": "cnt",
    "cnt": "cnt",
}


class BatchConfigurationError(ValueError):
    """Raised when the requested batch matrix is invalid."""


class BatchExecutionError(RuntimeError):
    """Raised when a child process cannot start or exits unsuccessfully."""


@dataclass(frozen=True)
class GeneratedRun:
    run_id: str
    bits: int
    fingerprint_format: str
    level: int
    growth_ratio: float
    config_path: Path
    config_digest: str


@dataclass(frozen=True)
class RunResult:
    run_id: str
    config_path: str
    log_path: str
    command: list[str]
    returncode: int
    duration_seconds: float
    config_digest: str
    status: str


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="hipplinteractomics-multiple-run",
        description=(
            "Generate a Cartesian matrix of JSON configs, then run each one "
            "serially through the headless terminal command."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--base-config",
        type=Path,
        required=True,
        metavar="PATH",
        help="Base JSON config inherited by every generated run.",
    )
    parser.add_argument(
        "--bits",
        nargs="+",
        required=True,
        metavar="BITS",
        help='Fingerprint lengths, e.g. 1024 2048 or "[1024, 2048]".',
    )
    parser.add_argument(
        "--formats",
        nargs="+",
        required=True,
        metavar="FORMAT",
        help="Fingerprint formats: bin/binary and cnt/count.",
    )
    parser.add_argument(
        "--levels-growth",
        nargs="+",
        required=True,
        metavar="LEVEL:GROWTH",
        help=(
            'Pairs such as 2:10 3:5 6:2 or '
            '"[(2,10), (3,5), (6,2)]".'
        ),
    )
    parser.add_argument(
        "--config-dir",
        type=Path,
        default=Path("generated_configs"),
        metavar="PATH",
        help="Directory receiving JSON files and the run summary.",
    )
    parser.add_argument(
        "--workdir-root",
        type=Path,
        metavar="PATH",
        help=(
            "Root for per-combination workdirs. The base config workdir is "
            "used when this option is omitted."
        ),
    )
    parser.add_argument(
        "--log-dir",
        type=Path,
        metavar="PATH",
        help="Log directory; defaults to CONFIG_DIR/logs.",
    )
    parser.add_argument(
        "--terminal-executable",
        type=Path,
        metavar="PATH",
        help=(
            "Terminal executable or Python script. When omitted, resolve the "
            "installed command, a sibling executable, or the source script."
        ),
    )
    parser.add_argument(
        "--continue-on-error",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="Continue with later configs after a failed child process.",
    )
    parser.add_argument(
        "--preserve-output-paths",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="Keep explicit output paths inherited from the base config.",
    )
    return parser


def _read_json_object(path: Path) -> dict[str, Any]:
    path = path.expanduser().resolve()
    if path.suffix.lower() != ".json":
        raise BatchConfigurationError(
            f"The base configuration must be a .json file: {path}"
        )
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except FileNotFoundError as exc:
        raise BatchConfigurationError(f"Base configuration not found: {path}") from exc
    except json.JSONDecodeError as exc:
        raise BatchConfigurationError(
            f"Invalid JSON in {path} at line {exc.lineno}, column {exc.colno}: "
            f"{exc.msg}"
        ) from exc
    if not isinstance(payload, dict):
        raise BatchConfigurationError(
            f"The base configuration root must be a JSON object: {path}"
        )
    return payload


def _literal_or_tokens(values: Sequence[str]) -> list[Any]:
    text = " ".join(values).strip()
    try:
        parsed = ast.literal_eval(text)
    except (SyntaxError, ValueError):
        return list(values)
    if isinstance(parsed, (list, tuple, set)):
        return list(parsed)
    return [parsed]


def _deduplicate(values: Iterable[Any]) -> list[Any]:
    unique: list[Any] = []
    for value in values:
        if value not in unique:
            unique.append(value)
    return unique


def parse_bits(values: Sequence[str]) -> list[int]:
    bits: list[int] = []
    for raw in _literal_or_tokens(values):
        if isinstance(raw, bool):
            raise BatchConfigurationError("Fingerprint lengths must be integers.")
        try:
            value = int(raw)
        except (TypeError, ValueError) as exc:
            raise BatchConfigurationError(
                f"Invalid fingerprint length: {raw!r}"
            ) from exc
        if value <= 0:
            raise BatchConfigurationError(
                f"Fingerprint lengths must be positive: {value}"
            )
        bits.append(value)
    if not bits:
        raise BatchConfigurationError("At least one fingerprint length is required.")
    return _deduplicate(bits)


def parse_formats(values: Sequence[str]) -> list[str]:
    formats: list[str] = []
    for raw in _literal_or_tokens(values):
        normalized = str(raw).strip().lower()
        try:
            formats.append(_FORMAT_ALIASES[normalized])
        except KeyError as exc:
            raise BatchConfigurationError(
                f"Invalid format {raw!r}; use bin/binary or cnt/count."
            ) from exc
    if not formats:
        raise BatchConfigurationError("At least one fingerprint format is required.")
    return _deduplicate(formats)


def parse_levels_growth(values: Sequence[str]) -> list[tuple[int, float]]:
    text = " ".join(values).strip()
    try:
        parsed = ast.literal_eval(text)
    except (SyntaxError, ValueError):
        parsed_items: list[Any] = list(values)
    else:
        if (
            isinstance(parsed, tuple)
            and len(parsed) == 2
            and not isinstance(parsed[0], (list, tuple))
        ):
            parsed_items = [parsed]
        elif isinstance(parsed, (list, tuple)):
            parsed_items = list(parsed)
        else:
            parsed_items = [parsed]

    pairs: list[tuple[int, float]] = []
    for raw in parsed_items:
        if isinstance(raw, str):
            separator = ":" if ":" in raw else ","
            components = [item.strip() for item in raw.strip("()[]").split(separator)]
        elif isinstance(raw, (list, tuple)) and len(raw) == 2:
            components = list(raw)
        else:
            raise BatchConfigurationError(
                f"Invalid level/growth pair: {raw!r}; expected LEVEL:GROWTH."
            )
        if len(components) != 2:
            raise BatchConfigurationError(
                f"Invalid level/growth pair: {raw!r}; expected two values."
            )
        try:
            level = int(components[0])
            growth_ratio = float(components[1])
        except (TypeError, ValueError) as exc:
            raise BatchConfigurationError(
                f"Invalid numeric level/growth pair: {raw!r}"
            ) from exc
        if level <= 0 or growth_ratio <= 0:
            raise BatchConfigurationError(
                f"Level and growth ratio must be positive: {raw!r}"
            )
        pairs.append((level, growth_ratio))
    if not pairs:
        raise BatchConfigurationError(
            "At least one (level, growth_ratio) pair is required."
        )
    return _deduplicate(pairs)


def _project_section(document: dict[str, Any]) -> dict[str, Any]:
    project = document.get("project")
    if project is None:
        return document
    if not isinstance(project, dict):
        raise BatchConfigurationError(
            "The optional 'project' value in the base config must be an object."
        )
    return project


def _terminal_section(document: dict[str, Any]) -> dict[str, Any]:
    terminal = document.get("terminal")
    if terminal is None:
        return document
    if not isinstance(terminal, dict):
        raise BatchConfigurationError(
            "The optional 'terminal' value in the base config must be an object."
        )
    return terminal


def _configured_ligands(value: Any) -> list[str]:
    if value is None:
        return []
    if isinstance(value, str):
        values = [item.strip() for item in value.split(",") if item.strip()]
    elif isinstance(value, (list, tuple, set)):
        values = [str(item).strip() for item in value if str(item).strip()]
    else:
        values = [str(value).strip()]
    if len(values) == 1 and values[0].upper() == "ALL":
        return []
    return values


def prepare_static_inputs(base_config: dict[str, Any]) -> dict[str, Any]:
    """Materialize shared ligands and scheduler CPU settings exactly once."""

    prepared = copy.deepcopy(base_config)
    project = _project_section(prepared)
    terminal = _terminal_section(prepared)
    selected = _configured_ligands(project.get("selected_ligands"))
    entries_file = terminal.get("entries_file") or terminal.get(
        "selected_ligands_file"
    )

    if entries_file:
        entries_path = Path(str(entries_file)).expanduser()
        try:
            selected = [
                line.strip()
                for line in entries_path.read_text(
                    encoding="utf-8", errors="replace"
                ).splitlines()
                if line.strip() and not line.lstrip().startswith("#")
            ]
        except OSError as exc:
            raise BatchConfigurationError(
                f"Could not read the shared entries file {entries_path}: {exc}"
            ) from exc
    elif not selected:
        ligand_file = str(project.get("ligand_file") or "").strip()
        ligand_path = Path(ligand_file).expanduser() if ligand_file else None
        if ligand_path is not None and ligand_path.exists():
            try:
                selected = ligand_io.parse_ligand_file(ligand_file)
            except (OSError, ValueError) as exc:
                raise BatchConfigurationError(
                    f"Could not parse the shared ligand source {ligand_path}: {exc}"
                ) from exc

    if selected:
        project["selected_ligands"] = selected
        terminal.pop("entries_file", None)
        terminal.pop("selected_ligands_file", None)

    allocation = detect_cpu_allocation()
    project["nproc"] = effective_nproc(project.get("nproc"))
    print(
        f"[multiple-run] nproc={project['nproc']} "
        f"from {allocation.source}"
    )
    return prepared


def _slug_number(value: float) -> str:
    if value.is_integer():
        return str(int(value))
    return format(value, "g").replace(".", "p")


def _write_json_atomic(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary_path: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            newline="",
            dir=path.parent,
            prefix=f".{path.name}.",
            suffix=".tmp",
            delete=False,
        ) as temporary:
            temporary_path = Path(temporary.name)
            json.dump(payload, temporary, indent=2, sort_keys=True, ensure_ascii=False)
            temporary.write("\n")
            temporary.flush()
            os.fsync(temporary.fileno())
        for attempt in range(20):
            try:
                os.replace(temporary_path, path)
                break
            except PermissionError:
                if attempt == 19:
                    raise
                # Windows can briefly lock the destination while another
                # orchestrator atomically replaces the same summary.
                time.sleep(0.01)
    finally:
        if temporary_path is not None and temporary_path.exists():
            temporary_path.unlink()


def _config_digest(payload: dict[str, Any]) -> str:
    canonical = json.dumps(
        payload,
        ensure_ascii=False,
        separators=(",", ":"),
        sort_keys=True,
    ).encode("utf-8")
    return hashlib.sha256(canonical).hexdigest()


def generate_configurations(
    *,
    base_config: dict[str, Any],
    bits: Sequence[int],
    formats: Sequence[str],
    levels_growth: Sequence[tuple[int, float]],
    config_dir: Path,
    workdir_root: Path | None = None,
    preserve_output_paths: bool = False,
) -> list[GeneratedRun]:
    """Write the full matrix before any child process is launched."""

    config_dir = config_dir.expanduser().resolve()
    base_project = _project_section(base_config)
    configured_workdir = str(base_project.get("workdir", "")).strip()
    if workdir_root is None:
        if not configured_workdir:
            raise BatchConfigurationError(
                "The base config must define workdir, or use --workdir-root."
            )
        workdir_root = Path(configured_workdir)
    workdir_root = workdir_root.expanduser().resolve()

    generated: list[GeneratedRun] = []
    for bit_count, fingerprint_format, pair in itertools.product(
        bits, formats, levels_growth
    ):
        level, growth_ratio = pair
        growth_slug = _slug_number(float(growth_ratio))
        run_id = f"B{bit_count}_L{level}_G{growth_slug}_{fingerprint_format}"
        document = copy.deepcopy(base_config)
        project = _project_section(document)
        project["ifp_length"] = int(bit_count)
        project["ifp_levels"] = int(level)
        project["ifp_radius"] = float(growth_ratio)
        project["ifp_bit"] = fingerprint_format == "bin"
        project["workdir"] = str(workdir_root / run_id)

        if not preserve_output_paths:
            project["ifp_output"] = ""
            project["sim_matrix_output"] = ""
            project["pse_path"] = ""

        config_path = config_dir / f"hipplinteractomics_{run_id}.json"
        _write_json_atomic(config_path, document)
        generated.append(
            GeneratedRun(
                run_id=run_id,
                bits=int(bit_count),
                fingerprint_format=fingerprint_format,
                level=int(level),
                growth_ratio=float(growth_ratio),
                config_path=config_path,
                config_digest=_config_digest(document),
            )
        )

    if not generated:
        raise BatchConfigurationError("The Cartesian product produced no configs.")
    return generated


def resolve_terminal_command(explicit: Path | None = None) -> list[str]:
    """Resolve an installed command, packaged executable, or source script."""

    if explicit is not None:
        candidate = explicit.expanduser().resolve()
        if not candidate.is_file():
            raise BatchConfigurationError(
                f"Terminal executable was not found: {candidate}"
            )
        if candidate.suffix.lower() == ".py":
            return [sys.executable, str(candidate)]
        return [str(candidate)]

    installed = shutil.which("hipplinteractomics-terminal")
    if installed:
        return [installed]

    source_script = Path(__file__).resolve().with_name(
        "hipplinteractomics_terminal.py"
    )
    for name in (
        "hipplinteractomics-terminal.exe",
        "hipplinteractomics-terminal",
    ):
        sibling = source_script.with_name(name)
        if sibling.is_file():
            return [str(sibling)]
    if source_script.is_file():
        return [sys.executable, str(source_script)]
    raise BatchConfigurationError(
        "Could not resolve hipplinteractomics-terminal. Install the package "
        "or pass --terminal-executable."
    )


def _format_command(command: Sequence[str]) -> str:
    if os.name == "nt":
        return subprocess.list2cmdline(list(command))
    return shlex.join(command)


def execute_serially(
    generated: Sequence[GeneratedRun],
    *,
    terminal_command: Sequence[str],
    log_dir: Path,
    continue_on_error: bool = False,
    on_result: Callable[[RunResult], None] | None = None,
    termination: TerminationController | None = None,
) -> list[RunResult]:
    """Run configs synchronously and capture combined output per run."""

    log_dir = log_dir.expanduser().resolve()
    log_dir.mkdir(parents=True, exist_ok=True)
    results: list[RunResult] = []
    termination = termination or TerminationController()
    total = len(generated)

    child_environment = os.environ.copy()
    child_environment.setdefault("MPLBACKEND", "Agg")
    child_environment.setdefault("QT_QPA_PLATFORM", "offscreen")
    child_environment.setdefault("PYTHONUNBUFFERED", "1")

    for position, generated_run in enumerate(generated, start=1):
        command = [*terminal_command, str(generated_run.config_path)]
        log_path = log_dir / f"{generated_run.run_id}.log"
        started = time.monotonic()
        print(
            f"[multiple-run] [{position}/{total}] starting "
            f"{generated_run.run_id}",
        )
        print(f"[multiple-run] command: {_format_command(command)}")

        try:
            with log_path.open(
                "w", encoding="utf-8", newline="", buffering=64 * 1024
            ) as log_handle:
                log_handle.write(f"$ {_format_command(command)}\n\n")
                process = subprocess.Popen(
                    command,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT,
                    text=True,
                    encoding="utf-8",
                    errors="replace",
                    env=child_environment,
                    bufsize=1,
                )
                termination.attach(process)
                try:
                    if process.stdout is None:
                        raise BatchExecutionError(
                            f"No output stream for {generated_run.run_id}."
                        )
                    for line in process.stdout:
                        print(line, end="")
                        log_handle.write(line)
                    returncode = process.wait()
                finally:
                    if process.stdout is not None:
                        process.stdout.close()
                    termination.detach(process)
        except OSError as exc:
            try:
                log_path.write_text(
                    f"Failed to start child process: {exc}\n",
                    encoding="utf-8",
                )
            except OSError:
                pass
            returncode = 127

        duration = time.monotonic() - started
        status = (
            "interrupted"
            if termination.received_signal is not None
            else "completed" if returncode == 0 else "failed"
        )
        result = RunResult(
            run_id=generated_run.run_id,
            config_path=str(generated_run.config_path),
            log_path=str(log_path),
            command=list(command),
            returncode=returncode,
            duration_seconds=round(duration, 3),
            config_digest=generated_run.config_digest,
            status=status,
        )
        results.append(result)
        if on_result is not None:
            on_result(result)

        if termination.received_signal is not None:
            print(
                f"[multiple-run] interrupted by signal "
                f"{termination.received_signal}",
                file=sys.stderr,
            )
            break

        if returncode == 0:
            print(
                f"[multiple-run] completed {generated_run.run_id} "
                f"in {duration:.1f}s",
            )
            continue

        message = (
            f"{generated_run.run_id} failed with exit code {returncode}. "
            f"See {log_path}"
        )
        print(f"[multiple-run] ERROR: {message}", file=sys.stderr)
        if not continue_on_error:
            break

    return results


def _summary_payload(
    generated: Sequence[GeneratedRun],
    results: Sequence[RunResult],
    *,
    status: str,
    error: str | None = None,
) -> dict[str, Any]:
    result_by_id = {result.run_id: result for result in results}
    runs: list[dict[str, Any]] = []
    for item in generated:
        row: dict[str, Any] = {
            "run_id": item.run_id,
            "bits": item.bits,
            "format": item.fingerprint_format,
            "level": item.level,
            "growth_ratio": item.growth_ratio,
            "config_path": str(item.config_path),
            "config_digest": item.config_digest,
            "status": "pending",
        }
        result = result_by_id.get(item.run_id)
        if result is not None:
            row.update(asdict(result))
        runs.append(row)
    payload: dict[str, Any] = {"status": status, "runs": runs}
    if error:
        payload["error"] = error
    return payload


def _load_completed_results(
    summary_path: Path,
    generated: Sequence[GeneratedRun],
) -> dict[str, RunResult]:
    if not summary_path.is_file():
        return {}
    try:
        payload = json.loads(summary_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return {}
    rows = payload.get("runs") if isinstance(payload, dict) else None
    if not isinstance(rows, list):
        return {}

    expected = {item.run_id: item for item in generated}
    completed: dict[str, RunResult] = {}
    for row in rows:
        if not isinstance(row, dict) or row.get("status") != "completed":
            continue
        run_id = str(row.get("run_id") or "")
        item = expected.get(run_id)
        if item is None or row.get("config_digest") != item.config_digest:
            continue
        try:
            completed[run_id] = RunResult(
                run_id=run_id,
                config_path=str(row.get("config_path") or item.config_path),
                log_path=str(row.get("log_path") or ""),
                command=[str(value) for value in row.get("command", [])],
                returncode=0,
                duration_seconds=float(row.get("duration_seconds") or 0.0),
                config_digest=item.config_digest,
                status="completed",
            )
        except (TypeError, ValueError):
            continue
    return completed


def main(argv: Sequence[str] | None = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)

    generated: list[GeneratedRun] = []
    results: list[RunResult] = []
    config_dir = args.config_dir.expanduser().resolve()
    summary_path = config_dir / "pipeline_summary.json"

    try:
        base_config = prepare_static_inputs(_read_json_object(args.base_config))
        bits = parse_bits(args.bits)
        formats = parse_formats(args.formats)
        levels_growth = parse_levels_growth(args.levels_growth)

        generated = generate_configurations(
            base_config=base_config,
            bits=bits,
            formats=formats,
            levels_growth=levels_growth,
            config_dir=config_dir,
            workdir_root=args.workdir_root,
            preserve_output_paths=args.preserve_output_paths,
        )
        print(
            f"[multiple-run] generated {len(generated)} JSON configs in "
            f"{config_dir}",
        )

        terminal_command = resolve_terminal_command(args.terminal_executable)
        log_dir = (
            args.log_dir.expanduser().resolve()
            if args.log_dir
            else config_dir / "logs"
        )
        result_by_id = _load_completed_results(summary_path, generated)
        if result_by_id:
            print(
                f"[multiple-run] resume: skipping {len(result_by_id)} "
                "completed matching configuration(s)"
            )

        def persist(status: str, error: str | None = None) -> None:
            _write_json_atomic(
                summary_path,
                _summary_payload(
                    generated,
                    list(result_by_id.values()),
                    status=status,
                    error=error,
                ),
            )

        persist("running")
        pending = [item for item in generated if item.run_id not in result_by_id]
        termination = TerminationController()

        def checkpoint(result: RunResult) -> None:
            result_by_id[result.run_id] = result
            persist("interrupted" if result.status == "interrupted" else "running")

        with termination.installed():
            new_results = execute_serially(
                pending,
                terminal_command=terminal_command,
                log_dir=log_dir,
                continue_on_error=args.continue_on_error,
                on_result=checkpoint,
                termination=termination,
            )
        results = list(result_by_id.values())
        if termination.received_signal is not None:
            persist("interrupted")
            return signal_exit_code(termination.received_signal)

        failed = [result for result in new_results if result.returncode != 0]
        all_completed = len(result_by_id) == len(generated) and not failed
        status = "completed" if all_completed else "failed"
        _write_json_atomic(
            summary_path,
            _summary_payload(generated, results, status=status),
        )
        if failed:
            print(
                f"[multiple-run] {len(failed)} run(s) failed; summary: "
                f"{summary_path}",
                file=sys.stderr,
            )
            return 1
        if all_completed:
            print(f"[multiple-run] all runs completed; summary: {summary_path}")
            return 0
        print(
            f"[multiple-run] batch stopped before all configurations completed; "
            f"summary: {summary_path}",
            file=sys.stderr,
        )
        return 1

    except (BatchConfigurationError, BatchExecutionError) as exc:
        if generated:
            _write_json_atomic(
                summary_path,
                _summary_payload(
                    generated,
                    results,
                    status="failed",
                    error=str(exc),
                ),
            )
        print(f"[multiple-run] ERROR: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
