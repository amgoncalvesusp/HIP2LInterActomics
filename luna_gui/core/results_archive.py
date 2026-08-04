"""Safe extraction and workdir discovery for portable Results archives."""
from __future__ import annotations

import shutil
import stat
import tarfile
import zipfile
from pathlib import Path, PurePosixPath


class ResultsArchiveError(RuntimeError):
    pass


def _safe_relative_path(name: str) -> Path:
    normalized = str(name or "").replace("\\", "/")
    candidate = PurePosixPath(normalized)
    if not normalized or candidate.is_absolute() or ".." in candidate.parts:
        raise ResultsArchiveError(f"Unsafe archive path: {name!r}")
    if candidate.parts and ":" in candidate.parts[0]:
        raise ResultsArchiveError(f"Unsafe archive drive path: {name!r}")
    return Path(*candidate.parts)


def _copy_stream(source, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open("wb") as target:
        shutil.copyfileobj(source, target, length=1024 * 1024)


def _extract_zip(path: Path, destination: Path) -> None:
    with zipfile.ZipFile(path) as archive:
        for member in archive.infolist():
            relative = _safe_relative_path(member.filename)
            mode = member.external_attr >> 16
            if stat.S_ISLNK(mode):
                raise ResultsArchiveError(f"Archive links are not allowed: {member.filename}")
            target = destination / relative
            if member.is_dir():
                target.mkdir(parents=True, exist_ok=True)
                continue
            with archive.open(member, "r") as source:
                _copy_stream(source, target)


def _extract_tar(path: Path, destination: Path) -> None:
    with tarfile.open(path, mode="r:*") as archive:
        for member in archive.getmembers():
            relative = _safe_relative_path(member.name)
            if member.issym() or member.islnk() or member.isdev():
                raise ResultsArchiveError(f"Archive links/devices are not allowed: {member.name}")
            target = destination / relative
            if member.isdir():
                target.mkdir(parents=True, exist_ok=True)
                continue
            if not member.isfile():
                continue
            source = archive.extractfile(member)
            if source is None:
                raise ResultsArchiveError(f"Could not read archive member: {member.name}")
            with source:
                _copy_stream(source, target)


def find_results_workdir(root: str | Path) -> Path:
    base = Path(root)
    project_files = sorted(base.rglob(".luna_gui.json"), key=lambda item: (len(item.parts), str(item)))
    if project_files:
        return project_files[0].parent
    candidates = [base]
    candidates.extend(path.parent for path in base.rglob("results") if path.is_dir())
    for candidate in sorted(set(candidates), key=lambda item: (len(item.parts), str(item))):
        if (candidate / "results").is_dir():
            return candidate
    raise ResultsArchiveError("The archive does not contain a HIP2LInterActomics results workdir.")


def extract_results_archive(archive_path: str | Path, destination: str | Path) -> Path:
    source = Path(archive_path)
    output = Path(destination)
    if not source.is_file():
        raise ResultsArchiveError(f"Archive not found: {source}")
    output.mkdir(parents=True, exist_ok=False)
    try:
        suffixes = [suffix.casefold() for suffix in source.suffixes]
        if source.suffix.casefold() == ".zip":
            _extract_zip(source, output)
        elif suffixes[-2:] in ([".tar", ".gz"], [".tar", ".bz2"], [".tar", ".xz"]) or source.suffix.casefold() in {".tgz", ".tbz2", ".txz", ".tar"}:
            _extract_tar(source, output)
        else:
            raise ResultsArchiveError("Supported formats: .zip, .tar, .tar.gz, .tgz, .tar.bz2, .tar.xz")
        return find_results_workdir(output)
    except BaseException:
        shutil.rmtree(output, ignore_errors=True)
        raise
