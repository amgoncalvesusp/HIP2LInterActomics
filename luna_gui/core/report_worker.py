"""Dedicated PDF renderer entry point used by the Qt application."""
from __future__ import annotations

import sys

from .report_export import execute_pdf_render_job


def main(argv: list[str] | None = None) -> int:
    args = list(sys.argv[1:] if argv is None else argv)
    if len(args) != 2:
        print("usage: python -m luna_gui.core.report_worker JOB_PATH STATUS_PATH", file=sys.stderr)
        return 2
    return execute_pdf_render_job(args[0], args[1])


if __name__ == "__main__":
    raise SystemExit(main())
