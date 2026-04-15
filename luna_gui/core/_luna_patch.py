"""Post-install patch for LUNA — run INSIDE the luna-env.

Fixes the Windows-specific `np.int_` overflow bug in LUNA's fingerprint
code. On Windows, `np.int_` maps to a 32-bit C long and overflows with
the 64-bit hashes mmh3 produces. Rewrites it to np.int64.

Idempotent — safe to run multiple times.
"""
from __future__ import annotations

import os
import sys


def patch_fingerprint() -> int:
    try:
        import luna  # type: ignore
    except ImportError:
        print("[patch] luna not importable — aborting", file=sys.stderr)
        return 1
    path = os.path.join(
        os.path.dirname(luna.__file__), "interaction", "fp", "fingerprint.py"
    )
    if not os.path.exists(path):
        print(f"[patch] file not found: {path}", file=sys.stderr)
        return 1
    with open(path, "r", encoding="utf-8") as f:
        src = f.read()
    count = src.count("np.int_")
    if count == 0:
        print("[patch] fingerprint.py: already patched, skipping")
        return 0
    patched = src.replace("np.int_", "np.int64")
    with open(path, "w", encoding="utf-8") as f:
        f.write(patched)
    print(f"[patch] fingerprint.py: np.int_ -> np.int64 ({count} replacements)")
    return 0


if __name__ == "__main__":
    sys.exit(patch_fingerprint())
