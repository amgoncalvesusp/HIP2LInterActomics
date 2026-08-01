from __future__ import annotations

import ast
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_translation_catalog_has_no_duplicate_literal_keys() -> None:
    source = (ROOT / "luna_gui" / "i18n.py").read_text(encoding="utf-8")
    tree = ast.parse(source)
    catalog = next(
        node.value
        for node in tree.body
        if isinstance(node, ast.AnnAssign)
        and isinstance(node.target, ast.Name)
        and node.target.id == "TRANSLATIONS"
    )
    assert isinstance(catalog, ast.Dict)

    seen: set[str] = set()
    duplicates: set[str] = set()
    for key in catalog.keys:
        if isinstance(key, ast.Constant) and isinstance(key.value, str):
            if key.value in seen:
                duplicates.add(key.value)
            seen.add(key.value)

    assert duplicates == set()


def test_windowed_entry_point_installs_exception_hook() -> None:
    source = (ROOT / "luna_gui" / "main.py").read_text(encoding="utf-8")
    assert "sys.excepthook = _exception_hook" in source
    assert "QMessageBox.critical" in source
