from __future__ import annotations

import ast
import json
from pathlib import Path

from luna_gui.i18n import language, language_notifier, set_language, t


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


def test_locale_json_catalogs_have_identical_keys_and_english_fallback() -> None:
    catalogs = {
        code: json.loads((ROOT / "luna_gui" / "locales" / f"{code}.json").read_text(encoding="utf-8"))
        for code in ("en", "pt", "es")
    }
    assert set(catalogs["en"]) == set(catalogs["pt"]) == set(catalogs["es"])
    set_language("unsupported")
    assert language() == "en"
    assert t("Carregar resultados") == "Load results"


def test_literal_widget_text_is_registered_in_locale_catalogs() -> None:
    catalog = json.loads(
        (ROOT / "luna_gui" / "locales" / "en.json").read_text(encoding="utf-8")
    )
    constructors = {"QLabel", "QPushButton", "QGroupBox", "QAction"}
    methods = {"setTitle", "setPlaceholderText"}
    missing: set[str] = set()
    for path in (ROOT / "luna_gui" / "ui").glob("*.py"):
        tree = ast.parse(path.read_text(encoding="utf-8"))
        for node in ast.walk(tree):
            if not isinstance(node, ast.Call):
                continue
            name = node.func.id if isinstance(node.func, ast.Name) else ""
            attribute = node.func.attr if isinstance(node.func, ast.Attribute) else ""
            if name in constructors and node.args:
                candidate = node.args[0]
            elif attribute in methods and node.args:
                candidate = node.args[0]
            elif attribute in {"addTab", "insertTab"} and node.args:
                candidate = node.args[-1]
            else:
                continue
            if isinstance(candidate, ast.Constant) and isinstance(candidate.value, str):
                text = candidate.value.strip()
                if text and any(char.isalpha() for char in text) and text not in catalog:
                    missing.add(text)
    assert missing == set(), json.dumps(sorted(missing), ensure_ascii=False, indent=2)


def test_language_change_emits_qt_signal() -> None:
    notifier = language_notifier()
    if notifier is None:
        return
    received: list[str] = []
    notifier.language_changed.connect(received.append)
    set_language("en")
    set_language("pt")
    set_language("en")
    assert received[-2:] == ["pt", "en"]
