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


def test_critical_workflow_strings_translate_in_all_three_languages() -> None:
    expected = {
        "&Aparência": ("&Appearance", "&Apariencia"),
        "Fechar aplicativo": ("Close application", "Cerrar aplicación"),
        "Usar bit fingerprints (padrão: count)": (
            "Use bit fingerprints (default: count)",
            "Usar bit fingerprints (predeterminado: count)",
        ),
        "Seed da importância:": ("Importance seed:", "Semilla de importancia:"),
        "Selecionar/desselecionar todas": ("Select/deselect all", "Seleccionar/deseleccionar todas"),
        "add_atom_atom (interações atômicas genéricas)": (
            "add_atom_atom (generic atom-atom interactions)",
            "add_atom_atom (interacciones atómicas genéricas)",
        ),
        "Gerar relatório PDF": ("Generate PDF report", "Generar reporte PDF"),
        "plots generation": ("plots generation", "generación de gráficos"),
        "All Frames": ("All Frames", "Todos los frames"),
    }
    for source, (english, spanish) in expected.items():
        assert t(source, lang="pt") in {source, "geração de gráficos", "Todos os frames"}
        assert t(source, lang="en") == english
        assert t(source, lang="es") == spanish


def test_dynamic_setup_and_count_patterns_are_translated() -> None:
    assert t("Prefixo do env: C:/Users/test/luna-env", lang="en") == "Environment prefix: C:/Users/test/luna-env"
    assert t("Python do env: /opt/luna/bin/python", lang="es") == "Python del entorno: /opt/luna/bin/python"
    assert t("8 de 10 ligantes selecionados", lang="en") == "8 of 10 selected ligands"
    assert t("8 de 10 ligantes selecionados", lang="es") == "8 de 10 ligandos seleccionados"


def test_fp_method_description_templates_are_translatable() -> None:
    source = "Features confiáveis por classe: {count}/{total}"
    assert t(source, lang="en").format(count=8, total=10) == "Class-reliable features: 8/10"
    assert t(source, lang="es").format(count=8, total=10) == "Features confiables por clase: 8/10"


def test_every_information_circle_tooltip_has_english_and_spanish_text() -> None:
    for path in (ROOT / "luna_gui" / "ui").glob("*.py"):
        tree = ast.parse(path.read_text(encoding="utf-8"))
        for node in ast.walk(tree):
            if not (
                isinstance(node, ast.Call)
                and isinstance(node.func, ast.Name)
                and node.func.id == "InfoButton"
                and node.args
                and isinstance(node.args[0], ast.Constant)
                and isinstance(node.args[0].value, str)
            ):
                continue
            source = node.args[0].value
            assert t(source, lang="en") != source, f"Missing English InfoButton tooltip: {source}"
            assert t(source, lang="es") != source, f"Missing Spanish InfoButton tooltip: {source}"
