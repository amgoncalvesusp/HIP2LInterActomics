"""Visual editor for LUNA's binding_modes.cfg.

The .cfg has one section per interaction type with two keys:
  accept_all = True/False
  accept_only = ["chain/resname/resnum/atom", ...]

A wildcard section [*] applies to all unlisted interaction types.
Selectors use "*" as wildcard for any field.
"""
from __future__ import annotations

import ast
from pathlib import Path

from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import (
    QDialog, QVBoxLayout, QHBoxLayout, QPushButton, QTableWidget,
    QTableWidgetItem, QHeaderView, QFileDialog, QMessageBox, QLabel,
    QComboBox, QInputDialog,
)


# Standard LUNA interaction types (from manual + InteractionCalculator).
DEFAULT_TYPES = [
    "Hydrogen bond",
    "Weak hydrogen bond",
    "Hydrophobic",
    "Aromatic stacking",
    "Edge-to-face pi-stacking",
    "Face-to-edge pi-stacking",
    "Parallel-displaced pi-stacking",
    "Cation-pi",
    "Amide-pi stacking",
    "Halogen bond",
    "Chalcogen bond",
    "Ionic",
    "Salt bridge",
    "Multipolar",
    "Ion-multipole",
    "Repulsive",
    "Van der Waals",
    "Van der Waals clash",
    "Proximal",
    "Covalent bond",
    "*",
]


class BindingModesEditor(QDialog):
    def __init__(self, parent=None, initial_path: str = "") -> None:
        super().__init__(parent)
        self.setWindowTitle("Editor de Binding Modes")
        self.resize(820, 520)
        self.path: str = initial_path

        layout = QVBoxLayout(self)

        layout.addWidget(QLabel(
            "Seletores no formato <CHAIN>/<RES>/<NUM>/<ATOM>. Use * como curinga.\n"
            "Ex.: */HIS/*/* (todos os átomos de histidinas), A/LYS/245/N* (nitrogênios da Lys245 cadeia A)."
        ))

        self.table = QTableWidget(0, 3)
        self.table.setHorizontalHeaderLabels(["Tipo de interação", "accept_all", "accept_only (lista)"])
        self.table.horizontalHeader().setSectionResizeMode(0, QHeaderView.ResizeMode.Stretch)
        self.table.horizontalHeader().setSectionResizeMode(1, QHeaderView.ResizeMode.ResizeToContents)
        self.table.horizontalHeader().setSectionResizeMode(2, QHeaderView.ResizeMode.Stretch)
        layout.addWidget(self.table, 1)

        btn_row = QHBoxLayout()
        btn_add = QPushButton("+ Adicionar tipo")
        btn_add.clicked.connect(self._add_row_dialog)
        btn_del = QPushButton("− Remover linha")
        btn_del.clicked.connect(self._del_row)
        btn_load = QPushButton("Carregar .cfg")
        btn_load.clicked.connect(self._load)
        btn_save = QPushButton("Salvar .cfg")
        btn_save.clicked.connect(self._save)
        btn_ok = QPushButton("OK")
        btn_ok.clicked.connect(self._save_and_close)
        btn_cancel = QPushButton("Cancelar")
        btn_cancel.clicked.connect(self.reject)
        for b in (btn_add, btn_del, btn_load, btn_save):
            btn_row.addWidget(b)
        btn_row.addStretch()
        btn_row.addWidget(btn_ok)
        btn_row.addWidget(btn_cancel)
        layout.addLayout(btn_row)

        if initial_path and Path(initial_path).exists():
            self._load_from(initial_path)
        else:
            # seed with the most common types
            for t in ["Hydrogen bond", "Hydrophobic", "Aromatic stacking", "*"]:
                self._add_row(t, False, "")

    # ---- table helpers ----
    def _add_row(self, type_name: str, accept_all: bool, accept_only: str) -> None:
        r = self.table.rowCount()
        self.table.insertRow(r)
        self.table.setItem(r, 0, QTableWidgetItem(type_name))
        cb = QComboBox(); cb.addItems(["False", "True"])
        cb.setCurrentText(str(accept_all))
        self.table.setCellWidget(r, 1, cb)
        self.table.setItem(r, 2, QTableWidgetItem(accept_only))

    def _add_row_dialog(self) -> None:
        type_name, ok = QInputDialog.getItem(
            self, "Tipo de interação", "Selecione ou digite:", DEFAULT_TYPES, 0, True
        )
        if ok and type_name:
            self._add_row(type_name, False, "")

    def _del_row(self) -> None:
        r = self.table.currentRow()
        if r >= 0:
            self.table.removeRow(r)

    # ---- file IO ----
    def _load(self) -> None:
        f, _ = QFileDialog.getOpenFileName(self, "Abrir binding_modes.cfg", "",
                                           "Config (*.cfg);;Todos (*)")
        if f:
            self._load_from(f)

    def _load_from(self, path: str) -> None:
        self.path = path
        self.table.setRowCount(0)
        current_section: str | None = None
        accept_all = False
        accept_only = ""
        try:
            for line in Path(path).read_text(encoding="utf-8", errors="replace").splitlines():
                s = line.strip()
                if not s or s.startswith(";") or s.startswith("#"):
                    continue
                if s.startswith("[") and s.endswith("]"):
                    if current_section is not None:
                        self._add_row(current_section, accept_all, accept_only)
                    current_section = s[1:-1]
                    accept_all = False
                    accept_only = ""
                elif "=" in s:
                    k, v = s.split("=", 1)
                    k = k.strip().lower(); v = v.strip()
                    if k == "accept_all":
                        accept_all = v.lower() == "true"
                    elif k == "accept_only":
                        accept_only = v
            if current_section is not None:
                self._add_row(current_section, accept_all, accept_only)
        except Exception as e:
            QMessageBox.critical(self, "Erro ao ler", str(e))

    def _serialize(self) -> str:
        lines = ["; Generated by HIP²LInterActomics"]
        for r in range(self.table.rowCount()):
            name_item = self.table.item(r, 0)
            if not name_item or not name_item.text().strip():
                continue
            name = name_item.text().strip()
            cb: QComboBox = self.table.cellWidget(r, 1)  # type: ignore
            accept_all = cb.currentText() if cb else "False"
            accept_only = (self.table.item(r, 2).text() if self.table.item(r, 2) else "").strip()
            lines.append(f"\n[{name}]")
            lines.append(f"accept_all={accept_all}")
            if accept_only:
                # Allow user to type either a Python list or a comma-separated list
                if not accept_only.startswith("["):
                    parts = [p.strip().strip('"').strip("'") for p in accept_only.split(",") if p.strip()]
                    accept_only = "[" + ", ".join(f'"{p}"' for p in parts) + "]"
                else:
                    try:
                        ast.literal_eval(accept_only)
                    except Exception:
                        raise ValueError(f"accept_only inválido em [{name}]: {accept_only}")
                lines.append(f"accept_only={accept_only}")
        return "\n".join(lines) + "\n"

    def _save(self) -> str | None:
        try:
            content = self._serialize()
        except ValueError as e:
            QMessageBox.critical(self, "Erro de validação", str(e))
            return None
        path = self.path
        if not path:
            f, _ = QFileDialog.getSaveFileName(self, "Salvar binding_modes.cfg",
                                               "binding_modes.cfg", "Config (*.cfg)")
            if not f:
                return None
            path = f
        Path(path).write_text(content, encoding="utf-8")
        self.path = path
        return path

    def _save_and_close(self) -> None:
        if self._save():
            self.accept()
