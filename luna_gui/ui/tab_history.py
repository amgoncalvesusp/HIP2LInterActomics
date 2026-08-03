"""History tab — list past projects (workdirs) and reload them."""
from __future__ import annotations

from pathlib import Path

from PyQt6.QtCore import Qt, pyqtSignal
from PyQt6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QPushButton, QListWidget,
    QListWidgetItem, QMessageBox, QLabel,
)

from ..core.project import (
    ProjectConfig, clear_history, load_history, remove_from_history, PROJECT_FILENAME,
)


class HistoryTab(QWidget):
    project_loaded = pyqtSignal(object)  # ProjectConfig

    def __init__(self) -> None:
        super().__init__()
        layout = QVBoxLayout(self)
        layout.addWidget(QLabel("Projetos LUNA recentes (workdirs salvos):"))
        self.list = QListWidget()
        self.list.itemDoubleClicked.connect(lambda _: self.load_selected())
        layout.addWidget(self.list, 1)

        row = QHBoxLayout()
        btn_reload = QPushButton("Atualizar")
        btn_reload.clicked.connect(self.refresh)
        btn_load = QPushButton("Carregar projeto selecionado")
        btn_load.clicked.connect(self.load_selected)
        btn_remove = QPushButton("Remover da lista")
        btn_remove.clicked.connect(self.remove_selected)
        btn_clear = QPushButton("Remover toda a lista")
        btn_clear.setToolTip("Remove todos os projetos do histórico sem apagar seus arquivos.")
        btn_clear.clicked.connect(self.remove_all)
        row.addWidget(btn_reload)
        row.addWidget(btn_load)
        row.addWidget(btn_remove)
        row.addWidget(btn_clear)
        row.addStretch()
        layout.addLayout(row)

        self.refresh()

    def refresh(self) -> None:
        self.list.clear()
        for wd in load_history():
            exists = Path(wd).exists()
            label = wd if exists else f"{wd}  (não encontrado)"
            it = QListWidgetItem(label)
            it.setData(Qt.ItemDataRole.UserRole, wd)
            if not exists:
                it.setForeground(Qt.GlobalColor.gray)
            self.list.addItem(it)

    def load_selected(self) -> None:
        it = self.list.currentItem()
        if not it:
            return
        wd = it.data(Qt.ItemDataRole.UserRole)
        cfg_file = Path(wd) / PROJECT_FILENAME
        if not cfg_file.exists():
            QMessageBox.warning(self, "Sem config",
                                f"Nenhum {PROJECT_FILENAME} encontrado em {wd}")
            return
        try:
            cfg = ProjectConfig.load(cfg_file)
        except Exception as e:
            QMessageBox.critical(self, "Erro ao carregar", str(e))
            return
        self.project_loaded.emit(cfg)

    def remove_selected(self) -> None:
        it = self.list.currentItem()
        if not it:
            return
        wd = it.data(Qt.ItemDataRole.UserRole)
        remove_from_history(str(wd))
        self.refresh()

    def remove_all(self) -> None:
        if not load_history():
            return
        answer = QMessageBox.question(
            self,
            "Limpar histórico",
            "Remover todos os projetos da lista? Os arquivos dos projetos serão preservados.",
        )
        if answer != QMessageBox.StandardButton.Yes:
            return
        clear_history()
        self.refresh()
