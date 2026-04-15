"""Project tab — pick protein, ligand library, select ligands, choose workdir."""
from __future__ import annotations

from pathlib import Path

from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QFormLayout, QLabel, QLineEdit,
    QPushButton, QFileDialog, QListWidget, QListWidgetItem, QMessageBox,
    QGroupBox, QCheckBox,
)

from ..core.ligand_io import (
    parse_ligand_file, consolidate_folder, consolidate_folder_clean,
)
from ..core.project import ProjectConfig


class ProjectTab(QWidget):
    def __init__(self, cfg: ProjectConfig) -> None:
        super().__init__()
        self.cfg = cfg

        layout = QVBoxLayout(self)

        intro = QLabel(
            "Nesta aba você escolhe os arquivos de entrada do projeto. A proteína é o receptor em "
            "PDB, os ligantes são as moléculas a comparar, e o diretório de trabalho é onde a "
            "execução salvará entradas, resultados e relatórios."
        )
        intro.setWordWrap(True)
        intro.setProperty("muted", True)
        layout.addWidget(intro)

        # --- Inputs ---
        form_box = QGroupBox("Entradas")
        form = QFormLayout(form_box)

        self.protein_edit = QLineEdit()
        self.protein_edit.setToolTip("Arquivo PDB da proteína ou receptor que será analisado pelo LUNA.")
        btn_p = QPushButton("Procurar...")
        btn_p.setToolTip("Seleciona o arquivo PDB da proteína.")
        btn_p.clicked.connect(self._pick_protein)
        prow = QHBoxLayout(); prow.addWidget(self.protein_edit); prow.addWidget(btn_p)
        form.addRow("Proteína (PDB):", self._wrap(prow))

        self.ligand_edit = QLineEdit()
        self.ligand_edit.setToolTip(
            "Arquivo com os ligantes. Pode ser um MOL2/SDF único ou um arquivo consolidado com vários ligantes."
        )
        btn_l = QPushButton("Arquivo...")
        btn_l.setToolTip("Seleciona um arquivo de ligantes já pronto para leitura.")
        btn_l.clicked.connect(self._pick_ligand)
        btn_lf = QPushButton("Pasta MOL2...")
        btn_lf.setToolTip("Consolida todos os .mol2 da pasta em um único arquivo (modo robusto: remove LP e renumera)")
        btn_lf.clicked.connect(self._pick_ligand_folder)
        lrow = QHBoxLayout(); lrow.addWidget(self.ligand_edit); lrow.addWidget(btn_l); lrow.addWidget(btn_lf)
        form.addRow("Ligantes (MOL2/SDF):", self._wrap(lrow))

        self.workdir_edit = QLineEdit()
        self.workdir_edit.setToolTip(
            "Pasta do projeto. A GUI grava aqui o entries.txt, os resultados do LUNA e os arquivos exportados."
        )
        btn_w = QPushButton("Procurar...")
        btn_w.setToolTip("Escolhe a pasta onde o projeto será salvo e reaberto depois.")
        btn_w.clicked.connect(self._pick_workdir)
        wrow = QHBoxLayout(); wrow.addWidget(self.workdir_edit); wrow.addWidget(btn_w)
        form.addRow("Diretório de trabalho:", self._wrap(wrow))

        self.cb_waters = QCheckBox("Incluir águas (HOH) — análise hidratada")
        self.cb_waters.setToolTip(
            "Quando marcado, HOH é mantido no PDB e LUNA é chamado com "
            "ignore_any_h2o=False (inclui interações mediadas por água).\n"
            "Requer opções avançadas — usa o runner Python API."
        )
        self.cb_waters.toggled.connect(lambda v: setattr(self.cfg, "include_waters", v))
        form.addRow(self.cb_waters)

        prep_row = QHBoxLayout()
        btn_prep = QPushButton("Preparar arquivos de docking...")
        btn_prep.setToolTip("Divide MOL2 de docking (proteína+ligante) em PDB + MOL2 separados")
        btn_prep.clicked.connect(self._open_prep_wizard)
        prep_row.addWidget(btn_prep)
        prep_row.addStretch()
        form.addRow("Pré-processamento:", self._wrap(prep_row))

        layout.addWidget(form_box)

        # --- Ligand selection ---
        lig_box = QGroupBox("Ligantes detectados")
        lig_layout = QVBoxLayout(lig_box)

        lig_help = QLabel(
            "Depois de carregar os ligantes, marque apenas os que devem entrar na análise. "
            "O filtro ajuda a localizar nomes específicos sem apagar a lista."
        )
        lig_help.setWordWrap(True)
        lig_help.setProperty("muted", True)
        lig_layout.addWidget(lig_help)

        filter_row = QHBoxLayout()
        self.filter_edit = QLineEdit()
        self.filter_edit.setPlaceholderText("Filtrar por nome (texto livre)...")
        self.filter_edit.setToolTip("Mostra apenas os ligantes cujo nome contém o texto digitado.")
        self.filter_edit.textChanged.connect(self._apply_filter)
        btn_all = QPushButton("Selecionar tudo")
        btn_all.setToolTip("Marca todos os ligantes visíveis na lista atual.")
        btn_all.clicked.connect(lambda: self._set_all(True))
        btn_none = QPushButton("Nenhum")
        btn_none.setToolTip("Desmarca todos os ligantes visíveis na lista atual.")
        btn_none.clicked.connect(lambda: self._set_all(False))
        filter_row.addWidget(self.filter_edit, 1)
        filter_row.addWidget(btn_all)
        filter_row.addWidget(btn_none)
        lig_layout.addLayout(filter_row)

        self.lig_list = QListWidget()
        self.lig_list.setSelectionMode(QListWidget.SelectionMode.NoSelection)
        self.lig_list.setToolTip("Cada linha representa um ligante detectado no arquivo. Marque os que devem ser processados.")
        lig_layout.addWidget(self.lig_list, 1)

        self.count_label = QLabel("0 ligantes")
        self.count_label.setToolTip("Resumo de quantos ligantes estão marcados para entrar na execução.")
        lig_layout.addWidget(self.count_label)

        layout.addWidget(lig_box, 1)

    # ---- helpers ----
    def _wrap(self, sub_layout) -> QWidget:
        w = QWidget(); w.setLayout(sub_layout); return w

    def _pick_protein(self) -> None:
        f, _ = QFileDialog.getOpenFileName(self, "Selecionar proteína", "",
                                           "PDB (*.pdb *.ent);;Todos (*)")
        if f:
            self.protein_edit.setText(f)
            self.cfg.protein_file = f

    def _pick_ligand(self) -> None:
        f, _ = QFileDialog.getOpenFileName(self, "Selecionar ligantes", "",
                                           "Ligantes (*.mol2 *.sdf *.mol);;Todos (*)")
        if not f:
            return
        self.ligand_edit.setText(f)
        self.cfg.ligand_file = f
        self._load_ligands(f)

    def _pick_ligand_folder(self) -> None:
        d = QFileDialog.getExistingDirectory(self, "Pasta com arquivos .mol2")
        if not d:
            return
        out = Path(d) / "_consolidated_ligands.mol2"
        try:
            n, _names = consolidate_folder_clean(d, out, drop_lp=True)
        except Exception as e:
            # Fall back to the naïve concat if cleaning fails for any reason
            try:
                n = consolidate_folder(d, out)
            except Exception as e2:
                QMessageBox.critical(self, "Erro ao consolidar", f"{e}\n{e2}")
                return
        if n == 0:
            QMessageBox.warning(self, "Sem ligantes", "Nenhum arquivo .mol2 encontrado na pasta.")
            return
        self.ligand_edit.setText(str(out))
        self.cfg.ligand_file = str(out)
        self._load_ligands(str(out))
        QMessageBox.information(
            self, "Consolidado",
            f"{n} ligantes consolidados (LP removido, átomos renumerados) em\n{out}"
        )

    def _open_prep_wizard(self) -> None:
        from .dialog_prep import DockingPrepDialog
        dlg = DockingPrepDialog(self)
        dlg.exec()
        if dlg.cb_open_after.isChecked() and dlg.result_ligand_dir:
            # Consolidate the freshly generated ligands folder and load it
            out = Path(dlg.result_ligand_dir) / "_consolidated_ligands.mol2"
            try:
                n, _ = consolidate_folder_clean(dlg.result_ligand_dir, out, drop_lp=True)
            except Exception as e:
                QMessageBox.critical(self, "Erro ao consolidar", str(e))
                return
            if n > 0:
                self.ligand_edit.setText(str(out))
                self.cfg.ligand_file = str(out)
                self._load_ligands(str(out))
            # Point the user at the protein folder — they pick one
            if dlg.result_protein_dir:
                QMessageBox.information(
                    self, "Próximo passo",
                    "Agora selecione em 'Proteína (PDB)' um arquivo da pasta:\n"
                    + dlg.result_protein_dir,
                )

    def _pick_workdir(self) -> None:
        d = QFileDialog.getExistingDirectory(self, "Diretório de trabalho")
        if d:
            self.workdir_edit.setText(d)
            self.cfg.workdir = d

    def _load_ligands(self, path: str) -> None:
        try:
            names = parse_ligand_file(path)
        except Exception as e:
            QMessageBox.critical(self, "Erro ao ler ligantes", str(e))
            return
        self.lig_list.clear()
        for n in names:
            item = QListWidgetItem(n or "<sem nome>")
            item.setFlags(item.flags() | Qt.ItemFlag.ItemIsUserCheckable)
            item.setCheckState(Qt.CheckState.Checked)
            item.setData(Qt.ItemDataRole.UserRole, n)
            self.lig_list.addItem(item)
        self._update_count()

    def _set_all(self, checked: bool) -> None:
        state = Qt.CheckState.Checked if checked else Qt.CheckState.Unchecked
        for i in range(self.lig_list.count()):
            it = self.lig_list.item(i)
            if not it.isHidden():
                it.setCheckState(state)
        self._update_count()

    def _apply_filter(self, text: str) -> None:
        text = text.lower().strip()
        for i in range(self.lig_list.count()):
            it = self.lig_list.item(i)
            it.setHidden(bool(text) and text not in it.text().lower())

    def _update_count(self) -> None:
        total = self.lig_list.count()
        sel = sum(
            1 for i in range(total)
            if self.lig_list.item(i).checkState() == Qt.CheckState.Checked
        )
        self.count_label.setText(f"{sel} de {total} ligantes selecionados")

    def collect(self) -> None:
        """Push UI state into self.cfg."""
        self.cfg.protein_file = self.protein_edit.text().strip()
        self.cfg.ligand_file = self.ligand_edit.text().strip()
        self.cfg.workdir = self.workdir_edit.text().strip()
        self.cfg.include_waters = self.cb_waters.isChecked()
        self.cfg.selected_ligands = [
            self.lig_list.item(i).data(Qt.ItemDataRole.UserRole)
            for i in range(self.lig_list.count())
            if self.lig_list.item(i).checkState() == Qt.CheckState.Checked
        ]
        self._update_count()
