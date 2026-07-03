"""Project tab — pick protein, ligand library, select ligands, choose workdir."""
from __future__ import annotations

from pathlib import Path

from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QFormLayout, QLabel, QLineEdit,
    QPushButton, QFileDialog, QListWidget, QListWidgetItem, QMessageBox,
    QGroupBox, QCheckBox, QScrollArea, QFrame, QSizePolicy,
)

from ..core.ligand_io import (
    parse_ligand_file, consolidate_ligand_folder,
)
from ..core.mol2_prep import count_water_molecules_in_inputs
from ..core.project import ProjectConfig
from .info import InfoButton


class ProjectTab(QWidget):
    def __init__(self, cfg: ProjectConfig) -> None:
        super().__init__()
        self.cfg = cfg

        outer = QVBoxLayout(self)
        outer.setContentsMargins(0, 0, 0, 0)
        outer.setSpacing(0)

        self.scroll = QScrollArea()
        self.scroll.setWidgetResizable(True)
        self.scroll.setFrameShape(QFrame.Shape.NoFrame)
        outer.addWidget(self.scroll)

        content = QWidget()
        layout = QVBoxLayout(content)
        layout.setContentsMargins(12, 12, 12, 12)
        layout.setSpacing(10)
        self.scroll.setWidget(content)

        intro = QLabel(
            "Nesta aba você define o modo do projeto, escolhe o diretório de trabalho, pode "
            "pré-processar arquivos de complexos de docking ou frames de dinâmica e então informa proteína e ligantes para a análise."
        )
        intro.setWordWrap(True)
        intro.setProperty("muted", True)
        layout.addWidget(intro)

        # --- Inputs ---
        form_box = QGroupBox("Entradas")
        form = QFormLayout(form_box)
        form.setContentsMargins(18, 18, 18, 16)
        form.setHorizontalSpacing(12)
        form.setVerticalSpacing(10)
        form.setFieldGrowthPolicy(QFormLayout.FieldGrowthPolicy.AllNonFixedFieldsGrow)
        form.setRowWrapPolicy(QFormLayout.RowWrapPolicy.WrapLongRows)
        form.setLabelAlignment(Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter)

        self.cb_fork = QCheckBox("Fork de projeto existente (desmarcado = projeto novo)")
        self.cb_fork.setToolTip(
            "Marque para usar um projeto LUNA já existente como base.\n"
            "Os novos ligantes e proteínas serão acrescentados ao fork no diretório de trabalho escolhido."
        )
        self.cb_fork.toggled.connect(self._on_fork_toggled)
        form.addRow("Modo do projeto:", self.cb_fork)

        self.fork_edit = QLineEdit()
        self.fork_edit.setToolTip(
            "Pasta do projeto existente que servirá como fonte para o fork."
        )
        btn_fork = QPushButton("Procurar...")
        btn_fork.setMinimumWidth(118)
        btn_fork.setToolTip("Seleciona a pasta do projeto fonte.")
        btn_fork.clicked.connect(self._pick_fork)
        fork_row = QHBoxLayout(); fork_row.addWidget(self.fork_edit, 1); fork_row.addWidget(btn_fork)
        fork_row.addWidget(InfoButton("Use quando quiser reabrir um projeto LUNA existente e acrescentar novas entradas sem apagar resultados antigos."))
        self.fork_row_widget = self._wrap(fork_row)
        form.addRow("Projeto fonte:", self.fork_row_widget)

        self.workdir_edit = QLineEdit()
        self.workdir_edit.setToolTip(
            "Pasta do projeto. A GUI grava aqui o entries.txt, os resultados do LUNA e os arquivos exportados."
        )
        btn_w = QPushButton("Procurar...")
        btn_w.setMinimumWidth(118)
        btn_w.setToolTip("Escolhe a pasta onde o projeto será salvo e reaberto depois.")
        btn_w.clicked.connect(self._pick_workdir)
        wrow = QHBoxLayout(); wrow.addWidget(self.workdir_edit, 1); wrow.addWidget(btn_w)
        wrow.addWidget(InfoButton("Diretório raiz do projeto. A GUI salva configurações, entradas, logs, resultados, relatórios e sessões PyMOL nesta pasta."))
        form.addRow("Diretório de trabalho:", self._wrap(wrow))

        prep_row = QHBoxLayout()
        btn_prep = QPushButton("Preparar arquivos de complexos...")
        btn_prep.setMinimumWidth(260)
        btn_prep.setToolTip("Divide complexos MOL2/PDB (proteína+ligante+águas) em entradas separadas.")
        btn_prep.clicked.connect(self._open_prep_wizard)
        prep_row.addWidget(btn_prep)
        prep_row.addWidget(InfoButton("Divide arquivos de complexos em proteína, ligante e águas. Aceita MOL2/PDB e remove pares LP de águas para evitar confusão no LUNA."))
        prep_row.addStretch()
        form.addRow("Pré-processamento:", self._wrap(prep_row))

        self.protein_edit = QLineEdit()
        self.protein_edit.setToolTip(
            "Arquivo PDB da proteína para todos os ligantes, ou pasta com um PDB por ligante/complexo.\n"
            "Quando uma pasta é usada, cada ligante será pareado com o PDB de mesmo nome-base."
        )
        self.protein_edit.textChanged.connect(lambda _text: self._update_water_count())
        self.btn_protein_file = QPushButton("Arquivo...")
        self.btn_protein_file.setMinimumWidth(106)
        self.btn_protein_file.setToolTip("Seleciona um arquivo PDB de proteína para todos os ligantes.")
        self.btn_protein_file.clicked.connect(self._pick_protein_file)
        self.btn_protein_folder = QPushButton("Pasta...")
        self.btn_protein_folder.setMinimumWidth(106)
        self.btn_protein_folder.setToolTip("Seleciona uma pasta com um arquivo .pdb para cada ligante/complexo.")
        self.btn_protein_folder.clicked.connect(self._pick_protein_folder)
        prow = QHBoxLayout()
        prow.addWidget(self.protein_edit, 1)
        prow.addWidget(self.btn_protein_file)
        prow.addWidget(self.btn_protein_folder)
        prow.addWidget(InfoButton("Arquivo PDB único: uma proteína para todos os ligantes. Pasta: usa a proteína correspondente a cada ligante/frame pelo nome-base."))
        form.addRow("Proteína (PDB):", self._wrap(prow))

        self.ligand_edit = QLineEdit()
        self.ligand_edit.setToolTip(
            "Arquivo ou pasta com ligantes. Pode ser MOL2/SDF/PDB único, pasta de PDBs ou arquivo consolidado com vários ligantes."
        )
        btn_l = QPushButton("Arquivo...")
        btn_l.setMinimumWidth(106)
        btn_l.setToolTip("Seleciona um arquivo de ligantes já pronto para leitura.")
        btn_l.clicked.connect(self._pick_ligand)
        btn_lf = QPushButton("Pasta MOL2/SDF")
        btn_lf.setMinimumWidth(154)
        btn_lf.setToolTip(
            "Consolida todos os .mol2 ou todos os .sdf/.mol da pasta em um arquivo "
            "único. Para MOL2, remove LP e renumera átomos."
        )
        btn_lf.clicked.connect(self._pick_ligand_folder)
        lrow = QHBoxLayout(); lrow.addWidget(self.ligand_edit, 1); lrow.addWidget(btn_l); lrow.addWidget(btn_lf)
        lrow.addWidget(InfoButton("Entrada de ligantes ou frames. Pode ser arquivo MOL2/SDF/PDB ou pasta para consolidar/selecionar multiplas moleculas."))
        form.addRow("Ligantes (MOL2/SDF):", self._wrap(lrow))

        self.cb_waters = QCheckBox("Incluir águas (HOH) — análise hidratada")
        self.cb_waters.setToolTip(
            "Quando marcado, HOH é mantido no PDB e LUNA é chamado com "
            "ignore_any_h2o=False (inclui interações mediadas por água).\n"
            "Requer opções avançadas — usa o runner Python API."
        )
        self.cb_waters.setToolTip(
            "Quando marcado, HOH é mantido no PDB e LUNA é chamado com "
            "ignore_any_h2o=False (inclui interações mediadas por água).\n"
            "Se a pasta do PDB contiver um arquivo por complexo, cada ligante "
            "será pareado com a proteína de mesmo nome-base.\n"
            "Requer opções avançadas - usa o runner Python API."
        )
        self.cb_waters.toggled.connect(self._on_waters_toggled)
        water_row = QHBoxLayout()
        water_row.addWidget(self.cb_waters)
        water_row.addWidget(InfoButton("Mantém águas do complexo no PDB da proteína e executa LUNA permitindo interações mediadas por água."))
        water_row.addStretch()
        form.addRow(self._wrap(water_row))

        self.water_count_label = QLabel("Águas detectadas nos inputs: 0")
        self.water_count_label.setProperty("muted", True)
        self.water_count_label.setWordWrap(True)
        self.water_count_label.setToolTip(
            "Conta resíduos de água reconhecidos como HOH/WAT/H2O/OH2/DOD nos PDB/MOL2 informados."
        )
        form.addRow("", self.water_count_label)

        layout.addWidget(form_box)
        self._on_fork_toggled(False)
        self._on_waters_toggled(bool(getattr(self.cfg, "include_waters", False)))

        self.cb_trajectory = QCheckBox("Análise de trajetória de dinâmica molecular/poses de docking (entradas = frames/poses)")
        self.cb_trajectory.setToolTip(
            "Marque quando cada entrada do projeto representar um frame de dinâmica molecular "
            "ou uma pose de docking. A aba Resultados > Estatísticas mostrará gráficos por "
            "frame/pose e percentuais de entradas por resíduo/interação."
        )
        self.cb_trajectory.setChecked(bool(getattr(self.cfg, "trajectory_analysis", False)))
        self.cb_trajectory.toggled.connect(lambda v: setattr(self.cfg, "trajectory_analysis", bool(v)))
        trajectory_row = QHBoxLayout()
        trajectory_row.addWidget(self.cb_trajectory)
        trajectory_row.addWidget(InfoButton("Ative quando cada entrada representa um frame/pose. Os gráficos passam a usar percentuais de frames/entradas."))
        trajectory_row.addStretch()
        layout.addWidget(self._wrap(trajectory_row))

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
        btn_clear = QPushButton("Limpar detecção")
        btn_clear.setToolTip("Remove todos os ligantes detectados da lista atual sem apagar os arquivos de entrada.")
        btn_clear.clicked.connect(self._clear_detected_ligands)
        filter_row.addWidget(self.filter_edit, 1)
        filter_row.addWidget(btn_all)
        filter_row.addWidget(btn_none)
        filter_row.addWidget(btn_clear)
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
        sub_layout.setContentsMargins(0, 0, 0, 0)
        sub_layout.setSpacing(8)
        w = QWidget()
        w.setLayout(sub_layout)
        w.setMinimumHeight(38)
        w.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        return w

    def _pick_protein_file(self) -> None:
        f, _ = QFileDialog.getOpenFileName(self, "Selecionar proteína", "",
                                           "PDB (*.pdb *.ent);;Todos (*)")
        if f:
            self.protein_edit.setText(f)
            self.cfg.protein_file = f
            self._update_water_count()

    def _pick_protein_folder(self) -> None:
        d = QFileDialog.getExistingDirectory(
            self,
            "Selecionar pasta com proteínas PDB",
            self.protein_edit.text().strip(),
        )
        if d:
            self.protein_edit.setText(d)
            self.cfg.protein_file = d
            self._update_water_count()

    def _pick_fork(self) -> None:
        d = QFileDialog.getExistingDirectory(self, "Projeto fonte para fork")
        if d:
            self.fork_edit.setText(d)
            self.cb_fork.setChecked(True)

    def _on_fork_toggled(self, enabled: bool) -> None:
        self.fork_row_widget.setVisible(enabled)
        if not enabled:
            self.fork_edit.clear()
            self.cfg.fork_from = ""

    def _pick_ligand(self) -> None:
        f, _ = QFileDialog.getOpenFileName(self, "Selecionar ligantes", "",
                                           "Ligantes (*.mol2 *.sdf *.sd *.mol *.pdb *.ent);;Todos (*)")
        if not f:
            return
        self.ligand_edit.setText(f)
        self.cfg.ligand_file = f
        self._load_ligands(f)
        self._update_water_count()

    def _pick_ligand_folder(self) -> None:
        d = QFileDialog.getExistingDirectory(self, "Pasta com arquivos MOL2/SDF")
        if not d:
            return
        try:
            n, _names, out = consolidate_ligand_folder(d, use_file_stem_as_name=True)
        except Exception as e:
            QMessageBox.critical(self, "Erro ao consolidar", str(e))
            return
        if n == 0:
            QMessageBox.warning(
                self,
                "Sem ligantes",
                "Nenhum arquivo .mol2, .sdf ou .mol encontrado na pasta.",
            )
            return
        self.ligand_edit.setText(str(out))
        self.cfg.ligand_file = str(out)
        self._load_ligands(str(out))
        self._update_water_count()
        fmt = out.suffix.upper().lstrip(".")
        QMessageBox.information(
            self, "Consolidado",
            f"{n} ligantes {fmt} consolidados em\n{out}"
        )

    def _open_prep_wizard(self) -> None:
        from .dialog_prep import DockingPrepDialog
        dlg = DockingPrepDialog(self)
        dlg.exec()
        if dlg.cb_open_after.isChecked() and dlg.result_ligand_dir:
            ligand_dir = Path(dlg.result_ligand_dir)
            pdb_ligands = any(
                item.is_file() and item.suffix.lower() in {".pdb", ".ent"}
                for item in ligand_dir.iterdir()
            )
            if pdb_ligands:
                self.ligand_edit.setText(str(ligand_dir))
                self.cfg.ligand_file = str(ligand_dir)
                self._load_ligands(str(ligand_dir))
                n = self.lig_list.count()
            else:
                try:
                    n, _names, out = consolidate_ligand_folder(
                        ligand_dir,
                        use_file_stem_as_name=False,
                    )
                except Exception as e:
                    QMessageBox.critical(self, "Erro ao consolidar", str(e))
                    return
                if n > 0:
                    self.ligand_edit.setText(str(out))
                    self.cfg.ligand_file = str(out)
                    self._load_ligands(str(out))
            # Point the project at the protein folder generated for each complex.
            if dlg.result_protein_dir:
                self.protein_edit.setText(str(dlg.result_protein_dir))
                self.cfg.protein_file = str(dlg.result_protein_dir)
                self._update_water_count()
                water_msg = ""
                if getattr(dlg, "water_molecules_detected", 0):
                    water_msg = f"\nÁguas detectadas na preparação: {dlg.water_molecules_detected}."
                QMessageBox.information(
                    self, "Próximo passo",
                    "Os ligantes extraídos foram carregados automaticamente.\n"
                    "O campo de proteína usa a pasta de PDBs gerados para parear cada ligante/complexo:\n"
                    + dlg.result_protein_dir
                    + water_msg,
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
        self._update_water_count()

    def _set_all(self, checked: bool) -> None:
        state = Qt.CheckState.Checked if checked else Qt.CheckState.Unchecked
        for i in range(self.lig_list.count()):
            it = self.lig_list.item(i)
            if not it.isHidden():
                it.setCheckState(state)
        self._update_count()

    def _clear_detected_ligands(self) -> None:
        self.lig_list.clear()
        self.filter_edit.clear()
        self.cfg.selected_ligands = []
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

    def _on_waters_toggled(self, enabled: bool) -> None:
        self.cfg.include_waters = bool(enabled)
        self._update_protein_picker_mode()
        self._update_water_count()

    def _update_protein_picker_mode(self) -> None:
        self.protein_edit.setToolTip(
            "Arquivo PDB da proteína para todos os ligantes, ou pasta com um PDB por ligante/complexo.\n"
            "Quando uma pasta é usada, cada ligante será pareado com o PDB de mesmo nome-base."
        )

    def _update_water_count(self) -> None:
        if not hasattr(self, "water_count_label"):
            return
        hydrated = self.cb_waters.isChecked()
        self.water_count_label.setVisible(hydrated)
        if not hydrated:
            return
        n_waters = count_water_molecules_in_inputs(
            self.protein_edit.text().strip(),
            self.ligand_edit.text().strip(),
        )
        self.water_count_label.setText(f"Águas detectadas nos inputs: {n_waters}")

    def collect(self) -> None:
        """Push UI state into self.cfg."""
        self.cfg.fork_from = self.fork_edit.text().strip() if self.cb_fork.isChecked() else ""
        self.cfg.protein_file = self.protein_edit.text().strip()
        self.cfg.ligand_file = self.ligand_edit.text().strip()
        self.cfg.workdir = self.workdir_edit.text().strip()
        self.cfg.include_waters = self.cb_waters.isChecked()
        self.cfg.trajectory_analysis = self.cb_trajectory.isChecked()
        self.cfg.selected_ligands = [
            self.lig_list.item(i).data(Qt.ItemDataRole.UserRole)
            for i in range(self.lig_list.count())
            if self.lig_list.item(i).checkState() == Qt.CheckState.Checked
        ]
        self._update_count()
        self._update_water_count()
