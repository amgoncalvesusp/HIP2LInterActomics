"""Analyses tab — pick which LUNA outputs to generate and tune parameters."""
from __future__ import annotations

from PyQt6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QFormLayout, QGroupBox, QCheckBox,
    QComboBox, QSpinBox, QDoubleSpinBox, QLineEdit, QPushButton, QFileDialog,
    QListWidget, QListWidgetItem, QLabel, QScrollArea,
)
from PyQt6.QtCore import Qt

from ..core.project import ProjectConfig


# LUNA's built-in interaction types (see luna.interaction.type).
LUNA_INTERACTION_TYPES = [
    "Hydrogen bond", "Weak hydrogen bond", "Halogen bond",
    "Chalcogen bond", "Ionic", "Salt bridge", "Cation-pi",
    "Pi-stacking", "Edge-to-face", "Face-to-face", "Parallel",
    "T-shaped", "Hydrophobic", "Amide-aromatic stacking",
    "Water-bridged hydrogen bond", "Disulfide bond",
    "Metal coordination", "Van der Waals", "Proximal",
]


class AnalysesTab(QWidget):
    def __init__(self, cfg: ProjectConfig) -> None:
        super().__init__()
        self.cfg = cfg

        # Scroll area wraps everything so the advanced section doesn't clip
        outer = QVBoxLayout(self)
        scroll = QScrollArea(); scroll.setWidgetResizable(True)
        inner = QWidget(); layout = QVBoxLayout(inner)
        scroll.setWidget(inner)
        outer.addWidget(scroll)

        intro = QLabel(
            "Aqui você escolhe quais saídas o LUNA vai gerar. Se estiver em dúvida, deixe os "
            "valores padrão e ative apenas o que pretende inspecionar depois: fingerprints, "
            "similaridade, sessões PyMOL ou filtros específicos."
        )
        intro.setWordWrap(True)
        intro.setProperty("muted", True)
        layout.addWidget(intro)

        # ----- Fingerprints -----
        fp_box = QGroupBox("Fingerprints de interação (IFP)")
        fp_box.setCheckable(True)
        fp_box.setChecked(True)
        fp_box.setToolTip("Gera a representação vetorial das interações para cada ligante.")
        self.fp_box = fp_box
        fp_form = QFormLayout(fp_box)

        fp_help = QLabel(
            "Fingerprints transformam as interações em uma tabela comparável entre ligantes. "
            "É a saída mais útil para análise posterior e costuma ser a opção principal."
        )
        fp_help.setWordWrap(True)
        fp_help.setProperty("muted", True)
        fp_form.addRow(fp_help)

        self.cb_type = QComboBox()
        self.cb_type.addItems(["EIFP", "HIFP", "FIFP"])
        self.cb_type.setToolTip(
            "Escolhe o formato do fingerprint. Em geral, mantenha o padrão do seu fluxo; altere apenas se souber qual representação precisa."
        )
        fp_form.addRow("Tipo:", self.cb_type)

        self.sp_levels = QSpinBox(); self.sp_levels.setRange(1, 20); self.sp_levels.setValue(2)
        self.sp_levels.setToolTip("Define quantos níveis de vizinhança estrutural entram no fingerprint.")
        fp_form.addRow("Levels (-L):", self.sp_levels)

        self.sp_radius = QDoubleSpinBox(); self.sp_radius.setRange(0.1, 50.0)
        self.sp_radius.setDecimals(5); self.sp_radius.setValue(5.73171)
        self.sp_radius.setToolTip("Controla o passo radial usado na geração do fingerprint.")
        fp_form.addRow("Radius step (-R):", self.sp_radius)

        self.sp_length = QSpinBox(); self.sp_length.setRange(64, 65536)
        self.sp_length.setSingleStep(256); self.sp_length.setValue(4096)
        self.sp_length.setToolTip("Comprimento final do fingerprint. Valores maiores aumentam detalhamento e tamanho do arquivo.")
        fp_form.addRow("Length (-S):", self.sp_length)

        self.cb_bit = QCheckBox("Usar bit fingerprints (padrão: count)")
        self.cb_bit.setToolTip("Troca contagens por presença/ausência de interações. Útil quando você quer uma assinatura binária.")
        fp_form.addRow(self.cb_bit)

        self.fp_out_edit = QLineEdit()
        self.fp_out_edit.setPlaceholderText("(padrão: <workdir>/results/fingerprints/ifp.csv)")
        self.fp_out_edit.setToolTip("Permite salvar o CSV de fingerprints em um caminho específico.")
        btn_fp = QPushButton("...")
        btn_fp.clicked.connect(lambda: self._pick_save(self.fp_out_edit, "fingerprints.csv"))
        row = QHBoxLayout(); row.addWidget(self.fp_out_edit); row.addWidget(btn_fp)
        fp_form.addRow("Saída IFP:", self._wrap(row))

        layout.addWidget(fp_box)

        # ----- Similarity matrix -----
        self.sm_box = QGroupBox("Matriz de similaridade (Tanimoto)")
        self.sm_box.setCheckable(True)
        self.sm_box.setChecked(False)
        self.sm_box.setToolTip("Gera uma matriz que compara a semelhança entre os fingerprints dos ligantes.")
        sm_form = QFormLayout(self.sm_box)
        sm_help = QLabel(
            "A matriz de similaridade ajuda a ver quais ligantes se comportam de forma parecida. "
            "Ela também alimenta os gráficos e clusters da aba de resultados."
        )
        sm_help.setWordWrap(True)
        sm_help.setProperty("muted", True)
        sm_form.addRow(sm_help)
        self.sm_out_edit = QLineEdit()
        self.sm_out_edit.setPlaceholderText("(padrão: <workdir>/sim_matrix.csv)")
        self.sm_out_edit.setToolTip("Caminho opcional para salvar o CSV da matriz de similaridade.")
        btn_sm = QPushButton("...")
        btn_sm.clicked.connect(lambda: self._pick_save(self.sm_out_edit, "sim_matrix.csv"))
        row = QHBoxLayout(); row.addWidget(self.sm_out_edit); row.addWidget(btn_sm)
        sm_form.addRow("Saída:", self._wrap(row))
        layout.addWidget(self.sm_box)

        # ----- PyMOL sessions -----
        self.pse_box = QGroupBox("Exportar sessões PyMOL (.pse)")
        self.pse_box.setCheckable(True)
        self.pse_box.setChecked(False)
        self.pse_box.setToolTip("Gera sessões prontas para inspeção visual no PyMOL.")
        pse_form = QFormLayout(self.pse_box)
        pse_help = QLabel(
            "As sessões PyMOL permitem abrir a proteína e os ligantes já com as interações destacadas visualmente."
        )
        pse_help.setWordWrap(True)
        pse_help.setProperty("muted", True)
        pse_form.addRow(pse_help)
        self.pse_dir_edit = QLineEdit()
        self.pse_dir_edit.setPlaceholderText("(padrão: <workdir>/results/pse/)")
        self.pse_dir_edit.setToolTip("Pasta onde os arquivos .pse serão gravados.")
        btn_pse = QPushButton("...")
        btn_pse.clicked.connect(self._pick_pse_dir)
        row = QHBoxLayout(); row.addWidget(self.pse_dir_edit); row.addWidget(btn_pse)
        pse_form.addRow("Pasta de saída:", self._wrap(row))
        layout.addWidget(self.pse_box)

        # ----- Filter binding modes -----
        self.fbm_box = QGroupBox("Filtrar por binding modes (.cfg)")
        self.fbm_box.setCheckable(True)
        self.fbm_box.setChecked(False)
        self.fbm_box.setToolTip("Restringe a análise aos modos de ligação definidos em um arquivo .cfg.")
        fbm_form = QFormLayout(self.fbm_box)
        fbm_help = QLabel(
            "Use este filtro quando você já tem regras de binding modes e quer limitar quais poses ou interações entram na análise."
        )
        fbm_help.setWordWrap(True)
        fbm_help.setProperty("muted", True)
        fbm_form.addRow(fbm_help)
        self.fbm_edit = QLineEdit()
        self.fbm_edit.setToolTip("Arquivo .cfg com as definições de binding modes aceitas pelo LUNA.")
        btn_fbm = QPushButton("...")
        btn_fbm.clicked.connect(self._pick_fbm)
        btn_fbm_edit = QPushButton("Editor visual")
        btn_fbm_edit.setToolTip("Abre um editor simples para criar ou ajustar o arquivo de binding modes.")
        btn_fbm_edit.clicked.connect(self._open_fbm_editor)
        row = QHBoxLayout(); row.addWidget(self.fbm_edit); row.addWidget(btn_fbm); row.addWidget(btn_fbm_edit)
        fbm_form.addRow("Arquivo .cfg:", self._wrap(row))
        layout.addWidget(self.fbm_box)

        # ----- PSE filter by interaction type (T2.7) -----
        self.pse_filter_box = QGroupBox("Filtrar PSE por tipo de interação (opcional)")
        self.pse_filter_box.setCheckable(True)
        self.pse_filter_box.setChecked(False)
        self.pse_filter_box.setToolTip("Limita as sessões PyMOL aos tipos de interação selecionados abaixo.")
        pse_f_layout = QVBoxLayout(self.pse_filter_box)
        pse_f_layout.addWidget(QLabel(
            "Se marcado, LUNA gera sessões PyMOL contendo apenas os tipos selecionados."
        ))
        self.pse_types_list = QListWidget()
        self.pse_types_list.setSelectionMode(QListWidget.SelectionMode.NoSelection)
        self.pse_types_list.setMaximumHeight(140)
        self.pse_types_list.setToolTip("Marque quais classes de interação devem aparecer nas sessões PyMOL exportadas.")
        for t in LUNA_INTERACTION_TYPES:
            it = QListWidgetItem(t)
            it.setFlags(it.flags() | Qt.ItemFlag.ItemIsUserCheckable)
            it.setCheckState(Qt.CheckState.Unchecked)
            self.pse_types_list.addItem(it)
        pse_f_layout.addWidget(self.pse_types_list)
        layout.addWidget(self.pse_filter_box)

        # ----- Advanced (T2.5 + T2.6) -----
        self.adv_box = QGroupBox(
            "Opções avançadas — DefaultInteractionConfig + InteractionCalculator"
        )
        self.adv_box.setCheckable(True)
        self.adv_box.setChecked(False)
        self.adv_box.setToolTip(
            "Marcar esta caixa ou alterar qualquer opção aqui faz a GUI usar o "
            "runner Python-API em vez da CLI (mais lento, mais flexível)."
        )
        adv_form = QFormLayout(self.adv_box)
        adv_help = QLabel(
            "Altere estas opções apenas se precisar ajustar critérios finos de detecção de interação. "
            "Para uso geral, o padrão do LUNA costuma ser suficiente."
        )
        adv_help.setWordWrap(True)
        adv_help.setProperty("muted", True)
        adv_form.addRow(adv_help)

        # Distance thresholds (most-used subset of DefaultInteractionConfig)
        self.sp_hb_da = QDoubleSpinBox()
        self.sp_hb_da.setRange(0.1, 20.0); self.sp_hb_da.setDecimals(2); self.sp_hb_da.setSingleStep(0.1)
        self.sp_hb_da.setValue(0.0)
        self.sp_hb_da.setSpecialValueText("(padrão)")
        self.sp_hb_da.setToolTip("max_da_dist_hb_inter — distância máxima doador-aceitador (Å)")
        adv_form.addRow("HB max D–A dist (Å):", self.sp_hb_da)

        self.sp_hb_dha = QDoubleSpinBox()
        self.sp_hb_dha.setRange(0.0, 180.0); self.sp_hb_dha.setDecimals(1); self.sp_hb_dha.setValue(0.0)
        self.sp_hb_dha.setSpecialValueText("(padrão)")
        self.sp_hb_dha.setToolTip("min_dha_ang_hb_inter — ângulo mínimo D–H–A (graus)")
        adv_form.addRow("HB min D–H–A ângulo (°):", self.sp_hb_dha)

        self.sp_hydrop = QDoubleSpinBox()
        self.sp_hydrop.setRange(0.0, 20.0); self.sp_hydrop.setDecimals(2); self.sp_hydrop.setValue(0.0)
        self.sp_hydrop.setSpecialValueText("(padrão)")
        self.sp_hydrop.setToolTip("max_cc_dist_hydrop_inter — distância máxima para contato hidrofóbico (Å)")
        adv_form.addRow("Hidrofóbica max dist (Å):", self.sp_hydrop)

        self.sp_pistack = QDoubleSpinBox()
        self.sp_pistack.setRange(0.0, 20.0); self.sp_pistack.setDecimals(2); self.sp_pistack.setValue(0.0)
        self.sp_pistack.setSpecialValueText("(padrão)")
        self.sp_pistack.setToolTip("max_cc_dist_pi_pi_inter — distância máxima centro-centro para π-stacking (Å)")
        adv_form.addRow("π–π max dist (Å):", self.sp_pistack)

        adv_form.addRow(QLabel("<b>Flags do InteractionCalculator</b>"))
        self.cb_ic_proximal = QCheckBox("add_proximal (inclui contatos só por proximidade)")
        self.cb_ic_proximal.setToolTip("Inclui contatos próximos mesmo quando eles não se encaixam em uma classe química mais específica.")
        self.cb_ic_atom_atom = QCheckBox("add_atom_atom (interações atômicas genéricas)")
        self.cb_ic_atom_atom.setToolTip("Adiciona contatos genéricos átomo-átomo, úteis para inspeção mais detalhada.")
        self.cb_ic_dep = QCheckBox("add_dependent_inter (bridges de água etc.)")
        self.cb_ic_dep.setChecked(True)
        self.cb_ic_dep.setToolTip("Mantém interações dependentes de contexto, como pontes mediadas por água.")
        self.cb_ic_h2o_nt = QCheckBox("add_h2o_pairs_with_no_target (águas sem alvo)")
        self.cb_ic_h2o_nt.setChecked(True)
        self.cb_ic_h2o_nt.setToolTip("Inclui pares com água mesmo quando ela não está ligada diretamente ao alvo principal.")
        self.cb_ic_self = QCheckBox("ignore_self_inter (recomendado)")
        self.cb_ic_self.setChecked(True)
        self.cb_ic_self.setToolTip("Ignora interações internas do mesmo grupo, o que evita ruído na análise principal.")
        for cb in (self.cb_ic_proximal, self.cb_ic_atom_atom, self.cb_ic_dep,
                   self.cb_ic_h2o_nt, self.cb_ic_self):
            adv_form.addRow(cb)

        layout.addWidget(self.adv_box)

        # ----- Fork -----
        fork_box = QGroupBox("Fork de projeto existente (opcional)")
        fork_form = QFormLayout(fork_box)
        fork_help = QLabel(
            "Use esta opção para reaproveitar um projeto LUNA anterior como base, mantendo sua estrutura e resultados associados."
        )
        fork_help.setWordWrap(True)
        fork_help.setProperty("muted", True)
        fork_form.addRow(fork_help)
        self.fork_edit = QLineEdit()
        self.fork_edit.setToolTip("Pasta de um projeto LUNA já existente que servirá de base para o novo projeto.")
        btn_fork = QPushButton("...")
        btn_fork.clicked.connect(self._pick_fork)
        row = QHBoxLayout(); row.addWidget(self.fork_edit); row.addWidget(btn_fork)
        fork_form.addRow("Projeto fonte:", self._wrap(row))
        layout.addWidget(fork_box)

        layout.addStretch()

    def _wrap(self, sub) -> QWidget:
        w = QWidget(); w.setLayout(sub); return w

    def _pick_save(self, edit: QLineEdit, default_name: str) -> None:
        f, _ = QFileDialog.getSaveFileName(self, "Salvar como", default_name, "CSV (*.csv)")
        if f:
            edit.setText(f)

    def _pick_pse_dir(self) -> None:
        d = QFileDialog.getExistingDirectory(self, "Pasta para sessões PyMOL")
        if d:
            self.pse_dir_edit.setText(d)

    def _pick_fbm(self) -> None:
        f, _ = QFileDialog.getOpenFileName(self, "Arquivo de binding modes", "",
                                           "Config (*.cfg);;Todos (*)")
        if f:
            self.fbm_edit.setText(f)

    def _open_fbm_editor(self) -> None:
        from .binding_modes_editor import BindingModesEditor
        dlg = BindingModesEditor(self, initial_path=self.fbm_edit.text().strip())
        if dlg.exec() and dlg.path:
            self.fbm_edit.setText(dlg.path)
            self.fbm_box.setChecked(True)

    def _pick_fork(self) -> None:
        d = QFileDialog.getExistingDirectory(self, "Projeto fonte para fork")
        if d:
            self.fork_edit.setText(d)

    def collect(self) -> None:
        c = self.cfg
        c.out_ifp = self.fp_box.isChecked()
        c.ifp_type = self.cb_type.currentText()
        c.ifp_levels = self.sp_levels.value()
        c.ifp_radius = self.sp_radius.value()
        c.ifp_length = self.sp_length.value()
        c.ifp_bit = self.cb_bit.isChecked()
        c.ifp_output = self.fp_out_edit.text().strip()

        c.sim_matrix = self.sm_box.isChecked()
        c.sim_matrix_output = self.sm_out_edit.text().strip()

        c.out_pse = self.pse_box.isChecked()
        c.pse_path = self.pse_dir_edit.text().strip()

        c.filter_binding_modes = self.fbm_box.isChecked()
        c.binding_modes_cfg = self.fbm_edit.text().strip()

        c.fork_from = self.fork_edit.text().strip()

        # PSE interaction-type filter
        if self.pse_filter_box.isChecked():
            c.pse_interaction_types = [
                self.pse_types_list.item(i).text()
                for i in range(self.pse_types_list.count())
                if self.pse_types_list.item(i).checkState() == Qt.CheckState.Checked
            ]
        else:
            c.pse_interaction_types = []

        # InteractionCalculator flags
        c.ic_add_proximal = self.cb_ic_proximal.isChecked()
        c.ic_add_atom_atom = self.cb_ic_atom_atom.isChecked()
        c.ic_add_dependent_inter = self.cb_ic_dep.isChecked()
        c.ic_add_h2o_pairs_with_no_target = self.cb_ic_h2o_nt.isChecked()
        c.ic_ignore_self_inter = self.cb_ic_self.isChecked()

        # DefaultInteractionConfig overrides (only fields where user entered a non-zero value)
        overrides: dict = {}
        if self.sp_hb_da.value() > 0:
            overrides["max_da_dist_hb_inter"] = self.sp_hb_da.value()
        if self.sp_hb_dha.value() > 0:
            overrides["min_dha_ang_hb_inter"] = self.sp_hb_dha.value()
        if self.sp_hydrop.value() > 0:
            overrides["max_cc_dist_hydrop_inter"] = self.sp_hydrop.value()
        if self.sp_pistack.value() > 0:
            overrides["max_cc_dist_pi_pi_inter"] = self.sp_pistack.value()
        c.inter_config_overrides = overrides

        c.force_python_api = self.adv_box.isChecked()
