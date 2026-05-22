"""Analyses tab — pick which LUNA outputs to generate and tune parameters."""
from __future__ import annotations

from pathlib import Path

from PyQt6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QFormLayout, QGroupBox, QCheckBox,
    QComboBox, QSpinBox, QDoubleSpinBox, QLineEdit, QPushButton, QFileDialog,
    QListWidget, QListWidgetItem, QLabel, QScrollArea, QMessageBox,
)
from PyQt6.QtCore import Qt, QUrl
from PyQt6.QtGui import QDesktopServices

from ..core.project import IFP_ALL, ProjectConfig


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

        hyd_box = QGroupBox("Preparação de hidrogênios")
        hyd_form = QFormLayout(hyd_box)
        hyd_help = QLabel(
            "Define se o LUNA deve protonar as estruturas antes da análise. "
            "Quando ativado, o pH abaixo é usado na adição de hidrogênios."
        )
        hyd_help.setWordWrap(True)
        hyd_help.setProperty("muted", True)
        hyd_form.addRow(hyd_help)

        self.cb_add_h = QCheckBox("Adicionar hidrogênios antes da análise")
        self.cb_add_h.setChecked(True)
        self.cb_add_h.setToolTip(
            "Quando marcado, a GUI deixa o LUNA reprotonar o sistema no pH informado."
        )
        hyd_form.addRow(self.cb_add_h)

        self.sp_ph = QDoubleSpinBox()
        self.sp_ph.setRange(0.0, 14.0)
        self.sp_ph.setDecimals(2)
        self.sp_ph.setSingleStep(0.1)
        self.sp_ph.setValue(7.4)
        self.sp_ph.setToolTip("pH usado na adição de hidrogênios pelo LUNA.")
        hyd_form.addRow("pH:", self.sp_ph)
        self.cb_add_h.toggled.connect(self.sp_ph.setEnabled)
        self.sp_ph.setEnabled(self.cb_add_h.isChecked())

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
        self.cb_type.addItem("EIFP", "EIFP")
        self.cb_type.addItem("HIFP", "HIFP")
        self.cb_type.addItem("FIFP", "FIFP")
        self.cb_type.addItem("Todos (H + E + F)", IFP_ALL)
        self.cb_type.setToolTip(
            "Escolhe o formato do fingerprint.\n"
            "A opção 'Todos (H + E + F)' gera e salva os três tipos: ifp_H.csv, ifp_E.csv e ifp_F.csv."
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
        self.fp_out_edit.setToolTip(
            "Permite salvar o CSV de fingerprints em um caminho específico.\n"
            "Se 'Todos (H + E + F)' estiver selecionado, a pasta deste caminho será usada para salvar ifp_H.csv, ifp_E.csv e ifp_F.csv."
        )
        btn_fp = QPushButton("...")
        btn_fp.clicked.connect(lambda: self._pick_save(self.fp_out_edit, "fingerprints.csv"))
        row = QHBoxLayout(); row.addWidget(self.fp_out_edit); row.addWidget(btn_fp)
        self.ifp_seed_edit = QLineEdit()
        self.ifp_seed_edit.setPlaceholderText("(opcional: arquivo .txt com um inteiro; padrao = 0)")
        self.ifp_seed_edit.setToolTip(
            "Carrega um seed para as etapas estocasticas de importancia dos fingerprints. "
            "Se vazio, a GUI usa seed 0 e salva seed_ifp_H/E/F_importance.txt nos resultados."
        )
        btn_ifp_seed = QPushButton("...")
        btn_ifp_seed.clicked.connect(self._pick_ifp_seed)
        seed_row = QHBoxLayout(); seed_row.addWidget(self.ifp_seed_edit); seed_row.addWidget(btn_ifp_seed)
        fp_form.addRow("Saída IFP:", self._wrap(row))

        self.fp_labels_box = QGroupBox("Rótulos para importância de fingerprints")
        self.fp_labels_box.setCheckable(True)
        self.fp_labels_box.setChecked(False)
        self.fp_labels_box.setToolTip(
            "Configura rótulos supervisionados e seed para o cálculo de importância "
            "das features de fingerprint."
        )
        fp_labels_form = QFormLayout(self.fp_labels_box)
        fp_labels_help = QLabel(
            "Se um CSV for informado, a aba FP análises usa estes rótulos para treinar a importância "
            "das features. O seed controla apenas os modelos estocásticos da importância."
        )
        fp_labels_help.setWordWrap(True)
        fp_labels_help.setProperty("muted", True)
        fp_labels_form.addRow(fp_labels_help)
        fp_labels_form.addRow("Seed da importância:", self._wrap(seed_row))

        self.fp_labels_edit = QLineEdit()
        self.fp_labels_edit.setPlaceholderText("(opcional: CSV/TSV com ligand_id + coluna de rótulo)")
        self.fp_labels_edit.setToolTip(
            "Arquivo CSV ou TSV supervisionado com uma coluna de identificador da molécula e uma coluna de rótulo."
        )
        btn_fp_labels = QPushButton("...")
        btn_fp_labels.clicked.connect(self._pick_fp_labels)
        row = QHBoxLayout(); row.addWidget(self.fp_labels_edit); row.addWidget(btn_fp_labels)
        fp_labels_form.addRow("Arquivo CSV:", self._wrap(row))

        self.fp_label_id_column_edit = QLineEdit()
        self.fp_label_id_column_edit.setPlaceholderText("Ex.: ligand_id, molecule_chembl_id")
        self.fp_label_id_column_edit.setToolTip(
            "Nome da coluna do CSV/TSV que identifica cada molécula. "
            "Se vazio, a GUI tenta detectar automaticamente."
        )
        fp_labels_form.addRow("Coluna do ligand_id:", self.fp_label_id_column_edit)

        self.fp_label_column_edit = QLineEdit()
        self.fp_label_column_edit.setPlaceholderText("Ex.: label, class, activity")
        self.fp_label_column_edit.setToolTip(
            "Nome da coluna do CSV que traz os rótulos usados na análise de importância."
        )
        fp_labels_form.addRow("Coluna de rótulo:", self.fp_label_column_edit)
        self.cb_fp_label_task = QComboBox()
        self.cb_fp_label_task.addItem("Regressão", "regression")
        self.cb_fp_label_task.addItem("Classificação", "classification")
        self.cb_fp_label_task.setCurrentIndex(0)
        self.cb_fp_label_task.setToolTip(
            "Define como a coluna de rótulo será interpretada no cálculo de importância."
        )
        fp_labels_form.addRow("Tarefa:", self.cb_fp_label_task)
        self.cb_fp_use_otsu = QCheckBox("Aplicar Otsu tambem a interacoes/residuos")
        self.cb_fp_use_otsu.setToolTip(
            "A atribuicao de classe das features de FP usa Otsu automaticamente quando "
            "nenhum bit tem z-score > 1. Marque esta opcao para aplicar o mesmo fallback "
            "tambem aos limiares de interacao, residuo e par interacao-residuo."
        )
        fp_labels_form.addRow(self.cb_fp_use_otsu)
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
        self.sm_out_edit.setPlaceholderText("(padrao: <workdir>/sim_matrix_<tipo>.csv)")
        self.sm_out_edit.setToolTip("Caminho opcional para salvar o CSV da matriz de similaridade.")
        btn_sm = QPushButton("...")
        btn_sm.clicked.connect(lambda: self._pick_save(self.sm_out_edit, "sim_matrix_E.csv"))
        row = QHBoxLayout(); row.addWidget(self.sm_out_edit); row.addWidget(btn_sm)
        sm_form.addRow("Saída:", self._wrap(row))
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
        # ----- PSE filter by interaction type (T2.7) -----
        self.pse_filter_box = QGroupBox("Filtrar PSE por tipo de interação (opcional)")
        self.pse_filter_box.setCheckable(True)
        self.pse_filter_box.setChecked(False)
        self.pse_filter_box.setToolTip("Limita as sessões PyMOL aos tipos de interação selecionados abaixo.")
        pse_f_layout = QVBoxLayout(self.pse_filter_box)
        pse_f_layout.addWidget(QLabel(
            "Se marcado, LUNA gera sessões PyMOL contendo apenas os tipos selecionados."
        ))
        self.cb_pse_select_all = QCheckBox("Selecionar/desselecionar todas")
        self.cb_pse_select_all.setToolTip("Marca ou desmarca todos os tipos de interacao do filtro PSE.")
        self.cb_pse_select_all.toggled.connect(self._set_all_pse_interaction_types)
        pse_f_layout.addWidget(self.cb_pse_select_all)
        self.pse_types_list = QListWidget()
        self.pse_types_list.setSelectionMode(QListWidget.SelectionMode.NoSelection)
        self.pse_types_list.setMaximumHeight(140)
        self.pse_types_list.setToolTip("Marque quais classes de interação devem aparecer nas sessões PyMOL exportadas.")
        for t in LUNA_INTERACTION_TYPES:
            it = QListWidgetItem(t)
            it.setFlags(it.flags() | Qt.ItemFlag.ItemIsUserCheckable)
            it.setCheckState(Qt.CheckState.Unchecked)
            self.pse_types_list.addItem(it)
        self.pse_types_list.itemChanged.connect(self._sync_pse_select_all_checkbox)
        pse_f_layout.addWidget(self.pse_types_list)
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
            "Use este quadro apenas para limitar globalmente as distâncias máximas do cálculo "
            "ou ajustar flags do InteractionCalculator. Se quiser mudar parâmetros específicos "
            "de cada interação, use o quadro do arquivo completo .cfg logo abaixo."
        )
        adv_help.setWordWrap(True)
        adv_help.setProperty("muted", True)
        adv_form.addRow(adv_help)

        self.dist_cap_box = QGroupBox("Distância máxima global do cálculo de interações")
        dist_cap_form = QFormLayout(self.dist_cap_box)
        dist_cap_help = QLabel(
            "Se definido, qualquer distância máxima do LUNA acima deste valor será limitada "
            "para este projeto. Distâncias padrão menores do que o valor informado serão mantidas."
        )
        dist_cap_help.setWordWrap(True)
        dist_cap_help.setProperty("muted", True)
        dist_cap_form.addRow(dist_cap_help)

        self.sp_inter_max_cap = QDoubleSpinBox()
        self.sp_inter_max_cap.setRange(0.0, 50.0)
        self.sp_inter_max_cap.setDecimals(2)
        self.sp_inter_max_cap.setSingleStep(0.1)
        self.sp_inter_max_cap.setValue(0.0)
        self.sp_inter_max_cap.setSpecialValueText("(sem limite global)")
        self.sp_inter_max_cap.setToolTip(
            "Cap global por projeto para distâncias máximas do cálculo de interações."
        )
        dist_cap_form.addRow("Distância máxima global (A):", self.sp_inter_max_cap)
        adv_form.addRow(self.dist_cap_box)

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

        self.inter_cfg_box = QGroupBox("Arquivo completo de interações (.cfg)")
        self.inter_cfg_box.setCheckable(True)
        self.inter_cfg_box.setChecked(False)
        self.inter_cfg_box.setToolTip(
            "Permite carregar um arquivo .cfg com todas as interações e parâmetros "
            "de cálculo. Quando marcado, a GUI usa esse arquivo como base no runner Python."
        )
        inter_cfg_form = QFormLayout(self.inter_cfg_box)
        inter_cfg_help = QLabel(
            "Use este quadro quando quiser controlar, em um único arquivo .cfg, "
            "quais interações o LUNA deve calcular e todos os parâmetros de cada uma. "
            "Os ajustes do quadro avançado acima continuam sendo aplicados por cima deste arquivo."
        )
        inter_cfg_help.setWordWrap(True)
        inter_cfg_help.setProperty("muted", True)
        inter_cfg_form.addRow(inter_cfg_help)

        self.inter_cfg_edit = QLineEdit()
        self.inter_cfg_edit.setPlaceholderText("(opcional: config_interacoes.cfg)")
        self.inter_cfg_edit.setToolTip(
            "Arquivo .cfg completo no formato aceito por luna.interaction.config.InteractionConfig."
        )
        btn_inter_cfg = QPushButton("...")
        btn_inter_cfg.clicked.connect(self._pick_inter_cfg)
        row = QHBoxLayout(); row.addWidget(self.inter_cfg_edit); row.addWidget(btn_inter_cfg)
        inter_cfg_form.addRow("Arquivo .cfg:", self._wrap(row))

        self.inter_cfg_link = QLabel(
            "<a href='open_default_interaction_config'>Abrir exemplo padrão utilizado pelo LUNA</a>"
        )
        self.inter_cfg_link.setOpenExternalLinks(False)
        self.inter_cfg_link.linkActivated.connect(self._open_default_interaction_config)
        self.inter_cfg_link.setToolTip(
            "Abre o arquivo de exemplo padrão do LUNA para servir de base na edição."
        )
        inter_cfg_form.addRow(self.inter_cfg_link)

        layout.addWidget(hyd_box)
        layout.addWidget(self.adv_box)
        layout.addWidget(self.inter_cfg_box)
        layout.addWidget(self.pse_box)
        layout.addWidget(self.pse_filter_box)
        layout.addWidget(self.fbm_box)
        layout.addWidget(fp_box)
        layout.addWidget(self.fp_labels_box)
        layout.addWidget(self.sm_box)

        layout.addStretch()

    def _wrap(self, sub) -> QWidget:
        w = QWidget(); w.setLayout(sub); return w

    def _set_all_pse_interaction_types(self, checked: bool) -> None:
        if not hasattr(self, "pse_types_list"):
            return
        state = Qt.CheckState.Checked if checked else Qt.CheckState.Unchecked
        self.pse_types_list.blockSignals(True)
        try:
            for i in range(self.pse_types_list.count()):
                self.pse_types_list.item(i).setCheckState(state)
        finally:
            self.pse_types_list.blockSignals(False)

    def _sync_pse_select_all_checkbox(self) -> None:
        if not hasattr(self, "cb_pse_select_all"):
            return
        total = self.pse_types_list.count()
        checked = sum(
            1
            for i in range(total)
            if self.pse_types_list.item(i).checkState() == Qt.CheckState.Checked
        )
        self.cb_pse_select_all.blockSignals(True)
        self.cb_pse_select_all.setChecked(total > 0 and checked == total)
        self.cb_pse_select_all.blockSignals(False)

    def _pick_save(self, edit: QLineEdit, default_name: str) -> None:
        f, _ = QFileDialog.getSaveFileName(self, "Salvar como", default_name, "CSV (*.csv)")
        if f:
            edit.setText(f)

    def _pick_pse_dir(self) -> None:
        d = QFileDialog.getExistingDirectory(self, "Pasta para sessões PyMOL")
        if d:
            self.pse_dir_edit.setText(d)

    def _pick_fp_labels(self) -> None:
        f, _ = QFileDialog.getOpenFileName(
            self,
            "CSV/TSV de rótulos para FP análises",
            "",
            "Tabela (*.csv *.tsv);;CSV (*.csv);;TSV (*.tsv);;Todos (*)",
        )
        if f:
            self.fp_labels_edit.setText(f)

    def _pick_ifp_seed(self) -> None:
        f, _ = QFileDialog.getOpenFileName(
            self,
            "Arquivo seed IFP",
            "",
            "Texto (*.txt);;Todos (*)",
        )
        if f:
            self.ifp_seed_edit.setText(f)

    def _pick_fbm(self) -> None:
        f, _ = QFileDialog.getOpenFileName(self, "Arquivo de binding modes", "",
                                           "Config (*.cfg);;Todos (*)")
        if f:
            self.fbm_edit.setText(f)

    def _pick_inter_cfg(self) -> None:
        f, _ = QFileDialog.getOpenFileName(
            self,
            "Arquivo completo de interações",
            "",
            "Config (*.cfg);;Todos (*)",
        )
        if f:
            self.inter_cfg_edit.setText(f)
            self.inter_cfg_box.setChecked(True)

    def _default_interaction_config_example(self) -> Path:
        candidates = [
            Path.home() / ".conda" / "envs" / "luna-env" / "Lib" / "site-packages" / "luna" / "interaction" / "config.cfg",
            Path(__file__).resolve().parents[1] / "examples" / "luna_default_interaction_config.cfg",
        ]
        for candidate in candidates:
            if candidate.exists():
                return candidate
        return candidates[-1]

    def _open_default_interaction_config(self, _link: str) -> None:
        path = self._default_interaction_config_example()
        if not path.exists():
            QMessageBox.warning(
                self,
                "Exemplo não encontrado",
                f"Não foi possível localizar o exemplo padrão do LUNA.\n\nEsperado em:\n{path}",
            )
            return
        QDesktopServices.openUrl(QUrl.fromLocalFile(str(path)))

    def _open_fbm_editor(self) -> None:
        from .binding_modes_editor import BindingModesEditor
        dlg = BindingModesEditor(self, initial_path=self.fbm_edit.text().strip())
        if dlg.exec() and dlg.path:
            self.fbm_edit.setText(dlg.path)
            self.fbm_box.setChecked(True)

    def current_ifp_type(self) -> str:
        return self.cb_type.currentData() or self.cb_type.currentText()

    def set_ifp_type(self, value: str) -> None:
        idx = self.cb_type.findData(value)
        if idx < 0:
            idx = self.cb_type.findText(value)
        if idx >= 0:
            self.cb_type.setCurrentIndex(idx)

    def collect(self) -> None:
        c = self.cfg
        c.add_h = self.cb_add_h.isChecked()
        c.ph = self.sp_ph.value()

        c.out_ifp = self.fp_box.isChecked()
        c.ifp_type = self.current_ifp_type()
        c.ifp_levels = self.sp_levels.value()
        c.ifp_radius = self.sp_radius.value()
        c.ifp_length = self.sp_length.value()
        c.ifp_bit = self.cb_bit.isChecked()
        c.ifp_output = self.fp_out_edit.text().strip()
        c.ifp_seed_file = self.ifp_seed_edit.text().strip()
        c.fp_labels_csv = self.fp_labels_edit.text().strip() if self.fp_labels_box.isChecked() else ""
        c.fp_labels_id_column = self.fp_label_id_column_edit.text().strip() if self.fp_labels_box.isChecked() else ""
        c.fp_labels_column = self.fp_label_column_edit.text().strip() if self.fp_labels_box.isChecked() else ""
        c.fp_label_task = str(self.cb_fp_label_task.currentData() or "regression") if self.fp_labels_box.isChecked() else "regression"
        c.fp_use_otsu_threshold = bool(self.fp_labels_box.isChecked() and self.cb_fp_use_otsu.isChecked())

        c.sim_matrix = self.sm_box.isChecked()
        c.sim_matrix_output = self.sm_out_edit.text().strip()

        c.out_pse = self.pse_box.isChecked()
        c.pse_path = self.pse_dir_edit.text().strip()

        c.filter_binding_modes = self.fbm_box.isChecked()
        c.binding_modes_cfg = self.fbm_edit.text().strip()

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
        c.interaction_config_file = (
            self.inter_cfg_edit.text().strip() if self.inter_cfg_box.isChecked() else ""
        )
        c.inter_max_distance_cap = self.sp_inter_max_cap.value() if self.adv_box.isChecked() else 0.0

        c.inter_config_overrides = {}

        c.force_python_api = self.adv_box.isChecked()
