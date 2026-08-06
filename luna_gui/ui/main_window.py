"""Main window — assembles all tabs and wires shared state."""
from __future__ import annotations

import logging
from pathlib import Path
from typing import TYPE_CHECKING

from PyQt6.QtWidgets import (
    QApplication,
    QMainWindow,
    QTabWidget,
    QFileDialog,
    QMessageBox,
    QScrollArea,
    QFrame,
    QSizePolicy,
    QWidget,
    QLabel,
    QPushButton,
    QVBoxLayout,
)
from PyQt6.QtCore import QSettings, Qt, QTimer
from PyQt6.QtGui import QAction, QActionGroup, QGuiApplication, QKeySequence

from ..core.project import ProjectConfig, save_to_workdir
from ..i18n import (
    LANGUAGES,
    install_translation_hooks,
    language,
    language_notifier,
    retranslate_ui,
    set_language,
    t,
)
from .tab_setup import SetupTab
from .tab_project import ProjectTab
from .tab_analyses import AnalysesTab
from .tab_run import RunTab
from .tab_history import HistoryTab
from .theme import THEME_MODES, apply_theme

if TYPE_CHECKING:
    from .tab_results_advanced import ResultsTab


LOGGER = logging.getLogger(__name__)


class MainWindow(QMainWindow):
    def __init__(self) -> None:
        super().__init__()
        install_translation_hooks()
        self.setWindowTitle("HIP²LInterActomics")
        self._fit_to_screen()
        self._settings = QSettings()

        preferred_language = str(self._settings.value("interface/language", "en") or "en")
        self.cfg = ProjectConfig(language=preferred_language if preferred_language in LANGUAGES else "en")
        set_language(getattr(self.cfg, "language", "en"))

        self.tabs = QTabWidget()
        self.tabs.setUsesScrollButtons(True)
        self.setCentralWidget(self.tabs)

        self.tab_setup = SetupTab()
        self.tab_project = ProjectTab(self.cfg)
        self.tab_analyses = AnalysesTab(self.cfg)
        self.tab_run = RunTab(self.cfg)
        self.tab_results: ResultsTab | None = None
        self._results_python = ""
        self._results_placeholder = self._build_results_placeholder()
        self.tab_history = HistoryTab()
        self._tab_scroll_pages: dict[QWidget, QScrollArea] = {}

        self._add_scrollable_tab(self.tab_setup, "1. Inicio")
        self._add_scrollable_tab(self.tab_project, "2. Projeto")
        self._add_scrollable_tab(self.tab_analyses, "3. Análises")
        self._add_scrollable_tab(self.tab_run, "4. Executar")
        self._add_scrollable_tab(self._results_placeholder, "5. Resultados")
        self.tabs.currentChanged.connect(self._on_main_tab_changed)
        self._add_scrollable_tab(self.tab_history, "6. Histórico")

        # Wire: when LUNA is detected, hand the paths to RunTab and ResultsTab
        self.tab_setup.luna_ready.connect(self.tab_run.set_luna)
        self.tab_setup.luna_ready.connect(lambda py, _run: self._set_results_python(py))
        QTimer.singleShot(0, self.tab_setup.detect)

        self._build_menu()
        self._language_notifier = language_notifier()
        if self._language_notifier is not None:
            self._language_notifier.language_changed.connect(self._apply_language_change)
        retranslate_ui(self)
        self.statusBar().showMessage(t("Pronto"))

        # Wire: before running, collect state and save project
        self.tab_run.collect_callback = self._collect_and_save

        # Wire: when LUNA finishes, jump to results and load
        self.tab_run.finished_ok.connect(self._on_run_finished)

        # Wire: history reload
        self.tab_history.project_loaded.connect(self._apply_loaded_cfg)

    def _create_scroll_page(self, content: QWidget) -> QScrollArea:
        """Create a responsive scroll page for one main workflow step."""
        scroll = QScrollArea()
        scroll.setObjectName("mainTabScrollArea")
        scroll.setWidgetResizable(True)
        scroll.setAlignment(Qt.AlignmentFlag.AlignTop)
        scroll.setFrameShape(QFrame.Shape.NoFrame)
        scroll.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        content.setSizePolicy(
            QSizePolicy.Policy.Ignored,
            QSizePolicy.Policy.Preferred,
        )
        scroll.setWidget(content)
        self._tab_scroll_pages[content] = scroll
        return scroll

    def _add_scrollable_tab(self, content: QWidget, label: str) -> None:
        """Keep every main tab usable on short screens without shrinking its content."""
        scroll = self._create_scroll_page(content)
        self.tabs.addTab(scroll, label)

    def _build_results_placeholder(self) -> QWidget:
        page = QWidget()
        layout = QVBoxLayout(page)
        layout.setContentsMargins(28, 28, 28, 28)
        layout.addStretch()
        title = QLabel("Resultados sob demanda")
        title.setAlignment(Qt.AlignmentFlag.AlignCenter)
        title.setStyleSheet("font-size: 20px; font-weight: 650;")
        description = QLabel(
            "Os módulos científicos e gráficos serão carregados somente quando necessários, "
            "reduzindo a memória e acelerando a abertura do aplicativo."
        )
        description.setWordWrap(True)
        description.setAlignment(Qt.AlignmentFlag.AlignCenter)
        description.setProperty("muted", True)
        button = QPushButton("Abrir resultados")
        button.setMaximumWidth(220)
        button.clicked.connect(self._ensure_results_tab)
        layout.addWidget(title)
        layout.addWidget(description)
        layout.addWidget(button, alignment=Qt.AlignmentFlag.AlignHCenter)
        layout.addStretch()
        return page

    def _ensure_results_tab(self) -> ResultsTab:
        if self.tab_results is not None:
            return self.tab_results

        from .tab_results_advanced import ResultsTab

        old_scroll = self._tab_scroll_pages.pop(self._results_placeholder)
        index = self.tabs.indexOf(old_scroll)
        label = self.tabs.tabText(index) if index >= 0 else "5. Resultados"
        was_current = self.tabs.currentWidget() is old_scroll

        results = ResultsTab(self.cfg)
        if self._results_python:
            results.set_python(self._results_python)
        new_scroll = self._create_scroll_page(results)

        self.tabs.blockSignals(True)
        if index >= 0:
            self.tabs.removeTab(index)
            self.tabs.insertTab(index, new_scroll, label)
        else:
            index = self.tabs.addTab(new_scroll, label)
        self.tabs.blockSignals(False)
        old_scroll.deleteLater()
        self.tab_results = results
        # The Results tab is created lazily, so apply the active language now.
        if hasattr(results, "retranslate_dynamic"):
            results.retranslate_dynamic()
        if was_current:
            self.tabs.setCurrentIndex(index)
        return results

    def _set_results_python(self, python_path: str) -> None:
        self._results_python = python_path
        if self.tab_results is not None:
            self.tab_results.set_python(python_path)

    def _on_main_tab_changed(self, _index: int) -> None:
        placeholder_scroll = self._tab_scroll_pages.get(self._results_placeholder)
        if placeholder_scroll is not None and self.tabs.currentWidget() is placeholder_scroll:
            QTimer.singleShot(0, self._ensure_results_tab)
        self._sync_cfg_if_idle()

    def _set_current_tab(self, content: QWidget) -> None:
        self.tabs.setCurrentWidget(self._tab_scroll_pages.get(content, content))

    def closeEvent(self, event) -> None:
        answer = QMessageBox.question(
            self,
            "Fechar aplicativo",
            "Deseja realmente fechar o HIP²LInterActomics?",
            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
            QMessageBox.StandardButton.No,
        )
        if answer == QMessageBox.StandardButton.Yes:
            event.accept()
        else:
            event.ignore()

    def _build_menu(self) -> None:
        bar = self.menuBar()
        m_file = bar.addMenu("&Arquivo")

        act_open = QAction("Abrir projeto...", self)
        act_open.setShortcut(QKeySequence.StandardKey.Open)
        act_open.triggered.connect(self._open_project)
        m_file.addAction(act_open)

        act_save = QAction("Salvar projeto", self)
        act_save.setShortcut(QKeySequence.StandardKey.Save)
        act_save.triggered.connect(self._save_project)
        m_file.addAction(act_save)

        m_file.addSeparator()
        act_quit = QAction("Sair", self)
        act_quit.setShortcut(QKeySequence.StandardKey.Quit)
        act_quit.triggered.connect(self.close)
        m_file.addAction(act_quit)

        m_lang = bar.addMenu("&Idioma")
        self._language_actions: dict[str, QAction] = {}
        lang_group = QActionGroup(self)
        lang_group.setExclusive(True)
        for code, label in LANGUAGES.items():
            act_lang = QAction(label, self)
            act_lang.setCheckable(True)
            act_lang.setData(code)
            act_lang.setChecked(code == language())
            act_lang.triggered.connect(lambda _checked=False, c=code: self._set_language(c))
            lang_group.addAction(act_lang)
            m_lang.addAction(act_lang)
            self._language_actions[code] = act_lang

        m_theme = bar.addMenu("&Aparência")
        self._theme_actions: dict[str, QAction] = {}
        theme_group = QActionGroup(self)
        theme_group.setExclusive(True)
        preferred_theme = str(self._settings.value("appearance/theme", "system") or "system")
        theme_labels = {
            "system": "Tema do sistema",
            "light": "Claro",
            "dark": "Escuro",
        }
        for mode in THEME_MODES:
            action = QAction(theme_labels[mode], self)
            action.setCheckable(True)
            action.setChecked(mode == preferred_theme)
            action.triggered.connect(lambda _checked=False, value=mode: self._set_theme(value))
            theme_group.addAction(action)
            m_theme.addAction(action)
            self._theme_actions[mode] = action

        m_help = bar.addMenu("Aj&uda")
        act_about = QAction("Sobre", self)
        act_about.triggered.connect(self._about)
        m_help.addAction(act_about)
        act_docs = QAction("Documentação LUNA", self)
        act_docs.triggered.connect(
            lambda: __import__("webbrowser").open("https://luna-toolkit.readthedocs.io")
        )
        m_help.addAction(act_docs)

    def _set_language(self, code: str) -> None:
        previous = language()
        code = set_language(code)
        if code == previous:
            self._apply_language_change(code)

    def _apply_language_change(self, code: str) -> None:
        self.cfg.language = code
        self._settings.setValue("interface/language", code)
        for lang_code, action in getattr(self, "_language_actions", {}).items():
            action.setChecked(lang_code == code)
        retranslate_ui(self)
        for tab in (
            self.tab_setup,
            self.tab_project,
            self.tab_analyses,
            self.tab_run,
            self.tab_results,
            self.tab_history,
        ):
            if tab is not None and hasattr(tab, "retranslate_dynamic"):
                tab.retranslate_dynamic()
        if self.tab_results is not None and hasattr(self.tab_results, "on_language_changed"):
            self.tab_results.on_language_changed(code)
        self.statusBar().showMessage(t(f"Idioma alterado para {LANGUAGES[code]}"), 4000)

    def _set_theme(self, mode: str) -> None:
        if mode not in THEME_MODES:
            mode = "system"
        self._settings.setValue("appearance/theme", mode)
        app = QApplication.instance()
        if app is not None:
            resolved = apply_theme(app, mode)
            self.statusBar().showMessage(
                t("Tema atualizado") + f": {resolved}",
                3500,
            )
        for theme_mode, action in self._theme_actions.items():
            action.setChecked(theme_mode == mode)

    def _open_project(self) -> None:
        f, _ = QFileDialog.getOpenFileName(
            self, "Abrir projeto", "", "HIP²LInterActomics project (*.json);;Todos (*)"
        )
        if not f:
            return
        try:
            cfg = ProjectConfig.load(f)
        except Exception as e:
            QMessageBox.critical(self, "Erro", str(e)); return
        self._apply_loaded_cfg(cfg)

    def _save_project(self) -> None:
        self.tab_project.collect()
        self.tab_analyses.collect()
        self.cfg.language = language()
        if not self.cfg.workdir:
            QMessageBox.warning(self, "Workdir", t("Defina um workdir antes de salvar."))
            return
        save_to_workdir(self.cfg)
        self.statusBar().showMessage(f"{t('Projeto salvo em')} {self.cfg.workdir}", 5000)

    def _about(self) -> None:
        QMessageBox.about(self, t("Sobre o HIP²LInterActomics"), self._about_html())

    def _about_html(self) -> str:
        lang = language()
        if lang == "en":
            description = (
                "Graphical interface for LUNA focused on protein-ligand and protein-protein "
                "intermolecular interactomics, virtual screening evaluation, molecular dynamics "
                "trajectory analysis, machine-learning data construction and pharmacophoric "
                "feature extraction."
            )
            authors_label = "Authors"
            citations_label = "Citations"
        elif lang == "es":
            description = (
                "Interfaz gráfica para LUNA enfocada en interactómica intermolecular "
                "proteína-ligando y proteína-proteína, evaluación de virtual screening, "
                "análisis de trayectorias de dinámica molecular, construcción de datos para "
                "machine learning y extracción de pharmacophoric features."
            )
            authors_label = "Autores"
            citations_label = "Citas"
        else:
            description = (
                "Interface gráfica para o LUNA focada em interactômica intermolecular "
                "proteína-ligante e proteína-proteína, avaliação de virtual screening, "
                "trajetórias de dinâmica molecular, construção de dados para machine learning "
                "e extração de pharmacophoric features."
            )
            authors_label = "Authores"
            citations_label = "Citações"
        return (
            "<h3>HIP²LInterActomics</h3>"
            f"<p>{description}</p>"
            f"<p><b>{authors_label}:</b> Daniel Andrés Grajales Ruiz e Adriano Marques Gonçalves</p>"
            "<p>LUNA: <a href='https://luna-toolkit.readthedocs.io'>luna-toolkit.readthedocs.io</a></p>"
            f"<p><b>{citations_label}</b></p>"
            "<ol>"
            "<li>Keiser, M. J.; Hert, J. In <i>Chemogenomics: Methods and Applications</i>; "
            "Jacoby, E., Ed.; Methods in Molecular Biology; Humana Press: Totowa, NJ, 2009; "
            "pp 195-205.</li>"
            "<li><i>Prioritizing Virtual Screening with Interpretable Interaction Fingerprints</i><br>"
            "Alexandre V. Fassio, Laura Shub, Luca Ponzoni, Jessica McKinley, Matthew J. O’Meara, "
            "Rafaela S. Ferreira, Michael J. Keiser, and Raquel C. de Melo Minardi<br>"
            "Journal of Chemical Information and Modeling 2022 62 (18), 4300-4318<br>"
            "DOI: 10.1021/acs.jcim.2c00695</li>"
            "<li><i>Extended-Connectivity Fingerprints</i><br>"
            "David Rogers and Mathew Hahn<br>"
            "Journal of Chemical Information and Modeling 2010 50 (5), 742-754<br>"
            "DOI: 10.1021/ci100050t</li>"
            "<li><i>A Simple Representation of Three-Dimensional Molecular Structure</i><br>"
            "Seth D. Axen, Xi-Ping Huang, Elena L. Cáceres, Leo Gendelev, Bryan L. Roth, and Michael J. Keiser<br>"
            "Journal of Medicinal Chemistry 2017 60 (17), 7393-7409<br>"
            "DOI: 10.1021/acs.jmedchem.7b00696</li>"
            "<li>Liao, P-S., Chen, T-S. and Chung, P-C., “A fast algorithm for multilevel thresholding”, "
            "Journal of Information Science and Engineering 17 (5): 713-727, 2001. Available at: "
            "&lt;https://ftp.iis.sinica.edu.tw/JISE/2001/200109_01.pdf&gt;.</li>"
            "<li>Pedregosa, F., Varoquaux, G., Gramfort, A., Michel, V., Thirion, B., Grisel, O., "
            "... &amp; Duchesnay, E. (2011). Scikit-learn: Machine learning in Python. "
            "Journal of Machine Learning Research, 12(Oct), 2825-2830.</li>"
            "</ol>"
        )

    def _collect_and_save(self) -> None:
        self.tab_project.collect()
        self.tab_analyses.collect()
        self.cfg.language = language()
        save_to_workdir(self.cfg)

    def _sync_cfg_if_idle(self) -> None:
        proc = getattr(self.tab_run, "proc", None)
        if proc is not None and proc.state() != proc.ProcessState.NotRunning:
            return
        self.tab_project.collect()
        self.tab_analyses.collect()
        self.cfg.language = language()

    def _on_run_finished(self) -> None:
        results = self._ensure_results_tab()
        self._set_current_tab(results)
        results.wd_edit.setText(self.cfg.workdir)
        results.load_all()
        self.tab_history.refresh()

    def _apply_loaded_cfg(self, cfg: ProjectConfig) -> None:
        # Replace shared cfg fields and reflect in UI
        for k, v in vars(cfg).items():
            setattr(self.cfg, k, v)
        self._set_language(getattr(self.cfg, "language", "en"))
        # Push values back into widgets
        self.tab_project.protein_edit.setText(self.cfg.protein_file)
        self.tab_project.ligand_edit.setText(self.cfg.ligand_file)
        self.tab_project.workdir_edit.setText(self.cfg.workdir)
        self.tab_project.cb_trajectory.setChecked(bool(getattr(self.cfg, "trajectory_analysis", False)))
        if self.cfg.ligand_file:
            try:
                self.tab_project._load_ligands(self.cfg.ligand_file)
                # restore selection
                selected = set(self.cfg.selected_ligands)
                lst = self.tab_project.lig_list
                for i in range(lst.count()):
                    it = lst.item(i)
                    name = it.data(Qt.ItemDataRole.UserRole)
                    it.setCheckState(
                        Qt.CheckState.Checked if name in selected else Qt.CheckState.Unchecked
                    )
                self.tab_project._update_count()
            except Exception as exc:
                LOGGER.warning("Não foi possível restaurar a seleção de ligantes", exc_info=True)
                self.statusBar().showMessage(
                    t("Projeto aberto, mas a seleção de ligantes não pôde ser restaurada")
                    + f": {exc}",
                    7000,
                )
        a = self.tab_analyses
        a.cb_add_h.setChecked(self.cfg.add_h)
        a.sp_ph.setValue(self.cfg.ph)
        a.sp_ph.setEnabled(self.cfg.add_h)
        a.fp_box.setChecked(self.cfg.out_ifp)
        a.set_ifp_type(self.cfg.ifp_type)
        a.sp_levels.setValue(self.cfg.ifp_levels)
        a.sp_radius.setValue(self.cfg.ifp_radius)
        a.sp_length.setValue(self.cfg.ifp_length)
        a.cb_bit.setChecked(self.cfg.ifp_bit)
        a.fp_out_edit.setText(self.cfg.ifp_output)
        if hasattr(a, "ifp_seed_edit"):
            a.ifp_seed_edit.setText(getattr(self.cfg, "ifp_seed_file", "") or "")
        a.fp_labels_box.setChecked(bool(self.cfg.fp_labels_csv) or bool(getattr(self.cfg, "ifp_seed_file", "")))
        a.fp_labels_edit.setText(self.cfg.fp_labels_csv)
        a.fp_label_id_column_edit.setText(getattr(self.cfg, "fp_labels_id_column", "") or "")
        a.fp_label_column_edit.setText(self.cfg.fp_labels_column)
        task_value = str(getattr(self.cfg, "fp_label_task", "regression") or "regression")
        task_idx = a.cb_fp_label_task.findData(task_value)
        a.cb_fp_label_task.setCurrentIndex(task_idx if task_idx >= 0 else 0)
        a.cb_fp_use_otsu.setChecked(bool(getattr(self.cfg, "fp_use_otsu_threshold", False)))
        a.sm_box.setChecked(self.cfg.sim_matrix)
        a.sm_out_edit.setText(self.cfg.sim_matrix_output)
        a.pse_box.setChecked(self.cfg.out_pse)
        a.pse_dir_edit.setText(self.cfg.pse_path)
        a.fbm_box.setChecked(self.cfg.filter_binding_modes)
        a.fbm_edit.setText(self.cfg.binding_modes_cfg)

        # Restore advanced-option widgets
        self.tab_project.cb_fork.setChecked(bool(self.cfg.fork_from))
        self.tab_project.fork_edit.setText(self.cfg.fork_from)
        self.tab_project.cb_waters.setChecked(self.cfg.include_waters)
        self.tab_project.cb_protein_heteroatoms.setChecked(
            bool(getattr(self.cfg, "include_protein_heteroatoms", False))
        )
        a.adv_box.setChecked(
            self.cfg.force_python_api
            or bool(self.cfg.inter_config_overrides)
            or float(getattr(self.cfg, "inter_max_distance_cap", 0.0) or 0.0) > 0.0
        )
        a.inter_cfg_box.setChecked(bool(getattr(self.cfg, "interaction_config_file", "")))
        a.inter_cfg_edit.setText(getattr(self.cfg, "interaction_config_file", "") or "")
        a.cb_ic_proximal.setChecked(self.cfg.ic_add_proximal)
        a.cb_ic_atom_atom.setChecked(self.cfg.ic_add_atom_atom)
        a.cb_ic_dep.setChecked(self.cfg.ic_add_dependent_inter)
        a.cb_ic_h2o_nt.setChecked(self.cfg.ic_add_h2o_pairs_with_no_target)
        a.cb_ic_self.setChecked(self.cfg.ic_ignore_self_inter)
        a.sp_inter_max_cap.setValue(float(getattr(self.cfg, "inter_max_distance_cap", 0.0) or 0.0))
        if self.cfg.pse_interaction_types:
            a.pse_filter_box.setChecked(True)
            selected = set(self.cfg.pse_interaction_types)
            for i in range(a.pse_types_list.count()):
                it = a.pse_types_list.item(i)
                it.setCheckState(
                    Qt.CheckState.Checked if it.text() in selected else Qt.CheckState.Unchecked
                )

        # Auto-populate Results from the loaded workdir (T2.8)
        if self.cfg.workdir and Path(self.cfg.workdir).exists():
            results = self._ensure_results_tab()
            results.wd_edit.setText(self.cfg.workdir)
            try:
                results.load_all()
            except Exception as exc:
                LOGGER.warning("Não foi possível carregar os resultados do projeto", exc_info=True)
                self.statusBar().showMessage(
                    t("Projeto aberto, mas os resultados não puderam ser carregados")
                    + f": {exc}",
                    7000,
                )

        self._set_current_tab(self.tab_project)

    def _fit_to_screen(self) -> None:
        screen = QGuiApplication.primaryScreen()
        if screen is None:
            self.resize(1280, 820)
            self.setMinimumSize(900, 620)
            return
        available = screen.availableGeometry()
        width = max(900, min(1280, int(available.width() * 0.90)))
        height = max(620, min(820, int(available.height() * 0.88)))
        min_width = max(760, min(900, int(available.width() * 0.72)))
        min_height = max(540, min(620, int(available.height() * 0.66)))
        self.resize(width, height)
        self.setMinimumSize(min_width, min_height)
