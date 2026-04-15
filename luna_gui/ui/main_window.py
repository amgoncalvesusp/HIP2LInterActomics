"""Main window — assembles all tabs and wires shared state."""
from __future__ import annotations

from pathlib import Path

from PyQt6.QtWidgets import QMainWindow, QTabWidget, QFileDialog, QMessageBox
from PyQt6.QtGui import QAction

from ..core.project import ProjectConfig, save_to_workdir
from .tab_setup import SetupTab
from .tab_project import ProjectTab
from .tab_analyses import AnalysesTab
from .tab_run import RunTab
from .tab_results_enhanced import ResultsTab
from .tab_history import HistoryTab


class MainWindow(QMainWindow):
    def __init__(self) -> None:
        super().__init__()
        self.setWindowTitle("LUNA GUI")
        self.resize(1100, 760)

        self.cfg = ProjectConfig()

        self.tabs = QTabWidget()
        self.setCentralWidget(self.tabs)

        self.tab_setup = SetupTab()
        self.tab_project = ProjectTab(self.cfg)
        self.tab_analyses = AnalysesTab(self.cfg)
        self.tab_run = RunTab(self.cfg)
        self.tab_results = ResultsTab(self.cfg)
        self.tab_history = HistoryTab()

        self.tabs.addTab(self.tab_setup, "1. Setup")
        self.tabs.addTab(self.tab_project, "2. Projeto")
        self.tabs.addTab(self.tab_analyses, "3. Análises")
        self.tabs.addTab(self.tab_run, "4. Executar")
        self.tabs.addTab(self.tab_results, "5. Resultados")
        self.tabs.addTab(self.tab_history, "6. Histórico")

        # Wire: when LUNA is detected, hand the paths to RunTab and ResultsTab
        self.tab_setup.luna_ready.connect(self.tab_run.set_luna)
        self.tab_setup.luna_ready.connect(
            lambda py, _run: self.tab_results.set_python(py)
        )

        self._build_menu()
        self.statusBar().showMessage("Pronto")

        # Wire: before running, collect state and save project
        self.tab_run.collect_callback = self._collect_and_save

        # Wire: when LUNA finishes, jump to results and load
        self.tab_run.finished_ok.connect(self._on_run_finished)

        # Wire: history reload
        self.tab_history.project_loaded.connect(self._apply_loaded_cfg)

    def _build_menu(self) -> None:
        bar = self.menuBar()
        m_file = bar.addMenu("&Arquivo")

        act_open = QAction("Abrir projeto...", self)
        act_open.triggered.connect(self._open_project)
        m_file.addAction(act_open)

        act_save = QAction("Salvar projeto", self)
        act_save.triggered.connect(self._save_project)
        m_file.addAction(act_save)

        m_file.addSeparator()
        act_quit = QAction("Sair", self)
        act_quit.triggered.connect(self.close)
        m_file.addAction(act_quit)

        m_help = bar.addMenu("Aj&uda")
        act_about = QAction("Sobre", self)
        act_about.triggered.connect(self._about)
        m_help.addAction(act_about)
        act_docs = QAction("Documentação LUNA", self)
        act_docs.triggered.connect(
            lambda: __import__("webbrowser").open("https://luna-toolkit.readthedocs.io")
        )
        m_help.addAction(act_docs)

    def _open_project(self) -> None:
        f, _ = QFileDialog.getOpenFileName(
            self, "Abrir projeto", "", "LUNA GUI project (*.json);;Todos (*)"
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
        if not self.cfg.workdir:
            QMessageBox.warning(self, "Workdir", "Defina um workdir antes de salvar.")
            return
        save_to_workdir(self.cfg)
        self.statusBar().showMessage(f"Projeto salvo em {self.cfg.workdir}", 5000)

    def _about(self) -> None:
        QMessageBox.about(
            self, "Sobre o LUNA GUI",
            "<h3>LUNA GUI</h3>"
            "<p>Interface gráfica para o toolkit LUNA — análise de interações "
            "proteína-ligante em larga escala.</p>"
            "<p>LUNA: <a href='https://luna-toolkit.readthedocs.io'>luna-toolkit.readthedocs.io</a></p>"
        )

    def _collect_and_save(self) -> None:
        self.tab_project.collect()
        self.tab_analyses.collect()
        save_to_workdir(self.cfg)

    def _on_run_finished(self) -> None:
        self.tabs.setCurrentWidget(self.tab_results)
        self.tab_results.wd_edit.setText(self.cfg.workdir)
        self.tab_results.load_all()
        self.tab_history.refresh()

    def _apply_loaded_cfg(self, cfg: ProjectConfig) -> None:
        # Replace shared cfg fields and reflect in UI
        for k, v in vars(cfg).items():
            setattr(self.cfg, k, v)
        # Push values back into widgets
        self.tab_project.protein_edit.setText(self.cfg.protein_file)
        self.tab_project.ligand_edit.setText(self.cfg.ligand_file)
        self.tab_project.workdir_edit.setText(self.cfg.workdir)
        if self.cfg.ligand_file:
            try:
                self.tab_project._load_ligands(self.cfg.ligand_file)
                # restore selection
                from PyQt6.QtCore import Qt
                selected = set(self.cfg.selected_ligands)
                lst = self.tab_project.lig_list
                for i in range(lst.count()):
                    it = lst.item(i)
                    name = it.data(Qt.ItemDataRole.UserRole)
                    it.setCheckState(
                        Qt.CheckState.Checked if name in selected else Qt.CheckState.Unchecked
                    )
                self.tab_project._update_count()
            except Exception:
                pass
        a = self.tab_analyses
        a.fp_box.setChecked(self.cfg.out_ifp)
        a.cb_type.setCurrentText(self.cfg.ifp_type)
        a.sp_levels.setValue(self.cfg.ifp_levels)
        a.sp_radius.setValue(self.cfg.ifp_radius)
        a.sp_length.setValue(self.cfg.ifp_length)
        a.cb_bit.setChecked(self.cfg.ifp_bit)
        a.fp_out_edit.setText(self.cfg.ifp_output)
        a.sm_box.setChecked(self.cfg.sim_matrix)
        a.sm_out_edit.setText(self.cfg.sim_matrix_output)
        a.pse_box.setChecked(self.cfg.out_pse)
        a.pse_dir_edit.setText(self.cfg.pse_path)
        a.fbm_box.setChecked(self.cfg.filter_binding_modes)
        a.fbm_edit.setText(self.cfg.binding_modes_cfg)
        a.fork_edit.setText(self.cfg.fork_from)

        # Restore advanced-option widgets
        self.tab_project.cb_waters.setChecked(self.cfg.include_waters)
        a.adv_box.setChecked(self.cfg.force_python_api or bool(self.cfg.inter_config_overrides))
        a.cb_ic_proximal.setChecked(self.cfg.ic_add_proximal)
        a.cb_ic_atom_atom.setChecked(self.cfg.ic_add_atom_atom)
        a.cb_ic_dep.setChecked(self.cfg.ic_add_dependent_inter)
        a.cb_ic_h2o_nt.setChecked(self.cfg.ic_add_h2o_pairs_with_no_target)
        a.cb_ic_self.setChecked(self.cfg.ic_ignore_self_inter)
        ov = self.cfg.inter_config_overrides or {}
        a.sp_hb_da.setValue(ov.get("max_da_dist_hb_inter", 0.0) or 0.0)
        a.sp_hb_dha.setValue(ov.get("min_dha_ang_hb_inter", 0.0) or 0.0)
        a.sp_hydrop.setValue(ov.get("max_cc_dist_hydrop_inter", 0.0) or 0.0)
        a.sp_pistack.setValue(ov.get("max_cc_dist_pi_pi_inter", 0.0) or 0.0)
        if self.cfg.pse_interaction_types:
            from PyQt6.QtCore import Qt
            a.pse_filter_box.setChecked(True)
            selected = set(self.cfg.pse_interaction_types)
            for i in range(a.pse_types_list.count()):
                it = a.pse_types_list.item(i)
                it.setCheckState(
                    Qt.CheckState.Checked if it.text() in selected else Qt.CheckState.Unchecked
                )

        # Auto-populate Results from the loaded workdir (T2.8)
        if self.cfg.workdir and Path(self.cfg.workdir).exists():
            self.tab_results.wd_edit.setText(self.cfg.workdir)
            try:
                self.tab_results.load_all()
            except Exception:
                pass

        self.tabs.setCurrentWidget(self.tab_project)
