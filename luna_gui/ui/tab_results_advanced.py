"""Runtime-aware Results tab with cached analytics and richer fingerprint views."""
from __future__ import annotations

import math
import shutil
import textwrap
from pathlib import Path

import numpy as np

from PyQt6.QtCore import Qt, QTimer
from PyQt6.QtGui import QPixmap
from PyQt6.QtWidgets import (
    QApplication,
    QComboBox,
    QCompleter,
    QDialog,
    QFileDialog,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QListWidget,
    QListWidgetItem,
    QLineEdit,
    QMessageBox,
    QPushButton,
    QScrollArea,
    QTabWidget,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
    QWidget,
)

from ..core.analysis_runtime import (
    generate_filtered_pse_sessions,
    generate_fp_session,
    run_analysis,
    run_fp_dashboard_analysis,
    run_fp_detail_analysis,
    run_residue_matrix,
)
from ..core.pymol_launcher import launch_pse_session
from ..core.project import PROJECT_FILENAME, ProjectConfig
from ..core.report_export import save_pdf_report_isolated
from ..i18n import translate_figure
from ..core.results_analysis import (
    CLASS_UNRELIABLE,
    FP_CLASS_ORDER,
    INTERACTION_COLORS,
    build_complete_heatmap_layers,
    build_fp_analysis_dashboard,
    build_ligand_atom_entry_counts,
    build_ligand_atom_frame_percentages,
    build_trajectory_entry_counts,
    build_trajectory_frame_percentages,
    cluster_rows,
    format_residue_label,
    get_interaction_color,
    interaction_priority_key,
    is_pi_stacking_interaction,
    is_unfavorable_or_repulsive_interaction,
    load_external_fp_labels,
    load_fp_detail_artifact,
    load_analysis_summary,
    load_fp_analysis_artifacts,
    load_ifp_sparse_matrix,
    load_residue_matrix_artifact,
    normalize_fp_class_name,
    resolve_fp_labels_file,
    resolve_fp_random_seed,
    trajectory_frame_number,
    _resolve_external_label_value,
)
from .tab_results_enhanced import (
    HAS_MPL,
    ResultsTab as EnhancedResultsTab,
    _apply_tick_labels,
    _sortable_item,
)
from .info import InfoButton

if HAS_MPL:
    from matplotlib import colormaps
    from matplotlib.colors import BoundaryNorm, ListedColormap
    from matplotlib.patches import Patch, Rectangle


class ResultsTab(EnhancedResultsTab):
    _REPORT_FIGURE_WIDTH_IN = 12.0
    _REPORT_FIGURE_HEIGHT_IN = 6.6
    _REPORT_FIGURE_DPI = 180

    def __init__(self, cfg) -> None:
        super().__init__(cfg)
        self._fp_artifacts: dict[str, dict] = {}
        self._fp_dashboards: dict[tuple[str, str], dict] = {}
        self._generated_fp_session: str = ""
        self._stats_hidden_interactions: set[str] = set()
        self._complete_hidden_interactions: set[str] = set()
        self._stats_legend_pick_cid: int | None = None
        self._stats_legend_artist_labels: dict[object, str] = {}
        self._fp_interaction_hidden_types: set[str] = set()
        self._fp_interaction_legend_pick_cid: int | None = None
        self._fp_interaction_legend_artist_labels: dict[object, str] = {}
        self._fp_manual_feature_ids: dict[str, list[int]] = {}
        self._fullscreen_dialog: QDialog | None = None
        self._workdir_project_cache: dict[str, tuple[float, ProjectConfig | None]] = {}
        self._fp_plot_canvases_ready = False
        self._fp_plot_specs: list[tuple[str, str, QWidget, tuple[float, float], QLabel]] = []
        self._install_fullscreen_button()
        self._install_stats_scope_control()
        self._install_stats_scroll_area()
        self._install_residue_heatmap_scroll_area()
        self._install_similarity_scroll_area()
        self._install_cluster_scroll_area()
        self._install_pse_dynamic_filter()
        self._install_complete_heatmap_tab()
        self._install_fp_analysis_tab()
        self._install_fp_session_tab()
        self._reorder_tabs()

    def _install_fullscreen_button(self) -> None:
        root_layout = self.layout()
        if root_layout is None or root_layout.count() < 2:
            return
        row = root_layout.itemAt(1).layout()
        if row is None:
            return
        btn_fullscreen = QPushButton("Full screen")
        btn_fullscreen.setToolTip("Abre o gráfico atual em uma janela maior para inspeção visual.")
        btn_fullscreen.clicked.connect(self.show_current_chart_fullscreen)
        row.addWidget(btn_fullscreen)
        row.addWidget(InfoButton("Amplia o gráfico atualmente visível. Use para inspecionar mapas de calor e barras sem mudar os dados."))

    def show_current_chart_fullscreen(self) -> None:
        if not HAS_MPL:
            QMessageBox.information(self, "Full screen", "matplotlib não está instalado.")
            return
        figure, default_name = self._current_figure()
        if figure is None or not getattr(figure, "axes", None):
            QMessageBox.information(self, "Sem gráfico", "A aba atual não possui um gráfico para ampliar.")
            return
        wd = self._current_wd() or Path.cwd()
        preview_dir = Path(wd) / "results" / "_fullscreen_previews"
        preview_dir.mkdir(parents=True, exist_ok=True)
        out = preview_dir / f"{_safe_name(default_name or 'grafico')}.png"
        try:
            translate_figure(figure)
            figure.savefig(out, dpi=260, bbox_inches="tight", pad_inches=0.18)
        except Exception as exc:
            QMessageBox.critical(self, "Full screen", str(exc))
            return

        dlg = QDialog(self)
        dlg.setWindowTitle(default_name or "Grafico")
        dlg.resize(1280, 820)
        layout = QVBoxLayout(dlg)
        toolbar = QHBoxLayout()
        btn_zoom_out = QPushButton("Afastar")
        btn_zoom_in = QPushButton("Aproximar")
        btn_zoom_100 = QPushButton("100%")
        btn_fit = QPushButton("Ajustar a tela")
        btn_close = QPushButton("Fechar")
        zoom_label = QLabel("100%")
        zoom_label.setMinimumWidth(58)
        toolbar.addWidget(btn_zoom_out)
        toolbar.addWidget(btn_zoom_in)
        toolbar.addWidget(btn_zoom_100)
        toolbar.addWidget(btn_fit)
        toolbar.addWidget(zoom_label)
        toolbar.addStretch()
        toolbar.addWidget(btn_close)
        layout.addLayout(toolbar)

        original_pixmap = QPixmap(str(out))
        image_label = QLabel()
        image_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        scroll = QScrollArea()
        scroll.setWidget(image_label)
        scroll.setWidgetResizable(False)
        layout.addWidget(scroll, 1)

        zoom_state = {"factor": 1.0}

        def apply_zoom(factor: float) -> None:
            factor = max(0.1, min(6.0, float(factor)))
            zoom_state["factor"] = factor
            width = max(1, int(original_pixmap.width() * factor))
            height = max(1, int(original_pixmap.height() * factor))
            image_label.setPixmap(
                original_pixmap.scaled(
                    width,
                    height,
                    Qt.AspectRatioMode.KeepAspectRatio,
                    Qt.TransformationMode.SmoothTransformation,
                )
            )
            zoom_label.setText(f"{int(round(factor * 100))}%")
            image_label.adjustSize()

        def fit_to_window() -> None:
            viewport = scroll.viewport().size()
            if viewport.width() <= 0 or viewport.height() <= 0:
                apply_zoom(1.0)
                return
            factor = min(
                viewport.width() / max(1, original_pixmap.width()),
                viewport.height() / max(1, original_pixmap.height()),
            )
            apply_zoom(factor)

        btn_zoom_out.clicked.connect(lambda: apply_zoom(zoom_state["factor"] / 1.25))
        btn_zoom_in.clicked.connect(lambda: apply_zoom(zoom_state["factor"] * 1.25))
        btn_zoom_100.clicked.connect(lambda: apply_zoom(1.0))
        btn_fit.clicked.connect(fit_to_window)
        btn_close.clicked.connect(dlg.close)
        self._fullscreen_dialog = dlg
        dlg.showMaximized()
        QTimer.singleShot(0, fit_to_window)

    def _install_stats_scope_control(self) -> None:
        st_layout = self.stats_tab.layout()
        if st_layout is None or st_layout.count() < 2:
            return
        st_ctrl = st_layout.itemAt(1).layout()
        if st_ctrl is None:
            return

        self.cb_stats_scope = QComboBox()
        self.cb_stats_scope.addItem("Totais do projeto", "__all__")
        self.cb_stats_scope.setToolTip(
            "Alterna entre o resumo total do projeto e a distribuicao por um ligante especifico."
        )
        self.cb_stats_scope.currentIndexChanged.connect(self._render_cached_stats_chart)
        st_ctrl.insertWidget(1, QLabel("Visao:"))
        st_ctrl.insertWidget(2, self.cb_stats_scope)
        self.btn_stats_toggle_all = QPushButton("Ligar/desligar todas")
        self.btn_stats_toggle_all.setToolTip(
            "Oculta ou mostra todos os tipos de interação do gráfico de estatísticas."
        )
        self.btn_stats_toggle_all.clicked.connect(self._toggle_all_stats_interactions)
        self.btn_stats_toggle_unfavorable = QPushButton("Ocultar desfav./repulsivas")
        self.btn_stats_toggle_unfavorable.setToolTip(
            "Oculta ou mostra, de uma vez, todas as interações desfavoráveis e repulsivas no gráfico de estatísticas."
        )
        self.btn_stats_toggle_unfavorable.clicked.connect(
            lambda: self._toggle_stats_interaction_group("unfavorable")
        )
        self.btn_stats_toggle_stacking = QPushButton("Ocultar empilhamentos")
        self.btn_stats_toggle_stacking.setToolTip(
            "Oculta ou mostra todos os tipos de pi-stacking/empilhamento no gráfico de estatísticas."
        )
        self.btn_stats_toggle_stacking.clicked.connect(
            lambda: self._toggle_stats_interaction_group("stacking")
        )
        st_ctrl.insertWidget(3, self.btn_stats_toggle_all)
        st_ctrl.insertWidget(4, self.btn_stats_toggle_unfavorable)
        st_ctrl.insertWidget(5, self.btn_stats_toggle_stacking)
        st_ctrl.insertWidget(
            6,
            InfoButton(
                "Estes filtros são apenas visuais: ligam/desligam famílias de interações no gráfico de estatísticas sem alterar os arquivos do projeto."
            ),
        )

    def _install_residue_heatmap_scroll_area(self) -> None:
        if not HAS_MPL or not hasattr(self, "hm_canvas"):
            return
        layout = self.residue_tab.layout()
        if layout is None:
            return
        layout.removeWidget(self.hm_canvas)
        self.hm_scroll = self._build_canvas_scroll_area(self.hm_canvas)
        layout.addWidget(self.hm_scroll, 1)

    def _install_stats_scroll_area(self) -> None:
        if not HAS_MPL or not hasattr(self, "st_canvas"):
            return
        layout = self.stats_tab.layout()
        if layout is None:
            return
        layout.removeWidget(self.st_canvas)
        self.stats_scroll = self._build_canvas_scroll_area(self.st_canvas)
        layout.addWidget(self.stats_scroll, 1)

    def _install_similarity_scroll_area(self) -> None:
        if not HAS_MPL or not hasattr(self, "canvas"):
            return
        layout = self.sim_tab.layout()
        if layout is None:
            return
        layout.removeWidget(self.canvas)
        self.sim_scroll = self._build_canvas_scroll_area(self.canvas)
        layout.addWidget(self.sim_scroll, 1)

    def _install_cluster_scroll_area(self) -> None:
        if not (HAS_MPL and hasattr(self, "cluster_canvas")):
            return
        layout = self.cluster_tab.layout()
        if layout is None:
            return
        layout.removeWidget(self.cluster_canvas)
        self.cluster_scroll = self._build_canvas_scroll_area(self.cluster_canvas)
        layout.insertWidget(2, self.cluster_scroll, 1)

    def _build_canvas_scroll_area(self, canvas) -> QScrollArea:
        scroll = QScrollArea()
        scroll.setWidget(canvas)
        scroll.setWidgetResizable(False)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        scroll.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        return scroll

    def _set_layout_visible(self, layout, visible: bool) -> None:
        if layout is None:
            return
        for index in range(layout.count()):
            item = layout.itemAt(index)
            widget = item.widget()
            child_layout = item.layout()
            if widget is not None:
                widget.setVisible(visible)
            if child_layout is not None:
                self._set_layout_visible(child_layout, visible)

    def _make_collapsible(self, group: QGroupBox, checked: bool = False) -> None:
        group.setCheckable(True)

        def apply(visible: bool) -> None:
            self._set_layout_visible(group.layout(), visible)

        group.toggled.connect(apply)
        group.setChecked(bool(checked))
        apply(group.isChecked())

    def _make_searchable_combo(self, combo: QComboBox, placeholder: str = "") -> None:
        combo.setEditable(True)
        combo.setInsertPolicy(QComboBox.InsertPolicy.NoInsert)
        combo.setMaxVisibleItems(24)
        line_edit = combo.lineEdit()
        if line_edit is not None:
            line_edit.setClearButtonEnabled(True)
            if placeholder:
                line_edit.setPlaceholderText(placeholder)
        completer = QCompleter(combo.model(), combo)
        completer.setCaseSensitivity(Qt.CaseSensitivity.CaseInsensitive)
        completer.setFilterMode(Qt.MatchFlag.MatchContains)
        completer.setCompletionMode(QCompleter.CompletionMode.PopupCompletion)
        completer.activated.connect(lambda text, c=combo: self._select_combo_text(c, str(text)))
        combo.setCompleter(completer)

    def _select_combo_text(self, combo: QComboBox, text: str) -> None:
        idx = combo.findText(text, Qt.MatchFlag.MatchExactly)
        if idx >= 0:
            combo.setCurrentIndex(idx)

    def _combo_current_data(self, combo: QComboBox):
        text = combo.currentText().strip()
        if text:
            idx = combo.findText(text, Qt.MatchFlag.MatchExactly)
            if idx >= 0:
                return combo.itemData(idx)
        return combo.currentData() if combo.count() else None

    @staticmethod
    def _contrast_text_color(color: object) -> str:
        text = str(color or "").strip()
        if not text:
            return "#111111"
        named_dark = {"black", "navy", "blue", "purple", "brown", "darkred", "darkblue", "darkgreen"}
        if text.lower() in named_dark:
            return "#ffffff"
        if text.startswith("#"):
            hex_value = text[1:]
            if len(hex_value) == 3:
                hex_value = "".join(ch * 2 for ch in hex_value)
            if len(hex_value) >= 6:
                try:
                    r = int(hex_value[0:2], 16) / 255.0
                    g = int(hex_value[2:4], 16) / 255.0
                    b = int(hex_value[4:6], 16) / 255.0
                except ValueError:
                    return "#111111"
                luminance = (0.299 * r) + (0.587 * g) + (0.114 * b)
                return "#ffffff" if luminance < 0.46 else "#111111"
        return "#111111"

    def _install_pse_dynamic_filter(self) -> None:
        layout = self.pse_tab.layout()
        if layout is None:
            return

        box = QGroupBox("Filtragens dinâmicas por binding modes")
        box_layout = QVBoxLayout(box)
        help_label = QLabel(
            "Gera novas sessões PyMOL em uma subpasta separada usando regras .cfg de binding modes. "
            "As sessões originais não são alteradas; você pode criar, apagar e recriar filtros quantas vezes precisar."
        )
        help_label.setWordWrap(True)
        help_label.setProperty("muted", True)
        box_layout.addWidget(help_label)

        cfg_row = QHBoxLayout()
        cfg_row.addWidget(QLabel("Regras .cfg:"))
        self.pse_filter_cfg_edit = QLineEdit()
        self.pse_filter_cfg_edit.setPlaceholderText("(binding_modes.cfg)")
        cfg_row.addWidget(self.pse_filter_cfg_edit, 1)
        btn_pick = QPushButton("...")
        btn_pick.clicked.connect(self._pick_pse_filter_cfg)
        btn_editor = QPushButton("Editor visual")
        btn_editor.clicked.connect(self._open_pse_filter_editor)
        cfg_row.addWidget(btn_pick)
        cfg_row.addWidget(btn_editor)
        box_layout.addLayout(cfg_row)

        name_row = QHBoxLayout()
        name_row.addWidget(QLabel("Nome do filtro:"))
        self.pse_filter_name_edit = QLineEdit()
        self.pse_filter_name_edit.setPlaceholderText("ex.: hbonds_ligante_ativo")
        name_row.addWidget(self.pse_filter_name_edit, 1)
        btn_generate = QPushButton("Gerar filtragem")
        btn_generate.clicked.connect(self._generate_pse_filter)
        btn_delete = QPushButton("Apagar filtragem selecionada")
        btn_delete.clicked.connect(self._delete_selected_pse_filter)
        name_row.addWidget(btn_generate)
        name_row.addWidget(btn_delete)
        box_layout.addLayout(name_row)

        self.pse_filter_status = QLabel("-")
        self.pse_filter_status.setWordWrap(True)
        self.pse_filter_status.setProperty("muted", True)
        box_layout.addWidget(self.pse_filter_status)

        self.pse_filter_list = QListWidget()
        self.pse_filter_list.setToolTip("Subpastas de sessões PyMOL geradas por filtragem dinâmica.")
        self.pse_filter_list.currentItemChanged.connect(lambda _current, _previous: self._refresh_selected_pse_files())
        box_layout.addWidget(self.pse_filter_list, 1)

        filtered_files_label = QLabel("Arquivos .pse da filtragem selecionada:")
        filtered_files_label.setProperty("muted", True)
        box_layout.addWidget(filtered_files_label)
        self.pse_filter_files_status = QLabel("Selecione uma filtragem para listar os arquivos.")
        self.pse_filter_files_status.setWordWrap(True)
        self.pse_filter_files_status.setProperty("muted", True)
        box_layout.addWidget(self.pse_filter_files_status)
        self.pse_filter_files_list = QListWidget()
        self.pse_filter_files_list.setToolTip("Arquivos .pse que passaram pela filtragem selecionada.")
        self.pse_filter_files_list.itemDoubleClicked.connect(
            lambda item: self._open_pse_path(str(item.data(Qt.ItemDataRole.UserRole)))
        )
        box_layout.addWidget(self.pse_filter_files_list, 1)
        layout.insertWidget(2, box)

    def _pick_pse_filter_cfg(self) -> None:
        f, _ = QFileDialog.getOpenFileName(self, "Arquivo de binding modes", "", "Config (*.cfg);;Todos (*)")
        if f:
            self.pse_filter_cfg_edit.setText(f)

    def _open_pse_filter_editor(self) -> None:
        from .binding_modes_editor import BindingModesEditor

        dlg = BindingModesEditor(self, initial_path=self.pse_filter_cfg_edit.text().strip())
        if dlg.exec() and dlg.path:
            self.pse_filter_cfg_edit.setText(dlg.path)

    @staticmethod
    def _safe_filter_name(value: str) -> str:
        import re

        return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value or "")).strip("_") or "filtered"

    def _pse_filtered_root(self, wd: Path) -> Path:
        return wd / "results" / "pse_filtered"

    def _refresh_pse_filter_list(self, wd: Path | None = None, select_path: Path | None = None) -> None:
        if not hasattr(self, "pse_filter_list"):
            return
        if wd is None:
            wd = self._current_wd()
        selected_path = str(select_path) if select_path is not None else ""
        if not selected_path and self.pse_filter_list.currentItem():
            selected_path = str(self.pse_filter_list.currentItem().data(Qt.ItemDataRole.UserRole))
        self.pse_filter_list.clear()
        if hasattr(self, "pse_filter_files_list"):
            self.pse_filter_files_list.clear()
        if not wd:
            return
        root = self._pse_filtered_root(wd)
        if not root.exists():
            return
        selected_item = None
        for folder in sorted(path for path in root.iterdir() if path.is_dir()):
            pse_count = len(list(folder.glob("*.pse")))
            item = QListWidgetItem(f"{folder.name} ({pse_count} PSE)")
            item.setData(Qt.ItemDataRole.UserRole, str(folder))
            self.pse_filter_list.addItem(item)
            if selected_path and str(folder) == selected_path:
                selected_item = item
        if selected_item is not None:
            self.pse_filter_list.setCurrentItem(selected_item)
        elif self.pse_filter_list.count():
            self.pse_filter_list.setCurrentRow(self.pse_filter_list.count() - 1)
        self._refresh_selected_pse_files()

    def _refresh_selected_pse_files(self) -> None:
        if not hasattr(self, "pse_filter_files_list"):
            return
        self.pse_filter_files_list.clear()
        item = self.pse_filter_list.currentItem() if hasattr(self, "pse_filter_list") else None
        if not item:
            self.pse_filter_files_status.setText("Selecione uma filtragem para listar os arquivos.")
            return
        folder = Path(str(item.data(Qt.ItemDataRole.UserRole)))
        files = sorted(path for path in folder.glob("*.pse") if path.is_file())
        self.pse_filter_files_status.setText(f"{len(files)} arquivos .pse nesta filtragem.")
        for file_path in files:
            file_item = QListWidgetItem(file_path.name)
            file_item.setData(Qt.ItemDataRole.UserRole, str(file_path))
            self.pse_filter_files_list.addItem(file_item)

    def _open_pse_path(self, path: str) -> None:
        if not path:
            return
        try:
            launch_pse_session(path, self.py_exe)
        except Exception as exc:
            QMessageBox.critical(self, "Erro ao abrir PyMOL", str(exc))

    def _generate_pse_filter(self) -> None:
        wd = self._current_wd()
        if not wd:
            return
        if not self.py_exe:
            QMessageBox.warning(self, "luna-env", "LUNA não detectado. Verifique a aba Setup.")
            return
        cfg_path = self.pse_filter_cfg_edit.text().strip()
        if not cfg_path:
            QMessageBox.warning(self, "Binding modes", "Escolha ou crie um arquivo .cfg primeiro.")
            return
        filter_name = self._safe_filter_name(self.pse_filter_name_edit.text().strip() or Path(cfg_path).stem)
        filtered_root = self._pse_filtered_root(wd)
        output_dir = filtered_root / filter_name
        suffix = 2
        while output_dir.exists():
            output_dir = filtered_root / f"{filter_name}_{suffix}"
            suffix += 1
        self.pse_filter_status.setText("Gerando sessões filtradas...")
        self.pse_filter_status.repaint()
        result = generate_filtered_pse_sessions(self.py_exe, str(wd), cfg_path, str(output_dir))
        if "error" in result:
            self.pse_filter_status.setText("Falha ao gerar filtragem.")
            QMessageBox.critical(self, "Erro na filtragem PyMOL", result["error"])
            return
        if result.get("fallback"):
            self.pse_filter_status.setText(
                f"{result.get('created', 0)} arquivos PSE copiados por filtragem cacheada; "
                f"{result.get('matched_interactions', 0)} ocorrências compatíveis."
            )
        else:
            self.pse_filter_status.setText(
                f"{result.get('created', 0)} sessões geradas; "
                f"{result.get('matched_interactions', 0)} interações filtradas."
            )
        warnings = result.get("warnings") or []
        if warnings:
            self.pse_filter_status.setText(self.pse_filter_status.text() + f" Avisos: {len(warnings)}.")
        self._refresh_pse_filter_list(wd, output_dir)
        self._load_pse(wd)

    def _delete_selected_pse_filter(self) -> None:
        item = self.pse_filter_list.currentItem() if hasattr(self, "pse_filter_list") else None
        if not item:
            QMessageBox.information(self, "Filtragem PyMOL", "Selecione uma filtragem para apagar.")
            return
        wd = self._current_wd()
        if not wd:
            return
        folder = Path(str(item.data(Qt.ItemDataRole.UserRole)))
        root = self._pse_filtered_root(wd).resolve()
        try:
            target = folder.resolve()
        except Exception:
            QMessageBox.warning(self, "Filtragem PyMOL", "Caminho inválido.")
            return
        if root not in target.parents:
            QMessageBox.warning(self, "Filtragem PyMOL", "A GUI só pode apagar filtragens dentro do workdir atual.")
            return
        if QMessageBox.question(
            self,
            "Apagar filtragem",
            f"Apagar a pasta de filtragem?\n{target}",
        ) != QMessageBox.StandardButton.Yes:
            return
        shutil.rmtree(target, ignore_errors=True)
        self._refresh_pse_filter_list(wd)
        self._load_pse(wd)

    def _resize_canvas(self, fig, canvas, width_in: float, height_in: float) -> None:
        fig.set_dpi(120)
        fig.set_size_inches(width_in, height_in, forward=True)
        canvas.setMinimumSize(int(width_in * fig.dpi), int(height_in * fig.dpi))

    def _resize_heatmap_canvas(self, fig, canvas, n_residues: int, n_entries: int) -> None:
        width = max(9.0, 2.5 + (0.34 * max(1, n_residues)))
        height = max(5.4, 2.5 + (0.12 * self._render_entry_count(n_entries)))
        self._resize_canvas(fig, canvas, width, height)

    def _residue_xticklabels(self, residues: list[str]) -> list[str]:
        return [format_residue_label(residue) for residue in residues]

    @staticmethod
    def _display_ligand_name(entry_name: str) -> str:
        text = str(entry_name or "").strip()
        if ":" in text:
            text = text.split(":")[-1].strip()
        suffixes = (
            "_ligand",
            "-ligand",
            "_lig",
            "-lig",
            "_LIGAND",
            "-LIGAND",
            "_LIG",
            "-LIG",
        )
        changed = True
        while changed and text:
            changed = False
            for suffix in suffixes:
                if text.endswith(suffix):
                    text = text[: -len(suffix)]
                    changed = True
                    break
        return text or str(entry_name or "").strip()

    @staticmethod
    def _format_shell_levels(values, breakdown=None) -> str:
        source = []
        if isinstance(values, (list, tuple, set)):
            source.extend(values)
        elif values is not None and str(values).strip():
            source.append(values)
        if not source and isinstance(breakdown, dict):
            source.extend(breakdown.keys())
        keys = {str(value).strip() for value in source if str(value).strip()}
        ordered = sorted(
            keys,
            key=lambda value: (
                not str(value).lstrip("-").isdigit(),
                int(value) if str(value).lstrip("-").isdigit() else str(value),
            ),
        )
        return ", ".join(ordered)

    @staticmethod
    def _fp_feature_label(feature: dict) -> str:
        feature_key = str((feature or {}).get("feature_key", "") or "").strip()
        if feature_key:
            return feature_key
        feature_id = str((feature or {}).get("feature_id", "") or "").strip()
        assigned_level = str((feature or {}).get("assigned_level", "") or "").strip()
        if feature_id and assigned_level:
            return f"{feature_id}_{assigned_level}"
        return feature_id

    @staticmethod
    def _coerce_float(value: object) -> float | None:
        text = str(value or "").strip()
        if not text:
            return None
        try:
            return float(text.replace(",", "."))
        except Exception:
            return None

    @staticmethod
    def _label_sort_key(value: object) -> tuple[int, float | str]:
        text = str(value or "").strip()
        try:
            return (0, float(text.replace(",", ".")))
        except Exception:
            return (1, text.casefold())

    @staticmethod
    def _label_group_colors(groups: list[str]) -> dict[str, str]:
        palette = [
            "#0f766e",
            "#b45309",
            "#2563eb",
            "#be123c",
            "#7c3aed",
            "#15803d",
            "#c2410c",
            "#0369a1",
            "#a21caf",
            "#4d7c0f",
        ]
        unique = sorted({str(group) for group in groups if str(group).strip()}, key=lambda value: value.casefold())
        return {group: palette[idx % len(palette)] for idx, group in enumerate(unique)}

    def _workdir_path_no_warning(self) -> Path | None:
        wd = str(self.wd_edit.text().strip() or getattr(self.cfg, "workdir", "") or "").strip()
        if not wd:
            return None
        path = Path(wd)
        return path if path.exists() else None

    @staticmethod
    def _same_existing_path(left: Path | str | None, right: Path | str | None) -> bool:
        if not left or not right:
            return False
        try:
            return Path(left).resolve() == Path(right).resolve()
        except Exception:
            return str(left) == str(right)

    def _workdir_project_cfg(self) -> ProjectConfig | None:
        wd = self._workdir_path_no_warning()
        if wd is None:
            return None
        cfg_path = wd / PROJECT_FILENAME
        if not cfg_path.exists():
            return None
        cache_key = str(cfg_path)
        try:
            mtime = cfg_path.stat().st_mtime
        except OSError:
            return None
        cached = self._workdir_project_cache.get(cache_key)
        if cached and cached[0] == mtime:
            return cached[1]
        try:
            loaded = ProjectConfig.load(cfg_path)
        except Exception:
            loaded = None
        self._workdir_project_cache[cache_key] = (mtime, loaded)
        return loaded

    def _workdir_cfg_is_active(self) -> bool:
        wd = self._workdir_path_no_warning()
        saved_cfg = self._workdir_project_cfg()
        if wd is None or saved_cfg is None:
            return False
        current_wd = str(getattr(self.cfg, "workdir", "") or "").strip()
        return not current_wd or not self._same_existing_path(wd, current_wd)

    def _results_label_settings(self) -> tuple[str, str, str, str]:
        saved_cfg = self._workdir_project_cfg()
        if saved_cfg is not None and self._workdir_cfg_is_active():
            saved_labels_csv = str(getattr(saved_cfg, "fp_labels_csv", "") or "").strip()
            if saved_labels_csv:
                return (
                    saved_labels_csv,
                    str(getattr(saved_cfg, "fp_labels_id_column", "") or "").strip(),
                    str(getattr(saved_cfg, "fp_labels_column", "") or "").strip(),
                    str(getattr(saved_cfg, "fp_label_task", "regression") or "regression").strip().lower(),
                )

        labels_csv = str(getattr(self.cfg, "fp_labels_csv", "") or "").strip()
        labels_id_column = str(getattr(self.cfg, "fp_labels_id_column", "") or "").strip()
        labels_column = str(getattr(self.cfg, "fp_labels_column", "") or "").strip()
        task = str(getattr(self.cfg, "fp_label_task", "regression") or "regression").strip().lower()
        if labels_csv or saved_cfg is None:
            return labels_csv, labels_id_column, labels_column, task

        return (
            str(getattr(saved_cfg, "fp_labels_csv", "") or "").strip(),
            str(getattr(saved_cfg, "fp_labels_id_column", "") or "").strip(),
            str(getattr(saved_cfg, "fp_labels_column", "") or "").strip(),
            str(getattr(saved_cfg, "fp_label_task", "regression") or "regression").strip().lower(),
        )

    def _results_trajectory_mode(self) -> bool:
        saved_cfg = self._workdir_project_cfg()
        if saved_cfg is not None and self._workdir_cfg_is_active():
            return bool(getattr(saved_cfg, "trajectory_analysis", False))
        return bool(getattr(self.cfg, "trajectory_analysis", False))

    def _entry_label_context(self, entries: list[str]) -> dict:
        labels_csv, labels_id_column, labels_column, task = self._results_label_settings()
        wd = self._current_wd()
        empty = {"task": "", "values": {}, "numeric": {}, "groups": {}, "colors": {}}
        if not wd or not labels_csv:
            return empty

        label_path = resolve_fp_labels_file(wd, labels_csv)
        if label_path is None:
            return empty

        try:
            label_map = load_external_fp_labels(
                label_path,
                label_column=labels_column,
                id_column=labels_id_column,
            )
        except Exception:
            return empty

        values: dict[str, str] = {}
        numeric: dict[str, float] = {}
        groups: dict[str, str] = {}
        for entry_name in entries:
            raw_value = _resolve_external_label_value(label_map, entry_name)
            if raw_value == "":
                continue
            entry_key = str(entry_name)
            values[entry_key] = str(raw_value)
            parsed = self._coerce_float(raw_value)
            if parsed is not None:
                numeric[entry_key] = float(parsed)
            groups[entry_key] = str(raw_value)

        if not values:
            return empty
        if task == "regression" and len(numeric) != len(values):
            task = "classification"
        elif task not in {"regression", "classification"}:
            task = "regression" if len(numeric) == len(values) else "classification"
        return {
            "task": task,
            "values": values,
            "numeric": numeric,
            "groups": groups,
            "colors": self._label_group_colors(list(groups.values())),
        }

    def _entry_activity_values(self, entries: list[str]) -> dict[str, float]:
        return dict(self._entry_label_context(entries).get("numeric", {}))

    def _order_entries_by_activity(self, entries: list[str]) -> tuple[list[str], list[int]]:
        if self._is_trajectory_mode():
            ranked = sorted(
                range(len(entries)),
                key=lambda idx: (
                    0 if trajectory_frame_number(entries[idx]) is not None else 1,
                    -(trajectory_frame_number(entries[idx]) or -1),
                    str(entries[idx]).lower(),
                    idx,
                ),
            )
            return [entries[idx] for idx in ranked], ranked
        context = self._entry_label_context(entries)
        task = str(context.get("task") or "")
        numeric = dict(context.get("numeric") or {})
        groups = dict(context.get("groups") or {})
        if task == "regression" and numeric:
            ranked = sorted(
                range(len(entries)),
                key=lambda idx: (
                    0 if entries[idx] in numeric else 1,
                    -float(numeric.get(entries[idx], 0.0)),
                    idx,
                ),
            )
            return [entries[idx] for idx in ranked], ranked
        if task == "classification" and groups:
            group_order = sorted({str(value) for value in groups.values()}, key=self._label_sort_key)
            group_rank = {group: idx for idx, group in enumerate(group_order)}
            ranked = sorted(
                range(len(entries)),
                key=lambda idx: (
                    0 if entries[idx] in groups else 1,
                    group_rank.get(str(groups.get(entries[idx], "")), len(group_rank)),
                    idx,
                ),
            )
            return [entries[idx] for idx in ranked], ranked
        if not numeric:
            return list(entries), list(range(len(entries)))
        return list(entries), list(range(len(entries)))

    def _color_ticklabels_by_entry_group(self, ax, entries: list[str], axis: str = "y", offset: float = 0.0) -> None:
        context = self._entry_label_context(entries)
        if str(context.get("task") or "") != "classification":
            return
        groups = dict(context.get("groups") or {})
        colors = dict(context.get("colors") or {})
        if not groups:
            return

        def _entry_for_position(pos: float) -> str | None:
            idx = int(round(pos - offset))
            if 0 <= idx < len(entries):
                return entries[idx]
            return None

        if axis in {"x", "both"}:
            for tick in ax.get_xticklabels():
                entry = _entry_for_position(float(tick.get_position()[0]))
                group = groups.get(str(entry)) if entry is not None else None
                if group:
                    tick.set_color(colors.get(str(group), "#2e241f"))
                    tick.set_fontweight("bold")
        if axis in {"y", "both"}:
            for tick in ax.get_yticklabels():
                entry = _entry_for_position(float(tick.get_position()[1]))
                group = groups.get(str(entry)) if entry is not None else None
                if group:
                    tick.set_color(colors.get(str(group), "#2e241f"))
                    tick.set_fontweight("bold")

    def _current_hbond_override(self) -> float:
        overrides = getattr(self.cfg, "inter_config_overrides", {}) or {}
        try:
            return float(overrides.get("max_da_dist_hb_inter", 0.0) or 0.0)
        except Exception:
            return 0.0

    def _build_interaction_colors(self, present_types: list[str]) -> list[str]:
        colors = ["#f4efe6"]
        fallback_pool: list[str] = []
        for cmap_name in ("tab20", "tab20b", "tab20c", "Set3"):
            cmap = colormaps.get_cmap(cmap_name)
            fallback_pool.extend(
                [
                    "#{:02x}{:02x}{:02x}".format(
                        int(round(rgba[0] * 255)),
                        int(round(rgba[1] * 255)),
                        int(round(rgba[2] * 255)),
                    )
                    for rgba in getattr(cmap, "colors", [cmap(i / max(1, cmap.N - 1)) for i in range(cmap.N)])
                ]
            )
        used = {value.lower() for value in INTERACTION_COLORS.values()}
        fallback_index = 0
        for name in present_types:
            mapped = get_interaction_color(name)
            if mapped and mapped != "#2f7f83":
                colors.append(mapped)
                continue
            mapped = INTERACTION_COLORS.get(name)
            if mapped:
                colors.append(mapped)
                continue
            while fallback_index < len(fallback_pool) and fallback_pool[fallback_index].lower() in used:
                fallback_index += 1
            mapped = fallback_pool[fallback_index] if fallback_index < len(fallback_pool) else "#2f7f83"
            fallback_index += 1
            used.add(mapped.lower())
            colors.append(mapped)
        return colors

    def _hb_warning_message(self, interaction_types: list[str]) -> str:
        hb_override = self._current_hbond_override()
        if hb_override <= 0.0:
            return ""
        has_hbond = "Hydrogen bond" in interaction_types
        has_weak_hbond = "Weak hydrogen bond" in interaction_types
        if hb_override <= 0.5 and not has_hbond and has_weak_hbond:
            return (
                "Aviso: a configuracao atual usa max_da_dist_hb_inter="
                f"{hb_override:.2f} A, um valor muito restritivo para Hydrogen bond. "
                "Isso tende a eliminar H-bonds classicos e deixar apenas Weak hydrogen bond. "
                "Volte a 0.0/(padrão = 3.9 A) na aba 3.Análises para restaurar o critério padrão do LUNA."
            )
        return ""

    def _load_sim_matrix(self, wd: Path) -> None:
        super()._load_sim_matrix(wd)
        if not (HAS_MPL and getattr(self, "_sim_matrix", None) is not None and getattr(self, "_sim_labels", None)):
            return

        labels = list(self._sim_labels)
        ordered_labels, order = self._order_entries_by_activity(labels)
        if order and order != list(range(len(labels))):
            self._sim_labels = ordered_labels
            self._sim_matrix = self._sim_matrix[np.ix_(order, order)]

        label_count = max(1, len(self._sim_labels))
        rendered_count = self._render_entry_count(label_count)
        width_in = max(8.4, 3.2 + ((0.5 / 2.54) * rendered_count))
        height_in = max(7.0, 3.0 + ((0.5 / 2.54) * rendered_count))
        self._resize_canvas(self.fig, self.canvas, width_in=width_in, height_in=height_in)
        self.fig.clear()
        ax = self.fig.add_subplot(111)
        im = ax.imshow(self._sim_matrix, cmap="viridis", aspect="auto", vmin=0, vmax=1)
        ax.set_title("Similaridade de Tanimoto")
        _apply_tick_labels(
            ax,
            [self._display_ligand_name(label) for label in self._sim_labels],
            axis="both",
            ligand_axis=True,
        )
        self._color_ticklabels_by_entry_group(ax, list(self._sim_labels), axis="both")
        self.fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        self.fig.tight_layout()
        self.canvas.draw()
        self.refresh_clusters()

    def _render_cluster_chart(self, result) -> None:
        super()._render_cluster_chart(result)
        if not (HAS_MPL and result and hasattr(self, "cluster_fig")):
            return
        for ax in self.cluster_fig.axes:
            if ax.get_title() == "Matriz reordenada por cluster":
                self._color_ticklabels_by_entry_group(ax, list(result.ordered_labels), axis="both")
        self.cluster_canvas.draw()

    def _install_complete_heatmap_tab(self) -> None:
        self.complete_heatmap_tab = QWidget()
        layout = QVBoxLayout(self.complete_heatmap_tab)
        help_label = QLabel(
            "Mostra, em cada par ligante x resíduo, todas as interações presentes em faixas coloridas. "
            "As cores indicam a classe da interação, não a contagem."
        )
        help_label.setWordWrap(True)
        help_label.setProperty("muted", True)
        layout.addWidget(help_label)

        ctrl = QHBoxLayout()
        btn = QPushButton("Atualizar mapas de calor")
        btn.clicked.connect(self.compute_residue_matrix)
        ctrl.addWidget(btn)
        self.btn_complete_toggle_all = QPushButton("Ligar/desligar todas")
        self.btn_complete_toggle_all.setToolTip(
            "Oculta ou mostra todos os tipos de interação no mapa de calor completo."
        )
        self.btn_complete_toggle_all.clicked.connect(self._toggle_all_complete_heatmap_interactions)
        self.btn_complete_toggle_unfavorable = QPushButton("Ocultar desfav./repulsivas")
        self.btn_complete_toggle_unfavorable.setToolTip(
            "Oculta ou mostra todas as interações desfavoráveis e repulsivas no mapa de calor completo."
        )
        self.btn_complete_toggle_unfavorable.clicked.connect(
            lambda: self._toggle_complete_heatmap_group("unfavorable")
        )
        self.btn_complete_toggle_stacking = QPushButton("Ocultar empilhamentos")
        self.btn_complete_toggle_stacking.setToolTip(
            "Oculta ou mostra todos os tipos de pi-stacking/empilhamento no mapa de calor completo."
        )
        self.btn_complete_toggle_stacking.clicked.connect(
            lambda: self._toggle_complete_heatmap_group("stacking")
        )
        ctrl.addWidget(self.btn_complete_toggle_all)
        ctrl.addWidget(self.btn_complete_toggle_unfavorable)
        ctrl.addWidget(self.btn_complete_toggle_stacking)
        ctrl.addWidget(
            InfoButton(
                "Filtros visuais do mapa de calor completo. Eles ajudam a remover interações desfavoráveis/repulsivas ou empilhamentos para comparar padrões."
            )
        )
        self.hm_all_status = QLabel("-")
        self.hm_all_status.setProperty("muted", True)
        ctrl.addWidget(self.hm_all_status, 1)
        layout.addLayout(ctrl)

        if HAS_MPL:
            from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
            from matplotlib.figure import Figure

            self.hm_all_fig = Figure(figsize=(8.2, 5.4))
            self.hm_all_canvas = FigureCanvas(self.hm_all_fig)
            self.hm_all_scroll = self._build_canvas_scroll_area(self.hm_all_canvas)
            layout.addWidget(self.hm_all_scroll, 1)
        else:
            layout.addWidget(QLabel("matplotlib não está instalado."))

    def _install_fp_analysis_tab(self) -> None:
        self.fp_analysis_tab = QWidget()
        outer_layout = QVBoxLayout(self.fp_analysis_tab)
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        content = QWidget()
        layout = QVBoxLayout(content)
        help_label = QLabel(
            "Classifica os fingerprints com base nos shells reais do LUNA, "
            "atribui classes confiáveis via limiar por z-score e mostra as features mais relevantes."
        )
        help_label.setWordWrap(True)
        help_label.setProperty("muted", True)
        layout.addWidget(help_label)

        ctrl = QHBoxLayout()
        btn = QPushButton("Carregar análises de FP")
        btn.clicked.connect(self.compute_fp_analyses)
        ctrl.addWidget(btn)
        ctrl.addWidget(QLabel("Tipo:"))
        self.cb_fp_analysis_type = QComboBox()
        self.cb_fp_analysis_type.currentIndexChanged.connect(self._on_fp_analysis_controls_changed)
        ctrl.addWidget(self.cb_fp_analysis_type)
        ctrl.addWidget(QLabel("Algoritmo:"))
        self.cb_fp_algorithm = QComboBox()
        self.cb_fp_algorithm.addItem("GradientBoosting", "gradient_boosting")
        self.cb_fp_algorithm.addItem("ExtraTrees", "extra_trees")
        self.cb_fp_algorithm.setCurrentIndex(0)
        self.cb_fp_algorithm.currentIndexChanged.connect(self._on_fp_analysis_controls_changed)
        ctrl.addWidget(self.cb_fp_algorithm)
        ctrl.addWidget(QLabel("Corte p-value:"))
        self.cb_fp_pvalue_cutoff = QComboBox()
        self.cb_fp_pvalue_cutoff.addItem("0.05 Flexível", 0.05)
        self.cb_fp_pvalue_cutoff.addItem("0.02 Médio", 0.02)
        self.cb_fp_pvalue_cutoff.addItem("0.01 Conservador", 0.01)
        self.cb_fp_pvalue_cutoff.setCurrentIndex(2)
        self.cb_fp_pvalue_cutoff.currentIndexChanged.connect(self._on_fp_analysis_controls_changed)
        ctrl.addWidget(self.cb_fp_pvalue_cutoff)
        self.fp_analysis_status = QLabel("-")
        self.fp_analysis_status.setProperty("muted", True)
        ctrl.addWidget(self.fp_analysis_status, 1)
        layout.addLayout(ctrl)

        manual_row = QHBoxLayout()
        manual_row.addWidget(QLabel("Fingerprints relevantes:"))
        self.fp_manual_features_edit = QLineEdit()
        self.fp_manual_features_edit.setPlaceholderText("IDs separados por vírgula, espaço ou quebra de linha")
        self.fp_manual_features_edit.setToolTip(
            "Permite substituir a lista calculada pelo programa por fingerprints escolhidos manualmente. "
            "A escolha atualiza a tabela e os gráficos de análises FP."
        )
        manual_row.addWidget(self.fp_manual_features_edit, 1)
        btn_manual_apply = QPushButton("Aplicar seleção")
        btn_manual_apply.setToolTip("Usa os IDs informados como fingerprints relevantes para os gráficos abaixo.")
        btn_manual_apply.clicked.connect(self._apply_manual_fp_selection)
        btn_manual_from_table = QPushButton("Usar linhas selecionadas")
        btn_manual_from_table.setToolTip("Copia os IDs das linhas selecionadas da tabela para a seleção manual.")
        btn_manual_from_table.clicked.connect(self._use_selected_fp_rows_as_manual_selection)
        btn_manual_reset = QPushButton("Restabelecer calculado")
        btn_manual_reset.setToolTip("Volta ao critério automático de importância calculado pelo programa.")
        btn_manual_reset.clicked.connect(self._reset_manual_fp_selection)
        manual_row.addWidget(btn_manual_apply)
        manual_row.addWidget(btn_manual_from_table)
        manual_row.addWidget(btn_manual_reset)
        manual_row.addWidget(
            InfoButton(
                "Substitui temporariamente os fingerprints importantes calculados por IDs escolhidos pelo usuário. Todos os gráficos de análises FP passam a usar essa seleção."
            )
        )
        layout.addLayout(manual_row)

        self.fp_analysis_summary = QLabel("-")
        self.fp_analysis_summary.setWordWrap(True)
        self.fp_analysis_summary.setProperty("muted", True)

        self.fp_analysis_formula = QLabel("-")
        self.fp_analysis_formula.setWordWrap(True)
        self.fp_analysis_formula.setProperty("muted", True)

        self.fp_analysis_method = QLabel("-")
        self.fp_analysis_method.setWordWrap(True)
        self.fp_analysis_method.setProperty("muted", True)

        self.fp_analysis_active_context = QLabel("-")
        self.fp_analysis_active_context.setWordWrap(True)
        self.fp_analysis_active_context.setProperty("muted", True)

        self.fp_analysis_method_box = QGroupBox("Descrição do método: classificação e seleção")
        method_layout = QVBoxLayout(self.fp_analysis_method_box)
        method_layout.addWidget(self.fp_analysis_summary)
        method_layout.addWidget(self.fp_analysis_formula)
        method_layout.addWidget(self.fp_analysis_method)
        method_layout.addWidget(self.fp_analysis_active_context)
        layout.addWidget(self.fp_analysis_method_box)
        self._make_collapsible(self.fp_analysis_method_box, False)

        self.fp_analysis_table = QTableWidget()
        self.fp_analysis_table.setColumnCount(14)
        self.fp_analysis_table.setHorizontalHeaderLabels(
            [
                "Feature",
                "Moléculas",
                "Cobertura (%)",
                "Classe prevalente (%)",
                "Z-score classe",
                "Classe atribuída",
                "Importância",
                "Z-score Importance",
                "p-value",
                "Colisões",
                "Nível assinado",
                "Níveis shell",
                "Níveis colisão",
                "Perfil da base",
            ]
        )
        self.fp_analysis_table.setSelectionBehavior(QTableWidget.SelectionBehavior.SelectRows)
        self.fp_analysis_table.setSelectionMode(QTableWidget.SelectionMode.ExtendedSelection)
        self.fp_analysis_table.itemSelectionChanged.connect(self._sync_fp_session_from_table)
        self.fp_analysis_table.horizontalHeader().setStretchLastSection(True)
        self.fp_analysis_table.setSortingEnabled(True)
        layout.addWidget(self.fp_analysis_table, 1)

        if HAS_MPL:
            self.fp_plot_tabs = QTabWidget()
            plot_specs = [
                ("fp_class_tab", "fp_class_fig", "fp_class_canvas", "Classes importantes", (7.2, 4.4)),
                ("fp_assign_tab", "fp_assign_fig", "fp_assign_canvas", "Frequência por classe", (7.2, 5.2)),
                ("fp_cover_tab", "fp_cover_fig", "fp_cover_canvas", "Cobertura e importância", (7.2, 5.2)),
                ("fp_heatmap_tab", "fp_heatmap_fig", "fp_heatmap_canvas", "Mapa de calor de importância", (9.6, 7.2)),
                ("fp_interaction_assign_tab", "fp_interaction_assign_fig", "fp_interaction_assign_canvas", "Frequência de interações", (9.2, 5.4)),
                ("fp_interaction_tab", "fp_interaction_fig", "fp_interaction_canvas", "Interações prevalentes", (9.2, 5.4)),
                ("fp_interaction_heatmap_tab", "fp_interaction_heatmap_fig", "fp_interaction_heatmap_canvas", "Mapa de calor de interações", (9.6, 7.2)),
            ]
            for tab_attr, fig_attr, canvas_attr, label, size in plot_specs:
                tab = QWidget()
                tab_layout = QVBoxLayout(tab)
                placeholder = QLabel("O gráfico será criado quando uma análise FP for carregada.")
                placeholder.setAlignment(Qt.AlignmentFlag.AlignCenter)
                placeholder.setProperty("muted", True)
                tab_layout.addWidget(placeholder, 1)
                setattr(self, tab_attr, tab)
                self._fp_plot_specs.append((fig_attr, canvas_attr, tab, size, placeholder))
                self.fp_plot_tabs.addTab(tab, label)

            layout.addWidget(self.fp_plot_tabs, 1)
        else:
            self.fp_plot_tabs = None
            layout.addWidget(QLabel("matplotlib não está instalado."))
        scroll.setWidget(content)
        outer_layout.addWidget(scroll, 1)

    def _ensure_fp_plot_canvases(self) -> None:
        if not HAS_MPL or self.fp_plot_tabs is None or self._fp_plot_canvases_ready:
            return

        from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
        from matplotlib.figure import Figure

        for fig_attr, canvas_attr, tab, size, placeholder in self._fp_plot_specs:
            tab_layout = tab.layout()
            tab_layout.removeWidget(placeholder)
            placeholder.deleteLater()
            figure = Figure(figsize=size)
            canvas = FigureCanvas(figure)
            setattr(self, fig_attr, figure)
            setattr(self, canvas_attr, canvas)
            tab_layout.addWidget(self._build_canvas_scroll_area(canvas), 1)
        self._fp_plot_canvases_ready = True

    def _install_fp_session_tab(self) -> None:
        self.fp_session_tab = QWidget()
        layout = QVBoxLayout(self.fp_session_tab)
        help_label = QLabel(
            "Gera uma sessão PyMOL para o fingerprint escolhido, recuperando os shells "
            "que originaram aquele bit para um ligante específico."
        )
        help_label.setWordWrap(True)
        help_label.setProperty("muted", True)
        layout.addWidget(help_label)

        ctrl = QHBoxLayout()
        ctrl.addWidget(QLabel("Tipo:"))
        self.cb_fp_session_type = QComboBox()
        self._make_searchable_combo(self.cb_fp_session_type, "Digite para buscar o tipo")
        self.cb_fp_session_type.currentIndexChanged.connect(self._populate_fp_session_features)
        ctrl.addWidget(self.cb_fp_session_type)
        ctrl.addWidget(QLabel("Fingerprint:"))
        self.cb_fp_session_feature = QComboBox()
        self._make_searchable_combo(self.cb_fp_session_feature, "Digite para buscar o fingerprint")
        self.cb_fp_session_feature.currentIndexChanged.connect(self._populate_fp_session_entries)
        ctrl.addWidget(self.cb_fp_session_feature, 1)
        ctrl.addWidget(QLabel("Molécula:"))
        self.cb_fp_session_entry = QComboBox()
        self._make_searchable_combo(self.cb_fp_session_entry, "Digite para buscar a molécula")
        self.cb_fp_session_entry.currentIndexChanged.connect(self._update_fp_session_info)
        ctrl.addWidget(self.cb_fp_session_entry, 1)
        layout.addLayout(ctrl)

        self.fp_session_info = QLabel("-")
        self.fp_session_info.setWordWrap(True)
        self.fp_session_info.setProperty("muted", True)
        layout.addWidget(self.fp_session_info)

        btn_row = QHBoxLayout()
        btn_generate = QPushButton("Gerar sessão")
        btn_generate.clicked.connect(self._generate_fp_session)
        btn_open = QPushButton("Abrir sessão")
        btn_open.clicked.connect(self._open_selected_fp_session)
        btn_row.addWidget(btn_generate)
        btn_row.addWidget(btn_open)
        self.fp_session_status = QLabel("-")
        self.fp_session_status.setProperty("muted", True)
        btn_row.addWidget(self.fp_session_status, 1)
        layout.addLayout(btn_row)

        self.fp_session_list = QListWidget()
        self.fp_session_list.itemDoubleClicked.connect(lambda _item: self._open_selected_fp_session())
        layout.addWidget(self.fp_session_list, 1)

    def _reorder_tabs(self) -> None:
        ordered_tabs = [
            (self.stats_tab, "Estatísticas"),
            (self.residue_tab, "Mapa de calor por tipo"),
            (self.complete_heatmap_tab, "Mapa de calor completo ligantes x resíduos"),
            (self.pse_tab, "Triagem virtual racional"),
            (self.fp_tab, "Fingerprints"),
            (self.sim_tab, "Matriz de similaridade"),
            (self.cluster_tab, "Clusters"),
            (self.fp_analysis_tab, "Análises FP"),
            (self.fp_session_tab, "Sessão FP em PyMOL"),
        ]
        while self.inner.count():
            self.inner.removeTab(0)
        for widget, label in ordered_tabs:
            self.inner.addTab(widget, label)

    def load_all(self) -> None:
        super().load_all()
        wd = self._current_wd()
        if not wd:
            return
        self._load_cached_stats(wd)
        self._load_cached_residue_matrix(wd)
        self._load_cached_fp_analyses(wd)
        self._load_existing_fp_sessions(wd)
        self._refresh_pse_filter_list(wd)

    def compute_stats(self) -> None:
        wd = self._current_wd()
        if not wd:
            return

        result = load_analysis_summary(wd)
        if result is None:
            if not self.py_exe:
                QMessageBox.warning(self, "luna-env", "LUNA não detectado. Verifique a aba Setup.")
                return
            self.st_status.setText("Processando... (pode levar alguns minutos)")
            self.st_status.repaint()
            result = run_analysis(self.py_exe, str(wd))

        if "error" in result:
            self.st_status.setText("Erro")
            QMessageBox.critical(self, "Erro na analise", result["error"])
            return

        self._apply_stats_result(result)

    def compute_residue_matrix(self) -> None:
        wd = self._current_wd()
        if not wd:
            return

        result = load_residue_matrix_artifact(wd)
        if result is None or not self._residue_matrix_has_ligand_atoms(result):
            if not self.py_exe:
                QMessageBox.warning(self, "luna-env", "LUNA não detectado. Veja a aba Setup.")
                return
            self.hm_status.setText("Processando...")
            self.hm_status.repaint()
            result = run_residue_matrix(self.py_exe, str(wd), require_ligand_atoms=True)

        if "error" in result:
            self.hm_status.setText("Erro")
            self.hm_all_status.setText("Erro")
            QMessageBox.critical(self, "Erro na analise", result["error"])
            return

        self._apply_residue_result(result)

    def compute_fp_analyses(self) -> None:
        wd = self._current_wd()
        if not wd:
            return
        artifacts = load_fp_analysis_artifacts(wd)
        if not artifacts:
            self._fp_artifacts = {}
            self.fp_analysis_status.setText("Sem artefatos")
            self.fp_analysis_table.setRowCount(0)
            QMessageBox.information(
                self,
                "Análises FP",
                "Nenhuma analise de fingerprint foi encontrada neste workdir. "
                "Execute o LUNA novamente com fingerprints habilitados.",
            )
            return
        self._apply_fp_artifacts(artifacts)

    def _load_cached_stats(self, wd: Path) -> None:
        cached = load_analysis_summary(wd)
        if cached is not None:
            self._apply_stats_result(cached)

    def _load_cached_residue_matrix(self, wd: Path) -> None:
        cached = load_residue_matrix_artifact(wd)
        if cached is not None:
            self._apply_residue_result(cached)

    def _load_cached_fp_analyses(self, wd: Path) -> None:
        artifacts = load_fp_analysis_artifacts(wd)
        if artifacts:
            self._apply_fp_artifacts(artifacts)

    def _apply_stats_result(self, result: dict) -> None:
        self._last_analysis = result
        processed = len(result.get("entry_interaction_counts", {})) or result.get("entries", 0)
        self.st_status.setText(f"{processed} entradas processadas")
        self._populate_stats_scope(result)
        self._render_cached_stats_chart()

    def _load_fingerprints(self, wd: Path) -> None:
        file_path = self._selected_fingerprint_path(wd)
        if not file_path:
            self.fp_path_label.setText("Nenhum fingerprint calculado foi encontrado")
            self.fp_table.clear()
            self.fp_table.setRowCount(0)
            self.fp_table.setColumnCount(0)
            return

        self.fp_path_label.setText(f"{file_path} | index: ligand_id | colunas: bits")
        try:
            labels, bit_ids, matrix = load_ifp_sparse_matrix(file_path)
        except Exception as exc:
            QMessageBox.critical(self, "Erro ao ler fingerprints", str(exc))
            return

        limit = min(self.fp_rows.value(), len(labels))
        labels = labels[:limit]
        matrix = matrix[:limit, :]
        self.fp_table.clear()
        self.fp_table.setRowCount(len(labels))
        self.fp_table.setColumnCount(len(bit_ids))
        self.fp_table.setHorizontalHeaderLabels([str(bit_id) for bit_id in bit_ids])
        self.fp_table.setVerticalHeaderLabels(labels)

        self.fp_table.setSortingEnabled(False)
        for row_idx, row in enumerate(matrix):
            for col_idx, value in enumerate(row):
                numeric = float(value)
                text = str(int(numeric)) if numeric.is_integer() else f"{numeric:.4f}".rstrip("0").rstrip(".")
                self.fp_table.setItem(row_idx, col_idx, _sortable_item(text, numeric))
        self.fp_table.resizeRowsToContents()
        self.fp_table.setSortingEnabled(True)

    def _apply_residue_result(self, result: dict) -> None:
        self._residue_matrix = result
        interaction_types = result.get("interaction_types", [])
        summary = f"{len(result.get('entries', []))} entradas - {len(interaction_types)} tipos"
        hb_warning = self._hb_warning_message(list(interaction_types))
        status_text = summary + (f" | {hb_warning}" if hb_warning else "")
        self.hm_status.setText(status_text)
        self.hm_all_status.setText(status_text)

        self.cb_itype.blockSignals(True)
        current = self.cb_itype.currentText()
        self.cb_itype.clear()
        self.cb_itype.addItems(interaction_types)
        if current:
            idx = self.cb_itype.findText(current)
            if idx >= 0:
                self.cb_itype.setCurrentIndex(idx)
        self.cb_itype.blockSignals(False)

        if interaction_types:
            self._render_residue_heatmap()
            self._render_complete_heatmap()
        elif HAS_MPL:
            self.hm_fig.clear()
            ax = self.hm_fig.add_subplot(111)
            ax.text(0.5, 0.5, "Sem matriz de resíduos disponível", ha="center", va="center")
            self.hm_fig.tight_layout()
            self.hm_canvas.draw()
            self.hm_all_fig.clear()
            ax_all = self.hm_all_fig.add_subplot(111)
            ax_all.text(0.5, 0.5, "Sem matriz completa disponível", ha="center", va="center")
            self.hm_all_fig.tight_layout()
            self.hm_all_canvas.draw()

    def _apply_fp_artifacts(self, artifacts: dict[str, dict]) -> None:
        self._fp_artifacts = artifacts
        self._fp_dashboards = {}
        current = self.cb_fp_analysis_type.currentData() if self.cb_fp_analysis_type.count() else None
        self.cb_fp_analysis_type.blockSignals(True)
        self.cb_fp_analysis_type.clear()
        for ifp_type, artifact in sorted(artifacts.items()):
            label = artifact.get("ifp_label", ifp_type)
            self.cb_fp_analysis_type.addItem(f"{label} ({ifp_type})", ifp_type)
        idx = self.cb_fp_analysis_type.findData(current)
        self.cb_fp_analysis_type.setCurrentIndex(idx if idx >= 0 else 0)
        self.cb_fp_analysis_type.blockSignals(False)
        self._render_fp_analysis_table()
        self._sync_fp_session_types()

    def _populate_stats_scope(self, result: dict) -> None:
        if not hasattr(self, "cb_stats_scope"):
            return
        current = self.cb_stats_scope.currentData()
        self.cb_stats_scope.blockSignals(True)
        self.cb_stats_scope.clear()
        if self._is_trajectory_mode():
            self.cb_stats_scope.addItem("Todos os frames", "__all__")
            artifact = self._ensure_trajectory_matrix(allow_compute=False)
            entries = list((artifact or {}).get("entries", []) or [])
            if not entries:
                entries = sorted((result.get("entry_interaction_counts") or {}).keys())
            entries, _order = self._order_entries_by_activity(entries)
            for frame_name in entries:
                self.cb_stats_scope.addItem(f"Frame: {self._display_ligand_name(str(frame_name))}", frame_name)
        else:
            self.cb_stats_scope.addItem("Todos os ligantes", "__all__")
            artifact = self._ensure_trajectory_matrix(allow_compute=False)
            entries = list((artifact or {}).get("entries", []) or [])
            if not entries:
                entries = sorted((result.get("entry_interaction_counts") or {}).keys())
            entries, _order = self._order_entries_by_activity(entries)
            for ligand_name in entries:
                self.cb_stats_scope.addItem(f"Ligante: {self._display_ligand_name(str(ligand_name))}", ligand_name)
        idx = self.cb_stats_scope.findData(current)
        self.cb_stats_scope.setCurrentIndex(idx if idx >= 0 else 0)
        self.cb_stats_scope.blockSignals(False)

    def _render_cached_stats_chart(self) -> None:
        if HAS_MPL and self._last_analysis:
            self._render_stats_chart(self._last_analysis)

    def _known_stat_interactions(self) -> list[str]:
        names: set[str] = set()
        for key in (self._last_analysis.get("interaction_counts") or {}).keys():
            names.add(str(key))
        for counts in (self._last_analysis.get("entry_interaction_counts") or {}).values():
            for key in (counts or {}).keys():
                names.add(str(key))
        artifact = self._ensure_trajectory_matrix(allow_compute=False)
        for key in (artifact or {}).get("interaction_types", []) or []:
            names.add(str(key))
        return sorted(names, key=interaction_priority_key)

    def _toggle_stats_interaction_group(self, group: str) -> None:
        predicate = (
            is_unfavorable_or_repulsive_interaction
            if group == "unfavorable"
            else is_pi_stacking_interaction
        )
        labels = [name for name in self._known_stat_interactions() if predicate(name)]
        if not labels:
            QMessageBox.information(self, "Interações", "Nenhuma interação deste grupo foi encontrada no gráfico atual.")
            return
        if all(label in self._stats_hidden_interactions for label in labels):
            self._stats_hidden_interactions.difference_update(labels)
        else:
            self._stats_hidden_interactions.update(labels)
        self._render_cached_stats_chart()

    def _toggle_all_stats_interactions(self) -> None:
        labels = self._known_stat_interactions()
        if not labels:
            QMessageBox.information(self, "Interações", "Nenhuma interação foi encontrada no gráfico atual.")
            return
        if all(label in self._stats_hidden_interactions for label in labels):
            self._stats_hidden_interactions.difference_update(labels)
        else:
            self._stats_hidden_interactions.update(labels)
        self._render_cached_stats_chart()

    def _install_stats_legend_toggle(self, legend, labels: list[str]) -> None:
        if not HAS_MPL or legend is None or not hasattr(self, "st_canvas"):
            return
        if self._stats_legend_pick_cid is not None:
            try:
                self.st_canvas.mpl_disconnect(self._stats_legend_pick_cid)
            except Exception:
                pass
        self._stats_legend_artist_labels = {}
        handles = list(getattr(legend, "legend_handles", []) or [])
        if not handles and hasattr(legend, "legendHandles"):
            handles = list(getattr(legend, "legendHandles", []) or [])
        texts = list(legend.get_texts())
        for idx, label in enumerate(labels):
            artists = []
            if idx < len(handles):
                artists.append(handles[idx])
            if idx < len(texts):
                artists.append(texts[idx])
            for artist in artists:
                artist.set_picker(True)
                artist.set_alpha(0.28 if label in self._stats_hidden_interactions else 1.0)
                self._stats_legend_artist_labels[artist] = label
        self._stats_legend_pick_cid = self.st_canvas.mpl_connect(
            "pick_event",
            self._on_stats_legend_pick,
        )

    def _on_stats_legend_pick(self, event) -> None:
        label = self._stats_legend_artist_labels.get(event.artist)
        if not label:
            return
        if label in self._stats_hidden_interactions:
            self._stats_hidden_interactions.remove(label)
        else:
            self._stats_hidden_interactions.add(label)
        self._render_cached_stats_chart()

    def _install_fp_interaction_legend_toggle(self, legend, labels: list[str]) -> None:
        if not HAS_MPL or legend is None or not hasattr(self, "fp_interaction_canvas"):
            return
        if self._fp_interaction_legend_pick_cid is not None:
            try:
                self.fp_interaction_canvas.mpl_disconnect(self._fp_interaction_legend_pick_cid)
            except Exception:
                pass
        self._fp_interaction_legend_artist_labels = {}
        handles = list(getattr(legend, "legend_handles", []) or [])
        if not handles and hasattr(legend, "legendHandles"):
            handles = list(getattr(legend, "legendHandles", []) or [])
        texts = list(legend.get_texts())
        for idx, label in enumerate(labels):
            artists = []
            if idx < len(handles):
                artists.append(handles[idx])
            if idx < len(texts):
                artists.append(texts[idx])
            for artist in artists:
                artist.set_picker(True)
                artist.set_alpha(0.28 if label in self._fp_interaction_hidden_types else 1.0)
                self._fp_interaction_legend_artist_labels[artist] = label
        self._fp_interaction_legend_pick_cid = self.fp_interaction_canvas.mpl_connect(
            "pick_event",
            self._on_fp_interaction_legend_pick,
        )

    def _on_fp_interaction_legend_pick(self, event) -> None:
        label = self._fp_interaction_legend_artist_labels.get(event.artist)
        if not label:
            return
        if label in self._fp_interaction_hidden_types:
            self._fp_interaction_hidden_types.remove(label)
        else:
            self._fp_interaction_hidden_types.add(label)
        if hasattr(self, "cb_fp_analysis_type") and self.cb_fp_analysis_type.count():
            dashboard = self._ensure_fp_dashboard(self.cb_fp_analysis_type.currentData())
            if dashboard:
                self._render_fp_interaction_summary_plot(self._dashboard_with_selected_cutoff(dashboard))

    @staticmethod
    def _stats_bar_label_color(fill_color: str) -> str:
        color = str(fill_color or "").strip()
        if color.startswith("#"):
            hex_value = color[1:]
            if len(hex_value) == 3:
                hex_value = "".join(ch * 2 for ch in hex_value)
            if len(hex_value) >= 6:
                try:
                    r = int(hex_value[0:2], 16) / 255.0
                    g = int(hex_value[2:4], 16) / 255.0
                    b = int(hex_value[4:6], 16) / 255.0
                except ValueError:
                    return "#111111"
                luminance = (0.299 * r) + (0.587 * g) + (0.114 * b)
                return "#ffffff" if luminance < 0.35 else "#111111"
        return "#111111"

    def _render_stats_chart(self, result: dict) -> None:
        if self._is_trajectory_mode():
            self._render_trajectory_stats_chart(result)
            return

        self._render_ligand_stats_chart(result)

    def _render_ligand_stats_chart(self, result: dict) -> None:
        self.st_fig.clear()
        artifact = self._ensure_trajectory_matrix(allow_compute=True, require_ligand_atoms=False)
        if artifact is None:
            ax = self.st_fig.add_subplot(111)
            ax.text(
                0.5,
                0.5,
                "Sem matriz resíduo x interação.\nCalcule o mapa de calor por tipo ou reexecute a análise.",
                ha="center",
                va="center",
            )
            self.st_fig.tight_layout()
            self.st_canvas.draw()
            return

        scope = "__all__"
        if hasattr(self, "cb_stats_scope"):
            scope = self.cb_stats_scope.currentData() or "__all__"

        if scope == "__all__":
            residues, interaction_types, matrix = build_trajectory_frame_percentages(artifact)
            self._render_trajectory_stacked_bars(
                residues,
                interaction_types,
                matrix,
                title="Interações por aminoácido no conjunto de ligantes",
                ylabel="% de ligantes",
                percent_values=True,
            )
        else:
            residues, interaction_types, matrix = build_trajectory_entry_counts(artifact, str(scope))
            self._render_trajectory_stacked_bars(
                residues,
                interaction_types,
                matrix,
                title=f"Interações por aminoácido - ligante {self._display_ligand_name(str(scope))}",
                ylabel="Número de interações no ligante",
                percent_values=False,
            )

    def _is_trajectory_mode(self) -> bool:
        return self._results_trajectory_mode()

    @staticmethod
    def _residue_matrix_has_ligand_atoms(artifact: dict | None) -> bool:
        if not isinstance(artifact, dict):
            return False
        return bool(artifact.get("ligand_atoms")) and isinstance(artifact.get("ligand_atom_matrix"), dict)

    def _ligand_atom_map_path(self) -> Path | None:
        wd = self._current_wd()
        if wd is None:
            return None
        for relative in (
            "results/ligand_atom_map.png",
            "results/ligand_atom_map.jpg",
            "results/ligand_atom_map.jpeg",
        ):
            path = wd / relative
            if path.exists():
                return path
        return None

    def _ensure_trajectory_matrix(self, allow_compute: bool = True, require_ligand_atoms: bool = True) -> dict | None:
        if getattr(self, "_residue_matrix", None):
            if (
                not allow_compute
                or not require_ligand_atoms
                or self._residue_matrix_has_ligand_atoms(self._residue_matrix)
            ):
                return self._residue_matrix
        wd = self._current_wd()
        if not wd:
            return None
        cached = load_residue_matrix_artifact(wd)
        if cached is not None:
            if (
                allow_compute
                and require_ligand_atoms
                and self.py_exe
                and not self._residue_matrix_has_ligand_atoms(cached)
            ):
                result = run_residue_matrix(self.py_exe, str(wd), require_ligand_atoms=True)
                if "error" not in result:
                    self._residue_matrix = result
                    return result
            self._residue_matrix = cached
            return cached
        if not allow_compute or not self.py_exe:
            return None
        result = run_residue_matrix(self.py_exe, str(wd), require_ligand_atoms=require_ligand_atoms)
        if "error" in result:
            return None
        self._residue_matrix = result
        return result

    def _render_trajectory_stats_chart(self, result: dict) -> None:
        self.st_fig.clear()
        artifact = self._ensure_trajectory_matrix(allow_compute=True)
        if artifact is None:
            ax = self.st_fig.add_subplot(111)
            ax.text(
                0.5,
                0.5,
                "Sem matriz resíduo x interação para trajetória.\nCalcule o mapa de calor por tipo ou reexecute a análise.",
                ha="center",
                va="center",
            )
            self.st_fig.tight_layout()
            self.st_canvas.draw()
            return

        scope = "__all__"
        if hasattr(self, "cb_stats_scope"):
            scope = self.cb_stats_scope.currentData() or "__all__"

        if scope == "__all__":
            residues, interaction_types, matrix = build_trajectory_frame_percentages(artifact)
            atoms, atom_interaction_types, atom_matrix = build_ligand_atom_frame_percentages(artifact)
            self._render_trajectory_stacked_bars(
                residues,
                interaction_types,
                matrix,
                title="Interações por aminoácido ao longo da trajetória",
                ylabel="% de frames (entradas)",
                percent_values=True,
                ligand_atoms=atoms,
                ligand_interaction_types=atom_interaction_types,
                ligand_matrix=atom_matrix,
            )
        else:
            residues, interaction_types, matrix = build_trajectory_entry_counts(artifact, str(scope))
            atoms, atom_interaction_types, atom_matrix = build_ligand_atom_entry_counts(artifact, str(scope))
            self._render_trajectory_stacked_bars(
                residues,
                interaction_types,
                matrix,
                title=f"Interações por aminoácido - frame {self._display_ligand_name(str(scope))}",
                ylabel="Número de interações no frame",
                percent_values=False,
                ligand_atoms=atoms,
                ligand_interaction_types=atom_interaction_types,
                ligand_matrix=atom_matrix,
            )

    def _render_trajectory_dual_stacked_bars(
        self,
        residues: list[str],
        interaction_types: list[str],
        matrix: np.ndarray,
        title: str,
        ylabel: str,
        percent_values: bool,
        ligand_atoms: list[str],
        ligand_interaction_types: list[str],
        ligand_matrix: np.ndarray,
        ligand_title: str,
    ) -> None:
        residue_matrix = np.asarray(matrix, dtype=float) if matrix is not None else np.zeros((0, 0), dtype=float)
        atom_matrix = np.asarray(ligand_matrix, dtype=float)
        has_residue_axis = bool(residues and interaction_types and residue_matrix.size)
        has_ligand_axis = bool(ligand_atoms and ligand_interaction_types and atom_matrix.size)
        if not has_residue_axis and not has_ligand_axis:
            ax = self.st_fig.add_subplot(111)
            ax.text(0.5, 0.5, "Sem interações para esta visão.", ha="center", va="center")
            self.st_fig.tight_layout()
            self.st_canvas.draw()
            return

        display_residues = self._residue_xticklabels(residues) if has_residue_axis else []
        display_atoms = [str(atom) for atom in ligand_atoms] if has_ligand_axis else []
        atom_map_path = self._ligand_atom_map_path() if has_ligand_axis and self._is_trajectory_mode() else None
        atom_map_image = None
        if atom_map_path is not None:
            try:
                import matplotlib.image as mpimg

                atom_map_image = mpimg.imread(atom_map_path)
            except Exception:
                atom_map_image = None
        has_atom_map = atom_map_image is not None
        legend_labels: list[str] = []
        for label in list(interaction_types or []) + list(ligand_interaction_types or []):
            text = str(label)
            if text and text not in legend_labels:
                legend_labels.append(text)

        longest_interaction_label = max((len(str(label)) for label in legend_labels), default=1)
        longest_x_label = max((len(str(label)) for label in display_residues + display_atoms), default=1)
        legend_cols = max(1, min(3, len(legend_labels) or 1))
        legend_rows = max(1, math.ceil(max(1, len(legend_labels)) / legend_cols))
        legend_col_width = max(2.8, 0.105 * float(longest_interaction_label))
        total_x = len(display_residues) + len(display_atoms)
        matrices_for_scale = [
            np.asarray(candidate, dtype=float)
            for candidate in (residue_matrix if has_residue_axis else None, atom_matrix if has_ligand_axis else None)
            if candidate is not None and np.asarray(candidate).size
        ]
        positive_segments = [
            float(value)
            for candidate in matrices_for_scale
            for value in np.asarray(candidate, dtype=float).ravel()
            if float(value) > 0.0
        ]
        max_stack_value = max(
            [
                float(np.max(np.sum(candidate, axis=1)))
                for candidate in matrices_for_scale
                if candidate.ndim == 2 and candidate.size
            ]
            or [1.0]
        )
        min_label_segment = min(positive_segments) if positive_segments else 1.0
        readable_segment = max(0.5 if percent_values else 1.0, min_label_segment)
        target_label_px = 22.0 if percent_values else 18.0
        figure_dpi = float(getattr(self.st_fig, "dpi", 100.0) or 100.0)
        dynamic_label_height = (
            (max_stack_value * 1.18 / max(readable_segment, 1e-6)) * target_label_px / figure_dpi
        )
        legend_gap_in = 2.0 / 2.54
        legend_area_in = 0.42 * legend_rows + 0.38
        x_label_area_in = min(3.2, 1.35 + 0.018 * float(longest_x_label))
        width_in = max(
            17.0,
            3.8 + (0.48 * total_x),
            1.8 + (legend_cols * legend_col_width),
            6.5 + (0.04 * float(longest_x_label)),
            20.0 if has_atom_map else 0.0,
        )
        height_in = max(
            17.0,
            10.0 + (0.18 * max(len(display_residues), len(display_atoms))) + legend_area_in + legend_gap_in,
            dynamic_label_height + x_label_area_in + legend_area_in + legend_gap_in + 2.2,
        )
        self._resize_canvas(self.st_fig, self.st_canvas, width_in=width_in, height_in=height_in)
        self.st_fig.clear()

        def _draw_axis(
            ax,
            xlabels: list[str],
            type_labels: list[str],
            values_matrix: np.ndarray,
            panel_title: str,
            xlabel: str,
        ) -> None:
            if not xlabels or not type_labels or values_matrix.size == 0:
                ax.text(0.5, 0.5, "Sem dados.", ha="center", va="center")
                ax.set_axis_off()
                return
            x = np.arange(len(xlabels))
            bottoms = np.zeros(len(xlabels), dtype=float)
            for col_idx, interaction_type in enumerate(type_labels):
                if interaction_type in self._stats_hidden_interactions:
                    continue
                if col_idx >= values_matrix.shape[1]:
                    continue
                values = np.asarray(values_matrix[:, col_idx], dtype=float)
                if not np.any(values > 0.0):
                    continue
                bar_color = get_interaction_color(interaction_type)
                bars = ax.bar(
                    x,
                    values,
                    bottom=bottoms,
                    label=interaction_type,
                    color=bar_color,
                    edgecolor="#17324d",
                    linewidth=0.65,
                    width=0.82,
                )
                label_color = self._stats_bar_label_color(bar_color)
                label_fontsize = 7.0 if len(xlabels) > 32 else 7.8
                for bar, value, base in zip(bars, values, bottoms):
                    if value <= 0.0:
                        continue
                    if percent_values:
                        text = f"{value:.1f}%" if value < 10.0 else f"{int(round(value))}%"
                        min_visible = 0.5
                    else:
                        text = str(int(round(value)))
                        min_visible = 1.0
                    if value >= min_visible:
                        ax.text(
                            bar.get_x() + bar.get_width() / 2.0,
                            base + (value / 2.0),
                            text,
                            ha="center",
                            va="center",
                            fontsize=label_fontsize,
                            fontweight="bold",
                            color=label_color,
                            rotation=0,
                        )
                bottoms += values

            ax.set_xticks(x)
            x_label_fontsize = 9 if len(xlabels) <= 80 else 8
            ax.set_xticklabels(xlabels, rotation=90, fontsize=x_label_fontsize)
            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)
            ax.set_title(panel_title)
            ax.tick_params(axis="x", pad=14)
            ax.tick_params(axis="y", pad=14)
            ax.set_ylim(0, max(1.0, float(np.max(bottoms)) * 1.12))
            ax.grid(False)

        def _draw_atom_map_axis(ax) -> None:
            if atom_map_image is None:
                ax.set_axis_off()
                return
            ax.imshow(atom_map_image)
            ax.set_title("Estrutura 2D\nIDs dos átomos")
            ax.set_axis_off()

        if has_residue_axis and has_ligand_axis:
            if has_atom_map:
                axes = self.st_fig.subplots(
                    1,
                    3,
                    gridspec_kw={
                        "width_ratios": [
                            max(1, len(display_residues)),
                            max(1, len(display_atoms)),
                            max(8, int(0.55 * max(1, len(display_atoms)))),
                        ],
                        "wspace": 0.25,
                    },
                )
                ax_residue, ax_atom, ax_map = axes[0], axes[1], axes[2]
            else:
                axes = self.st_fig.subplots(
                    1,
                    2,
                    gridspec_kw={
                        "width_ratios": [max(1, len(display_residues)), max(1, len(display_atoms))],
                        "wspace": 0.22,
                    },
                )
                ax_residue, ax_atom = axes[0], axes[1]
                ax_map = None
            _draw_axis(ax_residue, display_residues, list(interaction_types), residue_matrix, title, "Aminoácidos")
            _draw_axis(ax_atom, display_atoms, ligand_interaction_types, atom_matrix, ligand_title, "Átomos do ligante")
            if ax_map is not None:
                _draw_atom_map_axis(ax_map)
        elif has_residue_axis:
            ax_residue = self.st_fig.add_subplot(111)
            _draw_axis(ax_residue, display_residues, list(interaction_types), residue_matrix, title, "Aminoácidos")
        else:
            if has_atom_map:
                ax_atom, ax_map = self.st_fig.subplots(
                    1,
                    2,
                    gridspec_kw={
                        "width_ratios": [max(1, len(display_atoms)), max(8, int(0.55 * max(1, len(display_atoms))))],
                        "wspace": 0.25,
                    },
                )
            else:
                ax_atom = self.st_fig.add_subplot(111)
                ax_map = None
            _draw_axis(ax_atom, display_atoms, ligand_interaction_types, atom_matrix, ligand_title, "Átomos do ligante")
            if ax_map is not None:
                _draw_atom_map_axis(ax_map)

        legend_handles = [
            Patch(
                facecolor=get_interaction_color(label),
                edgecolor="#17324d",
                label=label,
                alpha=0.28 if label in self._stats_hidden_interactions else 1.0,
            )
            for label in legend_labels
        ]
        legend = self.st_fig.legend(
            handles=legend_handles,
            loc="lower center",
            bbox_to_anchor=(0.5, 0.015),
            ncol=legend_cols,
            fontsize=8,
            frameon=False,
            title="Clique para ocultar/mostrar",
            columnspacing=1.4,
            handletextpad=0.7,
            borderaxespad=0.0,
        )
        self._install_stats_legend_toggle(legend, list(legend_labels))
        bottom_fraction = min(
            0.82,
            max(
                0.34,
                (x_label_area_in + legend_gap_in + legend_area_in) / max(height_in, 1.0),
            ),
        )
        self.st_fig.subplots_adjust(
            left=max(0.08, min(0.18, 0.065 + (0.004 * len(str(ylabel))))),
            right=0.985,
            top=0.90,
            bottom=bottom_fraction,
        )
        self.st_canvas.draw()

    def _render_trajectory_stacked_bars(
        self,
        residues: list[str],
        interaction_types: list[str],
        matrix: np.ndarray,
        title: str,
        ylabel: str,
        percent_values: bool,
        ligand_atoms: list[str] | None = None,
        ligand_interaction_types: list[str] | None = None,
        ligand_matrix: np.ndarray | None = None,
        ligand_title: str | None = None,
    ) -> None:
        self.st_fig.clear()
        if ligand_matrix is not None and ligand_atoms and ligand_interaction_types:
            self._render_trajectory_dual_stacked_bars(
                residues,
                interaction_types,
                matrix,
                title,
                ylabel,
                percent_values,
                list(ligand_atoms),
                list(ligand_interaction_types),
                np.asarray(ligand_matrix, dtype=float),
                ligand_title or "Interações por átomo do ligante",
            )
            return
        ax = self.st_fig.add_subplot(111)
        if matrix.size == 0 or not residues or not interaction_types:
            ax.text(0.5, 0.5, "Sem interações para esta visão.", ha="center", va="center")
            self.st_fig.tight_layout()
            self.st_canvas.draw()
            return

        display_residues = self._residue_xticklabels(residues)
        longest_interaction_label = max((len(str(label)) for label in interaction_types), default=1)
        longest_residue_label = max((len(str(label)) for label in display_residues), default=1)
        if longest_interaction_label >= 28:
            max_legend_cols = 2
        elif longest_interaction_label >= 18:
            max_legend_cols = 3
        else:
            max_legend_cols = 4
        legend_cols = max(1, min(max_legend_cols, len(interaction_types)))
        legend_rows = max(1, math.ceil(max(1, len(interaction_types)) / legend_cols))
        legend_col_width = max(2.8, 0.105 * float(longest_interaction_label))
        width_in = max(
            12.5,
            3.4 + (0.36 * len(display_residues)),
            1.8 + (legend_cols * legend_col_width),
            6.5 + (0.04 * float(longest_residue_label)),
        )
        height_in = max(15.0, 9.2 + (0.18 * len(display_residues)) + (0.78 * legend_rows))
        self._resize_canvas(self.st_fig, self.st_canvas, width_in=width_in, height_in=height_in)
        self.st_fig.clear()
        ax = self.st_fig.add_subplot(111)

        x = np.arange(len(display_residues))
        bottoms = np.zeros(len(display_residues), dtype=float)
        rendered_labels: list[str] = []
        for col_idx, interaction_type in enumerate(interaction_types):
            if interaction_type in self._stats_hidden_interactions:
                continue
            values = np.asarray(matrix[:, col_idx], dtype=float)
            if not np.any(values > 0.0):
                continue
            rendered_labels.append(interaction_type)
            bar_color = get_interaction_color(interaction_type)
            bars = ax.bar(
                x,
                values,
                bottom=bottoms,
                label=interaction_type,
                color=bar_color,
                edgecolor="#17324d",
                linewidth=0.65,
                width=0.82,
            )
            label_color = self._stats_bar_label_color(bar_color)
            label_fontsize = 6.0 if len(display_residues) > 32 else 6.5
            for bar, value, base in zip(bars, values, bottoms):
                if value <= 0.0:
                    continue
                if percent_values:
                    text = f"{int(round(value))}%"
                    min_visible = 6.0
                else:
                    text = str(int(round(value)))
                    min_visible = 1.2
                if value >= min_visible:
                    ax.text(
                        bar.get_x() + bar.get_width() / 2.0,
                        base + (value / 2.0),
                        text,
                        ha="center",
                        va="center",
                        fontsize=label_fontsize,
                        fontweight="bold",
                        color=label_color,
                        rotation=0,
                    )
            bottoms += values

        ax.set_xticks(x)
        x_label_fontsize = 9 if len(display_residues) <= 80 else 8
        ax.set_xticklabels(display_residues, rotation=90, fontsize=x_label_fontsize)
        ax.set_xlabel("Aminoácidos")
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.tick_params(axis="x", pad=14)
        ax.tick_params(axis="y", pad=14)
        ax.set_ylim(0, max(1.0, float(np.max(bottoms)) * 1.12))
        ax.grid(False)
        legend_handles = [
            Patch(
                facecolor=get_interaction_color(label),
                edgecolor="#17324d",
                label=label,
                alpha=0.28 if label in self._stats_hidden_interactions else 1.0,
            )
            for label in interaction_types
        ]
        legend = ax.legend(
            handles=legend_handles,
            loc="upper center",
            bbox_to_anchor=(0.5, -0.32),
            ncol=legend_cols,
            fontsize=8,
            frameon=False,
            title="Clique para ocultar/mostrar",
            columnspacing=1.4,
            handletextpad=0.7,
            borderaxespad=0.0,
        )
        self._install_stats_legend_toggle(legend, list(interaction_types))
        self.st_fig.subplots_adjust(
            left=max(0.08, min(0.18, 0.065 + (0.004 * len(str(ylabel))))),
            right=0.985,
            top=0.90,
            bottom=min(0.72, 0.28 + (0.055 * legend_rows)),
        )
        self.st_canvas.draw()

    def _render_residue_heatmap(self) -> None:
        if not HAS_MPL or not self._residue_matrix:
            return
        interaction_type = self.cb_itype.currentText()
        if not interaction_type:
            return

        residues = list(self._residue_matrix.get("residues", []) or [])
        data = self._residue_matrix.get("matrix", {}).get(interaction_type)
        if not data:
            return

        arr = np.asarray(data, dtype=float)
        ordered_entries, row_order = self._order_entries_by_activity(list(self._residue_matrix.get("entries", []) or []))
        if row_order:
            arr = arr[row_order, :]
        keep = arr.sum(axis=0) > 0
        arr = arr[:, keep]
        residues = [residue for residue, flag in zip(residues, keep) if flag]
        residue_labels = self._residue_xticklabels(residues)

        width_in = max(7.5 / 2.54, 2.2 + (0.5 / 2.54) * max(1, len(residue_labels)))
        title_width_in = 4.5 + (0.115 * len(str(interaction_type)))
        width_in = max(9.0, width_in, title_width_in)
        height_in = max(5.0, 2.5 + (0.16 * self._render_entry_count(arr.shape[0] if arr.ndim == 2 else 0)))
        self._resize_canvas(
            self.hm_fig,
            self.hm_canvas,
            width_in=width_in,
            height_in=height_in,
        )
        self.hm_fig.clear()
        ax = self.hm_fig.add_subplot(111)
        if arr.size == 0:
            ax.text(0.5, 0.5, f"Sem ocorrências de '{interaction_type}'", ha="center", va="center")
        else:
            im = ax.imshow(arr, cmap="viridis", aspect="auto")
            ax.set_xticks(list(range(len(residue_labels))))
            ax.set_xticklabels(residue_labels, rotation=90, fontsize=7)
            _apply_tick_labels(
                ax,
                [self._display_ligand_name(entry) for entry in ordered_entries],
                axis="y",
                max_labels=28,
                rotation=0,
                ligand_axis=True,
            )
            self._color_ticklabels_by_entry_group(ax, ordered_entries, axis="y")
            if len(ordered_entries) <= 220:
                ax.set_ylabel("Ligantes")
            ax.set_xlabel("Resíduos")
            title = "\n".join(textwrap.wrap(str(interaction_type), width=72, break_long_words=False)) or str(interaction_type)
            ax.set_title(title, pad=16)
            ax.tick_params(axis="y", pad=14)
            self.hm_fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
        self.hm_fig.tight_layout()
        self.hm_canvas.draw()

    def _toggle_complete_heatmap_group(self, group: str) -> None:
        if not self._residue_matrix:
            QMessageBox.information(self, "Mapa de calor", "Carregue ou calcule a matriz de resíduos primeiro.")
            return
        predicate = (
            is_unfavorable_or_repulsive_interaction
            if group == "unfavorable"
            else is_pi_stacking_interaction
        )
        labels = [
            str(name)
            for name in (self._residue_matrix.get("interaction_types", []) or [])
            if predicate(str(name))
        ]
        if not labels:
            QMessageBox.information(self, "Interações", "Nenhuma interação deste grupo foi encontrada no mapa atual.")
            return
        if all(label in self._complete_hidden_interactions for label in labels):
            self._complete_hidden_interactions.difference_update(labels)
        else:
            self._complete_hidden_interactions.update(labels)
        self._render_complete_heatmap()

    def _toggle_all_complete_heatmap_interactions(self) -> None:
        if not self._residue_matrix:
            QMessageBox.information(self, "Mapa de calor", "Carregue ou calcule a matriz de resíduos primeiro.")
            return
        labels = [str(name) for name in (self._residue_matrix.get("interaction_types", []) or [])]
        if not labels:
            QMessageBox.information(self, "Interações", "Nenhuma interação foi encontrada no mapa atual.")
            return
        if all(label in self._complete_hidden_interactions for label in labels):
            self._complete_hidden_interactions.difference_update(labels)
        else:
            self._complete_hidden_interactions.update(labels)
        self._render_complete_heatmap()

    def _render_complete_heatmap(self) -> None:
        if not HAS_MPL or not self._residue_matrix:
            return

        entries, residues, layered_cells, interaction_types = build_complete_heatmap_layers(self._residue_matrix)
        ordered_entries, row_order = self._order_entries_by_activity(entries)
        if row_order:
            layered_cells = [layered_cells[idx] for idx in row_order]
            entries = ordered_entries
        if self._complete_hidden_interactions:
            interaction_types = [
                name for name in interaction_types
                if name not in self._complete_hidden_interactions
            ]
            layered_cells = [
                [
                    [
                        interaction_name
                        for interaction_name in cell_types
                        if interaction_name not in self._complete_hidden_interactions
                    ]
                    for cell_types in row
                ]
                for row in layered_cells
            ]
        residue_labels = self._residue_xticklabels(residues)
        width_in = max(5.0 / 2.54, 2.2 + (0.5 / 2.54) * max(1, len(residue_labels)))
        legend_cols = max(1, min(3, len(interaction_types) or 1))
        legend_rows = max(1, math.ceil(max(1, len(interaction_types)) / legend_cols))
        height_in = max(5.8, 2.5 + (0.16 * self._render_entry_count(len(entries))) + (0.36 * legend_rows))
        self._resize_canvas(
            self.hm_all_fig,
            self.hm_all_canvas,
            width_in=width_in,
            height_in=height_in,
        )
        self.hm_all_fig.clear()
        ax = self.hm_all_fig.add_subplot(111)

        if not layered_cells or not residue_labels:
            ax.text(0.5, 0.5, "Sem dados de interação para o mapa de calor completo", ha="center", va="center")
            self.hm_all_fig.tight_layout()
            self.hm_all_canvas.draw()
            return

        ax.set_facecolor("#f4efe6")
        for row_idx, row in enumerate(layered_cells):
            for col_idx, cell_types in enumerate(row):
                if not cell_types:
                    ax.add_patch(
                        Rectangle(
                            (col_idx, row_idx),
                            1.0,
                            1.0,
                            facecolor="#f4efe6",
                            edgecolor="#ffffff",
                            linewidth=0.35,
                        )
                    )
                    continue
                ordered = sorted(cell_types, key=interaction_priority_key)
                stripe_width = 1.0 / max(1, len(ordered))
                for stripe_idx, interaction_name in enumerate(ordered):
                    ax.add_patch(
                        Rectangle(
                            (col_idx + stripe_idx * stripe_width, row_idx),
                            stripe_width,
                            1.0,
                            facecolor=get_interaction_color(interaction_name),
                            edgecolor="#ffffff",
                            linewidth=0.2,
                        )
                    )

        ax.set_xlim(0, len(residue_labels))
        ax.set_ylim(len(entries), 0)
        ax.set_xticks([idx + 0.5 for idx in range(len(residue_labels))])
        ax.set_xticklabels(residue_labels, rotation=90, fontsize=7)
        if len(entries) > 220:
            ax.set_yticks([])
            ax.set_ylabel("Todos os ligantes")
        else:
            ax.set_yticks([idx + 0.5 for idx in range(len(entries))])
            ax.set_yticklabels(
                [self._display_ligand_name(entry) for entry in entries],
                fontsize=7 * 0.85,
            )
            self._color_ticklabels_by_entry_group(ax, entries, axis="y", offset=0.5)
            ax.set_ylabel("Ligantes")
        ax.set_xlabel("Resíduos")
        ax.set_title("Tipos de interação por par ligante x resíduo")

        legend_types = sorted(interaction_types, key=interaction_priority_key)
        if legend_types:
            handles = [
                Patch(facecolor=get_interaction_color(name), edgecolor="none", label=name)
                for name in legend_types
            ]
            ax.legend(
                handles=handles,
                loc="upper center",
                bbox_to_anchor=(0.5, -0.09),
                ncol=legend_cols,
                fontsize=8,
                frameon=False,
            )
        displayed_entries = [self._display_ligand_name(entry) for entry in entries]
        longest_entry = max((len(label) for label in displayed_entries), default=1)
        left_margin = 0.08 if len(entries) > 220 else min(0.38, max(0.12, 0.06 + (0.006 * longest_entry)))
        self.hm_all_fig.subplots_adjust(
            left=left_margin,
            bottom=min(0.48, 0.12 + (0.045 * legend_rows)),
        )
        self.hm_all_canvas.draw()

    def _ensure_fp_dashboard(self, ifp_type: str | None) -> dict | None:
        if not ifp_type:
            return None
        labels_csv, labels_id_column, labels_column, task_kind_preference = self._results_label_settings()
        algorithm_preference = self._selected_fp_algorithm()
        use_otsu_threshold = bool(getattr(self.cfg, "fp_use_otsu_threshold", False))
        artifact = self._fp_artifacts.get(ifp_type)
        if artifact is None:
            return None

        wd = self._current_wd()
        if wd is None:
            return None

        seed_override = ""
        seed_file = str(getattr(self.cfg, "ifp_seed_file", "") or "").strip()
        if seed_file:
            try:
                seed_path = Path(seed_file)
                if seed_path.exists():
                    seed_override = seed_path.read_text(encoding="utf-8", errors="replace")
            except Exception:
                seed_override = ""
        random_seed = resolve_fp_random_seed(wd, artifact, seed_override)
        cache_key = (
            ifp_type,
            f"{labels_csv}|{labels_id_column}|{labels_column}|{task_kind_preference}|"
            f"{algorithm_preference}|otsu={int(use_otsu_threshold)}|seed={random_seed}",
        )
        if cache_key in self._fp_dashboards:
            cached = self._fp_dashboards[cache_key]
            cached_model = str(cached.get("model_name", "") or "")
            if cached_model != "Unavailable" and not cached_model.startswith("Fallback"):
                return cached

        detail_error = ""
        if load_fp_detail_artifact(wd, ifp_type) is None:
            detail_result = run_fp_detail_analysis(self.py_exe or "", str(wd), str(ifp_type))
            if isinstance(detail_result, dict) and detail_result.get("error"):
                detail_error = str(detail_result.get("error") or "")

        local_dashboard = build_fp_analysis_dashboard(
            wd,
            artifact,
            labels_csv=labels_csv,
            labels_id_column=labels_id_column,
            labels_column=labels_column,
            algorithm_preference=algorithm_preference,
            task_kind_preference=task_kind_preference,
            use_otsu_threshold=use_otsu_threshold,
            random_seed=random_seed,
        )
        dashboard = local_dashboard
        helper_dashboard = None
        helper_error = ""
        local_model = str(local_dashboard.get("model_name", "") or "")
        if self.py_exe and (local_model == "Unavailable" or local_model.startswith("Fallback")):
            helper_dashboard = run_fp_dashboard_analysis(
                self.py_exe,
                str(wd),
                str(ifp_type),
                labels_csv=labels_csv,
                labels_id_column=labels_id_column,
                labels_column=labels_column,
                algorithm_preference=algorithm_preference,
                task_kind_preference=task_kind_preference,
                use_otsu_threshold=use_otsu_threshold,
                random_seed=random_seed,
            )
            helper_model = str(helper_dashboard.get("model_name", "") or "")
            if "error" not in helper_dashboard and helper_model != "Unavailable" and not helper_model.startswith("Fallback"):
                dashboard = helper_dashboard
            else:
                helper_error = str(helper_dashboard.get("error", "") or "")
                if (
                    helper_dashboard
                    and "error" not in helper_dashboard
                    and str(helper_dashboard.get("model_name", "") or "") == "Unavailable"
                ):
                    dashboard = helper_dashboard
        if helper_error and str(dashboard.get("model_name", "") or "") == "Unavailable":
            dashboard["model_note"] = (
                str(dashboard.get("model_note", "") or "").strip()
                + (" " if str(dashboard.get("model_note", "") or "").strip() else "")
                + f" Falha no helper do luna-env: {helper_error}"
            ).strip()
        if detail_error and not dashboard.get("detail_available"):
            dashboard["detail_error"] = detail_error
        self._fp_dashboards[cache_key] = dashboard
        return dashboard

    def _selected_fp_pvalue_cutoff(self) -> float:
        if hasattr(self, "cb_fp_pvalue_cutoff") and self.cb_fp_pvalue_cutoff.count():
            try:
                return float(self.cb_fp_pvalue_cutoff.currentData() or 0.01)
            except Exception:
                return 0.01
        return 0.01

    def _selected_fp_algorithm(self) -> str:
        if hasattr(self, "cb_fp_algorithm") and self.cb_fp_algorithm.count():
            return str(self.cb_fp_algorithm.currentData() or "gradient_boosting")
        return "gradient_boosting"

    def _on_fp_analysis_controls_changed(self) -> None:
        self._render_fp_analysis_table()

    @staticmethod
    def _parse_fp_id_text(text: str) -> list[int]:
        ids: list[int] = []
        seen: set[int] = set()
        for token in str(text or "").replace("\n", ",").replace("\t", ",").replace(";", ",").split(","):
            for piece in token.split():
                cleaned = piece.strip()
                if not cleaned:
                    continue
                if "_" in cleaned:
                    cleaned = cleaned.split("_", 1)[0]
                try:
                    feature_id = int(float(cleaned))
                except Exception:
                    continue
                if feature_id not in seen:
                    seen.add(feature_id)
                    ids.append(feature_id)
        return ids

    def _current_fp_analysis_type(self) -> str:
        if hasattr(self, "cb_fp_analysis_type") and self.cb_fp_analysis_type.count():
            return str(self.cb_fp_analysis_type.currentData() or "")
        return ""

    def _apply_manual_fp_selection(self) -> None:
        ifp_type = self._current_fp_analysis_type()
        if not ifp_type:
            return
        ids = self._parse_fp_id_text(self.fp_manual_features_edit.text())
        if not ids:
            QMessageBox.information(
                self,
                "Fingerprints relevantes",
                "Informe pelo menos um ID de fingerprint para aplicar a seleção manual.",
            )
            return
        self._fp_manual_feature_ids[ifp_type] = ids
        self._render_fp_analysis_table()

    def _use_selected_fp_rows_as_manual_selection(self) -> None:
        if not hasattr(self, "fp_analysis_table"):
            return
        selected_rows = sorted({index.row() for index in self.fp_analysis_table.selectedIndexes()})
        ids: list[int] = []
        for row in selected_rows:
            item = self.fp_analysis_table.item(row, 0)
            if item is None:
                continue
            try:
                ids.append(int(float(item.text())))
            except Exception:
                continue
        if not ids:
            QMessageBox.information(self, "Fingerprints relevantes", "Selecione uma ou mais linhas na tabela.")
            return
        self.fp_manual_features_edit.setText(", ".join(str(feature_id) for feature_id in ids))
        self._apply_manual_fp_selection()

    def _reset_manual_fp_selection(self) -> None:
        ifp_type = self._current_fp_analysis_type()
        if ifp_type:
            self._fp_manual_feature_ids.pop(ifp_type, None)
        self.fp_manual_features_edit.clear()
        self._render_fp_analysis_table()

    def _dashboard_with_manual_selection(self, dashboard: dict, ifp_type: str | None) -> dict:
        ids = list(self._fp_manual_feature_ids.get(str(ifp_type or ""), []) or [])
        if not ids:
            return dashboard
        feature_by_id = {
            int(feature.get("feature_id", 0) or 0): feature
            for feature in list(dashboard.get("features", []) or [])
        }
        selected = [feature_by_id[feature_id] for feature_id in ids if feature_id in feature_by_id]
        class_counts: dict[str, int] = {}
        for feature in selected:
            class_name = str(feature.get("assigned_class", CLASS_UNRELIABLE))
            class_counts[class_name] = class_counts.get(class_name, 0) + 1
        class_share = {
            class_name: (100.0 * count / len(selected)) if selected else 0.0
            for class_name, count in class_counts.items()
        }
        filtered = dict(dashboard)
        filtered["important_features"] = selected
        filtered["class_counts"] = class_counts
        filtered["class_share"] = class_share
        filtered["important_selection"] = "manual"
        filtered["manual_feature_ids"] = ids
        filtered["manual_found_feature_ids"] = [int(feature.get("feature_id", 0) or 0) for feature in selected]
        return filtered

    def _dashboard_with_selected_cutoff(self, dashboard: dict) -> dict:
        cutoff = self._selected_fp_pvalue_cutoff()
        important = sorted(
            [
                feature
                for feature in list(dashboard.get("features", []) or [])
                if float(feature.get("importance_pvalue", 1.0) or 1.0) < cutoff
                or bool(feature.get("importance_selected", False))
            ],
            key=lambda row: (
                str(row.get("assigned_level", "") or ""),
                float(row.get("importance_pvalue", 1.0) or 1.0),
                -float(row.get("importance_score", 0.0) or 0.0),
                int(row.get("feature_id", 0) or 0),
            ),
        )
        class_counts: dict[str, int] = {}
        for feature in important:
            class_name = str(feature.get("assigned_class", CLASS_UNRELIABLE))
            class_counts[class_name] = class_counts.get(class_name, 0) + 1
        class_share = {
            class_name: (100.0 * count / len(important)) if important else 0.0
            for class_name, count in class_counts.items()
        }
        filtered = dict(dashboard)
        filtered["important_features"] = important
        filtered["class_counts"] = class_counts
        filtered["class_share"] = class_share
        filtered["important_selection"] = f"pvalue_lt_{cutoff:.2f}"
        filtered["pvalue_cutoff"] = cutoff
        return filtered

    def _render_fp_analysis_table(self) -> None:
        ifp_type = self.cb_fp_analysis_type.currentData() if self.cb_fp_analysis_type.count() else None
        if not ifp_type:
            self.fp_analysis_status.setText("Sem dados")
            self.fp_analysis_summary.setText("-")
            self.fp_analysis_formula.setText("-")
            self.fp_analysis_method.setText("-")
            self.fp_analysis_active_context.setText("-")
            self.fp_analysis_table.setRowCount(0)
            self._clear_fp_analysis_plots("Sem dados de fingerprints.")
            return

        self.fp_analysis_status.setText("Processando dashboard...")
        self.fp_analysis_status.repaint()
        try:
            dashboard = self._ensure_fp_dashboard(ifp_type)
        except Exception as exc:
            self.fp_analysis_status.setText("Erro")
            self.fp_analysis_summary.setText("-")
            self.fp_analysis_formula.setText("-")
            self.fp_analysis_method.setText("-")
            self.fp_analysis_active_context.setText("-")
            self.fp_analysis_table.setRowCount(0)
            self._clear_fp_analysis_plots("Erro ao calcular o dashboard.")
            QMessageBox.critical(self, "Erro na analise de FP", str(exc))
            return

        if not dashboard:
            self.fp_analysis_status.setText("Sem dados")
            self.fp_analysis_summary.setText("-")
            self.fp_analysis_formula.setText("-")
            self.fp_analysis_method.setText("-")
            self.fp_analysis_active_context.setText("-")
            self.fp_analysis_table.setRowCount(0)
            self._clear_fp_analysis_plots("Sem dados de fingerprints.")
            return

        manual_ids = self._fp_manual_feature_ids.get(str(ifp_type), [])
        if manual_ids:
            self.fp_manual_features_edit.setText(", ".join(str(feature_id) for feature_id in manual_ids))
        elif not self.fp_manual_features_edit.hasFocus():
            self.fp_manual_features_edit.clear()

        filtered_dashboard = self._dashboard_with_manual_selection(
            self._dashboard_with_selected_cutoff(dashboard),
            str(ifp_type),
        )
        features = list(filtered_dashboard.get("features", []) or [])
        important = list(filtered_dashboard.get("important_features", []) or [])
        pvalue_cutoff = float(filtered_dashboard.get("pvalue_cutoff", 0.01) or 0.01)
        threshold_pct = float(dashboard.get("threshold_pct", 0.0) or 0.0)
        model_name = str(dashboard.get("model_name", "Unavailable") or "Unavailable")
        model_note = str(dashboard.get("model_note", "") or "")
        label_source = str(dashboard.get("label_source", "derived_clusters") or "derived_clusters")
        label_kind = str(dashboard.get("label_kind", "classification") or "classification")
        labels_id_column = str(dashboard.get("labels_id_column", "") or "").strip()
        labels_column = str(dashboard.get("labels_column", "") or "").strip()
        matched_molecules = int(dashboard.get("matched_molecules", 0) or 0)
        reliable_count = sum(1 for row in features if row.get("reliable"))
        algorithm_preference = str(dashboard.get("algorithm_preference", "gradient_boosting") or "gradient_boosting")
        task_kind_preference = str(dashboard.get("task_kind_preference", label_kind) or label_kind)
        interaction_threshold_pct = float(dashboard.get("interaction_threshold_pct", 100.0) or 100.0)
        interaction_threshold_source = str(dashboard.get("interaction_threshold_source", "") or "")
        pair_threshold_pct = float(dashboard.get("pair_threshold_pct", interaction_threshold_pct) or interaction_threshold_pct)
        pair_threshold_source = str(dashboard.get("pair_threshold_source", interaction_threshold_source) or interaction_threshold_source)
        threshold_source = str(dashboard.get("threshold_source", "") or "")
        use_otsu_threshold = bool(dashboard.get("use_otsu_threshold", False))
        class_zscore_mean = float(dashboard.get("class_zscore_mean", 0.0) or 0.0)
        class_zscore_std = float(dashboard.get("class_zscore_std", 0.0) or 0.0)
        importance_zscore_mean = float(dashboard.get("importance_zscore_mean", 0.0) or 0.0)
        importance_zscore_std = float(dashboard.get("importance_zscore_std", 0.0) or 0.0)
        importance_eligible_count = int(
            dashboard.get(
                "importance_eligible_count",
                sum(1 for row in features if row.get("importance_eligible")),
            )
            or 0
        )
        random_seed = int(dashboard.get("random_seed", 0) or 0)
        detail_note = ""
        if not dashboard.get("detail_available"):
            detail_note = str(dashboard.get("detail_error", "") or "").strip()
        selection_mode = str(filtered_dashboard.get("important_selection", "") or "")
        if selection_mode == "manual":
            selection_text = (
                f"seleção manual ({len(important)} de {len(filtered_dashboard.get('manual_feature_ids', []) or [])} IDs encontrados)"
            )
        elif selection_mode == "per_level_pvalue_or_otsu":
            selection_text = f"selecionadas por modelo por nível (p < {pvalue_cutoff:.2f} ou Otsu)"
        else:
            selection_text = f"selecionadas por p < {pvalue_cutoff:.2f}"
        level_assignment = dashboard.get("level_assignment") or {}
        assigned_matrix = dashboard.get("assigned_matrix") or {}
        level_threshold_pct = float(level_assignment.get("threshold_pct", 100.0) or 100.0)
        level_threshold_source = str(level_assignment.get("threshold_source", "") or "")
        assigned_count = int(level_assignment.get("assigned_count", 0) or 0)
        undetermined_count = int(level_assignment.get("undetermined_count", 0) or 0)
        assigned_matrix_note = (
            " | matriz FP regravada com níveis"
            if bool(assigned_matrix.get("rewritten", False))
            else ""
        )
        self.fp_analysis_status.setText(
            f"{dashboard.get('total_molecules', 0)} moléculas - {len(features)} features"
        )
        self.fp_analysis_summary.setText(
            f"Limiar de atribuição da classe: {threshold_pct:.2f}% | "
            f"fonte: {threshold_source or 'z-score'}"
            + (" (Otsu habilitado)" if use_otsu_threshold else "")
            + " | "
            f"Limiar de atribuição do nível: {level_threshold_pct:.2f}% | "
            f"fonte nível: {level_threshold_source or 'z-score'} | "
            f"Limiar de atribuição do tipo de interação: {interaction_threshold_pct:.2f}% | "
            f"Limiar do par interação/resíduo: {pair_threshold_pct:.2f}% | "
            f"Features confiáveis por classe: {reliable_count}/{len(features)} | "
            f"Features confiáveis por nível: {assigned_count}/{len(features)}"
            + (f" ({undetermined_count} indeterminados)" if undetermined_count else "")
            + assigned_matrix_note
            + " | "
            f"Features elegíveis para importância: {importance_eligible_count}/{len(features)} | "
            f"Features importantes: {len(important)} | "
            f"Rótulos: {'CSV externo' if label_source == 'external_csv' else 'fallback automático'}"
            + (f" [ID: {labels_id_column}]" if label_source == 'external_csv' and labels_id_column else "")
            + (f" ({labels_column})" if label_source == 'external_csv' and labels_column else "")
            + (f" | pareadas: {matched_molecules}" if label_source == "external_csv" else "")
            + f" | tarefa: {'regressão' if label_kind == 'regression' else 'classificação'}"
            + f" | seed: {random_seed}"
            + (f" | detalhe FP: indisponível ({detail_note})" if detail_note else "")
            + f" | {selection_text} | "
            + f"Algoritmo: {algorithm_preference} | "
            + f"Modelo: {model_name}. {model_note}"
        )
        self.fp_analysis_formula.setText(
            "Equação de Keiser and Hert [1] para transformar os s-score dos coeficientes de importâncias em p-values: "
            "p = 1 - exp(-exp(((-z*pi)/sqrt(6)) - 0.577215665)), "
            "onde z é o Z-score Importance mostrado na tabela para a feature. "
            "A coluna 'Cobertura (%)' segue o bit na base inteira, mas os percentuais do perfil "
            "da base usam apenas as ocorrências classificadas do bit. "
            f"Os gráficos abaixo usam {selection_text}; os modelos estocásticos usam seed {random_seed}."
        )
        interaction_rule = (
            f"limiar do par interação/resíduo = {pair_threshold_pct:.2f}% "
            "definido pelo menor percentual entre pares com z-score > 1"
            if pair_threshold_source == "zscore_gt_1"
            else (
                f"limiar do par interação/resíduo = {pair_threshold_pct:.2f}% definido por Otsu's Thresholding"
                if pair_threshold_source in {"otsu", "otsu_single_value"}
                else "sem pares interação/resíduo com z-score > 1; apenas 100% foi aceito como par prevalente"
            )
        )
        active_ifp_label = str(filtered_dashboard.get("ifp_label", ifp_type) or ifp_type)
        self.fp_analysis_active_context.setText(
            f"Gráficos ativos: base = {active_ifp_label}; "
            f"{selection_text}; "
            f"algoritmo = {algorithm_preference}."
        )
        self.fp_analysis_method.setText(
            f"Configuração atual: cutoff p-value = {pvalue_cutoff:.2f}; "
            f"tarefa = {task_kind_preference}; algoritmo solicitado = {algorithm_preference}; modelo usado = {model_name}; seed = {random_seed}. "
            f"Z-score de classe: média = {class_zscore_mean:.4f}, desvio = {class_zscore_std:.4f}. "
            f"Z-score Importance: média = {importance_zscore_mean:.4f}, desvio = {importance_zscore_std:.4f}. "
            f"Os modelos de importância são ajustados separadamente por nível assinado; fingerprints sem nível "
            f"assinado ficam fora do treino e da importância. "
            f"Interação prevalente: para cada feature importante, calcula-se a frequência percentual "
            f"do par tipo de interação/resíduo nos shells reais do LUNA. O resíduo só aparece no gráfico "
            f"quando o par exato passa o {interaction_rule}."
            + (f" Sem fp_detail, estes dois gráficos não recebem contagens de interação/resíduo: {detail_note}" if detail_note else "")
        )
        self.fp_analysis_table.setSortingEnabled(False)
        self.fp_analysis_table.setRowCount(len(features))

        for row_index, feature in enumerate(features):
            feature_id = int(feature.get("feature_id", 0))
            breakdown = feature.get("class_breakdown", {}) or {}
            class_percentages = feature.get("class_percentages", {}) or {}
            profile_parts = [
                f"{label}: {count} ({float(class_percentages.get(label, 0.0)):.2f}%)"
                for label, count in sorted(breakdown.items(), key=lambda item: (-item[1], item[0]))
            ]
            missing_molecules = int(feature.get("missing_molecules", 0) or 0)
            if missing_molecules > 0:
                profile_parts.append(f"Sem ocorrencia na base: {missing_molecules}")
            profile = "; ".join(profile_parts)
            shell_levels = self._format_shell_levels(
                feature.get("shell_levels"),
                feature.get("shell_level_breakdown"),
            )
            collision_levels = self._format_shell_levels(
                feature.get("collision_shell_levels"),
                feature.get("collision_level_breakdown"),
            )
            assigned_level = str(feature.get("assigned_level", "") or "").strip()
            assigned_level_display = assigned_level or str(feature.get("assigned_level_label", "") or "").strip()
            items = [
                _sortable_item(str(feature_id), feature_id),
                _sortable_item(str(feature.get("molecule_hits", 0)), int(feature.get("molecule_hits", 0) or 0)),
                _sortable_item(f"{float(feature.get('coverage_pct', 0.0)):.2f}", float(feature.get("coverage_pct", 0.0) or 0.0)),
                _sortable_item(f"{float(feature.get('top_class_pct', 0.0)):.2f}", float(feature.get("top_class_pct", 0.0) or 0.0)),
                _sortable_item(f"{float(feature.get('zscore', 0.0)):.3f}", float(feature.get("zscore", 0.0) or 0.0)),
                _sortable_item(str(feature.get("assigned_class", "-")), str(feature.get("assigned_class", "-"))),
                _sortable_item(f"{float(feature.get('importance_score', 0.0)):.6f}", float(feature.get("importance_score", 0.0) or 0.0)),
                _sortable_item(f"{float(feature.get('importance_zscore', 0.0)):.3f}", float(feature.get("importance_zscore", 0.0) or 0.0)),
                _sortable_item(f"{float(feature.get('importance_pvalue', 1.0)):.6f}", float(feature.get("importance_pvalue", 1.0) or 1.0)),
                _sortable_item(str(feature.get("collision_hits", 0)), int(feature.get("collision_hits", 0) or 0)),
                _sortable_item(assigned_level_display or "-", assigned_level_display),
                _sortable_item(shell_levels or "-", shell_levels),
                _sortable_item(collision_levels or "-", collision_levels),
                _sortable_item(profile, profile),
            ]
            for col_index, item in enumerate(items):
                self.fp_analysis_table.setItem(row_index, col_index, item)
        self.fp_analysis_table.resizeRowsToContents()
        self.fp_analysis_table.setSortingEnabled(True)
        try:
            self._render_fp_analysis_plots(filtered_dashboard)
        except Exception as exc:
            self.fp_analysis_status.setText("Tabela carregada; erro ao desenhar gráficos.")
            self._clear_fp_analysis_plots(
                f"Tabela carregada, mas os gráficos de FP não puderam ser desenhados.\n{type(exc).__name__}: {exc}"
            )
        self._sync_fp_session_from_table()

    def _clear_fp_analysis_plots(self, message: str) -> None:
        if not HAS_MPL or self.fp_plot_tabs is None:
            return
        self._ensure_fp_plot_canvases()
        for fig, canvas in (
            (self.fp_class_fig, self.fp_class_canvas),
            (self.fp_assign_fig, self.fp_assign_canvas),
            (self.fp_cover_fig, self.fp_cover_canvas),
            (self.fp_heatmap_fig, self.fp_heatmap_canvas),
            (self.fp_interaction_assign_fig, self.fp_interaction_assign_canvas),
            (self.fp_interaction_fig, self.fp_interaction_canvas),
            (self.fp_interaction_heatmap_fig, self.fp_interaction_heatmap_canvas),
        ):
            fig.clear()
            ax = fig.add_subplot(111)
            ax.text(0.5, 0.5, message, ha="center", va="center")
            fig.tight_layout()
            canvas.draw()

    def _render_fp_analysis_plots(self, dashboard: dict) -> None:
        if not HAS_MPL or self.fp_plot_tabs is None:
            return
        self._ensure_fp_plot_canvases()
        self._render_fp_class_summary_plot(dashboard)
        self._render_fp_assignment_plot(dashboard)
        self._render_fp_coverage_plot(dashboard)
        self._render_fp_importance_heatmap(dashboard)
        self._render_fp_interaction_assignment_plot(dashboard)
        self._render_fp_interaction_summary_plot(dashboard)
        self._render_fp_interaction_heatmap(dashboard)

    def _render_fp_class_summary_plot(self, dashboard: dict) -> None:
        labels = [label for label in FP_CLASS_ORDER if (dashboard.get("class_share", {}) or {}).get(label, 0.0) > 0]
        self._resize_canvas(
            self.fp_class_fig,
            self.fp_class_canvas,
            width_in=max(10.2, 1.15 * max(1, len(labels)) + 4.2),
            height_in=6.2,
        )
        self.fp_class_fig.clear()
        ax = self.fp_class_fig.add_subplot(111)

        class_share = dashboard.get("class_share", {}) or {}
        class_counts = dashboard.get("class_counts", {}) or {}
        labels = [label for label in FP_CLASS_ORDER if class_share.get(label, 0.0) > 0]
        if not labels:
            ax.text(0.5, 0.5, "Sem features importantes para resumir.", ha="center", va="center")
            self.fp_class_fig.tight_layout()
            self.fp_class_canvas.draw()
            return

        values = [float(class_share.get(label, 0.0)) for label in labels]
        colors = [dashboard.get("class_colors", {}).get(label, "#6f9ec7") for label in labels]
        bars = ax.bar(range(len(labels)), values, color=colors)
        ax.set_ylabel("% de features importantes")
        ax.set_xticks(range(len(labels)))
        ax.set_xticklabels(labels, rotation=35, ha="right")
        ax.set_ylim(0, max(values) * 1.25 if values else 1.0)
        ax.set_title("Distribuição das classes entre as features mais importantes")
        for bar, label in zip(bars, labels):
            count = int(class_counts.get(label, 0))
            ax.text(bar.get_x() + bar.get_width() / 2.0, bar.get_height() + 0.8, str(count), ha="center", va="bottom")
        self.fp_class_fig.tight_layout()
        self.fp_class_canvas.draw()

    def _render_fp_assignment_plot(self, dashboard: dict) -> None:
        features = list(dashboard.get("important_features", []) or [])
        self._resize_canvas(
            self.fp_assign_fig,
            self.fp_assign_canvas,
            width_in=12.0,
            height_in=max(7.0, 0.35 * max(1, len(features)) + 3.2),
        )
        self.fp_assign_fig.clear()
        ax = self.fp_assign_fig.add_subplot(111)

        if not features:
            ax.text(0.5, 0.5, "Sem features importantes para plotar.", ha="center", va="center")
            self.fp_assign_fig.tight_layout()
            self.fp_assign_canvas.draw()
            return

        y_pos = list(range(len(features)))
        running = [0.0] * len(features)
        colors = dashboard.get("class_colors", {})
        for class_name in FP_CLASS_ORDER:
            widths = [float((feature.get("class_percentages", {}) or {}).get(class_name, 0.0)) for feature in features]
            if not any(width > 0.0 for width in widths):
                continue
            left_values = list(running)
            ax.barh(
                y_pos,
                widths,
                left=left_values,
                color=colors.get(class_name, "#cccccc"),
                edgecolor="white",
                label=class_name,
            )
            for row_idx, (left, width) in enumerate(zip(left_values, widths)):
                if width >= 4.0:
                    ax.text(
                        left + (width / 2.0),
                        row_idx,
                        f"{width:.1f}%",
                        ha="center",
                        va="center",
                        fontsize=7,
                        color=self._contrast_text_color(colors.get(class_name, "#cccccc")),
                    )
            running = [left + width for left, width in zip(running, widths)]

        threshold_pct = float(dashboard.get("threshold_pct", 0.0) or 0.0)
        level_assignment = dashboard.get("level_assignment") or {}
        level_threshold_pct = float(level_assignment.get("threshold_pct", 100.0) or 100.0)
        ax.axvline(threshold_pct, color="#444444", linestyle="--", linewidth=1)
        ax.text(
            min(103.0, max(1.0, threshold_pct + 1.0)),
            -0.72,
            f"classe {threshold_pct:.1f}%",
            fontsize=8,
            va="center",
            ha="left",
            color="#333333",
        )
        unreliable_rows = [
            idx for idx, feature in enumerate(features)
            if str(feature.get("assigned_class")) == CLASS_UNRELIABLE
        ]
        if unreliable_rows:
            ax.scatter(
                [102.0] * len(unreliable_rows),
                unreliable_rows,
                marker="x",
                color="black",
                s=36,
                zorder=5,
                clip_on=False,
            )

        ax.set_xlim(0, 105)
        ax.set_xlabel("Frequência de atribuição de cada classe (%)")
        ax.set_ylabel("ID da feature")
        ax.set_yticks(y_pos)
        ax.set_yticklabels([self._fp_feature_label(feature) for feature in features])
        ax.invert_yaxis()
        ax.set_title(
            "Frequência de atribuição de classes nas features importantes\n"
            f"limiar classe = {threshold_pct:.2f}% | limiar nível = {level_threshold_pct:.2f}%"
        )
        handles, labels_text = ax.get_legend_handles_labels()
        if handles:
            legend_cols = min(3, max(1, len(handles)))
            ax.legend(
                handles,
                labels_text,
                loc="upper center",
                bbox_to_anchor=(0.5, -0.18),
                ncol=legend_cols,
                fontsize=8,
                frameon=False,
            )
            self.fp_assign_fig.subplots_adjust(bottom=min(0.42, 0.16 + (0.06 * math.ceil(len(handles) / legend_cols))))
        else:
            self.fp_assign_fig.tight_layout()
        self.fp_assign_canvas.draw()

    def _render_fp_coverage_plot(self, dashboard: dict) -> None:
        features = list(dashboard.get("important_features", []) or [])
        self._resize_canvas(
            self.fp_cover_fig,
            self.fp_cover_canvas,
            width_in=12.0,
            height_in=max(7.0, 0.35 * max(1, len(features)) + 3.2),
        )
        self.fp_cover_fig.clear()
        ax = self.fp_cover_fig.add_subplot(111)

        if not features:
            ax.text(0.5, 0.5, "Sem features importantes para plotar.", ha="center", va="center")
            self.fp_cover_fig.tight_layout()
            self.fp_cover_canvas.draw()
            return

        colors = dashboard.get("class_colors", {})
        y_pos = list(range(len(features)))
        coverage = [float(feature.get("coverage_pct", 0.0)) for feature in features]
        bar_colors = [colors.get(str(feature.get("assigned_class")), "#6f9ec7") for feature in features]
        bars = ax.barh(y_pos, coverage, color=bar_colors)
        ax.set_xlim(0, 110)
        ax.set_xlabel("% de fingerprints contendo a feature")
        ax.set_ylabel("ID da feature")
        ax.set_yticks(y_pos)
        ax.set_yticklabels([self._fp_feature_label(feature) for feature in features])
        ax.invert_yaxis()
        ax.set_title("Cobertura das features importantes e importância do modelo")

        max_importance = max(float(feature.get("importance_pct", 0.0)) for feature in features) if features else 0.0
        star_positions = [
            float(feature.get("importance_pct", 0.0)) if max_importance > 0 else float(feature.get("coverage_pct", 0.0))
            for feature in features
        ]
        ax.scatter(star_positions, y_pos, marker="*", color="red", s=44, zorder=5)

        for idx, (bar, feature, star_x) in enumerate(zip(bars, features, star_positions)):
            coverage_value = float(feature.get("coverage_pct", 0.0))
            importance_value = float(feature.get("importance_score", 0.0))
            ax.text(max(1.0, coverage_value * 0.03), idx, f"{coverage_value:.1f}", va="center", ha="left", fontsize=8)
            if importance_value > 0:
                ax.text(min(108.0, star_x + 1.5), idx, f"{importance_value:.4f}", va="center", ha="left", fontsize=8, color="red")
        self.fp_cover_fig.tight_layout()
        self.fp_cover_canvas.draw()

    def _render_fp_importance_heatmap(self, dashboard: dict) -> None:
        features = list(dashboard.get("important_features", []) or [])
        entries = list(dashboard.get("entry_labels", []) or [])
        entries, _row_order = self._order_entries_by_activity(entries)
        self._resize_canvas(
            self.fp_heatmap_fig,
            self.fp_heatmap_canvas,
            width_in=max(9.5, 0.34 * max(1, len(features)) + 3.5),
            height_in=max(7.5, 0.12 * self._render_entry_count(len(entries)) + 3.0),
        )
        self.fp_heatmap_fig.clear()
        ax = self.fp_heatmap_fig.add_subplot(111)

        if not features or not entries:
            ax.text(0.5, 0.5, "Sem features importantes para gerar o mapa de calor.", ha="center", va="center")
            self.fp_heatmap_fig.tight_layout()
            self.fp_heatmap_canvas.draw()
            return

        class_order = [label for label in FP_CLASS_ORDER if label != CLASS_UNRELIABLE] + [CLASS_UNRELIABLE]
        class_to_id = {label: idx for idx, label in enumerate(class_order, start=1)}
        matrix = np.zeros((len(entries), len(features)), dtype=int)
        for col_idx, feature in enumerate(features):
            class_id = class_to_id.get(str(feature.get("assigned_class", CLASS_UNRELIABLE)), class_to_id[CLASS_UNRELIABLE])
            entry_counts = feature.get("entry_counts", {}) or {}
            for row_idx, entry_name in enumerate(entries):
                if float(entry_counts.get(entry_name, 0.0) or 0.0) > 0.0:
                    matrix[row_idx, col_idx] = class_id

        colors = ["#3f3f3f"] + [dashboard.get("class_colors", {}).get(label, "#7d7d7d") for label in class_order]
        cmap = ListedColormap(colors)
        norm = BoundaryNorm(range(len(colors) + 1), cmap.N)
        im = ax.imshow(matrix, cmap=cmap, norm=norm, aspect="auto", interpolation="nearest")
        ax.set_xticks(list(range(len(features))))
        ax.set_xticklabels([self._fp_feature_label(feature) for feature in features], rotation=90, fontsize=7)
        _apply_tick_labels(
            ax,
            [self._display_ligand_name(entry) for entry in entries],
            axis="y",
            max_labels=28,
            rotation=0,
            ligand_axis=True,
        )
        self._color_ticklabels_by_entry_group(ax, entries, axis="y")
        if len(entries) <= 220:
            ax.set_ylabel("Ligantes")
        ax.set_xlabel("ID da feature")
        ax.set_title("Mapa de presença das features importantes por classe")

        star_rows = [idx for idx, feature in enumerate(features) if float(feature.get("importance_pvalue", 1.0) or 1.0) < 0.01]
        if star_rows:
            ax.scatter(star_rows, [-1.1] * len(star_rows), marker="*", color="red", s=28, clip_on=False)

        cbar = self.fp_heatmap_fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
        cbar.set_ticks([idx + 0.5 for idx in range(len(class_order) + 1)])
        cbar.set_ticklabels(["Ausente"] + class_order)
        self.fp_heatmap_fig.tight_layout()
        self.fp_heatmap_canvas.draw()

    def _render_fp_interaction_assignment_plot(self, dashboard: dict) -> None:
        features = list(dashboard.get("important_features", []) or [])
        self._resize_canvas(
            self.fp_interaction_assign_fig,
            self.fp_interaction_assign_canvas,
            width_in=12.0,
            height_in=max(7.0, 0.35 * max(1, len(features)) + 3.2),
        )
        self.fp_interaction_assign_fig.clear()
        ax = self.fp_interaction_assign_fig.add_subplot(111)

        if not features:
            ax.text(0.5, 0.5, "Sem features importantes para plotar.", ha="center", va="center")
            self.fp_interaction_assign_fig.tight_layout()
            self.fp_interaction_assign_canvas.draw()
            return

        interaction_names: list[str] = []
        for feature in features:
            for interaction_name, count in sorted((feature.get("interaction_breakdown") or {}).items()):
                if int(count) > 0 and interaction_name not in interaction_names:
                    interaction_names.append(interaction_name)
        interaction_names = sorted(interaction_names, key=interaction_priority_key)

        y_pos = list(range(len(features)))
        running = [0.0] * len(features)
        for interaction_name in interaction_names:
            widths = []
            for feature in features:
                breakdown = feature.get("interaction_breakdown") or {}
                total = float(sum(breakdown.values()) or 0.0)
                widths.append((100.0 * float(breakdown.get(interaction_name, 0) or 0.0) / total) if total > 0 else 0.0)
            if not any(width > 0.0 for width in widths):
                continue
            color = get_interaction_color(interaction_name)
            left_values = list(running)
            ax.barh(
                y_pos,
                widths,
                left=left_values,
                color=color,
                edgecolor="white",
                label=interaction_name,
            )
            for row_idx, (left, width) in enumerate(zip(left_values, widths)):
                if width >= 4.0:
                    ax.text(
                        left + (width / 2.0),
                        row_idx,
                        f"{width:.1f}%",
                        ha="center",
                        va="center",
                        fontsize=7,
                        color=self._contrast_text_color(color),
                    )
            running = [left + width for left, width in zip(running, widths)]

        threshold_pct = float(dashboard.get("interaction_threshold_pct", 100.0) or 100.0)
        pair_threshold_pct = float(dashboard.get("pair_threshold_pct", threshold_pct) or threshold_pct)
        ax.axvline(threshold_pct, color="#444444", linestyle="--", linewidth=1)
        ax.text(
            min(103.0, max(1.0, threshold_pct + 1.0)),
            -0.72,
            f"tipo {threshold_pct:.1f}%",
            fontsize=8,
            va="center",
            ha="left",
            color="#333333",
        )
        unreliable_rows = [
            idx for idx, feature in enumerate(features)
            if str(feature.get("prevalent_interaction", "")) in {"", CLASS_UNRELIABLE}
        ]
        if unreliable_rows:
            ax.scatter(
                [102.0] * len(unreliable_rows),
                unreliable_rows,
                marker="x",
                color="black",
                s=36,
                zorder=5,
                clip_on=False,
            )

        ax.set_xlim(0, 105)
        ax.set_xlabel("Frequência de atribuição de cada interação (%)")
        ax.set_ylabel("ID da feature")
        ax.set_yticks(y_pos)
        ax.set_yticklabels([self._fp_feature_label(feature) for feature in features])
        ax.invert_yaxis()
        ax.set_title(
            "Frequência de atribuição da interação prevalente nas features importantes\n"
            f"limiar tipo de interação = {threshold_pct:.2f}% | limiar par interação/resíduo = {pair_threshold_pct:.2f}%"
        )
        handles, labels_text = ax.get_legend_handles_labels()
        if handles:
            legend_cols = min(3, max(1, len(handles)))
            ax.legend(
                handles,
                labels_text,
                loc="upper center",
                bbox_to_anchor=(0.5, -0.18),
                ncol=legend_cols,
                fontsize=8,
                frameon=False,
            )
            self.fp_interaction_assign_fig.subplots_adjust(
                bottom=min(0.48, 0.16 + (0.06 * math.ceil(len(handles) / legend_cols)))
            )
        else:
            self.fp_interaction_assign_fig.tight_layout()
        self.fp_interaction_assign_canvas.draw()

    def _render_fp_interaction_summary_plot(self, dashboard: dict) -> None:
        all_features = [
            feature
            for feature in (dashboard.get("important_features", []) or [])
            if str(feature.get("prevalent_pair", "")) not in {"", CLASS_UNRELIABLE}
        ]
        legend_types = []
        for feature in all_features:
            interaction_name = str(feature.get("prevalent_interaction", "")).strip()
            if interaction_name and interaction_name not in legend_types:
                legend_types.append(interaction_name)
        legend_types = sorted(legend_types, key=interaction_priority_key)
        features = [
            feature
            for feature in all_features
            if str(feature.get("prevalent_interaction", "")).strip() not in self._fp_interaction_hidden_types
        ]
        self._resize_canvas(
            self.fp_interaction_fig,
            self.fp_interaction_canvas,
            width_in=max(10.0, 0.6 * max(1, len(features)) + 3.5),
            height_in=6.2,
        )
        self.fp_interaction_fig.clear()
        ax = self.fp_interaction_fig.add_subplot(111)

        if not all_features:
            reason = str(dashboard.get("detail_error", "") or "").strip()
            if not reason and not dashboard.get("detail_available"):
                reason = "O artefato fp_detail não está disponível para este IFP."
            message = "Sem interações prevalentes confiáveis nas features importantes."
            if reason:
                message += "\n" + textwrap.fill(reason, width=80)
            ax.text(0.5, 0.5, message, ha="center", va="center")
            self.fp_interaction_fig.tight_layout()
            self.fp_interaction_canvas.draw()
            return
        if not features:
            ax.text(0.5, 0.5, "Todas as interações estão ocultas pela legenda.", ha="center", va="center")

        x = list(range(len(features)))
        heights = [len(feature.get("prevalent_pair_entries", []) or []) for feature in features]
        colors = [
            get_interaction_color(str(feature.get("prevalent_interaction", "")))
            for feature in features
        ]
        bars = ax.bar(x, heights, color=colors) if features else []
        ax.set_xticks(x)
        ax.set_xticklabels([self._fp_feature_label(feature) for feature in features], rotation=45, ha="right")
        ax.set_ylabel("Número de ligantes")
        ax.set_xlabel("ID da feature")
        ax.set_title("Interação e resíduo prevalentes nas features importantes")
        ax.set_ylim(0, max(heights) * 1.18 if heights else 1.0)

        for bar, feature in zip(bars, features):
            residue = str(feature.get("prevalent_residue", "")).strip()
            if residue and residue != CLASS_UNRELIABLE:
                residue = format_residue_label(residue)
            pair_pct = float(feature.get("prevalent_pair_pct", 0.0) or 0.0)
            label = residue or str(feature.get("prevalent_interaction", ""))
            ax.text(
                bar.get_x() + bar.get_width() / 2.0,
                bar.get_height() + 0.2,
                f"{label}\n{pair_pct:.2f}%",
                ha="center",
                va="bottom",
                fontsize=8,
            )

        if legend_types:
            handles = [
                Patch(
                    facecolor=get_interaction_color(name),
                    edgecolor="none",
                    label=name,
                    alpha=0.28 if name in self._fp_interaction_hidden_types else 1.0,
                )
                for name in legend_types
            ]
            legend_cols = min(3, max(1, len(handles)))
            legend = ax.legend(
                handles=handles,
                labels=legend_types,
                loc="upper center",
                bbox_to_anchor=(0.5, -0.18),
                fontsize=8,
                ncol=legend_cols,
                frameon=False,
                title="Clique para ocultar/mostrar",
            )
            self._install_fp_interaction_legend_toggle(legend, legend_types)
            self.fp_interaction_fig.subplots_adjust(
                bottom=min(0.44, 0.16 + (0.06 * math.ceil(len(handles) / legend_cols)))
            )
        else:
            self.fp_interaction_fig.tight_layout()
        self.fp_interaction_canvas.draw()

    def _render_fp_interaction_heatmap(self, dashboard: dict) -> None:
        features = [
            feature
            for feature in (dashboard.get("important_features", []) or [])
            if str(feature.get("prevalent_interaction", "")) not in {"", CLASS_UNRELIABLE}
        ]
        entries = list(dashboard.get("entry_labels", []) or [])
        entries, _row_order = self._order_entries_by_activity(entries)
        self._resize_canvas(
            self.fp_interaction_heatmap_fig,
            self.fp_interaction_heatmap_canvas,
            width_in=max(9.5, 0.34 * max(1, len(features)) + 3.5),
            height_in=max(7.5, 0.12 * self._render_entry_count(len(entries)) + 3.0),
        )
        self.fp_interaction_heatmap_fig.clear()
        ax = self.fp_interaction_heatmap_fig.add_subplot(111)

        if not features or not entries:
            reason = str(dashboard.get("detail_error", "") or "").strip()
            if not reason and not dashboard.get("detail_available"):
                reason = "O artefato fp_detail não está disponível para este IFP."
            message = "Sem interações prevalentes para gerar o mapa de calor."
            if reason:
                message += "\n" + textwrap.fill(reason, width=80)
            ax.text(0.5, 0.5, message, ha="center", va="center")
            self.fp_interaction_heatmap_fig.tight_layout()
            self.fp_interaction_heatmap_canvas.draw()
            return

        interaction_names: list[str] = []
        for feature in features:
            name = str(feature.get("prevalent_interaction", "")).strip()
            if name and name not in interaction_names:
                interaction_names.append(name)
        interaction_names = sorted(interaction_names, key=interaction_priority_key)
        interaction_to_id = {name: idx for idx, name in enumerate(interaction_names, start=1)}
        matrix = np.zeros((len(entries), len(features)), dtype=int)
        for col_idx, feature in enumerate(features):
            interaction_name = str(feature.get("prevalent_interaction", "")).strip()
            interaction_id = interaction_to_id.get(interaction_name, 0)
            entry_names = set(feature.get("prevalent_interaction_entries", []) or [])
            if not entry_names:
                entry_names = set((feature.get("entry_counts") or {}).keys())
            for row_idx, entry_name in enumerate(entries):
                if entry_name in entry_names:
                    matrix[row_idx, col_idx] = interaction_id

        colors = ["#3f3f3f"] + [get_interaction_color(name) for name in interaction_names]
        cmap = ListedColormap(colors)
        norm = BoundaryNorm(range(len(colors) + 1), cmap.N)
        im = ax.imshow(matrix, cmap=cmap, norm=norm, aspect="auto", interpolation="nearest")
        ax.set_xticks(list(range(len(features))))
        ax.set_xticklabels([self._fp_feature_label(feature) for feature in features], rotation=90, fontsize=7)
        _apply_tick_labels(
            ax,
            [self._display_ligand_name(entry) for entry in entries],
            axis="y",
            max_labels=28,
            rotation=0,
            ligand_axis=True,
        )
        self._color_ticklabels_by_entry_group(ax, entries, axis="y")
        if len(entries) <= 220:
            ax.set_ylabel("Ligantes")
        ax.set_xlabel("ID da feature")
        ax.set_title("Interações prevalentes das features importantes por ligante")

        cbar = self.fp_interaction_heatmap_fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
        cbar.set_ticks([idx + 0.5 for idx in range(len(interaction_names) + 1)])
        cbar.set_ticklabels(["Ausente"] + interaction_names)
        self.fp_interaction_heatmap_fig.tight_layout()
        self.fp_interaction_heatmap_canvas.draw()

    def _sync_fp_session_types(self) -> None:
        current = self._combo_current_data(self.cb_fp_session_type) if self.cb_fp_session_type.count() else None
        self.cb_fp_session_type.blockSignals(True)
        self.cb_fp_session_type.clear()
        for ifp_type, artifact in sorted(self._fp_artifacts.items()):
            label = artifact.get("ifp_label", ifp_type)
            self.cb_fp_session_type.addItem(f"{label} ({ifp_type})", ifp_type)
        idx = self.cb_fp_session_type.findData(current)
        self.cb_fp_session_type.setCurrentIndex(idx if idx >= 0 else 0)
        self.cb_fp_session_type.blockSignals(False)
        self._populate_fp_session_features()

    def _populate_fp_session_features(self) -> None:
        ifp_type = self._combo_current_data(self.cb_fp_session_type) if self.cb_fp_session_type.count() else None
        artifact = self._fp_artifacts.get(ifp_type) if ifp_type else None
        current = self._combo_current_data(self.cb_fp_session_feature) if self.cb_fp_session_feature.count() else None
        self.cb_fp_session_feature.blockSignals(True)
        self.cb_fp_session_feature.clear()
        if artifact:
            feature_rows: dict[int, dict] = {}
            for rows in (artifact.get("entry_index") or {}).values():
                for row in rows:
                    if "feature_id" not in row:
                        continue
                    feature_rows.setdefault(int(row["feature_id"]), row)
            for feature_id in sorted(feature_rows):
                row = feature_rows[feature_id]
                dominant_nature = normalize_fp_class_name(
                    row.get("dominant_nature")
                    or next(iter(row.get("nature_tags", []) or []), "feature")
                )
                levels = self._format_shell_levels(row.get("shell_levels"), row.get("shell_level_breakdown"))
                label = f"{feature_id} - {dominant_nature}" + (f" | L{levels}" if levels else "")
                self.cb_fp_session_feature.addItem(label, feature_id)
        idx = self.cb_fp_session_feature.findData(current)
        self.cb_fp_session_feature.setCurrentIndex(idx if idx >= 0 else 0)
        self.cb_fp_session_feature.blockSignals(False)
        self._populate_fp_session_entries()

    def _populate_fp_session_entries(self) -> None:
        ifp_type = self._combo_current_data(self.cb_fp_session_type) if self.cb_fp_session_type.count() else None
        artifact = self._fp_artifacts.get(ifp_type) if ifp_type else None
        feature_id = self._combo_current_data(self.cb_fp_session_feature) if self.cb_fp_session_feature.count() else None
        current = self._combo_current_data(self.cb_fp_session_entry) if self.cb_fp_session_entry.count() else None
        self.cb_fp_session_entry.blockSignals(True)
        self.cb_fp_session_entry.clear()
        if artifact and feature_id is not None:
            ranked_entries = []
            for entry_name, rows in sorted((artifact.get("entry_index") or {}).items()):
                for row in rows:
                    if int(row.get("feature_id", -1)) == int(feature_id):
                        ranked_entries.append((int(row.get("count", 0) or 0), entry_name))
                        break
            for count, entry_name in sorted(ranked_entries, key=lambda item: (-item[0], item[1])):
                self.cb_fp_session_entry.addItem(f"{entry_name} (count={count})", entry_name)
        idx = self.cb_fp_session_entry.findData(current)
        self.cb_fp_session_entry.setCurrentIndex(idx if idx >= 0 else 0)
        self.cb_fp_session_entry.blockSignals(False)
        self._update_fp_session_info()

    def _update_fp_session_info(self) -> None:
        ifp_type = self._combo_current_data(self.cb_fp_session_type) if self.cb_fp_session_type.count() else None
        artifact = self._fp_artifacts.get(ifp_type) if ifp_type else None
        entry_name = self._combo_current_data(self.cb_fp_session_entry) if self.cb_fp_session_entry.count() else None
        feature_id = self._combo_current_data(self.cb_fp_session_feature) if self.cb_fp_session_feature.count() else None
        if not artifact or not entry_name or feature_id is None:
            self.fp_session_info.setText("-")
            return

        row = None
        for item in (artifact.get("entry_index") or {}).get(entry_name, []):
            if int(item.get("feature_id", -1)) == int(feature_id):
                row = item
                break

        if row is None:
            self.fp_session_info.setText("-")
            return

        natures = ", ".join(row.get("nature_tags", []))
        original = ", ".join(str(value) for value in row.get("original_features", []))
        levels = self._format_shell_levels(row.get("shell_levels"), row.get("shell_level_breakdown"))
        collision_levels = self._format_shell_levels(
            row.get("collision_shell_levels"),
            row.get("collision_level_breakdown"),
        )
        self.fp_session_info.setText(
            f"Feature {feature_id} em {entry_name}. Natureza: {natures or '-'}"
            + (f" | Shells originais: {original}" if original else "")
            + (f" | Níveis shell: {levels}" if levels else "")
            + (f" | Níveis em colisão: {collision_levels}" if collision_levels else "")
        )

    def _sync_fp_session_from_table(self) -> None:
        selected = self.fp_analysis_table.selectedItems()
        if not selected:
            return
        feature_item = self.fp_analysis_table.item(selected[0].row(), 0)
        if feature_item is None:
            return
        feature_id = feature_item.data(Qt.ItemDataRole.UserRole)
        if feature_id is None:
            return

        ifp_type = self.cb_fp_analysis_type.currentData() if self.cb_fp_analysis_type.count() else None
        idx_type = self.cb_fp_session_type.findData(ifp_type)
        if idx_type >= 0:
            self.cb_fp_session_type.setCurrentIndex(idx_type)

        artifact = self._fp_artifacts.get(ifp_type) if ifp_type else None
        if not artifact:
            return
        entry_index = artifact.get("entry_index") or {}
        idx_feature = self.cb_fp_session_feature.findData(int(feature_id))
        if idx_feature >= 0:
            self.cb_fp_session_feature.setCurrentIndex(idx_feature)
        ranked_entries = []
        for candidate, rows in sorted(entry_index.items()):
            for row in rows:
                if int(row.get("feature_id", -1)) == int(feature_id):
                    ranked_entries.append((int(row.get("count", 0) or 0), candidate))
                    break
        if ranked_entries:
            best_entry = sorted(ranked_entries, key=lambda item: (-item[0], item[1]))[0][1]
            idx_entry = self.cb_fp_session_entry.findData(best_entry)
            if idx_entry >= 0:
                self.cb_fp_session_entry.setCurrentIndex(idx_entry)

    def _generate_fp_session(self) -> None:
        wd = self._current_wd()
        if not wd:
            return
        if not self.py_exe:
            QMessageBox.warning(self, "luna-env", "LUNA não detectado. Verifique a aba Setup.")
            return

        ifp_type = self._combo_current_data(self.cb_fp_session_type) if self.cb_fp_session_type.count() else None
        entry_name = self._combo_current_data(self.cb_fp_session_entry) if self.cb_fp_session_entry.count() else None
        feature_id = self._combo_current_data(self.cb_fp_session_feature) if self.cb_fp_session_feature.count() else None
        if not ifp_type or not entry_name or feature_id is None:
            QMessageBox.information(self, "Sessão FP", "Escolha um tipo, um ligante e um fingerprint.")
            return

        out_dir = wd / "results" / "fingerprints" / "pse" / str(ifp_type)
        out_dir.mkdir(parents=True, exist_ok=True)
        out_path = out_dir / f"{_safe_name(entry_name)}__feature_{int(feature_id)}.pse"
        self.fp_session_status.setText("Gerando sessão...")
        self.fp_session_status.repaint()
        result = generate_fp_session(
            self.py_exe,
            str(wd),
            str(ifp_type),
            str(entry_name),
            int(feature_id),
            str(out_path),
        )
        if "error" in result:
            self.fp_session_status.setText("Erro")
            QMessageBox.critical(self, "Erro ao gerar sessão", result["error"])
            return

        out_path = Path(str(result.get("output") or out_path))
        self._generated_fp_session = str(out_path)
        source = str(result.get("source") or "")
        shell_labels = int(result.get("shell_labels", 0) or 0)
        label_note = f"; {shell_labels} shells numerados" if shell_labels else ""
        if source.startswith("live_project"):
            self.fp_session_status.setText(f"Sessão salva em {out_path.name} (shells regenerados{label_note})")
        elif source.startswith("cached_payload"):
            self.fp_session_status.setText(f"Sessão salva em {out_path.name} (cache legado{label_note})")
        else:
            self.fp_session_status.setText(f"Sessão salva em {out_path.name}{label_note}")
        self._load_existing_fp_sessions(wd)
        self._select_fp_session(str(out_path))

    def _load_existing_fp_sessions(self, wd: Path) -> None:
        self.fp_session_list.clear()
        fp_pse_root = wd / "results" / "fingerprints" / "pse"
        if not fp_pse_root.exists():
            return
        paths = sorted(fp_pse_root.glob("**/*.pse"))
        for path in paths:
            item = QListWidgetItem(path.name)
            item.setData(Qt.ItemDataRole.UserRole, str(path))
            self.fp_session_list.addItem(item)

    def _select_fp_session(self, path: str) -> None:
        for index in range(self.fp_session_list.count()):
            item = self.fp_session_list.item(index)
            if item.data(Qt.ItemDataRole.UserRole) == path:
                self.fp_session_list.setCurrentItem(item)
                break

    def _open_selected_fp_session(self) -> None:
        item = self.fp_session_list.currentItem()
        if item is None and self._generated_fp_session:
            self._open_path(self._generated_fp_session)
            return
        if item is None:
            QMessageBox.information(self, "Sessão FP", "Nenhuma sessão de fingerprint foi selecionada.")
            return
        self._open_path(item.data(Qt.ItemDataRole.UserRole))

    def _open_path(self, path: str) -> None:
        try:
            launch_pse_session(path, self.py_exe)
        except Exception as exc:
            QMessageBox.critical(self, "Erro ao abrir PyMOL", str(exc))

    def _save_report_figure(self, fig, path: Path) -> Path | None:
        if not HAS_MPL or fig is None or not getattr(fig, "axes", None):
            return None
        old_size = tuple(fig.get_size_inches())
        old_dpi = fig.dpi
        try:
            translate_figure(fig)
            width_in = max(float(old_size[0]), self._REPORT_FIGURE_WIDTH_IN)
            height_in = max(float(old_size[1]), self._REPORT_FIGURE_HEIGHT_IN)
            fig.set_size_inches(width_in, height_in, forward=False)
            fig.savefig(
                path,
                dpi=self._REPORT_FIGURE_DPI,
                bbox_inches="tight",
                pad_inches=0.18,
                facecolor=fig.get_facecolor(),
            )
        except Exception:
            return None
        finally:
            try:
                fig.set_size_inches(old_size, forward=False)
                fig.set_dpi(old_dpi)
            except Exception:
                pass
        return path if path.exists() else None

    def _save_report_stats_overview(self, path: Path) -> Path | None:
        """Save the statistics chart in the project-wide scope for PDF reports."""
        if not HAS_MPL or not self._last_analysis:
            return None
        old_index = None
        old_hidden = set(getattr(self, "_stats_hidden_interactions", set()))
        combo = getattr(self, "cb_stats_scope", None)
        if combo is not None:
            old_index = combo.currentIndex()
        try:
            self._stats_hidden_interactions.clear()
            if combo is not None:
                idx = combo.findData("__all__")
                if idx >= 0:
                    combo.blockSignals(True)
                    combo.setCurrentIndex(idx)
                    combo.blockSignals(False)
            self._render_stats_chart(self._last_analysis)
            return self._save_report_figure(getattr(self, "st_fig", None), path)
        finally:
            self._stats_hidden_interactions = old_hidden
            if combo is not None and old_index is not None and 0 <= old_index < combo.count():
                combo.blockSignals(True)
                combo.setCurrentIndex(old_index)
                combo.blockSignals(False)
            if self._last_analysis:
                self._render_stats_chart(self._last_analysis)

    def export_pdf_report(self) -> None:
        wd = self._current_wd()
        if not wd:
            return
        out, _ = QFileDialog.getSaveFileName(
            self,
            "Salvar relatório PDF",
            str(wd / "luna_report.pdf"),
            "PDF (*.pdf)",
        )
        if not out:
            return

        if not self._last_analysis:
            self.compute_stats()
        if not self._last_analysis:
            return

        similarity_png = self._save_report_figure(getattr(self, "fig", None), wd / "_report_pdf_similarity.png")
        interactions_png = self._save_report_stats_overview(wd / "_report_pdf_interactions_all.png")
        cluster_png = None
        if getattr(self, "_cluster_result", None):
            cluster_png = self._save_report_figure(getattr(self, "cluster_fig", None), wd / "_report_pdf_clusters.png")

        extra_images: list[tuple[str, Path, str]] = []
        figure_specs = [
            ("hm_fig", "Mapa de calor por tipo", "Linhas representam ligantes e colunas representam resíduos. A intensidade da célula indica quantas vezes o tipo de interação selecionado aparece naquele par; colunas densas destacam resíduos recorrentes."),
            ("hm_all_fig", "Mapa de calor completo ligantes x resíduos", "Cada célula pode conter várias faixas de cor, uma por classe de interação. Esse gráfico mostra quando um mesmo resíduo participa de mecanismos químicos diferentes entre ligantes."),
            ("fp_class_fig", "FP - classes importantes", "Resume as classes atribuídas às features que passaram pelo corte de p-value. Classes dominantes indicam a natureza estrutural das features que mais influenciam o modelo."),
            ("fp_assign_fig", "FP - frequência por classe", "Mostra a composição de classe das features importantes. Barras mistas sugerem features com colisões ou natureza ambígua; barras concentradas sugerem assinatura mais interpretável."),
            ("fp_cover_fig", "FP - cobertura e importância", "Compara quantos ligantes possuem cada feature com a importância estimada pelo modelo. Features com alta importância e cobertura moderada costumam ser boas candidatas para inspeção."),
            ("fp_heatmap_fig", "FP - mapa de calor de importância", "Relaciona features importantes e classes de fingerprint. Tons mais intensos indicam maior peso relativo; use os IDs das features para gerar sessão PyMOL correspondente."),
            ("fp_interaction_assign_fig", "FP - frequência de interações", "Mostra a distribuição das interações prevalentes nos shells das features importantes. Ajuda a ligar uma feature abstrata a contatos químicos observáveis."),
            ("fp_interaction_fig", "FP - interações prevalentes", "Resume quais tipos de interação aparecem como dominantes nas features importantes, após aplicar o limiar configurado por z-score ou Otsu."),
            ("fp_interaction_heatmap_fig", "FP - mapa de calor de interações", "Cruza features importantes com interações prevalentes dos shells reais do LUNA. Ele revela se diferentes features importantes apontam para a mesma família de contatos."),
        ]
        for attr, title, caption in figure_specs:
            saved = self._save_report_figure(
                getattr(self, attr, None),
                wd / f"_report_pdf_{attr}.png",
            )
            if saved:
                extra_images.append((title, saved, caption))

        cluster_items = None
        if getattr(self, "_cluster_result", None):
            cluster_items = [
                (label, cluster_id)
                for label, cluster_id, _ in cluster_rows(self._cluster_result)
            ]

        try:
            save_pdf_report_isolated(
                out,
                cfg=self.cfg,
                analysis=self._last_analysis,
                heatmap_png=similarity_png,
                interactions_png=interactions_png,
                cluster_png=cluster_png,
                clusters=cluster_items,
                fp_dashboards=getattr(self, "_fp_dashboards", {}),
                extra_images=extra_images,
                progress_callback=QApplication.processEvents,
            )
        except Exception as exc:
            QMessageBox.critical(self, "Erro ao gerar PDF", str(exc))
            return
        QMessageBox.information(self, "Relatório PDF salvo", out)

    def _current_figure(self):
        current = self.inner.currentWidget()
        if HAS_MPL and current is self.complete_heatmap_tab:
            return self.hm_all_fig, "heatmap_completo"
        if (
            HAS_MPL
            and current is self.fp_analysis_tab
            and self.fp_plot_tabs is not None
            and self._fp_plot_canvases_ready
        ):
            current_plot = self.fp_plot_tabs.currentWidget()
            if current_plot is self.fp_class_tab:
                return self.fp_class_fig, "fp_classes_importantes"
            if current_plot is self.fp_assign_tab:
                return self.fp_assign_fig, "fp_frequencia_por_classe"
            if current_plot is self.fp_cover_tab:
                return self.fp_cover_fig, "fp_cobertura_importancia"
            if current_plot is self.fp_heatmap_tab:
                return self.fp_heatmap_fig, "fp_heatmap_importancia"
            if current_plot is self.fp_interaction_assign_tab:
                return self.fp_interaction_assign_fig, "fp_frequencia_interacoes"
            if current_plot is self.fp_interaction_tab:
                return self.fp_interaction_fig, "fp_interacoes_prevalentes"
            if current_plot is self.fp_interaction_heatmap_tab:
                return self.fp_interaction_heatmap_fig, "fp_heatmap_interacoes"
        return super()._current_figure()


def _safe_name(value: str) -> str:
    chars = []
    for char in str(value):
        chars.append(char if char.isalnum() or char in "._-" else "_")
    cleaned = "".join(chars).strip("_")
    return cleaned or "entry"
