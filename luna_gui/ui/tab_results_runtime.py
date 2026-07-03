"""Runtime-aware Results tab with cached analytics and richer fingerprint views."""
from __future__ import annotations

import os
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import (
    QComboBox,
    QLabel,
    QListWidget,
    QListWidgetItem,
    QMessageBox,
    QPushButton,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
    QWidget,
    QHBoxLayout,
)

from ..core.analysis_runtime import generate_fp_session, run_analysis, run_residue_matrix
from ..core.results_analysis import (
    build_complete_heatmap,
    load_analysis_summary,
    load_fp_analysis_artifacts,
    load_residue_matrix_artifact,
)
from .tab_results_enhanced import HAS_MPL, ResultsTab as EnhancedResultsTab, _apply_tick_labels

if HAS_MPL:
    from matplotlib.colors import BoundaryNorm, ListedColormap


class ResultsTab(EnhancedResultsTab):
    def __init__(self, cfg) -> None:
        super().__init__(cfg)
        self._install_stats_scope_control()

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
            "Alterna entre o resumo total do projeto e a distribuição por um ligante específico."
        )
        self.cb_stats_scope.currentIndexChanged.connect(self._render_cached_stats_chart)
        st_ctrl.insertWidget(1, QLabel("Visão:"))
        st_ctrl.insertWidget(2, self.cb_stats_scope)

    def load_all(self) -> None:
        super().load_all()
        wd = self._current_wd()
        if not wd:
            return
        self._load_cached_stats(wd)
        self._load_cached_residue_matrix(wd)

    def compute_stats(self) -> None:
        wd = self._current_wd()
        if not wd:
            return

        cached = load_analysis_summary(wd)
        result = cached
        if result is None:
            if not self.py_exe:
                QMessageBox.warning(self, "luna-env", "LUNA não detectado. Verifique a aba Setup.")
                return
            self.st_status.setText("Processando... (pode levar alguns minutos)")
            self.st_status.repaint()
            result = run_analysis(self.py_exe, str(wd))

        if "error" in result:
            self.st_status.setText("Erro")
            QMessageBox.critical(self, "Erro na análise", result["error"])
            return

        self._last_analysis = result
        processed = len(result.get("entry_interaction_counts", {})) or result.get("entries", 0)
        self.st_status.setText(f"{processed} entradas processadas")
        self._populate_stats_scope(result)
        self._render_cached_stats_chart()

    def compute_residue_matrix(self) -> None:
        wd = self._current_wd()
        if not wd:
            return

        cached = load_residue_matrix_artifact(wd)
        result = cached
        if result is None:
            if not self.py_exe:
                QMessageBox.warning(self, "luna-env", "LUNA não detectado. Veja a aba Setup.")
                return
            self.hm_status.setText("Processando...")
            self.hm_status.repaint()
            result = run_residue_matrix(self.py_exe, str(wd))

        if "error" in result:
            self.hm_status.setText("Erro")
            QMessageBox.critical(self, "Erro na análise", result["error"])
            return

        self._residue_matrix = result
        types = result.get("interaction_types", [])
        self.hm_status.setText(f"{len(result.get('entries', []))} entradas · {len(types)} tipos")
        self.cb_itype.blockSignals(True)
        self.cb_itype.clear()
        self.cb_itype.addItems(types)
        self.cb_itype.blockSignals(False)
        if types:
            self._render_residue_heatmap()
        elif HAS_MPL:
            self.hm_fig.clear()
            ax = self.hm_fig.add_subplot(111)
            ax.text(0.5, 0.5, "Sem matriz de resíduos disponível", ha="center", va="center")
            self.hm_fig.tight_layout()
            self.hm_canvas.draw()

    def _load_cached_stats(self, wd: Path) -> None:
        cached = load_analysis_summary(wd)
        if cached is None:
            return
        self._last_analysis = cached
        processed = len(cached.get("entry_interaction_counts", {})) or cached.get("entries", 0)
        self.st_status.setText(f"{processed} entradas processadas")
        self._populate_stats_scope(cached)
        self._render_cached_stats_chart()

    def _load_cached_residue_matrix(self, wd: Path) -> None:
        cached = load_residue_matrix_artifact(wd)
        if cached is None:
            return
        self._residue_matrix = cached
        types = cached.get("interaction_types", [])
        self.hm_status.setText(f"{len(cached.get('entries', []))} entradas · {len(types)} tipos")
        self.cb_itype.blockSignals(True)
        current = self.cb_itype.currentText()
        self.cb_itype.clear()
        self.cb_itype.addItems(types)
        if current:
            idx = self.cb_itype.findText(current)
            if idx >= 0:
                self.cb_itype.setCurrentIndex(idx)
        self.cb_itype.blockSignals(False)
        if types:
            self._render_residue_heatmap()

    def _populate_stats_scope(self, result: dict) -> None:
        if not hasattr(self, "cb_stats_scope"):
            return
        current = self.cb_stats_scope.currentData()
        self.cb_stats_scope.blockSignals(True)
        self.cb_stats_scope.clear()
        self.cb_stats_scope.addItem("Totais do projeto", "__all__")
        for ligand_name in sorted((result.get("entry_interaction_counts") or {}).keys()):
            self.cb_stats_scope.addItem(ligand_name, ligand_name)
        idx = self.cb_stats_scope.findData(current)
        self.cb_stats_scope.setCurrentIndex(idx if idx >= 0 else 0)
        self.cb_stats_scope.blockSignals(False)

    def _render_cached_stats_chart(self) -> None:
        if HAS_MPL and self._last_analysis:
            self._render_stats_chart(self._last_analysis)

    def _render_stats_chart(self, result: dict) -> None:
        self.st_fig.clear()
        ax = self.st_fig.add_subplot(111)

        scope = "__all__"
        if hasattr(self, "cb_stats_scope"):
            scope = self.cb_stats_scope.currentData() or "__all__"

        if scope == "__all__":
            counts = result.get("interaction_counts", {}) or {}
            title = "Contagem por tipo de interação"
            xlabel = "Total (todas as entradas)"
        else:
            counts = (result.get("entry_interaction_counts", {}) or {}).get(scope, {}) or {}
            title = f"Interações por tipo — {scope}"
            xlabel = "Total neste ligante"

        if not counts:
            ax.text(0.5, 0.5, "Sem dados de interação", ha="center", va="center")
            self.st_fig.tight_layout()
            self.st_canvas.draw()
            return

        items = sorted(counts.items(), key=lambda item: (-item[1], item[0]))
        labels = [label for label, _value in items]
        values = [value for _label, value in items]
        ax.barh(labels, values, color="#c8693a")
        ax.invert_yaxis()
        ax.set_xlabel(xlabel)
        ax.set_title(title)
        self.st_fig.tight_layout()
        self.st_canvas.draw()
