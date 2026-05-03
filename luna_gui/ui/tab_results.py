"""Results tab — preview fingerprint CSV, render similarity heatmap, open PSE."""
from __future__ import annotations

import csv
import shutil
import subprocess
import sys
from pathlib import Path

from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QPushButton, QLineEdit, QLabel,
    QTableWidget, QTableWidgetItem, QFileDialog, QListWidget, QListWidgetItem,
    QMessageBox, QTabWidget, QGroupBox, QSpinBox, QComboBox,
)

from ..core.project import ProjectConfig
from ..core.analysis_helper import run_analysis, run_residue_matrix
from ..core.report import save_report
from ..i18n import translate_figure

# matplotlib is optional — heatmap is disabled gracefully if missing.
try:
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
    HAS_MPL = True
except Exception:
    HAS_MPL = False


class ResultsTab(QWidget):
    def __init__(self, cfg: ProjectConfig) -> None:
        super().__init__()
        self.cfg = cfg
        self.py_exe: str = ""
        self._last_analysis: dict = {}

        layout = QVBoxLayout(self)

        # Workdir picker
        wd_row = QHBoxLayout()
        wd_row.addWidget(QLabel("Workdir:"))
        self.wd_edit = QLineEdit()
        self.wd_edit.setPlaceholderText("(usa o workdir do projeto atual)")
        btn_wd = QPushButton("Procurar...")
        btn_wd.clicked.connect(self._pick_wd)
        btn_load = QPushButton("Carregar resultados")
        btn_load.clicked.connect(self.load_all)
        btn_report = QPushButton("Exportar relatório HTML")
        btn_report.clicked.connect(self.export_report)
        wd_row.addWidget(self.wd_edit, 1)
        wd_row.addWidget(btn_wd)
        wd_row.addWidget(btn_load)
        wd_row.addWidget(btn_report)
        layout.addLayout(wd_row)

        # Inner tabs: Fingerprints | Similarity | PyMOL sessions
        self.inner = QTabWidget()
        layout.addWidget(self.inner, 1)

        # --- Fingerprints preview ---
        fp_w = QWidget(); fp_l = QVBoxLayout(fp_w)
        fp_ctrl = QHBoxLayout()
        fp_ctrl.addWidget(QLabel("Linhas a exibir:"))
        self.fp_rows = QSpinBox(); self.fp_rows.setRange(10, 10000); self.fp_rows.setValue(200)
        fp_ctrl.addWidget(self.fp_rows)
        self.fp_path_label = QLabel("—")
        self.fp_path_label.setStyleSheet("color:#666;")
        fp_ctrl.addWidget(self.fp_path_label, 1)
        fp_l.addLayout(fp_ctrl)
        self.fp_table = QTableWidget()
        fp_l.addWidget(self.fp_table, 1)
        self.inner.addTab(fp_w, "Fingerprints")

        # --- Similarity heatmap ---
        sm_w = QWidget(); sm_l = QVBoxLayout(sm_w)
        self.sm_path_label = QLabel("—")
        self.sm_path_label.setStyleSheet("color:#666;")
        sm_l.addWidget(self.sm_path_label)
        if HAS_MPL:
            self.fig = Figure(figsize=(6, 5))
            self.canvas = FigureCanvas(self.fig)
            sm_l.addWidget(self.canvas, 1)
        else:
            sm_l.addWidget(QLabel(
                "matplotlib não está instalado. Instale com: pip install matplotlib"
            ))
        self.inner.addTab(sm_w, "Matriz de similaridade")

        # --- Statistics ---
        st_w = QWidget(); st_l = QVBoxLayout(st_w)
        st_ctrl = QHBoxLayout()
        btn_st = QPushButton("Calcular estatísticas (usa luna-env)")
        btn_st.clicked.connect(self.compute_stats)
        self.st_status = QLabel("—")
        self.st_status.setStyleSheet("color:#666;")
        st_ctrl.addWidget(btn_st)
        st_ctrl.addWidget(self.st_status, 1)
        st_l.addLayout(st_ctrl)
        if HAS_MPL:
            self.st_fig = Figure(figsize=(6, 4))
            self.st_canvas = FigureCanvas(self.st_fig)
            st_l.addWidget(self.st_canvas, 1)
        else:
            st_l.addWidget(QLabel("matplotlib não está instalado."))
        self.inner.addTab(st_w, "Estatísticas")

        # --- Heatmap per interaction type (T1.4) ---
        hm_w = QWidget(); hm_l = QVBoxLayout(hm_w)
        hm_ctrl = QHBoxLayout()
        btn_hm = QPushButton("Calcular heatmap (usa luna-env)")
        btn_hm.clicked.connect(self.compute_residue_matrix)
        hm_ctrl.addWidget(btn_hm)
        hm_ctrl.addWidget(QLabel("Tipo:"))
        self.cb_itype = QComboBox()
        self.cb_itype.currentIndexChanged.connect(self._render_residue_heatmap)
        hm_ctrl.addWidget(self.cb_itype, 1)
        self.hm_status = QLabel("—"); self.hm_status.setStyleSheet("color:#666;")
        hm_ctrl.addWidget(self.hm_status)
        hm_l.addLayout(hm_ctrl)
        if HAS_MPL:
            self.hm_fig = Figure(figsize=(8, 5))
            self.hm_canvas = FigureCanvas(self.hm_fig)
            hm_l.addWidget(self.hm_canvas, 1)
        else:
            hm_l.addWidget(QLabel("matplotlib não está instalado."))
        self._residue_matrix: dict = {}
        self.inner.addTab(hm_w, "Heatmap por tipo")

        # --- PyMOL sessions ---
        pse_w = QWidget(); pse_l = QVBoxLayout(pse_w)
        pse_ctrl = QHBoxLayout()
        btn_open = QPushButton("Abrir no PyMOL")
        btn_open.clicked.connect(self._open_pse)
        pse_ctrl.addWidget(btn_open)
        pse_ctrl.addStretch()
        pse_l.addLayout(pse_ctrl)
        self.pse_list = QListWidget()
        self.pse_list.itemDoubleClicked.connect(lambda _: self._open_pse())
        pse_l.addWidget(self.pse_list, 1)
        self.inner.addTab(pse_w, "Sessões PyMOL")

    # ---- helpers ----
    def _pick_wd(self) -> None:
        d = QFileDialog.getExistingDirectory(self, "Workdir do projeto LUNA")
        if d:
            self.wd_edit.setText(d)

    def _current_wd(self) -> Path | None:
        wd = self.wd_edit.text().strip() or self.cfg.workdir
        if not wd:
            QMessageBox.warning(self, "Workdir", "Defina um workdir primeiro.")
            return None
        p = Path(wd)
        if not p.exists():
            QMessageBox.warning(self, "Workdir", f"Diretório não existe:\n{p}")
            return None
        return p

    def load_all(self) -> None:
        wd = self._current_wd()
        if not wd:
            return
        self._load_fingerprints(wd)
        self._load_sim_matrix(wd)
        self._load_pse(wd)

    # ---- fingerprints ----
    def _find_first(self, wd: Path, candidates: list[Path]) -> Path | None:
        for c in candidates:
            if c.exists():
                return c
        return None

    def _fingerprint_candidates(self, wd: Path) -> list[Path]:
        custom = Path(self.cfg.ifp_output) if self.cfg.ifp_output else None
        fp_dir = custom.parent if custom else wd / "results" / "fingerprints"
        return [
            c for c in [
                custom,
                wd / "results" / "fingerprints" / "ifp.csv",
                fp_dir / "ifp_E.csv",
                fp_dir / "ifp_H.csv",
                fp_dir / "ifp_F.csv",
            ] if c
        ]

    def _sim_matrix_candidates(self, wd: Path) -> list[Path]:
        custom = Path(self.cfg.sim_matrix_output) if self.cfg.sim_matrix_output else None
        sim_dir = custom.parent if custom else wd
        return [
            c for c in [
                custom,
                wd / "sim_matrix.csv",
                wd / "results" / "sim_matrix.csv",
                sim_dir / "sim_matrix_E.csv",
                sim_dir / "sim_matrix_H.csv",
                sim_dir / "sim_matrix_F.csv",
            ] if c
        ]

    def _load_fingerprints(self, wd: Path) -> None:
        f = self._find_first(wd, self._fingerprint_candidates(wd))
        if not f:
            self.fp_path_label.setText("ifp.csv não encontrado")
            self.fp_table.clear(); self.fp_table.setRowCount(0); self.fp_table.setColumnCount(0)
            return
        self.fp_path_label.setText(str(f))
        max_rows = self.fp_rows.value()
        try:
            with f.open("r", encoding="utf-8", errors="replace") as fh:
                reader = csv.reader(fh)
                rows = []
                for i, row in enumerate(reader):
                    if i > max_rows:
                        break
                    rows.append(row)
        except Exception as e:
            QMessageBox.critical(self, "Erro ao ler CSV", str(e))
            return
        if not rows:
            return
        header, data = rows[0], rows[1:]
        self.fp_table.clear()
        self.fp_table.setColumnCount(len(header))
        self.fp_table.setHorizontalHeaderLabels(header)
        self.fp_table.setRowCount(len(data))
        for r, row in enumerate(data):
            for c, val in enumerate(row):
                self.fp_table.setItem(r, c, QTableWidgetItem(val))

    # ---- similarity matrix ----
    def _load_sim_matrix(self, wd: Path) -> None:
        if not HAS_MPL:
            return
        f = self._find_first(wd, self._sim_matrix_candidates(wd))
        if not f:
            self.sm_path_label.setText("sim_matrix.csv não encontrado")
            self.fig.clear(); self.canvas.draw()
            return
        self.sm_path_label.setText(str(f))
        try:
            with f.open("r", encoding="utf-8", errors="replace") as fh:
                reader = list(csv.reader(fh))
        except Exception as e:
            QMessageBox.critical(self, "Erro", str(e))
            return
        if not reader:
            return
        # Try to detect a header row / index column
        first = reader[0]
        has_header = any(not _is_float(x) for x in first[1:])
        labels: list[str] = []
        data: list[list[float]] = []
        start_row = 1 if has_header else 0
        for row in reader[start_row:]:
            if not row:
                continue
            if not _is_float(row[0]):
                labels.append(row[0])
                vals = row[1:]
            else:
                vals = row
            data.append([_safe_float(v) for v in vals])
        self.fig.clear()
        ax = self.fig.add_subplot(111)
        if not data:
            self.canvas.draw(); return
        im = ax.imshow(data, cmap="viridis", aspect="auto", vmin=0, vmax=1)
        ax.set_title("Tanimoto similarity")
        if labels and len(labels) <= 40:
            ax.set_xticks(range(len(labels))); ax.set_xticklabels(labels, rotation=90, fontsize=7)
            ax.set_yticks(range(len(labels))); ax.set_yticklabels(labels, fontsize=7)
        self.fig.colorbar(im, ax=ax)
        self.fig.tight_layout()
        self.canvas.draw()

    # ---- PSE files ----
    def _load_pse(self, wd: Path) -> None:
        self.pse_list.clear()
        candidates = []
        if self.cfg.pse_path:
            candidates.append(Path(self.cfg.pse_path))
        candidates += [wd / "results" / "pse", wd / "pse"]
        found_dir: Path | None = None
        for c in candidates:
            if c.exists() and c.is_dir():
                found_dir = c; break
        if not found_dir:
            return
        for f in sorted(found_dir.glob("*.pse")):
            it = QListWidgetItem(f.name)
            it.setData(Qt.ItemDataRole.UserRole, str(f))
            self.pse_list.addItem(it)

    def _open_pse(self) -> None:
        it = self.pse_list.currentItem()
        if not it:
            QMessageBox.information(self, "PSE", "Selecione um arquivo .pse na lista.")
            return
        path = it.data(Qt.ItemDataRole.UserRole)
        pymol = shutil.which("pymol") or shutil.which("pymol.exe")
        try:
            if pymol:
                subprocess.Popen([pymol, path])
            else:
                # Fall back to OS file association
                if sys.platform == "win32":
                    import os
                    os.startfile(path)  # type: ignore
                elif sys.platform == "darwin":
                    subprocess.Popen(["open", path])
                else:
                    subprocess.Popen(["xdg-open", path])
        except Exception as e:
            QMessageBox.critical(self, "Erro ao abrir", str(e))


    # ---- statistics (post-analysis via luna-env helper) ----
    def set_python(self, py_exe: str) -> None:
        self.py_exe = py_exe

    def compute_stats(self) -> None:
        wd = self._current_wd()
        if not wd:
            return
        if not self.py_exe:
            QMessageBox.warning(self, "luna-env",
                                "LUNA não detectado. Verifique a aba Setup.")
            return
        self.st_status.setText("Processando... (pode levar alguns minutos)")
        self.st_status.repaint()
        result = run_analysis(self.py_exe, str(wd))
        if "error" in result:
            self.st_status.setText("Erro")
            QMessageBox.critical(self, "Erro na análise", result["error"])
            return
        self._last_analysis = result
        n = result.get("entries", 0)
        self.st_status.setText(f"{n} entradas processadas")
        if HAS_MPL:
            self._render_stats_chart(result)

    def _render_stats_chart(self, result: dict) -> None:
        counts = result.get("interaction_counts", {})
        self.st_fig.clear()
        ax = self.st_fig.add_subplot(111)
        if not counts:
            ax.text(0.5, 0.5, "Sem dados de interação", ha="center", va="center")
            self.st_canvas.draw(); return
        items = sorted(counts.items(), key=lambda x: -x[1])
        labels = [k for k, _ in items]
        values = [v for _, v in items]
        ax.barh(labels, values, color="#2c5aa0")
        ax.invert_yaxis()
        ax.set_xlabel("Total (todas as entradas)")
        ax.set_title("Contagem por tipo de interação")
        self.st_fig.tight_layout()
        self.st_canvas.draw()

    # ---- Residue matrix heatmap ----
    def compute_residue_matrix(self) -> None:
        wd = self._current_wd()
        if not wd:
            return
        if not self.py_exe:
            QMessageBox.warning(self, "luna-env", "LUNA não detectado. Veja a aba Setup.")
            return
        self.hm_status.setText("Processando...")
        self.hm_status.repaint()
        r = run_residue_matrix(self.py_exe, str(wd))
        if "error" in r:
            self.hm_status.setText("Erro")
            QMessageBox.critical(self, "Erro na análise", r["error"])
            return
        self._residue_matrix = r
        types = r.get("interaction_types", [])
        self.hm_status.setText(f"{len(r.get('entries', []))} entradas · {len(types)} tipos")
        self.cb_itype.blockSignals(True)
        self.cb_itype.clear()
        self.cb_itype.addItems(types)
        self.cb_itype.blockSignals(False)
        if types:
            self._render_residue_heatmap()

    def _render_residue_heatmap(self) -> None:
        if not HAS_MPL or not self._residue_matrix:
            return
        it = self.cb_itype.currentText()
        if not it:
            return
        m = self._residue_matrix
        residues = m.get("residues", [])
        entries = m.get("entries", [])
        data = m.get("matrix", {}).get(it)
        if not data:
            return
        # Drop all-zero residue columns so the plot stays readable
        import numpy as _np
        arr = _np.array(data, dtype=float)
        keep = arr.sum(axis=0) > 0
        arr = arr[:, keep]
        residues = [r for r, k in zip(residues, keep) if k]
        self.hm_fig.clear()
        ax = self.hm_fig.add_subplot(111)
        if arr.size == 0:
            ax.text(0.5, 0.5, f"Sem ocorrências de '{it}'", ha="center", va="center")
        else:
            im = ax.imshow(arr, cmap="viridis", aspect="auto")
            ax.set_xticks(range(len(residues)))
            ax.set_xticklabels(residues, rotation=90, fontsize=7)
            ax.set_yticks(range(len(entries)))
            ax.set_yticklabels(entries, fontsize=7)
            ax.set_title(f"Resíduos × ligantes — {it}")
            self.hm_fig.colorbar(im, ax=ax)
        self.hm_fig.tight_layout()
        self.hm_canvas.draw()

    # ---- HTML report ----
    def export_report(self) -> None:
        wd = self._current_wd()
        if not wd:
            return
        from pathlib import Path as _P
        out, _ = QFileDialog.getSaveFileName(
            self, "Salvar relatório", str(wd / "luna_report.html"), "HTML (*.html)"
        )
        if not out:
            return
        # Make sure analyses are loaded
        if not self._last_analysis:
            self.compute_stats()
        # Save current heatmap & stats charts as PNGs to embed
        heatmap_png = wd / "_report_heatmap.png"
        inter_png = wd / "_report_interactions.png"
        try:
            if HAS_MPL and self.fig.axes:
                translate_figure(self.fig)
                self.fig.savefig(heatmap_png, dpi=120, bbox_inches="tight")
            if HAS_MPL and self.st_fig.axes:
                translate_figure(self.st_fig)
                self.st_fig.savefig(inter_png, dpi=120, bbox_inches="tight")
        except Exception:
            pass
        try:
            save_report(
                out, cfg=self.cfg, analysis=self._last_analysis,
                heatmap_png=heatmap_png if heatmap_png.exists() else None,
                interactions_png=inter_png if inter_png.exists() else None,
            )
        except Exception as e:
            QMessageBox.critical(self, "Erro", str(e))
            return
        QMessageBox.information(self, "Relatório salvo", out)


def _is_float(s: str) -> bool:
    try:
        float(s); return True
    except Exception:
        return False


def _safe_float(s: str) -> float:
    try:
        return float(s)
    except Exception:
        return 0.0
