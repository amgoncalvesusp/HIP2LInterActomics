"""Results tab — preview fingerprints, inspect charts, export graphics, and cluster ligands."""
from __future__ import annotations

import csv
import math
import shutil
import subprocess
import sys
from pathlib import Path

from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QPushButton, QLineEdit, QLabel,
    QTableWidget, QTableWidgetItem, QFileDialog, QListWidget, QListWidgetItem,
    QMessageBox, QTabWidget, QSpinBox, QComboBox, QHeaderView,
)

from ..core.project import ProjectConfig
from ..core.analysis_runtime import run_analysis, run_residue_matrix
from ..core.pymol_launcher import launch_pse_session
from ..core.report_export import save_pdf_report, save_report
from ..i18n import translate_figure
from ..core.results_analysis import (
    count_tanimoto_similarity,
    load_analysis_summary,
    load_residue_matrix_artifact,
    load_ifp_sparse_matrix,
    load_similarity_matrix,
)

try:
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
    HAS_MPL = True
except Exception:
    HAS_MPL = False

try:
    from scipy.cluster.hierarchy import dendrogram
    from ..core.results_analysis import (
        cluster_rows,
        cluster_similarity_matrix,
        export_cluster_assignments,
    )
    HAS_CLUSTERING = True
except Exception:
    HAS_CLUSTERING = False


class SortableTableWidgetItem(QTableWidgetItem):
    def __lt__(self, other):  # type: ignore[override]
        if isinstance(other, QTableWidgetItem):
            left = self.data(Qt.ItemDataRole.UserRole)
            right = other.data(Qt.ItemDataRole.UserRole)
            if left is not None and right is not None:
                try:
                    return left < right
                except Exception:
                    pass
        return super().__lt__(other)


def _sortable_item(text: str, sort_value=None) -> QTableWidgetItem:
    item = SortableTableWidgetItem(text)
    if sort_value is not None:
        item.setData(Qt.ItemDataRole.UserRole, sort_value)
    return item


def _popen_detached(args: list[str]) -> subprocess.Popen:
    kwargs = {}
    if sys.platform == "win32":
        flags = 0
        flags |= getattr(subprocess, "CREATE_NEW_PROCESS_GROUP", 0)
        flags |= getattr(subprocess, "DETACHED_PROCESS", 0)
        if flags:
            kwargs["creationflags"] = flags
    else:
        kwargs["start_new_session"] = True
    return subprocess.Popen(args, **kwargs)


class ResultsTab(QWidget):
    def __init__(self, cfg: ProjectConfig) -> None:
        super().__init__()
        self.cfg = cfg
        self.py_exe: str = ""
        self._last_analysis: dict = {}
        self._sim_labels: list[str] = []
        self._sim_matrix = None
        self._cluster_result = None
        self._residue_matrix: dict = {}

        layout = QVBoxLayout(self)

        intro = QLabel(
            "Use esta aba para revisar o que o LUNA gerou. Você pode carregar um workdir já existente, "
            "visualizar tabelas e gráficos, exportar figuras e agrupar ligantes por similaridade."
        )
        intro.setWordWrap(True)
        intro.setProperty("muted", True)
        layout.addWidget(intro)

        wd_row = QHBoxLayout()
        wd_row.addWidget(QLabel("Workdir:"))
        self.wd_edit = QLineEdit()
        self.wd_edit.setPlaceholderText("(usa o workdir do projeto atual)")
        self.wd_edit.setToolTip("Pasta do projeto cujos resultados serão carregados nesta aba.")
        btn_wd = QPushButton("Procurar...")
        btn_wd.setToolTip("Escolhe manualmente um workdir já processado pelo LUNA.")
        btn_wd.clicked.connect(self._pick_wd)
        btn_load = QPushButton("Carregar resultados")
        btn_load.setToolTip("Lê fingerprints, matriz de similaridade e sessões PyMOL do workdir selecionado.")
        btn_load.clicked.connect(self.load_all)
        btn_export_chart = QPushButton("Exportar gráfico atual...")
        btn_export_chart.setToolTip("Salva a figura da aba atual em PNG, SVG ou PDF.")
        btn_export_chart.clicked.connect(self.export_current_chart)
        btn_report = QPushButton("Exportar relatório HTML")
        btn_report.setToolTip("Gera um relatório HTML com os principais gráficos e resumos disponíveis.")
        btn_report.clicked.connect(self.export_report)
        btn_report_pdf = QPushButton("Gerar relatório PDF")
        btn_report_pdf.setToolTip("Gera um PDF com parâmetros, explicações e interpretação básica das análises.")
        btn_report_pdf.clicked.connect(self.export_pdf_report)
        wd_row.addWidget(self.wd_edit, 1)
        wd_row.addWidget(btn_wd)
        wd_row.addWidget(btn_load)
        wd_row.addWidget(btn_export_chart)
        wd_row.addWidget(btn_report)
        wd_row.addWidget(btn_report_pdf)
        layout.addLayout(wd_row)

        self.inner = QTabWidget()
        self.inner.setUsesScrollButtons(True)
        layout.addWidget(self.inner, 1)

        self.fp_tab = QWidget()
        fp_layout = QVBoxLayout(self.fp_tab)
        fp_help = QLabel(
            "Pré-visualização do arquivo `ifp.csv`. Cada linha representa um ligante e cada coluna descreve partes do fingerprint gerado."
        )
        fp_help.setWordWrap(True)
        fp_help.setProperty("muted", True)
        fp_layout.addWidget(fp_help)
        fp_ctrl = QHBoxLayout()
        fp_ctrl.addWidget(QLabel("Linhas a exibir:"))
        self.fp_rows = QSpinBox()
        self.fp_rows.setRange(10, 10000)
        self.fp_rows.setValue(200)
        self.fp_rows.setToolTip("Limita quantas linhas do CSV serão mostradas na tabela para facilitar a inspeção.")
        self.fp_rows.valueChanged.connect(self._reload_fingerprints_preview)
        fp_ctrl.addWidget(self.fp_rows)
        fp_ctrl.addWidget(QLabel("Tipo:"))
        self.cb_fp_preview_type = QComboBox()
        self.cb_fp_preview_type.setToolTip("Escolhe qual fingerprint calculado será exibido na tabela.")
        self.cb_fp_preview_type.currentIndexChanged.connect(self._reload_fingerprints_preview)
        fp_ctrl.addWidget(self.cb_fp_preview_type)
        self.fp_path_label = QLabel("—")
        self.fp_path_label.setProperty("muted", True)
        fp_ctrl.addWidget(self.fp_path_label, 1)
        fp_layout.addLayout(fp_ctrl)
        self.fp_table = QTableWidget()
        self.fp_table.horizontalHeader().setStretchLastSection(True)
        self.fp_table.setSortingEnabled(True)
        fp_layout.addWidget(self.fp_table, 1)
        self.inner.addTab(self.fp_tab, "Fingerprints")

        self.sim_tab = QWidget()
        sim_layout = QVBoxLayout(self.sim_tab)
        sim_help = QLabel(
            "Mapa de calor da similaridade entre ligantes. Valores mais altos indicam fingerprints mais parecidos."
        )
        sim_help.setWordWrap(True)
        sim_help.setProperty("muted", True)
        sim_layout.addWidget(sim_help)
        sim_ctrl = QHBoxLayout()
        sim_ctrl.addWidget(QLabel("Tipo:"))
        self.cb_sim_type = QComboBox()
        self.cb_sim_type.setToolTip("Escolhe qual matriz de similaridade carregada será exibida.")
        self.cb_sim_type.currentIndexChanged.connect(self._reload_similarity_preview)
        sim_ctrl.addWidget(self.cb_sim_type)
        sim_ctrl.addStretch()
        sim_layout.addLayout(sim_ctrl)
        self.sm_path_label = QLabel("—")
        self.sm_path_label.setProperty("muted", True)
        sim_layout.addWidget(self.sm_path_label)
        if HAS_MPL:
            self.fig = Figure(figsize=(6.4, 5.2))
            self.canvas = FigureCanvas(self.fig)
            sim_layout.addWidget(self.canvas, 1)
        else:
            sim_layout.addWidget(QLabel("matplotlib não está instalado."))
        self.inner.addTab(self.sim_tab, "Matriz de similaridade")

        self.stats_tab = QWidget()
        st_layout = QVBoxLayout(self.stats_tab)
        stats_help = QLabel(
            "Resume quantas interações de cada tipo aparecem no conjunto analisado. É uma visão global do perfil químico observado."
        )
        stats_help.setWordWrap(True)
        stats_help.setProperty("muted", True)
        st_layout.addWidget(stats_help)
        st_ctrl = QHBoxLayout()
        btn_st = QPushButton("Calcular estatísticas (usa luna-env)")
        btn_st.setToolTip("Varre os resultados do projeto e gera um resumo por tipo de interação.")
        btn_st.clicked.connect(self.compute_stats)
        self.st_status = QLabel("—")
        self.st_status.setProperty("muted", True)
        st_ctrl.addWidget(btn_st)
        st_ctrl.addWidget(self.st_status, 1)
        st_layout.addLayout(st_ctrl)
        if HAS_MPL:
            self.st_fig = Figure(figsize=(6.2, 4.2))
            self.st_canvas = FigureCanvas(self.st_fig)
            st_layout.addWidget(self.st_canvas, 1)
        else:
            st_layout.addWidget(QLabel("matplotlib não está instalado."))
        self.inner.addTab(self.stats_tab, "Estatísticas")

        self.residue_tab = QWidget()
        hm_layout = QVBoxLayout(self.residue_tab)
        residue_help = QLabel(
            "Mostra, para um tipo de interação escolhido, quais resíduos da proteína aparecem associados aos ligantes."
        )
        residue_help.setWordWrap(True)
        residue_help.setProperty("muted", True)
        hm_layout.addWidget(residue_help)
        hm_ctrl = QHBoxLayout()
        btn_hm = QPushButton("Calcular heatmap (usa luna-env)")
        btn_hm.setToolTip("Calcula a matriz resíduo × ligante a partir dos resultados do projeto.")
        btn_hm.clicked.connect(self.compute_residue_matrix)
        hm_ctrl.addWidget(btn_hm)
        hm_ctrl.addWidget(QLabel("Tipo:"))
        self.cb_itype = QComboBox()
        self.cb_itype.setToolTip("Escolhe qual classe de interação será exibida no heatmap.")
        self.cb_itype.currentIndexChanged.connect(self._render_residue_heatmap)
        hm_ctrl.addWidget(self.cb_itype, 1)
        self.hm_status = QLabel("—")
        self.hm_status.setProperty("muted", True)
        hm_ctrl.addWidget(self.hm_status)
        hm_layout.addLayout(hm_ctrl)
        if HAS_MPL:
            self.hm_fig = Figure(figsize=(8, 5))
            self.hm_canvas = FigureCanvas(self.hm_fig)
            hm_layout.addWidget(self.hm_canvas, 1)
        else:
            hm_layout.addWidget(QLabel("matplotlib não está instalado."))
        self.inner.addTab(self.residue_tab, "Heatmap por tipo")

        self.cluster_tab = QWidget()
        cluster_layout = QVBoxLayout(self.cluster_tab)
        cluster_help = QLabel(
            "Agrupa ligantes com base na matriz de similaridade. O dendrograma mostra a relação hierárquica, e a tabela lista o cluster atribuído a cada ligante."
        )
        cluster_help.setWordWrap(True)
        cluster_help.setProperty("muted", True)
        cluster_layout.addWidget(cluster_help)
        cluster_ctrl = QHBoxLayout()
        cluster_ctrl.addWidget(QLabel("Método:"))
        self.cluster_method = QComboBox()
        self.cluster_method.addItems(["average", "complete", "single"])
        self.cluster_method.setToolTip("Escolhe o critério de ligação hierárquica usado para formar os clusters.")
        cluster_ctrl.addWidget(self.cluster_method)
        cluster_ctrl.addWidget(QLabel("Clusters:"))
        self.cluster_count = QSpinBox()
        self.cluster_count.setRange(2, 50)
        self.cluster_count.setValue(4)
        self.cluster_count.setToolTip("Define em quantos grupos finais os ligantes serão separados.")
        cluster_ctrl.addWidget(self.cluster_count)
        btn_cluster = QPushButton("Atualizar clusters")
        btn_cluster.setToolTip("Recalcula o dendrograma e as atribuições de cluster com as opções escolhidas.")
        btn_cluster.clicked.connect(self.refresh_clusters)
        btn_clusters_csv = QPushButton("Exportar clusters CSV")
        btn_clusters_csv.setToolTip("Salva uma tabela CSV com ligante, cluster e ordem hierárquica.")
        btn_clusters_csv.clicked.connect(self.export_clusters_csv)
        cluster_ctrl.addWidget(btn_cluster)
        cluster_ctrl.addWidget(btn_clusters_csv)
        self.cluster_status = QLabel("—")
        self.cluster_status.setProperty("muted", True)
        cluster_ctrl.addWidget(self.cluster_status, 1)
        cluster_layout.addLayout(cluster_ctrl)
        if HAS_MPL and HAS_CLUSTERING:
            self.cluster_fig = Figure(figsize=(8, 6.2))
            self.cluster_canvas = FigureCanvas(self.cluster_fig)
            cluster_layout.addWidget(self.cluster_canvas, 1)
        else:
            msg = "Instale scipy para habilitar clustering." if HAS_MPL else "matplotlib não está instalado."
            cluster_layout.addWidget(QLabel(msg))
        self.cluster_table = QTableWidget()
        self.cluster_table.setColumnCount(3)
        self.cluster_table.setHorizontalHeaderLabels(["Ligante", "Cluster", "Ordem"])
        self.cluster_table.setToolTip("Lista final dos agrupamentos calculados para cada ligante.")
        self.cluster_table.setSortingEnabled(True)
        self.cluster_table.horizontalHeader().setSectionResizeMode(0, QHeaderView.ResizeMode.Stretch)
        self.cluster_table.horizontalHeader().setSectionResizeMode(1, QHeaderView.ResizeMode.ResizeToContents)
        self.cluster_table.horizontalHeader().setSectionResizeMode(2, QHeaderView.ResizeMode.ResizeToContents)
        cluster_layout.addWidget(self.cluster_table, 1)
        self.inner.addTab(self.cluster_tab, "Clusters")

        self.pse_tab = QWidget()
        pse_layout = QVBoxLayout(self.pse_tab)
        pse_help = QLabel(
            "Lista as sessões PyMOL exportadas pelo LUNA. Abra um arquivo para inspecionar visualmente interações e poses."
        )
        pse_help.setWordWrap(True)
        pse_help.setProperty("muted", True)
        pse_layout.addWidget(pse_help)
        pse_ctrl = QHBoxLayout()
        btn_open = QPushButton("Abrir no PyMOL")
        btn_open.setToolTip("Abre a sessão selecionada no PyMOL ou no programa associado ao sistema.")
        btn_open.clicked.connect(self._open_pse)
        pse_ctrl.addWidget(btn_open)
        pse_ctrl.addStretch()
        pse_layout.addLayout(pse_ctrl)
        self.pse_list = QListWidget()
        self.pse_list.setToolTip("Arquivos .pse encontrados no workdir atual.")
        self.pse_list.itemDoubleClicked.connect(lambda _: self._open_pse())
        pse_layout.addWidget(self.pse_list, 1)
        self.inner.addTab(self.pse_tab, "Sessões PyMOL")

    def _pick_wd(self) -> None:
        d = QFileDialog.getExistingDirectory(self, "Workdir do projeto LUNA")
        if d:
            self.wd_edit.setText(d)

    def _current_wd(self) -> Path | None:
        wd = self.wd_edit.text().strip() or self.cfg.workdir
        if not wd:
            QMessageBox.warning(self, "Workdir", "Defina um workdir primeiro.")
            return None
        path = Path(wd)
        if not path.exists():
            QMessageBox.warning(self, "Workdir", f"Diretório não existe:\n{path}")
            return None
        return path

    def load_all(self) -> None:
        wd = self._current_wd()
        if not wd:
            return
        self._populate_fp_preview_sources(wd)
        self._populate_similarity_sources(wd)
        self._load_fingerprints(wd)
        self._load_sim_matrix(wd)
        self._load_pse(wd)

    def _reload_fingerprints_preview(self) -> None:
        wd = self._current_wd()
        if wd:
            self._populate_fp_preview_sources(wd)
            self._load_fingerprints(wd)

    def _reload_similarity_preview(self) -> None:
        wd = self._current_wd()
        if wd:
            self._populate_similarity_sources(wd)
            self._load_sim_matrix(wd)

    def _find_first(self, candidates: list[Path]) -> Path | None:
        for candidate in candidates:
            if candidate.exists():
                return candidate
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

    def _normalize_ifp_source_key(self, path: Path) -> tuple[str, str]:
        stem = path.stem.lower()
        if stem.endswith("_e"):
            return "EIFP", path.name
        if stem.endswith("_f"):
            return "FIFP", path.name
        if stem.endswith("_h"):
            return "HIFP", path.name
        return "DEFAULT", path.name

    def _populate_fp_preview_sources(self, wd: Path) -> None:
        current = self.cb_fp_preview_type.currentData() if self.cb_fp_preview_type.count() else None
        candidates = [candidate for candidate in self._fingerprint_candidates(wd) if candidate.exists()]
        typed = []
        default = []
        for candidate in candidates:
            key, label = self._normalize_ifp_source_key(candidate)
            bucket = typed if key in {"EIFP", "FIFP", "HIFP"} else default
            bucket.append((key, label, str(candidate)))
        source_rows = typed or default
        self.cb_fp_preview_type.blockSignals(True)
        self.cb_fp_preview_type.clear()
        seen_keys: set[str] = set()
        for key, label, path_str in source_rows:
            if key in seen_keys:
                continue
            seen_keys.add(key)
            self.cb_fp_preview_type.addItem(label, path_str)
        if self.cb_fp_preview_type.count():
            idx = self.cb_fp_preview_type.findData(current)
            self.cb_fp_preview_type.setCurrentIndex(idx if idx >= 0 else 0)
        self.cb_fp_preview_type.blockSignals(False)

    def _selected_fingerprint_path(self, wd: Path) -> Path | None:
        selected = self.cb_fp_preview_type.currentData() if self.cb_fp_preview_type.count() else None
        if selected:
            candidate = Path(str(selected))
            if candidate.exists():
                return candidate
        return self._find_first(self._fingerprint_candidates(wd))

    def _sim_matrix_candidates(self, wd: Path) -> list[Path]:
        custom = Path(self.cfg.sim_matrix_output) if self.cfg.sim_matrix_output else None
        sim_dir = custom.parent if custom else wd
        custom_square = custom.with_name(f"{custom.stem}_square.csv") if custom else None
        return [
            c for c in [
                custom_square,
                custom,
                wd / "sim_matrix.csv",
                wd / "results" / "sim_matrix.csv",
                wd / "sim_matrix_square.csv",
                wd / "results" / "sim_matrix_square.csv",
                sim_dir / "sim_matrix_E_square.csv",
                sim_dir / "sim_matrix_H_square.csv",
                sim_dir / "sim_matrix_F_square.csv",
                sim_dir / "sim_matrix_E.csv",
                sim_dir / "sim_matrix_H.csv",
                sim_dir / "sim_matrix_F.csv",
            ] if c
        ]

    def _normalize_similarity_source_key(self, path: Path) -> tuple[str, str]:
        stem = path.stem.lower()
        if "sim_matrix_e" in stem:
            return "EIFP", "E (Extended)"
        if "sim_matrix_f" in stem:
            return "FIFP", "F (Functional)"
        if "sim_matrix_h" in stem:
            return "HIFP", "H (Hybrid)"
        return "DEFAULT", path.name

    def _populate_similarity_sources(self, wd: Path) -> None:
        current = self.cb_sim_type.currentData() if self.cb_sim_type.count() else None
        candidates = [candidate for candidate in self._sim_matrix_candidates(wd) if candidate.exists()]
        typed = []
        default = []
        for candidate in candidates:
            key, label = self._normalize_similarity_source_key(candidate)
            row = (key, label, str(candidate))
            if key in {"EIFP", "FIFP", "HIFP"}:
                typed.append(row)
            else:
                default.append(row)
        source_rows = typed or default
        self.cb_sim_type.blockSignals(True)
        self.cb_sim_type.clear()
        seen_keys: set[str] = set()
        for _key, label, path_str in source_rows:
            if _key in seen_keys:
                continue
            seen_keys.add(_key)
            self.cb_sim_type.addItem(label, path_str)
        if self.cb_sim_type.count():
            idx = self.cb_sim_type.findData(current)
            self.cb_sim_type.setCurrentIndex(idx if idx >= 0 else 0)
        self.cb_sim_type.blockSignals(False)

    def _selected_similarity_path(self, wd: Path) -> Path | None:
        selected = self.cb_sim_type.currentData() if self.cb_sim_type.count() else None
        if selected:
            candidate = Path(str(selected))
            if candidate.exists():
                return candidate
        return self._find_first(self._sim_matrix_candidates(wd))

    def _load_fingerprints(self, wd: Path) -> None:
        file_path = self._selected_fingerprint_path(wd)
        if not file_path:
            self.fp_path_label.setText("Nenhum fingerprint calculado foi encontrado")
            self.fp_table.clear()
            self.fp_table.setRowCount(0)
            self.fp_table.setColumnCount(0)
            return

        self.fp_path_label.setText(str(file_path))
        try:
            with file_path.open("r", encoding="utf-8", errors="replace") as fh:
                reader = csv.reader(fh)
                rows = []
                limit = self.fp_rows.value() + 1
                for i, row in enumerate(reader):
                    if i >= limit:
                        break
                    rows.append(row)
        except Exception as exc:
            QMessageBox.critical(self, "Erro ao ler CSV", str(exc))
            return

        if not rows:
            return
        header, data = rows[0], rows[1:]
        self.fp_table.clear()
        self.fp_table.setColumnCount(len(header))
        self.fp_table.setHorizontalHeaderLabels(header)
        self.fp_table.setRowCount(len(data))
        self.fp_table.setSortingEnabled(False)
        for r, row in enumerate(data):
            for c, value in enumerate(row):
                try:
                    sort_value = float(value)
                except Exception:
                    sort_value = value
                self.fp_table.setItem(r, c, _sortable_item(value, sort_value))
        self.fp_table.setSortingEnabled(True)

    def _load_sim_matrix(self, wd: Path) -> None:
        file_path = self._selected_similarity_path(wd)
        if not file_path:
            fp_path = self._selected_fingerprint_path(wd)
            if fp_path and fp_path.exists():
                try:
                    labels, _feature_ids, fp_matrix = load_ifp_sparse_matrix(fp_path)
                    self._sim_labels = labels
                    self._sim_matrix = count_tanimoto_similarity(fp_matrix)
                    self.sm_path_label.setText(f"Reconstruida a partir de {fp_path.name}")
                except Exception as exc:
                    self._sim_labels = []
                    self._sim_matrix = None
                    self.sm_path_label.setText("Nenhuma matriz de similaridade calculada foi encontrada")
                    if HAS_MPL:
                        self.fig.clear()
                        self.canvas.draw()
                    self._clear_clusters("Matriz de similaridade não encontrada.")
                    QMessageBox.critical(self, "Erro ao reconstruir matriz", str(exc))
                    return
            else:
                self._sim_labels = []
                self._sim_matrix = None
                self.sm_path_label.setText("Nenhuma matriz de similaridade calculada foi encontrada")
                if HAS_MPL:
                    self.fig.clear()
                    self.canvas.draw()
                self._clear_clusters("Matriz de similaridade não encontrada.")
                return
        else:
            self.sm_path_label.setText(str(file_path))
            try:
                labels, matrix = load_similarity_matrix(file_path)
            except Exception as exc:
                self._sim_labels = []
                self._sim_matrix = None
                self._clear_clusters("Erro ao carregar a matriz.")
                QMessageBox.critical(self, "Erro ao ler matriz", str(exc))
                return

            self._sim_labels = labels
            self._sim_matrix = matrix

        if HAS_MPL:
            label_count = max(1, len(self._sim_labels))
            width_in = max(8.4, 3.2 + ((0.5 / 2.54) * label_count))
            height_in = max(7.0, 3.0 + ((0.5 / 2.54) * label_count))
            self.fig.set_dpi(120)
            self.fig.set_size_inches(width_in, height_in, forward=True)
            self.canvas.setMinimumSize(int(width_in * self.fig.dpi), int(height_in * self.fig.dpi))
            self.fig.clear()
            ax = self.fig.add_subplot(111)
            im = ax.imshow(self._sim_matrix, cmap="viridis", aspect="auto", vmin=0, vmax=1)
            ax.set_title("Similaridade de Tanimoto")
            _apply_tick_labels(ax, self._sim_labels, axis="both", ligand_axis=True)
            self.fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
            self.fig.tight_layout()
            self.canvas.draw()

        self.refresh_clusters()

    def _load_pse(self, wd: Path) -> None:
        self.pse_list.clear()
        candidates = []
        if self.cfg.pse_path:
            candidates.append(Path(self.cfg.pse_path))
        candidates.extend([wd / "results" / "pse", wd / "pse"])
        found_dir = None
        for candidate in candidates:
            if candidate.exists() and candidate.is_dir():
                found_dir = candidate
                break
        if not found_dir:
            return
        pse_files = list(found_dir.glob("*.pse"))
        filtered_root = wd / "results" / "pse_filtered"
        if filtered_root.exists():
            pse_files.extend(filtered_root.glob("*/*.pse"))
        for file_path in sorted(pse_files):
            label = file_path.name
            try:
                if filtered_root.exists() and filtered_root.resolve() in file_path.resolve().parents:
                    label = f"{file_path.parent.name}/{file_path.name}"
            except Exception:
                pass
            item = QListWidgetItem(label)
            item.setData(Qt.ItemDataRole.UserRole, str(file_path))
            self.pse_list.addItem(item)

    def _open_pse(self) -> None:
        item = self.pse_list.currentItem()
        if not item:
            QMessageBox.information(self, "PSE", "Selecione um arquivo .pse na lista.")
            return
        path = item.data(Qt.ItemDataRole.UserRole)
        try:
            launch_pse_session(path, self.py_exe)
        except Exception as exc:
            QMessageBox.critical(self, "Erro ao abrir PyMOL", str(exc))

    def set_python(self, py_exe: str) -> None:
        self.py_exe = py_exe

    def compute_stats(self) -> None:
        wd = self._current_wd()
        if not wd:
            return
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
        self.st_status.setText(f"{result.get('entries', 0)} entradas processadas")
        if HAS_MPL:
            self._render_stats_chart(result)

    def _render_stats_chart(self, result: dict) -> None:
        counts = result.get("interaction_counts", {})
        self.st_fig.clear()
        ax = self.st_fig.add_subplot(111)
        if not counts:
            ax.text(0.5, 0.5, "Sem dados de interação", ha="center", va="center")
            self.st_canvas.draw()
            return
        items = sorted(counts.items(), key=lambda x: -x[1])
        labels = [k for k, _ in items]
        values = [v for _, v in items]
        ax.barh(labels, values, color="#c8693a")
        ax.invert_yaxis()
        ax.set_xlabel("Total (todas as entradas)")
        ax.set_title("Contagem por tipo de interação")
        self.st_fig.tight_layout()
        self.st_canvas.draw()

    def compute_residue_matrix(self) -> None:
        wd = self._current_wd()
        if not wd:
            return
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

    def _render_residue_heatmap(self) -> None:
        if not HAS_MPL or not self._residue_matrix:
            return
        interaction_type = self.cb_itype.currentText()
        if not interaction_type:
            return

        residues = self._residue_matrix.get("residues", [])
        entries = self._residue_matrix.get("entries", [])
        data = self._residue_matrix.get("matrix", {}).get(interaction_type)
        if not data:
            return

        import numpy as _np

        arr = _np.array(data, dtype=float)
        keep = arr.sum(axis=0) > 0
        arr = arr[:, keep]
        residues = [residue for residue, flag in zip(residues, keep) if flag]

        self.hm_fig.clear()
        ax = self.hm_fig.add_subplot(111)
        if arr.size == 0:
            ax.text(0.5, 0.5, f"Sem ocorrências de '{interaction_type}'", ha="center", va="center")
        else:
            im = ax.imshow(arr, cmap="viridis", aspect="auto")
            _apply_tick_labels(ax, residues, axis="x")
            _apply_tick_labels(ax, entries, axis="y", max_labels=28, rotation=0, ligand_axis=True)
            ax.set_title(f"Resíduos × ligantes — {interaction_type}")
            self.hm_fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        self.hm_fig.tight_layout()
        self.hm_canvas.draw()

    def refresh_clusters(self) -> None:
        if not HAS_CLUSTERING:
            self._clear_clusters("Clustering indisponível: scipy não encontrado.")
            return
        if self._sim_matrix is None or not self._sim_labels:
            self._clear_clusters("Gere ou carregue uma matriz de similaridade para clusterizar.")
            return
        try:
            result = cluster_similarity_matrix(
                self._sim_labels,
                self._sim_matrix,
                method=self.cluster_method.currentText(),
                n_clusters=self.cluster_count.value(),
            )
        except Exception as exc:
            self._clear_clusters("Falha ao clusterizar.")
            QMessageBox.critical(self, "Erro de clustering", str(exc))
            return

        self._cluster_result = result
        self.cluster_status.setText(f"{result.n_clusters} clusters · método {result.method}")
        self._render_cluster_chart(result)
        self._populate_cluster_table(result)

    def _render_cluster_chart(self, result) -> None:
        if not (HAS_MPL and HAS_CLUSTERING):
            return
        label_count = max(1, len(result.labels))
        width_in = max(9.0, 4.4 + ((0.5 / 2.54) * label_count))
        height_in = max(8.0, 4.2 + ((0.5 / 2.54) * label_count))
        self.cluster_fig.set_dpi(120)
        self.cluster_fig.set_size_inches(width_in, height_in, forward=True)
        self.cluster_canvas.setMinimumSize(int(width_in * self.cluster_fig.dpi), int(height_in * self.cluster_fig.dpi))
        self.cluster_fig.clear()
        grid = self.cluster_fig.add_gridspec(2, 1, height_ratios=[1.3, 2.4], hspace=0.35)

        ax_tree = self.cluster_fig.add_subplot(grid[0])
        labels = result.labels if len(result.labels) <= 40 else None
        dendrogram(
            result.linkage_matrix,
            labels=labels,
            ax=ax_tree,
            leaf_rotation=90,
            leaf_font_size=7,
            color_threshold=None,
        )
        ax_tree.set_title("Clustering hierárquico")
        ax_tree.set_ylabel("Distância")
        if labels is None:
            ax_tree.set_xticks([])
            ax_tree.set_xlabel(f"{len(result.labels)} ligantes em ordem hierárquica")

        ax_heat = self.cluster_fig.add_subplot(grid[1])
        im = ax_heat.imshow(result.ordered_matrix, cmap="magma", aspect="auto", vmin=0, vmax=1)
        ax_heat.set_title("Matriz reordenada por cluster")
        _apply_tick_labels(ax_heat, result.ordered_labels, axis="both", ligand_axis=True)
        self.cluster_fig.colorbar(im, ax=ax_heat, fraction=0.046, pad=0.04)
        self.cluster_fig.tight_layout()
        self.cluster_canvas.draw()

    def _populate_cluster_table(self, result) -> None:
        rows = cluster_rows(result)
        self.cluster_table.setSortingEnabled(False)
        self.cluster_table.setRowCount(len(rows))
        for row_index, (label, cluster_id, leaf_order) in enumerate(rows):
            self.cluster_table.setItem(row_index, 0, _sortable_item(label, label))
            self.cluster_table.setItem(row_index, 1, _sortable_item(str(cluster_id), int(cluster_id)))
            self.cluster_table.setItem(row_index, 2, _sortable_item(str(leaf_order), int(leaf_order)))
        self.cluster_table.resizeRowsToContents()
        self.cluster_table.setSortingEnabled(True)

    def _clear_clusters(self, message: str) -> None:
        self._cluster_result = None
        self.cluster_status.setText(message)
        self.cluster_table.setRowCount(0)
        if HAS_MPL and HAS_CLUSTERING:
            self.cluster_fig.clear()
            self.cluster_canvas.draw()

    def export_current_chart(self) -> None:
        figure, default_name = self._current_figure()
        if figure is None or not figure.axes:
            QMessageBox.information(self, "Sem gráfico", "A aba atual não possui um gráfico exportável.")
            return

        wd = self._current_wd()
        default_dir = wd if wd else Path.cwd()
        out, _ = QFileDialog.getSaveFileName(
            self,
            "Exportar gráfico",
            str(default_dir / f"{default_name}.png"),
            "PNG (*.png);;SVG (*.svg);;PDF (*.pdf)",
        )
        if not out:
            return
        try:
            translate_figure(figure)
            figure.savefig(out, dpi=180, bbox_inches="tight")
        except Exception as exc:
            QMessageBox.critical(self, "Erro ao exportar gráfico", str(exc))

    def _current_figure(self):
        current = self.inner.currentWidget()
        if HAS_MPL and current is self.sim_tab:
            return self.fig, "similaridade"
        if HAS_MPL and current is self.stats_tab:
            return self.st_fig, "estatisticas_interacoes"
        if HAS_MPL and current is self.residue_tab:
            interaction_type = self.cb_itype.currentText().strip().replace(" ", "_").lower() or "heatmap"
            return self.hm_fig, f"heatmap_{interaction_type}"
        if HAS_MPL and HAS_CLUSTERING and current is self.cluster_tab:
            return self.cluster_fig, "clusters"
        return None, ""

    def export_clusters_csv(self) -> None:
        if not self._cluster_result:
            QMessageBox.information(self, "Clusters", "Nenhum cluster foi calculado ainda.")
            return
        wd = self._current_wd()
        default_dir = wd if wd else Path.cwd()
        out, _ = QFileDialog.getSaveFileName(
            self,
            "Exportar clusters",
            str(default_dir / "clusters.csv"),
            "CSV (*.csv)",
        )
        if not out:
            return
        try:
            export_cluster_assignments(out, self._cluster_result)
        except Exception as exc:
            QMessageBox.critical(self, "Erro ao exportar clusters", str(exc))

    def export_report(self) -> None:
        wd = self._current_wd()
        if not wd:
            return
        out, _ = QFileDialog.getSaveFileName(
            self,
            "Salvar relatório",
            str(wd / "luna_report.html"),
            "HTML (*.html)",
        )
        if not out:
            return

        if not self._last_analysis:
            self.compute_stats()
        if not self._last_analysis:
            return

        heatmap_png = wd / "_report_heatmap.png"
        inter_png = wd / "_report_interactions.png"
        cluster_png = wd / "_report_clusters.png"
        try:
            if HAS_MPL and self.fig.axes:
                translate_figure(self.fig)
                self.fig.savefig(heatmap_png, dpi=140, bbox_inches="tight")
            if HAS_MPL and self.st_fig.axes:
                translate_figure(self.st_fig)
                self.st_fig.savefig(inter_png, dpi=140, bbox_inches="tight")
            if HAS_MPL and HAS_CLUSTERING and self._cluster_result and self.cluster_fig.axes:
                translate_figure(self.cluster_fig)
                self.cluster_fig.savefig(cluster_png, dpi=140, bbox_inches="tight")
        except Exception:
            pass

        cluster_items = None
        if self._cluster_result:
            cluster_items = [
                (label, cluster_id)
                for label, cluster_id, _ in cluster_rows(self._cluster_result)
            ]

        try:
            save_report(
                out,
                cfg=self.cfg,
                analysis=self._last_analysis,
                heatmap_png=heatmap_png if heatmap_png.exists() else None,
                interactions_png=inter_png if inter_png.exists() else None,
                cluster_png=cluster_png if cluster_png.exists() else None,
                clusters=cluster_items,
            )
        except Exception as exc:
            QMessageBox.critical(self, "Erro", str(exc))
            return
        QMessageBox.information(self, "Relatório salvo", out)

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

        heatmap_png = wd / "_report_pdf_similarity.png"
        inter_png = wd / "_report_pdf_interactions.png"
        cluster_png = wd / "_report_pdf_clusters.png"
        try:
            if HAS_MPL and hasattr(self, "fig") and self.fig.axes:
                translate_figure(self.fig)
                self.fig.savefig(heatmap_png, dpi=600, bbox_inches="tight", pad_inches=0.18)
            if HAS_MPL and hasattr(self, "st_fig") and self.st_fig.axes:
                translate_figure(self.st_fig)
                self.st_fig.savefig(inter_png, dpi=600, bbox_inches="tight", pad_inches=0.18)
            if HAS_MPL and HAS_CLUSTERING and self._cluster_result and self.cluster_fig.axes:
                translate_figure(self.cluster_fig)
                self.cluster_fig.savefig(cluster_png, dpi=600, bbox_inches="tight", pad_inches=0.18)
        except Exception:
            pass

        cluster_items = None
        if self._cluster_result:
            cluster_items = [
                (label, cluster_id)
                for label, cluster_id, _ in cluster_rows(self._cluster_result)
            ]

        try:
            save_pdf_report(
                out,
                cfg=self.cfg,
                analysis=self._last_analysis,
                heatmap_png=heatmap_png if heatmap_png.exists() else None,
                interactions_png=inter_png if inter_png.exists() else None,
                cluster_png=cluster_png if cluster_png.exists() else None,
                clusters=cluster_items,
                fp_dashboards=getattr(self, "_fp_dashboards", {}),
            )
        except Exception as exc:
            QMessageBox.critical(self, "Erro ao gerar PDF", str(exc))
            return
        QMessageBox.information(self, "Relatório PDF salvo", out)


def _apply_tick_labels(
    ax,
    labels: list[str],
    axis: str = "both",
    max_labels: int = 32,
    rotation: int = 90,
    ligand_axis: bool = False,
    ligand_limit: int = 220,
) -> None:
    labels = [str(label) for label in labels]
    if ligand_axis and len(labels) > ligand_limit:
        if axis in ("x", "both"):
            ax.set_xticks([])
            ax.set_xlabel("Todos os ligantes")
        if axis in ("y", "both"):
            ax.set_yticks([])
            ax.set_ylabel("Todos os ligantes")
        return
    effective_max = len(labels) if ligand_axis and labels else max_labels
    fontsize = 7 * 0.85 if ligand_axis else 7
    positions = _tick_positions(len(labels), max_labels=effective_max)
    if axis in ("x", "both"):
        ax.set_xticks(positions)
        ax.set_xticklabels([labels[i] for i in positions], rotation=rotation, fontsize=fontsize)
    if axis in ("y", "both"):
        ax.set_yticks(positions)
        ax.set_yticklabels([labels[i] for i in positions], fontsize=fontsize)


def _tick_positions(total: int, max_labels: int = 32) -> list[int]:
    if total <= 0:
        return []
    step = max(1, math.ceil(total / max_labels))
    return list(range(0, total, step))
