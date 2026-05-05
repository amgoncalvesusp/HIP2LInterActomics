"""Wizard dialog - split complex files into protein + ligand inputs."""
from __future__ import annotations

from pathlib import Path

from PyQt6.QtWidgets import (
    QApplication,
    QDialog, QVBoxLayout, QHBoxLayout, QFormLayout, QLineEdit, QPushButton,
    QSpinBox, QFileDialog, QLabel, QPlainTextEdit, QMessageBox, QCheckBox,
    QGroupBox, QRadioButton, QButtonGroup, QWidget, QProgressBar,
)
from PyQt6.QtCore import Qt

from ..core import env_manager as em
from ..core.mol2_prep import split_complex_folder, detect_last_protein_atom


class DockingPrepDialog(QDialog):
    """Wizard to split combined complex files into separate receptor/ligand files."""

    def __init__(self, parent=None) -> None:
        super().__init__(parent)
        self.setWindowTitle("Preparar arquivos de complexos")
        self.resize(760, 560)

        self.result_protein_dir: str | None = None
        self.result_ligand_dir: str | None = None
        self.water_molecules_detected: int = 0
        self._detected_last_pa: int | None = None

        layout = QVBoxLayout(self)

        layout.addWidget(QLabel(
            "Selecione uma pasta com arquivos .mol2, .pdb ou .ent "
            "(proteína + ligante no mesmo arquivo).\n"
            "Serão geradas subpastas separadas para proteínas e ligantes compatíveis "
            "com LUNA: MOL2 quando a origem é MOL2, SDF quando a origem é PDB/ENT, "
            "com águas preservadas junto à proteína."
        ))

        form = QFormLayout()

        # --- Pasta de origem ---
        self.src_edit = QLineEdit()
        btn_src = QPushButton("...")
        btn_src.setFixedWidth(30)
        btn_src.clicked.connect(self._pick_src)
        row = QHBoxLayout()
        row.addWidget(self.src_edit)
        row.addWidget(btn_src)
        form.addRow("Pasta de origem:", self._wrap(row))

        # --- Pasta de saída ---
        self.out_edit = QLineEdit()
        self.out_edit.setPlaceholderText("(padrão: mesma pasta da origem)")
        btn_out = QPushButton("...")
        btn_out.setFixedWidth(30)
        btn_out.clicked.connect(self._pick_out)
        row = QHBoxLayout()
        row.addWidget(self.out_edit)
        row.addWidget(btn_out)
        form.addRow("Pasta de saída:", self._wrap(row))

        layout.addLayout(form)

        # --- Seção: último átomo da proteína ---
        last_box = QGroupBox("Último átomo da proteína")
        last_layout = QVBoxLayout(last_box)

        # Radio buttons
        self.rb_auto = QRadioButton("Detectar automaticamente")
        self.rb_manual = QRadioButton("Informar manualmente")
        self.rb_auto.setChecked(True)
        btn_grp = QButtonGroup(self)
        btn_grp.addButton(self.rb_auto)
        btn_grp.addButton(self.rb_manual)
        self.rb_auto.toggled.connect(self._on_mode_changed)

        radio_row = QHBoxLayout()
        radio_row.addWidget(self.rb_auto)
        radio_row.addWidget(self.rb_manual)
        radio_row.addStretch()
        last_layout.addLayout(radio_row)

        # Auto-detect panel
        self.auto_panel = QWidget()
        auto_layout = QHBoxLayout(self.auto_panel)
        auto_layout.setContentsMargins(0, 0, 0, 0)
        self.btn_detect = QPushButton("Detectar agora")
        self.btn_detect.clicked.connect(self._detect)
        self.detect_result_label = QLabel("(selecione a pasta de origem primeiro)")
        self.detect_result_label.setStyleSheet("color: #555;")
        auto_layout.addWidget(self.btn_detect)
        auto_layout.addWidget(self.detect_result_label, 1)
        last_layout.addWidget(self.auto_panel)

        # Manual panel
        self.manual_panel = QWidget()
        manual_layout = QHBoxLayout(self.manual_panel)
        manual_layout.setContentsMargins(0, 0, 0, 0)
        self.sp_last = QSpinBox()
        self.sp_last.setRange(1, 1_000_000)
        self.sp_last.setValue(4068)
        self.sp_last.setToolTip(
            "Número do último átomo da proteína no MOL2.\n"
            "Átomos com índice > este valor são tratados como ligante."
        )
        manual_layout.addWidget(QLabel("Último átomo da proteína:"))
        manual_layout.addWidget(self.sp_last)
        manual_layout.addStretch()
        last_layout.addWidget(self.manual_panel)
        self.manual_panel.setVisible(False)

        layout.addWidget(last_box)

        # --- Opções extras ---
        self.cb_open_after = QCheckBox(
            "Usar as pastas geradas como entradas do projeto ao fechar"
        )
        self.cb_open_after.setChecked(True)
        layout.addWidget(self.cb_open_after)

        # --- Botões de ação ---
        btn_row = QHBoxLayout()
        self.btn_run = QPushButton("Executar preparação")
        self.btn_run.setStyleSheet("padding: 6px 14px;")
        self.btn_run.clicked.connect(self._run)
        self.btn_close = QPushButton("Fechar")
        self.btn_close.clicked.connect(self.accept)
        btn_row.addWidget(self.btn_run)
        btn_row.addStretch()
        btn_row.addWidget(self.btn_close)
        layout.addLayout(btn_row)

        progress_row = QHBoxLayout()
        self.progress = QProgressBar()
        self.progress.setRange(0, 1)
        self.progress.setValue(0)
        self.progress.setFormat("%v / %m")
        progress_row.addWidget(self.progress, 1)
        self.progress_status = QLabel("Aguardando execução.")
        self.progress_status.setProperty("muted", True)
        progress_row.addWidget(self.progress_status, 1)
        layout.addLayout(progress_row)

        # --- Log ---
        self.log = QPlainTextEdit()
        self.log.setReadOnly(True)
        self.log.setStyleSheet("font-family: Consolas, monospace; font-size: 11px;")
        layout.addWidget(self.log, 1)

    # ---- helpers ----

    def _wrap(self, sub_layout) -> QWidget:
        w = QWidget(); w.setLayout(sub_layout); return w

    def _on_mode_changed(self) -> None:
        auto = self.rb_auto.isChecked()
        self.auto_panel.setVisible(auto)
        self.manual_panel.setVisible(not auto)

    # ---- pickers ----

    def _pick_src(self) -> None:
        d = QFileDialog.getExistingDirectory(self, "Pasta com arquivos de complexos")
        if d:
            self.src_edit.setText(d)
            # Trigger auto-detect whenever a new folder is chosen
            if self.rb_auto.isChecked():
                self._detect()

    def _pick_out(self) -> None:
        d = QFileDialog.getExistingDirectory(self, "Pasta de saída")
        if d:
            self.out_edit.setText(d)

    # ---- auto-detect ----

    def _detect(self) -> None:
        src = self.src_edit.text().strip()
        if not src or not Path(src).is_dir():
            self.detect_result_label.setText("Selecione a pasta de origem primeiro.")
            return

        # Pick the first .mol2 in the folder as representative sample.
        # PDB/ENT inputs do not need the last-protein-atom boundary.
        mol2_files = sorted(Path(src).glob("*.mol2"))
        if not mol2_files:
            pdb_files = [
                candidate for candidate in sorted(Path(src).iterdir())
                if candidate.is_file() and candidate.suffix.lower() in {".pdb", ".ent"}
            ]
            if pdb_files:
                self._detected_last_pa = None
                self.detect_result_label.setText(
                    f"{len(pdb_files)} arquivo(s) PDB/ENT detectado(s). "
                    "Não é necessário informar o último átomo da proteína."
                )
                self.detect_result_label.setStyleSheet("color: #1a7a1a;")
                return
            self.detect_result_label.setText("Nenhum arquivo .mol2, .pdb ou .ent encontrado.")
            return

        sample = mol2_files[0]
        self.detect_result_label.setText(f"Analisando {sample.name}…")
        self.detect_result_label.repaint()

        try:
            r = detect_last_protein_atom(sample)
        except Exception as e:
            self.detect_result_label.setText(f"Erro: {e}")
            self._detected_last_pa = None
            return

        self._detected_last_pa = r.last_pa
        lig_str = ", ".join(r.ligand_names[:4])
        if len(r.ligand_names) > 4:
            lig_str += f" … (+{len(r.ligand_names)-4})"

        msg = (
            f"Último átomo da proteína: {r.last_pa}  "
            f"(método: {r.method})  |  "
            f"Proteína: {r.n_protein_atoms} átomos  |  "
            f"Ligante: {r.n_ligand_atoms} átomos"
        )
        if lig_str:
            msg += f"\nLigantes detectados: {lig_str}"
        self.detect_result_label.setText(msg)
        self.detect_result_label.setStyleSheet("color: #1a7a1a;")

        self.log.appendPlainText(f"[auto-detect] {sample.name}")
        self.log.appendPlainText(f"  last_pa={r.last_pa}  método={r.method}")
        self.log.appendPlainText(f"  proteína: {r.n_protein_atoms} átomos")
        self.log.appendPlainText(f"  ligante(s): {r.n_ligand_atoms} átomos — {lig_str or 'nenhum'}")

    # ---- main action ----

    def _run(self) -> None:
        src = self.src_edit.text().strip()
        if not src or not Path(src).is_dir():
            QMessageBox.warning(self, "Pasta inválida", "Selecione uma pasta de origem válida.")
            return

        mol2_files = sorted(Path(src).glob("*.mol2"))
        pdb_files = [
            candidate for candidate in sorted(Path(src).iterdir())
            if candidate.is_file() and candidate.suffix.lower() in {".pdb", ".ent"}
        ]

        if mol2_files and pdb_files:
            QMessageBox.warning(
                self,
                "Pasta mista",
                "Use uma pasta com apenas arquivos MOL2 ou apenas arquivos PDB/ENT por preparação.",
            )
            return

        # Resolve last_pa only for MOL2 inputs.
        if self.rb_auto.isChecked():
            if mol2_files and self._detected_last_pa is None:
                # Try to detect right now if user forgot to click "Detectar"
                self._detect()
            if mol2_files and self._detected_last_pa is None:
                QMessageBox.warning(
                    self, "Detecção falhou",
                    "Não foi possível detectar o último átomo automaticamente.\n"
                    "Troque para 'Informar manualmente' e insira o valor."
                )
                return
            last_pa = self._detected_last_pa
        else:
            last_pa = self.sp_last.value() if mol2_files else None

        out = self.out_edit.text().strip() or None

        self.log.appendPlainText(
            f"\nIniciando preparação — last_pa={last_pa}, src={src}"
        )
        total_files = len(mol2_files) + len(pdb_files)
        self.progress.setRange(0, max(1, total_files))
        self.progress.setValue(0)
        self.progress_status.setText(f"0 / {total_files} arquivos processados")
        self.btn_run.setEnabled(False)
        self.btn_close.setEnabled(False)

        def _on_progress(processed: int, total: int, filename: str, ok: bool, error_message: str) -> None:
            self.progress.setRange(0, max(1, total))
            self.progress.setValue(min(processed, max(1, total)))
            if filename:
                status = "ok" if ok else f"erro: {error_message}"
                self.progress_status.setText(f"{processed} / {total} — {filename} ({status})")
            else:
                self.progress_status.setText(f"{processed} / {total} arquivos processados")
            QApplication.processEvents()

        try:
            chemistry_python = self._chemistry_python() if pdb_files else None
            if pdb_files and chemistry_python:
                self.log.appendPlainText(
                    f"Conversao quimica PDB -> SDF usando luna-env: {chemistry_python}"
                )
            elif pdb_files:
                self.log.appendPlainText(
                    "Aviso: luna-env nao foi localizado; tentando converter PDB -> SDF apenas com o Python atual."
                )
            r = split_complex_folder(
                src,
                last_pa,
                out,
                progress_cb=_on_progress,
                chemistry_python=chemistry_python,
            )
        except Exception as e:
            self.btn_run.setEnabled(True)
            self.btn_close.setEnabled(True)
            QMessageBox.critical(self, "Erro", str(e))
            return
        finally:
            self.btn_run.setEnabled(True)
            self.btn_close.setEnabled(True)

        self.log.appendPlainText(
            f"Arquivos lidos: {r.files_processed} | "
            f"Proteínas: {r.proteins_written} | Ligantes: {r.ligands_written} | "
            f"Águas detectadas: {r.water_molecules_detected}"
        )
        self.log.appendPlainText(f"Proteínas → {r.protein_dir}")
        self.log.appendPlainText(f"Ligantes  → {r.ligand_dir}")
        for err in r.errors:
            self.log.appendPlainText(f"[erro] {err}")

        if r.ligands_written == 0:
            QMessageBox.warning(
                self, "Nenhum ligante gerado",
                f"Nenhum ligante foi extraído com last_pa={last_pa}.\n\n"
                "Se o valor foi detectado automaticamente, tente trocar para "
                "'Informar manualmente' e ajustar o número."
            )
            return

        self.result_protein_dir = r.protein_dir
        self.result_ligand_dir = r.ligand_dir
        self.water_molecules_detected = r.water_molecules_detected
        self.progress.setValue(self.progress.maximum())
        self.progress_status.setText(
            f"Concluído — {r.files_processed} arquivos processados | "
            f"{r.water_molecules_detected} águas detectadas"
        )
        self.log.appendPlainText("\nConcluido com sucesso.")

    def _chemistry_python(self) -> str | None:
        try:
            conda = em.find_conda()
            if not conda:
                return None
            py = em.env_python(conda)
            if py and py.exists():
                return str(py)
        except Exception:
            return None
        return None
