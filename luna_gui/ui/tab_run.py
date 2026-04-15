"""Run tab — execute LUNA via QProcess and stream the log."""
from __future__ import annotations

import os
import sys
from pathlib import Path

from PyQt6.QtCore import QProcess, QProcessEnvironment, pyqtSignal
from PyQt6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QFormLayout, QLabel, QPushButton,
    QSpinBox, QCheckBox, QPlainTextEdit, QMessageBox, QGroupBox,
)

from ..core.project import ProjectConfig
from ..core import luna_runner, ligand_io, luna_api_runner


class RunTab(QWidget):
    finished_ok = pyqtSignal()

    def __init__(self, cfg: ProjectConfig) -> None:
        super().__init__()
        self.cfg = cfg
        self.py_exe: str = ""
        self.run_py: str = ""
        self.proc: QProcess | None = None
        self.collect_callback = None  # set by MainWindow

        layout = QVBoxLayout(self)

        intro = QLabel(
            "Esta aba executa o projeto com as opções já escolhidas. O número de núcleos controla "
            "o paralelismo, e o log abaixo mostra exatamente o que o LUNA está fazendo."
        )
        intro.setWordWrap(True)
        intro.setProperty("muted", True)
        layout.addWidget(intro)

        opt_box = QGroupBox("Opções de execução")
        form = QFormLayout(opt_box)
        self.sp_nproc = QSpinBox()
        self.sp_nproc.setRange(1, max(1, (os.cpu_count() or 1)))
        # Default to 1 on Windows: LUNA's multiprocessing mode fails on spawn
        # ("cannot pickle '_thread.lock' object"). Linux/Mac use fork and are fine.
        default_np = 1 if sys.platform == "win32" else max(1, (os.cpu_count() or 2) // 2)
        self.sp_nproc.setValue(default_np)
        if sys.platform == "win32":
            self.sp_nproc.setToolTip(
                "No Windows, recomenda-se 1 — LUNA tem um bug de pickling "
                "com multiprocessing (spawn) que quebra com nproc > 1."
            )
        else:
            self.sp_nproc.setToolTip("Define quantos núcleos de CPU o LUNA pode usar na execução.")
        form.addRow("Núcleos (--nproc):", self.sp_nproc)
        self.cb_overwrite = QCheckBox("Sobrescrever projeto existente (--overwrite)")
        self.cb_overwrite.setToolTip("Permite reutilizar o mesmo diretório de projeto, substituindo saídas antigas.")
        form.addRow(self.cb_overwrite)
        layout.addWidget(opt_box)

        btn_row = QHBoxLayout()
        self.btn_run = QPushButton("▶ Executar LUNA")
        self.btn_run.setStyleSheet("font-size: 14px; padding: 8px;")
        self.btn_run.setToolTip("Valida as entradas, grava o entries.txt e inicia a execução do LUNA.")
        self.btn_run.clicked.connect(self.run)
        self.btn_cancel = QPushButton("■ Cancelar")
        self.btn_cancel.setEnabled(False)
        self.btn_cancel.setToolTip("Interrompe a execução em andamento.")
        self.btn_cancel.clicked.connect(self.cancel)
        btn_row.addWidget(self.btn_run)
        btn_row.addWidget(self.btn_cancel)
        btn_row.addStretch()
        layout.addLayout(btn_row)

        self.cmd_label = QLabel("")
        self.cmd_label.setWordWrap(True)
        self.cmd_label.setToolTip("Mostra o comando efetivo usado para executar o LUNA nesta rodada.")
        self.cmd_label.setStyleSheet("color: #555; font-family: Consolas, monospace; font-size: 10px;")
        layout.addWidget(self.cmd_label)

        self.log = QPlainTextEdit()
        self.log.setReadOnly(True)
        self.log.setToolTip("Exibe o log completo da execução, incluindo mensagens do LUNA e erros.")
        self.log.setStyleSheet("font-family: Consolas, monospace; font-size: 11px;")
        layout.addWidget(self.log, 1)

    # ---- public API ----
    def set_luna(self, py_exe: str, run_py: str) -> None:
        self.py_exe = py_exe
        self.run_py = run_py

    def run(self) -> None:
        if not self.py_exe or not self.run_py:
            QMessageBox.warning(self, "LUNA não pronto",
                                "Vá à aba 'Setup' e instale/verifique o LUNA primeiro.")
            return

        if self.collect_callback:
            self.collect_callback()

        self.cfg.nproc = self.sp_nproc.value()
        self.cfg.overwrite = self.cb_overwrite.isChecked()

        errs = luna_runner.validate(self.cfg)
        if errs:
            QMessageBox.warning(self, "Configuração inválida", "\n".join(errs))
            return

        # Write entries.txt to workdir
        wd = Path(self.cfg.workdir)
        wd.mkdir(parents=True, exist_ok=True)
        entries_file = wd / "entries.txt"
        ligand_io.write_entries_file(entries_file, self.cfg.selected_ligands)

        # Pick CLI or Python-API runner based on advanced-option flags
        if self.cfg.uses_python_api():
            cmd = luna_api_runner.build_api_command(
                self.py_exe, self.cfg, self.cfg.selected_ligands
            )
            runner_label = "Python API"
        else:
            cmd = luna_runner.build_command(
                self.py_exe, self.run_py, self.cfg, str(entries_file)
            )
            runner_label = "CLI"
        self.cmd_label.setText("$ " + " ".join(cmd))
        self.log.clear()
        self.log.appendPlainText(f"=== Iniciando LUNA ({runner_label}) ===\n")

        self.proc = QProcess(self)
        self.proc.setProcessChannelMode(QProcess.ProcessChannelMode.MergedChannels)
        # Force UTF-8 so LUNA's unicode status characters don't crash on cp1252 (Windows)
        env = QProcessEnvironment.systemEnvironment()
        env.insert("PYTHONIOENCODING", "utf-8")
        env.insert("PYTHONUTF8", "1")
        self.proc.setProcessEnvironment(env)
        self.proc.readyReadStandardOutput.connect(self._on_stdout)
        self.proc.finished.connect(self._on_finished)
        self.btn_run.setEnabled(False)
        self.btn_cancel.setEnabled(True)
        self.proc.start(cmd[0], cmd[1:])

    def cancel(self) -> None:
        if self.proc and self.proc.state() != QProcess.ProcessState.NotRunning:
            self.proc.kill()
            self.log.appendPlainText("\n[cancelado pelo usuário]")

    def _on_stdout(self) -> None:
        data = bytes(self.proc.readAllStandardOutput()).decode("utf-8", "replace")
        self.log.insertPlainText(data)
        sb = self.log.verticalScrollBar()
        sb.setValue(sb.maximum())

    def _on_finished(self, code: int, _status) -> None:
        self.log.appendPlainText(f"\n=== LUNA finalizou (exit code {code}) ===")
        self.btn_run.setEnabled(True)
        self.btn_cancel.setEnabled(False)
        if code == 0:
            self.finished_ok.emit()
            QMessageBox.information(self, "Concluído",
                                    f"Análise concluída.\nResultados em: {self.cfg.workdir}")
