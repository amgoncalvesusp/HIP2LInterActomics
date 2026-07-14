"""Run tab — execute LUNA via QProcess and stream the log."""
from __future__ import annotations

import os
import re
import shutil
import signal
import sys
from pathlib import Path

from PyQt6.QtCore import QProcess, QProcessEnvironment, QTimer, pyqtSignal
from PyQt6.QtGui import QFontDatabase
from PyQt6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QFormLayout, QLabel, QPushButton,
    QSpinBox, QCheckBox, QPlainTextEdit, QMessageBox, QGroupBox,
    QProgressBar,
)

from ..core.project import ProjectConfig
from ..core import env_manager as em, luna_runner, ligand_io, luna_api_runner


class RunTab(QWidget):
    finished_ok = pyqtSignal()

    def __init__(self, cfg: ProjectConfig) -> None:
        super().__init__()
        self.cfg = cfg
        self.py_exe: str = ""
        self.run_py: str = ""
        self.proc: QProcess | None = None
        self._process_group_owner: QProcess | None = None
        self._cancel_helper: QProcess | None = None
        self.collect_callback = None  # set by MainWindow
        self._last_progress_pct = 0

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
        max_nproc = 1 if sys.platform == "win32" else max(1, (os.cpu_count() or 1))
        self.sp_nproc.setRange(1, max_nproc)
        # Default to 1 on Windows: LUNA's multiprocessing mode fails on spawn
        # ("cannot pickle '_thread.lock' object"). Linux/Mac use fork and are fine.
        default_np = 1 if sys.platform == "win32" else max(1, (os.cpu_count() or 2) // 2)
        self.sp_nproc.setValue(default_np)
        if sys.platform == "win32":
            self.sp_nproc.setToolTip(
                "No Windows nativo, o LUNA roda com nproc=1. Para paralelismo real, "
                "use Linux ou WSL2, onde o multiprocessing usa fork."
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

        self.progress = QProgressBar()
        self.progress.setRange(0, 100)
        self.progress.setValue(0)
        self.progress.setTextVisible(True)
        self.progress.setFormat("Aguardando execução")
        self.progress.setToolTip("Mostra o progresso informado pelo LUNA durante a etapa atual.")
        layout.addWidget(self.progress)

        self.cmd_label = QLabel("")
        self.cmd_label.setWordWrap(True)
        self.cmd_label.setToolTip("Mostra o comando efetivo usado para executar o LUNA nesta rodada.")
        self.cmd_label.setFont(QFontDatabase.systemFont(QFontDatabase.SystemFont.FixedFont))
        self.cmd_label.setStyleSheet("color: #555; font-size: 10px;")
        layout.addWidget(self.cmd_label)

        self.log = QPlainTextEdit()
        self.log.setReadOnly(True)
        self.log.setToolTip("Exibe o log completo da execução, incluindo mensagens do LUNA e erros.")
        self.log.setFont(QFontDatabase.systemFont(QFontDatabase.SystemFont.FixedFont))
        self.log.setStyleSheet("font-size: 11px;")
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
            if self.collect_callback() is False:
                return

        self.cfg.nproc = luna_runner.safe_nproc(self.sp_nproc.value())
        self.cfg.overwrite = self.cb_overwrite.isChecked()

        errs = luna_runner.validate(self.cfg)
        if errs:
            QMessageBox.warning(self, "Configuração inválida", "\n".join(errs))
            return

        add_h_errs = luna_api_runner.validate_hydrogen_inputs(self.cfg)
        if add_h_errs:
            QMessageBox.warning(self, "Alerta de hidrogênios (Add_H)", "\n".join(add_h_errs))
            return

        # Write entries.txt to workdir
        wd = Path(self.cfg.workdir)
        wd.mkdir(parents=True, exist_ok=True)
        entries_file = wd / "entries.txt"
        ligand_io.write_entries_file(entries_file, self.cfg.selected_ligands)
        self.cfg.force_python_api = True

        # Pick CLI or Python-API runner based on advanced options and receptor handling.
        if luna_api_runner.should_use_api_runner(self.cfg):
            extra_errs = luna_api_runner.validate_entry_specs(
                self.cfg, self.cfg.selected_ligands
            )
            if extra_errs:
                QMessageBox.warning(self, "Configuração inválida", "\n".join(extra_errs))
                return
            try:
                cmd = luna_api_runner.build_api_command(
                    self.py_exe, self.cfg, self.cfg.selected_ligands
                )
            except Exception as e:
                QMessageBox.warning(self, "Configuração inválida", str(e))
                return
            runner_label = "Python API"
        else:
            cmd = luna_runner.build_command(
                self.py_exe, self.run_py, self.cfg, str(entries_file)
            )
            runner_label = "CLI"
        self.cmd_label.setText("$ " + " ".join(cmd))
        self.log.clear()
        self.log.appendPlainText(f"=== Iniciando LUNA ({runner_label}) ===\n")
        self.progress.setValue(0)
        self.progress.setFormat("0% - iniciando")
        self._last_progress_pct = 0

        self.proc = QProcess(self)
        proc = self.proc
        proc.setProcessChannelMode(QProcess.ProcessChannelMode.MergedChannels)
        program, arguments = cmd[0], cmd[1:]
        if sys.platform != "win32" and hasattr(proc, "setChildProcessModifier"):
            proc.setChildProcessModifier(os.setsid)
            self._process_group_owner = proc
        elif sys.platform != "win32" and (setsid := shutil.which("setsid")):
            program, arguments = setsid, cmd
            self._process_group_owner = proc
        # Force UTF-8 so LUNA's unicode status characters don't crash on cp1252 (Windows)
        proc_env = QProcessEnvironment()
        for key, value in em.python_process_env(self.py_exe).items():
            proc_env.insert(key, value)
        proc_env.insert("PYTHONIOENCODING", "utf-8")
        proc_env.insert("PYTHONUTF8", "1")
        proc.setProcessEnvironment(proc_env)
        proc.readyReadStandardOutput.connect(self._on_stdout)
        proc.finished.connect(self._on_finished)
        proc.errorOccurred.connect(self._on_process_error)
        self.btn_run.setEnabled(False)
        self.btn_cancel.setEnabled(True)
        proc.start(program, arguments)

    def is_running(self) -> bool:
        return self.proc is not None and self.proc.state() != QProcess.ProcessState.NotRunning

    def cancel(self) -> None:
        proc = self.proc
        if not proc or proc.state() == QProcess.ProcessState.NotRunning:
            return
        pid = int(proc.processId())
        self.btn_cancel.setEnabled(False)
        try:
            if pid > 0 and sys.platform == "win32":
                helper = QProcess(self)
                self._cancel_helper = helper
                helper.finished.connect(
                    lambda code, _status, original=proc, owner=helper: self._taskkill_finished(
                        original, owner, code
                    )
                )
                helper.errorOccurred.connect(
                    lambda _error, original=proc, owner=helper: self._taskkill_failed(
                        original, owner
                    )
                )
                helper.start("taskkill", ["/PID", str(pid), "/T", "/F"])
            elif pid > 0 and self._process_group_owner is proc:
                os.killpg(pid, signal.SIGTERM)
                QTimer.singleShot(2000, lambda p=pid, original=proc: self._kill_process_group(p, original))
            else:
                proc.terminate()
                QTimer.singleShot(2000, lambda original=proc: self._kill_if_running(original))
        except OSError:
            proc.kill()
        self.log.appendPlainText("\n[cancelado pelo usuário]")

    def _kill_process_group(self, pid: int, proc: QProcess) -> None:
        if self.proc is not proc or proc.state() == QProcess.ProcessState.NotRunning:
            return
        try:
            os.killpg(pid, signal.SIGKILL)
        except OSError:
            proc.kill()

    def _kill_if_running(self, proc: QProcess) -> None:
        if self.proc is proc and proc.state() != QProcess.ProcessState.NotRunning:
            proc.kill()

    def _taskkill_finished(self, proc: QProcess, helper: QProcess, code: int) -> None:
        if self._cancel_helper is helper:
            self._cancel_helper = None
        helper.deleteLater()
        if code != 0:
            self._kill_if_running(proc)

    def _taskkill_failed(self, proc: QProcess, helper: QProcess) -> None:
        self._taskkill_finished(proc, helper, 1)

    def _on_stdout(self) -> None:
        sender = self.sender()
        if isinstance(sender, QProcess) and sender is not self.proc:
            return
        if self.proc is None:
            return
        data = bytes(self.proc.readAllStandardOutput()).decode("utf-8", "replace")
        self.log.insertPlainText(data)
        self._update_progress_from_text(data)
        sb = self.log.verticalScrollBar()
        sb.setValue(sb.maximum())

    def _update_progress_from_text(self, text: str) -> None:
        clean = re.sub(r"\x1b\[[0-?]*[ -/]*[@-~]", "", text)
        api_progress = list(re.finditer(r"\[luna-api-progress\]\s*(\d{1,3})%\s*-\s*([^\r\n]+)", clean))
        if api_progress:
            match = api_progress[-1]
            pct = max(0, min(100, int(match.group(1))))
            stage = match.group(2).strip()
            self._set_progress(pct, stage)
            return

        pattern = re.compile(r"(\d{1,3})%\s+\[[^\]]*\]\s+(\d+)/(\d+).*?-\s*([^\r\n.]+)")
        matches = list(pattern.finditer(clean))
        if matches:
            match = matches[-1]
            pct = max(0, min(100, int(match.group(1))))
            done = match.group(2)
            total = match.group(3)
            stage = match.group(4).strip()
            self._set_progress(pct, f"{stage} ({done}/{total})")
            return

        stage_markers = [
            (r"novas entries carregadas|total de entries", 10, "entradas carregadas"),
            (r"iniciando proj\.run", 20, "calculando interacoes"),
            (r"resumo salvo|matriz de res", 72, "consolidando resultados"),
            (r"gerando fingerprints", 78, "gerando fingerprints"),
            (r"IFP .* salvo", 86, "fingerprints salvos"),
            (r"Similaridade", 91, "calculando similaridade"),
            (r"PSE salvos", 96, "gerando sessoes PyMOL"),
            (r"conclu", 100, "concluido"),
        ]
        for pattern_text, pct, stage in stage_markers:
            if re.search(pattern_text, clean, re.IGNORECASE):
                self._set_progress(pct, stage)

    def _set_progress(self, pct: int, stage: str) -> None:
        pct = max(0, min(100, int(pct)))
        if pct < self._last_progress_pct and pct != 0:
            pct = self._last_progress_pct
        self._last_progress_pct = pct
        self.progress.setValue(pct)
        self.progress.setFormat(f"{pct}% - {stage}")

    def _on_finished(self, code: int, _status) -> None:
        sender = self.sender()
        if isinstance(sender, QProcess) and sender is not self.proc:
            return
        proc = self.proc
        self.proc = None
        if self._process_group_owner is proc:
            self._process_group_owner = None
        if proc is not None:
            proc.deleteLater()
        self.log.appendPlainText(f"\n=== LUNA finalizou (exit code {code}) ===")
        self.btn_run.setEnabled(True)
        self.btn_cancel.setEnabled(False)
        if code == 0:
            self.progress.setValue(100)
            self.progress.setFormat("100% - concluído")
            self.finished_ok.emit()
            QMessageBox.information(self, "Concluído",
                                    f"Análise concluída.\nResultados em: {self.cfg.workdir}")
        else:
            self.progress.setFormat(f"Falhou (exit code {code})")

    def _on_process_error(self, error: QProcess.ProcessError) -> None:
        sender = self.sender()
        if isinstance(sender, QProcess) and sender is not self.proc:
            return
        if error != QProcess.ProcessError.FailedToStart:
            return
        proc = self.proc
        detail = proc.errorString() if proc is not None else "processo não iniciado"
        self.log.appendPlainText(f"\n[erro ao iniciar LUNA] {detail}")
        self.progress.setFormat("Falha ao iniciar")
        self.btn_run.setEnabled(True)
        self.btn_cancel.setEnabled(False)
        self.proc = None
        if self._process_group_owner is proc:
            self._process_group_owner = None
        if proc is not None:
            proc.deleteLater()
