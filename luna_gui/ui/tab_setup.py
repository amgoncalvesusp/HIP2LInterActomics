"""Setup tab — detect/install conda + luna-env, manage LUNA installation."""
from __future__ import annotations

import webbrowser
from pathlib import Path

from PyQt6.QtCore import QProcess, pyqtSignal
from PyQt6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QPushButton, QPlainTextEdit,
    QMessageBox,
)

from ..core import env_manager as em


class SetupTab(QWidget):
    """Detects conda, creates luna-env, installs LUNA."""

    luna_ready = pyqtSignal(str, str)  # (env_python_path, run_py_path)

    def __init__(self) -> None:
        super().__init__()
        self.proc: QProcess | None = None
        self._cmd_queue: list[list[str]] = []

        layout = QVBoxLayout(self)

        self.status_label = QLabel("Verificando ambiente...")
        self.status_label.setStyleSheet("font-size: 14px; padding: 8px;")
        layout.addWidget(self.status_label)

        help_label = QLabel(
            "Use esta aba para preparar o ambiente da GUI. Primeiro confira se o Conda foi "
            "encontrado; se não existir, baixe o Miniconda. Depois instale o LUNA no ambiente "
            "separado 'luna-env'."
        )
        help_label.setWordWrap(True)
        help_label.setProperty("muted", True)
        layout.addWidget(help_label)

        btn_row = QHBoxLayout()
        self.btn_refresh = QPushButton("Verificar novamente")
        self.btn_refresh.setToolTip("Refaz a detecção do Conda, do ambiente 'luna-env' e da instalação do LUNA.")
        self.btn_refresh.clicked.connect(self.detect)
        self.btn_install_miniconda = QPushButton("Baixar Miniconda")
        self.btn_install_miniconda.setToolTip("Abre a página oficial do Miniconda para instalar o Conda no sistema.")
        self.btn_install_miniconda.clicked.connect(
            lambda: webbrowser.open(em.miniconda_download_url())
        )
        self.btn_install_luna = QPushButton("Instalar LUNA (cria luna-env)")
        self.btn_install_luna.setToolTip(
            "Cria ou atualiza o ambiente 'luna-env' e instala nele as dependências usadas pela GUI."
        )
        self.btn_install_luna.clicked.connect(self.install_luna)
        btn_row.addWidget(self.btn_refresh)
        btn_row.addWidget(self.btn_install_miniconda)
        btn_row.addWidget(self.btn_install_luna)
        btn_row.addStretch()
        layout.addLayout(btn_row)

        self.log = QPlainTextEdit()
        self.log.setReadOnly(True)
        self.log.setToolTip("Mostra o passo a passo da verificação e da instalação do ambiente.")
        self.log.setStyleSheet("font-family: Consolas, monospace; font-size: 11px;")
        layout.addWidget(self.log, 1)

        self.detect()

    # ---- detection ----
    def detect(self) -> None:
        self.log.appendPlainText("=== Verificação de ambiente ===")
        conda = em.find_conda()
        if not conda:
            self.status_label.setText("❌ Conda não encontrado. Instale o Miniconda primeiro.")
            self.btn_install_luna.setEnabled(False)
            self.log.appendPlainText("Conda: NÃO encontrado.")
            return
        self.log.appendPlainText(f"Conda: {conda}")
        self.conda = conda

        if not em.env_exists(conda):
            self.status_label.setText("⚠️  Conda OK. Ambiente 'luna-env' não existe — clique em 'Instalar LUNA'.")
            self.btn_install_luna.setEnabled(True)
            return

        py = em.env_python(conda)
        if not py:
            self.status_label.setText("⚠️  luna-env existe mas python não foi localizado.")
            return
        self.log.appendPlainText(f"Python do env: {py}")

        if not em.luna_installed(py):
            self.status_label.setText("⚠️  luna-env existe mas LUNA não está instalado.")
            self.btn_install_luna.setEnabled(True)
            return

        run_py = em.luna_run_py_path(py)
        if not run_py:
            self.status_label.setText("⚠️  LUNA importado mas run.py não foi localizado.")
            return

        self.status_label.setText(f"✅ LUNA pronto. run.py: {run_py}")
        self.log.appendPlainText(f"run.py: {run_py}")
        self.luna_ready.emit(str(py), str(run_py))

    # ---- install ----
    def install_luna(self) -> None:
        conda = em.find_conda()
        if not conda:
            QMessageBox.warning(self, "Conda ausente", "Instale o Miniconda primeiro.")
            return
        self._cmd_queue = em.install_commands(conda)
        self.btn_install_luna.setEnabled(False)
        self.log.appendPlainText("\n=== Iniciando instalação do LUNA ===")
        self.log.appendPlainText("Isso pode levar vários minutos.\n")
        self._run_next()

    def _run_next(self) -> None:
        if not self._cmd_queue:
            self.log.appendPlainText("\n=== Instalação concluída ===")
            self.btn_install_luna.setEnabled(True)
            self.detect()
            return
        cmd = self._cmd_queue.pop(0)
        self.log.appendPlainText(f"\n$ {' '.join(cmd)}\n")
        self.proc = QProcess(self)
        self.proc.setProcessChannelMode(QProcess.ProcessChannelMode.MergedChannels)
        self.proc.readyReadStandardOutput.connect(self._on_stdout)
        self.proc.finished.connect(self._on_finished)
        self.proc.start(cmd[0], cmd[1:])

    def _on_stdout(self) -> None:
        data = bytes(self.proc.readAllStandardOutput()).decode("utf-8", "replace")
        self.log.insertPlainText(data)
        sb = self.log.verticalScrollBar()
        sb.setValue(sb.maximum())

    def _on_finished(self, code: int, _status) -> None:
        self.log.appendPlainText(f"[exit code: {code}]")
        if code != 0:
            self.log.appendPlainText("⚠️  Comando falhou — abortando instalação.")
            self.btn_install_luna.setEnabled(True)
            self._cmd_queue.clear()
            return
        self._run_next()
