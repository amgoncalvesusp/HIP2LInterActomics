"""Setup tab — detect/install conda + luna-env, manage LUNA installation."""
from __future__ import annotations

import webbrowser
from pathlib import Path

from PyQt6.QtCore import Qt, QSize, QProcess, QProcessEnvironment, pyqtSignal
from PyQt6.QtGui import QFontDatabase, QPixmap
from PyQt6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QGridLayout, QLabel, QPushButton, QPlainTextEdit,
    QMessageBox, QGroupBox,
)

from ..core import env_manager as em
from .async_task import AsyncTask
from .info import InfoButton


class SetupTab(QWidget):
    """Detects conda, creates luna-env, installs LUNA."""

    luna_ready = pyqtSignal(str, str)  # (env_python_path, run_py_path)

    def __init__(self) -> None:
        super().__init__()
        self.proc: QProcess | None = None
        self._cmd_queue: list[list[str]] = []
        self._detect_task: AsyncTask | None = None

        layout = QVBoxLayout(self)
        layout.setContentsMargins(12, 12, 12, 12)
        layout.setSpacing(10)

        self.status_label = QLabel("Verificando ambiente...")
        self.status_label.setStyleSheet("font-size: 14px; padding: 8px;")

        about_box = QGroupBox("O que é o HIP²LInterActomics")
        about_layout = QVBoxLayout(about_box)
        about_layout.setSpacing(14)
        icon_label = QLabel()
        icon_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        icon_path = Path(__file__).resolve().parents[1] / "assets" / "hip2l_interactomics_icon.png"
        if icon_path.exists():
            pixmap = QPixmap(str(icon_path))
            icon_label.setPixmap(
                pixmap.scaled(
                    QSize(190, 190),
                    Qt.AspectRatioMode.KeepAspectRatio,
                    Qt.TransformationMode.SmoothTransformation,
                )
            )
        icon_label.setFixedHeight(200)
        about_layout.addWidget(icon_label)

        help_label = QLabel(
            "HIP²LInterActomics é uma interface gráfica para estabelecer um fluxo de trabalho "
            "eficiente na análise de interações intermoleculares protein-ligand e protein-protein.\n\n"
            "O software ajuda a avaliar virtual screening, inspecionar poses de docking ou frames "
            "de dinâmica molecular, calcular fingerprints interpretáveis, organizar dados para "
            "modelos de machine learning e extrair pharmacophoric features de um alvo específico.\n\n"
            "Autores do software: Daniel Andrés Grajales Ruiz e Adriano Marques Gonçalves."
        )
        help_label.setWordWrap(True)
        help_label.setProperty("muted", True)
        help_label.setStyleSheet("font-size: 15px; line-height: 1.65; padding: 2px 8px 6px 8px;")
        about_layout.addWidget(help_label)
        workflow = QLabel(
            "Entradas estruturais -> Pré-processamento de complexos -> LUNA -> "
            "Fingerprints e matriz resíduo x interação -> Estatísticas, mapas de calor, modelos e sessões PyMOL"
        )
        workflow.setWordWrap(True)
        workflow.setAlignment(Qt.AlignmentFlag.AlignCenter)
        workflow.setStyleSheet(
            "padding: 10px; border: 1px solid #9db2c8; border-radius: 8px; "
            "background: #f4f8fc; color: #0b2b45; font-weight: 600;"
        )
        about_layout.addWidget(workflow)
        layout.addWidget(about_box, 2)

        env_box = QGroupBox("Configuração do ambiente e instalação de pacotes")
        env_layout = QVBoxLayout(env_box)
        env_layout.addWidget(self.status_label)

        btn_row = QGridLayout()
        self.btn_refresh = QPushButton("Verificar novamente")
        self.btn_refresh.setToolTip("Refaz a detecção do Conda, do ambiente 'luna-env' e da instalação do LUNA.")
        self.btn_refresh.clicked.connect(self.detect)
        self.btn_install_miniconda = QPushButton("Baixar Miniconda")
        self.btn_install_miniconda.setToolTip("Abre a página oficial do Miniconda para instalar o Conda no sistema.")
        self.btn_install_miniconda.clicked.connect(self._open_miniconda_download)
        self.btn_install_luna = QPushButton("Instalar LUNA (cria luna-env)")
        self.btn_install_luna.setToolTip(
            "Cria ou atualiza o ambiente 'luna-env' e instala nele as dependências usadas pela GUI."
        )
        self.btn_install_luna.clicked.connect(self.install_luna)
        btn_row.addWidget(self.btn_refresh, 0, 0)
        btn_row.addWidget(self.btn_install_miniconda, 0, 1)
        btn_row.addWidget(self.btn_install_luna, 1, 0, 1, 2)
        btn_row.addWidget(InfoButton("Verifica e prepara o ambiente luna-env. Esta etapa garante LUNA, scikit-learn e dependencias de analise para executar o fluxo."), 1, 2)
        btn_row.setColumnStretch(2, 1)
        env_layout.addLayout(btn_row)
        layout.addWidget(env_box)

        self.log = QPlainTextEdit()
        self.log.setReadOnly(True)
        self.log.setToolTip("Mostra o passo a passo da verificação e da instalação do ambiente.")
        self.log.setFont(QFontDatabase.systemFont(QFontDatabase.SystemFont.FixedFont))
        self.log.setStyleSheet("font-size: 11px;")
        layout.addWidget(self.log, 1)

        self.detect()

    # ---- detection ----
    def _open_miniconda_download(self) -> None:
        try:
            url = em.miniconda_download_url()
        except RuntimeError as exc:
            QMessageBox.critical(self, "Arquitetura não suportada", str(exc))
            return
        webbrowser.open(url)

    def detect(self) -> None:
        self.log.appendPlainText("=== Verificação de ambiente ===")
        if self._detect_task is not None and self._detect_task.is_running:
            return
        self.btn_refresh.setEnabled(False)
        self.btn_install_luna.setEnabled(False)
        self.status_label.setText("Verificando ambiente...")
        task = AsyncTask(self._probe_environment, self)
        self._detect_task = task
        task.succeeded.connect(self._apply_environment_status)
        task.failed.connect(self._detection_failed)
        task.finished.connect(self._detection_finished)
        task.start()

    @staticmethod
    def _probe_environment() -> dict[str, object]:
        logs: list[str] = []
        conda = em.find_conda()
        if not conda:
            return {
                "status": "❌ Conda não encontrado. Instale o Miniconda primeiro.",
                "install_enabled": False,
                "logs": ["Conda: NÃO encontrado."],
            }
        logs.append(f"Conda: {conda}")
        prefix = em.env_prefix(conda)
        logs.append(f"Prefixo do env: {prefix}")

        if em.env_is_partial(conda):
            return {"conda": conda, "status": "⚠️  Conda OK. Ambiente 'luna-env' está incompleto e será recriado ao clicar em 'Instalar LUNA'.", "install_enabled": True, "logs": logs}

        if not em.env_exists(conda):
            return {"conda": conda, "status": "⚠️  Conda OK. Ambiente 'luna-env' não existe — clique em 'Instalar LUNA'.", "install_enabled": True, "logs": logs}

        py = em.env_python(conda)
        if not py:
            return {"conda": conda, "status": "⚠️  luna-env existe mas python não foi localizado.", "install_enabled": True, "logs": logs}
        logs.append(f"Python do env: {py}")

        if not em.luna_installed(py):
            return {"conda": conda, "status": "⚠️  luna-env existe mas LUNA não está instalado.", "install_enabled": True, "logs": logs}

        run_py = em.luna_run_py_path(py)
        if not run_py:
            return {"conda": conda, "status": "⚠️  LUNA importado mas run.py não foi localizado.", "install_enabled": True, "logs": logs}

        missing = em.missing_runtime_packages(py)
        if missing:
            logs.append(
                "Dependências ausentes no luna-env: "
                + ", ".join(missing)
                + "\nUse 'Instalar LUNA' para executar conda install/update com scikit-learn."
            )
            return {"conda": conda, "status": "⚠️  LUNA pronto, mas faltam dependências para análises avançadas: " + ", ".join(missing) + ". Clique em 'Instalar LUNA' para atualizar o luna-env.", "install_enabled": True, "logs": logs, "ready": (str(py), str(run_py))}

        logs.append(f"run.py: {run_py}")
        return {"conda": conda, "status": f"✅ LUNA pronto. run.py: {run_py}", "install_enabled": True, "logs": logs, "ready": (str(py), str(run_py))}

    def _apply_environment_status(self, result: dict[str, object]) -> None:
        self.status_label.setText(str(result.get("status", "")))
        self.btn_install_luna.setEnabled(bool(result.get("install_enabled", False)))
        for line in result.get("logs", []):
            self.log.appendPlainText(str(line))
        conda = result.get("conda")
        if conda:
            self.conda = str(conda)
        ready = result.get("ready")
        if isinstance(ready, tuple) and len(ready) == 2:
            self.luna_ready.emit(str(ready[0]), str(ready[1]))

    def _detection_failed(self, message: str) -> None:
        self.status_label.setText("❌ Falha ao verificar o ambiente.")
        self.log.appendPlainText(f"Erro: {message}")

    def _detection_finished(self) -> None:
        self.btn_refresh.setEnabled(True)
        self._detect_task = None

    def has_running_tasks(self) -> bool:
        detecting = self._detect_task is not None and self._detect_task.is_running
        installing = self.proc is not None and self.proc.state() != QProcess.ProcessState.NotRunning
        return detecting or installing

    # ---- install ----
    def install_luna(self) -> None:
        if self._detect_task is not None and self._detect_task.is_running:
            return
        conda = em.find_conda()
        if not conda:
            QMessageBox.warning(self, "Conda ausente", "Instale o Miniconda primeiro.")
            return
        prefix = em.env_prefix(conda)
        try:
            removed = em.cleanup_partial_env(conda)
        except Exception as exc:
            QMessageBox.critical(
                self,
                "Erro ao preparar ambiente",
                f"Não foi possível limpar o ambiente parcial em:\n{prefix}\n\n{exc}",
            )
            return
        self._cmd_queue = em.install_commands(conda)
        self.btn_install_luna.setEnabled(False)
        self.btn_refresh.setEnabled(False)
        self.log.appendPlainText("\n=== Iniciando instalação do LUNA ===")
        self.log.appendPlainText("Isso pode levar vários minutos.\n")
        self.log.appendPlainText(f"Prefixo alvo do ambiente: {prefix}")
        if removed:
            self.log.appendPlainText(f"Ambiente parcial removido: {removed}")
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
        proc_env = QProcessEnvironment()
        for key, value in em.conda_process_env(cmd[0]).items():
            proc_env.insert(key, value)
        self.proc.setProcessEnvironment(proc_env)
        self.proc.readyReadStandardOutput.connect(self._on_stdout)
        self.proc.finished.connect(self._on_finished)
        self.proc.errorOccurred.connect(self._on_process_error)
        self.proc.start(cmd[0], cmd[1:])

    def _on_stdout(self) -> None:
        sender = self.sender()
        if isinstance(sender, QProcess) and sender is not self.proc:
            return
        if self.proc is None:
            return
        data = bytes(self.proc.readAllStandardOutput()).decode("utf-8", "replace")
        self.log.insertPlainText(data)
        sb = self.log.verticalScrollBar()
        sb.setValue(sb.maximum())

    def _on_finished(self, code: int, _status) -> None:
        sender = self.sender()
        if isinstance(sender, QProcess) and sender is not self.proc:
            return
        proc = self.proc
        self.proc = None
        if proc is not None:
            proc.deleteLater()
        self.log.appendPlainText(f"[exit code: {code}]")
        if code != 0:
            self.log.appendPlainText("⚠️  Comando falhou — abortando instalação.")
            self.btn_install_luna.setEnabled(True)
            self.btn_refresh.setEnabled(True)
            self._cmd_queue.clear()
            return
        self._run_next()

    def _on_process_error(self, error: QProcess.ProcessError) -> None:
        sender = self.sender()
        if isinstance(sender, QProcess) and sender is not self.proc:
            return
        if error != QProcess.ProcessError.FailedToStart:
            return
        proc = self.proc
        detail = proc.errorString() if proc is not None else "processo não iniciado"
        self.log.appendPlainText(f"[erro ao iniciar comando] {detail}")
        self.btn_install_luna.setEnabled(True)
        self.btn_refresh.setEnabled(True)
        self._cmd_queue.clear()
        self.proc = None
        if proc is not None:
            proc.deleteLater()
