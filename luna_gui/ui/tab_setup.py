"""Setup tab — detect/install conda + luna-env, manage LUNA installation."""
from __future__ import annotations

import webbrowser
from pathlib import Path

from PyQt6.QtCore import Qt, QSize, QProcess, QProcessEnvironment, QThread, pyqtSignal
from PyQt6.QtGui import QFontDatabase, QImageReader, QPixmap
from PyQt6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QPushButton, QPlainTextEdit,
    QMessageBox, QGroupBox,
)

from ..core import env_manager as em
from .info import InfoButton


def _probe_environment() -> dict[str, object]:
    """Inspect the LUNA runtime without touching Qt widgets."""
    logs = ["=== Verificação de ambiente ==="]
    result: dict[str, object] = {
        "logs": logs,
        "install_enabled": False,
        "ready": None,
    }
    conda = em.find_conda()
    if not conda:
        result["status"] = "Erro — Conda não encontrado. Instale o Miniconda primeiro."
        logs.append("Conda: NÃO encontrado.")
        return result

    result["conda"] = conda
    result["install_enabled"] = True
    logs.append(f"Conda: {conda}")
    runtime = em.find_luna_runtime(conda)
    if runtime is not None:
        python_path, run_py = runtime
        prefix = em.python_prefix(python_path)
        logs.append(f"Prefixo do env: {prefix}")
        logs.append(f"Python do env: {python_path}")
        logs.append(f"run.py: {run_py}")
        missing = em.missing_runtime_packages(python_path)
        result["ready"] = (str(python_path), str(run_py))
        if missing:
            result["status"] = (
                "Atenção — LUNA pronto, mas faltam dependências para análises avançadas: "
                + ", ".join(missing)
                + ". Clique em 'Instalar LUNA' para atualizar o luna-env."
            )
            logs.append(
                "Dependências ausentes no luna-env: "
                + ", ".join(missing)
                + "\nUse 'Instalar LUNA' para executar conda install/update com scikit-learn."
            )
            return result
        result["status"] = f"LUNA pronto — run.py: {run_py}"
        return result

    prefix = em.env_prefix(conda)
    logs.append(f"Prefixo do env: {prefix}")

    if em.env_is_partial(conda):
        result["status"] = (
            "Atenção — Conda OK. Ambiente 'luna-env' está incompleto e será recriado ao "
            "clicar em 'Instalar LUNA'."
        )
        return result
    if not em.env_exists(conda):
        result["status"] = (
            "Atenção — Conda OK. Ambiente 'luna-env' não existe — clique em 'Instalar LUNA'."
        )
        return result

    python_path = em.env_python(conda)
    if not python_path:
        result["status"] = "Atenção — luna-env existe mas python não foi localizado."
        return result
    logs.append(f"Python do env: {python_path}")

    if not em.luna_installed(python_path):
        result["status"] = "Atenção — luna-env existe mas LUNA não está instalado."
        return result
    run_py = em.luna_run_py_path(python_path)
    if not run_py:
        result["status"] = "Atenção — LUNA importado mas run.py não foi localizado."
        return result

    missing = em.missing_runtime_packages(python_path)
    result["ready"] = (str(python_path), str(run_py))
    if missing:
        result["status"] = (
            "Atenção — LUNA pronto, mas faltam dependências para análises avançadas: "
            + ", ".join(missing)
            + ". Clique em 'Instalar LUNA' para atualizar o luna-env."
        )
        logs.append(
            "Dependências ausentes no luna-env: "
            + ", ".join(missing)
            + "\nUse 'Instalar LUNA' para executar conda install/update com scikit-learn."
        )
        return result

    result["status"] = f"LUNA pronto — run.py: {run_py}"
    logs.append(f"run.py: {run_py}")
    return result


class _EnvironmentProbeThread(QThread):
    detected = pyqtSignal(object)

    def run(self) -> None:
        try:
            result = _probe_environment()
        except Exception as exc:
            result = {
                "status": f"Erro — Falha ao verificar o ambiente: {exc}",
                "logs": ["=== Verificação de ambiente ===", repr(exc)],
                "install_enabled": True,
                "ready": None,
            }
        self.detected.emit(result)


class SetupTab(QWidget):
    """Detects conda, creates luna-env, installs LUNA."""

    luna_ready = pyqtSignal(str, str)  # (env_python_path, run_py_path)

    def __init__(self) -> None:
        super().__init__()
        self.proc: QProcess | None = None
        self._probe_thread: _EnvironmentProbeThread | None = None
        self._cmd_queue: list[list[str]] = []
        self._last_ready: tuple[str, str] | None = None

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
            reader = QImageReader(str(icon_path))
            reader.setAutoTransform(True)
            reader.setScaledSize(QSize(190, 190))
            image = reader.read()
            if not image.isNull():
                icon_label.setPixmap(QPixmap.fromImage(image))
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
        workflow.setObjectName("workflowSummary")
        about_layout.addWidget(workflow)
        env_box = QGroupBox("Configuração do ambiente e instalação de pacotes")
        env_layout = QVBoxLayout(env_box)
        env_layout.addWidget(self.status_label)

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
        btn_row.addWidget(InfoButton("Verifica e prepara o ambiente luna-env. Esta etapa garante LUNA, scikit-learn e dependencias de analise para executar o fluxo."))
        btn_row.addStretch()
        env_layout.addLayout(btn_row)
        layout.addWidget(env_box)
        layout.addWidget(about_box, 2)

        self.log = QPlainTextEdit()
        self.log.setReadOnly(True)
        self.log.setToolTip("Mostra o passo a passo da verificação e da instalação do ambiente.")
        log_font = QFontDatabase.systemFont(QFontDatabase.SystemFont.FixedFont)
        log_font.setPointSize(10)
        self.log.setFont(log_font)
        layout.addWidget(self.log, 1)

    # ---- detection ----
    def detect(self) -> None:
        if self._probe_thread is not None and self._probe_thread.isRunning():
            return
        self.status_label.setText("Verificando ambiente em segundo plano...")
        self.btn_refresh.setEnabled(False)
        thread = _EnvironmentProbeThread(self)
        thread.detected.connect(self._apply_probe_result)
        thread.finished.connect(self._probe_finished)
        self._probe_thread = thread
        thread.start()

    def _apply_probe_result(self, result: dict[str, object]) -> None:
        for line in result.get("logs", []):
            self.log.appendPlainText(str(line))
        self.status_label.setText(str(result.get("status", "Verificação concluída.")))
        self.btn_install_luna.setEnabled(bool(result.get("install_enabled", False)))
        conda = result.get("conda")
        if conda:
            self.conda = str(conda)
        ready = result.get("ready")
        if isinstance(ready, tuple) and len(ready) == 2:
            self._last_ready = (str(ready[0]), str(ready[1]))
            self.luna_ready.emit(*self._last_ready)
        else:
            self._last_ready = None

    def ready_runtime(self) -> tuple[str, str] | None:
        """Return the most recent verified runtime, if one is available."""
        return self._last_ready

    def _probe_finished(self) -> None:
        thread = self._probe_thread
        self._probe_thread = None
        self.btn_refresh.setEnabled(True)
        if thread is not None:
            thread.deleteLater()

    # ---- install ----
    def install_luna(self) -> None:
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
        self.proc.start(cmd[0], cmd[1:])

    def _on_stdout(self) -> None:
        data = bytes(self.proc.readAllStandardOutput()).decode("utf-8", "replace")
        self.log.insertPlainText(data)
        sb = self.log.verticalScrollBar()
        sb.setValue(sb.maximum())

    def _on_finished(self, code: int, _status) -> None:
        self.log.appendPlainText(f"[exit code: {code}]")
        if code != 0:
            self.log.appendPlainText("Atenção — comando falhou; instalação interrompida.")
            self.btn_install_luna.setEnabled(True)
            self._cmd_queue.clear()
            return
        self._run_next()
