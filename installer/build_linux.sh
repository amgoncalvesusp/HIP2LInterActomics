#!/usr/bin/env bash
# Build a native Linux PyInstaller bundle and a portable tar.gz archive.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
OUTPUT_ROOT="${1:-${REPO_ROOT}/build/release-linux}"
PYTHON_BIN="${PYTHON:-python3}"

if [[ "$(uname -s)" != "Linux" ]]; then
    echo "[ERRO] O bundle Linux precisa ser gerado em um host Linux." >&2
    exit 1
fi

VENV_DIR="${OUTPUT_ROOT}/build-venv"
SITE_DIR="${OUTPUT_ROOT}/build-site"
WORK_DIR="${OUTPUT_ROOT}/pyinstaller-work"
BUNDLE_ROOT="${OUTPUT_ROOT}/bundle"
ARCHIVE_ROOT="${OUTPUT_ROOT}/installer"

mkdir -p "${OUTPUT_ROOT}" "${WORK_DIR}" "${BUNDLE_ROOT}" "${ARCHIVE_ROOT}"
if "${PYTHON_BIN}" -c "import ensurepip" >/dev/null 2>&1 \
    && "${PYTHON_BIN}" -m venv "${VENV_DIR}"; then
    BUILD_PYTHON="${VENV_DIR}/bin/python"
    "${BUILD_PYTHON}" -m pip install --upgrade pip
    "${BUILD_PYTHON}" -m pip install -r "${REPO_ROOT}/requirements.txt" "pyinstaller>=6.10,<7"
    PYTHONPATH_PREFIX=""
else
    echo "[AVISO] python3-venv indisponivel; usando dependencias privadas em ${SITE_DIR}." >&2
    mkdir -p "${SITE_DIR}"
    if ! PYTHONPATH="${SITE_DIR}" "${PYTHON_BIN}" -c \
        "import PyQt6, matplotlib, numpy, scipy, sklearn, PyInstaller" >/dev/null 2>&1; then
        "${PYTHON_BIN}" -m pip install --upgrade --target "${SITE_DIR}" \
            -r "${REPO_ROOT}/requirements.txt" "pyinstaller>=6.10,<7"
    fi
    BUILD_PYTHON="${PYTHON_BIN}"
    PYTHONPATH_PREFIX="${SITE_DIR}"
fi

cd "${REPO_ROOT}"
PYTHONPATH="${PYTHONPATH_PREFIX}${PYTHONPATH_PREFIX:+:}${PYTHONPATH:-}" "${BUILD_PYTHON}" -m PyInstaller \
    --noconfirm \
    --clean \
    --workpath "${WORK_DIR}" \
    --distpath "${BUNDLE_ROOT}" \
    HIP2LInterActomics.spec

EXECUTABLE="${BUNDLE_ROOT}/HIP2LInterActomics/HIP2LInterActomics"
if [[ ! -x "${EXECUTABLE}" ]]; then
    echo "[ERRO] Binario Linux nao gerado: ${EXECUTABLE}" >&2
    exit 1
fi

ENVIRONMENT_FILE="${BUNDLE_ROOT}/HIP2LInterActomics/environment.yml"
if [[ ! -f "${ENVIRONMENT_FILE}" ]]; then
    echo "[ERRO] environment.yml nao foi incluido no bundle: ${ENVIRONMENT_FILE}" >&2
    exit 1
fi
for cli_file in hipplinteractomics_terminal.py hipplinteractomics_multiple_run.py; do
    cli_path="${BUNDLE_ROOT}/HIP2LInterActomics/${cli_file}"
    if [[ ! -f "${cli_path}" ]]; then
        echo "[ERRO] ${cli_file} nao foi incluido no bundle: ${cli_path}" >&2
        exit 1
    fi
done

ARCH="$(uname -m)"
ARCHIVE="${ARCHIVE_ROOT}/HIP2LInterActomics-1.0.0-linux-${ARCH}.tar.gz"
tar -C "${BUNDLE_ROOT}" -czf "${ARCHIVE}" HIP2LInterActomics

echo "Binario: ${EXECUTABLE}"
echo "Ambiente terminal-only: ${ENVIRONMENT_FILE}"
echo "Pacote: ${ARCHIVE}"
