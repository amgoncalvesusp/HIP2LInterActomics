#!/usr/bin/env bash
# ============================================================
#  HIP2LInterActomics_GUI launcher — Linux
# ============================================================
#  Main distribution path for Linux:
#    - locates conda robustly
#    - activates the GUI env
#    - validates the required plotting/clustering modules
#    - launches run.py from the repo root
# ============================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"

ENV_NAME="${HIP2LINTERACTOMICS_GUI_ENV:-${LUNA_GUI_ENV:-luna-gui}}"
GUI_PYTHON="${HIP2LINTERACTOMICS_GUI_PYTHON:-${LUNA_GUI_PYTHON:-}}"

if [[ -n "${GUI_PYTHON}" ]]; then
    if [[ ! -x "${GUI_PYTHON}" ]]; then
        echo "[ERRO] HIP2LINTERACTOMICS_GUI_PYTHON aponta para um executável inexistente:"
        echo "       ${GUI_PYTHON}"
        exit 1
    fi
    exec "${GUI_PYTHON}" "${REPO_ROOT}/run.py" "$@"
fi

if command -v conda >/dev/null 2>&1; then
    CONDA_BIN="$(command -v conda)"
else
    for candidate in \
        "$HOME/miniconda3/bin/conda" \
        "$HOME/anaconda3/bin/conda" \
        "/opt/miniconda3/bin/conda" \
        "/opt/anaconda3/bin/conda"; do
        if [[ -x "$candidate" ]]; then
            CONDA_BIN="$candidate"
            break
        fi
    done
fi

if [[ -z "${CONDA_BIN:-}" ]]; then
    echo "[ERRO] conda não encontrado. Instale Miniconda primeiro:"
    echo "       https://docs.conda.io/en/latest/miniconda.html"
    exit 1
fi

CONDA_BASE="$("${CONDA_BIN}" info --base)"
CONDA_SH="${CONDA_BASE}/etc/profile.d/conda.sh"
if [[ ! -f "${CONDA_SH}" ]]; then
    echo "[ERRO] não foi possível localizar conda.sh em:"
    echo "       ${CONDA_SH}"
    exit 1
fi

# shellcheck disable=SC1090
source "${CONDA_SH}"

if ! conda activate "${ENV_NAME}" >/dev/null 2>&1; then
    echo "[ERRO] env '${ENV_NAME}' não encontrado."
    echo "       Crie com: conda create -n ${ENV_NAME} python=3.11 -y"
    echo "                 conda activate ${ENV_NAME}"
    echo "                 pip install -r ${REPO_ROOT}/requirements.txt"
    exit 1
fi

if ! python - <<'PY'
import importlib.util
missing = [mod for mod in ("PyQt6", "matplotlib", "numpy", "scipy") if importlib.util.find_spec(mod) is None]
raise SystemExit(0 if not missing else 1)
PY
then
    echo "[ERRO] O env '${ENV_NAME}' está ativo, mas faltam dependências da GUI."
    echo "       Rode: pip install -r ${REPO_ROOT}/requirements.txt"
    exit 1
fi

exec python "${REPO_ROOT}/run.py" "$@"
