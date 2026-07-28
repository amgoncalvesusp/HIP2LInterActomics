#!/usr/bin/env bash
set -euo pipefail

APP_DIR="${SNAP}/app"
PYTHON_BIN="${SNAP}/usr/bin/python3"

export PYTHONNOUSERSITE=1
export PYTHONUTF8=1
export PYTHONIOENCODING=utf-8
export MPLCONFIGDIR="${SNAP_USER_COMMON}/matplotlib"
export HIP2LINTERACTOMICS_SNAP=1

QT_PLUGIN_DIR="$(find "${SNAP}/usr/lib" -type d -path '*/qt6/plugins' -print -quit 2>/dev/null || true)"
if [[ -n "${QT_PLUGIN_DIR}" ]]; then
    export QT_PLUGIN_PATH="${QT_PLUGIN_DIR}"
    QT_LIB_DIR="$(dirname "$(dirname "${QT_PLUGIN_DIR}")")"
    export LD_LIBRARY_PATH="${QT_LIB_DIR}:${LD_LIBRARY_PATH:-}"
fi

mkdir -p "${MPLCONFIGDIR}"
cd "${APP_DIR}"

exec "${PYTHON_BIN}" "${APP_DIR}/run.py" "$@"
