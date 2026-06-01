#!/usr/bin/env bash
set -euo pipefail

APP_DIR="${SNAP}/app"
PYTHON_BIN="${SNAP}/usr/bin/python3"

export PYTHONNOUSERSITE=1
export PYTHONUTF8=1
export PYTHONIOENCODING=utf-8
export MPLCONFIGDIR="${SNAP_USER_COMMON}/matplotlib"
export HIP2LINTERACTOMICS_SNAP=1

mkdir -p "${MPLCONFIGDIR}"
cd "${APP_DIR}"

exec "${PYTHON_BIN}" "${APP_DIR}/run.py" "$@"
