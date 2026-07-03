#!/usr/bin/env bash
# Install HIP2LInterActomics on a fresh Linux workstation.
#
# The script installs Miniforge when conda is not available, creates the
# luna-gui and luna-env environments, installs the GUI requirements and installs
# LUNA with the compatibility patch used by the application.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"

GUI_ENV="${HIP2L_GUI_ENV:-luna-gui}"
LUNA_ENV="${HIP2L_LUNA_ENV:-luna-env}"
CONDA_ROOT="${HIP2L_CONDA_ROOT:-${HOME}/.hip2linteractomics/miniforge3}"
INSTALL_LUNA=1
CREATE_DESKTOP=1

usage() {
    cat <<'USAGE'
HIP2LInterActomics Linux installer

Usage:
  ./dist/linux/install_hip2linteractomics.sh [options]

Options:
  --gui-only       Install only luna-gui. Skip luna-env/LUNA.
  --skip-luna      Same as --gui-only.
  --no-desktop     Do not create a desktop launcher.
  --conda-root DIR Install Miniforge in DIR when conda is missing.
  -h, --help       Show this help.

Environment overrides:
  HIP2L_GUI_ENV      GUI env name. Default: luna-gui
  HIP2L_LUNA_ENV     LUNA env name. Default: luna-env
  HIP2L_CONDA_ROOT   Miniforge install path when conda is missing.
USAGE
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --gui-only|--skip-luna)
            INSTALL_LUNA=0
            shift
            ;;
        --no-desktop)
            CREATE_DESKTOP=0
            shift
            ;;
        --conda-root)
            CONDA_ROOT="$2"
            shift 2
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "[ERRO] Opcao desconhecida: $1" >&2
            usage
            exit 2
            ;;
    esac
done

if [[ ! -f "${REPO_ROOT}/requirements.txt" || ! -f "${REPO_ROOT}/run.py" ]]; then
    echo "[ERRO] Rode este instalador a partir de uma copia completa do repositorio." >&2
    exit 1
fi

log() {
    printf '\n[HIP2L] %s\n' "$*"
}

find_conda() {
    local candidate
    if [[ -n "${HIP2LINTERACTOMICS_GUI_CONDA:-}" && -x "${HIP2LINTERACTOMICS_GUI_CONDA}" ]]; then
        printf '%s\n' "${HIP2LINTERACTOMICS_GUI_CONDA}"
        return 0
    fi
    if command -v conda >/dev/null 2>&1; then
        command -v conda
        return 0
    fi
    for candidate in \
        "${CONDA_ROOT}/bin/conda" \
        "${HOME}/miniforge3/bin/conda" \
        "${HOME}/miniconda3/bin/conda" \
        "${HOME}/anaconda3/bin/conda" \
        "${HOME}/mambaforge/bin/conda" \
        "/opt/miniforge3/bin/conda" \
        "/opt/miniconda3/bin/conda" \
        "/opt/anaconda3/bin/conda"; do
        if [[ -x "${candidate}" ]]; then
            printf '%s\n' "${candidate}"
            return 0
        fi
    done
    return 1
}

download_file() {
    local url="$1"
    local output="$2"
    if command -v curl >/dev/null 2>&1; then
        curl -L --fail --output "${output}" "${url}"
    elif command -v wget >/dev/null 2>&1; then
        wget -O "${output}" "${url}"
    else
        echo "[ERRO] Instale curl ou wget para baixar o Miniforge automaticamente." >&2
        exit 1
    fi
}

install_miniforge() {
    local arch installer url
    arch="$(uname -m)"
    case "${arch}" in
        x86_64|amd64) arch="x86_64" ;;
        aarch64|arm64) arch="aarch64" ;;
        *)
            echo "[ERRO] Arquitetura Linux nao suportada automaticamente: ${arch}" >&2
            exit 1
            ;;
    esac

    mkdir -p "$(dirname "${CONDA_ROOT}")"
    installer="$(mktemp -t miniforge-installer.XXXXXX.sh)"
    url="https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-${arch}.sh"

    log "Conda nao encontrado. Baixando Miniforge: ${url}"
    download_file "${url}" "${installer}"
    bash "${installer}" -b -p "${CONDA_ROOT}"
    rm -f "${installer}"
}

env_exists() {
    "${CONDA_BIN}" run -n "$1" python -V >/dev/null 2>&1
}

conda_create_or_install() {
    local env_name="$1"
    shift
    if env_exists "${env_name}"; then
        "${CONDA_BIN}" install -n "${env_name}" --override-channels -c conda-forge -y "$@"
    else
        "${CONDA_BIN}" create -n "${env_name}" --override-channels -c conda-forge -y "$@"
    fi
}

install_gui_env() {
    log "Criando/atualizando ambiente ${GUI_ENV}"
    conda_create_or_install "${GUI_ENV}" python=3.11 pip
    "${CONDA_BIN}" run -n "${GUI_ENV}" python -m pip install --upgrade pip
    "${CONDA_BIN}" run -n "${GUI_ENV}" python -m pip install -r "${REPO_ROOT}/requirements.txt"
}

install_luna_env() {
    log "Criando/atualizando ambiente ${LUNA_ENV}"
    conda_create_or_install "${LUNA_ENV}" \
        python=3.9 pip rdkit openbabel pymol-open-source biopython=1.79 \
        numpy pandas scipy scikit-learn matplotlib seaborn networkx

    "${CONDA_BIN}" run -n "${LUNA_ENV}" python -s -m pip install pdbecif 'mmh3<4' xopen colorlog
    "${CONDA_BIN}" run -n "${LUNA_ENV}" python -s -m pip install --no-build-isolation -U luna
    "${CONDA_BIN}" run -n "${LUNA_ENV}" python -s "${REPO_ROOT}/luna_gui/core/_luna_patch.py"
}

create_desktop_launcher() {
    local desktop_dir desktop_file icon_path
    desktop_dir="${XDG_DATA_HOME:-${HOME}/.local/share}/applications"
    desktop_file="${desktop_dir}/hip2linteractomics.desktop"
    icon_path="${REPO_ROOT}/luna_gui/assets/hip2l_interactomics_icon.png"
    mkdir -p "${desktop_dir}"
    cat > "${desktop_file}" <<EOF
[Desktop Entry]
Type=Application
Name=HIP2LInterActomics
Comment=Protein-ligand interaction analysis with LUNA
Exec=env HIP2LINTERACTOMICS_GUI_ENV=${GUI_ENV} HIP2LINTERACTOMICS_GUI_CONDA=${CONDA_BIN} ${REPO_ROOT}/dist/linux/run_gui.sh
Icon=${icon_path}
Terminal=false
Categories=Science;Education;
EOF
    chmod +x "${desktop_file}"
    log "Atalho criado em ${desktop_file}"
}

if ! CONDA_BIN="$(find_conda)"; then
    install_miniforge
    CONDA_BIN="${CONDA_ROOT}/bin/conda"
fi

if [[ ! -x "${CONDA_BIN}" ]]; then
    echo "[ERRO] Conda nao foi encontrado apos a instalacao: ${CONDA_BIN}" >&2
    exit 1
fi

log "Usando conda: ${CONDA_BIN}"
"${CONDA_BIN}" config --set channel_priority flexible >/dev/null 2>&1 || true

install_gui_env
if [[ "${INSTALL_LUNA}" -eq 1 ]]; then
    install_luna_env
else
    log "Instalacao do luna-env pulada por opcao do usuario."
fi

chmod +x "${REPO_ROOT}/dist/linux/run_gui.sh"
if [[ "${CREATE_DESKTOP}" -eq 1 ]]; then
    create_desktop_launcher
fi

log "Instalacao concluida."
echo "Abra com:"
echo "  ${REPO_ROOT}/dist/linux/run_gui.sh"
echo
echo "Paralelismo: no Linux voce pode usar varios nucleos em nproc."
