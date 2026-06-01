#!/usr/bin/env bash
# Build the Linux .snap package from snap/snapcraft.yaml.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"

cd "${REPO_ROOT}"

if [[ "$(uname -s)" != "Linux" ]]; then
    echo "[ERRO] Snapcraft gera .snap apenas em Linux. Use uma maquina Linux ou WSL2 com snapcraft/LXD." >&2
    exit 1
fi

if ! command -v snapcraft >/dev/null 2>&1; then
    echo "[ERRO] snapcraft nao encontrado." >&2
    echo "Instale em Ubuntu/Debian com: sudo snap install snapcraft --classic" >&2
    exit 1
fi

if [[ ! -f "snap/snapcraft.yaml" ]]; then
    echo "[ERRO] snap/snapcraft.yaml nao encontrado." >&2
    exit 1
fi

snapcraft "$@"
