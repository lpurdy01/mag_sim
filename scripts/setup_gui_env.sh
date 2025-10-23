#!/usr/bin/env bash
set -euo pipefail

# Bootstrap dependencies required by the Flask GUI prototype.
#
# Usage:
#   ./scripts/setup_gui_env.sh [python-executable]
#   PYTHON=python3.11 ./scripts/setup_gui_env.sh
#
# When no interpreter is supplied the script defaults to the current
# environment's "python3". The command installs Flask, Matplotlib and
# NumPy so the GUI can render previews and field maps. It also prepares the
# runtime directories used to store uploads and generated artefacts.

ROOT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
PYTHON_BIN=${1:-${PYTHON:-python3}}

"${PYTHON_BIN}" -m pip install --upgrade pip
"${PYTHON_BIN}" -m pip install --upgrade \
  flask \
  matplotlib \
  numpy

mkdir -p "${ROOT_DIR}/python/gui/uploads" "${ROOT_DIR}/python/gui/results"

cat <<SETUP_MSG
[setup_gui_env] Installed Flask GUI runtime dependencies using ${PYTHON_BIN}.
[setup_gui_env] Upload/result directories are ready under python/gui/.
SETUP_MSG
