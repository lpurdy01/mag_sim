#!/usr/bin/env bash
set -euo pipefail

# Bootstrap dependencies required to run the Flask GUI end-to-end Playwright suite.
# Usage:
#   ./scripts/setup_gui_e2e_env.sh [python-executable]
#   PYTHON=python3.11 ./scripts/setup_gui_e2e_env.sh
#
# The script installs Flask GUI runtime deps plus pytest-playwright and the
# Chromium browser binary used by the E2E test.

ROOT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
PYTHON_BIN=${1:-${PYTHON:-python3}}

"${PYTHON_BIN}" -m pip install --upgrade pip
"${PYTHON_BIN}" -m pip install --upgrade \
  flask \
  matplotlib \
  numpy \
  ezdxf \
  pytest \
  pytest-playwright \
  playwright

"${PYTHON_BIN}" -m playwright install chromium

mkdir -p "${ROOT_DIR}/python/gui/uploads" "${ROOT_DIR}/python/gui/results"

cat <<SETUP_MSG
[setup_gui_e2e_env] Installed Flask GUI + Playwright dependencies using ${PYTHON_BIN}.
[setup_gui_e2e_env] Chromium browser binaries are available for headless E2E runs.
[setup_gui_e2e_env] Upload/result directories are ready under python/gui/.
SETUP_MSG
