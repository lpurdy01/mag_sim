#!/usr/bin/env bash
set -euo pipefail

# Light-weight housekeeping for the Flask GUI runtime assets.
#
# The script prunes stale uploads/results and surfaces dependency
# inconsistencies so environments stay clean between simulation runs.
#
# Usage:
#   ./scripts/maintain_gui_env.sh [days-to-keep]
#   KEEP_DAYS=3 ./scripts/maintain_gui_env.sh

ROOT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
KEEP_DAYS=${1:-${KEEP_DAYS:-7}}
UPLOAD_DIR="${ROOT_DIR}/python/gui/uploads"
RESULTS_DIR="${ROOT_DIR}/python/gui/results"
PYTHON_BIN=${PYTHON:-python3}

prune_dir() {
  local path="$1"
  if [[ -d "$path" ]]; then
    find "$path" -type f -mtime "+${KEEP_DAYS}" -print -delete
  fi
}

prune_dir "$UPLOAD_DIR"
prune_dir "$RESULTS_DIR"

"${PYTHON_BIN}" -m pip check >/dev/null

echo "[maintain_gui_env] Removed artefacts older than ${KEEP_DAYS} day(s)."

outdated=$("${PYTHON_BIN}" -m pip list --outdated)
if [[ $(echo "$outdated" | awk 'END {print NR}') -gt 2 ]]; then
  echo "[maintain_gui_env] The following packages can be updated:" >&2
  echo "$outdated" >&2
else
  echo "[maintain_gui_env] All installed packages are up to date." >&2
fi
