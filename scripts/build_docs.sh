#!/usr/bin/env bash
set -euo pipefail

# Determine repository root (the parent directory of this script).
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${SCRIPT_DIR}/.."


print_usage() {
  cat <<'EOF'
Usage: scripts/build_docs.sh [--skip-install]

Builds the MkDocs documentation with strict validation enabled. By default the
script installs the documentation dependencies listed in requirements-docs.txt.
Pass --skip-install when you already have the dependencies available (for
example after creating a virtual environment).
EOF
}

INSTALL_DEPS=1

while [[ $# -gt 0 ]]; do
  case "$1" in
    --skip-install)
      INSTALL_DEPS=0
      shift
      ;;
    -h|--help)
      print_usage
      exit 0
      ;;
    *)
      echo "Unknown argument: $1" >&2
      print_usage >&2
      exit 1
      ;;
  esac
done

if [[ $INSTALL_DEPS -eq 1 ]]; then
  python3 -m pip install --requirement requirements-docs.txt
fi

python3 tools/check_docs_math_blocks.py

mkdocs build --strict
