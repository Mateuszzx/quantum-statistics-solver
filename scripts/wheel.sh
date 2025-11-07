#!/usr/bin/env bash
set -euo pipefail

VENV_DIR="${VENV_DIR:-.venv}"
PY="${VENV_DIR}/bin/python"
echo "[info] Building wheel"
"${PY}" -m build --wheel
