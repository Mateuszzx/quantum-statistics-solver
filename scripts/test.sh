#!/usr/bin/env bash
set -euo pipefail

# Precedence: env BUILD_TYPE > first arg > Release
BUILD_TYPE="${BUILD_TYPE:-${1:-Release}}"
VENV_DIR="${VENV_DIR:-.venv}"
PY="${VENV_DIR}/bin/python"
BUILD_DIR="build/${BUILD_TYPE}"

# Ensure venv Python exists
if [[ ! -x "${PY}" ]]; then
  echo "[err] Missing venv Python at ${PY}. Create it first: python -m venv .venv"
  exit 1
fi

# Configure if build dir missing
if [[ ! -d "${BUILD_DIR}" || ! -f "${BUILD_DIR}/CMakeCache.txt" ]]; then
  echo "[info] Configuring (${BUILD_TYPE})"
  ./scripts/configure.sh "${BUILD_TYPE}"
fi

echo "[info] Building target 'qss' (${BUILD_TYPE})"
cmake --build "${BUILD_DIR}" --target qss -j

# Install into the venv's platlib (site-packages for extensions), not /usr/local
PLATLIB="$("${PY}" -c 'import sysconfig; print(sysconfig.get_paths()["platlib"])')"
echo "[info] Installing into: ${PLATLIB}"
cmake --install "${BUILD_DIR}" --prefix "${PLATLIB}"

# Ensure runtime deps (numpy) before import
if ! "${PY}" -c "import numpy" >/dev/null 2>&1; then
  echo "[info] Installing numpy (runtime dependency)"
  "${PY}" -m pip install -U pip
  "${PY}" -m pip install numpy
fi

echo "[info] Import smoke test"
if ! "${PY}" -c "import qss; print('Imported qss from', qss.__file__)"; then
  echo "[warn] Import failed from site-packages; trying build dir"
  export PYTHONPATH="$PWD/${BUILD_DIR}/python:${PYTHONPATH:-}"
  "${PY}" -c "import qss; print('Imported qss from (fallback)', qss.__file__)"
fi
