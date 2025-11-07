#!/usr/bin/env bash
set -euo pipefail

BUILD_TYPE="${1:-Release}"
VENV_DIR="${VENV_DIR:-.venv}"

if [[ ! -x "${VENV_DIR}/bin/python" ]]; then
  echo "[err] Python not found at ${VENV_DIR}/bin/python"; exit 1
fi

PYEXE="${VENV_DIR}/bin/python"
NANOBIND_DIR="$("${PYEXE}" -c 'import nanobind; print(nanobind.cmake_dir())')"

GEN=""
if command -v ninja >/dev/null 2>&1; then
  GEN="-G Ninja"
  echo "[info] Using Ninja generator"
else
  echo "[info] Using default CMake generator"
fi

echo "[info] Configuring (BUILD_TYPE=${BUILD_TYPE})"
cmake -S . -B "build/${BUILD_TYPE}" ${GEN} \
  -DCMAKE_BUILD_TYPE="${BUILD_TYPE}" \
  -DENABLE_PY_MODULE=ON \
  -DQSS_INSTALL_TO_VENV=ON \
  -DPython_EXECUTABLE="${PYEXE}" \
  -Dnanobind_DIR="${NANOBIND_DIR}"
