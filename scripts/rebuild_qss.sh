#!/usr/bin/env bash
set -euo pipefail

# -----------------------------------------------------------------------------
# rebuild_qss.sh
# Rebuilds the nanobind-based Python extension (qss) for the current project.
#
# Features:
#  - Uses/creates .venv (unless --system-python is used)
#  - Ensures nanobind is installed and exposes nanobind_DIR to CMake
#  - Configures CMake for Debug/Release (Ninja if available)
#  - Builds to build/<type>/ and installs into the venv's site-packages (on by default)
#  - Optionally cleans build dir, skips install, or uses system Python
#
# Usage:
#   ./rebuild_qss.sh [--clean] [--type Debug|Release|RelWithDebInfo] [--no-install] [--system-python]
#
# Examples:
#   ./rebuild_qss.sh
#   ./rebuild_qss.sh --type Release
#   ./rebuild_qss.sh --clean --type RelWithDebInfo
#   ./rebuild_qss.sh --no-install
#   ./rebuild_qss.sh --system-python  # uses current python without creating/activating .venv
# -----------------------------------------------------------------------------

BUILD_TYPE="Debug"
DO_CLEAN=0
DO_INSTALL=1
USE_SYSTEM_PY=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --clean) DO_CLEAN=1; shift ;;
    --no-install) DO_INSTALL=0; shift ;;
    --type)
      [[ $# -ge 2 ]] || { echo "Missing arg for --type"; exit 1; }
      BUILD_TYPE="$2"; shift 2 ;;
    --system-python) USE_SYSTEM_PY=1; shift ;;
    -h|--help)
      grep '^# ' "$0" | sed 's/^# //'; exit 0 ;;
    *)
      echo "Unknown arg: $1"; exit 1 ;;
  esac
done

# 1) Ensure we have a Python to use
if (( USE_SYSTEM_PY )); then
  PYBIN="$(command -v python)"
  if [[ -z "${PYBIN}" ]]; then
    echo "python not found in PATH"; exit 1
  fi
  echo "[info] Using system Python: ${PYBIN}"
else
  # Create .venv if missing; don’t force-activate for fish compatibility — we call it explicitly
  if [[ ! -d ".venv" ]]; then
    echo "[info] Creating .venv ..."
    python -m venv .venv
  fi
  PYBIN="$(pwd)/.venv/bin/python"
  if [[ ! -x "${PYBIN}" ]]; then
    echo "Expected venv python at ${PYBIN} not found or not executable"; exit 1
  fi
  echo "[info] Using venv Python: ${PYBIN}"

  # Upgrade pip/wheel and ensure nanobind is present
  "${PYBIN}" -m pip install -q -U pip wheel
fi

# 2) Ensure nanobind is installed in the chosen interpreter
if ! "${PYBIN}" - <<'PY' >/dev/null 2>&1
import importlib; importlib.import_module("nanobind")
PY
then
  echo "[info] Installing nanobind into $( "${PYBIN}" -c 'import sys; print(sys.executable)' ) ..."
  "${PYBIN}" -m pip install -U nanobind
fi

# 3) Gather Python/CMake hints
PYEXE="$("${PYBIN}" -c 'import sys; print(sys.executable)')"
NANOBIND_DIR="$("${PYBIN}" -c 'import nanobind; print(nanobind.cmake_dir())')"

# 4) Select generator
if command -v ninja >/dev/null 2>&1; then
  GENERATOR="-G Ninja"
  echo "[info] Using Ninja generator"
else
  GENERATOR=""
  echo "[info] Using default CMake generator"
fi

# 5) Build directory per config
BUILD_DIR="build/${BUILD_TYPE}"
if (( DO_CLEAN )) && [[ -d "${BUILD_DIR}" ]]; then
  echo "[info] Cleaning ${BUILD_DIR}"
  rm -rf "${BUILD_DIR}"
fi
mkdir -p "${BUILD_DIR}"

# 6) Configure
echo "[info] Configuring CMake in ${BUILD_DIR} (type: ${BUILD_TYPE})"
cmake -S . -B "${BUILD_DIR}" ${GENERATOR} \
  -DCMAKE_BUILD_TYPE="${BUILD_TYPE}" \
  -DENABLE_PY_MODULE=ON \
  -DQSS_INSTALL_TO_VENV="$(( DO_INSTALL ? 1 : 0 ))" \
  -DPython_EXECUTABLE="${PYEXE}" \
  -Dnanobind_DIR="${NANOBIND_DIR}"

# 7) Build
echo "[info] Building ..."
cmake --build "${BUILD_DIR}" -j

# 8) Optionally install into site-packages
if (( DO_INSTALL )); then
  echo "[info] Installing into Python site-packages ..."
  cmake --install "${BUILD_DIR}"
else
  echo "[info] Skipping install (artifacts in ${BUILD_DIR}/python)"
fi

# 9) Smoke test import
echo "[info] Import test ..."
"${PYBIN}" - <<'PY'
import sys
try:
    import qss
    print("OK: imported qss from:", qss.__file__)
except Exception as e:
    print("ERROR: failed to import qss:", e)
    sys.exit(1)
PY

echo "[done] Rebuild complete."
