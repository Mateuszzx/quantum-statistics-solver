#!/usr/bin/env bash
set -euo pipefail

# -----------------------------------------------------------------------------
# setup_nanobind.sh
#
# One-shot setup for nanobind-based C++/Python builds on Ubuntu or Arch Linux.
# - Installs system packages (CMake, Ninja, compilers, Python dev)
# - Creates/uses .venv (unless --system-python given)
# - Installs nanobind into the chosen Python
# - Configures CMake to use venv Python + nanobind_DIR
# - Builds, and optionally installs the module into the venv's site-packages
#
# Usage:
#   ./setup_nanobind.sh [--release] [--relwithdebinfo] [--no-install] [--clean]
#                       [--system-python] [--build-dir <dir>] [--target <tgt>]
#
# Examples:
#   ./setup_nanobind.sh
#   ./setup_nanobind.sh --release
#   ./setup_nanobind.sh --no-install --build-dir build/dev
#   ./setup_nanobind.sh --system-python --target qss
# -----------------------------------------------------------------------------

BUILD_TYPE="Debug"
DO_INSTALL=1
DO_CLEAN=0
USE_SYSTEM_PY=0
BUILD_DIR="build/${BUILD_TYPE}"
TARGET=""   # default: build all

while [[ $# -gt 0 ]]; do
  case "$1" in
    --release) BUILD_TYPE="Release"; shift ;;
    --relwithdebinfo) BUILD_TYPE="RelWithDebInfo"; shift ;;
    --no-install) DO_INSTALL=0; shift ;;
    --clean) DO_CLEAN=1; shift ;;
    --system-python) USE_SYSTEM_PY=1; shift ;;
    --build-dir) BUILD_DIR="$2"; shift 2 ;;
    --target) TARGET="$2"; shift 2 ;;
    -h|--help)
      grep '^# ' "$0" | sed 's/^# //'; exit 0 ;;
    *) echo "Unknown arg: $1"; exit 1 ;;
  esac
done

detect_distro() {
  if [[ -r /etc/os-release ]]; then
    . /etc/os-release
    echo "${ID:-unknown}"
  else
    echo "unknown"
  fi
}

need_sudo() {
  command -v sudo >/dev/null 2>&1 && echo "sudo" || echo ""
}

install_deps() {
  local id="$1"
  local sudocmd
  sudocmd="$(need_sudo)"

  case "$id" in
    ubuntu|debian)
      echo "[info] Installing deps for Ubuntu/Debian..."
      ${sudocmd} apt-get update -y
      ${sudocmd} apt-get install -y \
        build-essential cmake ninja-build \
        python3 python3-venv python3-dev python3-pip
      ;;
    arch)
      echo "[info] Installing deps for Arch..."
      ${sudocmd} pacman -Sy --needed --noconfirm \
        base-devel cmake ninja python python-pip
      # python venv is in the 'python' package on Arch
      ;;
    *)
      echo "[warn] Unknown distro '$id' — skipping system package install."
      echo "       Please ensure: CMake, Ninja, a C++ compiler, Python dev headers are installed."
      ;;
  esac
}

ensure_python() {
  if (( USE_SYSTEM_PY )); then
    PYBIN="$(command -v python || true)"
    [[ -n "${PYBIN}" ]] || { echo "[err] 'python' not found in PATH"; exit 1; }
    echo "[info] Using system Python: ${PYBIN}"
  else
    if [[ ! -d ".venv" ]]; then
      echo "[info] Creating .venv ..."
      python -m venv .venv
    fi
    PYBIN="$(pwd)/.venv/bin/python"
    [[ -x "${PYBIN}" ]] || { echo "[err] venv python not found at ${PYBIN}"; exit 1; }
    echo "[info] Using venv Python: ${PYBIN}"
    "${PYBIN}" -m pip install -U pip wheel
  fi
}

ensure_nanobind() {
  if ! "${PYBIN}" - <<'PY' >/dev/null 2>&1
import importlib; importlib.import_module("nanobind")
PY
  then
    echo "[info] Installing nanobind into $("${PYBIN}" -c 'import sys; print(sys.executable)') ..."
    "${PYBIN}" -m pip install -U nanobind
  else
    echo "[info] nanobind already installed."
  fi
}

configure_and_build() {
  local pyexe
  local nbdir
  pyexe="$("${PYBIN}" -c 'import sys; print(sys.executable)')"
  nbdir="$("${PYBIN}" -c 'import nanobind; print(nanobind.cmake_dir())')"

  # Select generator
  local gen=""
  if command -v ninja >/dev/null 2>&1; then
    gen="-G Ninja"
    echo "[info] Using Ninja generator"
  else
    echo "[info] Using default CMake generator"
  fi

  # Clean/build dir
  BUILD_DIR="${BUILD_DIR%/}"  # trim trailing slash
  if (( DO_CLEAN )) && [[ -d "${BUILD_DIR}" ]]; then
    echo "[info] Cleaning ${BUILD_DIR}"
    rm -rf "${BUILD_DIR}"
  fi
  mkdir -p "${BUILD_DIR}"

  echo "[info] Configuring CMake:"
  echo "       BUILD_TYPE=${BUILD_TYPE}"
  echo "       BUILD_DIR=${BUILD_DIR}"
  echo "       Python_EXECUTABLE=${pyexe}"
  echo "       nanobind_DIR=${nbdir}"
  echo "       QSS_INSTALL_TO_VENV=$(( DO_INSTALL ? 1 : 0 ))"

  cmake -S . -B "${BUILD_DIR}" ${gen} \
    -DCMAKE_BUILD_TYPE="${BUILD_TYPE}" \
    -DENABLE_PY_MODULE=ON \
    -DQSS_INSTALL_TO_VENV="$(( DO_INSTALL ? 1 : 0 ))" \
    -DPython_EXECUTABLE="${pyexe}" \
    -Dnanobind_DIR="${nbdir}"

  echo "[info] Building..."
  if [[ -n "${TARGET}" ]]; then
    cmake --build "${BUILD_DIR}" --target "${TARGET}" -j
  else
    cmake --build "${BUILD_DIR}" -j
  fi

  if (( DO_INSTALL )); then
    echo "[info] Installing into Python site-packages ..."
    cmake --install "${BUILD_DIR}"
  else
    echo "[info] Skipping install (artifacts in ${BUILD_DIR}/python)"
  fi

  echo "[info] Import smoke test:"
  "${PYBIN}" - <<'PY'
try:
    import qss
    print("OK: imported qss from:", qss.__file__)
except Exception as e:
    import sys
    print("ERROR: failed to import qss:", e)
    sys.exit(1)
PY
}

main() {
  local distro
  distro="$(detect_distro)"
  echo "[info] Detected distro: ${distro}"

  install_deps "${distro}"
  ensure_python
  ensure_nanobind
  configure_and_build

  echo "[done] nanobind setup/build complete."
  echo "       - Build dir: ${BUILD_DIR}"
  if (( DO_INSTALL )); then
    echo "       - Module installed into: $(${PYBIN} -c 'import sysconfig; print(sysconfig.get_paths().get(\"platlib\") or sysconfig.get_paths().get(\"purelib\"))')"
  else
    echo "       - To import without install: export PYTHONPATH=\"$PWD/${BUILD_DIR}/python:\$PYTHONPATH\""
  fi
}

main "$@"
