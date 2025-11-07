#!/usr/bin/env bash
set -euo pipefail

BUILD_TYPE="${1:-Release}"
echo "[info] Building target 'qss' (${BUILD_TYPE})"
cmake --build "build/${BUILD_TYPE}" --target qss -j
