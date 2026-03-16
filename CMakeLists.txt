# ==============================================================================
# Quantum Statistics Solver - Main CMake Configuration
# ==============================================================================
# This is the top-level CMakeLists.txt that orchestrates the entire build.
#
# Build types:
#   - Standalone C++ library/CLI (default)
#   - Python extension module (with -DENABLE_PY_MODULE=ON)
#
# Key options:
#   ENABLE_PY_MODULE   - Build Python bindings (default: OFF)
#   ENABLE_TESTING     - Build tests (default: OFF)
#   ENABLE_COVERAGE    - Enable coverage instrumentation (default: OFF)
#   ENABLE_WARNINGS    - Enable strict compiler warnings (default: OFF)
#   ENABLE_OPENMP      - Enable OpenMP parallelization (default: OFF)
#   ENABLE_IPO         - Enable link-time optimization (default: OFF)
#   BUILD_DOCS         - Build documentation (default: OFF)
#
# Examples:
#   cmake -S . -B build
#   cmake -S . -B build -DENABLE_PY_MODULE=ON -DENABLE_TESTING=ON
#   cmake --build build
# ==============================================================================

cmake_minimum_required(VERSION 3.22...3.28)

# ------------------------------------------------------------------------------
# Project Definition
# ------------------------------------------------------------------------------

project(
  quantum-statistics-solver
  VERSION 1.0.0
  DESCRIPTION "High-performance quantum statistics calculations"
  HOMEPAGE_URL "https://github.com/your-org/quantum-statistics-solver"
  LANGUAGES CXX
)

# Detect if this is the top-level project (not included via add_subdirectory)
if(CMAKE_SOURCE_DIR STREQUAL PROJECT_SOURCE_DIR)
  set(PROJECT_IS_TOP_LEVEL TRUE)
else()
  set(PROJECT_IS_TOP_LEVEL FALSE)
endif()

# ------------------------------------------------------------------------------
# CMake Module Path Configuration
# ------------------------------------------------------------------------------

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")

# Add nanobind cmake path (for Python bindings)
if(EXISTS "${PROJECT_SOURCE_DIR}/ext/nanobind/cmake")
  list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/ext/nanobind/cmake")
endif()

# ------------------------------------------------------------------------------
# Build Options
# ------------------------------------------------------------------------------

option(ENABLE_PY_MODULE "Build Python extension module (_qss)" ON)
option(ENABLE_TESTING "Build C++ and Python tests" OFF)
option(ENABLE_COVERAGE "Enable code coverage instrumentation" OFF)
option(ENABLE_WARNINGS "Enable strict compiler warnings" OFF)
option(ENABLE_IPO "Enable interprocedural optimization (LTO)" OFF)
option(BUILD_DOCS "Build documentation (Doxygen + Sphinx)" OFF)

# Support standard CMake variable names
if(BUILD_TESTING AND NOT ENABLE_TESTING)
    set(ENABLE_TESTING ON CACHE BOOL "Build tests" FORCE)
endif()

# Print configuration summary
message(STATUS "========================================")
message(STATUS "Quantum Statistics Solver Configuration")
message(STATUS "========================================")
message(STATUS "Version:         ${PROJECT_VERSION}")
message(STATUS "Build type:      ${CMAKE_BUILD_TYPE}")
message(STATUS "C++ compiler:    ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")
message(STATUS "Install prefix:  ${CMAKE_INSTALL_PREFIX}")
message(STATUS "")
message(STATUS "Options:")
message(STATUS "  Python module: ${ENABLE_PY_MODULE}")
message(STATUS "  Testing:       ${ENABLE_TESTING}")
message(STATUS "  Coverage:      ${ENABLE_COVERAGE}")
message(STATUS "  Warnings:      ${ENABLE_WARNINGS}")
message(STATUS "  IPO/LTO:       ${ENABLE_IPO}")
message(STATUS "  Documentation: ${BUILD_DOCS}")
message(STATUS "========================================")

# ------------------------------------------------------------------------------
# Include CMake Modules
# ------------------------------------------------------------------------------

# Code coverage must be included early (before add_subdirectory)
# so coverage flags are applied to all targets
if(ENABLE_COVERAGE)
  include(CodeCoverage)
  message(STATUS "Coverage instrumentation enabled")
endif()

# Compiler warnings
if(ENABLE_WARNINGS)
  include(Warnings)
  message(STATUS "Strict compiler warnings enabled")
endif()

# ------------------------------------------------------------------------------
# Python Configuration (if building Python module)
# ------------------------------------------------------------------------------

if(ENABLE_PY_MODULE)
  # Prefer virtual environment Python
  if(EXISTS "${CMAKE_SOURCE_DIR}/.venv/bin/python")
    set(Python_EXECUTABLE "${CMAKE_SOURCE_DIR}/.venv/bin/python" CACHE FILEPATH "venv Python" FORCE)
  endif()

  set(Python_FIND_VIRTUALENV ONLY)

  # CMake 3.18+ component name
  if(CMAKE_VERSION VERSION_LESS 3.18)
    set(DEV_MODULE Development)
  else()
    set(DEV_MODULE Development.Module)
  endif()

  # Find Python with Module target (provides Python::Module)
  find_package(Python 3.10 REQUIRED COMPONENTS Interpreter ${DEV_MODULE})

  message(STATUS "Python interpreter: ${Python_EXECUTABLE}")
  message(STATUS "Python version:     ${Python_VERSION}")
  message(STATUS "Python include dir: ${Python_INCLUDE_DIRS}")

  # Find nanobind
  execute_process(
    COMMAND "${Python_EXECUTABLE}" -m nanobind --cmake_dir
    OUTPUT_VARIABLE nanobind_DIR
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET
  )

  find_package(nanobind CONFIG REQUIRED)
  message(STATUS "Using nanobind from: ${nanobind_DIR}")
endif()

# ------------------------------------------------------------------------------
# Add Project Subdirectories
# ------------------------------------------------------------------------------

# 1. Core C++ library (libqss) - only add once
if(NOT TARGET qss_library)
  add_subdirectory(src)
endif()

# 2. Command-line application (qss_cli)
if(NOT TARGET qss_cli)
  add_subdirectory(app)
endif()

# 3. Configured headers (version.hpp, config.hpp)
if(EXISTS "${PROJECT_SOURCE_DIR}/configured")
  add_subdirectory(configured)
endif()

# 4. Tests (Google Test for C++, pytest for Python)
if(ENABLE_TESTING AND NOT TARGET qss_tests)
  include(CTest)
  enable_testing()

  # Fetch GoogleTest
  message(STATUS "Configuring GoogleTest...")
  set(CPM_LOCAL "${CMAKE_SOURCE_DIR}/cmake/CPM.cmake")
  if(EXISTS "${CPM_LOCAL}")
    include("${CPM_LOCAL}")
    message(STATUS "Using local CPM.cmake")
  else()
    message(STATUS "Downloading CPM.cmake")
    file(DOWNLOAD
      https://github.com/cpm-cmake/CPM.cmake/releases/latest/download/CPM.cmake
      "${CPM_LOCAL}"
    )
    include("${CPM_LOCAL}")
  endif()

  CPMAddPackage("gh:google/googletest@1.14.0")
  message(STATUS "GoogleTest configured")

  # Add tests subdirectory (ONLY ONCE!)
  add_subdirectory(tests)
endif()

# 5. Documentation (Doxygen + Sphinx)
if(BUILD_DOCS)
  include(Docs)
  message(STATUS "Documentation build enabled")
endif()

# ------------------------------------------------------------------------------
# Code Coverage Configuration
# ------------------------------------------------------------------------------
# Coverage must be enabled AFTER all targets are defined

if(ENABLE_COVERAGE)
  # Add coverage instrumentation to the main library
  # This generates coverage data when tests are run
  add_coverage(QSS::library)
  message(STATUS "Coverage instrumentation added to QSS::library")
endif()

# ------------------------------------------------------------------------------
# Link-Time Optimization (LTO/IPO)
# ------------------------------------------------------------------------------
# LTO can significantly improve performance but increases build time

if(ENABLE_IPO)
  include(LTO)

  # Enable LTO for main targets
  enable_lto(qss_library)
  enable_lto(qss_cli)

  # Enable LTO for Python module if it exists
  if(TARGET _qss)
    enable_lto(_qss)
  endif()

  message(STATUS "Link-time optimization enabled")
endif()

# ------------------------------------------------------------------------------
# Installation Rules
# ------------------------------------------------------------------------------
# Install Python module and source files (for scikit-build-core)

if(ENABLE_PY_MODULE)
  # Get platform-specific library directory from scikit-build-core
  # SKBUILD_PLATLIB_DIR is set by scikit-build during pip install
  if(NOT DEFINED SKBUILD_PLATLIB_DIR)
    set(SKBUILD_PLATLIB_DIR "${CMAKE_INSTALL_PREFIX}/lib/python${Python_VERSION_MAJOR}.${Python_VERSION_MINOR}/site-packages")
  endif()

  # Install Python source files to the `qss` package namespace.
  # The compiled `_qss` extension is installed in src/bindings/CMakeLists.txt.
  install(
    DIRECTORY ${PROJECT_SOURCE_DIR}/qss/
    DESTINATION ${SKBUILD_PLATLIB_DIR}/qss
    FILES_MATCHING
    PATTERN "*.py"
    PATTERN "__pycache__" EXCLUDE
    PATTERN "*.pyc" EXCLUDE
    PATTERN ".pytest_cache" EXCLUDE
  )
endif()

# ------------------------------------------------------------------------------
# Export Package Configuration
# ------------------------------------------------------------------------------
# This allows other CMake projects to find this package with find_package()

if(NOT ENABLE_PY_MODULE)
  include(GNUInstallDirs)
  include(CMakePackageConfigHelpers)

  # Install C++ library and headers
  install(
    TARGETS qss_library
    EXPORT qss-targets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  )

  # Install public headers
  install(
    DIRECTORY ${PROJECT_SOURCE_DIR}/include/
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    FILES_MATCHING PATTERN "*.hpp" PATTERN "*.h"
  )

  # Generate and install package configuration files
  write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/qss-config-version.cmake"
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
  )

  install(
    EXPORT qss-targets
    FILE qss-targets.cmake
    NAMESPACE QSS::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/qss
  )

  configure_file(
    cmake/qss-config.cmake.in
    "${CMAKE_CURRENT_BINARY_DIR}/qss-config.cmake"
    @ONLY
  )

  install(
    FILES
      "${CMAKE_CURRENT_BINARY_DIR}/qss-config.cmake"
      "${CMAKE_CURRENT_BINARY_DIR}/qss-config-version.cmake"
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/qss
  )
endif()

# ==============================================================================
# End of Main CMakeLists.txt
# ==============================================================================
