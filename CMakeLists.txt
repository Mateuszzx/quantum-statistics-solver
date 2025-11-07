cmake_minimum_required(VERSION 3.22)

project(quantum-statistics-solver
  VERSION 1.0.0
  LANGUAGES CXX)

# mark top-level project so options are defined
if (CMAKE_SOURCE_DIR STREQUAL PROJECT_SOURCE_DIR)
  set(PROJECT_IS_TOP_LEVEL TRUE)
endif()

# make coverage option available at configure time
option(ENABLE_COVERAGE "Enable coverage report target" OFF)

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/ext/nanobind/cmake")

# If coverage requested, include helper early so flags apply to all targets
if (ENABLE_COVERAGE)
  include(${CMAKE_SOURCE_DIR}/cmake/CodeCoverage.cmake)
endif()

if (ENABLE_WARNINGS)
  include(Warnings)
endif()

if (ENABLE_OPENMP)
  find_package(OpenMP REQUIRED)
endif()

if (ENABLE_PY_MODULE)
  if (EXISTS "${CMAKE_SOURCE_DIR}/.venv/bin/python")
    set(Python_EXECUTABLE "${CMAKE_SOURCE_DIR}/.venv/bin/python" CACHE FILEPATH "Preferred venv Python" FORCE)
  endif()

  set(Python_FIND_VIRTUALENV ONLY)
  if (CMAKE_VERSION VERSION_LESS 3.18)
    set(DEV_MODULE Development)
  else()
    set(DEV_MODULE Development.Module)
  endif()

  # Add Development to obtain Python::Python
  find_package(Python 3.8 REQUIRED COMPONENTS Interpreter Development ${DEV_MODULE})
  message(STATUS "Using Python interpreter: ${Python_EXECUTABLE}")

  cmake_policy(SET CMP0148 NEW)
  execute_process(
    COMMAND "${Python_EXECUTABLE}" -m nanobind --cmake_dir
    OUTPUT_VARIABLE nanobind_DIR
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET
  )
  find_package(nanobind CONFIG QUIET)
  if (NOT nanobind_FOUND)
    list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/ext/nanobind/cmake")
    find_package(nanobind CONFIG REQUIRED)
  endif()
endif()

# Prefer building library and headers first
add_subdirectory(src)
add_subdirectory(app)
add_subdirectory(configured)

if (ENABLE_TESTING)
  include(CTest)
  enable_testing()

  # Prefer a local CPM.cmake if available (good for offline/CI). Otherwise bootstrap CPM.
  set(CPM_LOCAL "${CMAKE_SOURCE_DIR}/cmake/CPM.cmake")
  if (EXISTS "${CPM_LOCAL}")
    include(${CPM_LOCAL})
  else()
    include(FetchContent)
    FetchContent_Declare(
      CPM
      URL https://github.com/cpm-cmake/CPM.cmake/releases/latest/download/CPM.cmake
    )
    FetchContent_MakeAvailable(CPM)
  endif()

  # Add GoogleTest via CPM (pinned). Use the gh: shorthand or full URL.
  CPMAddPackage("gh:google/googletest@1.14.0")

  add_subdirectory(tests)
endif()

if (ENABLE_COVERAGE)
  add_coverage(QSS::library)   # instrument the code under test
endif()



if (ENABLE_IPO)
  include(LTO)
  enable_lto(qss_library)
  enable_lto(qss_cli)
  if (TARGET _qss)
    enable_lto(_qss)
  endif()
endif()

# Install Python module (only if target exists and not editable install)
if(TARGET _qss AND NOT SKBUILD_STATE STREQUAL "editable")
  install(TARGETS _qss LIBRARY DESTINATION ${SKBUILD_PLATLIB_DIR}/qss)
endif()

# Always install Python source files
install(
  DIRECTORY ${PROJECT_SOURCE_DIR}/qss/
  DESTINATION ${SKBUILD_PLATLIB_DIR}/qss
  FILES_MATCHING
  PATTERN "*.py"
  PATTERN "__pycache__" EXCLUDE
  PATTERN "*.pyc" EXCLUDE
)
