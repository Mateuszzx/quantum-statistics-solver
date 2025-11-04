cmake_minimum_required(VERSION 3.22)

project(quantum-statistics-solver
  VERSION 1.0.0
  LANGUAGES CXX)


if (PROJECT_IS_TOP_LEVEL)
  option(ENABLE_TESTING "Enable unit tests" ON)
  option(ENABLE_WARNINGS "Enable compiler warnings" ON)
  option(ENABLE_WARNINGS_AS_ERRORS "Treat warnings as errors" OFF)
  option(ENABLE_OPENMP "Enable OpenMP acceleration" ON)
  option(ENABLE_PY_MODULE "Build Python extension (requires Python dev)" ON)
  option(ENABLE_IPO "Enable interprocedural optimization (LTO) for Release" OFF)
endif()


list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/ext/nanobind/cmake")
if (ENABLE_WARNINGS)
  include(Warnings)
endif()

if (ENABLE_OPENMP)
  find_package(OpenMP REQUIRED)
endif()


if (ENABLE_PY_MODULE)

  if (CMAKE_VERSION VERSION_LESS 3.18)
    set(DEV_MODULE Development)
  else()
    set(DEV_MODULE Development.Module)
  endif()

  find_package(Python 3.8 COMPONENTS Interpreter ${DEV_MODULE} REQUIRED)

  cmake_policy(SET CMP0148 NEW)
  find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)
  add_subdirectory(ext/nanobind)

  # Detect the installed nanobind package and import it into CMake
  execute_process(
    COMMAND "${Python_EXECUTABLE}" -m nanobind --cmake_dir
    OUTPUT_STRIP_TRAILING_WHITESPACE OUTPUT_VARIABLE nanobind_ROOT)

    find_package(nanobind CONFIG REQUIRED)

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

if (ENABLE_IPO)
  include(LTO)
  enable_lto(qss_library)
  enable_lto(qss_cli)
  if (TARGET qss)  # Python module target exists only when ENABLE_PY_MODULE
    enable_lto(qss)
  endif()
endif()
