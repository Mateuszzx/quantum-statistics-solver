cmake_minimum_required(VERSION 3.22)

project(
        quantum-statics-solver
        VERSION 1.0.0
        LANGUAGES CXX)


# set(CMAKE_CXX_STANDARD 17)
# set(CMAKE_CXX_STANDARD_REQUIRED ON)
# set(CMAKE_CXX_EXTENSIONS OFF)


# ---- Options (top-level only) ----
if (PROJECT_IS_TOP_LEVEL)
  option(ENABLE_TESTING "Enable unit tests" ON)
  option(ENABLE_WARNINGS "Enable compiler warnings" ON)
  option(ENABLE_WARNINGS_AS_ERRORS "Treat warnings as errors" OFF)
  option(ENABLE_OPENMP "Enable OpenMP acceleration" ON)
  option(ENABLE_PY_MODULE "Build Python extension (requires Python dev)" OFF)
  option(ENABLE_IPO "Enable interprocedural optimization (LTO) for Release" OFF)
endif()


# ---- Module path & helpers ----
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")

if(ENABLE_WARNINGS)
    include(Warnings)
endif()

if (ENABLE_PY_MODULE)
  # Only if you actually build a Python extension (pybind11/nanobind)
  find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)
  # find_package(nanobind CONFIG REQUIRED) # if you use it
endif()

if (ENABLE_OPENMP)
  find_package(OpenMP REQUIRED)
endif()


# ---- Library / app subdirs ----
add_subdirectory(src)       # defines your library target(s); set cxx_std and includes there
add_subdirectory(app)       # defines your executable(s) that link the library
add_subdirectory(configured)  # if this generates headers/sources via configure_file


# ---- Tests ----
if (PROJECT_IS_TOP_LEVEL AND ENABLE_TESTING)
  include(CTest)
  enable_testing()

  include(CPM)
  CPMAddPackage("gh:google/googletest@1.14.0")

  add_subdirectory(tests)
endif()

# ---- IPO / LTO ----
if(PROJECT_IS_TOP_LEVEL AND ENABLE_IPO)
  include(LTO)
  enable_lto(qss_library)
  enable_lto(qss_cli)
endif()