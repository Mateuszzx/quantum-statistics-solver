# ==============================================================================
# Core C++ Library (libqss)
# ==============================================================================
# This file builds the main quantum statistics solver library as a static lib.
# The library is later linked by:
#   - The command-line application (qss_cli)
#   - The Python extension module (_qss)
#   - The test suite (qss_tests)
# ==============================================================================

# Prevent multiple inclusions if this directory is added multiple times
include_guard(GLOBAL)

# ------------------------------------------------------------------------------
# Collect Source Files from Subdirectories
# ------------------------------------------------------------------------------

# Core utility sources (in this directory)
set(LIBRARY_SOURCES
    math_utils.cpp
)

# Add module subdirectories - each sets a variable with its sources
# These must be added before we can append their sources to LIBRARY_SOURCES
add_subdirectory(distributions)       # Sets DISTRIBUTIONS_SOURCES
add_subdirectory(chemical_potential)  # Sets CHEMICAL_POTENTIAL_SOURCES
add_subdirectory(solver)              # Sets SOLVER_SOURCES

# Combine all module sources into one list
list(APPEND LIBRARY_SOURCES
    ${DISTRIBUTIONS_SOURCES}
    ${CHEMICAL_POTENTIAL_SOURCES}
    ${SOLVER_SOURCES}
)

# ------------------------------------------------------------------------------
# Create Static Library Target
# ------------------------------------------------------------------------------

# Build as static library (could be SHARED if needed)
add_library(qss_library STATIC ${LIBRARY_SOURCES})

# Create alias for clean namespacing: allows "QSS::library" in target_link_libraries
add_library(QSS::library ALIAS qss_library)

# ------------------------------------------------------------------------------
# Configure Library Target
# ------------------------------------------------------------------------------

# Require C++17 for all consumers of this library
target_compile_features(qss_library PUBLIC cxx_std_17)

# Set include directories
# PUBLIC: these paths are used both when building AND when linking this library
target_include_directories(qss_library
  PUBLIC
    # When building: use source tree headers
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    # When installed: use installed headers
    $<INSTALL_INTERFACE:include>
  PRIVATE
    # Private: only needed when building this library
    ${CMAKE_CURRENT_SOURCE_DIR}
)

# Link optional dependencies
if(ENABLE_OPENMP)
  # OpenMP for parallel computation
  target_link_libraries(qss_library PUBLIC OpenMP::OpenMP_CXX)
endif()

# Add compiler warnings (if enabled via -DENABLE_WARNINGS=ON)
if(ENABLE_WARNINGS)
  target_compile_options(qss_library PRIVATE
    # GCC/Clang warnings
    $<$<CXX_COMPILER_ID:GNU,Clang>:-Wall -Wextra -Wpedantic>
    # MSVC warnings
    $<$<CXX_COMPILER_ID:MSVC>:/W4>
  )
endif()

# Position-independent code (required for shared libraries and Python extensions)
set_target_properties(qss_library PROPERTIES
  POSITION_INDEPENDENT_CODE ON
  OUTPUT_NAME "qss"  # Library will be named libqss.a (or qss.lib on Windows)
)

# ------------------------------------------------------------------------------
# Python Bindings (Conditional)
# ------------------------------------------------------------------------------

# Only build Python bindings if explicitly enabled
if(ENABLE_PY_MODULE)
  add_subdirectory(bindings)
endif()

# ==============================================================================
# End of src/CMakeLists.txt
# ==============================================================================
