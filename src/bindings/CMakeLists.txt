# src/bindings/CMakeLists.txt

# ==============================================================================
# Python Bindings (nanobind)
# ==============================================================================
# Creates the Python extension module '_qss' that wraps the C++ library.
# Uses nanobind (https://github.com/wjakob/nanobind) for Python<->C++ bindings.
#
# The module structure:
#   - module_main.cpp:               Module initialization, version info
#   - module_distributions.cpp:      Distribution functions (fd, nflb, etc.)
#   - module_chemical_potential.cpp: Chemical potential solvers
# ==============================================================================

# Prevent multiple inclusions (critical for binding modules!)
include_guard(GLOBAL)

# Early exit if Python module building is disabled
if(NOT ENABLE_PY_MODULE)
  return()
endif()

# Sanity check: Python should have been found in main CMakeLists.txt
if(NOT Python_FOUND)
  message(FATAL_ERROR "Python not found but ENABLE_PY_MODULE=ON")
endif()

# Ensure nanobind is available (should be found in main CMakeLists.txt)
if(NOT TARGET nanobind::nanobind)
  find_package(nanobind REQUIRED CONFIG)
endif()

# ------------------------------------------------------------------------------
# Collect Binding Source Files
# ------------------------------------------------------------------------------

set(BINDING_SOURCES
    module_main.cpp               # Module entry point (NB_MODULE macro)
    module_distributions.cpp      # Wraps qss::fd, qss::nflb, etc.
    module_chemical_potential.cpp # Wraps qss::solve_mu_grid, options classes
    module_solver.cpp
)

# Verify all sources exist (helps catch typos)
foreach(src IN LISTS BINDING_SOURCES)
  if(NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${src}")
    message(WARNING "Binding source not found: ${src}")
  endif()
endforeach()

# ------------------------------------------------------------------------------
# Create Python Extension Module
# ------------------------------------------------------------------------------

# nanobind_add_module creates a Python extension module target
# STABLE_ABI: Use Python stable ABI (forward compatibility)
# NB_STATIC: Embed nanobind statically (no runtime dependency)
nanobind_add_module(
  _qss          # Target name (will produce _qss.cpython-312-x86_64-linux-gnu.so)
  STABLE_ABI    # Use Python's stable ABI for forward compatibility
  NB_STATIC     # Link nanobind statically
  ${BINDING_SOURCES}
)

# ------------------------------------------------------------------------------
# Link Libraries
# ------------------------------------------------------------------------------

target_link_libraries(_qss PRIVATE
    QSS::library     # Our C++ library with all the quantum statistics code
    Python::Module   # Python C API (provides _Py_NoneStruct, PyFloat_FromDouble, etc.)
)

# ------------------------------------------------------------------------------
# Include Directories
# ------------------------------------------------------------------------------

target_include_directories(_qss PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}  # For bindings.hpp
    ${PROJECT_SOURCE_DIR}/include # For qss/*.hpp headers
)

# ------------------------------------------------------------------------------
# Target Properties
# ------------------------------------------------------------------------------

set_target_properties(_qss PROPERTIES
    OUTPUT_NAME "_qss"              # Output file will be _qss.so (not lib_qss.so)
    PREFIX ""                       # No 'lib' prefix on Unix
    CXX_STANDARD 17                 # C++17 required
    CXX_STANDARD_REQUIRED YES
    POSITION_INDEPENDENT_CODE ON    # Required for shared libraries
    CXX_VISIBILITY_PRESET hidden    # Hide symbols by default (smaller binary)
    VISIBILITY_INLINES_HIDDEN YES
)

# Set output directory (for development builds)
# Final location: build/python/qss/_qss.cpython-312-x86_64-linux-gnu.so
set_target_properties(_qss PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/python/qss"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/python/qss"
)

# ------------------------------------------------------------------------------
# Optional: Code Coverage Support
# ------------------------------------------------------------------------------

if(ENABLE_COVERAGE)
  target_compile_options(_qss PRIVATE --coverage)
  target_link_options(_qss PRIVATE --coverage)
endif()

# ------------------------------------------------------------------------------
# Optional: Strip Symbols in Release Mode
# ------------------------------------------------------------------------------

# Strip debug symbols to reduce binary size (Release builds only)
if(CMAKE_BUILD_TYPE STREQUAL "Release")
  if(UNIX AND NOT APPLE)
    add_custom_command(TARGET _qss POST_BUILD
      COMMAND ${CMAKE_STRIP} --strip-unneeded $<TARGET_FILE:_qss>
      COMMENT "Stripping symbols from _qss module"
    )
  endif()
endif()

# ------------------------------------------------------------------------------
# Installation Rules
# ------------------------------------------------------------------------------

# Install the Python extension to the appropriate site-packages directory.
# Keep this active for editable installs as well, otherwise `pip/uv install -e .`
# will expose Python sources but miss the compiled `_qss` module.
if(DEFINED SKBUILD_PLATLIB_DIR)
  install(TARGETS _qss
    LIBRARY DESTINATION "${SKBUILD_PLATLIB_DIR}/qss"
    RUNTIME DESTINATION "${SKBUILD_PLATLIB_DIR}/qss"
  )
else()
  # Fallback: use Python's platlib location (for 'make install')
  execute_process(
    COMMAND "${Python_EXECUTABLE}" -c
      "import sysconfig; print(sysconfig.get_path('platlib'))"
    OUTPUT_VARIABLE Python_SITEARCH
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )
  install(TARGETS _qss
    LIBRARY DESTINATION "${Python_SITEARCH}/qss"
    RUNTIME DESTINATION "${Python_SITEARCH}/qss"
  )
endif()

# ------------------------------------------------------------------------------
# Status Messages
# ------------------------------------------------------------------------------

message(STATUS "Python module '_qss' configured")
message(STATUS "  Output: ${CMAKE_BINARY_DIR}/python/qss/_qss${CMAKE_SHARED_MODULE_SUFFIX}")

# ==============================================================================
# End of src/bindings/CMakeLists.txt
# ==============================================================================
