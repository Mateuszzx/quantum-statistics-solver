# src/bindings/CMakeLists.txt

# Only build the Python extension when requested
if (NOT ENABLE_PY_MODULE)
  return()
endif()

# At the root you've already tried find_package(Python ...) and nanobind (or fetched it).
# Still, be robust if this directory is configured standalone.
if (NOT TARGET nanobind::nanobind)
  find_package(nanobind REQUIRED CONFIG)
endif()

# Some older CMake setups may not bring in Python::Module; ensure it's there.
if (NOT TARGET Python::Module)
  find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)
endif()

# ---- Sources ----
# List your binding sources here. Only existing files will be added.
set(QSS_PY_SOURCES
  module_main.cpp
  module_chemical_potential.cpp
  module_distributions.cpp
)

set(_existing)
foreach(f IN LISTS QSS_PY_SOURCES)
  if (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${f}")
    list(APPEND _existing "${f}")
  endif()
endforeach()

# Option to enable CPython limited (stable) ABI
option(QSS_STABLE_ABI "Build Python module with CPython limited/stable ABI" OFF)

set(_nb_args MODULE NB_STATIC)         # Static self-contained extension
if (QSS_STABLE_ABI)
  list(APPEND _nb_args STABLE_ABI)
endif()

# ---- Build the nanobind module ----
nanobind_add_module(qss
  ${_nb_args}
  ${_existing}
)

# Language / warnings / includes
set_target_properties(qss PROPERTIES
  CXX_STANDARD 17
  CXX_STANDARD_REQUIRED YES
  POSITION_INDEPENDENT_CODE ON
  CXX_VISIBILITY_PRESET hidden
  VISIBILITY_INLINES_HIDDEN YES
)

target_include_directories(qss PRIVATE
  ${PROJECT_SOURCE_DIR}/include
)

# Link against your core library (which already pulls OpenMP if enabled)
target_link_libraries(qss PRIVATE QSS::library)

# ----------------------------------------------------------------------------
# Developer convenience: drop the built extension into a predictable folder
# so it's easy to import from the build tree (add build/python to PYTHONPATH).
set_target_properties(qss PROPERTIES
  LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/python
  RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/python
)

# ----------------------------------------------------------------------------
# Optional: install the extension directly into the current Python's site-packages
# (useful when building inside a virtualenv). Toggle with -DQSS_INSTALL_TO_VENV=ON/OFF
option(QSS_INSTALL_TO_VENV "Install qss Python extension into current Python site-packages" ON)
if (QSS_INSTALL_TO_VENV)
  # Determine the Python executable used for the build
  if (NOT Python_EXECUTABLE)
    if (TARGET Python::Interpreter)
      get_target_property(_qss_pyexe Python::Interpreter IMPORTED_LOCATION)
      set(Python_EXECUTABLE "${_qss_pyexe}")
    endif()
  endif()

  if (Python_EXECUTABLE)
    # Query site-packages path (prefer platlib for binary extensions; fall back to purelib)
    execute_process(
      COMMAND ${Python_EXECUTABLE} -c "import sysconfig,sys; p=sysconfig.get_paths(); sys.stdout.write(p.get('platlib') or p.get('purelib') or '')"
      OUTPUT_VARIABLE QSS_PY_SITE
      OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    if (QSS_PY_SITE)
      message(STATUS "qss: will install Python extension to ${QSS_PY_SITE}")
      # Absolute DESTINATION is allowed; handle both ELF/.so (LIBRARY) and Windows .pyd (RUNTIME)
      install(TARGETS qss
        LIBRARY DESTINATION "${QSS_PY_SITE}"
        RUNTIME DESTINATION "${QSS_PY_SITE}"
      )
    else()
      message(WARNING "qss: could not determine Python site-packages; skipping install to venv. Set QSS_INSTALL_TO_VENV=OFF or provide Python_EXECUTABLE.")
    endif()
  else()
    message(WARNING "qss: Python_EXECUTABLE not known; skipping install to venv. Set -DPython_EXECUTABLE or enable Python::Interpreter.")
  endif()
endif()

