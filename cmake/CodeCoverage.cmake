# Guard against multiple inclusion
if (DEFINED CODECOVERAGE_INCLUDED)
  return()
endif()
set(CODECOVERAGE_INCLUDED TRUE)

if (NOT DEFINED ENABLE_COVERAGE)
  set(ENABLE_COVERAGE OFF)
endif()

if (NOT ENABLE_COVERAGE)
  return()
endif()

message(STATUS "CodeCoverage: ENABLE_COVERAGE=ON — enabling coverage flags and target")

# Apply safe coverage-friendly compile/link flags for GNU/Clang
if (CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  add_compile_options(--coverage -O0 -g -fno-inline -fno-omit-frame-pointer)
  # ensure link flags propagate
  string(APPEND CMAKE_EXE_LINKER_FLAGS " --coverage")
  string(APPEND CMAKE_SHARED_LINKER_FLAGS " --coverage")
else()
  message(WARNING "CodeCoverage: coverage flags not applied for this compiler (${CMAKE_CXX_COMPILER_ID})")
endif()

# Prepare output dir
set(COVERAGE_OUTPUT_DIR "${CMAKE_BINARY_DIR}/coverage")
file(MAKE_DIRECTORY "${COVERAGE_OUTPUT_DIR}")

# If a coverage target already exists (maybe from another include), skip creating it here.
if (TARGET coverage)
  message(STATUS "CodeCoverage: 'coverage' target already exists — skipping target creation")
  # still define helper function so add_coverage() works
else()
  # Coverage target: run tests, then try gcovr, otherwise lcov/genhtml fallback
  add_custom_target(coverage
    COMMENT "Run tests and generate coverage report (html in ${COVERAGE_OUTPUT_DIR})"
  )

  # The coverage command uses a small shell snippet so we can pick gcovr or lcov at build time.
  add_custom_command(TARGET coverage
    PRE_BUILD
      COMMAND ${CMAKE_COMMAND} -E echo "Running tests (ctest)..."
      COMMAND /bin/sh -c "LD_LIBRARY_PATH='${CMAKE_BINARY_DIR}':\\$LD_LIBRARY_PATH ${CMAKE_CTEST_COMMAND} --output-on-failure || true"
      COMMAND ${CMAKE_COMMAND} -E echo "Generating coverage report (src/ only)..."
      COMMAND /bin/sh -c " \
        if command -v gcovr >/dev/null 2>&1; then \
          gcovr -r '${CMAKE_SOURCE_DIR}' \
            --object-directory '${CMAKE_BINARY_DIR}' \
            --filter '${CMAKE_SOURCE_DIR}/src/.*' \
            --exclude '${CMAKE_SOURCE_DIR}/src/bindings/.*' \
            --exclude '.*nanobind.*' \
            --exclude '.*build/.*' \
            --exclude '.*tests/.*' \
            --html --html-details \
            -o '${COVERAGE_OUTPUT_DIR}/index.html'; \
        else \
          if ! command -v lcov >/dev/null 2>&1; then echo 'ERROR: install gcovr or lcov'; exit 1; fi; \
          lcov --rc geninfo_unexecuted_blocks=1 --directory '${CMAKE_BINARY_DIR}' --capture --output-file '${COVERAGE_OUTPUT_DIR}/coverage.info' || true; \
          lcov --remove '${COVERAGE_OUTPUT_DIR}/coverage.info' \
            '/usr/*' \
            '*/external/*' \
            '*/build/*' \
            '*/tests/*' \
            '*/src/bindings/*' \
            '*nanobind*' \
            '*/app/main.cpp' \
            --output-file '${COVERAGE_OUTPUT_DIR}/coverage.cleaned.info' || true; \
          genhtml '${COVERAGE_OUTPUT_DIR}/coverage.cleaned.info' --output-directory '${COVERAGE_OUTPUT_DIR}/html' || true; \
        fi"
    WORKING_DIRECTORY "${CMAKE_BINARY_DIR}"
    VERBATIM
  )
endif()

# Convenience function: ensure coverage target depends on build/test targets (so they build first)
function(add_coverage tgt_in)
  if (NOT ENABLE_COVERAGE)
    return()
  endif()

  if (NOT TARGET ${tgt_in})
    message(FATAL_ERROR "add_coverage: target '${tgt_in}' does not exist")
  endif()

  # Resolve alias to its real target, if needed
  get_target_property(_aliased ${tgt_in} ALIASED_TARGET)
  if (_aliased)
    set(tgt ${_aliased})
  else()
    set(tgt ${tgt_in})
  endif()

  add_dependencies(coverage ${tgt})
endfunction()
