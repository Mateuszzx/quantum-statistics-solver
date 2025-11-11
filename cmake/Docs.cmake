# ==============================================================================
# Documentation Configuration
# ==============================================================================
# Builds documentation using Doxygen (C++) and Sphinx (Python).
# Output: build/docs/html/index.html
# ==============================================================================

find_package(Doxygen REQUIRED dot)
find_package(Sphinx REQUIRED)

if(DOXYGEN_FOUND)
    # Configure Doxygen
    set(DOXYGEN_IN ${CMAKE_SOURCE_DIR}/docs/Doxyfile.in)
    set(DOXYGEN_OUT ${CMAKE_BINARY_DIR}/Doxyfile)

    configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)

    # Add documentation target
    add_custom_target(docs
        COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
        COMMENT "Generating API documentation with Doxygen"
        VERBATIM
    )

    message(STATUS "Documentation target 'docs' configured")
    message(STATUS "  Build with: cmake --build build --target docs")
    message(STATUS "  Output: ${CMAKE_BINARY_DIR}/docs/html/index.html")
endif()
