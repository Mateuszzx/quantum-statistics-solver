# ==============================================================================
# Distribution Functions Module
# ==============================================================================

set(DISTRIBUTIONS_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/distributions.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/fd.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/flb.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/nfl.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/nfl_reworked.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/nflb.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/ssl.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/hubbard.cpp
    PARENT_SCOPE  # Make visible to parent CMakeLists.txt
)
