# ==============================================================================
# Distribution Functions Module
# ==============================================================================

set(DISTRIBUTIONS_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/distributions.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/fd.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/flb.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/nfl.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/nflb.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/nflb_i.cpp
    PARENT_SCOPE  # Make visible to parent CMakeLists.txt
)
