# ==============================================================================
# Chemical Potential Solver Module
# ==============================================================================
# Algorithms for computing chemical potential from particle density.
# ==============================================================================

set(CHEMICAL_POTENTIAL_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/chemical_potential.cpp
    PARENT_SCOPE  # Make visible to parent CMakeLists.txt
)
