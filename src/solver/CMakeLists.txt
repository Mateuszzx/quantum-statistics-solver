# ==============================================================================
# Root-Finding Solver Module
# ==============================================================================
# Generic root-finding algorithms:
#   - Bracketing methods (bisection, false position)
#   - Brent's method (hybrid root finder)
#   - Simulated annealing (global optimization)
# ==============================================================================

set(SOLVER_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/solver.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/bracket.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/solver_sa.cpp
    PARENT_SCOPE  # Make visible to parent CMakeLists.txt
)
