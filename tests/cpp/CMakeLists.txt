# tests/cpp/CMakeLists.txt
# CMakeLists.txt for C++ tests (Google Test)

# Find or fetch Google Test
find_package(GTest QUIET)
if(NOT GTest_FOUND)
  include(FetchContent)
  FetchContent_Declare(
    googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG        v1.14.0
  )
  # Enable test discovery and reuse
  include(GoogleTest)
  set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
  FetchContent_MakeAvailable(googletest)
endif()

# One test runner; list test files explicitly (preferred over GLOB)
add_executable(qss_tests
  test_main.cpp
  test_distributions.cpp      # add more test_*.cpp files as you create them
  test_solver.cpp
  test_math_utils.cpp
  test_chemical_potential.cpp
  test_brent_root.cpp
  test_distribution_total.cpp
)

# C++ standard and warnings (re-use your warnings target if present)
target_compile_features(qss_tests PRIVATE cxx_std_17)
if (TARGET project_warnings)
  target_link_libraries(qss_tests PRIVATE project_warnings)
endif()

# Link your library and GoogleTest
target_link_libraries(qss_tests
  PRIVATE
    QSS::library
    GTest::gtest_main
)

# If tests include generated headers, uncomment:
# target_include_directories(qss_tests PRIVATE ${PROJECT_SOURCE_DIR}/configured)

# If tests use OpenMP directly (usually not needed), uncomment:
# if (ENABLE_OPENMP AND TARGET OpenMP::OpenMP_CXX)
#   target_link_libraries(qss_tests PRIVATE OpenMP::OpenMP_CXX)
# endif()

# Let CTest discover and register tests from TEST() macros
include(GoogleTest)
# Enable test discovery and reuse
gtest_discover_tests(qss_tests
  DISCOVERY_TIMEOUT 30
)

if (ENABLE_COVERAGE)
  include(${CMAKE_SOURCE_DIR}/cmake/CodeCoverage.cmake)
  add_coverage(qss_tests)      # instrument test runner too
endif()
