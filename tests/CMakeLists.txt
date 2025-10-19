# tests/CMakeLists.txt

# One test runner; list test files explicitly (preferred over GLOB)
add_executable(qss_tests
  test_main.cpp
  test_distributions.cpp      # add more test_*.cpp files as you create them
  # test_chemical_potential.cpp
  # test_solver.cpp
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
gtest_discover_tests(qss_tests
  # WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}
  DISCOVERY_TIMEOUT 30
)
