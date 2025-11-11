# ==============================================================================
# Test Suite
# ==============================================================================

include_guard(GLOBAL)

message(STATUS "Processing tests/CMakeLists.txt")

# Add C++ tests subdirectory
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/cpp/CMakeLists.txt")
  message(STATUS "Adding C++ tests from: ${CMAKE_CURRENT_SOURCE_DIR}/cpp")
  add_subdirectory(cpp)
else()
  message(WARNING "C++ tests CMakeLists.txt not found")
endif()

message(STATUS "Test suite configured")
