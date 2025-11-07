add_library(qss_library)
add_library(QSS::library ALIAS qss_library)

target_compile_features(qss_library PUBLIC cxx_std_17)

target_include_directories(qss_library
  PUBLIC  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
          $<INSTALL_INTERFACE:include>
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}
)

if (TARGET project_warnings)
  target_link_libraries(qss_library PRIVATE project_warnings)
endif()

if (ENABLE_OPENMP)
  find_package(OpenMP REQUIRED)
  target_link_libraries(qss_library PUBLIC OpenMP::OpenMP_CXX)
endif()

# Optional but recommended for good ABI hygiene
set_target_properties(qss_library PROPERTIES
  POSITION_INDEPENDENT_CODE ON
  CXX_VISIBILITY_PRESET hidden
  VISIBILITY_INLINES_HIDDEN YES
)

# Public headers (FILE_SET)
target_sources(qss_library PUBLIC
  FILE_SET public_headers TYPE HEADERS
  BASE_DIRS ${PROJECT_SOURCE_DIR}/include
  FILES
    ${PROJECT_SOURCE_DIR}/include/qss/chemical_potential.hpp
    ${PROJECT_SOURCE_DIR}/include/qss/distributions.hpp
    ${PROJECT_SOURCE_DIR}/include/qss/solver.hpp
)

# Private sources shared across modules
target_sources(qss_library PRIVATE
  math_utils.cpp
)

# Modules
add_subdirectory(chemical_potential)
add_subdirectory(distributions)
add_subdirectory(solver)
add_subdirectory(bindings)
