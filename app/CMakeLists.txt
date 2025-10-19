# Executable
add_executable(qss_cli main.cpp)
add_executable(QSS::cli ALIAS qss_cli)

# C++ standard (target-scoped)
target_compile_features(qss_cli PRIVATE cxx_std_17)

# Link the core library (defined in src/)
target_link_libraries(qss_cli PRIVATE QSS::library)

# Project-wide warnings (optional helper from cmake/Warnings.cmake)
if (TARGET project_warnings)
  target_link_libraries(qss_cli PRIVATE project_warnings)
endif()


# Install the binary when building this project standalone
if (PROJECT_IS_TOP_LEVEL)
  include(GNUInstallDirs)
  install(TARGETS qss_cli
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  )
endif()
