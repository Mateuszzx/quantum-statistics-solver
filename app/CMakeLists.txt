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


# Install the CLI only when explicitly requested. By default OFF to avoid
# writing into system locations like /usr/local/bin when using `cmake --install .`.
option(QSS_INSTALL_CLI "Install the qss_cli executable" OFF)
if (PROJECT_IS_TOP_LEVEL AND QSS_INSTALL_CLI)
  include(GNUInstallDirs)
  install(TARGETS qss_cli
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  )
endif()
