# LTO.cmake — minimal, safe defaults
function(enable_lto target)
  if (NOT TARGET ${target})
    message(FATAL_ERROR "enable_lto: target '${target}' does not exist")
  endif()

  # Run IPO support check once and cache it
  if (NOT DEFINED QSS_IPO_SUPPORTED)
    include(CheckIPOSupported)
    check_ipo_supported(RESULT QSS_IPO_SUPPORTED OUTPUT QSS_IPO_OUTPUT)
    set(QSS_IPO_SUPPORTED "${QSS_IPO_SUPPORTED}" CACHE BOOL   "IPO/LTO support")
    set(QSS_IPO_OUTPUT    "${QSS_IPO_OUTPUT}"    CACHE STRING "IPO/LTO check output")
  endif()

  if (NOT QSS_IPO_SUPPORTED)
    # Quietly skip when not supported
    return()
  endif()

  # Prefer per-config enabling; keep Debug fast
  set_property(TARGET ${target} PROPERTY INTERPROCEDURAL_OPTIMIZATION FALSE)
  foreach(cfg IN ITEMS RELEASE RELWITHDEBINFO)
    set_property(TARGET ${target} PROPERTY INTERPROCEDURAL_OPTIMIZATION_${cfg} TRUE)
  endforeach()
endfunction()
