# Append this module’s sources to the main library
target_sources(qss_library PRIVATE
  distributions.cpp
  impl.hpp
  fd.cpp
  flb.cpp
  nflb.cpp
  nflb_i.cpp
)
