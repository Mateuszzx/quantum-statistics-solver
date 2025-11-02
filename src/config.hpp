#pragma once
#include <cstddef>

#ifndef QSS_OMP_THRESHOLD_OVERRIDE
  // Default: ~16k elements — good starting point for exp()-heavy loops
  inline constexpr std::size_t QSS_OMP_THRESHOLD = (1u << 14);
#else
  inline constexpr std::size_t QSS_OMP_THRESHOLD = QSS_OMP_THRESHOLD_OVERRIDE;
#endif
