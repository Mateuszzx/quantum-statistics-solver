#include <cmath>
#include <tuple>
#include "qss/distributions.hpp"
#include "qss/math_utils.hpp"

namespace qss {

// Dispatcher (scalar)
std::optional<double> distribution(double e_k, double mu, double beta, DistType dist) {
  switch (dist) {
    case DistType::FD:    return fd_distribution(e_k, mu, beta);
    case DistType::FLB:   return flb_distribution(e_k, mu, beta);
    case DistType::NFL:   return nfl_distribution(e_k, mu, beta);
    case DistType::NFL_i: return nfl_i_distribution(e_k, mu, beta);
    default:              return std::nullopt;
  }
}

// Dispatcher (vector)
std::optional<Vec> distribution(const Vec& e_k, double mu, double beta, DistType dist) {
  switch (dist) {
    case DistType::FD:    return fd_distribution(e_k, mu, beta);
    case DistType::FLB:   return flb_distribution(e_k, mu, beta);
    case DistType::NFL:   return nfl_distribution(e_k, mu, beta);
    case DistType::NFL_i: return nfl_i_distribution(e_k, mu, beta);
    default:              return std::nullopt;
  }
}

// 2D symmetry: (1/(2π)) * |e_k| * f(e_k)
Mat distribution_symmetry_2d(const Mat& e_k, double mu, double beta, DistType dist) {
  const std::size_t rows = e_k.size();
  const std::size_t cols = rows ? e_k[0].size() : 0;
  Mat out = math::make_empty_mat(rows, cols);
  const double pref = 1.0 / (2.0 * PI);

  for (std::size_t j = 0; j < rows; ++j) {
    for (std::size_t i = 0; i < cols; ++i) {
      const double ek = e_k[j][i];
      double val = 0.0;
      switch (dist) {
        case DistType::FD:    val = fd_distribution(ek, mu, beta); break;
        case DistType::FLB:   val = flb_distribution(ek, mu, beta); break;
        case DistType::NFL:   val = nfl_distribution(ek, mu, beta); break;
        case DistType::NFL_i: val = nfl_i_distribution(ek, mu, beta); break;
        default:              val = 0.0; break;
      }
      out[j][i] = pref * std::abs(ek) * val;
    }
  }
  return out;
}

// 2D lattice preparation
std::tuple<Mat, Mat, Mat> prepare_2D_lattice_data(const Vec& kx, const Vec& ky,
                                                  double t, double t_prim) {
  const std::size_t nx = kx.size(), ny = ky.size();
  Mat KX = math::make_empty_mat(ny, nx);
  Mat KY = math::make_empty_mat(ny, nx);
  Mat e_k = math::make_empty_mat(ny, nx);

  for (std::size_t j = 0; j < ny; ++j) {
    for (std::size_t i = 0; i < nx; ++i) {
      const double kxi = kx[i], kyj = ky[j];
      KX[j][i] = kxi;
      KY[j][i] = kyj;
      e_k[j][i] = 2.0 * t * (std::cos(kxi) + std::cos(kyj))
                + 2.0 * t_prim * (std::cos(kxi - kyj) + std::cos(kxi + kyj));
    }
  }
  return {e_k, KX, KY};
}

} // namespace qss
