#include <cmath>
#include <tuple>
#include "qss/distributions.hpp"
#include "qss/math_utils.hpp"
#include "impl.hpp"  // <- private kernels

namespace qss {

// Scalar dispatcher
std::optional<double> distribution(double e_k, double mu, double beta, DistType dist) {
  using namespace qss::detail;
  switch (dist) {
    case DistType::FD:      return fd_distribution   (e_k, mu, beta);
    case DistType::FLB:     return flb_distribution  (e_k, mu, beta);
    case DistType::NFLB:    return nflb_distribution (e_k, mu, beta);
    case DistType::NFLB_i:  return nflb_i_distribution(e_k, mu, beta);
    default:                return std::nullopt;
  }
}

// Vector dispatcher (returns Vec)
std::optional<Vec> distribution(const Vec& e_k, double mu, double beta, DistType dist) {
  using namespace qss::detail;
  const std::size_t n = e_k.size();
  Vec out(n);

  switch (dist) {
    case DistType::FD:      fd_distribution    (e_k.data(), n, mu, beta, out.data()); break;
    case DistType::FLB:     flb_distribution   (e_k.data(), n, mu, beta, out.data()); break;
    case DistType::NFLB:    nflb_distribution  (e_k.data(), n, mu, beta, out.data()); break;
    case DistType::NFLB_i:  nflb_i_distribution(e_k.data(), n, mu, beta, out.data()); break;
    default: return std::nullopt;
  }
  return out;
}

// 2D lattice preparation (unchanged)
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
