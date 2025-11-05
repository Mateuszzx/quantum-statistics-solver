#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/vector.h>
#include <cstring>
#ifdef _OPENMP
  #include <omp.h>
#endif
#include "qss/types.hpp"
#include "qss/distributions.hpp"
#include "qss/solver.hpp"
#include "qss/solver_params.hpp"
#include "qss/chemical_potential.hpp"
#include "bindings.hpp"

namespace nb = nanobind;
using nb::arg;

using Vec1 = nb::ndarray<nb::numpy, const double, nb::ndim<1>, nb::c_contig>;
using Mat2 = nb::ndarray<nb::numpy, const double, nb::ndim<2>, nb::c_contig>;
using Out4D = nb::ndarray<nb::numpy, double, nb::ndim<4>>;

static qss::SolverParams make_params(double tol, int max_iters,
                                     double expand_factor, int max_expands,
                                     double a0, double b0) {
  qss::SolverParams p{};
  p.tol = tol; p.max_iters = max_iters;
  p.expand_factor = expand_factor; p.max_expands = max_expands;
  p.a0 = a0; p.b0 = b0;
  return p;
}

static Out4D make_out_4d(size_t d0, size_t d1, size_t d2, size_t d3) {
  const size_t N = d0*d1*d2*d3;
  double* data = new double[N];
  nb::capsule owner(data, [](void* p) noexcept { delete[] (double*)p; });
  return Out4D(data, {d0,d1,d2,d3}, owner, { (int64_t)(d1*d2*d3),
                                             (int64_t)(d2*d3),
                                             (int64_t)d3,
                                             1 });
}

namespace qss_bindings {

void bind_chemical_potential(nb::module_& m) {
  m.def("chemical_potential_B_grid",
    [](std::vector<qss::DistType> distros,
       Vec1 n_assumed, Vec1 beta, Vec1 H_a_T,
       Vec1 kx, Vec1 ky, Mat2 E,
       double g, double muB_eV_per_T,
       double tol, int max_iters,
       double expand_factor, int max_expands,
       double a0, double b0,
       int threads) {

      const size_t nx = (size_t)kx.shape(0);
      const size_t ny = (size_t)ky.shape(0);
      if ((size_t)E.shape(0) != ny || (size_t)E.shape(1) != nx)
        throw std::runtime_error("E must have shape (len(ky), len(kx))");

      const size_t ND = distros.size();
      const size_t NN = (size_t)n_assumed.shape(0);
      const size_t NB = (size_t)beta.shape(0);
      const size_t NH = (size_t)H_a_T.shape(0);
      if (!ND || !NN || !NB || !NH)
        throw std::runtime_error("distros, n_assumed, beta, H_a_T must be non-empty");

      // Convert grid arrays once
      qss::Vec vx(nx), vy(ny);
      std::memcpy(vx.data(), kx.data(), nx * sizeof(double));
      std::memcpy(vy.data(), ky.data(), ny * sizeof(double));
      qss::Mat Em(ny, qss::Vec(nx));
      for (size_t j = 0; j < ny; ++j) {
        const double* row = E.data() + j * nx;
        std::memcpy(Em[j].data(), row, nx * sizeof(double));
      }

      auto out = make_out_4d(ND, NN, NB, NH);
      double* __restrict P = out.data();

      qss::SolverParams params =
        make_params(tol, max_iters, expand_factor, max_expands, a0, b0);

      const size_t NBxNH = NB * NH;

      nb::gil_scoped_release release;
      #ifdef _OPENMP
        if (threads > 0) omp_set_num_threads(threads);
        #pragma omp parallel for collapse(3) schedule(static)
      #endif
      for (ptrdiff_t id = 0; id < (ptrdiff_t)ND; ++id)
      for (ptrdiff_t inx= 0; inx < (ptrdiff_t)NN; ++inx)
      for (ptrdiff_t ib = 0; ib < (ptrdiff_t)NB; ++ib) {
        const double beta_v = beta.data()[ib];
        const double n_v    = n_assumed.data()[inx];
        const size_t base   = (size_t)id * (NN*NBxNH)
                            + (size_t)inx *  NBxNH
                            + (size_t)ib  *  NH;

        for (ptrdiff_t ih = 0; ih < (ptrdiff_t)NH; ++ih) {
          const double mu = qss::chemical_potential_B(
              distros[(size_t)id], beta_v, n_v,
              vx, vy, Em,
              H_a_T.data()[ih], g, muB_eV_per_T, params);
          P[base + (size_t)ih] = mu;
        }
      }
      return out;
    },
    arg("distros"), arg("n_assumed"), arg("beta"), arg("H_a_T"),
    arg("kx"), arg("ky"), arg("E"),
    arg("g") = 2.0,
    arg("muB_eV_per_T") = qss::QSS_MU_B_EV_PER_T,
    arg("tol") = 1e-6, arg("max_iters") = 300,
    arg("expand_factor") = 1.5, arg("max_expands") = 200,
    arg("a0") = -10.0, arg("b0") = 10.0,
    arg("threads") = 0
  );
}

} // namespace qss_bindings
