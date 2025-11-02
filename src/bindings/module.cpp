#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include "qss/types.hpp"
#include "qss/distributions.hpp"
#include "qss/solver.hpp"               // for factories (Brent + expansion)
#include "qss/math_utils.hpp"           // mean_abs_2d if needed in future
// chemical_potential_2d_default is implemented in your chemical_potential.cpp
// (per earlier step)
namespace nb = nanobind;
using namespace nb::literals;

NB_MODULE(qss, m) {
  // ---- Enum DistType ----
  nb::enum_<qss::DistType>(m, "DistType")
      .value("FD",       qss::DistType::FD)
      .value("FLB",      qss::DistType::FLB)
      .value("NFLB",     qss::DistType::NFLB)
      .value("NFLB_i",   qss::DistType::NFLB_i);

  // ---- prepare_2d_lattice ----
  m.def("prepare_2d_lattice",
        [](nb::ndarray<const double, nb::c_contig> kx,
           nb::ndarray<const double, nb::c_contig> ky,
           double t, double t_prim) {
          // Copy kx/ky into std::vector (cheap vs sizes)
          qss::Vec vx(kx.size()), vy(ky.size());
          std::memcpy(vx.data(), kx.data(), kx.size() * sizeof(double));
          std::memcpy(vy.data(), ky.data(), ky.size() * sizeof(double));

          auto [E, KX, KY] = qss::prepare_2D_lattice_data(vx, vy, t, t_prim);

          // Allocate NumPy arrays and copy out
          const std::size_t ny = KY.size();
          const std::size_t nx = (ny ? KY[0].size() : 0);

          nb::ndarray<double, nb::c_contig> e_py({ny, nx});
          nb::ndarray<double, nb::c_contig> kx_py({ny, nx});
          nb::ndarray<double, nb::c_contig> ky_py({ny, nx});

          for (std::size_t j = 0; j < ny; ++j) {
            std::memcpy(e_py.data()  + j*nx, E [j].data(),  nx * sizeof(double));
            std::memcpy(kx_py.data() + j*nx, KX[j].data(), nx * sizeof(double));
            std::memcpy(ky_py.data() + j*nx, KY[j].data(), nx * sizeof(double));
          }
          return nb::make_tuple(e_py, kx_py, ky_py);
        },
        "kx"_a, "ky"_a, "t"_a, "t_prim"_a,
        "Return (E, KX, KY) as NumPy arrays with shape (len(ky), len(kx)).");

  // ---- distribution_df ----
  // Build a pandas.DataFrame from E and f(E; mu, beta, dist)
  m.def("distribution_df",
        [](nb::ndarray<const double, nb::c_contig> E,
           double mu, double beta, qss::DistType dist,
           nb::object kx_opt = nb::none(), nb::object ky_opt = nb::none()) {
          // E is (ny, nx)
          if (E.ndim() != 2)
            throw std::runtime_error("distribution_df: E must be 2D array");

          const std::size_t ny = E.shape(0);
          const std::size_t nx = E.shape(1);

          // Compute F row-by-row using public dispatcher
          nb::ndarray<double, nb::c_contig> F({ny, nx});
          for (std::size_t j = 0; j < ny; ++j) {
            const double* rowE = E.data() + j * nx;
            qss::Vec row_vec(rowE, rowE + nx);
            auto frow_opt = qss::distribution(row_vec, mu, beta, dist);
            if (!frow_opt) throw std::runtime_error("distribution() failed");
            std::memcpy(F.data() + j * nx, frow_opt->data(), nx * sizeof(double));
          }

          // Build pandas.DataFrame(F, index=ky?, columns=kx?)
          nb::module_ pd = nb::module_::import_("pandas");
          if (!kx_opt.is_none() && !ky_opt.is_none()) {
            return pd.attr("DataFrame")(F,
                                        "index"_a=ky_opt,
                                        "columns"_a=kx_opt);
          } else {
            return pd.attr("DataFrame")(F);
          }
        },
        "E"_a, "mu"_a, "beta"_a, "dist"_a,
        "kx"_a = nb::none(), "ky"_a = nb::none(),
        "Return pandas.DataFrame of f(E). Optionally label index/columns with ky/kx.");

  // ---- chemical_potential_df ----
  // Accept beta scalar or 1D array. Return DataFrame with columns ['beta','mu'].
  m.def("chemical_potential_df",
        [](qss::DistType dist,
           nb::object beta_obj,  // float or 1D array
           double n_assumed,
           nb::ndarray<const double, nb::c_contig> kx,
           nb::ndarray<const double, nb::c_contig> ky,
           nb::ndarray<const double, nb::c_contig> E,
           // optional knobs
           double a0, double b0,
           double tol, int max_iter,
           double expand_factor, int max_expands) {
          // Convert kx, ky to vectors for size checks / labels
          qss::Vec vx(kx.size()), vy(ky.size());
          std::memcpy(vx.data(), kx.data(), kx.size() * sizeof(double));
          std::memcpy(vy.data(), ky.data(), ky.size() * sizeof(double));

          // Sanity on E shape
          if (E.ndim() != 2)
            throw std::runtime_error("chemical_potential_df: E must be 2D array");
          const std::size_t ny = E.shape(0);
          const std::size_t nx = E.shape(1);
          if (ny != vy.size() || nx != vx.size())
            throw std::runtime_error("chemical_potential_df: shape(E) must be (len(ky), len(kx))");

          // Build Mat E_mat (row-major copy)
          qss::Mat E_mat(ny, qss::Vec(nx));
          for (std::size_t j = 0; j < ny; ++j)
            std::memcpy(E_mat[j].data(), E.data() + j*nx, nx * sizeof(double));

          // Prepare solver pieces (Expansion bracketer + Brent)
          auto br = qss::make_bracketer(qss::BracketerType::Expansion);
          auto sl = qss::make_solver(qss::SolverType::Brent);
          qss::BracketerOptions bopt{expand_factor, max_expands};
          qss::RootOptions      ropt{tol, max_iter};

          // Gather betas
          std::vector<double> betas;
          if (nb::isinstance<nb::float_>(beta_obj)) {
            betas.push_back(nb::cast<double>(beta_obj));
          } else if (nb::isinstance<nb::ndarray<nb::c_contig>>(beta_obj)) {
            nb::ndarray<const double, nb::c_contig> b_arr = nb::cast<nb::ndarray<const double, nb::c_contig>>(beta_obj);
            if (b_arr.ndim() != 1)
              throw std::runtime_error("chemical_potential_df: beta must be scalar or 1D array");
            betas.resize(b_arr.size());
            std::memcpy(betas.data(), b_arr.data(), b_arr.size() * sizeof(double));
          } else {
            throw std::runtime_error("chemical_potential_df: beta must be float or 1D float array");
          }

          // Compute mus
          std::vector<double> mus; mus.reserve(betas.size());
          for (double beta : betas) {
            double mu = qss::chemical_potential_2d(
                dist, beta, n_assumed, vx, vy, E_mat,
                a0, b0, *br, bopt, *sl, ropt);
            mus.push_back(mu);
          }

          // Build pandas.DataFrame({'beta': betas, 'mu': mus})
          nb::module_ pd = nb::module_::import_("pandas");
          nb::dict cols;
          // Convert vectors to NumPy arrays
          nb::ndarray<double> beta_np(betas.size());
          std::memcpy(beta_np.data(), betas.data(), betas.size()*sizeof(double));
          nb::ndarray<double> mu_np(mus.size());
          std::memcpy(mu_np.data(), mus.data(), mus.size()*sizeof(double));

          cols["beta"] = beta_np;
          cols["mu"]   = mu_np;
          return pd.attr("DataFrame")(cols);
        },
        "dist"_a,
        "beta"_a,                   // float or 1D np.ndarray
        "n_assumed"_a,
        "kx"_a, "ky"_a, "E"_a,      // numpy arrays
        "a0"_a = -5000.0, "b0"_a = 1000.0,
        "tol"_a = 1e-6, "max_iter"_a = 300,
        "expand_factor"_a = 1.5, "max_expands"_a = 200,
        "Return pandas.DataFrame with columns ['beta','mu'].");
}
